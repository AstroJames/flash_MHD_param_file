#!/usr/bin/env python

"""
	written by James Beattie, 2021.
	Used for automatically setting up MHD mean-field turbulence simulations in FLASH
	> sets the dump times for the desired Mach
	> sets the mean magnetic field for the desired Alfven Mach number
	> adds the forcing file to the flash.par file

	############################################################################
	Instructions:
	add to ~/bin/
	chmod +x prep_MHD

	Exmaple useage:
	prep_MHD -m 10 -ma 1.0 -forcing forcingfile_M10.dat

	for mach = 10, m_A0 = 1.0, forcing file = forcingfile_M10.dat

	############################################################################

	based on prep_restart.py, written by Christoph Federrath
"""

import argparse
import numpy as np
import os
import sys
import time
import array
import re
import subprocess
import fnmatch
from tempfile import mkstemp
from shutil import move, copyfile
from os import remove, close, path

ap = argparse.ArgumentParser(description='command line inputs')
ap.add_argument('-f', '--file_name',default='flash.par',help='the name of the flash parameter file', type=str)
ap.add_argument('-m', '--mach',default=None,required=True,help='the rms mach number used to set the runtime', type=float)
ap.add_argument('-d', '--driving_scale',default=2.0,required=False,help='the driving scale of the turbulence', type=float)
ap.add_argument('-ma', '--mach_alfven',default=None,required=True,help='the mach Alfven w.r.t the mean B-field', type=float)
ap.add_argument('-n', '--dumps_per_turnover',default=10.0,required=False,help='the number of dumps per turnover time', type=float)
ap.add_argument('-T', '--total_turn_overs',default=10.0,required=False,help='the total number of turnover times', type=float)
ap.add_argument('-forcing', '--forcing_file',default=None,required=False,help='the turbulent forcing file', type=str)
args = vars(ap.parse_args())

class UpdateFlashParameterFile:

	def __init__(self,
				 mach_in,
				 mach_alfven_in,
				 driving_scale_in,
				 dumps_per_turnover_in,
				 total_turn_overs_in,
				 forcing_file_in,
				 file_name_in):
		self.mach 						= mach_in
		self.mach_alfven 				= mach_alfven_in
		self.driving_scale 				= driving_scale_in
		self.dumps_per_turnover			= dumps_per_turnover_in
		self.total_turn_overs 			= total_turn_overs_in
		self.file_name 					= file_name_in
		self.forcing_file 				= forcing_file_in
		self.debug 						= True
		self.B0 						= []
		self.tmax 						= []
		self.checkpointFileIntervalTime = []
		self.plotFileIntervalTime 		= []
		self.movie_dt_dump 				= []
		self.particleFileIntervalTime 	= []
		self.dtmax 						= []

	def compute_B0(self):
		rho_0 = 1 	# assumes rho_0 = 1
		c_s = 1		# assumes the sound speed = 1
		self.B0 = 2*c_s*np.sqrt(rho_0*np.pi) * self.mach / self.mach_alfven

	def compute_turnover_time_parameters(self):
		self.tmax = self.total_turn_overs / (self.driving_scale * self.mach) 							# the maxmimum time domain
		self.checkpointFileIntervalTime = self.tmax / (self.total_turn_overs * self.dumps_per_turnover)	# the dump time for checkpoints
		self.plotFileIntervalTime = self.checkpointFileIntervalTime										# the dump time for plot files
		self.movie_dt_dump = self.plotFileIntervalTime / 10												# the dump time for movie files
		self.particleFileIntervalTime = self.movie_dt_dump												# the dump time for particles files
		self.dtmax = self.movie_dt_dump / 2																# the maximum integration time step

	def replace_magnetic_field(self,line):
		# compute the desired B0 for mach, ma0
		self.compute_B0()
		# replace magnetic field
		if line.find("MagField_z")==0 or line.find("st_MPzBmeanTarget")==0:
			if self.debug == True:
				print(f"{self.file_name}: found line:\t\t{line.rstrip()}")
			i = line.find("=")
			newline = line[0:i+1]+ f" {self.B0}\n"
			line = newline
			if self.debug == True:
				print(f"{self.file_name}: replaced with:\t{line.rstrip()}")
		return line

	def replace_forcing_file(self,line):
		# replace forcing file
		if line.find("st_infilename")==0:
			if self.debug == True:
				print(f"{self.file_name}: found line:\t\t{line.rstrip()}")
			i = line.find("=")
			newline = line[0:i+1]+ f" {self.forcing_file}\n"
			line = newline
			if self.debug == True:
				print(f"{self.file_name}: replaced with:\t{line.rstrip()}")
		return line

	def replace_turnover_times(self,line):
		# compute all of the turnover things
		self.compute_turnover_time_parameters()
		# replace turnover time parameters
		if line.find("movie_dt_dump")==0:
			if self.debug == True:
				print(f"{self.file_name}: found line:\t\t{line.rstrip()}")
			i = line.find("=")
			newline = line[0:i+1]+ f" {self.movie_dt_dump}\n"
			line = newline
			if self.debug == True:
				print(f"{self.file_name}: replaced with:\t{line.rstrip()}")
		elif line.find("checkpointFileIntervalTime")==0:
			if self.debug == True:
				print(f"{self.file_name}: found line:\t\t{line.rstrip()}")
			i = line.find("=")
			newline = line[0:i+1]+ f" {self.checkpointFileIntervalTime}\n"
			line = newline
			if self.debug == True:
				print(f"{self.file_name}: replaced with:\t{line.rstrip()}")
		elif line.find("plotFileIntervalTime")==0:
			if self.debug == True:
				print(f"{self.file_name}: found line:\t\t{line.rstrip()}")
			i = line.find("=")
			newline = line[0:i+1]+ f" {self.plotFileIntervalTime}\n"
			line = newline
			if self.debug == True:
				print(f"{self.file_name}: replaced with:\t{line.rstrip()}")
		elif line.find("particleFileIntervalTime")==0:
			if self.debug == True:
				print(f"{self.file_name}: found line:\t\t{line.rstrip()}")
			i = line.find("=")
			newline = line[0:i+1]+ f" {self.particleFileIntervalTime}\n"
			line = newline
			if self.debug == True:
				print(f"{self.file_name}: replaced with:\t{line.rstrip()}")
		elif line.find("dtmax")==0:
			if self.debug == True:
				print(f"{self.file_name}: found line:\t\t{line.rstrip()}")
			i = line.find("=")
			newline = line[0:i+1]+ f" {self.dtmax}\n"
			line = newline
			if self.debug == True:
				print(f"{self.file_name}: replaced with:\t{line.rstrip()}")
		elif line.find("tmax")==0:
			if self.debug == True:
				print(f"{self.file_name}: found line:\t\t{line.rstrip()}")
			i = line.find("=")
			newline = line[0:i+1]+ f" {self.tmax}\n"
			line = newline
			if self.debug == True:
				print(f"{self.file_name}: replaced with:\t{line.rstrip()}")
		return line

	def update_flash_par(self):
		fh, tempfile = mkstemp()
		ftemp = open(tempfile, 'w')
		f = open(self.file_name, 'r')
		for line in f:
			# replace magnetic field
			line = self.replace_magnetic_field(line)
			# replace turnover time data based on Mach
			line = self.replace_turnover_times(line)
			# replace turbulent driving
			if self.forcing_file is not None:
				line = self.replace_forcing_file(line)
			# write the new lines
			ftemp.write(line)
		ftemp.close()
		close(fh)
		f.close()
		remove(self.file_name)
		move(tempfile, self.file_name)
		os.chmod(self.file_name, 644)

def main():

	# extract command line arguments
	file_name 			= args["file_name"]
	mach 				= args["mach"]
	mach_alfven 		= args["mach_alfven"]
	driving_scale 		= args["driving_scale"]
	dumps_per_T 		= args["dumps_per_turnover"]
	total_T 			= args["total_turn_overs"]
	forcing_file 		= args["forcing_file"]

	# make a backup copy of flash.par
	print("Copying 'flash.par' to 'flash.par_restart_backup' as backup.")
	copyfile(file_name,"flash.par_restart_backup")

	# read in the flash parameters
	flash_paramters = UpdateFlashParameterFile(mach,
											   mach_alfven,
											   driving_scale,
											   dumps_per_T,
											   total_T,
											   forcing_file,
											   file_name)

	# udpate the flash.par parameter file
	flash_paramters.update_flash_par()

if __name__ == '__main__':
	main()
