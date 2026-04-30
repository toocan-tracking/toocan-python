import sys
import numpy as np 
import netCDF4
import pandas as pd
import xarray as xr
import scipy as sp
import glob
import random
import itertools
import datetime
import re
import os
from os import listdir
from os.path import isfile, join
import math
import time
import ctypes
import argparse
import subprocess
import gzip
from tqdm import tqdm


########################################################################################################################################

# python3 Write_Filelist.py

########################################################################################################################################

class launchParameters(object): 
  
    def __init__(self):
        self.region             = 0
        self.lonmin             = 0
        self.lonmax             = 0
        self.latmin             = 0
        self.latmax             = 0
        self.path_OUT           = 0
        self.path_imagesgeo     = 0
        self.nomenclature_image = 0
        self.reso_tempo			= 0

def load_launchparam(file):
	"""
	Load launch parameters from a text file.
	"""
    
    launchparam=launchParameters()

    lunit=open(file,'rt')
    print ("extract launch parameters from ",file)
    lines = lunit.readlines()
    
    launchparam.region               = lines[0].split()[-1]
    launchparam.lonmin               = lines[1].split()[-1]
    launchparam.lonmax               = lines[2].split()[-1]
    launchparam.latmin               = lines[3].split()[-1]
    launchparam.latmax               = lines[4].split()[-1]
    launchparam.path_OUT             = lines[5].split()[-1]
    launchparam.path_imagesgeo       = lines[6].split()[-1]
    launchparam.nomenclature_image   = lines[7].split()[-1]
    launchparam.reso_tempo			 = lines[8].split()[-1]

    return (launchparam)

######################################################################################################
######################################################################################################
# print("\n--- Writing filelist TOOCAN ---")

params_launch = sys.argv[1]
launchparam   = load_launchparam(params_launch)

year_input   = sys.argv[2]
month_input  = sys.argv[3]
day_input    = sys.argv[4]
hour_input   = sys.argv[5]
minute_input = sys.argv[6]

year_input_begin   = sys.argv[7]
month_input_begin  = sys.argv[8]
day_input_begin    = sys.argv[9]
hour_input_begin   = sys.argv[10]
minute_input_begin = sys.argv[11]

output_path = sys.argv[12]

reso_tempo = int(launchparam.reso_tempo)
# reso_tempo = 15  ##### pour MAESTRO
# reso_tempo = 30  ##### pour test

######################################################################################################
######################################################################################################


dateEND   = datetime.datetime(int(year_input), int(month_input), int(day_input), int(hour_input), int(minute_input))

dateBEGIN   = datetime.datetime(int(year_input_begin), int(month_input_begin), int(day_input_begin), int(hour_input_begin), int(minute_input_begin))

delta_resotempo = datetime.timedelta(minutes=float(reso_tempo))
nb_step = ((dateEND - dateBEGIN) // delta_resotempo) + 1 # number of times reso_tempo between the two intervals

time = dateEND


# --------------------------------------------------------------------------------------------------------------- #

prev_month = int(month_input)-1
if int(month_input)==1 : prev_month=12
if len(str(prev_month))==1: prev_month = '0' + str(prev_month)


path_IR = f'{launchparam.path_imagesgeo}/'
files_IR = sorted(glob.glob(path_IR+f'*/{launchparam.nomenclature_image}*.nc'))    # adjust according to the IR file paths
print("files_IR created")

for ifile_IR in files_IR :	# l.139 need ifile_IR
	continue


year       = np.zeros(nb_step,dtype=('i4'))
month      = np.zeros(nb_step,dtype=('i4'))
day        = np.zeros(nb_step,dtype=('i4'))
hour       = np.zeros(nb_step,dtype=('i4'))
minute     = np.zeros(nb_step,dtype=('i4'))
julianday  = np.zeros(nb_step,dtype=('i4'))
slot       = np.zeros(nb_step,dtype=('i4'))
exists     = np.zeros(nb_step,dtype=('i4'))
missing    = np.zeros(nb_step,dtype=('i8'))
state      = np.zeros(nb_step,dtype=('i4'))
mode       = np.zeros(nb_step,dtype=('i4'))
file       = ["" for i in range(nb_step)]
path       = ["" for i in range(nb_step)]



_missing   = 0
flagdeb    = 0
nb_issue   = 0

# invalid_files_path = '/home/jbounan/toocan-python/FileINPUT/invalid_files_2008_2024.txt'   # spatially invalid files
# with open(invalid_files_path) as f:
#     invalid_files = [line.split() for line in f.readlines()]
invalid_files = []


file_by_time = {}

pattern = re.compile(r'_(\d{4})-(\d{2})-(\d{2})T(\d{2})-(\d{2})-(\d{2})')

for f in files_IR:
    fname = f.split('/')[-1]
    m = pattern.search(fname)
    if m is None:
        continue

    year_i, month_i, day_i, HH_i, MM_i, _ = m.groups()
    t = datetime.datetime(
        int(year_i), int(month_i), int(day_i),
        int(HH_i), int(MM_i)
    )

    file_by_time[t] = f


for itime in tqdm(range(0, nb_step, 1)):

	year[itime]      = time.year
	month[itime]     = time.month
	day[itime]       = time.day
	hour[itime]      = time.hour
	minute[itime]    = time.minute
	julianday[itime] = (time - datetime.datetime(time.year, 1, 1)).days + 1
	slot[itime]      = time.hour * (60 // reso_tempo) + 1 + time.minute // reso_tempo

	path[itime] = 'N/A'
	file[itime] = 'N/A'
	flagOK = 0
	
	if time in file_by_time:
		f = file_by_time[time]
		fname = f.split('/')[-1]

		if f not in invalid_files:

			iexist = 1
			imissing = 0

			path[itime] = f[:ifile_IR.rfind('/')+1]
			
			if flagdeb == 0:
				state[itime] = 1
				flagdeb = 1

				flagOK         = 1 
				exists[itime]  = iexist       # images exists
				missing[itime] = imissing     
				mode[itime]    = 0            # for the GEOS mode
				file[itime]    = fname   # file name
				
				
			else:                                   # find the image corresponding to the time, and it is not the last image for the run
				state[itime]       = 0
				exists[itime]      = iexist
				missing[itime]     = imissing
				_missing           = 0
				flagOK             = 1
				mode[itime]        = 0
				file[itime]        = fname

			
	if(flagOK==0):
		print(f"missing {time}")
		_missing = _missing +1
		missing[itime] = _missing
		iexist = 0
		exists[itime] = iexist
		# flagdeb = 0
		
	time -= delta_resotempo


# Start of the run: state = 2 for the oldest non-missing image in the list
for i in np.arange(1, nb_step+1, 1):
	if exists[-i] == 1 :
		state[-i] = 2	
		break
	else:
		continue


# set state = 1 before a batch of 5 or more missing images
for i in np.arange(1, nb_step, 1):
	if missing[i-1] >= 5 and missing[i]==0 :
		state[i] = 1


# Handle batches of 5 or more missing images to determine the value of state after the gap; start of a new run, so state = 2
flag_run_cut = 0
flag_reprise = nb_step-1
for i in np.arange(nb_step-1, -1, -1):
	if missing[i] >= 5 :
		flag_run_cut = 1
	if flag_run_cut == 0 :
		continue
	if missing[i] == 0 :
		state[i] = 2
		flag_reprise = i # use i as a flag to write to the file only from the point where it resumes
		flag_run_cut = 0


# fill in missing file and path values with the nearest existing value
for itime in np.arange(nb_step-2, -1, -1):
	if missing[itime] > 0 :
		if missing[itime + 1] == 0:
			nb_missing = missing[itime]
			paths_toto = [path[itime + 1], path[itime - nb_missing]]
			files_toto = [file[itime + 1], file[itime - nb_missing]]
			path[itime - (nb_missing + 1) // 2 + 1 : itime + 1] = [paths_toto[0]] * ((nb_missing + 1) // 2)
			file[itime - (nb_missing + 1) // 2 + 1 : itime + 1] = [files_toto[0]] * ((nb_missing + 1) // 2)
			path[itime - nb_missing + 1 : itime - (nb_missing + 1) // 2 + 1] = [paths_toto[1]] * (nb_missing - (nb_missing + 1) // 2)
			file[itime - nb_missing + 1 : itime - (nb_missing + 1) // 2 + 1] = [files_toto[1]] * (nb_missing - (nb_missing + 1) // 2)

		else:
			continue



df = pd.DataFrame({
    "year": year,
    "month": month,
    "day": day,
    "hour": hour,
    "minute": minute,
    "julian": julianday,
    "slot": slot,
    "exists": exists,
    "missing": missing,
    "state": state,
    "mode": mode,
    "path_ir": path,
    "file_ir": file,
})

df["time"] = pd.to_datetime(
    dict(
        year=df.year,
        month=df.month,
        day=df.day,
        hour=df.hour,
        minute=df.minute
    )
)
df = df.set_index("time")

# Save in NetCDF format
ds = xr.Dataset.from_dataframe(df)

ds.to_netcdf(
    output_path
)
