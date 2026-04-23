import sys
import numpy as np 
import netCDF4
import pandas as pd
import xarray as xr
import scipy as sp
import glob
import random
import itertools
# from datetime import date, time, datetime, timedelta
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

reso_tempo = int(launchparam.reso_tempo)
# reso_tempo = 15  ##### pour MAESTRO
# reso_tempo = 30  ##### pour test

######################################################################################################
######################################################################################################


dateEND   = datetime.datetime(int(year_input), int(month_input), int(day_input), int(hour_input), int(minute_input))

# dateBEGIN = dateEND - datetime.timedelta(days=365) + datetime.timedelta(minutes=reso_tempo)
dateBEGIN   = datetime.datetime(int(year_input_begin), int(month_input_begin), int(day_input_begin), int(hour_input_begin), int(minute_input_begin))

delta_resotempo = datetime.timedelta(minutes=float(reso_tempo))
nombre_de_pas = ((dateEND - dateBEGIN) // delta_resotempo) + 1 # nb de fois reso_tempo entre les deux intervales

time = dateEND


# --------------------------------------------------------------------------------------------------------------- #

prev_month = int(month_input)-1
if int(month_input)==1 : prev_month=12
if len(str(prev_month))==1: prev_month = '0' + str(prev_month)

# list_years = [y for y in range(int(year_input_begin), int(year_input) + 1)]

path_IR = f'{launchparam.path_imagesgeo}/'
# files_IR = sorted(glob.glob(path_IR+f'/{year_input}_{prev_month}*/{launchparam.nomenclature_image}*.nc')) + sorted(glob.glob(path_IR+f'/{year_input}_{month_input}*/{launchparam.nomenclature_image}*.nc'))
# files_IR = sorted(glob.glob(path_IR+f'/*/*/{launchparam.nomenclature_image}*.nc')) + sorted(glob.glob(path_IR+f'/{year_input}_{month_input}*/{launchparam.nomenclature_image}*.nc'))
# files_IR = sorted(glob.glob(path_IR+f'{year_input}/*/*/{launchparam.nomenclature_image}*.nc'))
files_IR = sorted(glob.glob(path_IR+f'*/{launchparam.nomenclature_image}*.nc'))
print("files_IR created")
##### pour test
# path_IR    = f'/bdd/MT_WORKSPACE/lgouttes/MEGHA_TROPIQUES/MSG+0000/GRID004.v2.00/'
# files_IR    = sorted(glob.glob(path_IR+f'/*/*.nc'))

for ifile_IR in files_IR :	# l.139 need ifile_IR
	# print(ifile_IR)
	continue


year       = np.zeros(nombre_de_pas,dtype=('i4'))
month      = np.zeros(nombre_de_pas,dtype=('i4'))
day        = np.zeros(nombre_de_pas,dtype=('i4'))
hour       = np.zeros(nombre_de_pas,dtype=('i4'))
minute     = np.zeros(nombre_de_pas,dtype=('i4'))
julianday  = np.zeros(nombre_de_pas,dtype=('i4'))
slot       = np.zeros(nombre_de_pas,dtype=('i4'))
exists     = np.zeros(nombre_de_pas,dtype=('i4'))
missing    = np.zeros(nombre_de_pas,dtype=('i8'))
state      = np.zeros(nombre_de_pas,dtype=('i4'))
mode       = np.zeros(nombre_de_pas,dtype=('i4'))
file       = ["" for i in range(nombre_de_pas)]
path       = ["" for i in range(nombre_de_pas)]



_missing   = 0
flagdeb    = 0
nb_issue   = 0

# invalid_files_path = '/home/jbounan/toocan-python/FileINPUT/invalid_files_2008_2024.txt'   # fichiers invalides spatialement
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


for itime in tqdm(range(0, nombre_de_pas, 1)):

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

			# print(time_i, '=', time, flagdeb)
			iexist = 1
			imissing = 0

			path[itime] = f[:ifile_IR.rfind('/')+1]
			
			if flagdeb == 0:
				state[itime] = 1
				flagdeb = 1

				flagOK         = 1 
				exists[itime]  = iexist       # l'image existe
				missing[itime] = imissing     
				mode[itime]    = 0            # pour les mode de GEOS
				file[itime]    = fname   # nom du fichier
				
				
			else:                                   # on trouve l'image qui correspond au temps et il ne s'agit pas de la derniere image pour le run
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


# for itime in tqdm(np.arange(0, nombre_de_pas, 1)):
	

# 	flagOK  = 0
	
# 	year[itime]      = time.year
# 	month[itime]     = time.month
# 	day[itime]       = time.day
# 	hour[itime]      = time.hour
# 	minute[itime]    = time.minute
# 	julianday[itime] = (time-datetime.datetime(year[itime] , 1, 1,0,)).days+1
# 	slot[itime]      = time.hour*(int(60/reso_tempo)) + 1 + int(time.minute/reso_tempo) # reso 15min

# 	path[itime]      ='N/A'
# 	file[itime]      ='N/A'

# 	for itime_files in range(len(files_IR)):
		
# 		ifile_name = files_IR[itime_files].split('/')[-1]
# 		#print(ifile_name)
		
# 		##### exemple ifile_name pour MAESTRO : '/bdd/msg/images/netcdf/fulldisk_00/2024/2024_01_08/Mmultic3kmNC4_msg03_202401080000.nc'
# 		ismatch = re.search(r'_(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})\.nc', ifile_name)
# 		##### exemple ifile_name pour test    : '/bdd/MT_WORKSPACE/lgouttes/MEGHA_TROPIQUES/MSG+0000/GRID004.v2.00/2020_01_31/GEO_L1C-MSG4_2020-01-31T23-30-00_G_IR108_004_V1.1.nc'
# 		#ismatch = re.search(r'_(\d{4})-(\d{2})-(\d{2})T(\d{2})-(\d{2})-(\d{2})', ifile_name)
		
# 		year_i  = ismatch.group(1)
# 		month_i = ismatch.group(2)
# 		day_i   = ismatch.group(3)
# 		HH_i    = ismatch.group(4)
# 		MM_i    = ismatch.group(5)
		
# 		time_i = datetime.datetime(int(year_i), int(month_i), int(day_i), int(HH_i), int(MM_i))

# 		if(time_i == time):
			
# 			# print(time_i, '=', time, flagdeb)
# 			iexist = 1
# 			imissing = 0

# 			path[itime] = files_IR[itime_files][:ifile_IR.rfind('/')+1]

# 			if(flagdeb == 0):				  # on trouve l'image qui correspond au temps et il s'agit de la dernière image pour le run (on boucle dans le sens time = time-delta_reso)
# 				state[itime]   = 1
# 				flagdeb        = 1

# 				flagOK         = 1 
# 				exists[itime]  = iexist       # l'image existe
# 				missing[itime] = imissing     
# 				mode[itime]    = 0            # pour les mode de GEOS
# 				file[itime]    = ifile_name   # nom du fichier
				
# 				# print (file[itime])

# 			else :							# on trouve l'image qui correspond au temps et il ne s'agit pas de la derniere image pour le run
# 				state[itime]       = 0
# 				exists[itime]      = iexist
# 				missing[itime]     = imissing
# 				_missing           = 0
# 				flagOK             = 1
# 				mode[itime]        = 0
# 				file[itime]        = ifile_name
				
# 				# print (file[itime])

# 				continue


# 	if(flagOK==0):
# 		print(f"missing {time}")
# 		_missing = _missing +1
# 		missing[itime] = _missing
# 		iexist = 0
# 		exists[itime] = iexist
# 		# flagdeb = 0
		
# 	time = time - delta_resotempo



# Début du run : state=2 pour la plus ancienne image non manquante dans la liste
for i in np.arange(1, nombre_de_pas+1, 1):
	if exists[-i] == 1 :
		state[-i] = 2	
		break
	else:
		continue


# definir state=1 avant un lot de 5 images manquantes ou plus
for i in np.arange(1, nombre_de_pas, 1):
	if missing[i-1] >= 5 and missing[i]==0 :
		state[i] = 1


# Gérer les lots de 5 images manquantes ou plus pour déterminer la valeur de state après ce trou, début nv run donc state=2
flag_run_cut = 0
flag_reprise = nombre_de_pas-1
for i in np.arange(nombre_de_pas-1, -1, -1):
	if missing[i] >= 5 :
		flag_run_cut = 1
	if flag_run_cut == 0 :
		continue
	if missing[i] == 0 :
		state[i] = 2
		flag_reprise = i # recupérer i en flag pour n'écrire dans le fichier que les lignes à partir du moment ou ça reprend ///// !!!!!!!!!!!!!!!!!!!!!!! \\\\\\
		flag_run_cut = 0


# print(file)


for itime in np.arange(nombre_de_pas-2, -1, -1):
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



######################################################################################
## Ecriture du fichier														   ####

# f = open('../../../FileINPUT/file_list_19810101-20231231_v2.08_15min_v3.csv','w')
# start = 0

# headers=['year'.rjust(8),'month'.rjust(8),'day'.rjust(8),'hour'.rjust(8),'minute'.rjust(8),'julian'.rjust(8),
# 		 'slot'.rjust(8),'exists'.rjust(8),'missing'.rjust(8),'state'.rjust(8),'mode'.rjust(8),
# 		 'path_ir'.rjust(100),'file_ir'.rjust(100)]



# for header in headers: f.write(header)
# f.write("\n")

# for itime in np.arange(nombre_de_pas-1, -1, -1):

# 	#if itime > flag_reprise : continue

# 	print(f'{year[itime]:5d}{month[itime]:3d}{day[itime]:3d}{hour[itime]:3d}{minute[itime]:3d}{julianday[itime]:8d}{slot[itime]:4d}{exists[itime]:3d}{missing[itime]:3d}{state[itime]:3d}{mode[itime]:3d}{path[itime].rjust(80)}{file[itime].rjust(75)}')

# 	f.write(f'{year[itime]:8d}{month[itime]:8d}{day[itime]:8d}{hour[itime]:8d}{minute[itime]:8d}{julianday[itime]:8d}{slot[itime]:8d}{exists[itime]:8d}{missing[itime]:8d}{state[itime]:8d}{mode[itime]:8d}{path[itime].rjust(100)}{file[itime].rjust(100)}\n')


# f.close()


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

# Index temporel (optionnel mais fortement conseillé)
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

# Enregistrement en netCDF
ds = xr.Dataset.from_dataframe(df)

ds.to_netcdf(
    "../../../FileINPUT/total_file_list_19810101-20041231_v2.08_05min_v3.nc"
)
# 	if((start == 1) ):
		
# 		print (year[itime],month[itime],day[itime],hour[itime],minute[itime],slot[itime],exists[itime],missing[itime],mode[itime],path[itime],file[itime])		
# 		# f.write( "%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%100s%100s\n" %(year[itime],month[itime],day[itime],hour[itime],minute[itime],julianday[itime],slot[itime],exists[itime],missing[itime],state[itime],mode[itime],path[itime],file[itime]) )
# 		f.write(f'{year[itime]:8d}{month[itime]:8d}{day[itime]:8d}{hour[itime]:8d}{minute[itime]:8d}{julianday[itime]:8d}{slot[itime]:8d}{exists[itime]:8d}{missing[itime]:8d}{state[itime]:8d}{mode[itime]:8d}{path[itime].rjust(100)}{file[itime].rjust(100)}\n')
	
# 	if((file[itime] != 'N/A') and (start == 0)):
		
# 		print (year[itime],month[itime],day[itime],hour[itime],minute[itime],slot[itime],exists[itime],missing[itime],mode[itime],path[itime],file[itime])
# 		# f.write( "%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%100s%100s\n" %(year[itime],month[itime],day[itime],hour[itime],minute[itime],julianday[itime],slot[itime],exists[itime],missing[itime],state[itime],mode[itime],path[itime],file[itime]) )
# 		f.write(f'{year[itime]:8d}{month[itime]:8d}{day[itime]:8d}{hour[itime]:8d}{minute[itime]:8d}{julianday[itime]:8d}{slot[itime]:8d}{exists[itime]:8d}{missing[itime]:8d}{state[itime]:8d}{mode[itime]:8d}{path[itime].rjust(100)}{file[itime].rjust(100)}\n')
# 		start = 1
		
	# if(state[itime] == 1): 
	# 	break		


