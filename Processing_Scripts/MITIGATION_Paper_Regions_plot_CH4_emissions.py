#!/bin/python

# Python module to generate regional plot for Regional Mitigation paper

# All on one plot
# ./MITIGATION_Paper_Regions_plot.py N All 0:6:11:12:16 4 7 2

# Split into several plots
# ./MITIGATION_Paper_Regions_plot_CH4_emissions.py N RCP19 26:1:23:25:2:22:17:10:11:21:20:3 0:6:11:12:16 4 3 2
# ./MITIGATION_Paper_Regions_plot_CH4_emissions.py N RCP19 8:7:13:24:4:9:5:19:6:0:15:18 0:6:11:12:16 4 3 2
# ./MITIGATION_Paper_Regions_plot_CH4_emissions.py N RCP19 8:7:13:24:4:9:5:19:6:14:16:12 0:6:11:12:16 4 3 2

# ESD Paper, April 2020
# ./MITIGATION_Paper_Regions_plot_CH4_emissions.py N RCP19 26:1:23:9:17:3:0:9:5:15:18:12 0:6:11:12:16 4 3 2

# Latex version, January 2021
# ./MITIGATION_Paper_Regions_plot_CH4_emissions.py N RCP19 26:1:23:9:17:3:0:9:5:15:18:12 0:6:11:12:16 4 3 2 pdf 20210102

# Garry Hayman
# Centre for Ecology and Hydrology
# June 2018

# Contains

import os
import sys
import numpy as np
import datetime

import plot_functions


DEBUG          = sys.argv[1]
RCP_OPT        = sys.argv[2]
REGIONS        = sys.argv[3]
SECTORS        = sys.argv[4]
NROWS          = int(sys.argv[5])
NCOLUMNS       = int(sys.argv[6])
PLOT_OPT       = sys.argv[7]
FILE_EXT       = sys.argv[8]
sDATE          = sys.argv[9]

HOME_DIR       = '/users/eow/garr/Work/Projects/CLIFFTOP/Paper_Synthesis/PLOTS_Methane/'
MISS_DATA      = -999.9
LEGEND_POS     =  -1
if FILE_EXT not in ['png','jpg','eps','pdf']:
	print "Invalid file format selected - will default to jpg"
	FILE_EXT       = 'jpg'

# Original settings
#WIDTH0,HEIGHT0 =    4.0, 5.0
#XMIN,XMAX,XINC = 2000,2100,10
#FONTSIZES      = [ 8,8,8,10 ]
# Settings for Figures in Synthesis Paper - submitted (tight_layout commented out)
#WIDTH0,HEIGHT0 =    4.0, 5.0
#XMIN,XMAX,XINC = 2000,2100,20
#FONTSIZES      = [ 10,10,10,12 ]
#LINE_WIDTH     =   1.0
#DPI            = 600
# Settings for Figures in Synthesis Paper - revised: same size with increased font size & linewidth 
WIDTH0,HEIGHT0 =    4.0, 5.0
XMIN,XMAX,XINC = 2000,2100,20
FONTSIZES      = [ 12,12,14,14 ]
LINE_WIDTH     =   2.0
DPI            = 300
# Settings for Figures in Synthesis Paper - revised: test smaller figure setting
#WIDTH0,HEIGHT0 =    2.0, 2.5
#XMIN,XMAX,XINC = 2000,2100,20
#FONTSIZES      = [  6, 6, 6, 8 ]
#LINE_WIDTH     =   1.0


os.chdir(HOME_DIR)

# Information on JULES offline runs

if RCP_OPT == 'RCP26':
	FILES_CH4      = [ \
#			   [' 0:','RCP2.6:', 'SSP2_CH4emissions_Region_RCP26.csv' ], \
#			   [' 1:','Ref:',    'SSP2_CH4emissions_Region_Ref.csv'   ]  \
			   [' 0:','SSP2 RCP2.6:',   'SSP2_CH4emissions_Region_RCP26.csv' ], \
			   [' 1:','SSP2 Baseline:', 'SSP2_CH4emissions_Region_Ref.csv'   ]  \
			 ]
elif RCP_OPT == 'RCP19':
	FILES_CH4      = [ \
			   [' 0:','SSP2 RCP1.9:',   'SSP2_CH4emissions_Region_RCP19.csv' ], \
			   [' 1:','SSP2 Baseline:', 'SSP2_CH4emissions_Region_Ref.csv'   ]  \
			 ]

DATA_REGIONS   = dict([ \
		   (  '1', [ 'CAN',  'Canada',              10,   1 ]), \
		   ( '23', [ 'USA',  'USA',                 60,  10 ]), \
		   ( '10', [ 'MEX',  'Mexico',              10,   1 ]), \
		   ( '13', [ 'RCAM', 'Central America',     10,   1 ]), \
		   (  '0', [ 'BRA',  'Brazil',              40,   5 ]), \
		   ( '15', [ 'RSAM', 'Rest of S. America',  50,   5 ]), \
		   ( '11', [ 'NAF',  'Northern Africa',     10,   1 ]), \
		   ( '24', [ 'WAF',  'Western Africa',      50,   5 ]), \
		   (  '4', [ 'EAF',  'Eastern Africa',      50,   5 ]), \
		   ( '18', [ 'SAF',  'South Africa',        50,   5 ]), \
		   ( '14', [ 'RSAF', 'Rest of S. Africa',   10,   1 ]), \
		   ( '25', [ 'WEU',  'Western Europe',      50,   5 ]), \
		   (  '2', [ 'CEU',  'Central Europe',      10,   1 ]), \
		   ( '21', [ 'TUR',  'Turkey',              10,   1 ]), \
		   ( '22', [ 'UKR',  'Ukraine Region',      10,   1 ]), \
		   ( '20', [ 'STAN', 'Central Asia',        20,   5 ]), \
		   ( '17', [ 'RUS',  'Russia Region',       50,   5 ]), \
		   (  '9', [ 'ME',   'Middle East',         50,   5 ]), \
		   (  '5', [ 'INDIA','India',              100,  10 ]), \
		   (  '8', [ 'KOR',  'Korea Region',        10,   1 ]), \
		   (  '3', [ 'CHN',  'China',              100,  10 ]), \
		   ( '19', [ 'SEAS', 'Southeastern Asia',   25,   5 ]), \
		   (  '6', [ 'INDO', 'Indonesia',           20,   5 ]), \
		   (  '7', [ 'JAP',  'Japan',               10,   1 ]), \
		   ( '16', [ 'RSAS', 'Rest of S. Asia',     20,   5 ]), \
		   ( '12', [ 'OCE',  'Oceania',             30,   5 ]), \
		   ( '26', [ 'World','World',              500,  50 ])  \
		])

DATA_SECTORS  = [ \
		  [ ' 0', 'Total'       ,'Emissions|CH4' ], \
		  [ ' 1', 'Energy Ind'  ,'Emissions|CH4|Energy Demand|Industry'], \
		  [ ' 2', 'Energy Res'  ,'Emissions|CH4|Energy Demand|Residential and Commercial'], \
		  [ ' 3', 'Energy Tra'  ,'Emissions|CH4|Energy Demand|Transportation'], \
		  [ ' 4', 'Energy Tra'  ,'Emissions|CH4|Energy Demand|Transportation|Ground Transportation'], \
		  [ ' 5', 'Energy Sup'  ,'Emissions|CH4|Energy Supply'], \
		  [ ' 6', 'Energy'      ,'Emissions|CH4|Energy Supply and Demand'], \
		  [ ' 7', 'Land Use'    ,'Emissions|CH4|Land Use'], \
		  [ ' 8', 'Agric Burn'  ,'Emissions|CH4|Land Use|Agricultural Waste Burning'], \
		  [ ' 9', 'Agric'       ,'Emissions|CH4|Land Use|Agriculture'], \
		  [ '10', 'Agric Waste' ,'Emissions|CH4|Land Use|Agriculture|AWM'], \
		  [ '11', 'Agric Cattle','Emissions|CH4|Land Use|Agriculture|Enteric Fermentation'], \
		  [ '12', 'Agric Rice'  ,'Emissions|CH4|Land Use|Agriculture|Rice'], \
		  [ '13', 'Forest Burn' ,'Emissions|CH4|Land Use|Forest Burning'], \
		  [ '14', 'Savanna Burn','Emissions|CH4|Land Use|Savannah Burning'], \
		  [ '15', 'Other'       ,'Emissions|CH4|Other'], \
		  [ '16', 'Waste'       ,'Emissions|CH4|Waste']  \
		]

nTIMES         = 11
nFILES_CH4     = len(FILES_CH4)
nREGIONS_ALL   = len(DATA_REGIONS)
nSECTORS_ALL   = len(DATA_SECTORS)

if REGIONS == 'All':
	nREGIONS       = nREGIONS_ALL
#	REGIONS_NUM    = [ str(iREGION) for iREGION in range(nREGIONS) ]
	REGIONS_NUM    = [ '26','1','23','25','2','22','17','10','11','21','20','3','8', \
			   '7','13','24','4','9','5','19','6','0','15','18','14','16','12' ]
else:
	REGIONS_NUM    = REGIONS.split(':')
        nREGIONS       = len(REGIONS_NUM)
print('Region codes: ',REGIONS_NUM)

if SECTORS == 'All':
	nSECTORS       = nSECTORS_ALL
        SECTORS_NUM    = [ str(iSECTOR) for iSECTOR in range(nSECTORS) ]
else:
	SECTORS_NUM    = SECTORS.split(':')
        nSECTORS       = len(SECTORS_NUM)
print('Sector codes: ',SECTORS_NUM)

if nREGIONS > NROWS*NCOLUMNS:
	print('Invalid plot configuration: nREGIONS <= nROWS*nCOLUMNS')
	quit()

# Loop over and input from CH4 files

INPUT_CH4_ALL  = []
SCENARIO_CH4   = []

for iCH4 in range(nFILES_CH4):

	INPUT_CH4_FILE = []
	SCENARIO_CH4.append(FILES_CH4[iCH4][1])
	FILE_CSV       = FILES_CH4[iCH4][2]
	print('Opening file:  '+FILE_CSV)

	for LINE in open(FILE_CSV):
		LINE_SPLIT     = LINE.replace('\n','').replace('\r','').split(',')
		INPUT_CH4_FILE.append(LINE_SPLIT)
		if DEBUG == 'Y': print(len(LINE_SPLIT),LINE_SPLIT)

	INPUT_CH4_ALL.append(INPUT_CH4_FILE)

# Extract information - each file: nREGIONS_ALL * nSECTORS_ALL + 1
# Assume regions and sectors in order in REGIONS_DATA and SECTORS_DATA

WIDTH,HEIGHT   = WIDTH0*NCOLUMNS,HEIGHT0*NROWS

NXTICKS        = int((XMAX-XMIN)/XINC)+1
XTICKS         = XMIN+XINC*np.arange(NXTICKS)
XLABEL         = 'Year'

YLABEL         = 'CH$_4$ Emissions (Mt CH$_4$ yr$^{-1}$)'

XPLOT_A,YPLOT_A                 = [],[]
XMIN_A,XMAX_A,XTICKS_A,XLABEL_A = [],[],[],[]
YMIN_A,YMAX_A,YTICKS_A,YLABEL_A = [],[],[],[]
NDATASETS      = []
SUBTITLES      = []
LEGEND         = []
PLOT_CODES     = []
LINE_CODES     = []

PLOT_TEXT      = []
PLOT_TITLE     = 'CH$_4$ Emissions by IMAGE Region'
FILE_PLOT      = HOME_DIR+'Plot_Methane_Emissions_Image_Region_'+sDATE+'.'+FILE_EXT
LINE_CODES_BAS = [':','-']
PLOT_COLOURS   = ['black','darkorange','blue' ,'green' ,'magenta' ,'yellow' ,'cyan' ]

for sREGION in REGIONS_NUM:

	iREGION        = int(sREGION)
	xFIRST         = True
	XPLOT,YPLOT    = [],[]

	XMIN_A.append(XMIN)
	XMAX_A.append(XMAX)
	XTICKS_A.append(XTICKS)
	XLABEL_A.append(XLABEL)

	YMIN,YMAX,YINC = 0,int(DATA_REGIONS[sREGION][2]),int(DATA_REGIONS[sREGION][3])
	NYTICKS        = int((YMAX-YMIN)/YINC)+1
	YTICKS         = YMIN+YINC*np.arange(NYTICKS)

	YMIN_A.append(YMIN)
	YMAX_A.append(YMAX)
	YTICKS_A.append(YTICKS)
	YLABEL_A.append(YLABEL)
	SUBTITLES.append(DATA_REGIONS[sREGION][1])
	
	LINE_CODES_TMP = []
	PLOT_CODES_TMP = []
        iCODE          = 0

	for sSECTOR in SECTORS_NUM:
		iSECTOR        = int(sSECTOR)
		iLINE          = iSECTOR+iREGION*nSECTORS_ALL+1

		TEMP_X         = []
		TEMP_Y         = []

		for iTIME in range(nTIMES):
			TEMP_X.append(float(INPUT_CH4_ALL[0][0][iTIME+5]))
			TEMP_Y.append(float(INPUT_CH4_ALL[0][iLINE][iTIME+5]))
		if xFIRST:
			XPLOT.append(TEMP_X)
			xFIRST         = False
		YPLOT.append(TEMP_Y)	

		LINE_CODES_TMP.append(LINE_CODES_BAS[0])
		PLOT_CODES_TMP.append(PLOT_COLOURS[iCODE])
		iCODE += 1

        iCODE          = 0

	for sSECTOR in SECTORS_NUM:
		iSECTOR        = int(sSECTOR)
		iLINE          = iSECTOR+iREGION*nSECTORS_ALL+1

		TEMP_X         = []
		TEMP_Y         = []

		for iTIME in range(nTIMES):
			TEMP_X.append(float(INPUT_CH4_ALL[1][0][iTIME+5]))
			TEMP_Y.append(float(INPUT_CH4_ALL[1][iLINE][iTIME+5]))
		YPLOT.append(TEMP_Y)	

		LINE_CODES_TMP.append(LINE_CODES_BAS[1])
		PLOT_CODES_TMP.append(PLOT_COLOURS[iCODE])
		iCODE += 1

	LINE_CODES.append(LINE_CODES_TMP)
	PLOT_CODES.append(PLOT_CODES_TMP)

	XPLOT_A.append(np.array(XPLOT).squeeze())
	YPLOT_A.append(np.swapaxes(np.array(YPLOT),0,1))
	NDATASETS.append(nSECTORS*2)

# Names of sectors
for sSECTOR in SECTORS_NUM:
	iSECTOR        = int(sSECTOR)
	LEGEND.append(SCENARIO_CH4[0]+' '+DATA_SECTORS[iSECTOR][1])
for sSECTOR in SECTORS_NUM:
	iSECTOR        = int(sSECTOR)
	LEGEND.append(SCENARIO_CH4[1]+' '+DATA_SECTORS[iSECTOR][1])

#DEBUG = 'Y'
if DEBUG == 'Y':
	print('XPLOT_A:    ',XPLOT_A)
	print('YPLOT_A:    ',YPLOT_A)
	print('PLOT_CODES: ',PLOT_CODES)
	print('LEGEND:     ',LEGEND_POS,LEGEND)

plot_functions.Plot_General_MultiPlot_general2( \
	NROWS,NCOLUMNS,NDATASETS,XPLOT_A,YPLOT_A, \
	XMIN_A,XMAX_A,XTICKS_A,XTICKS_A,XLABEL_A, \
	YMIN_A,YMAX_A,YTICKS_A,YTICKS_A,YLABEL_A, \
	WIDTH,HEIGHT,FONTSIZES,PLOT_TITLE,SUBTITLES,PLOT_TEXT, \
	LEGEND,LEGEND_POS,PLOT_CODES,PLOT_OPT,FILE_PLOT,DEBUG,DPI,LINE_WIDTH,LINE_CODES)

# End of Program
