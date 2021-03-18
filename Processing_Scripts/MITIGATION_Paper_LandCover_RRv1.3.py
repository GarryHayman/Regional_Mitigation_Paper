#!/bin/env python2.7

#import matplotlib as mpl
#mpl.use('Agg')
#mpl.style.use('classic')

import numpy as np
import netCDF4 as nc
import sys,os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from copy import deepcopy
import datetime as dt

from imogen import data_info
import plot_functions

#import ipdb

def box_and_whisker(DS,xpos,color,ax,boxwidth=0.1,med_color=None,fill=False):
    if med_color==None: med_color=color
    datastats = DS.describe()
    patch = patches.Rectangle( (xpos-(boxwidth/2.),datastats['25%']),
                                boxwidth, datastats['75%']-datastats['25%'],
                                fill=fill, edgecolor=color, color=color, lw=1.5)
    ax.add_patch(patch)
    line = ax.plot( [xpos-(boxwidth/2.),xpos+(boxwidth/2.)],
                    [datastats['50%'],datastats['50%']], c=med_color, lw=2)
    ax.plot( [xpos,xpos],[datastats['75%'],datastats['max']], c=color, lw=1.5)
    ax.plot( [xpos,xpos],[datastats['25%'],datastats['min']], c=color, lw=1.5)
    return line

def optional_argparse(arg,default):
    if arg in sys.argv:
        temp_loc=sys.argv.index(arg)
        temp_arg=sys.argv.pop(temp_loc)
        value=sys.argv.pop(temp_loc)
    else:
        value=default
    return value

INTERACTIVE   = '-interactive' in sys.argv
kg_to_Gt      = 1e-12
kg_to_Mt      = 1e-9
kg_to_Tg      = kg_to_Mt
kg_to_Gg      = 1e-6
m2_to_Mha     = 1e-10
GtC_to_ppm    = 0.471
C_to_CH4      = 16.04/12.011
ppm_to_kgC    = 1e12/GtC_to_ppm
C_to_water    = 530.0*1.0E-03 # Convert from GtC yr-1 to Tm3 yr-1
kg_to_Tm3     = 1.0E-15       # Assume 1 kg of water is 10-3 m3
npft          = 13

# Optional input parameters
DEBUG         = optional_argparse('-debug','N')
sIDX_REGION   = optional_argparse('-regions','ALL')
subregions    = optional_argparse('-subregions','IMAGE').upper()
PLATFORM      = optional_argparse('-platform','JASMIN')
sDATE         = optional_argparse('-date',dt.datetime.strftime(dt.datetime.now(),'%Y%m%d'))
PLOT_TAG      = optional_argparse('-plottag', 'LandUse')
PLOT_OPT      = str(optional_argparse('-plot_opt','1'))
LAND_npy      = optional_argparse('-LAND_NPY','Use')
FILE_EXT      = optional_argparse('-ext','.jpg')

# Directories containing JULES output and plot output directories:
if PLATFORM == 'JASMIN':
    HOME_DIR      = '/gws/nopw/j04/clifftop/SYNTHESIS/PostReview_output/'
    DATA_DIR      = HOME_DIR
    ANCILS_DIR    ='/gws/nopw/j04/clifftop/COMMON_DATA/ANCILS/'
    PLOT_DIR      = optional_argparse('-plotdir',DATA_DIR+'plots/'+PLOT_TAG+'/')
    os.system('mkdir -p '+PLOT_DIR )

    Q10_exps      = [ 'lowQ10', 'highQ10' ]
    Ozone_exps    = [['L','lowO3'], ['H','highO3'], ]
elif PLATFORM == 'CEH':
    HOME_DIR      = '/prj/CLIFFTOP/SYNTHESIS/'
    PLOT_DIR      = '/data/grp/eow/garr/Projects/CLIFFTOP/Paper_Synthesis/PLOTS_LAND_COVER/CEH/'
    DATA_DIR      = HOME_DIR+'Review_Response_Check/GCM_Output/'
    ANCILS_DIR    = HOME_DIR+'Land_Cover/'

    Q10_exps      = [ 'highQ10' ] 
    Ozone_exps    = [['L','lowO3']]

DATA_DIR=optional_argparse('-data_dir', DATA_DIR)
print('DATA_DIR: '+DATA_DIR)
print('PLOT_DIR: '+PLOT_DIR)

Tile_names    = data_info.TILE_short_names()
Tile_colours  = data_info.TILE_colours()
nTiles        = len(Tile_names)
nQ10s         = len(Q10_exps)
nO3s          = len(Ozone_exps)

COMPS_keys    = [ 'CTL', 'LULUC_CCS', 'LULUC_Nat' ]
COMPS         = { 
                  'CTL': { 'config': 'highCH4_OZONESUB_LULUCBL', 'runid':'H_OZONESUB_BL' }
                , 'LULUC_CCS': { 'config':'highCH4_OZONESUB_LULUC1.9', 'runid':'H_OZONESUB_19' }
                , 'LULUC_Nat': { 'config':'highCH4_OZONESUB_LULUC1.9Nat', 'runid':'H_OZONESUB_19N' }
                }

COMPS_DIFF    = { 
                  'LULUC_CCS': { 'config':'highCH4_OZONESUB_LULUC1.9', 'runid':'H_OZONESUB_19' }
                , 'LULUC_Nat': { 'config':'highCH4_OZONESUB_LULUC1.9Nat', 'runid':'H_OZONESUB_19N' }
                }

START_YEAR    = 2000
END_YEAR      = 2100
nYEARS        = END_YEAR-START_YEAR # No data for 2100

# Scenarios to plot:
TEMPs         = ['1p5deg', '2deg' ] #'1p81p5deg',   # tag in the JULES output file directory 
TEMP_names    = ['1.5$^o$C (2100)','2.0$^o$C (2100)'] # '1.5$^o$C Overshoot (2100)',
                # Name to appear on plots etc.
nTEMPs        = len(TEMPs)

# File containing the pre-industrial conditions to use as a baseline:

LAND_TYPES    = [ 'Trees','Agriculture','Bioenergy','Grasses','Shrubs' ]
nLAND_TYPES   = len(LAND_TYPES)

# Indices of pfts for aggregation to trees, grasses, shrubs
TREE_INDICES  = [ 0,1,2,3,4 ]
GRASS_INDICES = [ 5,8 ] # natural grasses
SHRUB_INDICES = [ 11,12 ]

# IMAGE regions

META_IMAGE    =  [ \
                  [ '1','Canada'          ,[[ 0.10, 0.02 ], [ 0.40, 0.05 ]], [[ 150,  25 ], [ 500,  50 ]], [[  80,  20 ], [ 100,  20 ]] ], \
                  [ '2','USA'             ,[[ 0.50, 0.10 ], [ 0.30, 0.05 ]], [[ 500,  50 ], [ 300,  50 ]], [[  50,  10 ], [  40,  10 ]] ], \
                  [ '3','Mexico'          ,[[ 0.70, 0.10 ], [ 0.40, 0.05 ]], [[ 120,  20 ], [  80,  10 ]], [[  50,  10 ], [  20,   5 ]] ], \
                  [ '4','C. America'      ,[[ 0.60, 0.10 ], [ 0.50, 0.10 ]], [[  50,  10 ], [  40,   5 ]], [[  25,   5 ], [  10,   2 ]] ], \
                  [ '5','Brazil'          ,[[ 0.60, 0.10 ], [ 0.60, 0.10 ]], [[ 500,  50 ], [ 500,  50 ]], [[ 250,  50 ], [ 100,  20 ]] ], \
                  [ '6','Rest S. America' ,[[ 0.50, 0.10 ], [ 0.40, 0.10 ]], [[ 400,  50 ], [ 300,  50 ]], [[ 150,  25 ], [  80,  20 ]] ], \
                  [ '7','N. Africa'       ,[[ 0.20, 0.05 ], [ 1.00, 0.10 ]], [[ 100,  20 ], [  80,  10 ]], [[  20,   5 ], [  30,   5 ]] ], \
                  [ '8','W. Africa'       ,[[ 0.70, 0.10 ], [ 0.40, 0.05 ]], [[ 700, 100 ], [ 400,  50 ]], [[ 300,  50 ], [ 200,  50 ]] ], \
                  [ '9','E. Africa'       ,[[ 0.70, 0.10 ], [ 0.50, 0.10 ]], [[ 500,  50 ], [ 400,  50 ]], [[ 200,  50 ], [  60,  10 ]] ], \
                  ['10','S. Africa'       ,[[ 1.00, 0.10 ], [ 1.00, 0.10 ]], [[ 120,  20 ], [ 100,  10 ]], [[  40,  10 ], [  10,   2 ]] ], \
                  ['11','W. Europe'       ,[[ 0.40, 0.05 ], [ 0.40, 0.05 ]], [[ 160,  20 ], [ 120,  20 ]], [[  25,   5 ], [  25,   5 ]] ], \
                  ['12','C. Europe'       ,[[ 0.70, 0.10 ], [ 0.50, 0.10 ]], [[  80,  10 ], [  60,  10 ]], [[  20,   5 ], [  25,   5 ]] ], \
                  ['13','Turkey'          ,[[ 0.80, 0.10 ], [ 0.40, 0.05 ]], [[  60,  10 ], [  40,   5 ]], [[  20,   5 ], [  10,   2 ]] ], \
                  ['14','Ukraine'         ,[[ 0.90, 0.10 ], [ 0.70, 0.10 ]], [[  60,  10 ], [  50,  10 ]], [[  12,   2 ], [  15,   5 ]] ], \
                  ['15','C. Asia'         ,[[ 0.70, 0.10 ], [ 0.50, 0.10 ]], [[ 250,  25 ], [ 220,  20 ]], [[ 100,  20 ], [  30,   5 ]] ], \
                  ['16','Russia'          ,[[ 0.25, 0.05 ], [ 0.30, 0.05 ]], [[ 500,  50 ], [ 600,  50 ]], [[ 150,  50 ], [ 150,  50 ]] ], \
                  ['17','Middle East'     ,[[ 0.40, 0.05 ], [ 1.00, 0.10 ]], [[ 200,  20 ], [ 140,  20 ]], [[  50,  10 ], [  10,   2 ]] ], \
                  ['18','India'           ,[[ 0.80, 0.10 ], [ 0.60, 0.10 ]], [[ 250,  25 ], [ 200,  20 ]], [[  50,  10 ], [  40,  10 ]] ], \
                  ['19','Korea'           ,[[ 0.50, 0.10 ], [ 0.70, 0.10 ]], [[  10,   1 ], [  15,   5 ]], [[   5,   1 ], [   5,   1 ]] ], \
                  ['20','China'           ,[[ 0.60, 0.10 ], [ 0.50, 0.10 ]], [[ 700, 100 ], [ 500,  50 ]], [[ 150,  50 ], [  60,  10 ]] ], \
                  ['21','S.E. Asia'       ,[[ 0.50, 0.10 ], [ 0.70, 0.10 ]], [[ 140,  20 ], [ 200,  20 ]], [[  50,  10 ], [  50,  10 ]] ], \
                  ['22','Indonesia'       ,[[ 0.50, 0.10 ], [ 0.70, 0.10 ]], [[ 100,  10 ], [ 120,  10 ]], [[  60,  10 ], [  10,   2 ]] ], \
                  ['23','Japan'           ,[[ 0.25, 0.05 ], [ 0.70, 0.10 ]], [[  12,   1 ], [  30,   5 ]], [[   5,   1 ], [   5,   1 ]] ], \
                  ['24','Oceania'         ,[[ 0.70, 0.10 ], [ 0.80, 0.10 ]], [[ 600, 100 ], [ 450,  50 ]], [[ 100,  20 ], [ 100,  20 ]] ], \
                  ['25','Rest S. Asia'    ,[[ 0.50, 0.10 ], [ 0.70, 0.10 ]], [[  70,  10 ], [  40,  10 ]], [[  15,   5 ], [   5,   1 ]] ], \
                  ['26','Rest S. Africa'  ,[[ 0.70, 0.10 ], [ 0.50, 0.10 ]], [[ 350,  50 ], [ 240,  40 ]], [[ 150,  25 ], [ 100,  20 ]] ]  \
                 ]
# Select subregions (IMAGE, TRANSCOM)
nREGIONs_ALL  = len(META_IMAGE)
REGION_dict   = data_info.REGION_DICTIONARIES()[subregions]

if sIDX_REGION == 'ALL':
        nREGIONs       = nREGIONs_ALL
        IDX_REGION     = [ str(iREGION+1) for iREGION in range(nREGIONs) ]
else:
        IDX_REGION     = sIDX_REGION.split(':')
        nREGIONs       = len(IDX_REGION)
print('Number of regions: ',nREGIONs)
print('REGION codes:      ',IDX_REGION)

# Directory of ancillary data:
# Grid File (My index for converting the 1D jules output to a 2D grid)
GRID_file= ANCILS_DIR+'grid_info.nc'
grinf=nc.Dataset(GRID_file,'r')
grindex=grinf.variables['land_index'][:]
lats_2d = grinf.variables['latitude'][:]
lons_2d = grinf.variables['longitude'][:]
#Area_2d = grinf.variables['Area'][:]       # I don't actually use this but it's here
land_index = grinf.variables['land_index'][:]
grinf.close()

# 1Dimension grid cell area data for calculating totals etc.
AREA_file=ANCILS_DIR+'Area_in_iris_format.nc'
Ainf=nc.Dataset(AREA_file,'r')
lats_1D = Ainf.variables['latitude'][:].squeeze()
lons_1D = Ainf.variables['longitude'][:].squeeze()
AREA_1D = Ainf.variables['area'][:].squeeze()
REGIONS_1D=Ainf.variables[REGION_dict['NCvarname']][:].squeeze()
Ainf.close()
#print(AREA_file)

####################################
# select GCMs:
GCMs=data_info.GCMs()
nGCMs=len(GCMs)
nGCMs_ALL=len(GCMs)
for igcm in range(nGCMs):
    print('%3i: '%igcm+GCMs[igcm])

if INTERACTIVE==True:
    GCM_index=raw_input('Select GCMs to plot (Press Enter to select all): ')
    if GCM_index.replace(' ','')!='':
        GCM_index = [ int(i) for i in GCM_index.split(':') ]
        GCMs=[ GCMs[i] for i in GCM_index ] 
        nGCMs=len(GCMs)
    else:
        GCM_index=[ i for i in range(nGCMs) ]
else:
    GCM_index=optional_argparse('-GCMs','ALL')
    if GCM_index=='ALL':
        GCM_index=[ i for i in range(nGCMs) ]
    else:
        GCM_index=[ int(i) for i in GCM_index.split(':') ]
        GCMs=[ GCMs[i] for i in GCM_index ] 
        nGCMs=len(GCMs)

cmip5_runs = [ data_info.cmip5_runs()[i] for i in GCM_index ]
print(' -GCMs ',GCM_index)
for igcm in range(nGCMs):
    print(igcm,GCM_index[igcm],GCMs[igcm])

# Select GCMs to plot:
###################################################################################################
# Read in the data

if LAND_npy == "Save":
    for comp in COMPS_keys:
        for itemp in range(nTEMPs): 
            temp=TEMPs[itemp]
            temptag=temp
            COMPS[comp][temp]={}
            for iQ10 in range(nQ10s):
                Q10 = Q10_exps[iQ10]
                COMPS[comp][temp][Q10]={}
                for iO3 in range(nO3s):
                    O3 = Ozone_exps[iO3]
                    COMPS[comp][temp][Q10][O3[1]]={ land_type:{ GCMs[igcm]: {} for igcm in range(nGCMs) } for land_type in LAND_TYPES  }
                    config=COMPS[comp]['config'].replace('OZONESUB',O3[1])
                    runid=COMPS[comp]['runid'].replace('OZONESUB',O3[0])
                    print(comp,temp,Q10,config,runid)
                    for igcm in range(nGCMs):
                        gcm=GCMs[igcm]
                        gcm_index=GCM_index[igcm]
                        print(gcm)
                        first_year = True
                        for iyear in range(nYEARS):
                            YEAR = START_YEAR+iyear

                            DUMP_FILE=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+temptag+'.dump.'+str(YEAR)+'0101.0.nc'
                            ANN_FILE=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+temptag+'.Annual_carbon.'+str(YEAR)+'.nc'
                            if iyear == 0 or iyear == nYEARS-1: print(YEAR,ANN_FILE,DUMP_FILE)
                            # Open current year files:
                            Dinf = nc.Dataset(DUMP_FILE,'r')
                            Ainf = nc.Dataset(ANN_FILE,'r')
                            
                            # Read in data from current year    
                            FRAC_PFTs  = Dinf.variables['frac'][:].squeeze()

                            if first_year:
                                first_year   = False
                                nLAND_PTS    = FRAC_PFTs.shape[1]
                                FRAC_TREES   = np.zeros((nYEARS,nLAND_PTS))
                                FRAC_GRASSES = np.zeros((nYEARS,nLAND_PTS))
                                FRAC_SHRUBS  = np.zeros((nYEARS,nLAND_PTS))
                                FRAC_AGRIC   = np.zeros((nYEARS,nLAND_PTS))
                                FRAC_BECCS   = np.zeros((nYEARS,nLAND_PTS))

                            # Trees:
                            for index in TREE_INDICES:  FRAC_TREES[iyear,:]   += FRAC_PFTs[index,:]  
                            # Natural grasses:
                            for index in GRASS_INDICES: FRAC_GRASSES[iyear,:] += FRAC_PFTs[index,:]  
                            # Shrubs:
                            for index in SHRUB_INDICES: FRAC_SHRUBS[iyear,:]  += FRAC_PFTs[index,:]

                            # All agriculture
                            # Issue with some of the annual files
                            TEMP_AGRIC          = Ainf.variables['frac_agr'][:].squeeze()
                            TEMP_PAST           = Ainf.variables['frac_past'][:].squeeze()
                            TEMP_BECCS          = Ainf.variables['frac_biocrop'][:].squeeze()

                            if len(TEMP_AGRIC) == 0 or len(TEMP_PAST) == 0:
                                print('No data for agriculture    for '+str(YEAR))
                                FRAC_AGRIC[iyear,:] = float('nan')
                            else:
                                FRAC_AGRIC[iyear,:] = ( TEMP_AGRIC + TEMP_PAST )

                            if len(TEMP_BECCS) == 0:
                                print('No data for bioenery crops for '+str(YEAR))
                                FRAC_BECCS[iyear,:] = float('nan')
                            else:
                                FRAC_BECCS[iyear,:] =   TEMP_BECCS

                            Dinf.close(); Ainf.close()

                        COMPS[comp][temp][Q10][O3[1]]['Trees'][gcm]        = FRAC_TREES  *AREA_1D*m2_to_Mha
                        COMPS[comp][temp][Q10][O3[1]]['Grasses'][gcm]      = FRAC_GRASSES*AREA_1D*m2_to_Mha
                        COMPS[comp][temp][Q10][O3[1]]['Shrubs'][gcm]       = FRAC_SHRUBS *AREA_1D*m2_to_Mha
                        COMPS[comp][temp][Q10][O3[1]]['Agriculture'][gcm]  = FRAC_AGRIC  *AREA_1D*m2_to_Mha
                        COMPS[comp][temp][Q10][O3[1]]['Bioenergy'][gcm]    = FRAC_BECCS  *AREA_1D*m2_to_Mha 

                        del FRAC_TREES, FRAC_GRASSES, FRAC_SHRUBS, FRAC_AGRIC, FRAC_BECCS

# Subtract baseline run

    for comp in ['LULUC_CCS', 'LULUC_Nat']:
        COMPS_DIFF[comp] = {}
        COMPS_DIFF[comp] = deepcopy(COMPS[comp])
        for itemp in range(nTEMPs): 
            temp=TEMPs[itemp]
            for iQ10 in range(nQ10s):
              Q10 = Q10_exps[iQ10]
              for iO3 in range(nO3s):
                O3 = Ozone_exps[iO3]
                for igcm in range(nGCMs):
                    gcm=GCMs[igcm]
                    #print(gcm)
                    for land_type in LAND_TYPES:
                        COMPS_DIFF[comp][temp][Q10][O3[1]][land_type][gcm] = \
                            COMPS[comp][temp][Q10][O3[1]][land_type][gcm] - \
                            COMPS['CTL'][temp][Q10][O3[1]][land_type][gcm]

################################################################################
# Get time series by IMAGE region
nREGIONs_DICT = REGION_dict['Nregions']

nSCENARIOs    = 2
print(nSCENARIOs,nTEMPs,nREGIONs_DICT,nLAND_TYPES,nYEARS,nGCMs,nQ10s,nO3s)

# Save numpy array DATA_REGION
if len(GCMs) == nGCMs_ALL:
        FILE_PART      = 'Plot_LandSurface_Area_Diff_Ctl_GCMs_All_'
else:
    FILE_PART      = 'Plot_LandSurface_Area_Diff_Ctl_GCMs_'
    for GCM in GCM_index: FILE_PART      = FILE_PART+('%02d_' % GCM)

FILE_NUMPY     = PLOT_DIR+FILE_PART+sDATE+'.npy'

if LAND_npy == "Save":

    DATA_REGION   = np.zeros((nSCENARIOs,nTEMPs,nREGIONs_DICT,nLAND_TYPES,nYEARS,nGCMs,nQ10s,nO3s))

    for itemp in range(nTEMPs):
        temp= TEMPs[itemp]
        for iQ10 in range(nQ10s):
            Q10 = Q10_exps[iQ10]
            for iO3 in range(nO3s):
                O3 = Ozone_exps[iO3]
                for igcm in range(nGCMs):
                    gcm=GCMs[igcm]
                    for itype in range(nLAND_TYPES):
                        land_type=LAND_TYPES[itype]
               
                        # Regional Breakdown of the area
                        for iregion in range(nREGIONs_DICT):
                            region =REGION_dict['Name'][iregion]
                            region_mask=REGIONS_1D==(iregion+1)
                            if region=='Global': region_mask[:]=True
                            if region=='International Transportation': region_mask[:]=False

                            for iyear in range(nYEARS):
                                DATA_REGION[0,itemp,iregion,itype,iyear,igcm,iQ10,iO3]   = \
                                     COMPS_DIFF['LULUC_CCS'][temp][Q10][O3[1]][land_type][gcm][iyear,region_mask].sum()

                                DATA_REGION[1,itemp,iregion,itype,iyear,igcm,iQ10,iO3] = \
                                     COMPS_DIFF['LULUC_Nat'][temp][Q10][O3[1]][land_type][gcm][iyear,region_mask].sum()

    np.save(FILE_NUMPY,DATA_REGION)

elif LAND_npy == 'Use':

    print('Reading from: '+FILE_NUMPY)
    LAND_ID      = open(FILE_NUMPY,'rb')
    DATA_REGION  = np.load(LAND_ID)
    LAND_ID.close()
 
###################################################################################################
# Plot of land cover

if True:
    Atmoscolor    = '#ffff99'
    Oceancolor    = '#a6cee3'
    CCScolor      = '#e6ab02'
    Landcolor     = '#1b9e77'
    totLULUCcolor ='#66a61e'
    CH4color      = '#7570b3'
    totcolor      = '#d95f02'
    totcolor_lin  = '#a6761d'
    CTLcolor      = '#666666' 
    bnwcolor      = '#e7298a'

WIDTH0,HEIGHT0 =    4.0, 5.0
nROWS          =    1
nCOLUMNS       =  nLAND_TYPES
WIDTH,HEIGHT   = WIDTH0*nCOLUMNS,HEIGHT0*nROWS

LINE_WIDTH     =    2.0
DPI            =  300
LEGEND_POS     =   -1
FONTSIZES      = [ 12,12,14,14 ]
PLOT_CODES_ALL = [ 'orange' ,'green' ,'black', 'blue' ,'magenta' ,'yellow' ,'cyan' ]
LINE_CODES_ALL = [ '--',    ':', '',  '',  '',  '',  ''  ]

xMIN,xMAX,xINC = START_YEAR,END_YEAR,25
nxTICKS        = int((xMAX-xMIN)/xINC)+1
xTICKS         = xMIN+xINC*np.arange(nxTICKS)
xPLOT_BASE     = START_YEAR+np.arange(nYEARS)
xLABEL         = 'Year'
yLABEL         = 'Area (Mha)'
iINDEX         = 1
LEGEND         = ['"BECCS" Scenario','"Natural" Scenario'] 

for iREGION in range(nREGIONs):
   
    jREGION        = int(IDX_REGION[iREGION]) 
    NAME_IMAGE     = META_IMAGE[jREGION-1][1]
    print('Land Cover - area:      '+NAME_IMAGE,REGION_dict['Name'][jREGION-1],iREGION,jREGION)
    FILE_IMAGE     = NAME_IMAGE.replace('.','').replace(' ','_')
    yMIN,yMAX,yINC = -int(META_IMAGE[jREGION-1][4][iINDEX][0]),int(META_IMAGE[jREGION-1][4][iINDEX][0]), \
                      int(META_IMAGE[jREGION-1][4][iINDEX][1])
    nyTICKS        = int((yMAX-yMIN)/yINC)+1
    yTICKS         = yMIN+yINC*np.arange(nyTICKS)

    xPLOT_A,yPLOT_A,yPLOT_AL,yPLOT_AH = [],[],[],[]
    xMIN_A,xMAX_A,xTICKS_A,xLABEL_A   = [],[],[],[]
    yMIN_A,yMAX_A,yTICKS_A,yLABEL_A   = [],[],[],[]

    nDATASETS      = []
    SUBTITLES      = []
    LINE_CODES     = []
    PLOT_CODES     = []
    PLOT_TEXT      = []
    PLOT_TITLE     = 'Land Cover by IMAGE Region'
    FILE_PLOT      = PLOT_DIR+FILE_PART+('R%02d_' % (jREGION))+FILE_IMAGE+'_'+sDATE+FILE_EXT

    for itype in range(nLAND_TYPES):

        xPLOT,yPLOT,yPLOT_L,yPLOT_H = [],[],[],[]

        xMIN_A.append(xMIN)
        xMAX_A.append(xMAX)
        xTICKS_A.append(xTICKS)
        xLABEL_A.append(xLABEL)

        yMIN_A.append(yMIN)
        yMAX_A.append(yMAX)
        yTICKS_A.append(yTICKS)
        yLABEL_A.append(yLABEL)

        SUBTITLES.append(LAND_TYPES[itype].replace('_',' '))
        PLOT_CODES.append(PLOT_CODES_ALL[0:2])
        LINE_CODES.append(LINE_CODES_ALL[0:2])
        nDATASETS.append(2)
    
        for iSCENARIO in range(nSCENARIOs):
            xPLOT_DATA         = deepcopy(xPLOT_BASE)
            yPLOT_DATA         = np.zeros(nYEARS)
            yPLOT_DATA_L       = np.zeros(nYEARS)
            yPLOT_DATA_H       = np.zeros(nYEARS)
            for iYEAR in range(nYEARS):
                 yPLOT_TEMP         = DATA_REGION[iSCENARIO,0,jREGION-1,itype,iYEAR,:,:,:]
                 yPLOT_TEMP         = yPLOT_TEMP[~np.isnan(yPLOT_TEMP)].flatten()
                 if len(yPLOT_TEMP) == 0:
                     yPLOT_DATA[iYEAR]  = float('nan')
                     yPLOT_DATA_L[iYEAR]= float('nan')
                     yPLOT_DATA_H[iYEAR]= float('nan')
                 else:
                     yPLOT_DATA[iYEAR]  = np.median(yPLOT_TEMP[~np.isnan(yPLOT_TEMP)])
                     yPLOT_DATA_L[iYEAR]= np.percentile(yPLOT_TEMP[~np.isnan(yPLOT_TEMP)],25.0)
                     yPLOT_DATA_H[iYEAR]= np.percentile(yPLOT_TEMP[~np.isnan(yPLOT_TEMP)],75.0)

            xPLOT.append(xPLOT_DATA)
            yPLOT.append(yPLOT_DATA)
            yPLOT_L.append(yPLOT_DATA_L)
            yPLOT_H.append(yPLOT_DATA_H)
            
#           print(iREGION,itype,iSCENARIO,yPLOT_DATA,yPLOT_DATA_L,yPLOT_DATA_H)
            del xPLOT_DATA,yPLOT_DATA,yPLOT_DATA_L,yPLOT_DATA_H

        xPLOT_A.append(np.swapaxes(np.array(xPLOT),0,1))
        yPLOT_A.append(np.swapaxes(np.array(yPLOT),0,1))
        yPLOT_AL.append(np.swapaxes(np.array(yPLOT_L),0,1))
        yPLOT_AH.append(np.swapaxes(np.array(yPLOT_H),0,1))

    plot_functions.Plot_General_MultiPlot_general2_range( \
        nROWS,nCOLUMNS,nDATASETS,xPLOT_A,yPLOT_A,yPLOT_AL,yPLOT_AH, \
        xMIN_A,xMAX_A,xTICKS_A,xTICKS_A,xLABEL_A, \
        yMIN_A,yMAX_A,yTICKS_A,yTICKS_A,yLABEL_A, \
        WIDTH,HEIGHT,FONTSIZES,PLOT_TITLE,SUBTITLES,PLOT_TEXT, \
        LEGEND,LEGEND_POS,PLOT_CODES,PLOT_OPT,FILE_PLOT,DEBUG, \
        DPI=DPI,LW=LINE_WIDTH,LINE_CODE=LINE_CODES)
#       DPI=DPI,LINE_CODE=LINE_CODES)
