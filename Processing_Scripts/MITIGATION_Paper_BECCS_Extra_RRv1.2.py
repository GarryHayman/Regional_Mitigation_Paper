#!/bin/env python2.7

import matplotlib as mpl
#mpl.use('Agg')
#mpl.style.use('classic')

import numpy as np
import netCDF4 as nc
import sys,os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gs
from copy import deepcopy
import datetime as dt

from imogen import data_info
from PlotTools import plot_tools as PTs

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

INTERACTIVE  = '-interactive' in sys.argv
kg_to_t      = 1e-3
kg_to_Gt     = 1e-12
kg_to_Mt     = 1e-9
kg_to_Tg     = kg_to_Mt
kg_to_Gg     = 1e-6
m2_to_Mha    = 1e-10
m2_to_Ha     = 1e-4
GtC_to_ppm   = 0.471
C_to_CH4     = 16.04/12.011
ppm_to_kgC   = 1e12/GtC_to_ppm
C_to_water   = 530.0*1.0E-03            # Convert from GtC yr-1 to Tm3 yr-1
kg_to_Tm3    = 1.0E-15                  # Assume 1 kg of water is 10-3 m3
sec_to_year  = 3600.0*24.0*360.0        # 360 day year
                                        # * kgC to biomass * 90% dry matter
BECCS_harvest_conv        = (1./m2_to_Ha) * kg_to_t * 2.0  * 0.9 * (365./360.)
kappa_efficiency_ratio    = 0.6/0.87
 
# Optional input parameters
DEBUG         = optional_argparse('-debug','N')
subregions    = optional_argparse('-subregions','IMAGE').upper()
PLATFORM      = optional_argparse('-platform','JASMIN')
sDATE         = optional_argparse('-date',dt.datetime.strftime(dt.datetime.now(),'%Y%m%d'))
PLOT_TAG      = optional_argparse('-plottag', 'BECCS_Extra')
VERSION       = optional_argparse('-BECCS','max_bioenergy')
BECCS_NPP     = optional_argparse('-BECCS_NPP','N')
PROC_STEP     = optional_argparse('-BECCS_NPY','Save')
CCS_minimum_threshold     = float(optional_argparse('CCS_min','1e-4'))
BIOFRAC_minimum_threshold = float(optional_argparse('CCSfrac_min','1e-2'))
FILE_EXT      = optional_argparse('-ext','.jpg')
PLOT_FIGURE   = False

# Directories containing JULES output and plot output directories:
if PLATFORM == 'JASMIN':
    HOME_DIR      = '/gws/nopw/j04/clifftop/SYNTHESIS/PostReview_output/'
    DATA_DIR      = HOME_DIR
    ANCILS_DIR    = '/gws/nopw/j04/clifftop/COMMON_DATA/ANCILS/'
    PLOT_DIR      = optional_argparse('-plotdir',DATA_DIR+'plots/'+PLOT_TAG+'/')
    BECCS_npy_DIR = '/gws/nopw/j04/jules/aharper/PYTHON/SYNTHESIS/npy_files/'+VERSION+'/'
    BECCS_npy_NPP = '/gws/nopw/j04/jules/ghayman/PYTHON/SYNTHESIS/npy_files/max_bioenergy_NPP/'
    SCNPP_npy_DIR = '/gws/nopw/j04/jules/ghayman/PYTHON/SYNTHESIS/npy_files/Carbon/'
    COMPS_npy_DIR = '/gws/nopw/j04/jules/ghayman/PYTHON/SYNTHESIS/npy_files/processed_output/'
    LAND_npy_DIR  = '/gws/nopw/j04/jules/ghayman/PYTHON/SYNTHESIS/npy_files/land_cover/'
#   BECCS_npy_NPP = '/work/scratch/ghayman/SYNTHESIS/npy_files/max_bioenergy_NPP/'
#   SCNPP_npy_DIR = '/work/scratch/ghayman/SYNTHESIS/npy_files/Carbon/'
#   COMPS_npy_DIR = '/work/scratch/ghayman/SYNTHESIS/npy_files/processed_output/'
#   LAND_npy_DIR  = '/work/scratch/ghayman/SYNTHESIS/npy_files/land_cover/'

    Q10_exps      = [ 'lowQ10', 'highQ10' ]
    Ozone_exps    = [['L','lowO3'], ['H','highO3'], ]

    COMPS_keys    = [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]

elif PLATFORM == 'CEH':
    HOME_DIR      = '/prj/CLIFFTOP/SYNTHESIS/'
    PLOT_DIR      = HOME_DIR+'Review_Response_Check2/plots/'
    DATA_DIR      = HOME_DIR+'Review_Response_Check2/GCM_Output/'
    ANCILS_DIR    = HOME_DIR+'Land_Cover/'
    BECCS_npy_DIR = HOME_DIR+'Review_Response_Check2/'+VERSION+'/'
    BECCS_npy_NPP = HOME_DIR+'Review_Response_Check2/max_bioenergy_NPP/'
    SCNPP_npy_DIR = HOME_DIR+'Review_Response_Check2/Carbon/'
    COMPS_npy_DIR = HOME_DIR+'Review_Response_Check2/processed_output/'
    LAND_npy_DIR  = HOME_DIR+'Review_Response_Check2/land_cover/'

    if PROC_STEP == "Save":
        Q10_exps      = [ 'highQ10' ]
        Ozone_exps    = [['L','lowO3']]
    else:
        Q10_exps      = [ 'lowQ10', 'highQ10' ]
        Ozone_exps    = [['L','lowO3'], ['H','highO3'], ]

    COMPS_keys    = [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat' ]

COMPS_opt     = [ 'LULUC_opt','Coupled_opt' ]
COMPS_keys_all= [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]+[ 'LULUC_opt','Coupled_opt' ]
COMPS         = {
                  'CTL': { 'config': 'highCH4_OZONESUB_LULUCBL', 'runid':'H_OZONESUB_BL' }
                , 'CH4': { 'config':'lowCH4_OZONESUB_LULUCBL', 'runid':'L_OZONESUB_BL' }
                , 'LULUC_CCS': { 'config':'highCH4_OZONESUB_LULUC1.9', 'runid':'H_OZONESUB_19' }
                , 'LULUC_Nat': { 'config':'highCH4_OZONESUB_LULUC1.9Nat', 'runid':'H_OZONESUB_19N' }
                , 'Coupled_CCS': { 'config':'lowCH4_OZONESUB_LULUC1.9', 'runid':'L_OZONESUB_19' }
                , 'Coupled_Nat': { 'config':'lowCH4_OZONESUB_LULUC1.9Nat', 'runid':'L_OZONESUB_19N' }
                }

COMPS_DIFF    = {
                  'LULUC_CCS': { 'config':'highCH4_OZONESUB_LULUC1.9', 'runid':'H_OZONESUB_19' }
                , 'LULUC_Nat': { 'config':'highCH4_OZONESUB_LULUC1.9Nat', 'runid':'H_OZONESUB_19N' }
                }

# Directories containing JULES output and plot output directories:
DATA_DIR      = optional_argparse('-data_dir', DATA_DIR)
PLOT_DIR      = optional_argparse('-plotdir',  PLOT_DIR+PLOT_TAG+'/')
print('DATA_DIR: '+DATA_DIR)
print('PLOT_DIR: '+PLOT_DIR)
os.system('mkdir -p '+PLOT_DIR )
os.system('mkdir -p '+PLOT_DIR+'npy_files' )

PRESENT_DAY_YEAR = 2015

ALPHABET      = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o']

START_YEAR    = PRESENT_DAY_YEAR
END_YEAR      = 2100
nYEARS        = END_YEAR-START_YEAR # No data for 2100

# Scenarios to plot:
TEMPs         = ['1p5deg', '2deg' ] #'1p81p5deg',   # tag in the JULES output file directory
TEMP_years    = [2099,2099,2099]    # tag in the JULES output file directory 
TEMP_names    = ['1.5$^o$C (2100)','2.0$^o$C (2100)'] # '1.5$^o$C Overshoot (2100)',
                # Name to appear on plots etc.

GLOBAL_pools  = [ 'Total','Atmos','Ocean' ]
LAND_pools    = [ 'Land','CS','CV','WP','CCS','BE_Harvest' ]
AREA_pools    = [ 'BECCS_Area' ]
BECCS_pools   = [ 'BECCS_productivity', 'BECCS_ScaleFactor', \
                  'BECCS_productivity_JULES','BECCS_productivity_NPP' ]
              #+[ 'BECCS_productivity','Harvest','kappaHarvest','kappamefficHarvest']
SOIL_pools    = [ 'Soil_Carbon', 'Soil_Carbon_gb', 'NPP_gb' ]

dzsoil        = np.array([ 0.05,0.08408964,0.11397535,0.14142136,0.16718508,0.19168293, \
                  0.21517585,0.23784142,0.25980762,0.28117066,0.30200527, \
                  0.32237098,0.34231625,0.36188121 ])

nLAYERs       = 7
soil_depth    = dzsoil[0:nLAYERs].sum()
 
LAND_TYPES    = [ 'Trees','Agriculture','Bioenergy','Grasses','Shrubs' ]
nLAND_TYPES   = len(LAND_TYPES)
nSCENARIOs    = 2
pools         = GLOBAL_pools+LAND_pools+BECCS_pools+SOIL_pools+LAND_TYPES+AREA_pools

Tile_names    = data_info.TILE_short_names()
Tile_colours  = data_info.TILE_colours()
nTiles        = len(Tile_names)
nTEMPs        = len(TEMPs)
nQ10s         = len(Q10_exps)
nO3s          = len(Ozone_exps)
npools        = len(pools)
nPFTs         = 13
nSOIL_L       = 14
nSOIL_P       =  4
nLAND_pts     = 1631

# Indices of pfts for aggregation to trees, grasses, shrubs
TREE_INDICES  = [ 0,1,2,3,4 ]
GRASS_INDICES = [ 5,8 ] # natural grasses
CROP_INDICES  = [ 6,9 ]
PAST_INDICES  = [ 7,10 ]
SHRUB_INDICES = [ 11,12 ]

#  BECCS_multiplier not required for BECCS sensitivty routine
### BECCS_multiplier = float(optional_argparse('-beccs_multiplier','1'))
### print(BECCS_multiplier)

# Select subregions (IMAGE, TRANSCOM)
REGION_dict   = data_info.REGION_DICTIONARIES()[subregions]
nREGIONs_DICT = REGION_dict['Nregions']
REGION_idx    = [ idx for idx in range(11) ] + [ 25 ] + [ idx for idx in range(11,23) ] + [ 24, 23, 27 ]
nREGIONS      = len(REGION_idx)

# Directory of Ocean Uptake data:
OCEAN_UPTAKE_DIR = optional_argparse('-ocean_uptake_dir',DATA_DIR) 
OCEAN_START_YEAR = 1850
print("Ocean Uptake data from: " + OCEAN_UPTAKE_DIR)

# Directory of ancillary data:
# Grid File (My index for converting the 1D jules output to a 2D grid)
GRID_file        = ANCILS_DIR+'grid_info.nc'
grinf            = nc.Dataset(GRID_file,'r')
grindex          = grinf.variables['land_index'][:]
lats_2d          = grinf.variables['latitude'][:]
lons_2d          = grinf.variables['longitude'][:]
#Area_2d         = grinf.variables['Area'][:]       # I don't actually use this but it's here
land_index       = grinf.variables['land_index'][:]
grinf.close()

# 1Dimension grid cell area data for calculating totals etc.
AREA_file        = ANCILS_DIR+'Area_in_iris_format.nc'
Ainf             = nc.Dataset(AREA_file,'r')
lats_1D          = Ainf.variables['latitude'][:].squeeze()
lons_1D          = Ainf.variables['longitude'][:].squeeze()
AREA_1D          = Ainf.variables['area'][:].squeeze()
REGIONS_1D       = Ainf.variables[REGION_dict['NCvarname']][:].squeeze()
MAXBIOFRAC_1D    = Ainf.variables['maxbiofrac_19'][:].squeeze()
Ainf.close()
#print(AREA_file)

MAXBIOFRAC_2D    = np.ma.masked_array(MAXBIOFRAC_1D[grindex], mask=grindex.mask)

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
    GCM_index  = optional_argparse('-GCMs','ALL')
    sGCM_index = GCM_index
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

# Arrays and partial filenames for land-use files
LAND_REGION    = np.zeros((nSCENARIOs,nTEMPs,nREGIONs_DICT,nLAND_TYPES,nYEARS,nGCMs,nQ10s,nO3s))
LAND_LANDPTS   = np.zeros((nSCENARIOs,nTEMPs,nLAND_TYPES,nYEARS,nGCMs,nQ10s,nO3s,nLAND_pts))
nFACTORIAL     = nGCMs*nQ10s*nO3s

FILE_PART      = 'Land_Area_Max_'
FILE_PART2     = 'Land_Area_Grid_'

if sGCM_index == 'ALL':
    FILE_PART3     = 'GCM_All_'+sDATE
else:
    FILE_PART3     = 'GCM_'+sGCM_index.replace(':','_')+'_'+sDATE

###################################################################################################

# Read in the Control data
OCEAN_START_YEAR=1850

if PROC_STEP == 'Save':
    for comp in COMPS_keys:
        comp_Ocean_in = { Q10 : { O3[1]:  np.load(OCEAN_UPTAKE_DIR+Q10+'/' 
                                        + COMPS[comp]['config'].replace('OZONESUB',O3[1])+'/'
                                        + COMPS[comp]['runid'].replace('OZONESUB',O3[0])
                                        +'_ocean_uptake_accum.npy' ) * -1. / kg_to_Gt
                         for O3 in Ozone_exps} for Q10 in Q10_exps }
        
        comp_Ocean = { '1p5deg':   { Q10:  { O3[1]: comp_Ocean_in[Q10][O3[1]][:,2,:] for O3 in Ozone_exps } for Q10 in Q10_exps},
                       '1p81p5deg':{ Q10:  { O3[1]: comp_Ocean_in[Q10][O3[1]][:,1,:] for O3 in Ozone_exps } for Q10 in Q10_exps},
                       '2deg':     { Q10:  { O3[1]: comp_Ocean_in[Q10][O3[1]][:,0,:] for O3 in Ozone_exps } for Q10 in Q10_exps}
                       }
        CH4 = COMPS[comp]['config'].split('_')[0]
        for iTEMP in range(nTEMPs): 
            TEMP=TEMPs[iTEMP]
            TEMP_year=TEMP_years[iTEMP]
            TEMP_tag=TEMP
            COMPS[comp][TEMP]={}
            for iQ10 in range(nQ10s):
                Q10 = Q10_exps[iQ10]
                COMPS[comp][TEMP][Q10]={}
                for iO3 in range(nO3s):
                    O3 = Ozone_exps[iO3]
                    COMPS[comp][TEMP][Q10][O3[1]]={ pool:{} for pool in pools  }  #+land_opt_pools }
                    config=COMPS[comp]['config'].replace('OZONESUB',O3[1])
                    runid=COMPS[comp]['runid'].replace('OZONESUB',O3[0])
                    print('GCM Input:',comp,TEMP,TEMP_year,Q10,config,runid)
                    for igcm in range(nGCMs):
                        gcm=GCMs[igcm]
                        gcm_index=GCM_index[igcm]
                        print(gcm)
                        DUMP_FILE=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+TEMP_tag+'.dump.YYYY0101.0.nc'
                        Ann_File=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+TEMP_tag+'.Annual_carbon.YYYY.nc'
                        #BECCSFile=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+TEMP_tag+'.be_harvest_gb.MAX.nc'
                        #DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+TEMP_tag+'.be_harvest_gb.MAX.nc'
                        BECCS_numpy_File=BECCS_npy_DIR + CH4+'_'+O3[1]+'_'+Q10+'_'+gcm+'_'+TEMP_tag+'.npy'
                        print(BECCS_numpy_File)
                        #quit()
                        #import  ipdb
                        #ipdb.set_trace()
                        #print(Ann_File)
                        # Open end of sim files:
                        Dinf = nc.Dataset(DUMP_FILE.replace('YYYY',str(TEMP_year+1)),'r')
                        Ainf = nc.Dataset(Ann_File.replace('YYYY',str(TEMP_year)),'r')
                        # Open Present day files:
                        Dinf_0 = nc.Dataset(DUMP_FILE.replace('YYYY',str(START_YEAR+1)),'r')
                        Ainf_0 = nc.Dataset(Ann_File.replace('YYYY',str(START_YEAR)),'r')
                        # and Max BECCS file:
                        #print(BECCSFile)
                        #Binf = nc.Dataset(BECCSFile,'r')
                        # Read in data, find delta from present day,    
                        CV = ( Ainf.variables['cv'][:]-Ainf_0.variables['cv'][:] ).squeeze() * AREA_1D
                        COMPS[comp][TEMP][Q10][O3[1]]['CV'][gcm] = CV *kg_to_Gt
                        
                        CS = ( Ainf.variables['cs_gb'][:]-Ainf_0.variables['cs_gb'][:] ).squeeze() * AREA_1D
                        COMPS[comp][TEMP][Q10][O3[1]]['CS'][gcm] = CS*kg_to_Gt
                        
                        # Correction for bug in test runs, can remove this with final runs, although will not affect results
                        #if '_Nat' in comp:
                        #    CCS = np.zeros_like(CS)
                        #else:
                        CCS = ( Ainf.variables['ccs_gb'][:]-Ainf_0.variables['ccs_gb'][:] ).squeeze() * AREA_1D 
                                    #  * BECCS_multiplier (Not Required for BECCS sensitivity script)
                        CCS[np.isfinite(CCS)==False]=0.
                        COMPS[comp][TEMP][Q10][O3[1]]['CCS'][gcm] = CCS*kg_to_Gt
                        
                        # Read in the maximum BECCS harvest rate, to scale with kappa_threshold for Figure2
                        #BE_H = ( Binf.variables['be_harvest_gb'][:] ).squeeze()
                        #BE_H[ np.isinf(BE_H) ] = 0.
                        BE_H = np.load(BECCS_numpy_File)  
                        ##COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity'][gcm] = BE_H * m2_to_Ha * kg_to_t 
                        COMPS[comp][TEMP][Q10][O3[1]]['BE_Harvest'][gcm] = BE_H * BECCS_harvest_conv 
                        COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_JULES'][gcm] = ( 
                                COMPS[comp][TEMP][Q10][O3[1]]['BE_Harvest'][gcm] )
                        COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity'][gcm] = ( 
                                COMPS[comp][TEMP][Q10][O3[1]]['BE_Harvest'][gcm] * kappa_efficiency_ratio )

                        data_temp = COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_JULES'][gcm]
                        data_temp = data_temp[MAXBIOFRAC_1D>BIOFRAC_minimum_threshold]
                        print('BE mean: ',comp,TEMP,Q10,O3[1],gcm,data_temp.shape,data_temp.mean())

                        #if comp=='LULUC_CCS':
                        #    ipdb.set_trace() 
                        WP =(( Dinf.variables['wood_prod_fast'][:]+Dinf.variables['wood_prod_med'][:] +
                               Dinf.variables['wood_prod_slow'][:]  )
                           - ( Dinf_0.variables['wood_prod_fast'][:]+Dinf_0.variables['wood_prod_med'][:] +
                               Dinf_0.variables['wood_prod_slow'][:]  ) )  *  AREA_1D
                        WP[np.isfinite(WP)==False]=0.
                        COMPS[comp][TEMP][Q10][O3[1]]['WP'][gcm] = WP*kg_to_Gt

                        AtmCO2_ppm = Dinf.variables['co2_ppmv'][0]-Dinf_0.variables['co2_ppmv'][0]
                        AtmCO2_kg = AtmCO2_ppm*ppm_to_kgC
                        COMPS[comp][TEMP][Q10][O3[1]]['Atmos'][gcm] = AtmCO2_kg *kg_to_Gt
                        
                        Ocean =  comp_Ocean[TEMP][Q10][O3[1]][gcm_index,TEMP_year-OCEAN_START_YEAR]   \
                               - comp_Ocean[TEMP][Q10][O3[1]][gcm_index,START_YEAR-OCEAN_START_YEAR]
                        COMPS[comp][TEMP][Q10][O3[1]]['Ocean'][gcm] = Ocean *kg_to_Gt
                        
                        COMPS[comp][TEMP][Q10][O3[1]]['Land'][gcm] = (  
                                                                       COMPS[comp][TEMP][Q10][O3[1]]['CV'][gcm]
                                                                     + COMPS[comp][TEMP][Q10][O3[1]]['CS'][gcm]
                                                                     + COMPS[comp][TEMP][Q10][O3[1]]['WP'][gcm] ) 
                        COMPS[comp][TEMP][Q10][O3[1]]['Total'][gcm] = (  
                                                                        COMPS[comp][TEMP][Q10][O3[1]]['Land'][gcm].sum()
                                                                      + COMPS[comp][TEMP][Q10][O3[1]]['CCS'][gcm].sum()
                                                                      + COMPS[comp][TEMP][Q10][O3[1]]['Atmos'][gcm]
                                                                      + COMPS[comp][TEMP][Q10][O3[1]]['Ocean'][gcm] )
                        
                        Dinf.close(); Ainf.close(); Dinf_0.close(); Ainf_0.close(); # Binf.close()

                        # Productivity based on NPP
                        first_year = True
                        if ('CTL' in comp or 'CCS' in comp or 'Nat' in comp) and (BECCS_NPP == 'Y'):

                            # Skip if files already exist for bioenergy, NPP and soil carbon
                            BECCS_NPP_FILE = BECCS_npy_NPP+comp+'_'+CH4+'_'+O3[1]+'_'+Q10+'_'+gcm+'_'+TEMP_tag+'_NPP.npy'
                            SC_NPP_FILE    = SCNPP_npy_DIR+comp+'_'+CH4+'_'+O3[1]+'_'+Q10+'_'+gcm+'_'+TEMP_tag+'_SCNPP.pkl'

                            print(not os.path.exists(BECCS_NPP_FILE))
                            print(not os.path.exists(SC_NPP_FILE))

                            if not os.path.exists(BECCS_NPP_FILE) and not os.path.exists(SC_NPP_FILE):

                                for iyear in range(nYEARS):
                                     YEAR = START_YEAR+iyear
                                     DUMP_FILE=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+TEMP_tag+'.dump.'+str(YEAR)+'0101.0.nc'
                                     ANN_FILE=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+TEMP_tag+'.Annual_carbon.'+str(YEAR)+'.nc'
                                     if iyear == 0 or iyear == nYEARS-1: print(YEAR,ANN_FILE,DUMP_FILE)
                                     # Open current year files:
                                     Dinf = nc.Dataset(DUMP_FILE,'r')
                                     Ainf = nc.Dataset(ANN_FILE,'r')

                                     # Read in data from current year
                                     FRAC       = Dinf.variables['frac'][:].squeeze()
                                     NPP        = Ainf.variables['npp'][:].squeeze()
                                     NPP_GB     = Ainf.variables['npp_gb'][:].squeeze()
                                     SOILC      = Ainf.variables['cs'][:].squeeze()
                                     SOILC_GB   = Ainf.variables['cs_gb'][:].squeeze()

                                     if first_year:
                                         first_year   = False
                                         nLAND_PTS    = FRAC.shape[1]
                                         NPP_ALL      = np.zeros((nYEARS,nPFTs,nLAND_PTS))
                                         NPP_GB_ALL   = np.zeros((nYEARS,nLAND_PTS))
                                         SOILC_ALL    = np.zeros((nYEARS,nSOIL_P,nSOIL_L,nLAND_PTS))
                                         SOILC_GB_ALL = np.zeros((nYEARS,nLAND_PTS))
                                         FRAC_ALL     = np.zeros((nYEARS,nPFTs,nLAND_PTS))
                                         FRAC_TREES   = np.zeros((nYEARS,nLAND_PTS))
                                         FRAC_GRASSES = np.zeros((nYEARS,nLAND_PTS))
                                         FRAC_SHRUBS  = np.zeros((nYEARS,nLAND_PTS))
                                         FRAC_CROPS   = np.zeros((nYEARS,nLAND_PTS))
                                         FRAC_PAST    = np.zeros((nYEARS,nLAND_PTS))
                                         FRAC_AGRIC   = np.zeros((nYEARS,nLAND_PTS))
                                         FRAC_BECCS   = np.zeros((nYEARS,nLAND_PTS))
                                         NPP_BECCS    = np.zeros((nYEARS,nLAND_PTS))

                                     if len(FRAC) == 0:
                                         print('No data for fractions      for '+str(YEAR))
                                         FRAC_ALL[iyear,:,:]     = float('nan')
                                     else:
                                         FRAC_ALL[iyear,:,:]     = FRAC[0:nPFTs,:]

                                     if len(NPP) == 0:
                                         print('No data for npp            for '+str(YEAR))
                                         NPP_ALL[iyear,:,:]      = float('nan')
                                     else:
                                         NPP_ALL[iyear,:,:]      = NPP

                                     if len(NPP_GB) == 0:
                                         print('No data for npp_gb         for '+str(YEAR))
                                         NPP_GB_ALL[iyear,:]     = float('nan')
                                     else:
                                         NPP_GB_ALL[iyear,:]     = NPP_GB

                                     if len(SOILC) == 0:
                                         print('No data for soil carbon    for '+str(YEAR))
                                         SOILC_ALL[iyear,:,:,:]  = float('nan')
                                     else:
                                         SOILC_ALL[iyear,:,:,:]  = SOILC

                                     if len(SOILC_GB) == 0:
                                         print('No data for soil carbon gb for '+str(YEAR))
                                         SOILC_GB_ALL[iyear,:]   = float('nan')
                                     else:
                                         SOILC_GB_ALL[iyear,:]   = SOILC_GB

                                     # Trees:
                                     for index in TREE_INDICES:  FRAC_TREES[iyear,:]   += FRAC_ALL[iyear,index,:]
                                     # Natural grasses:
                                     for index in GRASS_INDICES: FRAC_GRASSES[iyear,:] += FRAC_ALL[iyear,index,:]
                                     # Shrubs:
                                     for index in SHRUB_INDICES: FRAC_SHRUBS[iyear,:]  += FRAC_ALL[iyear,index,:]
                                     # Agriculture - Crops
                                     for index in CROP_INDICES:  FRAC_CROPS[iyear,:]   += FRAC_ALL[iyear,index,:]
                                     # Agriculture - Pasture
                                     for index in PAST_INDICES:  FRAC_PAST[iyear,:]    += FRAC_ALL[iyear,index,:]

                                     # All agriculture
                                     # Issue with some of the annual files
                                     TEMP_AGRIC          = Ainf.variables['frac_agr'][:].squeeze()
                                     TEMP_BECCS          = Ainf.variables['frac_biocrop'][:].squeeze()

                                     if len(TEMP_AGRIC) == 0:
                                         print('No data for agriculture    for '+str(YEAR))
                                         FRAC_AGRIC[iyear,:] = float('nan')
                                     else:
                                         FRAC_AGRIC[iyear,:] = TEMP_AGRIC

                                     if len(TEMP_BECCS) == 0:
                                         print('No data for bioenery crops for '+str(YEAR))
                                         FRAC_BECCS[iyear,:] = float('nan')
                                     else:
                                         FRAC_BECCS[iyear,:] = TEMP_BECCS

                                     Dinf.close(); Ainf.close()

                                     for index in CROP_INDICES:
                                          NPP_BECCS[iyear,:]  += FRAC_ALL[iyear,index,:]*NPP_ALL[iyear,index,:] \
                                                                *0.5*BECCS_harvest_conv*sec_to_year \
                                                                /FRAC_CROPS[iyear,:]
                                                                 
                                     if iyear == 0 or iyear == nYEARS-1:
                                         NPP_TEMP            = NPP_BECCS[iyear,:]
                                         print(YEAR,FRAC_CROPS[iyear,:].min(),FRAC_CROPS[iyear,:].max(), \
                                                    FRAC_PAST[iyear,:].min(),FRAC_PAST[iyear,:].max(),   \
                                                    FRAC_AGRIC[iyear,:].min(),FRAC_AGRIC[iyear,:].max(), \
                                                    FRAC_BECCS[iyear,:].min(),FRAC_BECCS[iyear,:].max() )

                                COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_NPP'][gcm] = NPP_BECCS[0,:]
                                COMPS[comp][TEMP][Q10][O3[1]]['Trees'][gcm]        = FRAC_TREES  *AREA_1D*m2_to_Mha
                                COMPS[comp][TEMP][Q10][O3[1]]['Grasses'][gcm]      = FRAC_GRASSES*AREA_1D*m2_to_Mha
                                COMPS[comp][TEMP][Q10][O3[1]]['Shrubs'][gcm]       = FRAC_SHRUBS *AREA_1D*m2_to_Mha
                                COMPS[comp][TEMP][Q10][O3[1]]['Agriculture'][gcm]  = FRAC_AGRIC  *AREA_1D*m2_to_Mha
                                COMPS[comp][TEMP][Q10][O3[1]]['Bioenergy'][gcm]    = FRAC_BECCS  *AREA_1D*m2_to_Mha

                                for iLAND in range(nLAND_PTS):
                                     NPP_TEMP            = NPP_BECCS[:,iLAND]
                                     if len(NPP_TEMP[~np.isnan(NPP_TEMP)]) == 0:
                                         COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_NPP'][gcm][iLAND] = float('nan') 
                                     else: 
                                         COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_NPP'][gcm][iLAND] = NPP_TEMP[~np.isnan(NPP_TEMP)].max() 

                                if ('CCS' in comp):
                                    # Save as numpy file: filename defined earlier (l. 357)
                                    np.save(BECCS_NPP_FILE,COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_NPP'][gcm])
                                    print('Writing to: '+BECCS_NPP_FILE)

                                if ('CCS' in comp or 'Nat' in comp):
                                    # Save as pickle file: filename defined earlier (l. 358)
                                    SC_NPP_FILE_ID = open(SC_NPP_FILE,'wb')
                                    print('Writing to: '+SC_NPP_FILE)
                                    pickle.dump( [ NPP_GB_ALL, SOILC_ALL, SOILC_GB_ALL ], SC_NPP_FILE_ID )
                                    SC_NPP_FILE_ID.close()

                                    COMPS[comp][TEMP][Q10][O3[1]]['NPP_gb'][gcm]         = NPP_GB_ALL
                                    COMPS[comp][TEMP][Q10][O3[1]]['Soil_Carbon'][gcm]    = SOILC_ALL
                                    COMPS[comp][TEMP][Q10][O3[1]]['Soil_Carbon_gb'][gcm] = SOILC_GB_ALL

                                del FRAC,FRAC_ALL,FRAC_TREES,FRAC_GRASSES,FRAC_SHRUBS, \
                                    FRAC_CROPS,FRAC_PAST,FRAC_AGRIC,FRAC_BECCS,NPP,NPP_ALL,NPP_TEMP, \
                                    NPP_GB, SOILC_ALL, SOILC_GB_ALL

                            else:

                                if ('CCS' in comp):
                                    # Read from numpy file
                                    COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_NPP'][gcm] = \
                                        deepcopy(COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity'][gcm])
                                    print('Reading from: '+BECCS_NPP_FILE)
                                    COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_NPP'][gcm] = \
                                        np.load(BECCS_NPP_FILE)

                                # Read from pickle file
                                SC_NPP_FILE_ID = open(SC_NPP_FILE,'rb')
                                print('Reading from: '+SC_NPP_FILE)
                                PKL_IN         = pickle.load(SC_NPP_FILE_ID)
                                COMPS[comp][TEMP][Q10][O3[1]]['NPP_gb'][gcm]         = PKL_IN[0]
                                COMPS[comp][TEMP][Q10][O3[1]]['Soil_Carbon'][gcm]    = PKL_IN[1]
                                COMPS[comp][TEMP][Q10][O3[1]]['Soil_Carbon_gb'][gcm] = PKL_IN[2]

                                SC_NPP_FILE_ID.close()
                                del PKL_IN

                        if not gcm in COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_NPP'].keys():
                            # Set to zero
                            COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_NPP'][gcm] = \
                                deepcopy(COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity'][gcm])
                            COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_NPP'][gcm][:] = 0.0

                        data_temp = COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity_NPP'][gcm]
                        data_temp = data_temp[MAXBIOFRAC_1D>BIOFRAC_minimum_threshold]
                        print('NPP mean: ',comp,TEMP,Q10,O3[1],gcm,data_temp.shape,data_temp.mean())
                        print

                        if ('CCS' in comp or 'Nat' in comp):
                            # Modify NPP and Soil Carbon:
                            # Change in variable over time series
                            
                            # (1) NPP_gb and convert to tC per hectare per year
                            var       = 'NPP_gb'
                            data_temp = COMPS[comp][TEMP][Q10][O3[1]][var][gcm]
                            print(var,data_temp.shape)
                            COMPS[comp][TEMP][Q10][O3[1]][var][gcm] = (data_temp[-1,:]-data_temp[0,:]) * \
                                sec_to_year
                            print(var+': ',COMPS[comp][TEMP][Q10][O3[1]][var][gcm][:].min(), \
                                           COMPS[comp][TEMP][Q10][O3[1]][var][gcm][:].max() )

                            # (2) Grid-box mean soil carbon and convert to tC per hectare
                            var       = 'Soil_Carbon_gb'
                            data_temp = COMPS[comp][TEMP][Q10][O3[1]][var][gcm]
                            COMPS[comp][TEMP][Q10][O3[1]][var][gcm] = (data_temp[-1,:]-data_temp[0,:])
                            print(var+': ',COMPS[comp][TEMP][Q10][O3[1]][var][gcm][:].min(), \
                                           COMPS[comp][TEMP][Q10][O3[1]][var][gcm][:].max() )

                            # (3) Soil carbon and convert to tC per hectare
                            var       = 'Soil_Carbon'
                            data_temp = np.zeros((COMPS[comp][TEMP][Q10][O3[1]][var][gcm].shape[0], \
                                                  COMPS[comp][TEMP][Q10][O3[1]][var][gcm].shape[3]))
                            for iLAYER in range(nLAYERs):
                                for iPOOL in range(2):
                                    data_temp[:,:] = data_temp[:,:] + \
                                        COMPS[comp][TEMP][Q10][O3[1]][var][gcm][:,iPOOL,iLAYER,:] * \
                                        dzsoil[iLAYER]/soil_depth
                            COMPS[comp][TEMP][Q10][O3[1]][var][gcm] = (data_temp[-1,:]-data_temp[0,:])
                            print(var+': ',COMPS[comp][TEMP][Q10][O3[1]][var][gcm][:].min(), \
                                           COMPS[comp][TEMP][Q10][O3[1]][var][gcm][:].max() )


    # Need to create dummy datasets for CH4, 
    # Full set: COMPS_keys    = [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]
    # CEH  set: COMPS_keys    = [ 'CTL', 'LULUC_CCS', 'LULUC_Nat' ]
    if PLATFORM == 'CEH':
        COMPS['CH4']         = deepcopy(COMPS['CTL'])
        COMPS['Coupled_CCS'] = deepcopy(COMPS['LULUC_CCS'])
        COMPS['Coupled_Nat'] = deepcopy(COMPS['LULUC_Nat'])

    #ipdb.set_trace()
    # Create optimised LULUC mitigation option by choosing CCS or return to Natural Vegetation
    for comp in COMPS_opt:
        compCCS = comp.replace('opt','CCS')
        compNat = comp.replace('opt','Nat')
        # copy _CCS dictionay to _opt
        COMPS[comp] = deepcopy(COMPS[compCCS])
        for iTEMP in range(nTEMPs): 
            TEMP=TEMPs[iTEMP]
            TEMP_tag=TEMP
            for iQ10 in range(nQ10s):
                Q10 = Q10_exps[iQ10]
                for iO3 in range(nO3s):
                    O3 = Ozone_exps[iO3]
                    for igcm in range(nGCMs):
                        gcm=GCMs[igcm]
                        print('Optimisation: ',comp,TEMP,Q10,O3,gcm)
                        # create mask, CCS>Nat = 1; CCS==Nat = 0; CCS<Nat = -1 
                        difference = ( ( COMPS[compCCS][TEMP][Q10][O3[1]]['Land'][gcm]
                                       + COMPS[compCCS][TEMP][Q10][O3[1]]['CCS'][gcm] )
                                     - ( COMPS[compNat][TEMP][Q10][O3[1]]['Land'][gcm]
                                       + COMPS[compNat][TEMP][Q10][O3[1]]['CCS'][gcm] )  )
                        flag_mask = (difference/np.abs(difference)).astype(int)
                        flag_mask[difference==0.] = 0
                        
                        # Substitute Nat for CCS in land/CCS arrays where Nat is the prefered choice:
                        for pool in LAND_pools:
                            COMPS[comp][TEMP][Q10][O3[1]][pool][gcm][flag_mask==-1] = \
                                    COMPS[compNat][TEMP][Q10][O3[1]][pool][gcm][flag_mask==-1]
                        
                        # Recalculate the Total emission budget:
                        COMPS[comp][TEMP][Q10][O3[1]]['Total'][gcm] = (
                                                 COMPS[comp][TEMP][Q10][O3[1]]['Land'][gcm].sum()
                                               + COMPS[comp][TEMP][Q10][O3[1]]['CCS'][gcm].sum()
                                               + COMPS[comp][TEMP][Q10][O3[1]]['Atmos'][gcm]
                                               + COMPS[comp][TEMP][Q10][O3[1]]['Ocean'][gcm] )

                        # Calculate the required scale factor for BECCS to become viable:
                        #   What scale factor is required for the CCS to be greater than the 
                        #   benefit of returning land to natural vegetaiton?
                        req_scale_factor = ( COMPS[compNat][TEMP][Q10][O3[1]]['Land'][gcm]
                                           - COMPS[compCCS][TEMP][Q10][O3[1]]['Land'][gcm])
                        #ccs_mask = COMPS[compCCS][TEMP][Q10][O3[1]]['CCS'][gcm]>CCS_minimum_threshold
                        ccs_mask = (MAXBIOFRAC_1D>BIOFRAC_minimum_threshold) \
                                & (COMPS[compCCS][TEMP][Q10][O3[1]]['CCS'][gcm]>CCS_minimum_threshold)
                        req_scale_factor[ccs_mask] /= COMPS[compCCS][TEMP][Q10][O3[1]]['CCS'][gcm][ccs_mask]
                        req_scale_factor[~ccs_mask] = -1e20
                        COMPS[comp][TEMP][Q10][O3[1]]['BECCS_ScaleFactor'][gcm] = \
                                np.ma.masked_equal(req_scale_factor,-1e20)

                        # Calculate the required productivity for BECCS to become viable:
                        #   including an efficiency increase of the farm to final store (kappa_efficiency_ratio)
                        req_BECCS_product = deepcopy( COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity'][gcm])
                        work_scale_factor = deepcopy(req_scale_factor)
                        #ipdb.set_trace()
                        work_scale_factor[ccs_mask][work_scale_factor[ccs_mask]<1.0] = 1.0
                        req_BECCS_product[ccs_mask] = req_BECCS_product[ccs_mask] * work_scale_factor[ccs_mask]
                        COMPS[comp][TEMP][Q10][O3[1]]['BECCS_productivity'][gcm] = \
                                np.ma.masked_equal(req_BECCS_product,-1e20)

                        # Process NPP_gb and Soil Carbon
                        if comp == 'LULUC_opt':
                            pool = 'NPP_gb'
                            COMPS[comp][TEMP][Q10][O3[1]][pool][gcm] = \
                                (COMPS[compCCS][TEMP][Q10][O3[1]][pool][gcm] - \
                                 COMPS[compNat][TEMP][Q10][O3[1]][pool][gcm])*ccs_mask  

                            pool = 'Soil_Carbon_gb'
                            COMPS[comp][TEMP][Q10][O3[1]][pool][gcm] = \
                                (COMPS[compCCS][TEMP][Q10][O3[1]][pool][gcm] - \
                                 COMPS[compNat][TEMP][Q10][O3[1]][pool][gcm])*ccs_mask  

                            pool = 'Soil_Carbon'
                            COMPS[comp][TEMP][Q10][O3[1]][pool][gcm] = \
                                (COMPS[compCCS][TEMP][Q10][O3[1]][pool][gcm] - \
                                 COMPS[compNat][TEMP][Q10][O3[1]][pool][gcm])*ccs_mask  

    for comp in COMPS_keys_all:
        CH4 = COMPS[comp]['config'].split('_')[0]

        pools_pkl = deepcopy(pools)

        if comp in [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]:
            del pools_pkl[pools_pkl.index('BECCS_ScaleFactor')]
        if comp in [ 'CTL', 'CH4' ]:
            for pool in SOIL_pools: del pools_pkl[pools_pkl.index(pool)]

        for pool in LAND_TYPES+AREA_pools: del pools_pkl[pools_pkl.index(pool)]

        for iTEMP in range(nTEMPs): 
            TEMP=TEMPs[iTEMP]
            TEMP_tag=TEMP
            for iQ10 in range(nQ10s):
                Q10 = Q10_exps[iQ10]
                for iO3 in range(nO3s):
                    O3 = Ozone_exps[iO3]
                    for igcm in range(nGCMs):
                        gcm=GCMs[igcm]

                        # Save COMPS[comp][TEMP][Q10][O3[1]][pool][gcm] as pickle file
                        COMPS_FILE     = COMPS_npy_DIR+comp+'_'+CH4+'_'+O3[1]+'_'+Q10+'_'+gcm+'_'+TEMP_tag+'_processed.pkl'
                        COMPS_FILE_ID  = open(COMPS_FILE,'wb')
                        print('Writing to: '+COMPS_FILE)
                        pickle.dump( [ COMPS[comp][TEMP][Q10][O3[1]][pool][gcm] for pool in pools_pkl ], COMPS_FILE_ID )
                        COMPS_FILE_ID.close()

    # Land use: subtract baseline run

    for comp in ['LULUC_CCS', 'LULUC_Nat']:
        COMPS_DIFF[comp] = deepcopy(COMPS['CTL'])
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
                            print(land_type,COMPS_DIFF[comp][temp][Q10][O3[1]][land_type], \
                                            COMPS[comp][temp][Q10][O3[1]][land_type], \
                                            COMPS['CTL'][temp][Q10][O3[1]][land_type])
                            COMPS_DIFF[comp][temp][Q10][O3[1]][land_type][gcm] = \
                                COMPS[comp][temp][Q10][O3[1]][land_type][gcm] - \
                                COMPS['CTL'][temp][Q10][O3[1]][land_type][gcm]


    # Save land use data into numpy file

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
                                LAND_REGION[0,itemp,iregion,itype,iyear,igcm,iQ10,iO3]   = \
                                     COMPS_DIFF['LULUC_CCS'][temp][Q10][O3[1]][land_type][gcm][iyear,region_mask].sum()

                                LAND_REGION[1,itemp,iregion,itype,iyear,igcm,iQ10,iO3]   = \
                                     COMPS_DIFF['LULUC_Nat'][temp][Q10][O3[1]][land_type][gcm][iyear,region_mask].sum()

                                LAND_LANDPTS[0,itemp,itype,iyear,igcm,iQ10,iO3,:]        = \
                                     COMPS_DIFF['LULUC_CCS'][temp][Q10][O3[1]][land_type][gcm][iyear,:]

                                LAND_LANDPTS[1,itemp,itype,iyear,igcm,iQ10,iO3,:]        = \
                                     COMPS_DIFF['LULUC_Nat'][temp][Q10][O3[1]][land_type][gcm][iyear,:]

    # Save by GCM
    for igcm in range(nGCMs):
        gcm=GCMs[igcm]

        FILE_LAND_NUMPY    = LAND_npy_DIR+FILE_PART+gcm+'.npy'
        FILE_LAND_NUMPY2   = LAND_npy_DIR+FILE_PART2+gcm+'.npy'

        np.save(FILE_LAND_NUMPY,LAND_REGION[:,:,:,:,:,igcm,:,:])
        np.save(FILE_LAND_NUMPY2,LAND_LANDPTS[:,:,:,:,igcm,:,:,:])

    # No plotting if saving to NPP numpy files
    print('Stopping without plots') 
    quit()

elif PROC_STEP == 'Use':

    # Need to create dummy datasets for CH4, 
    # Full set: COMPS_keys    = [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]
    # CEH  set: COMPS_keys    = [ 'CTL', 'LULUC_CCS', 'LULUC_Nat' ]
    if PLATFORM == 'CEH':
        # COMPS['CH4']         = deepcopy(COMPS['CTL'])
        # COMPS['Coupled_CCS'] = deepcopy(COMPS['LULUC_CCS'])
        # COMPS['Coupled_Nat'] = deepcopy(COMPS['LULUC_Nat'])
        COMPS['Coupled_CCS'] = deepcopy(COMPS['CH4'])
        COMPS['Coupled_Nat'] = deepcopy(COMPS['CH4'])

    COMPS['LULUC_opt']   = deepcopy(COMPS['LULUC_CCS'])
    COMPS['Coupled_opt'] = deepcopy(COMPS['Coupled_CCS'])

    for comp in COMPS_keys_all:
        CH4 = COMPS[comp]['config'].split('_')[0]

        pools_pkl = deepcopy(pools)
        if comp in [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]:
            del pools_pkl[pools_pkl.index('BECCS_ScaleFactor')]
        if comp in [ 'CTL', 'CH4' ]:
            for pool in SOIL_pools: del pools_pkl[pools_pkl.index(pool)]

        for pool in LAND_TYPES+AREA_pools: del pools_pkl[pools_pkl.index(pool)]

        for iTEMP in range(nTEMPs): 
            TEMP=TEMPs[iTEMP]
            TEMP_year=TEMP_years[iTEMP]
            TEMP_tag=TEMP
            COMPS[comp][TEMP]={}
            for iQ10 in range(nQ10s):
                Q10 = Q10_exps[iQ10]
                COMPS[comp][TEMP][Q10]={}
                for iO3 in range(nO3s):
                    O3 = Ozone_exps[iO3]
                    COMPS[comp][TEMP][Q10][O3[1]]={ pool:{} for pool in pools  }  #+land_opt_pools }
                    config=COMPS[comp]['config'].replace('OZONESUB',O3[1])
                    runid=COMPS[comp]['runid'].replace('OZONESUB',O3[0])
                    print('GCM Input:',comp,TEMP,TEMP_year,Q10,config,runid)
                    for igcm in range(nGCMs):
                        gcm=GCMs[igcm]

                        # Read from pickle file
                        COMPS_FILE     = COMPS_npy_DIR+comp+'_'+CH4+'_'+O3[1]+'_'+Q10+'_'+gcm+'_'+TEMP_tag+'_processed.pkl'
                        COMPS_FILE_ID  = open(COMPS_FILE,'rb')
                        print('Reading from:'+COMPS_FILE)
                        PKL_LOAD       = pickle.load(COMPS_FILE_ID)

                        for ipool, pool in enumerate(pools_pkl):
                            COMPS[comp][TEMP][Q10][O3[1]][pool][gcm] = PKL_LOAD[ipool]
    
                            if (pool in [ 'BECCS_productivity_JULES','BECCS_productivity' ]) and DEBUG == "Y":
                                for iregion in REGION_idx:
                                    region            = REGION_dict['Name'][iregion]
                                    region_mask       = REGIONS_1D==(iregion+1)
                                    if region=='Global': region_mask[:]=True
                                    if region=='International Transportation': region_mask[:]=False
    
                                    if len(COMPS[comp][TEMP][Q10][O3[1]][pool][gcm][region_mask]) > 0:
                                        print(comp,TEMP,Q10,O3[1],pool,gcm,region, \
                                            np.min(COMPS[comp][TEMP][Q10][O3[1]][pool][gcm][region_mask]), \
                                            np.max(COMPS[comp][TEMP][Q10][O3[1]][pool][gcm][region_mask]) )
                                    else:
                                        print(comp,TEMP,Q10,O3[1],pool,gcm,region)

                        COMPS_FILE_ID.close()
                        del PKL_LOAD

    # Input land use data from numpy file
    # Saved by GCM
    for igcm in range(nGCMs):
        gcm=GCMs[igcm]

        FILE_LAND_NUMPY    = LAND_npy_DIR+FILE_PART+gcm+'.npy'
        print('Reading from: '+FILE_LAND_NUMPY)
        LAND_ID      = open(FILE_LAND_NUMPY,'rb')
        LAND_FILE    = np.load(LAND_ID)
        LAND_ID.close()
        LAND_REGION[:,:,:,:,:,igcm,:,:] = LAND_FILE

        FILE_LAND_NUMPY2   = LAND_npy_DIR+FILE_PART2+gcm+'.npy'
        print('Reading from: '+FILE_LAND_NUMPY2)
        LAND_ID      = open(FILE_LAND_NUMPY2,'rb')
        LAND_FILE    = np.load(LAND_ID)
        LAND_ID.close()
        LAND_LANDPTS[:,:,:,:,igcm,:,:,:] = LAND_FILE
        
    del LAND_FILE

# For each IMAGE region, identify year of maximum bioenergy land use
BECCS_year   = np.zeros((nTEMPs,nREGIONs_DICT))
itype        = LAND_TYPES.index('Bioenergy')
      
# Need to assign BECCS_AREA from LAND_LANDPTS
for iCOMP,COMP in enumerate(['LULUC_CCS','LULUC_Nat']):
    for itemp in range(nTEMPs): 
        TEMP=TEMPs[itemp]
        for iQ10 in range(nQ10s):
            Q10 = Q10_exps[iQ10]
            for iO3 in range(nO3s):
                O3 = Ozone_exps[iO3]
                COMPS[COMP]['config'].replace('OZONESUB',O3[1])
                COMPS[COMP]['runid'].replace('OZONESUB',O3[0])
                for igcm in range(nGCMs):
                    gcm=GCMs[igcm]

                    iyear_max = np.where(LAND_REGION[0,itemp,-1,itype,:,igcm,iQ10,iO3] == \
                        LAND_REGION[0,itemp,-1,itype,:,igcm,iQ10,iO3].max())[0]

                    COMPS[COMP][TEMP][Q10][O3[1]]['BECCS_Area'][gcm] = \
                            LAND_LANDPTS[iCOMP,itemp,itype,iyear_max,igcm,iQ10,iO3,:].squeeze()

                    if DEBUG == 'Y':
                        for iregion in REGION_idx:
                            region            = REGION_dict['Name'][iregion]
                            region_mask       = REGIONS_1D==(iregion+1)
                            if region=='Global': region_mask[:]=True
                            if region=='International Transportation': region_mask[:]=False

                            print(iCOMP,COMP,TEMP,Q10,O3,region, \
                                 COMPS[COMP][TEMP][Q10][O3[1]]['BECCS_Area'][gcm][region_mask].sum())

                if COMP == 'LULUC_CCS':
                    print('LULUC_CCS')
                    COMPS['Coupled_CCS'][TEMP][Q10][O3[1]]['BECCS_Area'] = deepcopy(COMPS[COMP][TEMP][Q10][O3[1]]['BECCS_Area'])
                    COMPS['LULUC_opt'][TEMP][Q10][O3[1]]['BECCS_Area']   = deepcopy(COMPS[COMP][TEMP][Q10][O3[1]]['BECCS_Area'])
                    COMPS['Coupled_opt'][TEMP][Q10][O3[1]]['BECCS_Area'] = deepcopy(COMPS[COMP][TEMP][Q10][O3[1]]['BECCS_Area'])
                if COMP == 'LULUC_Nat':
                    print('LULUC_Nat')
                    COMPS['Coupled_Nat'][TEMP][Q10][O3[1]]['BECCS_Area'] = deepcopy(COMPS[COMP][TEMP][Q10][O3[1]]['BECCS_Area'])

# From PLOT_MitigationOptions_synthesis_RRv1.4.py
################################################################################

delCStores              = { TEMP: {} for TEMP in TEMPs }
delCStores_unc          = { TEMP: {} for TEMP in TEMPs }

delCStores_regions      =  { TEMP: { region : { Q10: {} for Q10 in Q10_exps } 
                                     for region in REGION_dict['Name'] } 
                                     for TEMP in TEMPs }
delCStores_regions_unc  =  { TEMP: { region : [] for region in REGION_dict['Name'] } 
                             for TEMP in TEMPs }

SummaryFactExps         = ['CTL','CH4','LULUC_opt','Coupled_opt']
FactExps                = ['CTL','CH4',
                           'LULUC_CCS','LULUC_Nat','LULUC_opt',
                           'Coupled_CCS','Coupled_Nat','Coupled_opt']

regional_DF_columns     = [ 'CH4_mit', 'NatLand_mit', 'CCS_mit', 'LULUC_mit' ]
nREGcolumns             = len(regional_DF_columns)

StoreSinks              = ['Atmos','Ocean','Land','CCS']
delCStores_SINKS        = { TEMP: { } for TEMP in TEMPs }
delCStores_SINKS_unc    = { TEMP: { exp:[] for exp in FactExps } for TEMP in TEMPs }
delCStores_mitSINKS_unc = { TEMP: { exp:[] for exp in FactExps } for TEMP in TEMPs }

os.system('mkdir -p '+PLOT_DIR+'CSVoutput/')
print('Storing basic stats in: '+PLOT_DIR+'CSVoutput/Mitigation_Summary_Region.csv')

outf       = open(PLOT_DIR+'CSVoutput/Mitigation_Summary_Global.csv','w')
outfreg    = open(PLOT_DIR+'CSVoutput/Mitigation_Summary_Region.csv','w')
outfunc    = open(PLOT_DIR+'CSVoutput/Mitigation_Summary_FullUnc.csv','w')
outfuncreg = open(PLOT_DIR+'CSVoutput/Mitigation_Summary_FullUnc_Reg.csv','w')
outfsink   = open(PLOT_DIR+'CSVoutput/Mitigation_Summary_Sink_FullUnc.csv','w')

#ipdb.set_trace()
for itemp in range(nTEMPs):
    TEMP= TEMPs[itemp]
    scen_list = []
    for iQ10 in range(nQ10s):
        Q10 = Q10_exps[iQ10]
        delCStores[TEMP][Q10]={}
        delCStores_SINKS[TEMP][Q10] = {}
        for iO3 in range(nO3s):
            O3 = Ozone_exps[iO3]
            outf.write('%-30s'%(TEMP+', '+Q10+', '+O3[1]+':\n'))
            outfreg.write('%-30s'%(TEMP+', '+Q10+', '+O3[1]+':\n'))
            outfreg.write('%-30s:'%('Region') + nREGcolumns*'%10s,' % tuple(regional_DF_columns)+'\n')
            
            # Global Emission budgets and mitigation potentials
            # Budgets:
            DataFrame_Tots         = pd.concat( [pd.Series(COMPS[factexp][TEMP][Q10][O3[1]]['Total'], index=GCMs)
                                                           for factexp in FactExps ], axis=1 )
            DataFrame_Tots.columns = [ factexp+'_tot' for factexp in FactExps ]

            # Mitigation (difference from control)
            DataFrame_Mits         = deepcopy(DataFrame_Tots)
            for col in DataFrame_Tots.columns: DataFrame_Mits[col] -= DataFrame_Tots['CTL_tot']
            DataFrame_Mits.columns = [ factexp+'_mit' for factexp in FactExps ]
            
            DF                     = pd.concat([DataFrame_Tots,DataFrame_Mits],axis=1)

            delCStores[TEMP][Q10][O3[1]] = DF
            DF.describe().to_csv(outf,float_format='%10.2f')
            scen_list.append(DF) 
            
            # Global Sinks:
            delCStores_SINKS[TEMP][Q10][O3[1]] = {}
            for exp in FactExps:
                #print(exp)
                delCStores_SINKS[TEMP][Q10][O3[1]][exp] = pd.DataFrame( 
                            { sink: pd.Series({gcm:COMPS[exp][TEMP][Q10][O3[1]][sink][gcm].sum() 
                                for gcm in GCMs } ) for sink in StoreSinks }  ) 
                # Append to full uncertatinty lists:
                delCStores_SINKS_unc[TEMP][exp].append(delCStores_SINKS[TEMP][Q10][O3[1]][exp])
                delCStores_mitSINKS_unc[TEMP][exp].append( delCStores_SINKS[TEMP][Q10][O3[1]][exp] 
                                                          -delCStores_SINKS[TEMP][Q10][O3[1]]['CTL'] )

            # Land uptake on land points
            LULUC_mit_landpts = { gcm: COMPS['LULUC_opt'][TEMP][Q10][O3[1]]['Land'][gcm]
                                     - COMPS['CTL'][TEMP][Q10][O3[1]]['Land'][gcm] for gcm in GCMs}
            CCS_mit_landpts   = {gcm:  COMPS['LULUC_opt'][TEMP][Q10][O3[1]]['CCS'][gcm]
                                     - COMPS['CTL'][TEMP][Q10][O3[1]]['CCS'][gcm] for gcm in GCMs}
            Total_Land_C_med  = (  np.median( [LULUC_mit_landpts[gcm] for gcm in GCMs], axis=0 )
                                 + np.median( [CCS_mit_landpts[gcm] for gcm in GCMs], axis=0 ) )
                         
            # Regional Breakdown of the mitigaiton options
            for iregion in range(REGION_dict['Nregions']):
                region            = REGION_dict['Name'][iregion]
                region_anthFrac   = REGION_dict['AnthroFraction'][iregion]
                #region_map_index = REGION_dict['Index'][iregion]  #Not required as stored in oreder of map index
                region_mask       = REGIONS_1D==(iregion+1)

                if region=='Global':                       region_mask[:]=True
                if region=='International Transportation': region_mask[:]=False

                regional_CH4_mit   = DF['CH4_mit']*region_anthFrac
                regional_Land_mit  = pd.Series({gcm:LULUC_mit_landpts[gcm][region_mask].sum() for gcm in GCMs} ) 
                regional_CCS_mit   = pd.Series({gcm:CCS_mit_landpts[gcm][region_mask].sum() for gcm in GCMs } ) 
                regional_LULUC_mit = regional_CCS_mit+regional_Land_mit
                DFreg              = pd.concat([regional_CH4_mit,regional_Land_mit,regional_CCS_mit,
                                                regional_LULUC_mit], axis=1) 
                DFreg.columns      = regional_DF_columns
                delCStores_regions[TEMP][region][Q10][O3[1]] = DFreg
                delCStores_regions_unc[TEMP][region].append(
                                       delCStores_regions[TEMP][region][Q10][O3[1]] )
                outfreg.write('%-30s: '%(region))
                DFreg.describe()[5:6].to_csv(outfreg,header=False,float_format='%10.2f',index=False)
                #print(region, DFreg.describe()[5:6])
    
    # compend Global Carbon Stores, full uncertainty to DataFrame and output to csv
    delCStores_unc[TEMP] = pd.concat(scen_list)
    outfunc.write(TEMP+': \n')
    outfunc.write('Global: \n')
    delCStores_unc[TEMP].describe().to_csv(outfunc,float_format='%10.2f')
    
    # Compend the global carbon stores by pool to DataFrame and output to csv 
    delCStores_SINKS_unc[TEMP] = pd.concat( { EXP: pd.concat(delCStores_SINKS_unc[TEMP][EXP]) 
                                                  for EXP in FactExps }, axis=1)
    delCStores_mitSINKS_unc[TEMP] = pd.concat( { EXP: pd.concat(delCStores_mitSINKS_unc[TEMP][EXP]) 
                                                     for EXP in FactExps }, axis=1)
    outfsink.write(TEMP+': \n')
    outfsink.write('Total Sink: \n')
    for EXP in FactExps: 
        outfsink.write('%16s'%(EXP+': \n'))
        delCStores_SINKS_unc[TEMP][EXP].describe().to_csv(outfsink,float_format='%10.2f')
    outfsink.write('Mitigation Potential: \n')
    for EXP in FactExps: 
        outfsink.write('%16s'%(EXP+': \n'))
        delCStores_mitSINKS_unc[TEMP][EXP].describe().to_csv(outfsink,float_format='%10.2f')

    # Compend the regional breakdown into DataFrame: 
    delCStores_regions_unc[TEMP] = pd.concat({ region: pd.concat(delCStores_regions_unc[TEMP][region])
                                                    for region in REGION_dict['Name'] }, axis=1 )
    # Output regional breakdown to csv file 
    outfuncreg.write(TEMP+': \nMedian of GCMs:\n')
    outfuncreg.write('%-30s '%('Region,')+nREGcolumns*'%9s,' % tuple(regional_DF_columns)+'\n')
    for region in REGION_dict['Name']: 
        outfuncreg.write('%-30s '%(region+','))
        delCStores_regions_unc[TEMP][region].describe()[5:6].to_csv(outfuncreg,float_format='%10.2f',
                                                                         header=False,index=False)
    outfuncreg.write('\n\n\n25% of GCMs: \n')
    outfuncreg.write('%-30s '%('Region,')+nREGcolumns*'%9s,' % tuple(regional_DF_columns)+'\n')
    for region in REGION_dict['Name']: 
        outfuncreg.write('%-30s '%(region+','))
        delCStores_regions_unc[TEMP][region].describe()[4:5].to_csv(outfuncreg,float_format='%10.2f',
                                                                         header=False,index=False)
    outfuncreg.write('\n\n\n75% of GCMs: \n')
    outfuncreg.write('%-30s '%('Region,')+nREGcolumns*'%9s,' % tuple(regional_DF_columns)+'\n')
    for region in REGION_dict['Name']: 
        outfuncreg.write('%-30s '%(region+','))
        delCStores_regions_unc[TEMP][region].describe()[6:7].to_csv(outfuncreg,float_format='%10.2f',
                                                                         header=False,index=False)

outfunc.close()
outfuncreg.close()
outfsink.close()
outf.close()
outfreg.close()

# BECCS additional analysis
# Produce output for BECCS scale factors (BECCS_multiplier) from 1 to 4

EXTRA_pools_1   = [ 'BECCS_productivity_JULES','BECCS_productivity',    'BECCS_ScaleFactor' ]
EXTRA_pools_tx1 = [ 'BECCS Productivity',      'Required Productivity', 'Required Scale Factor' ]
EXTRA_pools_2   = [ 'BECCS_Area',    'CCS',     'Land' ]
EXTRA_pools_tx2 = [ 'Area of BECCS', 'BECCS C', 'Land C' ]
EXTRA_pools     = EXTRA_pools_1   + EXTRA_pools_2
EXTRA_pools_txt = EXTRA_pools_tx1 + EXTRA_pools_tx2
#EXTRA_pools    = [ 'BECCS_Area' ]
nEXTRA          = len(EXTRA_pools)

BECCS_SFs   = [ 1, 2, 3, 4 ]
#BECCS_SFs   = [ 1, 3 ] 
nBECCS_SFs  = len(BECCS_SFs)

EXTRA_global  = { TEMP: { BECCS_SF: { Q10: {} for Q10 in Q10_exps } 
                          for BECCS_SF in BECCS_SFs }
                          for TEMP in TEMPs }

EXTRA_region  = { TEMP: { BECCS_SF: { pool : { region: []
                          for region in REGION_dict['Name'] }
                          for pool in EXTRA_pools } 
                          for BECCS_SF in BECCS_SFs }
                          for TEMP in TEMPs }

EXTRA_output  = { TEMP: { BECCS_SF: []
                          for BECCS_SF in BECCS_SFs }
                          for TEMP in TEMPs }


comp    = 'LULUC_opt'
compCCS = comp.replace('opt','CCS')
compNat = comp.replace('opt','Nat')

for itemp in range(nTEMPs): 
    TEMP    = TEMPs[itemp]
    for BECCS_SF in BECCS_SFs:
        BECCS_multiplier = float(BECCS_SF)
        for iQ10 in range(nQ10s):
            Q10 = Q10_exps[iQ10]
            EXTRA_global[TEMP][BECCS_SF][Q10] = {Ozone_exps[iO3][1]: {} for iO3 in range(nO3s) }
            for iO3 in range(nO3s):
                O3 = Ozone_exps[iO3]
                EXTRA_global[TEMP][BECCS_SF][Q10][O3[1]] = {pool: {} for pool in EXTRA_pools}
                for igcm in range(nGCMs):
                    gcm=GCMs[igcm]

                    # create mask, CCS>Nat = 1; CCS==Nat = 0; CCS<Nat = -1 
                    difference = ( ( COMPS[compCCS][TEMP][Q10][O3[1]]['Land'][gcm]
                                   + COMPS[compCCS][TEMP][Q10][O3[1]]['CCS'][gcm]*BECCS_multiplier )
                                 - ( COMPS[compNat][TEMP][Q10][O3[1]]['Land'][gcm]
                                   + COMPS[compNat][TEMP][Q10][O3[1]]['CCS'][gcm] )  )

                    flag_mask = (difference/np.abs(difference)).astype(int)
                    flag_mask[difference==0.] = 0
                    
                    # Calculate the required scale factor for BECCS to become viable:
                    # What scale factor is required for the CCS to be greater than the 
                    # benefit of returning land to natural vegetaiton?
                    req_scale_factor = ( COMPS[compNat][TEMP][Q10][O3[1]]['Land'][gcm]
                                       - COMPS[compCCS][TEMP][Q10][O3[1]]['Land'][gcm])
                    ccs_mask         = (MAXBIOFRAC_1D>BIOFRAC_minimum_threshold) \
                         & (COMPS[compCCS][TEMP][Q10][O3[1]]['CCS'][gcm]>CCS_minimum_threshold)
                    req_scale_factor[ccs_mask] /= COMPS[compCCS][TEMP][Q10][O3[1]]['CCS'][gcm][ccs_mask]
                    req_scale_factor[~ccs_mask] = -1e20
                    COMPS[comp][TEMP][Q10][O3[1]]['BECCS_ScaleFactor'][gcm] = \
                         np.ma.masked_equal(req_scale_factor,-1e20)
                    if DEBUG == "Y": print(TEMP,BECCS_SF,Q10,O3,gcm,'BECCS_ScaleFactor', \
                        np.min(COMPS[comp][TEMP][Q10][O3[1]]['BECCS_ScaleFactor'][gcm]),
                        np.min(COMPS[comp][TEMP][Q10][O3[1]]['BECCS_ScaleFactor'][gcm]))

                    # Substitute Nat for CCS in land/CCS arrays where Nat is the prefered choice
                    for pool in EXTRA_pools:
                        if DEBUG == "Y": print(comp,TEMP,Q10,O3,pool,gcm)

                        if pool in EXTRA_pools_2:
                            COMPS[comp][TEMP][Q10][O3[1]][pool][gcm][flag_mask!=-1] = \
                                deepcopy(COMPS[compCCS][TEMP][Q10][O3[1]][pool][gcm][flag_mask!=-1])
                            COMPS[comp][TEMP][Q10][O3[1]][pool][gcm][flag_mask==-1] = \
                                deepcopy(COMPS[compNat][TEMP][Q10][O3[1]][pool][gcm][flag_mask==-1])

                        # Subtract "CTL" for 'Land' pool to get mitigation potential
                        if pool == 'Land':
                            EXTRA_global[TEMP][BECCS_SF][Q10][O3[1]][pool][gcm] = \
                                COMPS[comp][TEMP][Q10][O3[1]][pool][gcm]-COMPS['CTL'][TEMP][Q10][O3[1]][pool][gcm]
                            print('Land',np.sum(EXTRA_global[TEMP][BECCS_SF][Q10][O3[1]][pool][gcm]), \
                                         np.sum(COMPS[comp][TEMP][Q10][O3[1]][pool][gcm]), \
                                         np.sum(COMPS['CTL'][TEMP][Q10][O3[1]][pool][gcm]))
                        elif pool == 'CCS':
                            EXTRA_global[TEMP][BECCS_SF][Q10][O3[1]][pool][gcm] = deepcopy( \
                               COMPS[comp][TEMP][Q10][O3[1]][pool][gcm]*BECCS_multiplier )
                        else:
                            EXTRA_global[TEMP][BECCS_SF][Q10][O3[1]][pool][gcm] = deepcopy( \
                               COMPS[comp][TEMP][Q10][O3[1]][pool][gcm] )

                        if pool == 'BECCS_Area' and DEBUG == 'Y':
                            print(TEMP,BECCS_SF,Q10,O3,pool,len(flag_mask[flag_mask ==-1]),
                                  len(flag_mask[flag_mask ==0]),len(flag_mask[flag_mask ==1]), \
                                  COMPS[compCCS][TEMP][Q10][O3[1]]['Land'][gcm].sum(), \
                                  COMPS[compCCS][TEMP][Q10][O3[1]]['CCS'][gcm].sum()*BECCS_multiplier, \
                                  COMPS[compNat][TEMP][Q10][O3[1]]['Land'][gcm].sum(), \
                                  COMPS[compNat][TEMP][Q10][O3[1]]['CCS'][gcm].sum())

                            for i in range(len(COMPS[comp][TEMP][Q10][O3[1]][pool][gcm][:])):
                                iregion = REGIONS_1D[i]
                                print(i,iregion,REGION_dict['Name'][iregion-1],flag_mask[i],difference[i], \
                                   COMPS[comp][TEMP][Q10][O3[1]][pool][gcm][i], \
                                   COMPS[compCCS][TEMP][Q10][O3[1]][pool][gcm][i], \
                                   COMPS[compNat][TEMP][Q10][O3[1]][pool][gcm][i],
                                   LAND_LANDPTS[0,itemp,itype,:,igcm,iQ10,iO3,i].max())

                    if DEBUG == "Y":
                        for pool in EXTRA_pools:
                            for iregion in REGION_idx:
                                region            = REGION_dict['Name'][iregion]
                                region_mask       = REGIONS_1D==(iregion+1)
                                if region=='Global': region_mask[:]=True
                                if region=='International Transportation': region_mask[:]=False

                                print(comp,TEMP,BECCS_SF,Q10,O3,pool,region)
                                if pool == 'BECCS_ScaleFactor':
                                    print(np.max(EXTRA_global[TEMP][BECCS_SF][Q10][O3[1]][pool][gcm][region_mask]))
                                elif pool in EXTRA_pools_1:
                                    print(np.max(EXTRA_global[TEMP][BECCS_SF][Q10][O3[1]][pool][gcm][region_mask]), \
                                       np.max(COMPS[compCCS][TEMP][Q10][O3[1]][pool][gcm][region_mask]), \
                                       np.max(COMPS[compNat][TEMP][Q10][O3[1]][pool][gcm][region_mask]))
                                elif pool in EXTRA_pools_2:
                                    print(np.sum(EXTRA_global[TEMP][BECCS_SF][Q10][O3[1]][pool][gcm][region_mask]), \
                                       np.sum(COMPS[compCCS][TEMP][Q10][O3[1]][pool][gcm][region_mask]), \
                                       np.sum(COMPS[compNat][TEMP][Q10][O3[1]][pool][gcm][region_mask]))

print("\n")
print('EXTRA_global: loaded')

for itemp in range(nTEMPs): 
    TEMP    = TEMPs[itemp]
    for BECCS_SF in BECCS_SFs:
        BECCS_multiplier = float(BECCS_SF)
        for pool in EXTRA_pools:
       
            # Regional Breakdown
            for iregion in REGION_idx:
                region            = REGION_dict['Name'][iregion]
                region_mask       = REGIONS_1D==(iregion+1)

                if region=='Global':                       region_mask[:]=True
                if region=='International Transportation': region_mask[:]=False

                for iQ10 in range(nQ10s):
                    Q10 = Q10_exps[iQ10]
                    for iO3 in range(nO3s):
                        O3 = Ozone_exps[iO3]
                        for igcm in range(nGCMs):
                            gcm=GCMs[igcm]

                            if pool in EXTRA_pools_1:
                                REGION_DATA  = EXTRA_global[TEMP][BECCS_SF][Q10][O3[1]][pool][gcm][region_mask]
                                EXTRA_region[TEMP][BECCS_SF][pool][region] = np.append( \
                                   EXTRA_region[TEMP][BECCS_SF][pool][region], \
                                   REGION_DATA[REGION_DATA > 0.0])
                                if DEBUG == "Y":
                                    print(TEMP,BECCS_SF,region,Q10,O3[1],pool,gcm, \
                                       len(REGION_DATA),REGION_DATA.min(),REGION_DATA.max())

                            if pool in EXTRA_pools_2:
                                EXTRA_region[TEMP][BECCS_SF][pool][region] = np.append( \
                                   EXTRA_region[TEMP][BECCS_SF][pool][region], \
                                   np.sum(EXTRA_global[TEMP][BECCS_SF][Q10][O3[1]][pool][gcm][region_mask]))
                                if DEBUG == "Y":
                                    print(TEMP,BECCS_SF,region,Q10,O3[1],pool,gcm, \
                                       EXTRA_global[TEMP][BECCS_SF][Q10][O3[1]][pool][gcm][region_mask].sum())

print("\n")
print('EXTRA_regional: loaded')
FILE_CSV    = PLOT_DIR+'BECCS_Extra_'+FILE_PART3+'.csv'
print('Writing to: '+FILE_CSV)
OUT_FID     = open(FILE_CSV,'w')

for itemp in range(nTEMPs): 
    TEMP        = TEMPs[itemp]
    OUT_TEXT1   = ('%-10s') % TEMP
    OUT_FID.write(OUT_TEXT1+'\n')

    OUT_TEXT1   = ',,,,,'
    OUT_TEXT2   = 'Region,Max BECCS Area,BECCS_Productivity,Required_Productivity,Required Scale Factor,'
    for BECCS_SF in BECCS_SFs:
        OUT_TEXT1   = OUT_TEXT1 + '%s' % ('BECCS scale factor = '+str(BECCS_SF))
        for pool in EXTRA_pools_tx2:
            OUT_TEXT1   = OUT_TEXT1 + ','
            OUT_TEXT2   = OUT_TEXT2 + '%s,' % pool
    OUT_FID.write(OUT_TEXT1+'\n')    
    OUT_FID.write(OUT_TEXT2+'\n')    
    
    # Regional Breakdown
    for iregion in REGION_idx:
        region            = REGION_dict['Name'][iregion]
        region_mask       = REGIONS_1D==(iregion+1)
        OUT_TEXT1   = ('%s,%.1f,') % (region, \
            np.sum(LAND_LANDPTS[0,itemp,itype,iyear_max,:,:,:,region_mask])/nFACTORIAL)
        for BECCS_SF in BECCS_SFs:
            
            if (BECCS_SF == 1):
                for pool in EXTRA_pools_1:
                    print(region,BECCS_SF,pool,len(EXTRA_region[TEMP][BECCS_SF][pool][region]))
                    if len(EXTRA_region[TEMP][BECCS_SF][pool][region]) > 0:
                        PERCENT     = np.percentile(np.array((EXTRA_region[TEMP][BECCS_SF][pool][region])),[50.0,10.0,90.0])
                        OUT_TEXT1   = OUT_TEXT1 + '%.2f (%.2f-%.2f),' % \
                            (PERCENT[0],PERCENT[1],PERCENT[2])
                    else:
                        OUT_TEXT1   = OUT_TEXT1 + '-,'
            
            for pool in EXTRA_pools_2:
                print(region,BECCS_SF,pool,len(EXTRA_region[TEMP][BECCS_SF][pool][region]))
                if len(EXTRA_region[TEMP][BECCS_SF][pool][region]) > 0:
                    PERCENT     = np.percentile(np.array((EXTRA_region[TEMP][BECCS_SF][pool][region])),[50.0,25.0,75.0])
                    OUT_TEXT1   = OUT_TEXT1 + '%.2f (%.2f-%.2f),' % \
                        (PERCENT[0],PERCENT[1],PERCENT[2])
                else:
                    OUT_TEXT1   = OUT_TEXT1 + '-,'

        OUT_FID.write(OUT_TEXT1+'\n')

OUT_FID.close()
               
print("\n")
print('EXTRA_output')
