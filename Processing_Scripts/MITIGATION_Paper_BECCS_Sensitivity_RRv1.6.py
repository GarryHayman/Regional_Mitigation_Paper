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
PLOT_TAG      = optional_argparse('-plottag', 'BECCS_Sensitivity')
VERSION       = optional_argparse('-BECCS','max_bioenergy')
BECCS_NPP     = optional_argparse('-BECCS_NPP','N')
PROC_STEP     = optional_argparse('-BECCS_NPY','Save')
CCS_minimum_threshold     = float(optional_argparse('CCS_min','1e-4'))
BIOFRAC_minimum_threshold = float(optional_argparse('CCSfrac_min','1e-2'))
PLOT_FIGURE   = False

# Directories containing JULES output and plot output directories:
if PLATFORM == 'JASMIN':
    HOME_DIR      = '/gws/nopw/j04/clifftop/SYNTHESIS/PostReview_output/'
    DATA_DIR      = HOME_DIR
    ANCILS_DIR    = '/gws/nopw/j04/clifftop/COMMON_DATA/ANCILS/'
    PLOT_DIR      = optional_argparse('-plotdir',DATA_DIR+'plots/'+PLOT_TAG+'/')
    BECCS_npy_DIR = '/gws/nopw/j04/jules/aharper/PYTHON/SYNTHESIS/npy_files/'+VERSION+'/'
#   BECCS_npy_NPP = '/gws/nopw/j04/jules/ghayman/PYTHON/SYNTHESIS/npy_files/max_bioenergy_NPP/'
#   SCNPP_npy_DIR = '/gws/nopw/j04/jules/ghayman/PYTHON/SYNTHESIS/npy_files/Carbon/'
    BECCS_npy_NPP = '/work/scratch/garr_output/SYNTHESIS/npy_files/max_bioenergy_NPP/'
    SCNPP_npy_DIR = '/work/scratch/garr_output/SYNTHESIS/npy_files/Carbon/'
    COMPS_npy_DIR = '/work/scratch/garr_output/SYNTHESIS/npy_files/processed_output/'

    Q10_exps      = [ 'lowQ10', 'highQ10' ]
    Ozone_exps    = [['L','lowO3'], ['H','highO3'], ]

    COMPS_keys    = [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]

elif PLATFORM == 'CEH':
    HOME_DIR      = '/prj/CLIFFTOP/SYNTHESIS/'
    PLOT_DIR      = '/data/grp/eow/garr/Projects/CLIFFTOP/Paper_Synthesis/PLOTS_LAND_COVER/CEH/'
    DATA_DIR      = HOME_DIR+'Review_Response_Check/GCM_Output/'
    ANCILS_DIR    = HOME_DIR+'Land_Cover/'
    BECCS_npy_DIR = HOME_DIR+'Review_Response_Check/'+VERSION+'/'
    BECCS_npy_NPP = HOME_DIR+'Review_Response_Check/max_bioenergy_NPP/'
    SCNPP_npy_DIR = HOME_DIR+'Review_Response_Check/Carbon/'
    COMPS_npy_DIR = HOME_DIR+'Review_Response_Check/processed_output/'

    Q10_exps      = [ 'highQ10' ]
    Ozone_exps    = [['L','lowO3']]

    COMPS_keys    = [ 'CTL', 'LULUC_CCS', 'LULUC_Nat' ]

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

# Directories containing JULES output and plot output directories:
DATA_DIR      = optional_argparse('-data_dir', DATA_DIR)
PLOT_DIR      = optional_argparse('-plotdir',DATA_DIR+'plots/'+PLOT_TAG+'/')
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

land_pools    = [ 'Land','CS','CV','WP','CCS', 'BE_Harvest' ]
pools         = ['Total','Atmos','Ocean',]+land_pools \
               +['Soil_Carbon', 'Soil_Carbon_gb', 'NPP_gb' ] \
               +['BECCS_productivity', 'BECCS_ScaleFactor', \
                 'BECCS_productivity_JULES','BECCS_productivity_NPP' ]
              #+['BECCS_productivity','Harvest','kappaHarvest','kappamefficHarvest']


dzsoil        = np.array([ 0.05,0.08408964,0.11397535,0.14142136,0.16718508,0.19168293, \
                  0.21517585,0.23784142,0.25980762,0.28117066,0.30200527, \
                  0.32237098,0.34231625,0.36188121 ])

nLAYERs       = 7
soil_depth    = dzsoil[0:nLAYERs].sum()
 
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


# Indices of pfts for aggregation to trees, grasses, shrubs
TREE_INDICES  = [ 0,1,2,3,4 ]
GRASS_INDICES = [ 5,8 ] # natural grasses
CROP_INDICES  = [ 6,9 ]
PAST_INDICES  = [ 7,10 ]
SHRUB_INDICES = [ 11,12 ]

#  BECCS_multiplier not required for BECCS sensitivty routine
### BECCS_multiplier = float(optional_argparse('-beccs_multiplier','1'))
### print(BECCS_multiplier)

# Directory of Ocean Uptake data:
OCEAN_UPTAKE_DIR = optional_argparse('-ocean_uptake_dir',DATA_DIR) 
OCEAN_START_YEAR = 1850
print("Ocean Uptake data from: " + OCEAN_UPTAKE_DIR)

REGION_dict      = data_info.REGION_DICTIONARIES()[subregions]

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
                        if ('CCS' in comp or 'Nat' in comp) and (BECCS_NPP == 'Y'):

                            # Skip if files already exist for bioenergy, NPP and soil carbon
                            BECCS_NPP_FILE = BECCS_npy_NPP+comp+'_'+CH4+'_'+O3[1]+'_'+Q10+'_'+gcm+'_'+TEMP_tag+'_NPP.npy'
                            SC_NPP_FILE    = SCNPP_npy_DIR+comp+'_'+CH4+'_'+O3[1]+'_'+Q10+'_'+gcm+'_'+TEMP_tag+'_SCNPP.pkl'

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
                        for pool in land_pools:
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
            del pools_pkl[pools_pkl.index('NPP_gb')]
            del pools_pkl[pools_pkl.index('Soil_Carbon')]
            del pools_pkl[pools_pkl.index('Soil_Carbon_gb')]

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

    # No plotting if saving to NPP numpy files
    print('Stopping without plots') 
    quit()

elif PROC_STEP == 'Use':

    # Need to create dummy datasets for CH4, 
    # Full set: COMPS_keys    = [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]
    # CEH  set: COMPS_keys    = [ 'CTL', 'LULUC_CCS', 'LULUC_Nat' ]
    if PLATFORM == 'CEH':
        COMPS['CH4']         = deepcopy(COMPS['CTL'])
        COMPS['Coupled_CCS'] = deepcopy(COMPS['LULUC_CCS'])
        COMPS['Coupled_Nat'] = deepcopy(COMPS['LULUC_Nat'])

    COMPS['LULUC_opt']   = deepcopy(COMPS['LULUC_CCS'])
    COMPS['Coupled_opt'] = deepcopy(COMPS['Coupled_CCS'])

    for comp in COMPS_keys_all:
        CH4 = COMPS[comp]['config'].split('_')[0]

        pools_pkl = deepcopy(pools)
        if comp in [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]:
            del pools_pkl[pools_pkl.index('BECCS_ScaleFactor')]
        if comp in [ 'CTL', 'CH4' ]:
            del pools_pkl[pools_pkl.index('NPP_gb')]
            del pools_pkl[pools_pkl.index('Soil_Carbon')]
            del pools_pkl[pools_pkl.index('Soil_Carbon_gb')]

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
                        COMPS_FILE_ID.close()
                        del PKL_LOAD
                
# BECCS productivity map
comp = 'LULUC_opt'

BECCS_colours = ['#bf812d','#f6e8c3','#35978f','#01665e',]  # Switch colours at 2 
BECCS_CLEVELS = [0,1.,2.0,3.0,5.]
BECCS_CLEVELS = [0,0.5,1.,1.5,2.0,2.5,3.0,5.,10.]

harvest_label = 'Required bioenergy crop yield (tonnesDM ha$^{-1}$ yr$^{-1}$)'
#harvest_CLEVELS = [0.,5.,10.,20.,30.]
harvest_CLEVELS  = [0.1,1.0,2.5,5.,7.5,10.,12.5,15.,20.,]

plotvarnames  = ['BECCS_ScaleFactor','BECCS_productivity']
plotvars = {
    'BECCS_ScaleFactor':  {'cbar_lab':'BECCS Scale Factor','ftag':'BECCS_Scale_Factor','clevels':BECCS_CLEVELS},
    'BECCS_productivity': {'cbar_lab':'BECCS Productivity','ftag':'BECCS_Productivity','clevels':harvest_CLEVELS},
    'BECCS_productivity_JULES': {'cbar_lab':'BECCS Productivity','ftag':'BECCS_Productivity','clevels':harvest_CLEVELS},
            }
npvars = len(plotvarnames)

#ipdb.set_trace()
if PLOT_FIGURE and PROC_STEP == 'Use':

# Figure Block (1)
# Plot maps of the median threshold BECCS scalefactor and productivity.
# Single map plots in sub-folders:
# (a) IndividualConfig_Maps: plot for each uncertainty configuration
# (b) Ozone Maps:            plot median for each ozone uncertainty configuration           
# (c) CH4_Q10_Maps           plot median for each CH4 Q10 uncertainty configuration
# (d) Summary_Maps           plot median of entire uncertainty configuration           
# pdf and png versions of each

    for ivar in range(npvars):
        var = plotvarnames[ivar]
        # Plot for each uncertainty configuration
        os.system('mkdir -p '+PLOT_DIR+'IndividualConfig_Maps/'+var)
        for iTEMP in range(nTEMPs):
            TEMP=TEMPs[iTEMP]
            for iQ10 in range(nQ10s):
                Q10 = Q10_exps[iQ10]
                for iO3 in range(nO3s):
                    #ipdb.set_trace()
                    O3 = Ozone_exps[iO3]
                    scen_median_DF = pd.DataFrame(COMPS[comp][TEMP][Q10][O3[1]][var]).median(axis=1)
                    scen_median_plotdata =  scen_median_DF.values[land_index]
                    GREY_MASK=np.ones_like(scen_median_plotdata)
                    GREY_MASK[(scen_median_plotdata<=0)&(land_index.mask==False)] = 0
                    GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
                    scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata, 
                        mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01)|np.isnan(scen_median_plotdata) )
                    print('Mean Productivity (l.  781): ',var,TEMP,Q10,O3[1],np.mean(scen_median_plotdata))
                    FILENAME = PLOT_DIR+'IndividualConfig_Maps/'+var+'/'+plotvars[var]['ftag']+'_Map_'+TEMP+O3[1]+Q10
                    PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
                        COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,CBAR_LABEL=plotvars[var]['cbar_lab'],
                        CLEVELS=plotvars[var]['clevels'], SET_OVER='#003c30',RESOLUTION='c', MAP_TYPE='Mesh',
                        FONTSIZES=[12,12,14,18],GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7',FILE_PLOT=FILENAME+'.png',)
                        #DISPLAY='N')
                    plt.savefig(FILENAME+'.pdf')
                    plt.close()

        # Plot median for each ozone uncertainty configuration
        os.system('mkdir -p '+PLOT_DIR+'Ozone_Maps/'+var)
        for iTEMP in range(nTEMPs):
            TEMP=TEMPs[iTEMP]
            for iO3 in range(nO3s):
                O3 = Ozone_exps[iO3]
                scen_median_DF = pd.concat(
                    [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[iQ10]][O3[1]][var] )
                      for iQ10 in range(nQ10s)  ],axis=1  ).median(axis=1)
                scen_median_plotdata =  scen_median_DF.values[land_index]
                GREY_MASK=np.ones_like(scen_median_plotdata)
                GREY_MASK[(scen_median_plotdata<=0)&(land_index.mask==False)] = 0
                GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
                scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata, 
                    mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01)|np.isnan(scen_median_plotdata) )
                FILENAME = PLOT_DIR+'Ozone_Maps/'+var+'/'+plotvars[var]['ftag']+'_Map_'+TEMP+O3[1]
                print('Mean Productivity (l.  807): ',var,TEMP,O3[1],np.mean(scen_median_plotdata))
                PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
                    COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,CBAR_LABEL=plotvars[var]['cbar_lab'],
                    CLEVELS=plotvars[var]['clevels'], SET_OVER='#003c30',RESOLUTION='c', MAP_TYPE='Mesh',
                    GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7',FONTSIZES=[12,12,14,18],FILE_PLOT=FILENAME+'.png')
                plt.savefig(FILENAME+'.pdf')
                plt.close()

        # Plot median for each CH4 Q10 uncertainty configuration
        os.system('mkdir -p '+PLOT_DIR+'CH4_Q10_Maps/'+var)
        for iTEMP in range(nTEMPs):
            TEMP=TEMPs[iTEMP]
            for iQ10 in range(nQ10s):
                Q10 = Q10_exps[iQ10]
                scen_median_DF = pd.concat(
                    [ pd.DataFrame( COMPS[comp][TEMP][Q10][Ozone_exps[iO3][1]][var] )
                      for iO3 in range(nO3s)  ],axis=1  ).median(axis=1)
                scen_median_plotdata =  scen_median_DF.values[land_index]
                GREY_MASK=np.ones_like(scen_median_plotdata)
                GREY_MASK[(scen_median_plotdata<=0)&(land_index.mask==False)] = 0
                GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
                scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata, 
                    mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01)|np.isnan(scen_median_plotdata) )
                print('Mean Productivity (l.  830): ',var,TEMP,Q10,np.mean(scen_median_plotdata))
                FILENAME = PLOT_DIR+'CH4_Q10_Maps/'+var+'/'+plotvars[var]['ftag']+'_Map_'+TEMP+Q10
                PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
                    COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,CBAR_LABEL=plotvars[var]['cbar_lab'],
                    CLEVELS=plotvars[var]['clevels'], SET_OVER='#003c30',RESOLUTION='c', MAP_TYPE='Mesh',
                    GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7',FONTSIZES=[12,12,14,18],FILE_PLOT=FILENAME+'.png')
                plt.savefig(FILENAME+'.pdf')
                plt.close()

        # Plot median of entire uncertainty configuration
        # Figure 7 in ESD paper (var = BECCS_ScaleFactor)
        # Note not printed in this version as PLOT_FIGURE is currently False
        os.system('mkdir -p '+PLOT_DIR+'Summary_Maps/'+var)
        for iTEMP in range(nTEMPs):
            TEMP=TEMPs[iTEMP]

            scen_median_DF = pd.concat(
                [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[0]][Ozone_exps[iO3][1]][var] )
                  for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  ).median(axis=1)
            scen_median_plotdata =  scen_median_DF.values[land_index]
            GREY_MASK=np.ones_like(scen_median_plotdata)
            GREY_MASK[(scen_median_plotdata<=0)&(land_index.mask==False)] = 0
            GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
            scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata, 
                mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01)|np.isnan(scen_median_plotdata) )
            print('Mean Productivity (l.  853): ',var,TEMP,np.mean(scen_median_plotdata))

            scen_median_DF = pd.concat(
                [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[iQ10]][Ozone_exps[iO3][1]][var] )
                  for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  ).median(axis=1)
            scen_median_plotdata =  scen_median_DF.values[land_index]
            GREY_MASK=np.ones_like(scen_median_plotdata)
            GREY_MASK[(scen_median_plotdata<=0)&(land_index.mask==False)] = 0
            GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
            scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata, 
                mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01)|np.isnan(scen_median_plotdata) )
            print('Mean Productivity (l.  864): ',var,TEMP,np.mean(scen_median_plotdata))

            FILENAME = PLOT_DIR+'Summary_Maps/'+var+'/'+plotvars[var]['ftag']+'_Map_'+TEMP
            PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
                COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,CBAR_LABEL=plotvars[var]['cbar_lab'],
                CLEVELS=plotvars[var]['clevels'], SET_OVER='#003c30',RESOLUTION='c', MAP_TYPE='Mesh',
                GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7',FONTSIZES=[12,12,14,18],FILE_PLOT=FILENAME+'.png')
            plt.savefig(FILENAME+'.pdf')
            plt.close()

# BECCS scale factor sensitivity # BECCS vs Rest, total contribution, sensitivity
# Calculate the CCS, Natural land, total land and total earth system contributions to AFFEB and MP
#    for a range of values of kappa (BECCS scale factor)
# Store the 25th and 75th percentiles for constructing the shaded regions on the graphs
if True:
    comps = ['LULUC','Coupled']
    BECCS_factors = np.arange(0.4,8.1,0.1)
    nBECCSf = len(BECCS_factors)
    CCS_tot_lims    = { comp: { TEMP: [ np.zeros_like(BECCS_factors),np.zeros_like(BECCS_factors) ] 
                            for TEMP in TEMPs} for comp in comps }
    NatLand_tot_lims = { comp: { TEMP: [ np.zeros_like(BECCS_factors),np.zeros_like(BECCS_factors) ] 
                            for TEMP in TEMPs} for comp in comps }
    Cland_tot_lims  = { comp: { TEMP: [ np.zeros_like(BECCS_factors),np.zeros_like(BECCS_factors) ] 
                            for TEMP in TEMPs} for comp in comps }
    CCS_mit_lims    = { comp: { TEMP: [ np.zeros_like(BECCS_factors),np.zeros_like(BECCS_factors) ] 
                            for TEMP in TEMPs} for comp in comps }
    NatLand_mit_lims = { comp: { TEMP: [ np.zeros_like(BECCS_factors),np.zeros_like(BECCS_factors) ] 
                            for TEMP in TEMPs} for comp in comps }
    Cland_mit_lims  = { comp: { TEMP: [ np.zeros_like(BECCS_factors),np.zeros_like(BECCS_factors) ] 
                            for TEMP in TEMPs} for comp in comps }
    Total_lims  = { comp: { TEMP: [ np.zeros_like(BECCS_factors),np.zeros_like(BECCS_factors) ]
                            for TEMP in TEMPs} for comp in comps+['Linear'] }
    Mitig_lims  = { comp: { TEMP: [ np.zeros_like(BECCS_factors),np.zeros_like(BECCS_factors) ]
                            for TEMP in TEMPs} for comp in comps+['Linear'] }

CTL=COMPS['CTL']
CH4=COMPS['CH4']
#ipdb.set_trace()
for comp in comps:  #['LULUC','Coupled']
    compCCS=comp+'_CCS';  compNat=comp+'_Nat'
    for iBECCSfactor in range(nBECCSf):
        BECCSfactor=BECCS_factors[iBECCSfactor]
        for iTEMP in range(nTEMPs):
            TEMP=TEMPs[iTEMP]
            Total_all = []; Cland_all = []; CCS_all = []; NatLand_all = []
            Mitig_all = []; Clandmit_all = []; CCSmit_all = []; NatLandmit_all = []
            if comp == 'LULUC': Mitig_all_Linear=[]; Total_all_Linear = []
            for iQ10 in range(nQ10s):
                Q10 = Q10_exps[iQ10]
                for iO3 in range(nO3s):
                    O3 = Ozone_exps[iO3]
                    # CCS Dataframes
                    NatLand_CCS    = pd.DataFrame(COMPS[compCCS][TEMP][Q10][O3[1]]['Land'])
                    CCS_scaled_CCS = pd.DataFrame(COMPS[compCCS][TEMP][Q10][O3[1]]['CCS'])* BECCSfactor
                    Cland_CCS      = NatLand_CCS+CCS_scaled_CCS
                    NonLand_CCS = ( pd.Series(COMPS[compCCS][TEMP][Q10][O3[1]]['Atmos']) 
                                   + pd.Series(COMPS[compCCS][TEMP][Q10][O3[1]]['Ocean']) )
                    #Total_CCS   = pd.Series(COMPS[compCCS][TEMP][Q10][O3[1]]['Total'])  
                    # Natural Dataframes
                    NatLand_Nat    = pd.DataFrame(COMPS[compNat][TEMP][Q10][O3[1]]['Land'])
                    CCS_scaled_Nat = pd.DataFrame(COMPS[compNat][TEMP][Q10][O3[1]]['CCS']) * BECCSfactor # This should equal zero
                    Cland_Nat      = NatLand_Nat + CCS_scaled_Nat  
                    #NonLand_Nat = ( pd.Series(COMPS[compNat][TEMP][Q10][O3[1]]['Atmos']) 
                    #               + pd.Series(COMPS[compNat][TEMP][Q10][O3[1]]['Ocean']) )
                    #Total_Nat   = pd.Series(COMPS[compNat][TEMP][Q10][O3[1]]['Total'])  
                    # Control Dataframes:
                    NatLand_ctl    = pd.DataFrame(CTL[TEMP][Q10][O3[1]]['Land'])
                    CCS_scaled_ctl = pd.DataFrame(CTL[TEMP][Q10][O3[1]]['CCS']) * BECCSfactor # This should = zero
                    Cland_ctl      = NatLand_ctl + CCS_scaled_ctl
                    NonLand_ctl = ( pd.Series(CTL[TEMP][Q10][O3[1]]['Atmos']) 
                                   + pd.Series(CTL[TEMP][Q10][O3[1]]['Ocean']) )
                    Total_ctl = pd.Series(CTL[TEMP][Q10][O3[1]]['Total'])
                    #NonLand_mit = NonLand_comp-NonLand_ctl
                    # create mask, CCS>Nat = 1; CCS==Nat = 0; CCS<Nat = -1 
                    flag_mask = (Cland_CCS-Cland_Nat)/np.abs(Cland_CCS-Cland_Nat)
                    flag_mask[Cland_CCS-Cland_Nat==0.] = 0
                    # Optimised/minimum Dataframes, start with CCS then fill with natural where neccessary:
                    NatLand_opt    = deepcopy(NatLand_CCS)
                    CCS_scaled_opt = deepcopy(CCS_scaled_CCS)
                    for gcm in GCMs:
                        NatLand_opt[gcm][flag_mask[gcm]==-1] = NatLand_Nat[gcm][flag_mask[gcm]==-1]
                        CCS_scaled_opt[gcm][flag_mask[gcm]==-1] = CCS_scaled_Nat[gcm][flag_mask[gcm]==-1]
                    Cland_opt = CCS_scaled_opt+NatLand_opt
                    # For non-land components, take the CCS component, this is slightly untrue as there are differences
                    # in the natural CH4 emissions caused by natural/managed land. However, the differences are <1GtC
                    # for the time being at least, we ignore this difference.
                    Total_opt = NonLand_CCS+Cland_opt.sum()
                    CCS_all.append(CCS_scaled_opt.sum())
                    NatLand_all.append(NatLand_opt.sum())
                    Cland_all.append(Cland_opt.sum())
                    Total_all.append(Total_opt)
                    CCSmit_all.append( CCS_scaled_opt.sum()-CCS_scaled_ctl.sum() )
                    NatLandmit_all.append( NatLand_opt.sum()-NatLand_ctl.sum() )
                    Clandmit_all.append( Cland_opt.sum()-Cland_ctl.sum() )
                    Mitig_all.append( Total_opt - Total_ctl )# NonLand_mit + Cland_opt.sum() - Cland_ctl.sum())

                    if comp=='LULUC':
                        Mitig_CH4 = ( pd.Series(CH4[TEMP][Q10][O3[1]]['Total'])
                                     -pd.Series(CTL[TEMP][Q10][O3[1]]['Total']))
                        Total_all_Linear.append( Mitig_CH4 + Total_opt )
                        Mitig_all_Linear.append( Mitig_CH4 + Total_opt - Total_ctl)# NonLand_mit+Cland_opt.sum()-Cland_min.sum())

            CCS_all = pd.concat(CCS_all)
            CCS_tot_lims[comp][TEMP][0][iBECCSfactor]=CCS_all.describe()['25%']
            CCS_tot_lims[comp][TEMP][1][iBECCSfactor]=CCS_all.describe()['75%']
            NatLand_all = pd.concat(NatLand_all)
            NatLand_tot_lims[comp][TEMP][0][iBECCSfactor]=NatLand_all.describe()['25%']
            NatLand_tot_lims[comp][TEMP][1][iBECCSfactor]=NatLand_all.describe()['75%']
            Cland_all = pd.concat(Cland_all)
            Cland_tot_lims[comp][TEMP][0][iBECCSfactor]=Cland_all.describe()['25%']
            Cland_tot_lims[comp][TEMP][1][iBECCSfactor]=Cland_all.describe()['75%']

            Total_all = pd.concat(Total_all)
            Total_lims[comp][TEMP][0][iBECCSfactor]=Total_all.describe()['25%']
            Total_lims[comp][TEMP][1][iBECCSfactor]=Total_all.describe()['75%']
            if comp=='LULUC':
                Total_all_Linear = pd.concat(Total_all_Linear)
                Total_lims['Linear'][TEMP][0][iBECCSfactor] = Total_all_Linear.describe()['25%']
                Total_lims['Linear'][TEMP][1][iBECCSfactor] = Total_all_Linear.describe()['75%']

            CCSmit_all = pd.concat(CCSmit_all)
            CCS_mit_lims[comp][TEMP][0][iBECCSfactor]=CCSmit_all.describe()['25%']
            CCS_mit_lims[comp][TEMP][1][iBECCSfactor]=CCSmit_all.describe()['75%']
            NatLandmit_all = pd.concat(NatLandmit_all)
            NatLand_mit_lims[comp][TEMP][0][iBECCSfactor]=NatLandmit_all.describe()['25%']
            NatLand_mit_lims[comp][TEMP][1][iBECCSfactor]=NatLandmit_all.describe()['75%']
            Clandmit_all = pd.concat(Clandmit_all)
            Cland_mit_lims[comp][TEMP][0][iBECCSfactor]=Clandmit_all.describe()['25%']
            Cland_mit_lims[comp][TEMP][1][iBECCSfactor]=Clandmit_all.describe()['75%']

            Mitig_all = pd.concat(Mitig_all)
            Mitig_lims[comp][TEMP][0][iBECCSfactor]=Mitig_all.describe()['25%']
            Mitig_lims[comp][TEMP][1][iBECCSfactor]=Mitig_all.describe()['75%']
            if comp=='LULUC':
                Mitig_all_Linear = pd.concat(Mitig_all_Linear)
                Mitig_lims['Linear'][TEMP][0][iBECCSfactor]=Mitig_all_Linear.describe()['25%']
                Mitig_lims['Linear'][TEMP][1][iBECCSfactor]=Mitig_all_Linear.describe()['75%']

           #ipdb.set_trace()

#ipdb.set_trace()

if True:
    Atmoscolor = '#ffff99'
    Oceancolor = '#a6cee3'
    CCScolor   = '#e6ab02'
    Landcolor  = '#1b9e77'
    totLULUCcolor='#66a61e'
    CH4color = '#7570b3'
    totcolor = '#d95f02'
    totcolor_lin = '#a6761d'
    CTLcolor = '#666666' # '#80b1d3'

xticks = [0.5,0.6,0.8,1,1.5,2,3,4,5,6] #,6,7,8]
xlims=[0.5,6]

os.system('mkdir -p '+PLOT_DIR+'Sensitivity_Plots/ ')
if PLOT_FIGURE and PROC_STEP == 'Use':

# Figure Block (2)
# All plots in sub-folder: Sensitivity_Plots
# Plots of sensitivity of mitigation potential for different scenario combinations 
# Separate plots for each temperature pathway
# pdf and png versions of each

    # Natural Land uptake and Total Land contribution
    for comp,ylims in zip(comps,([-50,150],[0,250])):
        for TEMP in TEMPs:
            fig,axis = plt.subplots(ncols=1,nrows=1,figsize=[8,5])
            axis.fill_between(BECCS_factors,Cland_tot_lims[comp][TEMP][0],Cland_tot_lims[comp][TEMP][1],
                              color=totLULUCcolor,alpha=0.5,label='Total Land Contribution',zorder=2)
            axis.fill_between(BECCS_factors,NatLand_tot_lims[comp][TEMP][0],NatLand_tot_lims[comp][TEMP][1],
                              color=Landcolor,alpha=0.7,label='Natural Land Uptake',zorder=3)
            axis.fill_between(BECCS_factors,CCS_tot_lims[comp][TEMP][0],CCS_tot_lims[comp][TEMP][1],
                              color=CCScolor,alpha=1.,label='BECCS',zorder=4)
            axis.plot(xlims,[0,0],c='k')
            axis.set_xscale('log')
            axis.set_ylim(ylims)
            axis.set_xlim(xlims)
            axis.set_xticks(xticks)
            axis.set_xticklabels(xticks)
            axis.tick_params(axis='both',labelsize=11)
            axis.set_xlabel('BECCS Scale Factor (-)', fontsize=12)
            axis.set_ylabel('Contribution to AFFEB, 2015-2100 (GtC)', fontsize=12)
            axis.legend(loc='upper left', ncol=2, fontsize=10)
            fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Uptake_'+comp+'_'+TEMP+'.png')
            fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Uptake_'+comp+'_'+TEMP+'.pdf')
            plt.close()

    # Natural Land uptake and Total Land mitigation
    for comp,ylims in zip(comps,([0,200],[0,250])):
        for TEMP in TEMPs:
            fig,axis = plt.subplots(ncols=1,nrows=1,figsize=[8,5])
            axis.fill_between(BECCS_factors,Cland_mit_lims[comp][TEMP][0],Cland_mit_lims[comp][TEMP][1],
                                color=totLULUCcolor,alpha=0.5,label='Total Land Mitigation',zorder=2)
            axis.fill_between(BECCS_factors,NatLand_mit_lims[comp][TEMP][0],NatLand_mit_lims[comp][TEMP][1],
                                color=Landcolor,alpha=0.7,label='Natural Land Uptake',zorder=3)
            axis.fill_between(BECCS_factors,CCS_mit_lims[comp][TEMP][0],CCS_mit_lims[comp][TEMP][1],
                                color=CCScolor,alpha=1.,label='BECCS',zorder=4)
            axis.plot(xlims,[0,0],c='k')
            #axis.spines['bottom'].set_position(('data',0.0))
            axis.set_xscale('log')
            axis.set_ylim(ylims)
            axis.set_xlim(xlims)
            axis.set_xticks(xticks)
            axis.set_xticklabels(xticks)
            axis.tick_params(axis='both',labelsize=11)
            axis.set_xlabel('BECCS Scale Factor (-)', fontsize=12)
            axis.set_ylabel('Mitigation Potential, 2015-2100 (GtC)', fontsize=12)
            axis.legend(loc='upper left', ncol=2, fontsize=10)
            fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Mitigation_'+comp+'_'+TEMP+'.png')
            fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Mitigation_'+comp+'_'+TEMP+'.pdf')
            plt.close()

    # Total AFFEB, Natural Land uptake and Total Land contribution
    for comp,ylims in zip(comps,([-50,500],[0,600])):
        for TEMP in TEMPs:
            fig,axis = plt.subplots(ncols=1,nrows=1,figsize=[8,5])
            axis.fill_between(BECCS_factors,Total_lims[comp][TEMP][0],Total_lims[comp][TEMP][1],
                color=totcolor,alpha=0.5,label='Total AFFEB',zorder=1)
            axis.fill_between(BECCS_factors,Cland_tot_lims[comp][TEMP][0],Cland_tot_lims[comp][TEMP][1],
                color=totLULUCcolor,alpha=0.5,label='Total Land Contribution',zorder=2)
            axis.fill_between(BECCS_factors,NatLand_tot_lims[comp][TEMP][0],NatLand_tot_lims[comp][TEMP][1],
                color=Landcolor,alpha=0.7,label='Natural Land Uptake',zorder=3)
            axis.fill_between(BECCS_factors,CCS_tot_lims[comp][TEMP][0],CCS_tot_lims[comp][TEMP][1],
                color=CCScolor,alpha=1.,label='BECCS',zorder=4)
            axis.set_xscale('log')
            #axis.set_ylim([0,axis.get_ylim()[1]])
            axis.set_ylim(ylims) 
            axis.set_xlim(xlims)
            axis.set_xticks(xticks)
            axis.set_xticklabels(xticks)
            axis.tick_params(axis='both', labelsize=11)
            axis.set_xlabel('BECCS Scale Factor (-)', fontsize=12)
            axis.set_ylabel('Contribution to AFFEB, 2015-2100 (GtC)', fontsize=12)
            axis.legend(loc='upper left', ncol=2, fontsize=10)
            fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Uptake_wTot_'+comp+'_'+TEMP+'.png')
            fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Uptake_wTot_'+comp+'_'+TEMP+'.pdf')
            #plt.show()
            plt.close()

    # Total, Natural Land and Total Land mitigation potentials
    for comp in comps:
        for TEMP in TEMPs:
          fig,axis = plt.subplots(ncols=1,nrows=1,figsize=[8,5])
          axis.fill_between(BECCS_factors,Mitig_lims[comp][TEMP][0],Mitig_lims[comp][TEMP][1],
              color=totcolor,alpha=0.5,label='Total Mitigation Potential',zorder=1)
          axis.fill_between(BECCS_factors,Cland_mit_lims[comp][TEMP][0],Cland_mit_lims[comp][TEMP][1],
              color=totLULUCcolor,alpha=0.5,label='Total Land Mitigation',zorder=2)
          axis.fill_between(BECCS_factors,NatLand_mit_lims[comp][TEMP][0],NatLand_mit_lims[comp][TEMP][1],
              color=Landcolor,alpha=0.7,label='Natural Land Uptake',zorder=3)
          axis.fill_between(BECCS_factors,CCS_mit_lims[comp][TEMP][0],CCS_mit_lims[comp][TEMP][1],
              color=CCScolor,alpha=1.,label='BECCS',zorder=4)
          axis.set_xscale('log')
          axis.set_ylim([0,axis.get_ylim()[1]])
          axis.set_xlim(xlims)
          axis.set_xticks(xticks)
          axis.set_xticklabels(xticks)
          axis.tick_params(axis='both',labelsize=15)
          axis.set_xlabel('BECCS Scale Factor (-)', fontsize=16)
          axis.set_ylabel('Mitigation Potential, 2015-2100 (GtC)', fontsize=16)
          axis.legend(loc='upper left')
          fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Mitigation_wTot_'+comp+'_'+TEMP+'.png')
          fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Mitigation_wTot_'+comp+'_'+TEMP+'.pdf')
          #plt.show()
          plt.close()

#ipdb.set_trace()
CH4tot_min = {}
CH4tot_max = {}
CH4tot_median = {}
CH4mit_min = {}
CH4mit_max = {}
CH4mit_median = {}

for TEMP in TEMPs:
    CH4_tot_all  = []
    CH4_mit_all  = []
    for iQ10 in range(nQ10s):
        Q10 = Q10_exps[iQ10]
        for iO3 in range(nO3s):
            O3 = Ozone_exps[iO3]
            CH4_tot_all.append(  pd.Series(CH4[TEMP][Q10][O3[1]]['Total']) )
            CH4_mit_all.append(  pd.Series(CH4[TEMP][Q10][O3[1]]['Total'])
                               - pd.Series(CTL[TEMP][Q10][O3[1]]['Total']) )
      
    CH4tot_ps = pd.concat( CH4_tot_all )
    CH4tot_min[TEMP] = CH4tot_ps.describe()['25%']
    CH4tot_median[TEMP] = CH4tot_ps.describe()['50%']
    CH4tot_max[TEMP] = CH4tot_ps.describe()['75%']
    CH4mit_ps = pd.concat( CH4_mit_all ) 
    CH4mit_min[TEMP] = CH4mit_ps.describe()['25%']
    CH4mit_median[TEMP] = CH4mit_ps.describe()['50%']
    CH4mit_max[TEMP] = CH4mit_ps.describe()['75%']

if True:
    xticks = [0.5,0.6,0.8,1,1.5,2,3,4,5,6]
    xlims = [min(xticks),max(xticks)]

#ipdb.set_trace()
if PLOT_FIGURE and PROC_STEP == 'Use':

# Figure Block (3)
# All plots in sub-folder: Sensitivity_Plots
# Plots of sensitivity of mitigation potential of everything
# Separate plots for each temperature pathway
# pdf and png versions of each

  # Mitigation potentials
  for TEMP in TEMPs:
      fig,axis = plt.subplots(ncols=1,nrows=1,figsize=[8,5])
      axis.fill_between(BECCS_factors,Mitig_lims['Linear'][TEMP][0],Mitig_lims['Linear'][TEMP][1],
          color=totcolor_lin,alpha=1.0,label='Total Mitigation Potential (Linear)',zorder=1)
      axis.fill_between(BECCS_factors,Mitig_lims['Coupled'][TEMP][0],Mitig_lims['Coupled'][TEMP][1],
          color=totcolor,alpha=0.7,label='Total Mitigation Potential (Coupled)',zorder=2)
      axis.fill_between(BECCS_factors[[0,-1]],
          [CH4mit_min[TEMP],CH4mit_min[TEMP]],[CH4mit_max[TEMP],CH4mit_max[TEMP]],
          color=CH4color,alpha=1.0,label='Methane Mitigation',zorder=3)
      axis.fill_between(BECCS_factors,Cland_mit_lims['LULUC'][TEMP][0],Cland_mit_lims['LULUC'][TEMP][1],
          color=totLULUCcolor,alpha=0.8,label='Total Land Mitigation',zorder=4)
      axis.fill_between(BECCS_factors,NatLand_mit_lims['LULUC'][TEMP][0],NatLand_mit_lims['LULUC'][TEMP][1],
          color=Landcolor,alpha=1.0,label='Natural Land Uptake',zorder=5)
      axis.fill_between(BECCS_factors,CCS_mit_lims['LULUC'][TEMP][0],CCS_mit_lims['LULUC'][TEMP][1],
          color=CCScolor,alpha=0.8,label='BECCS',zorder=6)

      axis.set_ylim([0,axis.get_ylim()[1]])
      axis.set_xscale('log')
      axis.set_xlim(xlims)
      axis.set_xticks(xticks)
      axis.set_xticklabels(xticks)
      
      axis.tick_params(axis='both',labelsize=11)
      axis.set_xlabel('BECCS Scale Factor (-)', fontsize=12)
      axis.set_ylabel('Mitigation Potential, 2015-2100 (GtC)', fontsize=12)
      axis.legend(loc='upper left', ncol=2, fontsize=10)
      fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Mitigation_everything_'+TEMP+'.png')
      fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Mitigation_everything_'+TEMP+'.pdf')
      #plt.show()
      plt.close()

#ipdb.set_trace()
  # Carbon budgets
  for TEMP in TEMPs:
      fig,axis = plt.subplots(ncols=1,nrows=1,figsize=[8,5])
      axis.fill_between(BECCS_factors,Total_lims['Linear'][TEMP][0],Total_lims['Linear'][TEMP][1],
          color=totcolor_lin,alpha=1.0,label='Total AFFEB (Linear)',zorder=1)
      axis.fill_between(BECCS_factors,Total_lims['Coupled'][TEMP][0],Total_lims['Coupled'][TEMP][1],
          color=totcolor,alpha=0.7,label='Total AFFEB (Coupled)',zorder=2)
      axis.fill_between(BECCS_factors[[0,-1]],
          [CH4tot_min[TEMP],CH4tot_min[TEMP]],[CH4tot_max[TEMP],CH4tot_max[TEMP]],
          color=CH4color,alpha=1.0,label='Methane AFFEB',zorder=3)
      axis.fill_between(BECCS_factors,Cland_tot_lims['LULUC'][TEMP][0],Cland_tot_lims['LULUC'][TEMP][1],
          color=totLULUCcolor,alpha=0.8,label='Total Land AFFEB',zorder=4)
      axis.fill_between(BECCS_factors,NatLand_tot_lims['LULUC'][TEMP][0],NatLand_tot_lims['LULUC'][TEMP][1],
          color=Landcolor,alpha=1.0,label='Natural Land Uptake',zorder=5)
      axis.fill_between(BECCS_factors,CCS_tot_lims['LULUC'][TEMP][0],CCS_tot_lims['LULUC'][TEMP][1],
          color=CCScolor,alpha=0.8,label='BECCS',zorder=6)

      #axis.set_ylim(ylims) 
      #axis.set_ylim([0,axis.get_ylim()[1]])
      axis.set_xscale('log')
      axis.set_xlim(xlims)
      axis.set_xticks(xticks)
      axis.set_xticklabels(xticks)
      
      axis.tick_params(axis='both', labelsize=11)
      axis.set_xlabel('BECCS Scale Factor (-)', fontsize=12)
      axis.set_ylabel('Contribution to AFFEB, 2015-2100 (GtC)', fontsize=12)
      axis.legend(loc='upper left', ncol=2, fontsize=10)

      fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Uptake_everything_'+TEMP+'.png')
      fig.savefig(PLOT_DIR+'Sensitivity_Plots/BECCSsensitivity_Uptake_everything_'+TEMP+'.pdf')
      #plt.show()
      plt.close()


os.system('mkdir -p '+PLOT_DIR+'Paper_MultiPlots/')
if PLOT_FIGURE and PROC_STEP == 'Use':

# Figure Block (4)
# Multi plot Figures for BECCS_ScaleFactor
# All plots in sub-folder: Paper_MultiPlots
# Panel (a) BECCS sensitivity; Panel (b) and (c) Ozone sensitivity
# Separate plots for each temperature profile
# eps, pdf and png versions of each

    for iTEMP in range(nTEMPs):
        fig,axes=plt.subplots(ncols=1,nrows=3,figsize=[6,10])
        fig.subplots_adjust(top=0.97,left=0.1,right=0.95,bottom=0.05,hspace=0.2)
        # At the top, BECCS sensitivity x-y plot
        axis=axes[0]
        axis.fill_between(BECCS_factors,Mitig_lims['Linear'][TEMP][0],Mitig_lims['Linear'][TEMP][1],
            color=totcolor_lin,alpha=1.0,label='Total Mitigation Potential (Linear)',zorder=1)
        axis.fill_between(BECCS_factors,Mitig_lims['Coupled'][TEMP][0],Mitig_lims['Coupled'][TEMP][1],
            color=totcolor,alpha=0.7,label='Total Mitigation Potential (Coupled)',zorder=2)
        axis.fill_between(BECCS_factors[[0,-1]],
            [CH4mit_min[TEMP],CH4mit_min[TEMP]],[CH4mit_max[TEMP],CH4mit_max[TEMP]],
            color=CH4color,alpha=1.0,label='Methane Mitigation',zorder=3)
        axis.fill_between(BECCS_factors,Cland_mit_lims['LULUC'][TEMP][0],Cland_mit_lims['LULUC'][TEMP][1],
            color=totLULUCcolor,alpha=0.8,label='Total Land Mitigation',zorder=4)
        axis.fill_between(BECCS_factors,NatLand_mit_lims['LULUC'][TEMP][0],NatLand_mit_lims['LULUC'][TEMP][1],
            color=Landcolor,alpha=1.0,label='Natural Land Uptake',zorder=5)
        axis.fill_between(BECCS_factors,CCS_mit_lims['LULUC'][TEMP][0],CCS_mit_lims['LULUC'][TEMP][1],
            color=CCScolor,alpha=0.8,label='BECCS',zorder=6)
        axis.set_ylim([0,axis.get_ylim()[1]])
        axis.set_xscale('log')
        axis.set_xlim(xlims)
        axis.set_xticks(xticks)
        axis.set_xticklabels(xticks)
        axis.tick_params(axis='both',labelsize=9)
        axis.set_xlabel('BECCS Scale Factor (-)', fontsize=10)
        axis.set_ylabel('Mitigation Potential, 2015-2100 (GtC)', fontsize=10)
        axis.legend(loc='upper left', ncol=2,fontsize=9)

        # Plot sensitivity maps for high-sensitivitiy and low-sensitivity to ozone
        comp='LULUC_opt'
        TEMP=TEMPs[iTEMP]
        BECCS_colours = ['#bf812d','#f6e8c3','#35978f','#01665e',]  # Switch colours at 2 
        for iO3 in range(nO3s):
            O3 = Ozone_exps[iO3]
            axis=axes[iO3+1]
            scen_median_plotdata = np.median(
               [COMPS[comp][TEMP][Q10_exps[iQ10]][O3[1]]['BECCS_ScaleFactor'][gcm] 
               for gcm in GCMs for iQ10 in range(nQ10s)], axis=0)
            scen_median_plotdata =  scen_median_plotdata[land_index]
            GREY_MASK=np.ones_like(scen_median_plotdata)
            GREY_MASK[(scen_median_plotdata<=0)&(land_index.mask==False)] = 0
            GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
            scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata, 
                mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
            print('Mean Productivity (l. 1294): ',TEMP,O3,np.mean(scen_median_plotdata))
            PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
                COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,
                CLEVELS=BECCS_CLEVELS, SET_OVER='#003c30',CBAR_LABEL='BECCS Scale Factor',
                RESOLUTION='c', MAP_TYPE='Mesh',FONTSIZES=[9,9,10,14],
                GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7',AXIS=axis,)

        axes[0].text(0.02,1.03,'(a)', transform=axes[0].transAxes, fontsize=14, fontweight='bold')
        axes[1].text(0.02,1.03,'(b)', transform=axes[1].transAxes, fontsize=14, fontweight='bold')
        axes[2].text(0.02,1.03,'(c)', transform=axes[2].transAxes, fontsize=14, fontweight='bold')

        fig.savefig(PLOT_DIR+'Paper_MultiPlots/BECCSsensitivity_xy_ScaleFactorMaps_'+TEMP+'.png')
        fig.savefig(PLOT_DIR+'Paper_MultiPlots/BECCSsensitivity_xy_ScaleFactorMaps_'+TEMP+'.pdf')
        fig.savefig(PLOT_DIR+'Paper_MultiPlots/BECCSsensitivity_xy_ScaleFactorMaps_'+TEMP+'.eps')
        #plt.show()
        plt.close()


if PLOT_FIGURE and PROC_STEP == 'Use':

# Figure Block (5)
# Multi plot Figures for BECCS_productivity
# All plots in sub-folder: Paper_MultiPlots
# Panel (a) BECCS sensitivity; Panel (b) and (c) Ozone sensitivity
# Separate plots for each temperature profile
# eps, pdf and png versions of each

    for iTEMP in range(nTEMPs):
        fig,axes=plt.subplots(ncols=1,nrows=3,figsize=[6,10])
        fig.subplots_adjust(top=0.97,left=0.1,right=0.95,bottom=0.05,hspace=0.2)
        # At the top, BECCS sensitivity x-y plot
        axis=axes[0]
        axis.fill_between(BECCS_factors,Mitig_lims['Linear'][TEMP][0],Mitig_lims['Linear'][TEMP][1],
                            color=totcolor_lin,alpha=1.0,label='Total Mitigation Potential (Linear)',zorder=1)
        axis.fill_between(BECCS_factors,Mitig_lims['Coupled'][TEMP][0],Mitig_lims['Coupled'][TEMP][1],
                            color=totcolor,alpha=0.7,label='Total Mitigation Potential (Coupled)',zorder=2)
        axis.fill_between(BECCS_factors[[0,-1]],
                        [CH4mit_min[TEMP],CH4mit_min[TEMP]],[CH4mit_max[TEMP],CH4mit_max[TEMP]],
                            color=CH4color,alpha=1.0,label='Methane Mitigation',zorder=3)
        axis.fill_between(BECCS_factors,Cland_mit_lims['LULUC'][TEMP][0],Cland_mit_lims['LULUC'][TEMP][1],
                            color=totLULUCcolor,alpha=0.8,label='Total Land Mitigation',zorder=4)
        axis.fill_between(BECCS_factors,NatLand_mit_lims['LULUC'][TEMP][0],NatLand_mit_lims['LULUC'][TEMP][1],
                            color=Landcolor,alpha=1.0,label='Natural Land Uptake',zorder=5)
        axis.fill_between(BECCS_factors,CCS_mit_lims['LULUC'][TEMP][0],CCS_mit_lims['LULUC'][TEMP][1],
                            color=CCScolor,alpha=0.8,label='BECCS',zorder=6)
        axis.set_ylim([0,axis.get_ylim()[1]])
        axis.set_xscale('log')
        axis.set_xlim(xlims)
        axis.set_xticks(xticks)
        axis.set_xticklabels(xticks)
        axis.tick_params(axis='both',labelsize=9)
        axis.set_xlabel('BECCS Scale Factor (-)', fontsize=10)
        axis.set_ylabel('Mitigation Potential, 2015-2100 (GtC)', fontsize=10)
        axis.legend(loc='upper left', ncol=2,fontsize=9)

        # Plot sensitivity maps for high-sensitivitiy and low-sensitivity to ozone
        comp='LULUC_opt'
        var='BECCS_productivity'
        #iTEMP=0
        TEMP=TEMPs[iTEMP]
        BECCS_colours = ['#bf812d','#f6e8c3','#35978f','#01665e',]  # Switch colours at 2 
        #ipdb.set_trace()
        for iO3 in range(nO3s):
            O3 = Ozone_exps[iO3]
            axis=axes[iO3+1]
            scen_median_DF = pd.concat(
                    [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[iQ10]][O3[1]][var] )
                      for iQ10 in range(nQ10s)  ],axis=1  ).median(axis=1)
            scen_median_plotdata =  scen_median_DF.values[land_index]
            GREY_MASK=np.ones_like(scen_median_plotdata)
            GREY_MASK[(MAXBIOFRAC_2D<0.01)&(land_index.mask==False)] = 0
            GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
            scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata, 
                    mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
            PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
                    COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,
                    CLEVELS=harvest_CLEVELS, SET_OVER='#003c30',CBAR_LABEL=harvest_label,
                    RESOLUTION='c', MAP_TYPE='Mesh', FONTSIZES=[12,12,14,18],
                    GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7', AXIS=axis, )
                    
        axes[0].text(0.02,1.03,'(a)', transform=axes[0].transAxes, fontsize=14, fontweight='bold')
        axes[1].text(0.02,1.03,'(b)', transform=axes[1].transAxes, fontsize=14, fontweight='bold')
        axes[2].text(0.02,1.03,'(c)', transform=axes[2].transAxes, fontsize=14, fontweight='bold')

        fig.savefig(PLOT_DIR+'Paper_MultiPlots/BECCSsensitivity_xy_ProductivityMaps_'+TEMP+'_byO3.png')
        fig.savefig(PLOT_DIR+'Paper_MultiPlots/BECCSsensitivity_xy_ProductivityMaps_'+TEMP+'_byO3.pdf')
        fig.savefig(PLOT_DIR+'Paper_MultiPlots/BECCSsensitivity_xy_ProductivityMaps_'+TEMP+'_byO3.eps')
        fig.savefig(PLOT_DIR+'Paper_MultiPlots/Figure2_'+TEMP+'_byO3.png')
        fig.savefig(PLOT_DIR+'Paper_MultiPlots/Figure2_'+TEMP+'_byO3.pdf')
        fig.savefig(PLOT_DIR+'Paper_MultiPlots/Figure2_'+TEMP+'_byO3.eps')
        plt.close()

#ipdb.set_trace()

if PLOT_FIGURE and PROC_STEP == 'Use':

# Figure Block (6)
# Multi plot Figures for BECCS_productivity
# All plots in top-level folder
# Panel (a) BECCS sensitivity; Panel (b) Map
# Separate plots for each temperature profile
# eps, pdf and png versions of each

    for iTEMP in range(nTEMPs):
        TEMP=TEMPs[iTEMP]
        fig,axes=plt.subplots(ncols=1,nrows=2,figsize=[6,8])
        fig.subplots_adjust(top=0.97,left=0.1,right=0.95,bottom=0.05,hspace=0.2)
        # At the top, BECCS sensitivity x-y plot
        axis=axes[0]
        axis.fill_between(BECCS_factors,Mitig_lims['Linear'][TEMP][0],Mitig_lims['Linear'][TEMP][1],
            color=totcolor_lin,alpha=1.0,label='Total Mitigation Potential (Linear)',zorder=1)
        axis.fill_between(BECCS_factors,Mitig_lims['Coupled'][TEMP][0],Mitig_lims['Coupled'][TEMP][1],
            color=totcolor,alpha=0.7,label='Total Mitigation Potential (Coupled)',zorder=2)
        axis.fill_between(BECCS_factors[[0,-1]],
            [CH4mit_min[TEMP],CH4mit_min[TEMP]],[CH4mit_max[TEMP],CH4mit_max[TEMP]],
            color=CH4color,alpha=1.0,label='Methane Mitigation',zorder=3)
        axis.fill_between(BECCS_factors,Cland_mit_lims['LULUC'][TEMP][0],Cland_mit_lims['LULUC'][TEMP][1],
            color=totLULUCcolor,alpha=0.8,label='Total Land Mitigation',zorder=4)
        axis.fill_between(BECCS_factors,NatLand_mit_lims['LULUC'][TEMP][0],NatLand_mit_lims['LULUC'][TEMP][1],
            color=Landcolor,alpha=1.0,label='Natural Land Uptake',zorder=5)
        axis.fill_between(BECCS_factors,CCS_mit_lims['LULUC'][TEMP][0],CCS_mit_lims['LULUC'][TEMP][1],
            color=CCScolor,alpha=0.8,label='BECCS',zorder=6)
        axis.set_ylim([0,axis.get_ylim()[1]])
        axis.set_xscale('log')
        axis.set_xlim(xlims)
        axis.set_xticks(xticks)
        axis.set_xticklabels(xticks)
        axis.tick_params(axis='both',labelsize=9)
        axis.set_xlabel('BECCS Scale Factor (-)', fontsize=10)
        axis.set_ylabel('Mitigation Potential, 2015-2100 (GtC)', fontsize=10)
        axis.legend(loc='upper left', ncol=2,fontsize=9)

        # Plot sensitivity maps for high-sensitivitiy and low-sensitivity to ozone
        comp='LULUC_opt'
        BECCS_colours = ['#bf812d','#f6e8c3','#35978f','#01665e',]  # Switch colours at 2 
        var  = 'BECCS_productivity' 
        ivar = 0
        axis=axes[1]

        scen_median_DF = pd.concat(
            [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[0]][Ozone_exps[iO3][1]][var] )
                for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  ).median(axis=1)
        scen_median_plotdata =  scen_median_DF.values[land_index]
        GREY_MASK=np.ones_like(scen_median_plotdata)
        GREY_MASK[(MAXBIOFRAC_2D<0.01)&(land_index.mask==False)] = 0
        GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
        scen_median_plotdata = np.ma.masked_array( scen_median_plotdata,mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
        print('Mean Productivity (l. 1441): ',TEMP,np.mean(scen_median_plotdata))

        scen_median_DF = pd.concat(
            [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[iQ10]][Ozone_exps[iO3][1]][var] )
                for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  ).median(axis=1)
        scen_median_plotdata =  scen_median_DF.values[land_index]
        GREY_MASK=np.ones_like(scen_median_plotdata)
        GREY_MASK[(MAXBIOFRAC_2D<0.01)&(land_index.mask==False)] = 0
        GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
        scen_median_plotdata = np.ma.masked_array( scen_median_plotdata,mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
        print('Mean Productivity (l. 1251): ',TEMP,np.mean(scen_median_plotdata))

        PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
                    COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,
                    CLEVELS=harvest_CLEVELS, SET_OVER='#003c30', SET_UNDER='#543005',
                    RESOLUTION='c', MAP_TYPE='Mesh',CBAR_LABEL=harvest_label,
                    GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7',FONTSIZES=[12,12,14,18], AXIS=axis )
                    
        axes[0].text(0.02,1.03,'(a)', transform=axes[0].transAxes, fontsize=14, fontweight='bold')
        axes[1].text(0.02,1.03,'(b)', transform=axes[1].transAxes, fontsize=14, fontweight='bold')

        FILENAME  = 'BECCSsensitivity_xy_ProductivityMaps_'+TEMP
        print('Writing to: '+FILENAME)
        fig.savefig(PLOT_DIR+FILENAME+'.png')
        fig.savefig(PLOT_DIR+FILENAME+'.pdf')
        fig.savefig(PLOT_DIR+FILENAME+'.eps')
        #plt.show()
        plt.close()


#ipdb.set_trace()

# Figure Block (7)
# Multi plot Figures for BECCS_productivity
# All plots in top-level folder
# Panel (a) BECCS sensitivity; Panel (b) BECCS_ScaleFactor; (c) BECCS_productivity
# Separate plots for each temperature profile
# eps, pdf and png versions of each

# Amended code as only using one Q10 in previous versions of code

for iTEMP in range(nTEMPs):
    fig,axes=plt.subplots(ncols=1,nrows=3,figsize=[6,10])
    fig.subplots_adjust(top=0.97,left=0.1,right=0.95,bottom=0.05,hspace=0.2)

    # At the top, BECCS sensitivity x-y plot
    axis=axes[0]
    axis.fill_between(BECCS_factors,Mitig_lims['Linear'][TEMP][0],Mitig_lims['Linear'][TEMP][1],
        color=totcolor_lin,alpha=1.0,label='Total Mitigation Potential (Linear)',zorder=1)
    axis.fill_between(BECCS_factors,Mitig_lims['Coupled'][TEMP][0],Mitig_lims['Coupled'][TEMP][1],
        color=totcolor,alpha=0.7,label='Total Mitigation Potential (Coupled)',zorder=2)
    axis.fill_between(BECCS_factors[[0,-1]],
        [CH4mit_min[TEMP],CH4mit_min[TEMP]],[CH4mit_max[TEMP],CH4mit_max[TEMP]],
        color=CH4color,alpha=1.0,label='Methane Mitigation',zorder=3)
    axis.fill_between(BECCS_factors,Cland_mit_lims['LULUC'][TEMP][0],Cland_mit_lims['LULUC'][TEMP][1],
        color=totLULUCcolor,alpha=0.8,label='Total Land Mitigation',zorder=4)
    axis.fill_between(BECCS_factors,NatLand_mit_lims['LULUC'][TEMP][0],NatLand_mit_lims['LULUC'][TEMP][1],
        color=Landcolor,alpha=1.0,label='Natural Land Uptake',zorder=5)
    axis.fill_between(BECCS_factors,CCS_mit_lims['LULUC'][TEMP][0],CCS_mit_lims['LULUC'][TEMP][1],
                        color=CCScolor,alpha=0.8,label='BECCS',zorder=6)
    axis.set_ylim([0,axis.get_ylim()[1]])
    axis.set_xscale('log')
    axis.set_xlim(xlims)
    axis.set_xticks(xticks)
    axis.set_xticklabels(xticks)
    axis.tick_params(axis='both',labelsize=9)
    axis.set_xlabel('BECCS Scale Factor (-)', fontsize=10)
    axis.set_ylabel('Mitigation Potential, 2015-2100 (GtC)', fontsize=10)
    axis.legend(loc='upper left', ncol=2,fontsize=9)

    # Maps of BECCS_ScaleFactor and BECCS_productivity
    comp='LULUC_opt'
    TEMP=TEMPs[iTEMP]
    BECCS_colours = ['#bf812d','#f6e8c3','#35978f','#01665e',]  # Switch colours at 2 
    plotvarnames = [ 'BECCS_ScaleFactor', 'BECCS_productivity']
    nvars = len(plotvarnames)
    for ivar in range(nvars):
        var = plotvarnames[ivar]
        axis=axes[ivar+1]

        scen_median_DF = pd.concat(
            [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[0]][Ozone_exps[iO3][1]][var] )
                for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  ).median(axis=1)
        scen_median_plotdata =  scen_median_DF.values[land_index]
        GREY_MASK=np.ones_like(scen_median_plotdata)
        GREY_MASK[(MAXBIOFRAC_2D<0.01)&(land_index.mask==False)] = 0
        GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
        scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata,mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
        print('Mean Productivity (l. 1529): ',var,TEMP,np.mean(scen_median_plotdata))

        scen_median_DF = pd.concat(
            [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[iQ10]][Ozone_exps[iO3][1]][var] )
                for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  ).median(axis=1)
        scen_median_plotdata =  scen_median_DF.values[land_index]
        GREY_MASK=np.ones_like(scen_median_plotdata)
        GREY_MASK[(MAXBIOFRAC_2D<0.01)&(land_index.mask==False)] = 0
        GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
        scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata,mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
        print('Mean Productivity (l. 1539): ',var,TEMP,np.mean(scen_median_plotdata))

        PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
            COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,
            CLEVELS=harvest_CLEVELS, SET_OVER='#003c30', SET_UNDER='#543005',
            RESOLUTION='c', MAP_TYPE='Mesh',
            CBAR_LABEL=harvest_label,
            FONTSIZES=[12,12,14,18],
            GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7',
            AXIS=axis  )

    axes[0].text(0.02,1.03,'(a)', transform=axes[0].transAxes, fontsize=14, fontweight='bold')
    axes[1].text(0.02,1.03,'(b)', transform=axes[1].transAxes, fontsize=14, fontweight='bold')
    axes[2].text(0.02,1.03,'(c)', transform=axes[2].transAxes, fontsize=14, fontweight='bold')

    FILENAME  = 'BECCSsensitivity_xy_ScaleFactor_Productivity_maps_'+TEMP
    print('Writing to: '+FILENAME)
    fig.savefig(PLOT_DIR+FILENAME+'.png')
    fig.savefig(PLOT_DIR+FILENAME+'.pdf')
    fig.savefig(PLOT_DIR+FILENAME+'.eps')
    #plt.show()
    plt.close()

#ipdb.set_trace()

# Figure Block (8)
# Multi plot Figures for BECCS_productivity
# Used in original NCC submission - NOTE: subsequent change to panel (b) label
# All plots in top-level folder
# Panel (a) BECCS sensitivity; Panel (b) BECCS_productivity_JULES; (c) BECCS_productivity
# Separate plots for each temperature profile
# eps, pdf and png versions of each

# Amended code as only using one Q10 in previous versions of code

for iTEMP in range(nTEMPs):
    fig,axes=plt.subplots(ncols=1,nrows=3,figsize=[6,10])
    fig.subplots_adjust(top=0.97,left=0.1,right=0.95,bottom=0.05,hspace=0.2)
    # At the top, BECCS sensitivity x-y plot
    axis=axes[0]
    axis.fill_between(BECCS_factors,Mitig_lims['Linear'][TEMP][0],Mitig_lims['Linear'][TEMP][1],
        color=totcolor_lin,alpha=1.0,label='Total Mitigation Potential (Linear)',zorder=1)
    axis.fill_between(BECCS_factors,Mitig_lims['Coupled'][TEMP][0],Mitig_lims['Coupled'][TEMP][1],
        color=totcolor,alpha=0.7,label='Total Mitigation Potential (Coupled)',zorder=2)
    axis.fill_between(BECCS_factors[[0,-1]],
        [CH4mit_min[TEMP],CH4mit_min[TEMP]],[CH4mit_max[TEMP],CH4mit_max[TEMP]],
        color=CH4color,alpha=1.0,label='Methane Mitigation',zorder=3)
    axis.fill_between(BECCS_factors,Cland_mit_lims['LULUC'][TEMP][0],Cland_mit_lims['LULUC'][TEMP][1],
        color=totLULUCcolor,alpha=0.8,label='Total Land Mitigation',zorder=4)
    axis.fill_between(BECCS_factors,NatLand_mit_lims['LULUC'][TEMP][0],NatLand_mit_lims['LULUC'][TEMP][1],
        color=Landcolor,alpha=1.0,label='Natural Land Uptake',zorder=5)
    axis.fill_between(BECCS_factors,CCS_mit_lims['LULUC'][TEMP][0],CCS_mit_lims['LULUC'][TEMP][1],
        color=CCScolor,alpha=0.8,label='BECCS',zorder=6)
    axis.set_ylim([0,axis.get_ylim()[1]])
    axis.set_xscale('log')
    axis.set_xlim(xlims)
    axis.set_xticks(xticks)
    axis.set_xticklabels(xticks)
    axis.tick_params(axis='both',labelsize=9)
    axis.set_xlabel('BECCS Scale Factor (-)', fontsize=10)
    axis.set_ylabel('Mitigation Potential, 2015-2100 (GtC)', fontsize=10)
    axis.legend(loc='upper left', ncol=2,fontsize=9)

    # Plot of BEECS productivity
    comp='LULUC_opt'
    TEMP=TEMPs[iTEMP]
    BECCS_colours = ['#bf812d','#f6e8c3','#35978f','#01665e',]  # Switch colours at 2 
    plotvarnames = [ 'BECCS_productivity_JULES', 'BECCS_productivity']
    nvars = len(plotvarnames)
    for ivar in range(nvars):
        var = plotvarnames[ivar]
        if var == 'BECCS_productivity_JULES':
            cbar_label = 'Modelled bioenergy crop yield (tonnesDM ha$^{-1}$ yr$^{-1}$)'
        else:
            cbar_label = harvest_label 
        axis=axes[ivar+1]

        scen_median_DF = pd.concat(
            [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[0]][Ozone_exps[iO3][1]][var] )
                for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  ).median(axis=1)
        scen_median_plotdata =  scen_median_DF.values[land_index]
        GREY_MASK=np.ones_like(scen_median_plotdata)
        GREY_MASK[(MAXBIOFRAC_2D<0.01)&(land_index.mask==False)] = 0
        GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
        scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata,mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
        print('Mean Productivity (l. 1624): ',var,TEMP,np.mean(scen_median_plotdata), \
            np.min(scen_median_plotdata),np.max(scen_median_plotdata))
        if DEBUG == 'Y': print(scen_median_plotdata[scen_median_plotdata.mask == False])

        scen_median_DF = pd.concat(
            [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[iQ10]][Ozone_exps[iO3][1]][var] )
                for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  ).median(axis=1)
        scen_median_plotdata =  scen_median_DF.values[land_index]
        GREY_MASK=np.ones_like(scen_median_plotdata)
        GREY_MASK[(MAXBIOFRAC_2D<0.01)&(land_index.mask==False)] = 0
        GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
        scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata,mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
        print('Mean Productivity (l. 1636): ',var,TEMP,np.mean(scen_median_plotdata), \
            np.min(scen_median_plotdata),np.max(scen_median_plotdata))
        if DEBUG == 'Y': print(scen_median_plotdata[scen_median_plotdata.mask == False])

        PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
            COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,
            CLEVELS=harvest_CLEVELS, SET_OVER='#003c30', SET_UNDER='#543005',
            RESOLUTION='c', MAP_TYPE='Mesh',
            CBAR_LABEL=cbar_label,
            FONTSIZES=[12,12,14,18],
            GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7',
            AXIS=axis  )

    axes[0].text(0.02,1.03,'(a)', transform=axes[0].transAxes, fontsize=14, fontweight='bold')
    axes[1].text(0.02,1.03,'(b)', transform=axes[1].transAxes, fontsize=14, fontweight='bold')
    axes[2].text(0.02,1.03,'(c)', transform=axes[2].transAxes, fontsize=14, fontweight='bold')

    FILENAME  = 'BECCSsensitivity_xy_Productivity_ProdJULES_maps_'+TEMP+'_paper'
    print('Writing to: '+FILENAME)
    fig.savefig(PLOT_DIR+FILENAME+'.png')
    fig.savefig(PLOT_DIR+FILENAME+'.pdf')
    fig.savefig(PLOT_DIR+FILENAME+'.eps')
    #plt.show()
    plt.close()

#ipdb.set_trace()

# Figure Block (9)
# Multi plot Figures for BECCS_productivity based on NPP
# Used in revised NCC submission
# All plots in top-level folder
# Panel (a) BECCS sensitivity; Panel (b) BECCS_productivity_JULES; (c) BECCS_productivity
# Separate plots for each temperature profile
# eps, pdf and png versions of each

# Amended code as only using one Q10 in previous versions of code

if BECCS_NPP == 'Y' and PROC_STEP == 'Use':

    for iTEMP in range(nTEMPs):
        fig,axes=plt.subplots(ncols=1,nrows=3,figsize=[6,10])
        fig.subplots_adjust(top=0.97,left=0.1,right=0.95,bottom=0.05,hspace=0.2)
        # At the top, BECCS sensitivity x-y plot
        axis=axes[0]
        axis.fill_between(BECCS_factors,Mitig_lims['Linear'][TEMP][0],Mitig_lims['Linear'][TEMP][1],
            color=totcolor_lin,alpha=1.0,label='Total Mitigation Potential (Linear)',zorder=1)
        axis.fill_between(BECCS_factors,Mitig_lims['Coupled'][TEMP][0],Mitig_lims['Coupled'][TEMP][1],
            color=totcolor,alpha=0.7,label='Total Mitigation Potential (Coupled)',zorder=2)
        axis.fill_between(BECCS_factors[[0,-1]],
            [CH4mit_min[TEMP],CH4mit_min[TEMP]],[CH4mit_max[TEMP],CH4mit_max[TEMP]],
            color=CH4color,alpha=1.0,label='Methane Mitigation',zorder=3)
        axis.fill_between(BECCS_factors,Cland_mit_lims['LULUC'][TEMP][0],Cland_mit_lims['LULUC'][TEMP][1],
            color=totLULUCcolor,alpha=0.8,label='Total Land Mitigation',zorder=4)
        axis.fill_between(BECCS_factors,NatLand_mit_lims['LULUC'][TEMP][0],NatLand_mit_lims['LULUC'][TEMP][1],
            color=Landcolor,alpha=1.0,label='Natural Land Uptake',zorder=5)
        axis.fill_between(BECCS_factors,CCS_mit_lims['LULUC'][TEMP][0],CCS_mit_lims['LULUC'][TEMP][1],
            color=CCScolor,alpha=0.8,label='BECCS',zorder=6)
        axis.set_ylim([0,axis.get_ylim()[1]])
        axis.set_xscale('log')
        axis.set_xlim(xlims)
        axis.set_xticks(xticks)
        axis.set_xticklabels(xticks)
        axis.tick_params(axis='both',labelsize=9)
        axis.set_xlabel('BECCS Scale Factor (-)', fontsize=10)
        axis.set_ylabel('Mitigation Potential, 2015-2100 (GtC)', fontsize=10)
        axis.legend(loc='upper left', ncol=2,fontsize=9)

        # BECCS_productivity_JULES and BECCS_productivity_NPP
        comp='LULUC_opt'
        TEMP=TEMPs[iTEMP]
        BECCS_colours = ['#bf812d','#f6e8c3','#35978f','#01665e',]  # Switch colours at 2 
        plotvarnames = [ 'BECCS_productivity_JULES', 'BECCS_productivity_NPP']
        nvars = len(plotvarnames)
        for ivar in range(nvars):
            var = plotvarnames[ivar]
            if var == 'BECCS_productivity_JULES':
                cbar_label = 'Modelled bioenergy crop yield, diagnostic(tonnesDM ha$^{-1}$ yr$^{-1}$)'
            else:
                cbar_label = 'Modelled bioenergy crop yield, NPP (tonnesDM ha$^{-1}$ yr$^{-1}$)'
            axis=axes[ivar+1]

            scen_median_DF = pd.concat(
                [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[0]][Ozone_exps[iO3][1]][var] )
                    for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  ).median(axis=1)
            scen_median_plotdata =  scen_median_DF.values[land_index]
            GREY_MASK=np.ones_like(scen_median_plotdata)
            GREY_MASK[(MAXBIOFRAC_2D<0.01)&(land_index.mask==False)] = 0
            GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
            scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata,mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
            print('Mean Productivity (l. 1725): ',var,TEMP,np.mean(scen_median_plotdata), \
                np.min(scen_median_plotdata),np.max(scen_median_plotdata))
            if DEBUG == 'Y': print(scen_median_plotdata[scen_median_plotdata.mask == False])

            scen_median_DF = pd.concat(
                [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[iQ10]][Ozone_exps[iO3][1]][var] )
                    for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  ).median(axis=1)
            scen_median_plotdata =  scen_median_DF.values[land_index]
            GREY_MASK=np.ones_like(scen_median_plotdata)
            GREY_MASK[(MAXBIOFRAC_2D<0.01)&(land_index.mask==False)] = 0
            GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
            scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata,mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
            print('Mean Productivity (l. 1737): ',var,TEMP,np.mean(scen_median_plotdata), \
                np.min(scen_median_plotdata),np.max(scen_median_plotdata))
            if DEBUG == 'Y': print(scen_median_plotdata[scen_median_plotdata.mask == False])

            PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
                COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,
                CLEVELS=harvest_CLEVELS, SET_OVER='#003c30', SET_UNDER='#543005',
                RESOLUTION='c', MAP_TYPE='Mesh',
                CBAR_LABEL=cbar_label,
                FONTSIZES=[12,12,14,18],
                GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7',
                AXIS=axis  )

        axes[0].text(0.02,1.03,'(a)', transform=axes[0].transAxes, fontsize=14, fontweight='bold')
        axes[1].text(0.02,1.03,'(b)', transform=axes[1].transAxes, fontsize=14, fontweight='bold')
        axes[2].text(0.02,1.03,'(c)', transform=axes[2].transAxes, fontsize=14, fontweight='bold')

        FILENAME  = 'BECCSsensitivity_xy_Productivity_ProdJULES_NPP_maps_'+TEMP+'_paper'
        print('Writing to: '+FILENAME)
        fig.savefig(PLOT_DIR+FILENAME+'.png')
        fig.savefig(PLOT_DIR+FILENAME+'.pdf')
        fig.savefig(PLOT_DIR+FILENAME+'.eps')
        #plt.show()
        plt.close()

#ipdb.set_trace()

# Figure Block (10)
# Figure 10 in ESD paper  
# Multi plot Figures for BECCS_productivity based on NPP
# All plots in top-level folder
# Panel (a) BECCS sensitivity; Panel (b) BECCS_productivity_JULES; (c) BECCS_productivity
# Separate plots for each temperature profile
# eps, pdf and png versions of each

# Amended code as only using one Q10 in previous versions of code

for iTEMP in range(nTEMPs):
    # Plot figure with subplots of different sizes
    fig,axes=plt.subplots(ncols=2,nrows=2,figsize=[12,6])
    fig.subplots_adjust(top=0.97,left=0.1,right=0.95,bottom=0.05,hspace=0.2)
    # set up subplot grid

    # At the top, BECCS sensitivity x-y plot
    axis=axes[0,0]
    axis.fill_between(BECCS_factors,Mitig_lims['Linear'][TEMP][0],Mitig_lims['Linear'][TEMP][1],
        color=totcolor_lin,alpha=1.0,label='Total Mitigation Potential (Linear)',zorder=1)
    axis.fill_between(BECCS_factors,Mitig_lims['Coupled'][TEMP][0],Mitig_lims['Coupled'][TEMP][1],
        color=totcolor,alpha=0.7,label='Total Mitigation Potential (Coupled)',zorder=2)
    axis.fill_between(BECCS_factors[[0,-1]],
        [CH4mit_min[TEMP],CH4mit_min[TEMP]],[CH4mit_max[TEMP],CH4mit_max[TEMP]],
        color=CH4color,alpha=1.0,label='Methane Mitigation',zorder=3)
    axis.fill_between(BECCS_factors,Cland_mit_lims['LULUC'][TEMP][0],Cland_mit_lims['LULUC'][TEMP][1],
        color=totLULUCcolor,alpha=0.8,label='Total Land Mitigation',zorder=4)
    axis.fill_between(BECCS_factors,NatLand_mit_lims['LULUC'][TEMP][0],NatLand_mit_lims['LULUC'][TEMP][1],
        color=Landcolor,alpha=1.0,label='Natural Land Uptake',zorder=5)
    axis.fill_between(BECCS_factors,CCS_mit_lims['LULUC'][TEMP][0],CCS_mit_lims['LULUC'][TEMP][1],
        color=CCScolor,alpha=0.8,label='BECCS',zorder=6)
    axis.set_ylim([0,axis.get_ylim()[1]])
    axis.set_xscale('log')
    axis.set_xlim(xlims)
    axis.set_xticks(xticks)
    axis.set_xticklabels(xticks)
    axis.tick_params(axis='both',labelsize=9)
    axis.set_xlabel('BECCS Scale Factor (-)', fontsize=10)
    axis.set_ylabel('Mitigation Potential, 2015-2100 (GtC)', fontsize=10)
    axis.legend(loc='upper left', ncol=2, fontsize=9)
    axis.text(0.02,1.03,'('+ALPHABET[0]+')', transform=axis.transAxes, fontsize=14, fontweight='bold')

    # Maps of BECCS_productivity_JULES, BECCS_productivity, NPP and Soil Carbon
    comp='LULUC_opt'
    TEMP=TEMPs[iTEMP]
    BECCS_colours = ['#bf812d','#f6e8c3','#35978f','#01665e',]  # Switch colours at 2 
    plotvarnames  = [ 'Soil_Carbon_gb', 'BECCS_productivity_JULES', 'BECCS_productivity' ]
    plotvarlabels = [ \
                      'Soil Carbon (kg-C m$^{-2}$)', \
                      'Modelled bioenergy crop yield (tonnesDM ha$^{-1}$ yr$^{-1}$)', \
                       harvest_label, \
                    ]
    nvars = len(plotvarnames)
    for ivar in range(nvars):
        var = plotvarnames[ivar]
        cbar_label = plotvarlabels[ivar]

        if var == 'BECCS_productivity_JULES':
            axis         = axes[1,0]
            plot_CLEVELS = harvest_CLEVELS 
        elif var == 'BECCS_productivity':
            axis         = axes[1,1]
            plot_CLEVELS = harvest_CLEVELS
        elif var == 'Soil_Carbon_gb':
            axis         = axes[0,1]
            # plot_CLEVELS = [ -20.0,-10.0,-5.0,-2.0,0.0,5.0,10.0,20.0 ] 
            plot_CLEVELS = [ -5.0,-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0,5.0 ]

        print('l. 1832: ',ivar,var)
        scen_median_DF = pd.concat(
            [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[iQ10]][Ozone_exps[iO3][1]][var] )
              for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  )
        # print(scen_median_DF)
        scen_median_DF = scen_median_DF.median(axis=1)
        scen_median_plotdata =  scen_median_DF.values[land_index]
        GREY_MASK=np.ones_like(scen_median_plotdata)
        GREY_MASK[(MAXBIOFRAC_2D<0.01)&(land_index.mask==False)] = 0
        GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
        scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata,mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
        print('Mean Productivity (l. 1843): ',var,TEMP,np.mean(scen_median_plotdata), \
            np.min(scen_median_plotdata),np.max(scen_median_plotdata))
        if DEBUG == 'Y': print(scen_median_plotdata[scen_median_plotdata.mask == False])

        PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
            COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,
            CLEVELS=plot_CLEVELS, SET_OVER='#003c30', SET_UNDER='#543005',
            RESOLUTION='c', MAP_TYPE='Mesh',
            CBAR_LABEL=cbar_label,
            FONTSIZES=[12,12,14,18],
            GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7',
            AXIS=axis )
        TEXT_label = '('+ALPHABET[ivar+1]+')'
        axis.text(0.02,1.03,TEXT_label, transform=axis.transAxes, fontsize=14, fontweight='bold')

    FILENAME  = 'BECCSsensitivity_xy_Productivity_ProdJULES_maps_Carbon_'+TEMP+'_paper'
    print('Writing to: '+FILENAME)
    fig.savefig(PLOT_DIR+FILENAME+'.png')
    fig.savefig(PLOT_DIR+FILENAME+'.pdf')
    fig.savefig(PLOT_DIR+FILENAME+'.eps')
    #plt.show()
    plt.close()

# Figure Block (11)
# Multi plot Figures for BECCS_productivity based on NPP
# All plots in top-level folder
# Panel (a) BECCS sensitivity; Panel (b) BECCS_productivity_JULES; (c) BECCS_productivity
# Separate plots for each temperature profile
# eps, pdf and png versions of each

# Amended code as only using one Q10 in previous versions of code

COMPS_list    = [ 'LULUC_CCS', 'LULUC_Nat', 'LULUC_opt']
COMPS_labels  = [ 'BECCS', 'Natural', 'Difference' ]
BECCS_colours = ['#bf812d','#f6e8c3','#35978f','#01665e',]  # Switch colours at 2 
plotvarnames  = [ 'Soil_Carbon_gb', 'Soil_Carbon', 'NPP_gb' ]
plotvarlabels = [ 'Total Soil Carbon', 'Soil Carbon', 'Net Primary Productivity' ]
plotvarunits  = ' (kg-C m$^{-2}$)'
nvars         = len(plotvarnames)

for iTEMP in range(nTEMPs):

    # Plot figure with subplots of different sizes
    fig,axes=plt.subplots(ncols=3,nrows=3,figsize=[18,9])
    fig.subplots_adjust(top=0.97,left=0.1,right=0.95,bottom=0.05,hspace=0.2)
    # set up subplot grid

    # Maps of NPP and Soil Carbon
    TEMP=TEMPs[iTEMP]
    iaxis = 0

    for icomp,comp in enumerate(COMPS_list):
        for ivar in range(nvars):
           var = plotvarnames[ivar]
           cbar_label = plotvarlabels[ivar]+plotvarunits

           if var == 'Soil_Carbon_gb':
               # plot_CLEVELS = [ -20.0,-10.0,-5.0,-2.0,0.0,5.0,10.0,20.0 ]
               plot_CLEVELS = [ -5.0,-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0,5.0 ] 
           else:
               plot_CLEVELS = [ -2.0,-1.0,-0.5,-0.2,0.0,0.2,0.5,1.0,2.0 ]

           axis=axes[icomp, ivar]
           print('l. 1906: ',icomp,ivar,comp,var)
           scen_median_DF = pd.concat(
               [ pd.DataFrame( COMPS[comp][TEMP][Q10_exps[iQ10]][Ozone_exps[iO3][1]][var] )
                 for iO3 in range(nO3s) for iQ10 in range(nQ10s) ],axis=1  )
           # print(scen_median_DF)
           scen_median_DF = scen_median_DF.median(axis=1)
           scen_median_plotdata =  scen_median_DF.values[land_index]
           GREY_MASK=np.ones_like(scen_median_plotdata)
           GREY_MASK[(MAXBIOFRAC_2D<0.01)&(land_index.mask==False)] = 0
           GREY_MASK = np.ma.masked_equal(GREY_MASK,1)
           scen_median_plotdata =  np.ma.masked_array( scen_median_plotdata,mask=(land_index.mask)|(MAXBIOFRAC_2D<0.01) )
           print('Mean Productivity (l. 1918): ',var,TEMP,np.mean(scen_median_plotdata), \
               np.min(scen_median_plotdata),np.max(scen_median_plotdata))
           if DEBUG == 'Y': print(scen_median_plotdata[scen_median_plotdata.mask == False])

           PTs.plot_map(scen_median_plotdata,lons_2d,lats_2d,
               COLOURS=BECCS_colours,INTERPOLATE_COLOURS=True,
               CLEVELS=plot_CLEVELS, SET_OVER='#003c30', SET_UNDER='#543005',
               RESOLUTION='c', MAP_TYPE='Mesh',
               CBAR_LABEL=cbar_label,
               FONTSIZES=[12,12,14,18],
               GREYMASK=GREY_MASK,MASKCOLOR='#d7d7d7',
               AXIS=axis )
           TEXT = '('+ALPHABET[iaxis]+')'+COMPS_labels[icomp]
           axis.text(0.02, 1.03, TEXT, transform=axis.transAxes, fontsize=14, fontweight='bold')

           iaxis += 1

    FILENAME  = 'BECCSsensitivity_maps_Carbon_'+TEMP+'_paper'
    print('Writing to: '+FILENAME)
    fig.savefig(PLOT_DIR+FILENAME+'.png')
    fig.savefig(PLOT_DIR+FILENAME+'.pdf')
    fig.savefig(PLOT_DIR+FILENAME+'.eps')
    #plt.show()
    plt.close()

