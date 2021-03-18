#!/bin/env python2.7

# Uncomment for the UK JASMIN platform
#import matplotlib as mpl
#mpl.use('Agg')
#mpl.style.use('classic')

import numpy as np
import netCDF4 as nc
import sys,os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from copy import deepcopy
import datetime as dt

from imogen import data_info
import plot_functions

#import ipdb

def optional_argparse(arg,default):
    if arg in sys.argv:
        temp_loc=sys.argv.index(arg)
        temp_arg=sys.argv.pop(temp_loc)
        value=sys.argv.pop(temp_loc)
    else:
        value=default
    return value

def get_ticks(MIN,MAX,INC):
    nTICKS         = int((MAX-MIN)/INC)+1
    TICKS          = MIN+INC*np.arange(nTICKS)
    return TICKS

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
PLOT_TAG      = optional_argparse('-plottag', 'Profiles')
PROFILE_npy   = optional_argparse('-PROF_NPY','Use')
CCS_minimum_threshold     = float(optional_argparse('CCS_min','1e-4'))
BIOFRAC_minimum_threshold = float(optional_argparse('CCSfrac_min','1e-2'))
PLOT_OPT      = str(optional_argparse('-plot_opt','1'))

# Directories containing JULES output and plot output directories:
if PLATFORM == 'JASMIN':
    HOME_DIR      = '/gws/nopw/j04/clifftop/SYNTHESIS/PostReview_output/'
    DATA_DIR      = HOME_DIR
    ANCILS_DIR    = '/gws/nopw/j04/clifftop/COMMON_DATA/ANCILS/'
    COMMON_DIR    = '/gws/nopw/j04/clifftop/COMMON_DATA/SCENARIOS/'
    GHG_DIR       = COMMON_DIR+'NCC_synthesis_files/'
    PLOT_DIR      = optional_argparse('-plotdir',DATA_DIR+'plots/'+PLOT_TAG+'/')
    PROF_npy_DIR  = '/group_workspaces/jasmin2/jules/ghayman/PYTHON/SYNTHESIS/npy_files/profiles/'

elif PLATFORM == 'CEH':
    HOME_DIR      = '/prj/CLIFFTOP/SYNTHESIS/'
    PLOT_DIR      = '/data/grp/eow/garr/Projects/CLIFFTOP/Paper_Synthesis/PLOTS_'+PLOT_TAG+'/'
    DATA_DIR      = HOME_DIR+'Review_Response_Check/GCM_Output/'
    COMMON_DIR    = HOME_DIR+'IMOGEN_data/'
    GHG_DIR       = COMMON_DIR+'z_RUNS_Revised/'
    ANCILS_DIR    = HOME_DIR+'Land_Cover/'
    PROF_npy_DIR  = HOME_DIR+'Review_Response_Check/profiles/'

if PLATFORM == 'CEH' and PROFILE_npy == 'Save':
    # For testing and development
    Q10_exps      = [ 'highQ10' ]
    Ozone_exps    = [['L','lowO3']]

    COMPS_opt     = [ 'LULUC_opt' ]
    COMPS_keys    = [ 'CTL', 'LULUC_CCS', 'LULUC_Nat' ]
    COMPS         = {
                      'CTL': { 'config': 'highCH4_OZONESUB_LULUCBL', 'runid':'H_OZONESUB_BL' }
                    , 'LULUC_CCS': { 'config':'highCH4_OZONESUB_LULUC1.9', 'runid':'H_OZONESUB_19' }
                    , 'LULUC_Nat': { 'config':'highCH4_OZONESUB_LULUC1.9Nat', 'runid':'H_OZONESUB_19N' }
                    }

else:
    # Full ensemble
    Q10_exps      = [ 'lowQ10', 'highQ10' ]
    Ozone_exps    = [['L','lowO3'], ['H','highO3'], ]

    COMPS_opt     = [ 'LULUC_opt','Coupled_opt' ]
    COMPS_keys    = [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]
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
PLOT_DIR      = optional_argparse('-plotdir',  PLOT_DIR)
print('DATA_DIR: '+DATA_DIR)
print('PLOT_DIR: '+PLOT_DIR)
os.system('mkdir -p '+PLOT_DIR )
os.system('mkdir -p '+PLOT_DIR+'npy_files' )

ALPHABET      = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o']

START_YEAR    = 2000
END_YEAR      = 2100
nYEARS        = END_YEAR-START_YEAR # No data for 2100

# Scenarios to plot:
TEMPs         = ['1p5deg', '2deg' ] #'1p81p5deg',   # tag in the JULES output file directory
TEMP_years    = [2099,2099,2099]    # tag in the JULES output file directory 
TEMP_names    = ['1.5$^o$C (2100)','2.0$^o$C (2100)'] # '1.5$^o$C Overshoot (2100)',
                # Name to appear on plots etc.

pools         =  ['CH4','CO2']
pools_gen     =  ['CH4','CH4_YEAR','CO2','CO2_YEAR','T','T_YEAR','Q','Q_YEAR']

Tile_names    = data_info.TILE_short_names()
Tile_colours  = data_info.TILE_colours()
nTiles        = len(Tile_names)
nTEMPs        = len(TEMPs)
nQ10s         = len(Q10_exps)
nO3s          = len(Ozone_exps)
npools        = len(pools)
npools_gen    = len(pools_gen)
nPFTs         = 13

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

for comp in COMPS_keys:
    CH4 = COMPS[comp]['config'].split('_')[0]

    for iTEMP in range(nTEMPs): 
        TEMP=TEMPs[iTEMP]
        TEMP_year=TEMP_years[iTEMP]
        TEMP_tag=TEMP
        COMPS[comp][TEMP]={}

        for igcm in range(nGCMs):
            gcm=GCMs[igcm]
            gcm_index=GCM_index[igcm]
            COMPS[comp][TEMP][gcm]={}

            # Time-series of atmospheric concentrations
            if PROFILE_npy != 'Use':

                for iQ10 in range(nQ10s):
                    Q10 = Q10_exps[iQ10]
                    COMPS[comp][TEMP][gcm][Q10]={}

                    for iO3 in range(nO3s):
                        O3 = Ozone_exps[iO3]
                        COMPS[comp][TEMP][gcm][Q10][O3[1]]={ pool:{} for pool in pools  }  #+land_opt_pools }
                        config=COMPS[comp]['config'].replace('OZONESUB',O3[1])
                        runid=COMPS[comp]['runid'].replace('OZONESUB',O3[0])

                        print; print('GCM Input:',comp,gcm,TEMP,TEMP_year,Q10,config,runid)

                        first_year = True

                        for iyear in range(nYEARS):
                             YEAR = START_YEAR+iyear
                             DUMP_FILE=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+TEMP_tag+'.dump.'+str(YEAR)+'0101.0.nc'
                             if iyear == 0 or iyear == nYEARS-1: print(YEAR,DUMP_FILE)
                             # Open current year files:
                             Dinf = nc.Dataset(DUMP_FILE,'r')

                             # Read in data from current year
                             ATM_CH4    = Dinf.variables['ch4_ppbv'][:].squeeze()
                             ATM_CO2    = Dinf.variables['co2_ppmv'][:].squeeze()

                             if first_year:
                                 first_year   = False
                                 ATM_CH4_ALL  = np.zeros((nYEARS))
                                 ATM_CO2_ALL  = np.zeros((nYEARS))

                             if ATM_CH4.size == 1:
                                 ATM_CH4_ALL[iyear] = ATM_CH4
                             else:
                                 print('No data for CH4 for '+str(YEAR))
                                 ATM_CH4_ALL[iyear] = float('nan')

                             if ATM_CO2.size == 1:
                                 ATM_CO2_ALL[iyear]  = ATM_CO2
                             else:
                                 print('No data for CO2 for '+str(YEAR))
                                 ATM_CO2_ALL[iyear]  = float('nan')

                             Dinf.close()

                             if iyear == 0 or iyear == nYEARS-1:
                                 print(YEAR,ATM_CH4_ALL.min(),ATM_CH4_ALL.max(), \
                                            ATM_CO2_ALL.min(),ATM_CO2_ALL.max() )

                        COMPS[comp][TEMP][gcm][Q10][O3[1]]['CH4'] = ATM_CH4_ALL[:]
                        COMPS[comp][TEMP][gcm][Q10][O3[1]]['CO2'] = ATM_CO2_ALL[:]

                        del ATM_CH4_ALL,ATM_CO2_ALL

            if PROFILE_npy != 'Use':

                # Save as pickle file
                PROFILE_FILE = PROF_npy_DIR+comp+'_'+gcm+'_'+TEMP_tag+'.pkl'
                print('Writing to: '+PROFILE_FILE)
                PROF_ID      = open(PROFILE_FILE,'wb')
                pickle.dump(COMPS[comp][TEMP][gcm],PROF_ID)
                PROF_ID.close()

            elif PROFILE_npy == 'Use':

                # Read from pickle file
                PROFILE_FILE = PROF_npy_DIR+comp+'_'+gcm+'_'+TEMP_tag+'.pkl'
                print('Reading from: '+PROFILE_FILE)
                PROF_ID      = open(PROFILE_FILE,'rb')
                COMPS[comp][TEMP][gcm] = pickle.load(PROF_ID)
                PROF_ID.close()

# No plotting if saving to NPP numpy files
if PROFILE_npy != 'Use':
    print('Stopping without plots') 
    quit()

DATA_PROF_COMMON = { pool:{} for pool in pools_gen }
# Get temperature profiles
for iTEMP in range(nTEMPs): 
    TEMP=TEMPs[iTEMP]
    FILE_COMMON = COMMON_DIR+TEMP+'_global_temp_anomaly.dat'
    DATA_INPUT  = pd.read_csv(FILE_COMMON, header=None, names=['Year','Temp'], delim_whitespace=True)
    DATA_PROF_COMMON['T_YEAR'][TEMP] = DATA_INPUT['Year'].values
    DATA_PROF_COMMON['T'][TEMP]      = DATA_INPUT['Temp'].values

# Get historic atmospheric concentrations: CO2
# Data in scenarios common to 2005
FILE_CO2    = GHG_DIR+'SSP2-baseline_IMAGE_concs_co2_vn3p0.txt'
DATA_INPUT  = pd.read_csv(FILE_CO2, header=None, names=['Year','CO2'], delim_whitespace=True)
DATA_PROF_COMMON['CO2_YEAR'] = DATA_INPUT['Year'][0:151].values 
DATA_PROF_COMMON['CO2']      = DATA_INPUT['CO2'][0:151].values

# Get historic atmospheric concentrations: CH4
# Data in scenarios common to 2005
FILE_CH4    = GHG_DIR+'SSP2-baseline_IMAGE_concs_ch4_n2o_vn3p0.txt'
DATA_INPUT  = pd.read_csv(FILE_CH4, header=None, names=['Year','CH4','N2O'], delim_whitespace=True)
DATA_PROF_COMMON['CH4_YEAR'] = DATA_INPUT['Year'][0:151].values 
DATA_PROF_COMMON['CH4']      = DATA_INPUT['CH4'][0:151].values

# Get non-CO2 radiative forcing for baseline scenario
FILE_Q      = GHG_DIR+'SSP2-1.9_IMAGE-baselineCH4_qnonco2_vn3p0.txt'
DATA_INPUT  = pd.read_csv(FILE_Q, header=None, names=['Year','Q'], delim_whitespace=True)
DATA_PROF_COMMON['Q_YEAR']['CTL'] = DATA_INPUT['Year'][:].values
DATA_PROF_COMMON['Q']['CTL']      = DATA_INPUT['Q'][:].values

FILE_Q      = GHG_DIR+'SSP2-1.9_IMAGE_qnonco2_vn3p0.txt'
DATA_INPUT  = pd.read_csv(FILE_Q, header=None, names=['Year','Q'], delim_whitespace=True)
DATA_PROF_COMMON['Q_YEAR']['CH4'] = DATA_INPUT['Year'][:].values
DATA_PROF_COMMON['Q']['CH4']      = DATA_INPUT['Q'][:].values

# Plot
ALPHA          =    0.2
LINE_WIDTH     =    1.0
DPI            =  300
LEGEND_POS     =   -1

xMIN,xMAX,xINC = 1850,2100,50
xTICKS         = get_ticks(xMIN,xMAX,xINC)
xPLOT_BASE     = START_YEAR+np.arange(nYEARS)
xLABEL         = 'Year'
TEMP_TEXT      = [ '1.5$^o$C','2.0$^o$C']
KEYS_SPECIES   = { 'CO2': ['CO2','CO2_YEAR','CO$_2$'], 'CH4': ['CH4','CH4_YEAR','CH$_4$'] }
PLOT_CODES_T   = [ 'black' , 'blue' , 'darkorange', 'green', 'purple' ]
PLOT_CODES_RF  = [ 'black' , 'green', 'purple', 'blue' , 'darkorange' ]
PLOT_CODES_ALL = { 'CTL': ['Control','black',], 'CH4': ['CH$_4$','blue'], \
                   'Hist': ['Historical', 'black'], \
                   'LULUC_CCS': ['BECCS', 'red'], 'LULUC_Nat': ['Natural', 'green'], \
                   'Coupled_CCS': ['BECCS+CH$_4$', 'pink'], 'Coupled_Nat': ['Natural+CH$_4$', 'purple'] }
PLOT_DATA_AXES = { 'T':   [  -0.25 ,   2.25,   0.25,'Temperature Anomaly ($^o$C)'], \
                   'Q':   [  -0.125,   1.125,  0.125,'non-CO$_2$ Radiative Forcing (W m$^{-2}$)'], \
                   'CH4': [ 500,    2500,    500,   'Atmospheric CH$_4$ (ppbv)'], \
                   'CO2': [ 250,     550,     50,   'Atmospheric CO$_2$ (ppmv)'] }

if PLATFORM == 'CEH' and PROFILE_npy == 'Save':
    KEYS_SCEN      = { '1': [ 'CTL vs CH$_4$', ['CTL']], \
                       '2': [ 'BECCS vs BECCS+CH$_4$', ['LULUC_CCS']], \
                       '3': [ 'Natural vs Natural+CH$_4$', ['LULUC_Nat']] }
else:
    KEYS_SCEN      = { '1': [ 'CTL vs CH$_4$', ['CTL','CH4']], \
                       '2': [ 'BECCS vs BECCS+CH$_4$', ['LULUC_CCS','Coupled_CCS']], \
                       '3': [ 'Natural vs Natural+CH$_4$', ['LULUC_Nat','Coupled_Nat']] }

# Multi plot profile figures for paper
# (1) Temperature and RF profiles
# Temperature
FIGURE,AXIS    = plt.subplots(ncols=1,nrows=1,figsize=[5,4])
yMIN,yMAX,yINC,yLABEL = \
    PLOT_DATA_AXES['T'][0],PLOT_DATA_AXES['T'][1],PLOT_DATA_AXES['T'][2],PLOT_DATA_AXES['T'][3]
yTICKS       = get_ticks(yMIN,yMAX,yINC)
LEGENDS      = [ 'Historical', TEMP_TEXT[0], TEMP_TEXT[1] ]
LINES        = []

for iTEMP in range(nTEMPs):

    TEMP         = TEMPs[iTEMP]
    xPLOT        = DATA_PROF_COMMON['T_YEAR'][TEMP]
    yPLOT        = DATA_PROF_COMMON['T'][TEMP]

    LINE,        = AXIS.plot(xPLOT,yPLOT,PLOT_CODES_T[iTEMP+1])
    LINES.append(LINE)
    LINE,        = AXIS.plot(xPLOT[xPLOT<=2000],yPLOT[xPLOT<=2000],PLOT_CODES_T[0])
    AXIS.plot([xMIN,xMAX],[0.0,0.0],'k--')

#AXIS.set_title('Temperature',fontsize=10)
AXIS.set_xlim([xMIN,xMAX])
AXIS.set_xticks(xTICKS)
AXIS.set_xticklabels(xTICKS, fontsize=10)
AXIS.set_xlabel(xLABEL, fontsize=10)
AXIS.set_ylim([yMIN,yMAX])
AXIS.set_yticks(yTICKS)
AXIS.set_yticklabels(yTICKS, fontsize=10)
AXIS.set_ylabel(yLABEL, fontsize=10)
LINES        = [ LINE ]+LINES
AXIS.legend((LINES),(LEGENDS),loc='upper left',frameon='False',fontsize=9)

FIGURE.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.10, \
               hspace=0.05, wspace=0.05)

if PLOT_OPT == '1':
    plt.show()
else:
    FILE_PLOT      = PLOT_DIR+'Profiles_T_'+PLATFORM+'_'+sDATE
    print(FILE_PLOT)
    FIGURE.savefig(FILE_PLOT+'.png',dpi=DPI)
    FIGURE.savefig(FILE_PLOT+'.pdf',dpi=DPI)
    FIGURE.savefig(FILE_PLOT+'.eps',dpi=DPI)

# Radiative Forcing
FIGURE,AXIS    = plt.subplots(ncols=1,nrows=1,figsize=[5,4])
yMIN,yMAX,yINC,yLABEL = \
    PLOT_DATA_AXES['Q'][0],PLOT_DATA_AXES['Q'][1],PLOT_DATA_AXES['Q'][2],PLOT_DATA_AXES['Q'][3]
yTICKS       = get_ticks(yMIN,yMAX,yINC)
LEGENDS      = [ 'Historical', 'Control (CTL)', 'CH$_4$ Mitigation' ]
LINES        = []

iKEY         = 0
for KEY in KEYS_SCEN['1'][1]:

    xPLOT        = DATA_PROF_COMMON['Q_YEAR'][KEY]
    yPLOT        = DATA_PROF_COMMON['Q'][KEY]

    LINE,        = AXIS.plot(xPLOT,yPLOT,PLOT_CODES_RF[iKEY+1])
    LINES.append(LINE)
    LINE,        = AXIS.plot(xPLOT[xPLOT<=2000],yPLOT[xPLOT<=2000],PLOT_CODES_T[0])
    AXIS.plot([xMIN,xMAX],[0.0,0.0],'k--')

    iKEY +=1

#AXIS.set_title('Radiative Forcing',fontsize=10)
AXIS.set_xlim([xMIN,xMAX])
AXIS.set_xticks(xTICKS)
AXIS.set_xticklabels(xTICKS, fontsize=10)
AXIS.set_xlabel(xLABEL, fontsize=10)
AXIS.set_ylim([yMIN,yMAX])
AXIS.set_yticks(yTICKS)
AXIS.set_yticklabels(yTICKS, fontsize=10)
AXIS.set_ylabel(yLABEL, fontsize=10)
LINES        = [ LINE ]+LINES
AXIS.legend((LINES),(LEGENDS),loc='upper left',frameon='False',fontsize=9)

FIGURE.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.10, \
               hspace=0.05, wspace=0.05)

if PLOT_OPT == '1':
    plt.show()
else:
    FILE_PLOT      = PLOT_DIR+'Profiles_RF_'+PLATFORM+'_'+sDATE
    print(FILE_PLOT)
    FIGURE.savefig(FILE_PLOT+'.png',dpi=DPI)
    FIGURE.savefig(FILE_PLOT+'.pdf',dpi=DPI)
    FIGURE.savefig(FILE_PLOT+'.eps',dpi=DPI)

# (2) CH4 and CO2 profiles
FIGURE,AXES  = plt.subplots(ncols=2,nrows=3,figsize=[10,10])
FIGURE.subplots_adjust(top=0.97,left=0.1,right=0.95,bottom=0.05,hspace=0.2)
iCOL         = 0
xPLOT        = xPLOT_BASE

for KEY_SP in ['CH4','CO2']:

    iROW          = 0
    print(KEY_SP,iROW,iCOL)
    yMIN,yMAX,yINC,yLABEL = \
        PLOT_DATA_AXES[KEY_SP][0],PLOT_DATA_AXES[KEY_SP][1],PLOT_DATA_AXES[KEY_SP][2],PLOT_DATA_AXES[KEY_SP][3]
    yTICKS        = get_ticks(yMIN,yMAX,yINC)

    # for KEY_SC in ['Methane','Methane+BECCS','Methane+Natural']:
    for KEY_SC in sorted(KEYS_SCEN.keys()):
        print(iROW,iCOL)
        AXIS          = AXES[iROW][iCOL]
        AXIS.set_title(KEYS_SCEN[KEY_SC][0],fontsize=10)
        AXIS.set_xlim([xMIN,xMAX])
        AXIS.set_xticks(xTICKS)
        AXIS.set_xticklabels(xTICKS, fontsize=10)
        AXIS.set_xlabel(xLABEL, fontsize=10)
        AXIS.set_ylim([yMIN,yMAX])
        AXIS.set_yticks(yTICKS)
        AXIS.set_yticklabels(yTICKS, fontsize=10)
        AXIS.set_ylabel(yLABEL, fontsize=10)

        iLINE         = 1
        LEGENDS,LINES = [],[]

        for iTEMP in range(nTEMPs):
            TEMP          = TEMPs[iTEMP]
            for comp in KEYS_SCEN[KEY_SC][1]:
                DATA_PLOT             = np.zeros((nYEARS,nGCMs,nQ10s,nO3s))
                yPLOT,yPLOT_L,yPLOT_H = np.zeros((nYEARS)),np.zeros((nYEARS)),np.zeros((nYEARS))
                LEGENDS.append(TEMP_TEXT[iTEMP]+': '+PLOT_CODES_ALL[comp][0])
          
                for igcm in range(nGCMs):
                    gcm = GCMs[igcm]
                    for iQ10 in range(nQ10s):
                        Q10 = Q10_exps[iQ10]
                        for iO3 in range(nO3s):
                            O3 = Ozone_exps[iO3]
                            DATA_PLOT[:,igcm,iQ10,iO3] = COMPS[comp][TEMP][gcm][Q10][O3[1]][KEY_SP][:]

                for iyear in range(nYEARS):
                    yPLOT_TEMP     = DATA_PLOT[iyear,:,:,:] 
                    yPLOT[iyear]   = np.median(yPLOT_TEMP[~np.isnan(yPLOT_TEMP)])
                    yPLOT_L[iyear] = np.percentile(yPLOT_TEMP[~np.isnan(yPLOT_TEMP)],25.0)
                    yPLOT_H[iyear] = np.percentile(yPLOT_TEMP[~np.isnan(yPLOT_TEMP)],75.0)

                print(comp,': ',yPLOT[0],yPLOT[-1])
                if DEBUG == 'Y':
                    print(DATA_PLOT.shape); print('x: ',xPLOT); print('y: ',yPLOT); \
                    print('yl: ',yPLOT_L); print('yh: ',yPLOT_H)

#               LINE,         = AXIS.plot(xPLOT,yPLOT,PLOT_CODES_ALL[comp][1],lw=LINE_WIDTH)
#               AXIS.fill_between(xPLOT,yPLOT_L,yPLOT_H,color=PLOT_CODES_ALL[comp][1],alpha=ALPHA)
                LINE,         = AXIS.plot(xPLOT,yPLOT,PLOT_CODES_T[iLINE],lw=LINE_WIDTH)
                AXIS.fill_between(xPLOT,yPLOT_L,yPLOT_H,color=PLOT_CODES_T[iLINE],alpha=ALPHA)
                LINES.append(LINE)

                del DATA_PLOT,yPLOT,yPLOT_L,yPLOT_H

                iLINE        += 1

        iROW         += 1

        xHIST         = DATA_PROF_COMMON[KEYS_SPECIES[KEY_SP][1]]+0.5
        yHIST         = DATA_PROF_COMMON[KEYS_SPECIES[KEY_SP][0]]
        LINE,         = AXIS.plot(xHIST[xHIST<=2000],yHIST[xHIST<=2000],PLOT_CODES_ALL['Hist'][1])
        LINES         = [ LINE ]+LINES; LEGENDS       = ['Historical']+LEGENDS
        AXIS.legend((LINES),(LEGENDS),loc='upper left',frameon='False',fontsize=9)
#       AXIS.legend((LEGENDS),loc='upper left',frameon='False',fontsize=9)

    iCOL         += 1

AXES[0][0].text(0.00,1.03,'(a)', transform=AXES[0][0].transAxes, fontsize=10, fontweight='bold')
AXES[0][1].text(0.00,1.03,'(b)', transform=AXES[0][1].transAxes, fontsize=10, fontweight='bold')
AXES[1][0].text(0.00,1.03,'(c)', transform=AXES[1][0].transAxes, fontsize=10, fontweight='bold')
AXES[1][1].text(0.00,1.03,'(d)', transform=AXES[1][1].transAxes, fontsize=10, fontweight='bold')
AXES[2][0].text(0.00,1.03,'(e)', transform=AXES[2][0].transAxes, fontsize=10, fontweight='bold')
AXES[2][1].text(0.00,1.03,'(f)', transform=AXES[2][1].transAxes, fontsize=10, fontweight='bold')

FIGURE.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05, \
                hspace=0.25, wspace=0.25)

if PLOT_OPT == '1':
    plt.show()
else:
    FILE_PLOT      = PLOT_DIR+'Profiles_CH4_CO2_'+PLATFORM+'_'+sDATE
    print(FILE_PLOT)
    FIGURE.savefig(FILE_PLOT+'.png',dpi=DPI)
    FIGURE.savefig(FILE_PLOT+'.pdf',dpi=DPI)
    FIGURE.savefig(FILE_PLOT+'.eps',dpi=DPI)

plt.close()
