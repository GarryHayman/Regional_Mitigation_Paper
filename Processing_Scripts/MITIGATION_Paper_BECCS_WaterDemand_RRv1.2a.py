#!/bin/env python2.7

import matplotlib as mpl
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

BECCS_WaterDem_Factor = 530  # * 1e-3 * 1e3   # 530 m3 t-1C => 530 kg

secs_to_year = 3600.*24.*360. 
kg_to_Gt     = 1e-12
kg_to_Mt     = 1e-9
kg_to_Tg     = kg_to_Mt
kg_to_Gg     = 1e-6
m2_to_Mha    = 1e-10
GtC_to_ppm   = 0.471
C_to_CH4     = 16.04/12.011
ppm_to_kgC   = 1e12/GtC_to_ppm
C_to_water   = 530.0*1.0E-03 # Convert from GtC yr-1 to Tm3 yr-1
kgH2O_to_m3  = 1e-3  
kgH2O_to_km3 = kgH2O_to_m3*1e-9
npft         = 13

# Optional input parameters
PLOT_FIGURE   = False
DEBUG         = optional_argparse('-debug','N')
subregions    = optional_argparse('-subregions','IMAGE').upper()
PLATFORM      = optional_argparse('-platform','JASMIN')
sDATE         = optional_argparse('-date',dt.datetime.strftime(dt.datetime.now(),'%Y%m%d'))
PLOT_TAG      = optional_argparse('-plottag', 'Water_Demand')
PLOT_OPT      = str(optional_argparse('-plot_opt','1'))
PROC_STEP     = optional_argparse('-BECCS_NPY','Save')
CCS_minimum_threshold     = float(optional_argparse('CCS_min','1e-4'))
BIOFRAC_minimum_threshold = float(optional_argparse('CCSfrac_min','1e-2'))
BECCSreg_minimum_threshold = float(optional_argparse('-BECCSreg_min','1'))  # GtC by 2100
sBECCS_multiplier          = optional_argparse('-beccs_multiplier','1')
BECCS_multiplier           = float(sBECCS_multiplier)
print('BECCS_multiplier: ', BECCS_multiplier)

# Directories containing JULES output and plot output directories:
if PLATFORM == 'JASMIN':
    HOME_DIR      = '/gws/nopw/j04/clifftop/SYNTHESIS/PostReview_output/'
    DATA_DIR      = HOME_DIR
    ANCILS_DIR    = '/gws/nopw/j04/clifftop/COMMON_DATA/ANCILS/'
    PLOT_DIR      = optional_argparse('-plotdir',DATA_DIR+'plots/'+PLOT_TAG+'_BECCS'+sBECCS_multiplier+'.0/')
    COMPS_npy_DIR = '/work/scratch/garr_output/SYNTHESIS/npy_files/processed_output/'

    Q10_exps      = [ 'lowQ10', 'highQ10' ]
    Ozone_exps    = [['L','lowO3'], ['H','highO3'], ]

    COMPS_keys    = [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]

elif PLATFORM == 'CEH':
    HOME_DIR      = '/prj/CLIFFTOP/SYNTHESIS/'
    PLOT_DIR      = '/data/grp/eow/garr/Projects/CLIFFTOP/Paper_Synthesis/PLOTS_WATER/CEH/'
    DATA_DIR      = HOME_DIR+'Review_Response_Check/GCM_Output/'
    ANCILS_DIR    = HOME_DIR+'Land_Cover/'
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
print('DATA_DIR: '+DATA_DIR)
print('PLOT_DIR: '+PLOT_DIR)

PRESENT_DAY_YEAR = 2015

ALPHABET       = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o']

START_YEAR     = PRESENT_DAY_YEAR
END_YEAR       = 2100
nYEARS         = END_YEAR-START_YEAR # No data for 2100

# Select Scenarios to plot:
START_SCEN       = '1p5deg'
START_YEAR       = PRESENT_DAY_YEAR

TEMPs            = ['1p5deg', '2deg' ] #'1p81p5deg',   # tag in the JULES output file directory
TEMP_years       = [2099,2099,2099]    # tag in the JULES output file directory
TEMP_names       = ['1.5$^o$C','2.0$^o$C'] # '1.5$^o$C Overshoot (2100)',
                # Name to appear on plots etc.

land_pools       = [ 'Land','CS','CV','WP','CCS' ]
pools            = ['Total','Atmos','Ocean',]+land_pools #\
#                 +['BECCS_productivity','Harvest','kappaHarvest','kappamefficHarvest']

# Water demand stuff:
water_years      = [2015,2060,2099]
water_SSP_years  = [2015,2060,2100]
water_tags       = [ '_t'+str(w_year) for w_year in water_years ]
water_meanlength = 3
n_watersteps     = len(water_years)
extra            = [ 'Runoff', 'BECCS', ]
#ipdb.set_trace()
extras           = [ ext+w_tag for ext in extra for w_tag in water_tags ]
#for w_tag in water_tags: extras += [ ext+w_tag for ext in extra ] 

nTEMPs           = len(TEMPs)
nQ10s            = len(Q10_exps)
nO3s             = len(Ozone_exps)
npools           = len(pools)


# Directory of Ocean Uptake data:
OCEAN_UPTAKE_DIR = optional_argparse('-ocean_uptake_dir',DATA_DIR)
OCEAN_START_YEAR = 1850
print("Ocean Uptake data from: " + OCEAN_UPTAKE_DIR)

REGION_dict      = data_info.REGION_DICTIONARIES()[subregions]
Nregions         = len(REGION_dict['Name'])

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
                    COMPS[comp][TEMP][Q10][O3[1]]={ pool:{} for pool in pools+extras } #+land_opt_pools }
                    config=COMPS[comp]['config'].replace('OZONESUB',O3[1])
                    runid=COMPS[comp]['runid'].replace('OZONESUB',O3[0])
                    print('GCM Input:',comp,TEMP,TEMP_year,Q10,config,runid)
                    for igcm in range(nGCMs):
                        gcm=GCMs[igcm]
                        gcm_index=GCM_index[igcm]
                        print(gcm)

                        DUMP_FILE=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+TEMP_tag+'.dump.YYYY0101.0.nc'
                        Ann_File=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+TEMP_tag+'.Annual_carbon.YYYY.nc'
                        H2O_File=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+TEMP_tag+'.Annual_h2o_t.YYYY.nc'

                        if DEBUG == 'Y':
                            print(DUMP_FILE)
                            print(Ann_File)
                            print(H2O_File)

                        # Open end of sim files:
                        Dinf   = nc.Dataset(DUMP_FILE.replace('YYYY',str(TEMP_year+1)),'r')
                        Ainf   = nc.Dataset(Ann_File.replace('YYYY',str(TEMP_year)),'r')
                        # Open Present day files:
                        Dinf_0 = nc.Dataset(DUMP_FILE.replace('YYYY',str(START_YEAR+1)),'r')
                        Ainf_0 = nc.Dataset(Ann_File.replace('YYYY',str(START_YEAR)),'r')
                        
                        # Read in data, find delta from present day,    
                        CV     = ( Ainf.variables['cv'][:]-Ainf_0.variables['cv'][:] ).squeeze() * AREA_1D
                        COMPS[comp][TEMP][Q10][O3[1]]['CV'][gcm] = CV *kg_to_Gt
                        
                        CS     = ( Ainf.variables['cs_gb'][:]-Ainf_0.variables['cs_gb'][:] ).squeeze() * AREA_1D
                        COMPS[comp][TEMP][Q10][O3[1]]['CS'][gcm] = CS*kg_to_Gt
                        
                        # Correction for bug in test runs, can remove this with final runs, although will not affect results
                        if '_Nat' in comp:
                            CCS = np.zeros_like(CS)
                        else:
                            CCS = (Ainf.variables['ccs_gb'][:]-Ainf_0.variables['ccs_gb'][:]).squeeze()*AREA_1D*BECCS_multiplier 
                            CCS[np.isfinite(CCS)==False]=0.
                        COMPS[comp][TEMP][Q10][O3[1]]['CCS'][gcm] = CCS*kg_to_Gt
                        
                        WP =(( Dinf.variables['wood_prod_fast'][:]+Dinf.variables['wood_prod_med'][:] +
                               Dinf.variables['wood_prod_slow'][:]  )
                           - ( Dinf_0.variables['wood_prod_fast'][:]+Dinf_0.variables['wood_prod_med'][:] +
                               Dinf_0.variables['wood_prod_slow'][:]  ) )  *  AREA_1D
                        WP[np.isfinite(WP)==False]=0.
                        COMPS[comp][TEMP][Q10][O3[1]]['WP'][gcm] = WP*kg_to_Gt

                        AtmCO2_ppm = Dinf.variables['co2_ppmv'][0]-Dinf_0.variables['co2_ppmv'][0]
                        AtmCO2_kg  = AtmCO2_ppm*ppm_to_kgC
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
                        
                        # Loop over years to calculate the water demand for.
                        for iyear in range(water_meanlength):
                            H2Oinfs = [ nc.Dataset(H2O_File.replace('YYYY',str(w_year-iyear)),'r') for w_year in water_years ]
                            CCSinfs = [ nc.Dataset(Ann_File.replace('YYYY',str(w_year-iyear)),'r') for w_year in water_years ]
                            if iyear==0:
                                RUNOFFS = [ H2Oinf.variables['runoff'][:].squeeze()*secs_to_year*kgH2O_to_km3 for H2Oinf in H2Oinfs ]
                                BECCS_0 = [ CCSinf.variables['ccs_gb'][:].squeeze()*BECCS_multiplier*kg_to_Gt for CCSinf in CCSinfs ]
                            else:
                                for itag in range(n_watersteps):
                                    RUNOFFS[itag] += H2Oinfs[itag].variables['runoff'][:].squeeze()*secs_to_year*kgH2O_to_km3

                        BECCS_1 = [  CCSinf.variables['ccs_gb'][:].squeeze()*BECCS_multiplier*kg_to_Gt  for CCSinf in CCSinfs ]
                        for ibeccs in range(len(CCSinfs)): 
                            BECCS_0[ibeccs][np.isfinite(BECCS_0[ibeccs])==False]=0.
                            BECCS_1[ibeccs][np.isfinite(BECCS_1[ibeccs])==False]=0.
                        BECCS = [(beccs_0 - beccs_1) for beccs_0,beccs_1 in zip(BECCS_0,BECCS_1) ]

                        for itag in range(n_watersteps):
                            COMPS[comp][TEMP][Q10][O3[1]]['Runoff'+water_tags[itag]][gcm] = RUNOFFS[itag]*AREA_1D/water_meanlength
                            COMPS[comp][TEMP][Q10][O3[1]]['BECCS'+water_tags[itag]][gcm] = BECCS[itag]*AREA_1D/(water_meanlength-1)

                        for H2Oinf in H2Oinfs: H2Oinf.close()
                        del H2Oinfs
                        del RUNOFFS 
                        for CCSinf in CCSinfs: CCSinf.close()
                        del CCSinfs
                        
                        Dinf.close(); Ainf.close(); Dinf_0.close(); Ainf_0.close()

    # Need to create dummy datasets for CH4,
    # Full set: COMPS_keys    = [ 'CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'Coupled_CCS', 'Coupled_Nat' ]
    # CEH  set: COMPS_keys    = [ 'CTL', 'LULUC_CCS', 'LULUC_Nat' ]
    if PLATFORM == 'CEH':
        COMPS['CH4']         = deepcopy(COMPS['CTL'])
        COMPS['Coupled_CCS'] = deepcopy(COMPS['LULUC_CCS'])
        COMPS['Coupled_Nat'] = deepcopy(COMPS['LULUC_Nat'])

    # ipdb.set_trace()
    # Create optimised LULUC mitigation option by choosing CCS or return to Natural Vegetation:
    # We use the atmosphere and ocean uptake components from the CCS simulation. There are slight
    # differences between Natural and CCS simualtions due to changes in the natural wetalnd contribution
    # but this is of the order 1-2 GtC so we ignore for the time being.

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
                        if DEBUG == 'Y': print('Optimisation: ',comp,TEMP,Q10,O3,gcm)
                        # create mask, CCS>Nat = 1; CCS==Nat = 0; CCS<Nat = -1
                        difference = ( ( COMPS[compCCS][TEMP][Q10][O3[1]]['Land'][gcm]
                                       + COMPS[compCCS][TEMP][Q10][O3[1]]['CCS'][gcm] )
                                     - ( COMPS[compNat][TEMP][Q10][O3[1]]['Land'][gcm]
                                       + COMPS[compNat][TEMP][Q10][O3[1]]['CCS'][gcm] )  )
                        flag_mask = (difference/np.abs(difference)).astype(int)
                        flag_mask[difference==0.] = 0

                        # Substitute Nat for CCS in land/CCS arrays where Nat is the prefered choice:
                        for pool in land_pools+extras:
                            COMPS[comp][TEMP][Q10][O3[1]][pool][gcm][flag_mask==-1] = \
                                    COMPS[compNat][TEMP][Q10][O3[1]][pool][gcm][flag_mask==-1]

                        # Recalculate the Total emission budget:
                        COMPS[comp][TEMP][Q10][O3[1]]['Total'][gcm] = (
                                                 COMPS[comp][TEMP][Q10][O3[1]]['Land'][gcm].sum()
                                               + COMPS[comp][TEMP][Q10][O3[1]]['CCS'][gcm].sum()
                                               + COMPS[comp][TEMP][Q10][O3[1]]['Atmos'][gcm]
                                               + COMPS[comp][TEMP][Q10][O3[1]]['Ocean'][gcm] )

                

# Create Linear Sum comparison
for Lin_comp in ['Nat','CCS','opt']:
    # Start with the CH4 comparison:
    COMPS['Linear_'+Lin_comp] = deepcopy(COMPS['CH4'])
    for iTEMP in range(nTEMPs): 
        TEMP=TEMPs[iTEMP]
        for iQ10 in range(nQ10s):
            Q10 = Q10_exps[iQ10]
            for iO3 in range(nO3s):
                O3 = Ozone_exps[iO3]
                for igcm in range(nGCMs):
                    gcm=GCMs[igcm]
                    #print(gcm)
                    # Then add on the difference between the LULUC simulation and the control simulation
                    for pool in pools+extras:
                        COMPS['Linear_'+Lin_comp][TEMP][Q10][O3[1]][pool][gcm] += (
                            COMPS['LULUC_'+Lin_comp][TEMP][Q10][O3[1]][pool][gcm] 
                          - COMPS['CTL'][TEMP][Q10][O3[1]][pool][gcm] )

if True:
    Demand_colours = ['#d6604d','#d7d7d7','#d1e5f0','#92c5de','#4393c3']   # 
    Demand_CLEVELS = [-50, -25, -10, 10, 25, 50, 100, 250, 500, 1000, 2000 ]
    in_WDvars      = ['BECCS','Runoff'] #,'Precip']
    WDvars         = ['SSP_WD','SSP_Irrig','BECCS_opt_WD','Water_Avail']+in_WDvars 
    comp           = 'LULUC_opt'
    
    #ensemble_data = { TEMP: { str(t_year): 
    #                    {var: { } for var in WDvars}  for t_year in water_years }
    #                  for TEMP in TEMPs } 
    ensemble_data  = { TEMP: { var: { region: { str(t_year): [] for t_year in water_years } 
                      for region in REGION_dict['Name'] }  for var in in_WDvars}  for TEMP in TEMPs } 
    BECCS_Region_index = {TEMP:[] for TEMP in TEMPs } 


for TEMP in TEMPs:
    for var in in_WDvars:
        for iregion in range(Nregions):
            region=REGION_dict['Name'][iregion]
            region_mask=REGIONS_1D==(iregion+1)
            if region=='Global': region_mask[:]=True
            if region=='International Transportation': region_mask[:]=False

            for itag in range(n_watersteps):
                tag = water_tags[itag]
                year = water_years[itag]
                yearSSP = water_SSP_years[itag] 

                #ensemble_data[TEMP][var][region][str(year)]=[]
                for iQ10 in range(nQ10s):
                    Q10 = Q10_exps[iQ10]
                    for iO3 in range(nO3s):
                        O3 = Ozone_exps[iO3]
                        ensemble_data[TEMP][var][region][str(year)].append( pd.Series(
                              [ COMPS[comp][TEMP][Q10][O3[1]][var+tag][gcm][region_mask].sum() for gcm in GCMs ], index=GCMs ) )
                  
                ensemble_data[TEMP][var][region][str(year)] = pd.concat( ensemble_data[TEMP][var][region][str(year)] ) 

            ensemble_data[TEMP][var][region] = pd.DataFrame(ensemble_data[TEMP][var][region])
            #print(iregion, ensemble_data[TEMP]['BECCS'][region].max().max())
            CCS_tot_median = np.median(  #pd.Series(
                    [COMPS[comp][TEMP][Q10][O3[1]]['CCS'][gcm][region_mask].sum() for gcm in GCMs]) #, index=GCMs ) 
            #if ( (ensemble_data[TEMP]['BECCS'][region].max().max()>BECCSreg_minimum_threshold) 
            if ( (CCS_tot_median >BECCSreg_minimum_threshold) 
               & (region!='Global') & (iregion not in BECCS_Region_index[TEMP]) ):
                print(region, CCS_tot_median)
                BECCS_Region_index[TEMP].append(iregion)

        ensemble_data[TEMP][var] = pd.concat(ensemble_data[TEMP][var], axis=1)
      
    # Now populate dictionary with derived variables
    ensemble_data[TEMP]['BECCS_WD']    = ensemble_data[TEMP]['BECCS'] * BECCS_WaterDem_Factor
        
    ensemble_data[TEMP]['Water_Avail'] = deepcopy(ensemble_data[TEMP]['Runoff'])
    for iregion in range(Nregions): 
        ensemble_data[TEMP]['Water_Avail'][REGION_dict['Name'][iregion]] *= REGION_dict['WaterAccessibility'][iregion]

    # Create panda dataframes for the SSP data read from data_info (provided by Jonathan Doelman)
    SSP_WD_dict    = {}
    SSP_Irrig_dict = {}
    for iregion in range(Nregions):
        region = REGION_dict['Name'][iregion]
        SSP_WD_dict[region] = pd.DataFrame(
                { str(water_years[itag]): [ REGION_dict['WaterWithdrawal'][str(water_SSP_years[itag])][iregion] ]
                       for itag in range(n_watersteps) }, index = ['SSP2-baseline'] )
        SSP_Irrig_dict[region] = pd.DataFrame(
                { str(water_years[itag]): [ REGION_dict['Irrigation'][str(water_SSP_years[itag])][iregion] ]
                       for itag in range(n_watersteps) }, index = ['SSP2-baseline'] )

    ensemble_data[TEMP]['SSP_WD']    = pd.concat(SSP_WD_dict, axis=1)
    ensemble_data[TEMP]['SSP_Irrig'] = pd.concat(SSP_Irrig_dict, axis=1)
    
    if DEBUG == 'Y': print( ensemble_data[TEMP]['SSP_Irrig'] )
    # ipdb.set_trace()
    # Don't convert into one big dataframe per TEMP as SSP data also included in these dictionaries
    # ensemble_data[TEMP] = pd.concat(ensemble_data[TEMP], axis=1)

if True:
    outf=open(PLOT_DIR+'Water_Demand_Data_'+sDATE+'.csv','w')
    
    for TEMP in TEMPs:
        output_pandas = []
        for var in ensemble_data[TEMP].keys():
            temp_pd       =  ensemble_data[TEMP][var].describe()[3:8]
            temp_pd.index = [ [var for i in range(5)],temp_pd.index ]
            output_pandas.append( temp_pd )
        
        output_panda = pd.concat(output_pandas)
        
        outf.write('# '+TEMP+':\n')
        for year in water_years:

            try:
                 outpd = output_panda.loc(axis=1)[:,str(year):str(year)]
                 outpd.transpose().to_csv(outf,float_format='%10.2f')
                 print('csv written as per original code')
            except:
                for iregion in range(Nregions):
                    outpd = output_panda[REGION_dict['Name'][iregion]][str(year)]
                    outpd.transpose().to_csv(outf,float_format='%10.2f')
                print('csv written as per modified code')

    outf.close()


############################################################################################################
# Now plot all regions on one plot:
# Bar plot of Lever impacts:
#ipdb.set_trace()
if True:
    Bar_list      = ['Water Availability', 'Water Demand']
    plotvar_list  = ['Water_Avail','SSP_WD','BECCS_WD']
    legend_names  = ['Water Availability','Water Demand: Other', \
                     'Water Demand: BECCS Additional'] 
    Bar_Colours   = ['#a6cee3','#1f78b4','#b2df8a']
    TEMP_space    = 0.8 # how much space all the bars should take up for a TEMP 
                        # 1 means TEMP will touch next TEMP
    nbars         = len(Bar_list)  # Number of bars per TEMP (i.e. Veg, Soil, Amos)
    bar_fraction  = 0.95 # fraction of bar space to fill with bar, 1 means bars touch
    bar_space     = (TEMP_space/nbars) # total space for a bar
    bar_width     = (bar_space)*(bar_fraction) # width of bar
    bar_gap       = (bar_space)*(1.-bar_fraction) # gap between bars
    bar_positions = [ -(TEMP_space/2.)+(i*bar_space)+(bar_gap/2.) for i in range(nbars) ]
    FONTSIZE      = 11
    
    #BECCS_Snames  = [REGION_dict['ShortName'][iregion] for iregion in BECCS_Region_index[TEMP]]

for iTEMP in range(nTEMPs):
    TEMP=TEMPs[iTEMP]
    print('Plotting Regional Water Demand (no irrigation): '+TEMP)
    nBECCS_regions=len(BECCS_Region_index[TEMP])
    BECCS_regions = [REGION_dict['Name'][iregion] for iregion in BECCS_Region_index[TEMP]]
    BECCS_Snames  = [REGION_dict['Name'][iregion].replace(' ','\n').replace('Rest\nof','Rest of')
                     for iregion in BECCS_Region_index[TEMP]]
    if DEBUG == 'Y':
       for i in range(nBECCS_regions): print(BECCS_Region_index[TEMP][i],BECCS_regions[i])
    for itag in range(n_watersteps):
        tag     = water_tags[itag]
        year    = water_years[itag]
        yearSSP = water_SSP_years[itag] 
        fig,ax=plt.subplots( ncols=1,nrows=1,figsize=(2+(1*nBECCS_regions),5) )
        fig.subplots_adjust(top=0.9,left=0.1,right=0.998,bottom=0.25)
        # Loop over regions:
        for iiregion in range(nBECCS_regions):#range(REGION_dict['Nregions']-2):
            iregion      = BECCS_Region_index[TEMP][iiregion] 
            region       = REGION_dict['Name'][iregion]
            scen_cen_pos = 0.5+iiregion
            bar_list=[]  # Append the bar objects to list for legend

            # Plot the run off as water availability
            ibar         = 0
            bar_base     = 0.0
            bar          = plotvar_list[ibar]
            xpos         = scen_cen_pos + bar_positions[ibar]
            bar_colour   = Bar_Colours[ibar]
            plotdata     = ensemble_data[TEMP][bar][region][str(year)]  
            bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
            ax.plot([xpos+bar_width/2 for i in range(len(plotdata))],plotdata,c='k',ls='',marker='.')

            # Plot the SSP availability data
            ibar         = 1
            bar          = plotvar_list[ibar]
            bar_base     = 0.0
            xpos         = scen_cen_pos + bar_positions[ibar]
            bar_colour   = Bar_Colours[ibar]
            bar_base     = 0.0
            plotdata     = ensemble_data[TEMP][bar][region][str(year)]
            bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))

            # plot the Irrigation bar on top of the Withdrawal bar
            ibar         = 2
            bar          = plotvar_list[ibar]
            bar_base     = bar_base + np.median(plotdata) # From previous bar
            bar_colour   = Bar_Colours[ibar]
            plotdata     = ensemble_data[TEMP][bar][region][str(year)]
            bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))

        ax.set_xlim([0,nBECCS_regions])
        ax.plot(ax.get_xlim(),[0,0],c='k',lw=2)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticks(np.arange(0.5,nBECCS_regions,1.))
        ax.set_xticklabels(BECCS_Snames,fontsize=FONTSIZE,rotation=60)
        ax.tick_params(axis='y',labelsize=FONTSIZE)
        ax.text(0.95,0.95,TEMP_names[iTEMP]+', '+str(yearSSP),transform=ax.transAxes,fontsize=FONTSIZE,ha='right')
        ax.set_ylim([0,ax.get_ylim()[1]])
        ax.set_ylabel('Water Availabilty/Demand (km$^3$ yr$^{-1}$)',fontsize=FONTSIZE)
        fig.legend(bar_list,legend_names, loc='upper center', fontsize=FONTSIZE, ncol=nbars+2)

        if PLOT_OPT == '1':
            plt.show()
        else:
            # Save as png
            fig.savefig(PLOT_DIR+'Water_Demand_incBECCS_ByRegion'+tag+'_'+TEMP+'_'+sDATE+'.png', bbox_inches='tight')
            # Save as eps
            fig.savefig(PLOT_DIR+'Water_Demand_incBECCS_ByRegion'+tag+'_'+TEMP+'_'+sDATE+'.eps', bbox_inches='tight')
            # Save as pdf
            fig.savefig(PLOT_DIR+'Water_Demand_incBECCS_ByRegion'+tag+'_'+TEMP+'_'+sDATE+'.pdf', bbox_inches='tight')
        plt.close()

#ipdb.set_trace()
############################################################################################################
# Now plot all regions on one plot:
# Bar plot of Lever impacts:
if True:
    Bar_list      = ['Water Availability', 'Water Demand']
    plotvar_list  = ['Water_Avail','SSP_Irrig','SSP_WD','BECCS_WD']
    legend_names  = ['Water Availability','Water Demand: Irrigation','Water Demand: Other', \
                     'Water Demand: BECCS Additional'] 
    Bar_Colours   = ['#a6cee3','#1f78b4','#ffff99','#b2df8a']
    TEMP_space    = 0.8 # how much space all the bars should take up for a TEMP 
                        # 1 means TEMP will touch next TEMP
    nbars         = len(Bar_list)  # Number of bars per TEMP (i.e. Veg, Soil, Amos)
    bar_fraction  = 0.95 # fraction of bar space to fill with bar, 1 means bars touch
    bar_space     = (TEMP_space/nbars) # total space for a bar
    bar_width     = (bar_space)*(bar_fraction) # width of bar
    bar_gap       = (bar_space)*(1.-bar_fraction) # gap between bars
    bar_positions = [ -(TEMP_space/2.)+(i*bar_space)+(bar_gap/2.) for i in range(nbars) ]
    FONTSIZE      = 11
    
    #BECCS_Snames  = [REGION_dict['ShortName'][iregion] for iregion in BECCS_Region_index[TEMP]]

for iTEMP in range(nTEMPs):
    TEMP=TEMPs[iTEMP]
    print('Plotting Plotting Regional Water Demand (irrigation): '+TEMP)
    nBECCS_regions=len(BECCS_Region_index[TEMP])
    BECCS_regions = [REGION_dict['Name'][iregion] for iregion in BECCS_Region_index[TEMP]]
    BECCS_Snames  = [REGION_dict['Name'][iregion].replace(' ','\n').replace('Rest\nof','Rest of')
                     for iregion in BECCS_Region_index[TEMP]]

    if DEBUG == 'Y':
        for i in range(nBECCS_regions): print(BECCS_Region_index[TEMP][i],BECCS_regions[i])

    for itag in range(n_watersteps):
        tag = water_tags[itag]
        year = water_years[itag]
        yearSSP = water_SSP_years[itag] 
        fig,ax=plt.subplots( ncols=1,nrows=1,figsize=(2+(1*nBECCS_regions),5) )
        fig.subplots_adjust(top=0.9,left=0.1,right=0.998,bottom=0.25)

        # Loop over regions:
        for iiregion in range(nBECCS_regions):#range(REGION_dict['Nregions']-2):
            iregion=BECCS_Region_index[TEMP][iiregion] 
            region = REGION_dict['Name'][iregion]
            scen_cen_pos=0.5+iiregion
            bar_list=[]  # Append the bar objects to list for legend

            # Plot the run off as water availability
            ibar = 0
            bar          = plotvar_list[ibar]
            xpos         = scen_cen_pos + bar_positions[ibar]
            bar_base     = 0.0
            bar_colour   = Bar_Colours[ibar]
            plotdata     = ensemble_data[TEMP][bar][region][str(year)]
            bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
            ax.plot([xpos+bar_width/2 for i in range(len(plotdata))],plotdata,c='k',ls='',marker='.')

            # Plot the SSP availability data
            ibar = 1
            bar          = plotvar_list[ibar]
            xpos         = scen_cen_pos + bar_positions[ibar]
            bar_base     = 0.0
            bar_colour   = Bar_Colours[ibar]
            plotdata     = ensemble_data[TEMP][bar][region][str(year)]  
            bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))

            # plot the SSP irrigation
            ibar = 2
            bar          = plotvar_list[ibar]
            bar_base     = bar_base + np.median(plotdata) # From previous bar
            bar_colour   = Bar_Colours[ibar]
            plotdata     = ensemble_data[TEMP][bar][region][str(year)] 
            bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))

            # plot the CCS bar on top of the SSP bar
            ibar = 3
            bar          = plotvar_list[ibar]
            bar_base     = bar_base + np.median(plotdata) # From previous bar
            bar_colour   = Bar_Colours[ibar]
            plotdata     = ensemble_data[TEMP][bar][region][str(year)]
            bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
            ax.plot([xpos+bar_width/2 for i in range(len(plotdata))],plotdata+bar_base,c='k',ls='',marker='.')

        ax.set_xlim([0,nBECCS_regions])
        ax.plot(ax.get_xlim(),[0,0],c='k',lw=2)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticks(np.arange(0.5,nBECCS_regions,1.))
        ax.set_xticklabels(BECCS_Snames,fontsize=FONTSIZE,rotation=60)
        ax.tick_params(axis='y',labelsize=FONTSIZE)
        ax.text(0.95,0.95,TEMP_names[iTEMP]+', '+str(yearSSP),transform=ax.transAxes, fontsize=FONTSIZE, ha='right')
        ax.set_ylim([0,ax.get_ylim()[1]])
        ax.set_ylabel('Water Availabilty/Demand (km$^3$ yr$^{-1}$)',fontsize=FONTSIZE)
        fig.legend(bar_list,legend_names, loc='upper center', fontsize=FONTSIZE, ncol=nbars+2)

        if PLOT_OPT == '1':
            plt.show()
        else:
            # Save as png
            fig.savefig(PLOT_DIR+'Water_Demand_Irrig_incBECCS_ByRegion'+tag+'_'+TEMP+'_'+sDATE+'.png', bbox_inches='tight')
            # Save as eps
            fig.savefig(PLOT_DIR+'Water_Demand_Irrig_incBECCS_ByRegion'+tag+'_'+TEMP+'_'+sDATE+'.eps', bbox_inches='tight')
            # Save as pdf
            fig.savefig(PLOT_DIR+'Water_Demand_Irrig_incBECCS_ByRegion'+tag+'_'+TEMP+'_'+sDATE+'.pdf', bbox_inches='tight')
        plt.close()

#ipdb.set_trace()
############################################################################################################
# Now plot all regions on one plot:
# Combined plot for 2060 and 2100
# Used for Figure 14 in published 
# Bar plot of Lever impacts:
if True:
    Bar_list      = ['Water Availability', 'Water Demand']
    plotvar_list  = ['Water_Avail','SSP_Irrig','SSP_WD','BECCS_WD']
    legend_names  = ['Water Availability','Water Demand: Irrigation','Water Demand: Other', \
                     'Water Demand: BECCS Additional'] 
    Bar_Colours   = ['#a6cee3','#1f78b4','#ffff99','#b2df8a']
    TEMP_space    = 0.8 # how much space all the bars should take up for a TEMP 
                        # 1 means TEMP will touch next TEMP
    nbars         = len(Bar_list)  # Number of bars per TEMP (i.e. Veg, Soil, Amos)
    bar_fraction  = 0.95 # fraction of bar space to fill with bar, 1 means bars touch
    bar_space     = (TEMP_space/nbars) # total space for a bar
    bar_width     = (bar_space)*(bar_fraction) # width of bar
    bar_gap       = (bar_space)*(1.-bar_fraction) # gap between bars
    bar_positions = [ -(TEMP_space/2.)+(i*bar_space)+(bar_gap/2.) for i in range(nbars) ]
    FONTSIZE      = 11
    
    #BECCS_Snames  = [REGION_dict['ShortName'][iregion] for iregion in BECCS_Region_index[TEMP]]

for iTEMP in range(nTEMPs):
    TEMP=TEMPs[iTEMP]
    print('Plotting Plotting Regional Water Demand (double plot): '+TEMP)
    nBECCS_regions=len(BECCS_Region_index[TEMP])
    BECCS_regions = [REGION_dict['Name'][iregion] for iregion in BECCS_Region_index[TEMP]]
    BECCS_Snames  = [REGION_dict['Name'][iregion].replace(' ','\n').replace('Rest\nof','Rest of')
                     for iregion in BECCS_Region_index[TEMP]]

    if DEBUG == 'Y':
        for i in range(nBECCS_regions): print(BECCS_Region_index[TEMP][i],BECCS_regions[i])

    
    fig,ax=plt.subplots( ncols=1,nrows=2,figsize=(2+(1*nBECCS_regions),12) )
#   fig.subplots_adjust(top=0.9,left=0.1,right=0.998,bottom=0.25)

    # 2060
    itag    = 1
    tag     = water_tags[itag]
    year    = water_years[itag]
    yearSSP = water_SSP_years[itag] 

    # Loop over regions:
    for iiregion in range(nBECCS_regions):#range(REGION_dict['Nregions']-2):
        iregion=BECCS_Region_index[TEMP][iiregion] 
        region = REGION_dict['Name'][iregion]
        scen_cen_pos=0.5+iiregion
        bar_list=[]  # Append the bar objects to list for legend

        # Plot the run off as water availability
        ibar = 0
        bar          = plotvar_list[ibar]
        xpos         = scen_cen_pos + bar_positions[ibar]
        bar_base     = 0.0
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax[0].bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
        ax[0].plot([xpos+bar_width/2 for i in range(len(plotdata))],plotdata,c='k',ls='',marker='.')

        # Plot the SSP availability data
        ibar = 1
        bar          = plotvar_list[ibar]
        xpos         = scen_cen_pos + bar_positions[ibar]
        bar_base     = 0.0
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax[0].bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))

        # plot the SSP irrigation
        ibar = 2
        bar          = plotvar_list[ibar]
        bar_base     = bar_base + np.median(plotdata) # From previous bar
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax[0].bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))

        # plot the CCS bar on top of the SSP bar
        ibar = 3
        bar          = plotvar_list[ibar]
        bar_base     = bar_base + np.median(plotdata) # From previous bar
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax[0].bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
        ax[0].plot([xpos+bar_width/2 for i in range(len(plotdata))],plotdata+bar_base,c='k',ls='',marker='.')

    ax[0].set_xlim([0,nBECCS_regions])
    ax[0].plot(ax[0].get_xlim(),[0,0],c='k',lw=2)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].set_xticks(np.arange(0.5,nBECCS_regions,1.))
    ax[0].set_xticklabels(BECCS_Snames,fontsize=FONTSIZE,rotation=0)
    ax[0].tick_params(axis='y',labelsize=FONTSIZE)
    ax[0].text(0.20,0.95,'(a) '+TEMP_names[iTEMP]+', '+str(yearSSP),transform=ax[0].transAxes, fontsize=FONTSIZE, ha='right')
    ax[0].set_ylabel('Water Availabilty/Demand (km$^3$ yr$^{-1}$)',fontsize=FONTSIZE)

    # 2100
    itag    = 2
    tag     = water_tags[itag]
    year    = water_years[itag]
    yearSSP = water_SSP_years[itag] 

    # Loop over regions:
    for iiregion in range(nBECCS_regions):#range(REGION_dict['Nregions']-2):
        iregion=BECCS_Region_index[TEMP][iiregion] 
        region = REGION_dict['Name'][iregion]
        scen_cen_pos=0.5+iiregion
        bar_list=[]  # Append the bar objects to list for legend

        # Plot the run off as water availability
        ibar = 0
        bar          = plotvar_list[ibar]
        bar_base     = 0.0
        xpos         = scen_cen_pos + bar_positions[ibar]
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax[1].bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
        ax[1].plot([xpos+bar_width/2 for i in range(len(plotdata))],plotdata,c='k',ls='',marker='.')

        # Plot the SSP availability data
        ibar = 1
        bar          = plotvar_list[ibar]
        bar_base     = 0.0
        xpos         = scen_cen_pos + bar_positions[ibar]
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax[1].bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))

        # plot the SSP irrigation
        ibar = 2
        bar          = plotvar_list[ibar]
        bar_base     = bar_base + np.median(plotdata) # From previous bar
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax[1].bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))

        # plot the CCS bar on top of the SSP bar
        ibar = 3
        bar          = plotvar_list[ibar]
        bar_base     = bar_base + np.median(plotdata) # From previous bar
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax[1].bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
        ax[1].plot([xpos+bar_width/2 for i in range(len(plotdata))],plotdata+bar_base,c='k',ls='',marker='.')

    ax[1].set_xlim([0,nBECCS_regions])
    ax[1].plot(ax[1].get_xlim(),[0,0],c='k',lw=2)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].set_xticks(np.arange(0.5,nBECCS_regions,1.))
    ax[1].set_xticklabels(BECCS_Snames,fontsize=FONTSIZE,rotation=0)
    ax[1].tick_params(axis='y',labelsize=FONTSIZE)
    ax[1].text(0.20,0.95,'(b) '+TEMP_names[iTEMP]+', '+str(yearSSP),transform=ax[1].transAxes, fontsize=FONTSIZE, ha='right')
    ax[1].set_ylabel('Water Availabilty/Demand (km$^3$ yr$^{-1}$)',fontsize=FONTSIZE)

    ymax = max(ax[0].get_ylim()[1],ax[1].get_ylim()[1])
    ax[0].set_ylim([0,ymax])
    ax[1].set_ylim([0,ymax])

    fig.legend(bar_list,legend_names, loc='lower center', fontsize=FONTSIZE, ncol=nbars+2)

    if PLOT_OPT == '1':
        plt.show()
    else:
        # Save as png
        fig.savefig(PLOT_DIR+'Water_Demand_Irrig_incBECCS_ByRegion_2060_2100_'+TEMP+'_'+sDATE+'.png', bbox_inches='tight')
        # Save as eps
        fig.savefig(PLOT_DIR+'Water_Demand_Irrig_incBECCS_ByRegion_2060_2100_'+TEMP+'_'+sDATE+'.eps', bbox_inches='tight')
        # Save as pdf
        fig.savefig(PLOT_DIR+'Water_Demand_Irrig_incBECCS_ByRegion_2060_2100_'+TEMP+'_'+sDATE+'.pdf', bbox_inches='tight')
    plt.close()

#ipdb.set_trace()
############################################################################################################
# Now plot Global Summary
if True:
    Bar_list      = ['Water Availability', 'Water Demand']
    plotvar_list  = ['Water_Avail','SSP_Irrig','SSP_WD','BECCS_WD']
    legend_names  = ['Water Availability','Water Demand: Irrigation','Water Demand: Other', \
                     'Water Demand: BECCS Additional'] 
    Bar_Colours   = ['#a6cee3','#1f78b4','#ffff99','#b2df8a']
    TEMP_space    = 0.8 # how much space all the bars should take up for a TEMP 
                        # 1 means TEMP will touch next TEMP
    nbars         = len(Bar_list)  # Number of bars per TEMP (i.e. Veg, Soil, Amos)
    bar_fraction  = 0.95 # fraction of bar space to fill with bar, 1 means bars touch
    bar_space     = (TEMP_space/nbars) # total space for a bar
    bar_width     = (bar_space)*(bar_fraction) # width of bar
    bar_gap       = (bar_space)*(1.-bar_fraction) # gap between bars
    bar_positions = [ -(TEMP_space/2.)+(i*bar_space)+(bar_gap/2.) for i in range(nbars) ]
    FONTSIZE      = 11

for iTEMP in range(nTEMPs):
    TEMP=TEMPs[iTEMP]
    print('Plotting Global Water Availability: '+TEMP)
    nBECCS_regions= 1
    iregion=0
    BECCS_regions = ['Global']
    region = 'Global'
    BECCS_Snames  = ['Global']

    if DEBUG == 'Y':
        for i in range(nBECCS_regions): print(BECCS_Region_index[TEMP][i],BECCS_regions[i])

    for itag in range(n_watersteps):
        tag          = water_tags[itag]
        year         = water_years[itag]
        yearSSP      = water_SSP_years[itag] 
        fig,ax       = plt.subplots( ncols=1,nrows=1,figsize=(2+(1*nBECCS_regions),6) )
        fig.subplots_adjust(top=0.75,left=0.15,right=0.98,bottom=0.15)
        # Loop over regions:
        scen_cen_pos = 0.5+iregion
        bar_list = []  # Append the bar objects to list for legend
        
        # Plot the run off and water availability
        ibar = 0
        bar          = plotvar_list[ibar]
        xpos         = scen_cen_pos + bar_positions[ibar]
        bar_base     = 0.0
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
        ax.plot([xpos+bar_width/2. for i in range(len(plotdata))],plotdata,c='k',ls='',marker='.')

        # Plot the water demand from other 
        ibar = 1
        bar          = plotvar_list[ibar]
        xpos         = scen_cen_pos + bar_positions[ibar]
        bar_base     = 0.0
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
        ax.plot([xpos+bar_width/2. for i in range(len(plotdata))],plotdata,c='k',ls='',marker='.')
        
        # plot the water demand from irrigation
        ibar = 2
        bar          = plotvar_list[ibar]
        bar_base     = bar_base + np.median(plotdata) # From previous bar
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
        ax.plot([xpos+bar_width/2. for i in range(len(plotdata))],plotdata+bar_base,c='k',ls='',marker='.')
        
        # plot the CCS bar on top of the SSP bar
        ibar = 3
        bar          = plotvar_list[ibar]
        bar_base     = bar_base + np.median(plotdata) # From previous bar
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)] 
        bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
        ax.plot([xpos+bar_width/2. for i in range(len(plotdata))],plotdata+bar_base,c='k',ls='',marker='.')
        
        ax.set_xlim([0,nBECCS_regions])
        ax.plot(ax.get_xlim(),[0,0],c='k',lw=2)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticks(np.arange(0.5,nBECCS_regions,1.))
        ax.set_xticklabels(BECCS_Snames,fontsize=FONTSIZE,rotation=0)
        ax.tick_params(axis='y',labelsize=FONTSIZE)
        ax.text(0.95,1.02,TEMP_names[iTEMP]+', '+str(yearSSP),transform=ax.transAxes, fontsize=FONTSIZE, ha='right')
        ax.set_ylim([0,ax.get_ylim()[1]])
        ax.set_ylabel('Water Availabilty/Demand (km$^3$ yr$^{-1}$)',fontsize=FONTSIZE)
        fig.legend(bar_list,legend_names, loc='upper center', fontsize=FONTSIZE, ncol=1) 

        if PLOT_OPT == '1':
            plt.show()
        else:
            # Save as png
            fig.savefig(PLOT_DIR+'Water_Demand_incBECCS_Global'+tag+'_'+TEMP+'_'+sDATE+'.png', bbox_inches='tight')
            # Save as eps
            fig.savefig(PLOT_DIR+'Water_Demand_incBECCS_Global'+tag+'_'+TEMP+'_'+sDATE+'.eps', bbox_inches='tight')
            # Save as pdf
            fig.savefig(PLOT_DIR+'Water_Demand_incBECCS_Global'+tag+'_'+TEMP+'_'+sDATE+'.pdf', bbox_inches='tight')
        plt.close()

############################################################################################################
# Now plot Global Summary for different years and warming
# Used for Figure 13 in paper (2015 and 2C output)
if True:
    Bar_list      = ['Water Availability', 'Water Demand']
    bar_pos_list  = [0,0,1,1,1]
    plotvar_list  = ['Runoff', 'Water_Avail','SSP_Irrig','SSP_WD','BECCS_WD']
    legend_names  = ['Total Runoff', 'Water Availability','Water Demand: Irrigation','Water Demand: Other', \
                     'Water Demand: BECCS Additional'] 
    Bar_Colours   = ['#a6cee3','#a6cee3','#1f78b4','#ffff99','#b2df8a']
    TEMP_space    = 0.8 # how much space all the bars should take up for a TEMP 
                        # 1 means TEMP will touch next TEMP
    nbars         = len(Bar_list)  # Number of bars per TEMP (i.e. Veg, Soil, Amos)
    bar_fraction  = 0.95 # fraction of bar space to fill with bar, 1 means bars touch
    bar_space     = (TEMP_space/nbars) # total space for a bar
    bar_width     = (bar_space)*(bar_fraction) # width of bar
    bar_gap       = (bar_space)*(1.-bar_fraction) # gap between bars
    bar_positions = [ -(TEMP_space/2.)+(i*bar_space)+(bar_gap/2.) for i in range(nbars) ]
    FONTSIZE      = 11

for iTEMP in range(nTEMPs):
    TEMP           =TEMPs[iTEMP]
    print('Plotting Global Water Availability & Runoff: '+TEMP)
    nBECCS_regions = 1
    iregion        = 0
    BECCS_regions  = ['Global']
    region         = 'Global'
    BECCS_Snames   = ['Global']

    if DEBUG == 'Y':
        for i in range(nBECCS_regions): print(BECCS_Region_index[TEMP][i],BECCS_regions[i])

    for itag in range(n_watersteps):
        tag          = water_tags[itag]
        year         = water_years[itag]
        yearSSP      = water_SSP_years[itag] 
        fig,ax       = plt.subplots( ncols=1,nrows=1,figsize=(2+(1*nBECCS_regions),6) )
        fig.subplots_adjust(top=0.75,left=0.2,right=0.98,bottom=0.1)

        # Loop over regions:
        if DEBUG == 'Y': print(region)
        scen_cen_pos = 0.5+iregion
       
        bar_list     = []  # Append the bar objects to list for legend

        # Plot the Total Runoff 
        ibar         = 0
        bar          = plotvar_list[ibar]
        xpos         = scen_cen_pos + bar_positions[bar_pos_list[ibar]]
        bar_base     = 0.0
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax.bar(xpos,np.median(plotdata),color='none',edgecolor=bar_colour,lw=3,width=bar_width,bottom=bar_base))
        ax.plot([xpos+bar_width/2. for i in range(len(plotdata))],plotdata,c='k',ls='',marker='.')

        # Plot the water availability
        ibar         = 1
        bar          = plotvar_list[ibar]
        xpos         = scen_cen_pos + bar_positions[bar_pos_list[ibar]]
        bar_colour   = Bar_Colours[ibar]
        #plotdata    = ensemble_data[TEMP][bar][region][str(year)] 
        plotdata     = np.sum([ensemble_data[TEMP][bar][region_name][str(year)] 
                                 for region_name in REGION_dict['Name'][:-1] ],axis=0 )
        bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,edgecolor='none',width=bar_width,bottom=bar_base))
        ax.plot([xpos+bar_width/2. for i in range(len(plotdata))],plotdata,c='k',ls='',marker='.')

        # Plot the SSP water demand irrigation
        ibar         = 2
        bar          = plotvar_list[ibar]
        xpos         = scen_cen_pos + bar_positions[bar_pos_list[ibar]]
        bar_base     = 0.0
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
        ax.plot([xpos+bar_width/2. for i in range(len(plotdata))],plotdata,c='k',ls='',marker='.')

        # Plot the SSP water demand other
        ibar         = 3
        bar          = plotvar_list[ibar]
        xpos         = scen_cen_pos + bar_positions[bar_pos_list[ibar]]
        bar_base     = bar_base + np.median(plotdata) # From previous bar
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
        ax.plot([xpos+bar_width/2. for i in range(len(plotdata))],plotdata+bar_base,c='k',ls='',marker='.')

        # plot the CCS bar on top of the SSP bar
        ibar         = 4
        bar          = plotvar_list[ibar]
        xpos         = scen_cen_pos + bar_positions[bar_pos_list[ibar]]
        bar_base     = bar_base + np.median(plotdata) # From previous bar
        bar_colour   = Bar_Colours[ibar]
        plotdata     = ensemble_data[TEMP][bar][region][str(year)]
        bar_list.append(ax.bar(xpos,np.median(plotdata),color=bar_colour,width=bar_width,bottom=bar_base))
        ax.plot([xpos+bar_width/2. for i in range(len(plotdata))],plotdata+bar_base,c='k',ls='',marker='.')

        ax.set_xlim([0,nBECCS_regions])
        ax.plot(ax.get_xlim(),[0,0],c='k',lw=2)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticks(np.arange(0.5,nBECCS_regions,1.))
        ax.set_xticklabels(BECCS_Snames,fontsize=FONTSIZE,rotation=0)
        ax.tick_params(axis='y',labelsize=FONTSIZE)
        ax.text(0.95,1.03,TEMP_names[iTEMP]+', '+str(yearSSP),transform=ax.transAxes, 
                fontsize=FONTSIZE, ha='right')
        ax.set_ylim([0,ax.get_ylim()[1]])
        ax.set_ylabel('Water Availabilty/Demand (km$^3$ yr$^{-1}$)',fontsize=FONTSIZE)
        fig.legend(bar_list,legend_names, 
                loc='upper center', fontsize=FONTSIZE, ncol=1) #nbars+1

        if PLOT_OPT == '1':
            plt.show()
        else:
            # Save as png
            fig.savefig(PLOT_DIR+'Water_Avail_incBECCS_Global'+tag+'_'+TEMP+'_'+sDATE+'.png', bbox_inches='tight')
            # Save as eps
            fig.savefig(PLOT_DIR+'Water_Avail_incBECCS_Global'+tag+'_'+TEMP+'_'+sDATE+'.eps', bbox_inches='tight')
            # Save as pdf
            fig.savefig(PLOT_DIR+'Water_Avail_incBECCS_Global'+tag+'_'+TEMP+'_'+sDATE+'.pdf', bbox_inches='tight')
        plt.close()

quit()
