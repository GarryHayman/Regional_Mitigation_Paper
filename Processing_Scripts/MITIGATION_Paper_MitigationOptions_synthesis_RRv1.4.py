#!/bin/env python2.7

import numpy as np
import netCDF4 as nc
import sys,os
import pandas as pd
from PlotTools import plot_tools as PTs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from copy import deepcopy

#import ipdb

import matplotlib as mpl
mpl.use('Agg')
mpl.style.use('classic')

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

# Select CH4 emission scenario for regional attribution
CH4_RCP = optional_argparse('-CH4_RCP','RCP26')
if CH4_RCP == 'RCP19':
   from imogen import data_info_CH4_RCP19 as data_info
elif CH4_RCP == 'RCP26':
   from imogen import data_info_CH4_RCP26 as data_info
else:
   from imogen import data_info

# Alphabet offset for plot labels
Alpha_Offset = int(optional_argparse('-alpha_offset','0'))

INTERACTIVE= '-interactive' in sys.argv
kg_to_Gt     = 1e-12
kg_to_Mt     = 1e-9
kg_to_Tg     = kg_to_Mt
kg_to_Gg     = 1e-6
m2_to_Mha    = 1e-10
GtC_to_ppm   = 0.471
C_to_CH4     = 16.04/12.011
ppm_to_kgC   = 1e12/GtC_to_ppm
C_to_water   = 530.0*1.0E-03 # Convert from GtC yr-1 to Tm3 yr-1
kg_to_Tm3    = 1.0E-15       # Assume 1 kg of water is 10-3 m3
npft         = 13

CCS_minimum_threshold = float(optional_argparse('CCS_min','1e-4'))
BIOFRAC_minimum_threshold = float(optional_argparse('CCSfrac_min','1e-2'))

PRESENT_DAY_YEAR=2015

ALPHABET=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o']

Tile_names=data_info.TILE_short_names()
Tile_colours = data_info.TILE_colours()
nTiles=len(Tile_names)

Q10_exps = [ 'lowQ10', 'highQ10' ]  # 'lowQ10' ] ## 
nQ10s = len(Q10_exps)

Ozone_exps = [['L','lowO3'], ['H','highO3'], ]
#Ozone_exps = [['H','highO3'], ]
nO3s = len(Ozone_exps)

COMPS_keys = ['CTL','CH4','LULUC_CCS','LULUC_Nat','Coupled_CCS','Coupled_Nat' ]
COMPS = { 
          'CTL': { 'config': 'highCH4_OZONESUB_LULUCBL', 'runid':'H_OZONESUB_BL' }
        , 'CH4': { 'config':'lowCH4_OZONESUB_LULUCBL', 'runid':'L_OZONESUB_BL' }
        , 'LULUC_CCS': { 'config':'highCH4_OZONESUB_LULUC1.9', 'runid':'H_OZONESUB_19' }
        , 'LULUC_Nat': { 'config':'highCH4_OZONESUB_LULUC1.9Nat', 'runid':'H_OZONESUB_19N' }
        , 'Coupled_CCS': { 'config':'lowCH4_OZONESUB_LULUC1.9', 'runid':'L_OZONESUB_19' }
        , 'Coupled_Nat': { 'config':'lowCH4_OZONESUB_LULUC1.9Nat', 'runid':'L_OZONESUB_19N' }
        }

#CONFIG_PLOTNAMES     = ['Standard','Permafrost Feedback','Methane Feedback','Both Feedbacks']

BECCS_multiplier = float(optional_argparse('-beccs_multiplier','1'))
print(BECCS_multiplier)
# Directories containing JULES output and plot output directories:
DATA_DIR=optional_argparse('-data_dir', 
                           '/gws/nopw/j04/clifftop/SYNTHESIS/PostReview_output/')
                           #'/gws/nopw/j04/clifftop/SYNTHESIS/ECP_output/')
                           #'/work/scratch/ecomynplatt/SYNTHESIS/PostReview_output/')
print('DATA_DIR: '+DATA_DIR)

PLOT_TAG = optional_argparse('-plottag', 'MitigationOptions_BECCS'+str(BECCS_multiplier))

PLOT_DIR = optional_argparse('-plotdir',DATA_DIR+'plots/'+PLOT_TAG+'/')
print('PLOT_DIR: '+PLOT_DIR)
os.system('mkdir -p '+PLOT_DIR )
os.system('mkdir -p '+PLOT_DIR+'Maps/')
os.system('mkdir -p '+PLOT_DIR+'BarCharts/')
#quit()
#ANNA_DIR = optional_argparse('-annadir','/gws/nopw/j04/jules/aharper/PYTHON/CLUES/npy_files/masked/')

# Directory of Ocean Uptake data:
OCEAN_UPTAKE_DIR = optional_argparse('-ocean_uptake_dir',DATA_DIR) 
OCEAN_START_YEAR=1850
print("Ocean Uptake data from: " + OCEAN_UPTAKE_DIR)

subregions = optional_argparse('-subregions','IMAGE').upper()
REGION_dict=data_info.REGION_DICTIONARIES()[subregions]

# Directory of ancillary data:
ANCILS_DIR='/gws/nopw/j04/clifftop/COMMON_DATA/ANCILS/'
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
MAXBIOFRAC_1D=Ainf.variables['maxbiofrac_19'][:].squeeze()
Ainf.close()
#print(AREA_file)

MAXBIOFRAC_2D=np.ma.masked_array(MAXBIOFRAC_1D[grindex], mask=grindex.mask)

####################################
# select GCMs:
GCMs=data_info.GCMs()
nGCMs=len(GCMs)
for igcm in range(nGCMs):
    print('%3i: '%igcm+GCMs[igcm])

if INTERACTIVE==True:
    GCM_index=raw_input('Select GCMs to plot (Press Enter to select all): ')
    if GCM_index.replace(' ','')!='':
        GCM_index = [ int(i) for i in GCM_index.split(',') ]
        GCMs=[ GCMs[i] for i in GCM_index ] 
        nGCMs=len(GCMs)
    else:
        GCM_index=[ i for i in range(nGCMs) ]
else:
    GCM_index=optional_argparse('-GCMs','ALL')
    if GCM_index=='ALL':
        GCM_index=[ i for i in range(nGCMs) ]
    else:
        GCM_index=[ int(i) for i in GCM_index.split(',') ]
        GCMs=[ GCMs[i] for i in GCM_index ] 
        nGCMs=len(GCMs)

cmip5_runs = [ data_info.cmip5_runs()[i] for i in GCM_index ]
print(' -GCMs ',GCM_index)
for igcm in range(nGCMs):
    print(igcm,GCM_index[igcm],GCMs[igcm])

# Select GCMs to plot:
#
# Select Scenarios to plot:
START_SCEN = '1p5deg'
START_YEAR = PRESENT_DAY_YEAR

SCENARIOs = ['1p5deg', '2deg' ] #'1p81p5deg',   # tag in the JULES output file directory 
SCENARIOyears = [2099,2099,2099]  # tag in the JULES output file directory 
                                   #  (see filename construction later)
SCENARIO_names=['1.5$^o$C (2100)','2.0$^o$C (2100)'] # '1.5$^o$C Overshoot (2100)',
                # Name to appear on plots etc.
nSCENARIOs=len(SCENARIOs)

# File containing the pre-industrial conditions to use as a baseline:

land_pools = [ 'Land','CS','CV','WP','CCS' ]
pools =  ['Total','Atmos','Ocean',]+land_pools #\
#        +['BECCS_productivity','Harvest','kappaHarvest','kappamefficHarvest']
npools=len(pools)

###################################################################################################
# Read in the Control data
OCEAN_START_YEAR=1850

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
    
    for iscenario in range(nSCENARIOs): 
        scenario=SCENARIOs[iscenario]
        scenyear=SCENARIOyears[iscenario]
        scentag=scenario
        COMPS[comp][scenario]={}
        for iQ10 in range(nQ10s):
          Q10 = Q10_exps[iQ10]
          COMPS[comp][scenario][Q10]={}
          for iO3 in range(nO3s):
            O3 = Ozone_exps[iO3]
            COMPS[comp][scenario][Q10][O3[1]]={ pool:{} for pool in pools  }  #+land_opt_pools }
            config=COMPS[comp]['config'].replace('OZONESUB',O3[1])
            runid=COMPS[comp]['runid'].replace('OZONESUB',O3[0])
            print(comp,scenario,scenyear,Q10,config,runid)
            for igcm in range(nGCMs):
                gcm=GCMs[igcm]
                gcm_index=GCM_index[igcm]
                print(gcm)
                DUMP_FILE=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+scentag+'.dump.YYYY0101.0.nc'
                Ann_File=DATA_DIR+Q10+'/'+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+scentag+'.Annual_carbon.YYYY.nc'
                #print(Ann_File)
                # Open end of sim files:
                Dinf = nc.Dataset(DUMP_FILE.replace('YYYY',str(scenyear+1)),'r')
                Ainf = nc.Dataset(Ann_File.replace('YYYY',str(scenyear)),'r')
                # Open Present day files:
                Dinf_0 = nc.Dataset(DUMP_FILE.replace('YYYY',str(START_YEAR+1)),'r')
                Ainf_0 = nc.Dataset(Ann_File.replace('YYYY',str(START_YEAR)),'r')
                
                # Read in data, find delta from present day,    
                CV = ( Ainf.variables['cv'][:]-Ainf_0.variables['cv'][:] ).squeeze() * AREA_1D
                COMPS[comp][scenario][Q10][O3[1]]['CV'][gcm] = CV *kg_to_Gt
                
                CS = ( Ainf.variables['cs_gb'][:]-Ainf_0.variables['cs_gb'][:] ).squeeze() * AREA_1D
                COMPS[comp][scenario][Q10][O3[1]]['CS'][gcm] = CS*kg_to_Gt
                
                # Correction for bug in test runs, can remove this with final runs, although will not affect results
                #if '_Nat' in comp:
                #    CCS = np.zeros_like(CS)
                #else:
                #    CCS = (Ainf.variables['ccs_gb'][:]-Ainf_0.variables['ccs_gb'][:]).squeeze()*AREA_1D*BECCS_multiplier 
                #    CCS[np.isfinite(CCS)==False]=0.

                CCS = (Ainf.variables['ccs_gb'][:]-Ainf_0.variables['ccs_gb'][:]).squeeze()*AREA_1D*BECCS_multiplier 
                CCS[np.isfinite(CCS)==False]=0.
                COMPS[comp][scenario][Q10][O3[1]]['CCS'][gcm] = CCS*kg_to_Gt
                
                WP =(( Dinf.variables['wood_prod_fast'][:]+Dinf.variables['wood_prod_med'][:] +
                       Dinf.variables['wood_prod_slow'][:]  )
                   - ( Dinf_0.variables['wood_prod_fast'][:]+Dinf_0.variables['wood_prod_med'][:] +
                       Dinf_0.variables['wood_prod_slow'][:]  ) )  *  AREA_1D
                WP[np.isfinite(WP)==False]=0.
                COMPS[comp][scenario][Q10][O3[1]]['WP'][gcm] = WP*kg_to_Gt

                AtmCO2_ppm = Dinf.variables['co2_ppmv'][0]-Dinf_0.variables['co2_ppmv'][0]
                AtmCO2_kg = AtmCO2_ppm*ppm_to_kgC
                COMPS[comp][scenario][Q10][O3[1]]['Atmos'][gcm] = AtmCO2_kg *kg_to_Gt
                
                Ocean =  comp_Ocean[scenario][Q10][O3[1]][gcm_index,scenyear-OCEAN_START_YEAR]   \
                       - comp_Ocean[scenario][Q10][O3[1]][gcm_index,START_YEAR-OCEAN_START_YEAR]
                COMPS[comp][scenario][Q10][O3[1]]['Ocean'][gcm] = Ocean *kg_to_Gt
                
                COMPS[comp][scenario][Q10][O3[1]]['Land'][gcm] = (  
                                                               COMPS[comp][scenario][Q10][O3[1]]['CV'][gcm]
                                                             + COMPS[comp][scenario][Q10][O3[1]]['CS'][gcm]
                                                             + COMPS[comp][scenario][Q10][O3[1]]['WP'][gcm] ) 
                COMPS[comp][scenario][Q10][O3[1]]['Total'][gcm] = (  
                                                                COMPS[comp][scenario][Q10][O3[1]]['Land'][gcm].sum()
                                                              + COMPS[comp][scenario][Q10][O3[1]]['CCS'][gcm].sum()
                                                              + COMPS[comp][scenario][Q10][O3[1]]['Atmos'][gcm]
                                                              + COMPS[comp][scenario][Q10][O3[1]]['Ocean'][gcm] )
                
                Dinf.close(); Ainf.close(); Dinf_0.close(); Ainf_0.close()

# Create optimised LULUC mitigation option by choosing CCS or return to Natural Vegetation:
#   We use the atmosphere and ocean uptake components from the CCS simulation. There are slight
#    differences between Natural and CCS simualtions due to changes in the natural wetalnd contribution
#    but this is of the order 1-2 GtC so we ignore for the time being.
for comp in ['LULUC_opt','Coupled_opt']:
    compCCS = comp.replace('opt','CCS')
    compNat = comp.replace('opt','Nat')
    # copy _CCS dictionay to _opt
    COMPS[comp] = deepcopy(COMPS[compCCS])
    for iscenario in range(nSCENARIOs): 
        scenario=SCENARIOs[iscenario]
        for iQ10 in range(nQ10s):
          Q10 = Q10_exps[iQ10]
          for iO3 in range(nO3s):
            O3 = Ozone_exps[iO3]
            for igcm in range(nGCMs):
                gcm=GCMs[igcm]
                #print(gcm)
                # create mask, CCS>Nat = 1; CCS==Nat = 0; CCS<Nat = -1 
                difference = ( ( COMPS[compCCS][scenario][Q10][O3[1]]['Land'][gcm]
                               + COMPS[compCCS][scenario][Q10][O3[1]]['CCS'][gcm] )
                             - ( COMPS[compNat][scenario][Q10][O3[1]]['Land'][gcm]
                               + COMPS[compNat][scenario][Q10][O3[1]]['CCS'][gcm] )  )
                flag_mask = (difference/np.abs(difference)).astype(int)
                flag_mask[difference==0.] = 0
                
                # Substitute Nat for CCS in land/CCS arrays where Nat is the prefered choice:
                for pool in land_pools:
                    COMPS[comp][scenario][Q10][O3[1]][pool][gcm][flag_mask==-1] = \
                            COMPS[compNat][scenario][Q10][O3[1]][pool][gcm][flag_mask==-1]
                
                # Recalculate the Total emission budget:
                COMPS[comp][scenario][Q10][O3[1]]['Total'][gcm] = (
                                         COMPS[comp][scenario][Q10][O3[1]]['Land'][gcm].sum()
                                       + COMPS[comp][scenario][Q10][O3[1]]['CCS'][gcm].sum()
                                       + COMPS[comp][scenario][Q10][O3[1]]['Atmos'][gcm]
                                       + COMPS[comp][scenario][Q10][O3[1]]['Ocean'][gcm] )


# Create Linear Sum comparison
for Lin_comp in ['Nat','CCS','opt']:
    # Start with the CH4 comparison:
    COMPS['Linear_'+Lin_comp] = deepcopy(COMPS['CH4'])
    for iscenario in range(nSCENARIOs): 
        scenario=SCENARIOs[iscenario]
        for iQ10 in range(nQ10s):
          Q10 = Q10_exps[iQ10]
          for iO3 in range(nO3s):
            O3 = Ozone_exps[iO3]
            for igcm in range(nGCMs):
                gcm=GCMs[igcm]
                #print(gcm)
                # Then add on the difference between the LULUC simulation and the control simulation
                for pool in pools:
                    COMPS['Linear_'+Lin_comp][scenario][Q10][O3[1]][pool][gcm] += (
                        COMPS['LULUC_'+Lin_comp][scenario][Q10][O3[1]][pool][gcm] 
                      - COMPS['CTL'][scenario][Q10][O3[1]][pool][gcm] )


#ipdb.set_trace()
################################################################################
delCStores     = { scenario: {} for scenario in SCENARIOs }
delCStores_unc = { scenario: {} for scenario in SCENARIOs }

delCStores_regions =  { scenario: { region : { Q10: {} for Q10 in Q10_exps } 
                                    for region in REGION_dict['Name'] } 
                        for scenario in SCENARIOs }
delCStores_regions_unc =  { scenario: { region : [] for region in REGION_dict['Name'] } 
                             for scenario in SCENARIOs }

SummaryFactExps = ['CTL','CH4','LULUC_opt','Linear_opt','Coupled_opt']
FactExps = ['CTL','CH4',
            'LULUC_CCS','LULUC_Nat','LULUC_opt',
            'Linear_CCS','Linear_Nat','Linear_opt',
            'Coupled_CCS','Coupled_Nat','Coupled_opt']

regional_DF_columns = [ 'CH4_mit', 'NatLand_mit', 'CCS_mit', 'LULUC_mit', 'Linear_mit' ]
nREGcolumns = len(regional_DF_columns)

StoreSinks = ['Atmos','Ocean','Land','CCS']
delCStores_SINKS     = { scenario: { } for scenario in SCENARIOs }
delCStores_SINKS_unc = { scenario: { exp:[] for exp in FactExps } for scenario in SCENARIOs }
delCStores_mitSINKS_unc = { scenario: { exp:[] for exp in FactExps } for scenario in SCENARIOs }

os.system('mkdir -p '+PLOT_DIR+'CSVoutput/')
print('Storing basic stats in: '+PLOT_DIR+'CSVoutput/Mitigation_Summary_Region.csv')
outf=open(PLOT_DIR+'CSVoutput/Mitigation_Summary_Global.csv','w')
outfreg=open(PLOT_DIR+'CSVoutput/Mitigation_Summary_Region.csv','w')

outfunc=open(PLOT_DIR+'CSVoutput/Mitigation_Summary_FullUnc.csv','w')

# Output in the following csv file used to generate data for Figure 11 in ESD paper
outfuncreg=open(PLOT_DIR+'CSVoutput/Mitigation_Summary_FullUnc_Reg.csv','w')

outfsink=open(PLOT_DIR+'CSVoutput/Mitigation_Summary_Sink_FullUnc.csv','w')
#
#ipdb.set_trace()
for iscenario in range(nSCENARIOs):
    scenario= SCENARIOs[iscenario]
    scen_list = []
    for iQ10 in range(nQ10s):
        Q10 = Q10_exps[iQ10]
        delCStores[scenario][Q10]={}
        delCStores_SINKS[scenario][Q10] = {}
        for iO3 in range(nO3s):
            O3 = Ozone_exps[iO3]
            outf.write('%22s'%(scenario+', '+Q10+', '+O3[1]+':\n'))
            outfreg.write('%22s'%(scenario+', '+Q10+', '+O3[1]+':\n'))
            outfreg.write('%15s'%('Region') + nREGcolumns*'%7s,' % tuple(regional_DF_columns)+'\n')
            
            # Global Emission budgets and mitigation potentials
            # Budgets:
            DataFrame_Tots = pd.concat( [pd.Series(COMPS[factexp][scenario][Q10][O3[1]]['Total'], index=GCMs)
                                       for factexp in FactExps ], axis=1 )
            DataFrame_Tots.columns = [ factexp+'_tot' for factexp in FactExps ]
            # Mitigattion (difference from control)
            DataFrame_Mits = deepcopy(DataFrame_Tots)
            for col in DataFrame_Tots.columns: DataFrame_Mits[col] -= DataFrame_Tots['CTL_tot']
            DataFrame_Mits.columns = [ factexp+'_mit' for factexp in FactExps ]
            
            DF = pd.concat([DataFrame_Tots,DataFrame_Mits],axis=1)

            delCStores[scenario][Q10][O3[1]] = DF
            DF.describe().to_csv(outf,float_format='%10.2f')
            scen_list.append(DF) 
            
            # Global Sinks:
            delCStores_SINKS[scenario][Q10][O3[1]]={}
            for exp in FactExps:
                #print(exp)
                delCStores_SINKS[scenario][Q10][O3[1]][exp] = pd.DataFrame( 
                            { sink: pd.Series({gcm:COMPS[exp][scenario][Q10][O3[1]][sink][gcm].sum() 
                                for gcm in GCMs } ) for sink in StoreSinks }  ) 
                # Append to full uncertatinty lists:
                delCStores_SINKS_unc[scenario][exp].append(delCStores_SINKS[scenario][Q10][O3[1]][exp])
                delCStores_mitSINKS_unc[scenario][exp].append( delCStores_SINKS[scenario][Q10][O3[1]][exp] 
                                                              -delCStores_SINKS[scenario][Q10][O3[1]]['CTL'] )

            # Land uptake on land points
            LULUC_mit_landpts = { gcm: COMPS['LULUC_opt'][scenario][Q10][O3[1]]['Land'][gcm]
                                     - COMPS['CTL'][scenario][Q10][O3[1]]['Land'][gcm] for gcm in GCMs}
            CCS_mit_landpts = {gcm:  COMPS['LULUC_opt'][scenario][Q10][O3[1]]['CCS'][gcm]
                                   - COMPS['CTL'][scenario][Q10][O3[1]]['CCS'][gcm] for gcm in GCMs}
            Total_Land_C_med = (  np.median( [LULUC_mit_landpts[gcm] for gcm in GCMs], axis=0 )
                                + np.median( [CCS_mit_landpts[gcm] for gcm in GCMs], axis=0 ) )
            plot_data = np.ma.masked_array( (Total_Land_C_med/AREA_1D)[grindex], mask = grindex.mask )/kg_to_Gt
            PTs.plot_map(plot_data,lons_2d,lats_2d,
                         DATA_RANGE=[0.1,10.], 
                         COLOURS=['#f6e8c3','#dfc27d','#d8b365','#bf812d','#8c510a'],INTERPOLATE_COLOURS=True,
                         TickLEVELS=[0.1,2,4,6,8,10],NLEVELS=50,SET_OVER='#8c510a',SET_UNDER='#f5f5f5',
                         PLOT_TITLE='Total Natural Land and BECCS carbon uptake by 2100',
                         CBAR_LABEL='kgC m$^{2}$',FONTSIZES=[12,12,14,18],RESOLUTION='c', MAP_TYPE='Mesh',
                         FILE_PLOT=PLOT_DIR+'Maps/TotalLULUC_Cuptake_Map_'+scenario+'_'+O3[1]+'_'+Q10+'.png',
                         iCLOSE='Y' )
                         
            # Regional Breakdown of the mitigaiton options
            for iregion in range(REGION_dict['Nregions']):
                region =REGION_dict['Name'][iregion]
                region_anthFrac=REGION_dict['AnthroFraction'][iregion]
                #region_map_index = REGION_dict['Index'][iregion]  #Not required as stored in oreder of map index
                region_mask=REGIONS_1D==(iregion+1)
                if region=='Global': region_mask[:]=True
                if region=='International Transportation': region_mask[:]=False
                regional_CH4_mit = DF['CH4_mit']*region_anthFrac
                regional_Land_mit = pd.Series({gcm:LULUC_mit_landpts[gcm][region_mask].sum() for gcm in GCMs} ) 
                regional_CCS_mit = pd.Series({gcm:CCS_mit_landpts[gcm][region_mask].sum() for gcm in GCMs } ) 
                regional_LULUC_mit = regional_CCS_mit+regional_Land_mit
                regional_Linear_mit = regional_CH4_mit+regional_LULUC_mit
                DFreg = pd.concat([regional_CH4_mit,regional_Land_mit,regional_CCS_mit,
                                   regional_LULUC_mit,regional_Linear_mit], axis=1) 
                DFreg.columns = regional_DF_columns
                delCStores_regions[scenario][region][Q10][O3[1]] = DFreg
                delCStores_regions_unc[scenario][region].append(
                                delCStores_regions[scenario][region][Q10][O3[1]] )
                outfreg.write('%15s: '%(region))
                DFreg.describe()[5:6].to_csv(outfreg,header=False,float_format='%10.2f',index=False)
                #print(region, DFreg.describe()[5:6])
    
    # compend Global Carbon Stores, full uncertainty to DataFrame and output to csv
    delCStores_unc[scenario] = pd.concat(scen_list)
    outfunc.write(scenario+': \n')
    outfunc.write('Global: \n')
    delCStores_unc[scenario].describe().to_csv(outfunc,float_format='%10.2f')
    
    # Compend the global carbon stores by pool to DataFrame and output to csv 
    delCStores_SINKS_unc[scenario] = pd.concat( { EXP: pd.concat(delCStores_SINKS_unc[scenario][EXP]) 
                                                  for EXP in FactExps }, axis=1)
    delCStores_mitSINKS_unc[scenario] = pd.concat( { EXP: pd.concat(delCStores_mitSINKS_unc[scenario][EXP]) 
                                                     for EXP in FactExps }, axis=1)
    outfsink.write(scenario+': \n')
    outfsink.write('Total Sink: \n')
    for EXP in FactExps: 
        outfsink.write('%16s'%(EXP+': \n'))
        delCStores_SINKS_unc[scenario][EXP].describe().to_csv(outfsink,float_format='%10.2f')
    outfsink.write('Mitigation Potential: \n')
    for EXP in FactExps: 
        outfsink.write('%16s'%(EXP+': \n'))
        delCStores_mitSINKS_unc[scenario][EXP].describe().to_csv(outfsink,float_format='%10.2f')

    # Compend the regional breakdown into DataFrame: 
    delCStores_regions_unc[scenario] = pd.concat({ region: pd.concat(delCStores_regions_unc[scenario][region])
                                                    for region in REGION_dict['Name'] }, axis=1 )
    # Output regional breakdown to csv file 
    outfuncreg.write(scenario+': \nMedian of GCMs:\n')
    outfuncreg.write('%25s '%('Region,')+nREGcolumns*'%9s,' % tuple(regional_DF_columns)+'\n')
    for region in REGION_dict['Name']: 
        outfuncreg.write('%25s '%(region+','))
        delCStores_regions_unc[scenario][region].describe()[5:6].to_csv(outfuncreg,float_format='%10.2f',
                                                                         header=False,index=False)
    outfuncreg.write('\n\n\n25% of GCMs: \n')
    outfuncreg.write('%25s '%('Region,')+nREGcolumns*'%9s,' % tuple(regional_DF_columns)+'\n')
    for region in REGION_dict['Name']: 
        outfuncreg.write('%25s '%(region+','))
        delCStores_regions_unc[scenario][region].describe()[4:5].to_csv(outfuncreg,float_format='%10.2f',
                                                                         header=False,index=False)
    outfuncreg.write('\n\n\n75% of GCMs: \n')
    outfuncreg.write('%25s '%('Region,')+nREGcolumns*'%9s,' % tuple(regional_DF_columns)+'\n')
    for region in REGION_dict['Name']: 
        outfuncreg.write('%25s '%(region+','))
        delCStores_regions_unc[scenario][region].describe()[6:7].to_csv(outfuncreg,float_format='%10.2f',
                                                                         header=False,index=False)

outfunc.close()
outfuncreg.close()
outfsink.close()
outf.close()
outfreg.close()
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
    CTLcolor = '#666666' 
    bnwcolor = '#e7298a'

###################################################################################################
# Used for Figure 9 in ESD paper (BECCS scale factors = 1 & 3)
# 2 Bar plots:
# 1a - AFFEB = Delta C stock from Present Day 
# 1b - Mitigation potential = AFFEB - Control_AFFEB
# Open fig with 2 wide rows for plotting
if True:
    # Open figure
    fig,axes = plt.subplots(ncols=1,nrows=2,figsize=[8,6])
    fig.subplots_adjust(top=0.94,left=0.07,right=0.9998)
    # Y limits for plots, manually selected
    AFFEB_limit = [-200,800]
    MP_limit = [0,350]
    # Plot variables, colours and labels
    plotvar_list_tots = ['CTL_tot', 'CH4_tot', 'LULUC_opt_tot', 'Coupled_opt_tot', 'Linear_opt_tot']  #
    plotvar_list_mits = ['CTL_mit', 'CH4_mit', 'LULUC_opt_mit', 'Coupled_opt_mit', 'Linear_opt_mit']  #
    plotvar_cols = [CTLcolor,CH4color,totLULUCcolor,totcolor,totcolor_lin]# 
    legend_names = ['Control', 'CH$_4$ Mitigation', 'Land Based Mitigation', 'Coupled', 'Linear']  #
    plotvar_list_pools= ['CTL', 'CH4', 'LULUC_opt', 'Coupled_opt', 'Linear_opt' ]  
    pool_list    = ['Atmos', 'Ocean', 'Land', 'CCS']
    poollabellist= ['Atmosphere', 'Ocean', 'Land', 'BECCS']
    pool_colours = [Atmoscolor,Oceancolor,Landcolor,CCScolor]
    npools = len(pool_list)
    # Bar spacing options
    scenario_space = 0.8 # how much space all the bars should take up for a scenario 
    nbars = len(plotvar_list_tots)  # Number of bars per scenario (i.e. Veg, Soil, Amos)
    bar_fraction=0.75 # fraction of bar space to fill with bar, 1 means bars touch
    bar_space = (scenario_space/nbars) # total space for a bar
    bar_width = (bar_space)*(bar_fraction) # width of bar
    bar_gap   = (bar_space)*(1.-bar_fraction) # gap between bars
    bar_positions = [ -(scenario_space/2.)+(i*bar_space)+(bar_gap/2.) for i in range(nbars) ]
    pool_bar_frac = 0.2 # fraction of bar taken up by pool breakdown
    pool_bar_width = pool_bar_frac *bar_width # fraction of bar taken up by pool breakdown
    pool_pos = bar_width-pool_bar_width
    pool_scat_inc = (pool_bar_width/npools)
    scatter_pos = pool_pos*0.5  # position of the scatter points for total bar
    # Loop and plot
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        scen_cen_pos = iscenario+0.5
        tot_bar_list=[]  # Append the bar objects to list for legend
        mit_bar_list=[]  # Append the bar objects to list for legend
        pool_bar_list=[]  # Append the bar objects to list for legend
        for ibar in range(nbars):
            totbar     = plotvar_list_tots[ibar]
            mitbar     = plotvar_list_mits[ibar]
            xpos       = scen_cen_pos + bar_positions[ibar] 
            barcolour  = plotvar_cols[ibar]
            # Plot AFFEB on the first axis:
            tot_plotdata = delCStores_unc[scenario][totbar]
            tot_bar_list.append(axes[0].bar(xpos,np.median(tot_plotdata),color=barcolour,width=bar_width,))
            axes[0].plot([xpos+scatter_pos for i in range(len(tot_plotdata))],tot_plotdata,c='k',ls='',marker='.')
            # Plot MP on the scedond axis:
            mit_plotdata = delCStores_unc[scenario][mitbar]
            mit_bar_list.append(axes[1].bar(xpos,np.median(mit_plotdata),color=barcolour,width=bar_width,))
            axes[1].plot([xpos+scatter_pos for i in range(len(mit_plotdata))],mit_plotdata,c='k',ls='',marker='.')
            # Overplot pools on the axes:
            poolbar    = plotvar_list_pools[ibar]
            for poolDF,ax in zip([delCStores_SINKS_unc,delCStores_mitSINKS_unc],axes):
                plotDF = poolDF[scenario][poolbar]
                negative_min,positive_max = 0,0  # vars for stacking bars
                if len(pool_bar_list)<npools+1:
                    pool_bar_list=[]  # Append the bar objects to list for legend
                else:
                    pool_bar_list=bar_list[:npools+1]
                for ipool in range(npools):
                    pool = pool_list[ipool]
                    colour = pool_colours[ipool]
                    plotdata = plotDF[pool]
                    median = np.median(plotdata)
                    if median>0:
                        pool_bar_list.append(ax.bar(xpos+pool_pos, median, 
                                            bottom=positive_max,color=colour,width=pool_bar_width,))
                        positive_max += median
                    elif median<0:
                        pool_bar_list.append(ax.bar(xpos+pool_pos, median, 
                                            bottom=negative_min,color=colour,width=pool_bar_width,))
                        negative_min += median
    # Axis options
    for iax in range(len(axes)):
        ax=axes[iax]
        ax.set_xlim([0,nSCENARIOs])
        ax.set_xticks(np.arange(0.5,nSCENARIOs,1.))
        ax.set_xticklabels(SCENARIO_names,fontsize=14)
        ax.tick_params(axis='y',labelsize=14)
        ax.plot([0,nSCENARIOs],[0,0],c='k',lw=1)
        # labels:
        ax.text(0.01,1.05,'('+ALPHABET[iax+Alpha_Offset]+')', 
                transform=ax.transAxes, fontsize=18, fontweight='bold')
    # AFFEB axis options:
    axes[0].set_ylabel('AFFEB 2015-2100 (GtC)',fontsize=15)
    axes[0].set_ylim(AFFEB_limit)
    # MP axis options
    axes[1].set_ylim(MP_limit)
    axes[1].set_ylabel('Mitigation Potential (GtC)',fontsize=15)

    # Legends
    if Alpha_Offset == 0:
        fig.legend(tot_bar_list,legend_names, loc='upper right', ncol=nbars+1, 
                fontsize=12, columnspacing=1.5, handlelength=1, handletextpad = 0.5 )
    else:
        fig.legend(pool_bar_list,poollabellist, loc='upper right', ncol=npools+1,
                fontsize=12, columnspacing=1.5, handlelength=1, handletextpad = 0.5 )

    #Save figs:
    fig.savefig(PLOT_DIR+'BarCharts/Figure_1.png',bbox_inches='tight')  # Store as png
    fig.savefig(PLOT_DIR+'BarCharts/Figure_1.eps',bbox_inches='tight')  # Store as eps
    fig.savefig(PLOT_DIR+'BarCharts/Figure_1.pdf',bbox_inches='tight')  # Store as pdf
    #plt.show()
    plt.close()

#ipdb.set_trace()

############################################################################################################
###################################################################################################
# Figure for synthesis paper.
# 2 Bar plots:
# 1a - AFFEB = Delta C stock from Present Day 
# 1b - Mitigation potential = AFFEB - Control_AFFEB
# Open fig with 2 wide rows for plotting
if True:
    # Open figure
    fig,axes = plt.subplots(ncols=1,nrows=2,figsize=[8,6])
    fig.subplots_adjust(top=0.95,left=0.07,right=0.92)
    
    # Plot variables, colours and labels
    plotvar_list_tots = ['CTL_tot', 'CH4_tot', 'LULUC_opt_tot', 'Coupled_opt_tot', 'Linear_opt_tot']  #
    plotvar_list_mits = ['CTL_mit', 'CH4_mit', 'LULUC_opt_mit', 'Coupled_opt_mit', 'Linear_opt_mit']  #
    plotvar_cols = [CTLcolor,CH4color,totLULUCcolor,totcolor,totcolor_lin]# 
    legend_names = ['Control', 'CH$_4$ Mitigaiton', 'Land Based Mitigation', 'Coupled', 'Linear']  #
    plotvar_list_pools= ['CTL', 'CH4', 'LULUC_opt', 'Coupled_opt', 'Linear_opt' ]  
    pool_list    = ['Atmos', 'Ocean', 'Land', 'CCS']
    poollabellist= ['Atmosphere', 'Ocean', 'Land', 'BECCS']
    pool_colours = [Atmoscolor,Oceancolor,Landcolor,CCScolor]
    npools = len(pool_list)
    
    # Bar spacing options
    scenario_space = 0.9 # how much space all the bars should take up for a scenario 
    nbars = len(plotvar_list_tots)  # Number of bars per scenario (i.e. Veg, Soil, Amos)
    bar_fraction=0.8  # fraction of bar space to fill with bar, 1 means bars touch
    bar_space = (scenario_space/nbars) # total space for a bar
    bar_width = (bar_space)*(bar_fraction) # width of bar
    bar_gap   = (bar_space)*(1.-bar_fraction) # gap between bars
    bar_positions = [ -(scenario_space/2.)+(i*bar_space)+(bar_gap/2.) for i in range(nbars) ]
    pool_bar_frac = 0.4 # fraction of bar taken up by pool breakdown
    pool_bar_width = pool_bar_frac *bar_width # fraction of bar taken up by pool breakdown
    pool_pos = bar_width-pool_bar_width
    pool_scat_inc = (pool_bar_width/npools)
    scatter_pos = pool_pos*0.5  # position of the scatter points for total bar
    # Loop and plot
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        scen_cen_pos = iscenario+0.5
        tot_bar_list=[]  # Append the bar objects to list for legend
        mit_bar_list=[]  # Append the bar objects to list for legend
        pool_bar_list=[]  # Append the bar objects to list for legend
        for ibar in range(nbars):
            totbar     = plotvar_list_tots[ibar]
            mitbar     = plotvar_list_mits[ibar]
            xpos       = scen_cen_pos + bar_positions[ibar] 
            barcolour  = plotvar_cols[ibar]
            # Plot AFFEB on the first axis:
            tot_plotdata = delCStores_unc[scenario][totbar]
            tot_bar_list.append(axes[0].bar(xpos,np.median(tot_plotdata),color=barcolour,width=bar_width,))
            axes[0].plot([xpos+scatter_pos for i in range(len(tot_plotdata))],tot_plotdata,c='k',ls='',marker='.')
            
            # Plot MP on the scedond axis:
            mit_plotdata = delCStores_unc[scenario][mitbar]
            mit_bar_list.append(axes[1].bar(xpos,np.median(mit_plotdata),color=barcolour,width=bar_width,))
            axes[1].plot([xpos+scatter_pos for i in range(len(mit_plotdata))],mit_plotdata,c='k',ls='',marker='.')

            # Overplot pools on the axes:
            poolbar    = plotvar_list_pools[ibar]
            for poolDF,ax in zip([delCStores_SINKS_unc,delCStores_mitSINKS_unc],axes):
                plotDF = poolDF[scenario][poolbar]
                negative_min,positive_max = 0,0  # vars for stacking bars
                if len(pool_bar_list)<npools+1:
                    pool_bar_list=[]  # Append the bar objects to list for legend
                else:
                    pool_bar_list=bar_list[:npools+1]
                for ipool in range(npools):
                    pool = pool_list[ipool]
                    colour = pool_colours[ipool]
                    plotdata = plotDF[pool]
                    median = np.median(plotdata)
                    if median>0:
                        pool_bar_list.append(ax.bar(xpos+pool_pos, median, 
                                            bottom=positive_max,color=colour,width=pool_bar_width,))
                        ax.scatter([xpos+pool_pos+pool_scat_inc/2.+ipool*pool_scat_inc for i in range(len(plotdata))],
                                    plotdata+positive_max,c=colour,edgecolor='k',marker='.',zorder=10,s=50)
                        positive_max += median
                    elif median<0:
                        pool_bar_list.append(ax.bar(xpos+pool_pos, median, 
                                            bottom=negative_min,color=colour,width=pool_bar_width,))
                        ax.scatter([xpos+pool_pos+pool_scat_inc/2.+ipool*pool_scat_inc for i in range(len(plotdata))],
                                    plotdata+negative_min,c=colour,edgecolor='k',marker='.',zorder=10,s=50)
                        negative_min += median
#if True:
    # Axis options
    for iax in range(len(axes)):
        ax=axes[iax]
        ax.set_xlim([0,nSCENARIOs])
        ax.set_xticks(np.arange(0.5,nSCENARIOs,1.))
        ax.set_xticklabels(SCENARIO_names,fontsize=14)
        ax.tick_params(axis='y',labelsize=14)
        ax.plot([0,nSCENARIOs],[0,0],c='k',lw=1)
        # labels:
        ax.text(0.02,1.03,'('+ALPHABET[iax+Alpha_Offset]+')', 
                transform=ax.transAxes, fontsize=18, fontweight='bold')
    # AFFEB axis options:
    axes[0].set_ylabel('AFFEB 2015-2100 (GtC)',fontsize=15)
    # MP axis options
    axes[1].set_ylim([0,axes[1].get_ylim()[1]])
    axes[1].set_ylabel('Mitigation Potential (GtC)',fontsize=15)

    # Legends 
    axes[0].legend(tot_bar_list,legend_names, loc='upper left', ncol=nbars+1, fontsize=10)
    axes[1].legend(pool_bar_list,poollabellist, loc='upper left', ncol=npools+1, fontsize=10)
    #Save figs:
    fig.savefig(PLOT_DIR+'BarCharts/Figure_1_withDOTS.png',bbox_inches='tight')  # Store as png
    fig.savefig(PLOT_DIR+'BarCharts/Figure_1_withDOTS.eps',bbox_inches='tight')  # Store as eps
    fig.savefig(PLOT_DIR+'BarCharts/Figure_1_withDOTS.pdf',bbox_inches='tight')  # Store as pdf
    #plt.show()
    plt.close()

#ipdb.set_trace()

############################################################################################################
# Used for Figure 12 in ESD paper (BECCS scale factors = 1, 2 & 3)
# Plot Regional Mitigation Potential as stacked bar plot
if True:
    plotvar_list = ['CH4_mit', 'NatLand_mit',         'CCS_mit']
    plot_colours = [CH4color,Landcolor,CCScolor]  #['#7570b3', '#66a61e',  '#e6ab02']
    total_colour = bnwcolor
    legend_names = ['CH4',     'Natural Land Uptake', 'BECCS'] 
    nstackbars = len(plotvar_list)  # Number of stacked bars 
    region_space = 0.9 # how much space all the bars should take up for a region 
    nbars = 1         # Number of bars per region 
    bar_fraction=0.9  # fraction of bar space to fill with bar, 1 means bars touch
    bar_space = (region_space/nbars) # total space for a bar
    bar_width = (bar_space)*(bar_fraction) # width of bar
    bar_gap   = (bar_space)*(1.-bar_fraction) # gap between bars
    bar_positions = [ -(region_space/2.)+(i*bar_space)+(bar_gap/2.) for i in range(nbars) ]

for iscenario in range(nSCENARIOs):
#    if True:
    scenario=SCENARIOs[iscenario]
    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(10,5))
    fig.subplots_adjust(top=0.9,left=0.1,right=0.99,bottom=0.1)
    # Loop over regions:
    bar_list=[]  # Append the bar objects to list for legend
    for iregion in range(REGION_dict['Nregions']-2):
        region =REGION_dict['Name'][iregion]
        scen_cen_pos=0.5+iregion
        if len(bar_list)<nstackbars+1:
            bar_list=[]  # Append the bar objects to list for legend
        else:
            bar_list=bar_list[:nstackbars+1]
        barbase = 0.
        medians    = np.zeros((nstackbars))
        for ibar in range(nstackbars):
            bar        = plotvar_list[ibar]
            xpos       = scen_cen_pos + bar_positions[0] 
            barcolour  = plot_colours[ibar]
            plotdata   = delCStores_regions_unc[scenario][region][bar]
            median     = np.median(plotdata)
            medians[ibar] = median
            if median>1e-1:
                median_tmp = medians[0:ibar]
                barbase    = median_tmp[median_tmp > 0].sum()
                bar_list.append(ax.bar(xpos,np.median(plotdata),color=barcolour,width=bar_width,bottom=barbase))
                #ax.scatter([xpos+bar_width/2 for i in range(len(plotdata))],plotdata+barbase,
                #              c=barcolour,edgecolor='k',marker='.',zorder=10,s=100)
            elif median<-1e-1:
                median_tmp = medians[0:ibar]
                barbase    = median_tmp[median_tmp < 0].sum()
                bar_list.append(ax.bar(xpos,np.median(plotdata),color=barcolour,width=bar_width,bottom=barbase))
            #print('Bar',iregion,ibar,median,barbase)

        xpos         = scen_cen_pos # + (bar_space/2.)
        mit_plotdata = delCStores_regions_unc[scenario][region]['Linear_mit']
        bar_list.append(box_and_whisker(mit_plotdata,xpos,total_colour,ax,boxwidth=bar_width*0.5)[0])

    ax.set_xlim([0,REGION_dict['Nregions']-2])
    ax.plot(ax.get_xlim(),[0,0],c='k',lw=2)
    ax.set_xticks(np.arange(0.5,REGION_dict['Nregions']-2,1.))
    ax.set_xticklabels(REGION_dict['ShortName'][:-2],fontsize=15,rotation=60)
    ax.tick_params(axis='y',labelsize=15)
    ax.text(0.03,0.91,SCENARIO_names[iscenario],transform=ax.transAxes, 
            fontsize=16)
    #ax.set_ylim([0,ax.get_ylim()[1]])
    ax.set_ylim([-10,45])
    #ax.set_ylabel('Anthropogenic Fossil Fuel Emission Budget\n 2015-2100 (GtC)',fontsize=17)
    ax.set_ylabel('AFFEB 2015-2100 (GtC)',fontsize=17)

    fig.legend(bar_list,legend_names+['Total'], 
            loc='upper center',ncol=nstackbars+1 ,fontsize=15)
    fig.savefig(PLOT_DIR+'BarCharts/Mitigation_Potential_Pools_Region_'+scenario+'.png', bbox_inches='tight')  # Store as png
    fig.savefig(PLOT_DIR+'BarCharts/Mitigation_Potential_Pools_Region_'+scenario+'.pdf', bbox_inches='tight')  # Store as pdf
    fig.savefig(PLOT_DIR+'BarCharts/Mitigation_Potential_Pools_Region_'+scenario+'.eps', bbox_inches='tight')  # Store as eps
    #plt.show()
    plt.close()

#quit()

#ipdb.set_trace()
###################################################################################################
# Bar plots of Mitigation Potential with Carbon Sinks
for iscenario in range(nSCENARIOs):
    scenario=SCENARIOs[iscenario]
    # Open fig with 1 wide row for plotting each sceanrio to seperate file
    fig,axes = plt.subplots(ncols=1,nrows=1,figsize=[10,5])
    fig.subplots_adjust(top=0.9,left=0.1,right=0.97)
    # Upper plot is for absolute AFFEB to 2100
    plotvar_list = ['CH4', 'LULUC_CCS', 'LULUC_Nat', 'LULUC_opt', 'Coupled_CCS', 'Coupled_Nat', 'Coupled_opt', 'Linear_opt' ]  
    barname_list = [ varname.replace('CH4','CH$_4$\nMitigation').replace('LULUC','Land Based\nMitigation').replace(
                                    '_CCS','\n(CCS)').replace('_Nat','\n(Natural Land)').replace('_opt','\n(Optimised)')
                          for varname in plotvar_list ]
    pool_list    = ['Atmos', 'Ocean', 'Land', 'CCS']
    pool_colours = [Atmoscolor,Oceancolor,Landcolor,CCScolor]
    total_colour = bnwcolor
    nexps = len(plotvar_list)
    npools = len(pool_list)
    exp_space = 0.8 # how much space all the bars should take up for a scenario 
    nbars = 2
    bar_fraction=0.9  # fraction of bar space to fill with bar, 1 means bars touch
    bar_space = (exp_space/nbars) # total space for a bar
    bar_width = (bar_space)*(bar_fraction) # width of bar
    bar_gap   = (bar_space)*(1.-bar_fraction) # gap between bars
    bar_positions = [ -(exp_space/2.)+(i*bar_space)+(bar_gap/2.) for i in range(nbars) ]
    bar_list=[]  # Append the bar objects to list for legend
    for iexp in range(nexps):
        exp = plotvar_list[iexp]
        exp_cen_pos = iexp+0.5
        plotDF = delCStores_mitSINKS_unc[scenario][exp]
        xpos   = exp_cen_pos + bar_positions[0] 
        negative_min,positive_max = 0,0
        if len(bar_list)<npools+1:
            bar_list=[]  # Append the bar objects to list for legend
        else:
            bar_list=bar_list[:npools+1]
        for ipool in range(npools):
            pool = pool_list[ipool]
            colour = pool_colours[ipool]
            plotdata = plotDF[pool]
            median = np.median(plotdata)
            if median>0:
                bar_list.append(axes.bar(xpos,median,bottom=positive_max,color=colour,width=bar_width,))
                axes.scatter([xpos+bar_width/4+((float(ipool)/npools)*(bar_width/2)) for i in range(len(plotdata))],
                                plotdata+positive_max,c=colour,edgecolor='k',marker='.',zorder=10,s=50)
                positive_max += median
            elif median<0:
                axes.bar(xpos,median,bottom=negative_min,color=colour,width=bar_width,)
                axes.scatter([xpos+bar_width/4+((float(ipool)/npools)*(bar_width/2)) for i in range(len(plotdata))],
                                plotdata+negative_min,c=colour,edgecolor='k',marker='.',zorder=10,s=50)
                negative_min += median
        xpos         = exp_cen_pos + (bar_space/2.)
        mit_plotdata = delCStores_unc[scenario][exp+'_mit']
        bar_list.append(box_and_whisker(mit_plotdata,xpos,total_colour,axes,boxwidth=bar_width*0.5)[0])

    axes.plot([0,nexps],[0,0],c='k',lw=1)
    axes.set_xlim([0,nexps])
    axes.set_xticks(np.arange(0.5,nexps,1.))
    axes.set_xticklabels(barname_list,fontsize=11)
    axes.tick_params(axis='y',labelsize=12)
    axes.set_ylabel('Mitigation Potential 2015-2100 (GtC)',fontsize=12)
    fig.legend(bar_list,pool_list+['Total'], loc='upper center', ncol=npools+1, fontsize=12)
    fig.savefig(PLOT_DIR+'BarCharts/MitigationPotential_Global_Pool_'+scenario+'.png',bbox_inches='tight')  # Store as png
    fig.savefig(PLOT_DIR+'BarCharts/MitigationPotential_Global_Pool_'+scenario+'.eps',bbox_inches='tight')  # Store as eps
    fig.savefig(PLOT_DIR+'BarCharts/MitigationPotential_Global_Pool_'+scenario+'.pdf',bbox_inches='tight')  # Store as pdf
    plt.close()
    #plt.show()


###################################################################################################
# Used for Figure 8 in ESD paper (BECCS scale factor = 1)
# Bar plots of Delta Stocks from Present Day with Carbon Sinks
for iscenario in range(nSCENARIOs):
    scenario=SCENARIOs[iscenario]
    # Open fig with 1 wide row for plotting each sceanrio to seperate file
    fig,axes = plt.subplots(ncols=1,nrows=1,figsize=[10,5])
    fig.subplots_adjust(top=0.9,left=0.1,right=0.97)
    # Upper plot is for absolute AFFEB to 2100
    plotvar_list = ['CTL', 'CH4', 'LULUC_CCS', 'LULUC_Nat', 'LULUC_opt', 
                    'Coupled_CCS', 'Coupled_Nat', 'Coupled_opt', 'Linear_opt' ]  
    barname_list = [ varname.replace('CH4','CH$_4$\nMitigation').replace('LULUC','Land Based\nMitigation').replace(
                  'CTL','Control').replace('_CCS','\n(CCS)').replace('_Nat','\n(Natural Land)').replace('_opt','\n(Optimised)')
                          for varname in plotvar_list ]
    pool_list    = ['Atmos', 'Ocean', 'Land', 'CCS']
    pool_colours = [Atmoscolor,Oceancolor,Landcolor,CCScolor]
    total_colour = bnwcolor 
    nexps = len(plotvar_list)
    npools = len(pool_list)
    exp_space = 0.8 # how much space all the bars should take up for a scenario 
    nbars = 2
    bar_fraction=0.9  # fraction of bar space to fill with bar, 1 means bars touch
    bar_space = (exp_space/nbars) # total space for a bar
    bar_width = (bar_space)*(bar_fraction) # width of bar
    bar_gap   = (bar_space)*(1.-bar_fraction) # gap between bars
    bar_positions = [ -(exp_space/2.)+(i*bar_space)+(bar_gap/2.) for i in range(nbars) ]
    bar_list=[]  # Append the bar objects to list for legend
    for iexp in range(nexps):
        exp = plotvar_list[iexp]
        exp_cen_pos = iexp+0.5
        plotDF = delCStores_SINKS_unc[scenario][exp]
        xpos       = exp_cen_pos + bar_positions[0] 
        negative_min,positive_max = 0,0
        if len(bar_list)<npools+1:
            bar_list=[]  # Append the bar objects to list for legend
        else:
            bar_list=bar_list[:npools+1]
        for ipool in range(npools):
            pool = pool_list[ipool]
            colour = pool_colours[ipool]
            plotdata = plotDF[pool]
            median = np.median(plotdata)
            if median>0:
                bar_list.append(axes.bar(xpos,median,bottom=positive_max,color=colour,width=bar_width,))
                axes.scatter([xpos+bar_width/4+((float(ipool)/npools)*(bar_width/2)) for i in range(len(plotdata))],
                                plotdata+positive_max,c=colour,edgecolor='k',marker='.',zorder=10,s=50)
                positive_max += median
            elif median<0:
                bar_list.append(axes.bar(xpos,median,bottom=negative_min,color=colour,width=bar_width,))
                axes.scatter([xpos+bar_width/4+((float(ipool)/npools)*(bar_width/2)) for i in range(len(plotdata))],
                                plotdata+negative_min,c=colour,edgecolor='k',marker='.',zorder=10,s=50)
                negative_min += median
        xpos         = exp_cen_pos + (bar_space/2.)
        tot_plotdata = delCStores_unc[scenario][exp+'_tot']
        bar_list.append(box_and_whisker(tot_plotdata,xpos,total_colour,axes,boxwidth=bar_width*0.5)[0])
    axes.plot([0,nexps],[0,0],c='k',lw=1)
    axes.set_xlim([0,nexps])
    axes.set_xticks(np.arange(0.5,nexps,1.))
    axes.set_xticklabels(barname_list,fontsize=12)
    axes.tick_params(axis='y',labelsize=12)
    axes.set_ylabel('Anthropogenic Fossil Fuel Emission Budget\n 2015-2100 (GtC)',fontsize=13)
    fig.legend(bar_list,pool_list+['Total'], loc='upper center', ncol=npools+1, fontsize=12)
    fig.savefig(PLOT_DIR+'BarCharts/AFFEB_Global_Pool_'+scenario+'.png',bbox_inches='tight')  # Store as png
    fig.savefig(PLOT_DIR+'BarCharts/AFFEB_Global_Pool_'+scenario+'.pdf',bbox_inches='tight')  # Store as png
    fig.savefig(PLOT_DIR+'BarCharts/AFFEB_Global_Pool_'+scenario+'.eps',bbox_inches='tight')  # Store as png
    plt.close()
    #plt.show()


#ipdb.set_trace()
