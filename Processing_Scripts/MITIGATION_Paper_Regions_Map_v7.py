#!/usr/bin/env python
# coding: utf-8

# In[1]:

import numpy as np
import os
import sys

from matplotlib import pylab, rc, rcParams, cm
import matplotlib.pyplot as plt
import matplotlib.colors as col
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import *


# In[2]:

# Configure for python (.py) or Jupyter notebook

if '-platform' in sys.argv:
    ARGLOC         = sys.argv.index('-platform')
    TEMP           = sys.argv.pop(ARGLOC)
    PLATFORM       = str(sys.argv.pop(ARGLOC))

    DEBUG          = sys.argv[1]
    iDISPLAY       = sys.argv[2]
    iFIX           = sys.argv[3]
    sTEMPER_OPT    = sys.argv[4]
    sBECCS_FACTOR  = sys.argv[5]
    sDATE          = sys.argv[6]
else:
    PLATFORM       = 'JUPYTER notebook'
    DEBUG          = 'Y'
    iDISPLAY       = 'Y'
    iFIX           = 'N'
    sTEMPER_OPT    = '1.5'
    sBECCS_FACTOR  = '3.0'
    sDATE          = '20191007'

print('Platform = '+PLATFORM)


# In[3]:

# Get data for map of IMAGE Regions

IMAGE_DIR      = '/prj/CLIFFTOP/CLUES/'
IMAGE_FILE     = IMAGE_DIR+'GREGION_IMAGE_5arcmin.asc'
IMAGE_FID      = open(IMAGE_FILE,'r')
IMAGE_LINES    = IMAGE_FID.readlines()
IMAGE_FID.close()

nCOLUMNS       = int(IMAGE_LINES.pop(0).split()[-1])
nROWS          = int(IMAGE_LINES.pop(0).split()[-1])
LON_MIN        = float(IMAGE_LINES.pop(0).split()[-1])
LAT_MIN        = float(IMAGE_LINES.pop(0).split()[-1])
RESOL          = float(IMAGE_LINES.pop(0).split()[-1])
FILL_VALUE     = int(IMAGE_LINES.pop(0).split()[-1])

LATS_1D        = (np.arange(nROWS)*RESOL)+LAT_MIN
LONGS_1D       = (np.arange(nCOLUMNS)*RESOL)+LON_MIN
LONGS_2D, LATS_2D = np.meshgrid(LONGS_1D, LATS_1D)
IMAGE_DATA     = np.zeros_like(LATS_2D)

for iLINE in range(nROWS):
    SPLIT          = IMAGE_LINES[iLINE].split()
    DATA_LINE      = np.array([ int(NUMBER) for NUMBER in SPLIT ])
    IMAGE_DATA[iLINE,:] = DATA_LINE

IMAGE_REG      = np.ma.masked_equal(IMAGE_DATA,FILL_VALUE)
DATA_IMAGE     = IMAGE_REG[::-1,:]


# In[4]:

# Define map grid
LAT_START      =  -90.0
LAT_END        =   90.0
DEL_LAT        =   30.0
DLAT           = LAT_END-LAT_START
LONG_START     = -180.0
LONG_END       =  180.0
DEL_LONG       =   30.0
DLONG          = LONG_END-LONG_START
if LONG_END == 360.0:
    LONG_PLOTS     = -180.0
    LONG_PLOTE     =  180.0
else:
    LONG_PLOTS     = LONG_START
    LONG_PLOTE     = LONG_END
    
LONG_GRID      = np.arange(LONG_PLOTS,LONG_PLOTE+DEL_LONG,DEL_LONG)
LAT_GRID       = np.arange(LAT_START,LAT_END+DEL_LAT,DEL_LAT)


# In[5]:

#Meta data for charts
sDATA_OPT      = '0'
META_CHART     = [ [ '0','MitigationOptions_BECCS'+sBECCS_FACTOR,['Reg']] ]

META_TEMP      = ['1.5','2.0']
iTEMPER_OPT    = META_TEMP.index(sTEMPER_OPT)

# IMAGE regions
META_IMAGE     = [ [ '1','Canada'          ,[  -90.0,  55.0,  10.0, 35.0], [ 10,  5, -5 ], [  30, 10, -10 ] ], \
                   [ '2','USA'             ,[ -140.0,  23.0,  10.0, 35.0], [ 20,  5,  0 ], [  80, 20,   0 ] ], \
                   [ '5','Mexico'          ,[ -120.0,  10.0,  10.0, 35.0], [  5,  5,  0 ], [  20, 10,   0 ] ], \
                   [ '4','C. America'      ,[  -60.0,  23.0,  10.0, 35.0], [  5,  5,  0 ], [  20, 10,   0 ] ], \
                   [ '5','Brazil'          ,[  -50.0, -55.0,  10.0, 35.0], [ 15,  5,  0 ], [  40, 10,   0 ] ], \
                   [ '6','Rest S. America' ,[  -90.0, -55.0,  10.0, 35.0], [ 15,  5,  0 ], [  40, 10,   0 ] ], \
                   [ '7','N. Africa'       ,[  -27.0,  23.0,  10.0, 35.0], [  5,  5,  0 ], [  20, 10,   0 ] ], \
                   [ '8','W. Africa'       ,[  -27.0, -20.0,  10.0, 35.0], [ 25,  5,  0 ], [  80, 20,   0 ] ], \
                   [ '9','E. Africa'       ,[   40.0,   0.0,  10.0, 35.0], [ 10,  5,  0 ], [  30, 10,   0 ] ], \
                   ['10','S. Africa'       ,[   35.0, -55.0,  10.0, 35.0], [ 15,  5, -5 ], [  60, 20, -10 ] ], \
                   ['11','W. Europe'       ,[   -5.0,  55.0,  10.0, 35.0], [ 15,  5, -5 ], [  40, 10, -10 ] ], \
                   ['12','C. Europe'       ,[   15.0,  55.0,  10.0, 35.0], [ 10,  5, -5 ], [  20, 10, -10 ] ], \
                   ['13','Turkey'          ,[   15.0,  35.0,  10.0, 35.0], [  5,  5,  0 ], [  20, 10,   0 ] ], \
                   ['14','Ukraine'         ,[   35.0,  55.0,  10.0, 35.0], [  5,  5, -5 ], [  20, 10, -10 ] ], \
                   ['15','C. Asia'         ,[   55.0,  55.0,  10.0, 35.0], [ 10,  5,  0 ], [  30, 10,   0 ] ], \
                   ['16','Russia'          ,[   90.0,  55.0,  10.0, 35.0], [ 15,  5,-10 ], [  60, 20, -20 ] ], \
                   ['17','Middle East'     ,[   60.0, -25.0,  10.0, 35.0], [ 15,  5,  0 ], [  50, 10,   0 ] ], \
                   ['18','India'           ,[   90.0,   0.0,  10.0, 35.0], [ 35,  5,  0 ], [ 120, 20,   0 ] ], \
                   ['19','Korea'           ,[  130.0,  43.0,  10.0, 35.0], [  5,  5, -5 ], [  20, 10, -10 ] ], \
                   ['20','China'           ,[  110.0,  43.0,  10.0, 35.0], [ 35,  5, -5 ], [ 120, 20, -20 ] ], \
                   ['21','S.E. Asia'       ,[  110.0,   0.0,  10.0, 35.0], [ 15,  5,  0 ], [  40, 10,   0 ] ], \
                   ['22','Indonesia'       ,[  130.0,   0.0,  10.0, 35.0], [ 10,  5,  0 ], [  30, 10,   0 ] ], \
                   ['23','Japan'           ,[  150.0,  43.0,  10.0, 35.0], [  5,  5, -5 ], [  20, 10, -10 ] ], \
                   ['24','Oceania'         ,[  157.0, -55.0,  10.0, 35.0], [ 10,  5, -5 ], [  40, 10, -10 ] ], \
                   ['25','Rest S. Asia'    ,[   70.0,   0.0,  10.0, 35.0], [ 10,  5,  0 ], [  30, 10,   0 ] ], \
                   ['26','Rest S. Africa'  ,[   -5.0, -40.0,  10.0, 35.0], [ 15,  5,  0 ], [  50, 10,   0 ] ]  \
                 ]

# Colour codes for IMAGE regions
#COLOURS        = [ '#afd7f9', '#0089e3', '#00287f', '#221d00', '#6e6a00', \
#                   '#aba266', '#8a005d', '#902977', '#d2a0c6', '#330013', \
#                   '#ddc73d', '#b2832f', '#6f3d45', '#290000', '#00003d', \
#                   '#b8aad1', '#3e1673', '#b74a00', '#69a43a', '#c9d0b7', \
#                   '#fff493', '#ffe600', '#003814', '#64746c', '#ed5100', \
#                   '#750039', '#c9c8d0', '#e4e3e7' ]

COLOURS        = [ '#dcdcdc', '#c0c0c0', '#a9a9a9', '#808080', '#dcdcdc', \
                   '#c9d0b7', '#dcdcdc', '#c0c0c0', '#a9a9a9', '#808080', \
                   '#c9d0b7', '#808080', '#c0c0c0', '#a9a9a9', '#808080', \
                   '#dcdcdc', '#c9d0b7', '#808080', '#a9a9a9', '#c9d0b7', \
                   '#dcdcdc', '#808080', '#808080', '#c9d0b7', '#c0c0c0', \
                   '#c9d0b7' ]

nREGIONS       = len(META_IMAGE)

print('Number of regions: ',nREGIONS)
print('Run option:        ',sDATA_OPT)


# In[6]:

# Get chart data
nFILES         =  1
nSTATS         =  3
nMITIGATION    =  5
nTEMPER        =  2
DATA_CHART_ALL = np.zeros((nSTATS,nREGIONS,nTEMPER,nMITIGATION))

iDATA_OPT      = int(sDATA_OPT)
print('Series:',iDATA_OPT)
sSCENARIO      = META_CHART[iDATA_OPT][1]
DATA_DIR       = '/prj/CLIFFTOP/SYNTHESIS/ReviewResponse/'+sSCENARIO+'/CSVoutput/'

for iFILE in range(nFILES):
    DATA_FILE      = DATA_DIR+'Mitigation_Summary_FullUnc_'+META_CHART[iDATA_OPT][2][iFILE]+'.csv'
    print('Opening: '+DATA_FILE)
    DATA_FID       = open(DATA_FILE,'r')
    for iTEMPER in range(nTEMPER):
        INPUT          = DATA_FID.readline().replace('\n','')
        if DEBUG == 'Y': print(INPUT)
        for iSTAT in range(nSTATS):
            INPUT          = DATA_FID.readline().replace('\n','')
            if DEBUG == 'Y': print(INPUT)
            HEADERS    = DATA_FID.readline().replace('\n','').split(',')
            for iREGION in range(nREGIONS):
                DATA_INPUT = DATA_FID.readline().replace('\n','').split(',')
                if DEBUG == 'Y': print(DATA_INPUT)
                for iMITIGATION in range(nMITIGATION):
                    DATA_CHART_ALL[iSTAT,iREGION,iTEMPER,iMITIGATION] = float(DATA_INPUT[iMITIGATION+1])
            if iSTAT == nSTATS-1:
                for iLINE in range(2): DATA_FID.readline()
            else:
                for iLINE in range(5): DATA_FID.readline()

if iFIX == 'Y':
    DATA_CHART     = DATA_CHART_ALL[:,:,:,0:3].squeeze()*44.0/12.0
else:
    DATA_CHART     = DATA_CHART_ALL[:,:,:,0:3].squeeze()
print(DATA_CHART_ALL.shape,DATA_CHART.shape)


# In[7]:



# In[8]:


#Plotting parameters
PLOT_TITLE     = 'IMAGE REGIONS'
PLOT_LABEL     = 'IMAGE REGION CODE'
SET_UNDER      = 'white'
SET_OVER       = 'white'
PLOT_DIR       = '/data/grp/eow/garr/Projects/CLIFFTOP/Paper_Synthesis/PLOTS_MITIGATION/'
PLOT_PART      = 'Synthesis_Paper_Mitigation_IMAGE_Region_'+META_TEMP[iTEMPER_OPT]+'_BECCS'+sBECCS_FACTOR+'_'+sDATE

if iFIX == 'Y':
    PLOT_FILE      = PLOT_PART+'_GtCO2.png'
    CSV_FILE       = PLOT_PART+'_GtCO2.csv'
else:
    PLOT_FILE      = PLOT_PART+'_GtC.png'
    CSV_FILE       = PLOT_PART+'_GtC.csv'

FILE_PLOT      = PLOT_DIR+PLOT_FILE
CLEVELS        = 1+np.arange(nREGIONS+1)
COLOUR_BARS    = ['#7570b3','#1b9e77','#e6ab02']
print('Plot file: '+FILE_PLOT)

#Plotting
FIGURE         = plt.figure(figsize=(16.0,12.0))
AX             = FIGURE.add_subplot(111)
FONT           = {'family' : 'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 10 }
matplotlib.rc('font', **FONT)
rc('text',usetex=True)

# Create polar stereographic Basemap instance.

M              = Basemap(projection='cyl', \
                 llcrnrlat =LAT_START,urcrnrlat=LAT_END,llcrnrlon =LONG_START,urcrnrlon=LONG_END, \
                 fix_aspect=False)
AX             = plt.gca()

# draw coastlines.
M.drawcoastlines(linewidth =0.2)

NCOLOURS       = len(CLEVELS)-1
COL_MAP        = col.ListedColormap(COLOURS[0:NCOLOURS],'indexed')
cm.register_cmap(cmap=COL_MAP)
PALETTE        = COL_MAP
NORM           = col.BoundaryNorm(CLEVELS, NCOLOURS, clip=False)
IMAGE          = M.imshow(DATA_IMAGE,cmap=PALETTE,interpolation='nearest',norm=NORM)


ERROR_BARS     = np.zeros((2,DATA_CHART.shape[3]))

for iREGION in range(nREGIONS):

    BOX_REGION     = META_IMAGE[iREGION][2][:]

    if iFIX == 'N':
        YMAX           = int(META_IMAGE[iREGION][3][0]) 
        YINC           = int(META_IMAGE[iREGION][3][1])
        YMIN           = int(META_IMAGE[iREGION][3][2])
        BOX_REGION[3]  = float(YMAX-YMIN)
        if YMIN < 0.0:
            BOX_REGION[1]  = BOX_REGION[1]+float(YMIN) 
    elif iFIX == 'Y':
        YMAX           = int(META_IMAGE[iREGION][4][0]) 
        YINC           = int(META_IMAGE[iREGION][4][1])
        YMIN           = int(META_IMAGE[iREGION][4][2]) 
        BOX_REGION[3]  = float(30*(YMAX-YMIN)/100)
        if YMIN < 0.0:
            BOX_REGION[1]  = BOX_REGION[1]+int(30*YMIN/100)
            
    if DEBUG == 'Y': print(iREGION,META_IMAGE[iREGION][2],BOX_REGION,YMIN,YMAX)

    AXIS           = inset_axes(AX, width="100%", height="100%",
                     bbox_to_anchor=BOX_REGION,
                     bbox_transform=AX.transData)

    AXIS.set_title(META_IMAGE[iREGION][1])
    AXIS.set_xlim(-0.5,1.5)
    AXIS.set_xticks([])
    AXIS.set_xticklabels([])
    AXIS.set_ylim(YMIN,YMAX)
    AXIS.set_yticks(np.arange(YMIN,YMAX+1,YINC))
    AXIS.set_yticklabels(np.arange(YMIN,YMAX+1,YINC))
    if YMIN < 0:
        AXIS.plot([-0.5,1.5],[0.0,0.0],'k:')
    
    ERROR_BARS[0,:]= DATA_CHART[2,iREGION,iTEMPER_OPT,:]-DATA_CHART[0,iREGION,iTEMPER_OPT,:]
    ERROR_BARS[1,:]= DATA_CHART[0,iREGION,iTEMPER_OPT,:]-DATA_CHART[1,iREGION,iTEMPER_OPT,:]
    
    if DATA_CHART[0,iREGION,iTEMPER_OPT,1] > 0:
        DATA_BOTTOM    = np.zeros((len(DATA_CHART[0,iREGION,iTEMPER_OPT,:])))
        DATA_BOTTOM[2] = DATA_CHART[0,iREGION,iTEMPER_OPT,1]
        plt.bar([0,1,1], DATA_CHART[0,iREGION,iTEMPER_OPT,:],         align='center', width=0.7, color=COLOUR_BARS, edgecolor=COLOUR_BARS, linewidth=1,
        bottom=DATA_BOTTOM, yerr=ERROR_BARS, ecolor='red')
    else:
        plt.bar([0,1,1], DATA_CHART[0,iREGION,iTEMPER_OPT,:],         align='center', width=0.7, color=COLOUR_BARS, edgecolor=COLOUR_BARS, linewidth=1,
        yerr=ERROR_BARS, ecolor='red')

plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, hspace=0.01, wspace=0.01)

if iDISPLAY =='Y':
    plt.show()
    # display onscreen.
else:
    # Write to file
    plt.savefig(FILE_PLOT)


# In[9]:

# Output chart data as CSV

HEADER         = "Region,CH4,Forest,BECCS,"+"CH4_error_pos,Forest_error_pos,BECCS_error_pos,"+ \
                 "CH4_error_neg,Forest_error_neg,BECCS_error_neg"

CSV_FID        = open(PLOT_DIR+CSV_FILE,"w")
print("Writing to CSV file: "+PLOT_DIR+CSV_FILE)
CSV_FID.write(HEADER+"\n")

for iREGION in range(nREGIONS):

    CSV_LINE       = ("%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f" % \
                      (META_IMAGE[iREGION][1], DATA_CHART[0,iREGION,iTEMPER_OPT,0], \
                       DATA_CHART[0,iREGION,iTEMPER_OPT,1],DATA_CHART[0,iREGION,iTEMPER_OPT,2], \
                       DATA_CHART[2,iREGION,iTEMPER_OPT,0],DATA_CHART[2,iREGION,iTEMPER_OPT,1], \
                       DATA_CHART[2,iREGION,iTEMPER_OPT,2],DATA_CHART[1,iREGION,iTEMPER_OPT,0], \
                       DATA_CHART[1,iREGION,iTEMPER_OPT,1],DATA_CHART[1,iREGION,iTEMPER_OPT,2] ))

    CSV_FID.write(CSV_LINE+"\n")

CSV_FID.close()

