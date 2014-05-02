__author__ = 'mpopovic'

from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com
import sqlite3 as lite
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pandas.io.sql as psql
import pandas.io.parsers as pp
import matplotlib.image as mpimg
import scipy.signal as signal
from matplotlib.colors import LogNorm

def fill_zeros(s,k):
    while len(s)< k:
        s = '0'+s
    return s

inPath = '/data/biophys/etournay/'
inName = 'WT_25deg_111102'

Rdata_path = '/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/WT_25deg_111102/shear_contrib/'
triList_name = 'triList.RData'
Ta_name = 'Ta_t.RData'
ro.r('load("'+Rdata_path+triList_name+'")')
triList_df = com.load_data('triList')
ro.r('load("'+Rdata_path+Ta_name+'")')
Ta_df = com.load_data('Ta_t')
ro.r('load("/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/WT_25deg_111102/roi_bt/lgRoiInclDead.RData")')
roi_df = com.load_data('lgRoiInclDead')
roi_balde = roi_df[roi_df['roi']=='blade']


inPath = '/data/biophys/etournay/'
inName = 'WT_25deg_111102'
inDB = inPath+'DB/'+inName+'/'+inName+'.sqlite'
con = lite.connect(inDB)
cells_df = psql.frame_query('SELECT * FROM cells WHERE cell_id>10000', con)
cells_df_blade = cells_df[cells_df['cell_id'].isin(roi_balde['cell_id'])]
time_data = psql.frame_query('SELECT * FROM timepoints', con)

frames = cells_df_blade['frame'].unique()
blade_area = np.array([np.sum(cells_df_blade[cells_df_blade['frame']==f]['area']) for f in frames])
blade_av_x = np.array([np.sum(cells_df_blade[cells_df_blade['frame']==f]['center_x']*cells_df_blade[cells_df_blade['frame']==f]['area']) for f in frames])
blade_av_y = np.array([np.sum(cells_df_blade[cells_df_blade['frame']==f]['center_y']*cells_df_blade[cells_df_blade['frame']==f]['area']) for f in frames])
blade_av_xx = np.array([np.sum(cells_df_blade[cells_df_blade['frame']==f]['center_x']*cells_df_blade[cells_df_blade['frame']==f]['center_x']*cells_df_blade[cells_df_blade['frame']==f]['area']) for f in frames])
blade_av_xy = np.array([np.sum(cells_df_blade[cells_df_blade['frame']==f]['center_x']*cells_df_blade[cells_df_blade['frame']==f]['center_y']*cells_df_blade[cells_df_blade['frame']==f]['area']) for f in frames])
blade_av_yy = np.array([np.sum(cells_df_blade[cells_df_blade['frame']==f]['center_y']*cells_df_blade[cells_df_blade['frame']==f]['center_y']*cells_df_blade[cells_df_blade['frame']==f]['area']) for f in frames])

m_xx = (blade_av_xx - blade_av_x*blade_av_x/blade_area)/blade_area**2
m_xy = (blade_av_xy - blade_av_x*blade_av_y/blade_area)/blade_area**2
m_yy = (blade_av_yy - blade_av_y*blade_av_y/blade_area)/blade_area**2

s = 0.5*np.log(m_xx*m_yy - m_xy**2)
Q = np.arcsinh(0.5*np.sqrt((m_xx-m_yy)**2.+(2*m_xy)**2.)/np.exp(s))
twophi = np.arctan2((2*m_xy),(m_xx-m_yy))
Q1 = Q*np.cos(twophi)
Q2 = Q*np.sin(twophi)



shear_path = '/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/movies_data/wild_type/111102/blade/dualMarginDeformation.dat'
shear_df = pp.read_csv(shear_path, sep='\t')
shear1 = np.cumsum(0.5*(shear_df['Delta u_{xx};']-shear_df['Delta u_{yy}']))
shear2 = np.cumsum(0.5*(shear_df['Delta u_{xy};']+shear_df['Delta u_{yx};']))

Q_time = 16+time_data['time_sec'][:201]/3600
shear_time = (np.array(Q_time[1:])+ np.array(Q_time[:-1]))/2.


plt.figure()
plt.plot(Q_time,0.5*(Q1-Q1[0]), label='0.5 (Q1-Q1[0])')
plt.plot(shear_time, shear1, label = 'shear1')
plt.xlabel('time[h]')
plt.grid()
plt.ylim(0,0.25)
plt.legend(loc = 'best')
plt.savefig('WT_25deg_111102_shear1_Q1.png')
plt.show()

plt.figure()
plt.plot(Q_time,0.5*(Q2-Q2[0]), label = '0.5 (Q2-Q2[0])')
plt.plot(shear_time,shear2, label = 'shear2')
plt.xlabel('time[h]')
plt.grid()
plt.legend(loc='best')
plt.savefig('WT_25deg_111102_shear2_Q2.png')
plt.show()