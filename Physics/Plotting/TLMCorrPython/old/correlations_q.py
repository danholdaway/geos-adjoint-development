#!/usr/bin/env python

###############################################################
# < next few lines under version control, D O  N O T  E D I T >
# $Date$
# $Revision$
# $Author$
# $Id$
###############################################################

__author__    = "Rahul Mahajan"
__email__     = "rahul.mahajan@nasa.gov"
__copyright__ = "Copyright 2013, NASA / GSFC / GMAO"
__license__   = "GPL"
__status__    = "Prototype"

'''
Computes and plots the comparisons and correlations between NLM and TLM perturbation evolution
'''

import os
import sys
import subprocess
import time
import numpy as np
from scipy.stats import t
from netCDF4 import Dataset
from matplotlib import pyplot, ticker, cm
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
from gmaopy.modules.obsstat import *
from lib_mapping import setProj, drawMap

def main():

    args = parse_args()

    if   ( args.calc_type == 'comp_map'  ): main_comparison_map(     args)
    elif ( args.calc_type == 'corr_map'  ): main_correlation_map(    args)
    elif ( args.calc_type == 'corr_prof' ): main_correlation_profile(args)
    else:
        print 'unrecognized calc_type %s' % args.calc_type
        sys.exit(1)

    if ( not args.savefigure ): pyplot.show()

    print '... all done.'

    sys.exit(0)

def main_comparison_map(args):

    grid = read_grid_data(args)

    adate = args.begin_date
    rdate = adate - 3
    vdate = adate + args.fhr

    root_dir = '%s/%s/prog/%s' % (args.root_tmpl,args.expid,rdate.format(format='Y%Y/M%m/D%d'))
    ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))

    filename = '%s/H21-NLM-%s/%s.prog.eta.%s.nc4' % (root_dir,'free',  args.expid,ftime_str)
    Free = read_state_2D(args,filename)

    filename = '%s/H21-NLM-%s/%s.prog.eta.%s.nc4' % (root_dir,'replay',args.expid,ftime_str)
    Replay = read_state_2D(args,filename)

    NLM_pert = XminusY(Replay,Free)

    filename = '%s/H21-TLM-GQ%d-%s/%s.fvpert.eta.%s.nc4' % (root_dir,1,args.norm,args.expid,ftime_str)
    TLM_gq1 = read_state_2D(args,filename)

    ftime_str = '%sz+%sz' % (adate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
    filename = '%s/H21-TLM-GQ%d-%s/%s.fvpert.eta.%s.nc4' % (root_dir,2,args.norm,args.expid,ftime_str)
    TLM_gq2 = read_state_2D(args,filename)

    show_NLMTLM_comparison_map(args,grid,NLM_pert,TLM_gq1,TLM2=TLM_gq2)

    return

def main_correlation_profile(args):

    check_all_files(args)

    grid = read_grid_data(args)

    print 'plotting correlation profile in %s region at %d hours' % (args.region, args.fhr)

    corr_gq1 = []
    corr_gq2 = []

    adates = Dates(args.begin_date.intvalue(),args.end_date.intvalue(),args.interval)
    for adate in adates:

        [NLM_pert, TLM_gq1, TLM_gq2] = read_NLM_TLM_data(args, adate)

        corr_gq1.append(compute_correlation_profile(args,grid,NLM_pert,TLM_gq1))
        corr_gq2.append(compute_correlation_profile(args,grid,NLM_pert,TLM_gq2))

    corr_gq1 = swap_variable_time_dimension(corr_gq1)
    corr_gq2 = swap_variable_time_dimension(corr_gq2)

    stats_corr_gq1 = compute_statistics(corr_gq1)
    stats_corr_gq2 = compute_statistics(corr_gq2)

    if ( args.begin_date == args.end_date ):
        show_NLMTLM_correlation_profile(args,grid,stats_corr_gq1.mean,stats_corr_gq2.mean)
    else:
        show_NLMTLM_correlation_profile(args,grid,stats_corr_gq1.mean,stats_corr_gq2.mean, \
                                 sgq1=stats_corr_gq1.std,sgq2=stats_corr_gq2.std)

    return

def main_correlation_map(args):

    check_all_files(args)

    grid = read_grid_data(args)

    print 'plotting correlation map at level %d hPa and %d hours' % (np.int(grid['plev'][args.level]), args.fhr)
    time.sleep(3.0)

    [NLM_pert, TLM_gq1, TLM_gq2] = read_NLM_TLM_level_data(args)

    corr_gq1 = compute_correlation_map(args,grid,NLM_pert,TLM_gq1)
    corr_gq2 = compute_correlation_map(args,grid,NLM_pert,TLM_gq2)

    show_NLMTLM_correlation_map(args,grid,corr_gq1,corr_gq2)

#    tcrit = t.ppf(confidence/100.0,Nt-1)

#    sig_corr_gq1 = compute_significant_correlation_map(corr_gq1)
#    sig_corr_gq2 = compute_significant_correlation_map(corr_gq2)

#    show_NLMTLM_correlation_map(sig_corr_gq1,sig_corr_gq2,ssig=True)

    return

def parse_args():

    parser = ArgumentParser(description='Compute NLM/TLM correlations', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--calc_type',help='type of correlation calculation',type=str,choices=['comp_map','corr_map','corr_prof'],default='corr_prof',required=False)
    parser.add_argument('--region',help='region to compute correlations',type=str,choices=['global','tropics','nhemp','shemp','midlat'],default='global',required=False)
    parser.add_argument('--expid',help='experiment ID',type=str,default='u001_C180',required=False)
    parser.add_argument('--begin_date',help='date to begin processing',type=str,default='2013010100',required=False)
    parser.add_argument('--end_date',help='date to end processing',type=str,default='2013013100',required=False)
    parser.add_argument('--interval',help='interval between dates',type=int,choices=range(6,30,6),default=24,required=False)
    parser.add_argument('--root_tmpl',help='location of the experiment',type=str,default=os.getenv('NOBACKUP'),required=False)
    parser.add_argument('--level',help='level to compute',type=int,choices=range(0,73),default=62,required=False)
    parser.add_argument('--fhr',help='forecast hour to process',type=int,choices=range(0,78,6),default=24,required=False)
    parser.add_argument('--confidence',help='confidence interval',type=float,default=0.0,required=False)
    parser.add_argument('--norm',help='norm to process',type=str,choices=['txe','twe','PH0','PH1','PH2'],default='txe',required=False)
    parser.add_argument('--variable',help='variable to plot',type=str,choices=['u','v','t','q','qi','ql','o3'],default='t',required=False)
    parser.add_argument('-s','--savefigure',help='save figures as eps and png',action='store_true',required=False)
    parser.add_argument('-v','--verbose',help='print correlations to the screen',action='store_true',required=False)

    args = parser.parse_args()
    args.begin_date = Date(args.begin_date)
    args.end_date   = Date(args.end_date)

    return args

def check_all_files(args):

    print 'check if all necessary files are present between %s and %s' % (args.begin_date.format(format='%Y%m%d_%H'), args.end_date.format(format='%Y%m%d_%H'))

    adates = Dates(args.begin_date.intvalue(),args.end_date.intvalue(),args.interval)
    for adate in adates:

        adate = Date(adate)
        rdate = adate - 3
        vdate = adate + args.fhr

        root_dir = '%s/%s/prog/%s' % (args.root_tmpl,args.expid,rdate.format(format='Y%Y/M%m/D%d'))

        ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
        filename = '%s/H21-NLM-%s/%s.prog.eta.%s.nc4' % (root_dir,'free',  args.expid,ftime_str)
        if ( not os.path.exists(filename) ): break

        ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
        filename = '%s/H21-NLM-%s/%s.prog.eta.%s.nc4' % (root_dir,'replay',args.expid,ftime_str)
        if ( not os.path.exists(filename) ): break

        ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
        filename = '%s/H21-TLM-GQ%d-%s/%s.fvpert.eta.%s.nc4' % (root_dir,1,args.norm,args.expid,ftime_str)
        if ( not os.path.exists(filename) ): break

        ftime_str = '%sz+%sz' % (adate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
        filename = '%s/H21-TLM-GQ%d-%s/%s.fvpert.eta.%s.nc4' % (root_dir,2,args.norm,args.expid,ftime_str)
        if ( not os.path.exists(filename) ): break

    if ( adate != args.end_date ):
            print 'Error occurred on analysis date %s' % (adate.format(format='%Y%m%d_%H'))
            print '%s does not exist' % filename
            sys.exit(0)
    else:
        print 'check done'

    return

def read_grid(filename):
    print 'reading grid  ... %s' % filename
    grid = {}
    nc = Dataset(filename,'r')
    if ( 'lon' in nc.variables ):
        grid['lon'] = np.squeeze(nc.variables['lon'][:])
        grid['Nx'] = len(grid['lon'])
    if ( 'lat' in nc.variables ):
        grid['lat'] = np.squeeze(nc.variables['lat'][:])
        grid['Ny'] = len(grid['lat'])
    if ( 'lev' in nc.variables ):
        grid['lev'] = np.squeeze(nc.variables['lev'][:])
        grid['Nz'] = len(grid['lev'])
    nc.close()
    grid['glon'],grid['glat'] = np.meshgrid(grid['lon'],grid['lat'])

    [ak,bk] = read_vertgrid72()
    grid['plev'] = np.zeros(grid['Nz']+1)
    for i in range(len(grid['plev'])):
        grid['plev'][i] = float(ak[i]*0.01 + bk[i]*1000.0)

    return grid

def read_vertgrid72():

    fh = open('VERTGRID72','r')
    lines = fh.readlines()
    fh.close()

    ak = np.zeros((72+1))
    bk = np.zeros((72+1))

    for line in lines:
        line = line.strip('\n')
        if ( line.startswith('#')  ): continue
        if ( line == 'VERTGRID72:' ): continue
        i = int(line.split()[0]) - 1
        ak[i] = float(line.split()[1])
        bk[i] = float(line.split()[2])

    return [ak, bk]

def read_grid_data(args):

    adate  = args.begin_date
    rdate  = adate - 3
    vdate  = adate + args.fhr

    ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
    root_dir = '%s/%s/prog/%s' % (args.root_tmpl,args.expid,rdate.format(format='Y%Y/M%m/D%d'))
    filename = '%s/H21-NLM-%s/%s.prog.eta.%s.nc4' % (root_dir,'free',args.expid,ftime_str)

    return read_grid(filename)

def read_state_3D(args,filename):
    print 'reading state ... %s' % filename
    state = {}
    nc = Dataset(filename,'r')
#    if ( 'u'     in nc.variables ): state['u']  = np.squeeze(nc.variables['u'][    :])
#    if ( 'v'     in nc.variables ): state['v']  = np.squeeze(nc.variables['v'][    :])
#    if ( 'tv'    in nc.variables ): state['t']  = np.squeeze(nc.variables['tv'][   :])
    if ( 'sphu'  in nc.variables ): state['q']  = np.squeeze(nc.variables['sphu'][ :])

#    if ( 'delp'  in nc.variables ): state['p']  = np.squeeze(nc.variables['delp'][ :])
#    if ( 'qitot' in nc.variables ): state['qi'] = np.squeeze(nc.variables['qitot'][:])
#    if ( 'qltot' in nc.variables ): state['ql'] = np.squeeze(nc.variables['qltot'][:])
#    if ( 'ozone' in nc.variables ): state['o3'] = np.squeeze(nc.variables['ozone'][:])
    nc.close()
    return state

def read_state_2D(args,filename):
    print 'reading state ... %s' % filename
    state = {}
    nc = Dataset(filename,'r')
#    if ( 'u'     in nc.variables ): state['u']  = np.squeeze(nc.variables['u'][    :,args.level])
#    if ( 'v'     in nc.variables ): state['v']  = np.squeeze(nc.variables['v'][    :,args.level])
#    if ( 'tv'    in nc.variables ): state['t']  = np.squeeze(nc.variables['tv'][   :,args.level])
    if ( 'sphu'  in nc.variables ): state['q']  = np.squeeze(nc.variables['sphu'][ :,args.level])

#    if ( 'delp'  in nc.variables ): state['p']  = np.squeeze(nc.variables['delp'][ :,args.level])
#    if ( 'qitot' in nc.variables ): state['qi'] = np.squeeze(nc.variables['qitot'][:,args.level])
#    if ( 'qltot' in nc.variables ): state['ql'] = np.squeeze(nc.variables['qltot'][:,args.level])
#    if ( 'ozone' in nc.variables ): state['o3'] = np.squeeze(nc.variables['ozone'][:,args.level])
    nc.close()
    return state

def read_NLM_TLM_data(args, adate):

    adate = Date(adate)
    rdate = adate - 3
    vdate = adate + args.fhr

    root_dir = '%s/%s/prog/%s' % (args.root_tmpl,args.expid,rdate.format(format='Y%Y/M%m/D%d'))

    ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
    filename = '%s/H21-NLM-%s/%s.prog.eta.%s.nc4' % (root_dir,'free',  args.expid,ftime_str)
    Free = read_state_3D(args,filename)

    ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
    filename = '%s/H21-NLM-%s/%s.prog.eta.%s.nc4' % (root_dir,'replay',args.expid,ftime_str)
    Replay = read_state_3D(args,filename)

    NLM_pert = XminusY(Replay,Free)

    ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
    filename = '%s/H21-TLM-GQ%d-%s/%s.fvpert.eta.%s.nc4' % (root_dir,1,args.norm,args.expid,ftime_str)
    TLM_gq1 = read_state_3D(args,filename)

    ftime_str = '%sz+%sz' % (adate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
    filename = '%s/H21-TLM-GQ%d-%s/%s.fvpert.eta.%s.nc4' % (root_dir,2,args.norm,args.expid,ftime_str)
    TLM_gq2 = read_state_3D(args,filename)

    return [NLM_pert, TLM_gq1, TLM_gq2]

def read_NLM_TLM_level_data(args):

    NLM_pert = []
    TLM_gq1  = []
    TLM_gq2  = []

    adates = Dates(args.begin_date.intvalue(),args.end_date.intvalue(),args.interval)
    Nt = len(adates)
    for adate in adates:

        adate = Date(adate)
        rdate = adate - 3
        vdate = adate + args.fhr

        root_dir = '%s/%s/prog/%s' % (args.root_tmpl,args.expid,rdate.format(format='Y%Y/M%m/D%d'))

        ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
        filename = '%s/H21-NLM-%s/%s.prog.eta.%s.nc4' % (root_dir,'free',  args.expid,ftime_str)
        Free = read_state_2D(args,filename)

        ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
        filename = '%s/H21-NLM-%s/%s.prog.eta.%s.nc4' % (root_dir,'replay',args.expid,ftime_str)
        Replay = read_state_2D(args,filename)

        NLM_pert.append(XminusY(Replay,Free))

        ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
        #filename = '%s/H21-TLM-GQ%d-%s/%s.fvpert.eta.%s.nc4' % (root_dir,1,args.norm,args.expid,ftime_str)
        filename = '%s/H21-TLM-GQ1-PH1/%s.fvpert.eta.%s.nc4' % (root_dir,args.expid,ftime_str)
        TLM_gq1.append(read_state_2D(args,filename))

        #ftime_str = '%sz+%sz' % (adate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
        #filename = '%s/H21-TLM-GQ%d-%s/%s.fvpert.eta.%s.nc4' % (root_dir,2,args.norm,args.expid,ftime_str)
        ftime_str = '%sz+%sz' % (rdate.format(format='%Y%m%d_%H'),vdate.format(format='%Y%m%d_%H'))
        filename = '%s/H21-TLM-GQ1-PH2/%s.fvpert.eta.%s.nc4' % (root_dir,args.expid,ftime_str)
        TLM_gq2.append(read_state_2D(args,filename))

    NLM_pert = swap_variable_time_dimension(NLM_pert)
    TLM_gq1  = swap_variable_time_dimension(TLM_gq1)
    TLM_gq2  = swap_variable_time_dimension(TLM_gq2)

    return [NLM_pert, TLM_gq1, TLM_gq2]

def XminusY(stateX, stateY):
    pert_state = {}
    for var in stateX:
        pert_state[var] = stateX[var] - stateY[var]
    return pert_state

def get_region_ij(args,grid):

    [Nx,Ny,Nz] = [grid['Nx'], grid['Ny'], grid['Nz']]

    if   ( args.region == 'global'  ):
        levels     = np.arange(29,Nz)
        latitudes  = np.arange(0,Ny)
        longitudes = np.arange(0,Nx)
    elif ( args.region == 'tropics' ):
        levels     = np.arange(39,Nz)
        latitudes  = np.arange(120,241)
        longitudes = np.arange(0,Nx)
    elif ( args.region == 'nhemp'   ):
        levels     = np.arange(29,Nz)
        latitudes  = np.arange(Ny/2,Ny)
        longitudes = np.arange(0,Nx)
    elif ( args.region == 'shemp'   ):
        levels     = np.arange(29,Nz)
        latitudes  = np.arange(0,Ny/2)
        longitudes = np.arange(0,Nx)
    elif ( args.region == 'midlat'  ):
        levels     = np.arange(29,Nz)
        latitudes  = np.concatenate([np.arange(0,120),np.arange(241,Ny)])
        longitudes = np.arange(0,Nx)

    return [levels,latitudes,longitudes]

def get_weights(grid,option='Cosine',debug=False):

    [Nx,Ny,Nz] = [grid['Nx'], grid['Ny'], grid['Nz']]

    r_earth = 6378.1
    lat_r = np.pi/180.0 * grid['lat']

    if ( option == 'DANHOLDAWAY' ):
        wght = 2.0 * np.pi**2 * r_earth**2 * np.cos(lat_r) / (Nx*Ny)
    elif ( option == 'SphericalCap' ):
        wght = 2.0 * np.pi * r_earth**2 * (np.sin(lat_r[1:]) - np.sin(lat_r[0:-1])) / Nx
    elif ( option == 'Circumference' ):
        wght = np.sqrt(2.0 * np.pi * r_earth * np.cos(lat_r) / Nx)
    elif ( option == 'Cosine' ):
        wght = np.sqrt(np.cos(lat_r))
    elif ( option == 'One' ):
        wght = np.ones(Ny)

    wght = np.tile(np.transpose(np.tile(wght,(Nx,1))),(Nz,1,1))

    if ( debug ):

        fig = pyplot.figure()

        level = 0

        cmap = cm.get_cmap('jet')
        cmax = np.max(wght)
        cntrs = np.arange(0.0,cmax+cmax/100,cmax/100)

        proj = setProj('mill',llcrnrlat=-89.0,urcrnrlat=89.0)
        map = drawMap(proj)
        map.drawcoastlines(color='0.8')
        x,y = map(grid['glon'],grid['glat'])
        cx = map.contourf(x,y,np.squeeze(wght[level,:,:]),cntrs,cmap=cmap)
        cb = pyplot.colorbar()
        cb.set_clim(vmin=0,vmax=cmax)
        pyplot.title('Weighting function at level = %d' % level)

    return wght

def compute_correlation_profile(args,grid,stateX,stateY):

    print 'computing correlations'

    [levs,lats,lons] = get_region_ij(args,grid)
    Nlevs = len(levs)
    weights = get_weights(grid,option='Cosine')

    corr_state = {}
    for var in stateX:

        corr_state[var] = np.zeros(Nlevs) * np.NaN

        x = np.squeeze(stateX[var][levs[0]:levs[-1]+1,lats[0]:lats[-1]+1,lons[0]:lons[-1]+1])
        y = np.squeeze(stateY[var][levs[0]:levs[-1]+1,lats[0]:lats[-1]+1,lons[0]:lons[-1]+1])
        w = np.squeeze(weights[    levs[0]:levs[-1]+1,lats[0]:lats[-1]+1,lons[0]:lons[-1]+1])

        for k in range(Nlevs):

            wk_vec = np.ravel(np.squeeze(w[k,]),order='F')
            xk_vec = np.ravel(np.squeeze(x[k,]),order='F')
            yk_vec = np.ravel(np.squeeze(y[k,]),order='F')

            cov_xy = np.sum((xk_vec*wk_vec)*(yk_vec*wk_vec))
            var_xx = np.sum((xk_vec*wk_vec)*(xk_vec*wk_vec))
            var_yy = np.sum((yk_vec*wk_vec)*(yk_vec*wk_vec))

            corr_state[var][k] = cov_xy / np.sqrt(var_xx*var_yy)

    return corr_state

def compute_correlation_map(args,grid,stateX,stateY):

    print 'computing correlation map'

    [Nx,Ny,Nz] = [grid['Nx'], grid['Ny'], grid['Nz']]
    Nt = len(Dates(args.begin_date.intvalue(),args.end_date.intvalue(),args.interval))

    corr_state = {}
    for var in stateX:

        corr = np.zeros(Ny*Nx) * np.NaN

        x = np.reshape(stateX[var][:],(Nt,Ny*Nx))
        y = np.reshape(stateY[var][:],(Nt,Ny*Nx))

        for ij in range(Ny*Nx):

            xvec = np.squeeze(x[:,ij])
            yvec = np.squeeze(y[:,ij])

            cov_xy = np.sum(xvec*yvec)
            var_xx = np.sum(xvec*xvec)
            var_yy = np.sum(yvec*yvec)

            corr[ij] = cov_xy / np.sqrt(var_xx*var_yy)

        corr_state[var] = np.reshape(corr,(Ny,Nx))

    return corr_state

def compute_significant_correlation_map(corr_state):

    print 'computing statistical significant correlation map'

    sig_corr_state = {}
    for var in corr_state:

        tstat = corr_state[var] * np.sqrt(Nt-1.0) / np.sqrt(1.0 - corr_state[var]**2.0)
        ssig = tstat > tcrit
        sig_corr_state[var] = np.ma.masked_where(ssig,corr_state[var])

    return sig_corr_state

def swap_variable_time_dimension(x):
    xswap = {}
    for var in x[0]:
        xswap[var] = []
        for q in range(len(x)):
            xswap[var].append(x[q][var])

    return xswap

def compute_statistics(x):
    stats = type('statistics',(),{})
    stats.mean = {}
    stats.std  = {}
    for var in x:
        stats.mean[var] = np.mean(x[var][:],axis=0)
        stats.std[ var] = np.std( x[var][:],axis=0,ddof=1)

    return stats


def show_NLMTLM_comparison_map(args,grid,NLM,TLM,TLM2=None):

    cmap = cm.get_cmap('jet')

    nrows = 2 if ( TLM2 == None ) else 3

    proj = setProj('mill',llcrnrlat=-85.0,urcrnrlat=85.0)

    for i, var in enumerate(TLM):

#        if ( var not in ['u','t'] ): continue
        if ( var not in ['q'] ): continue

        nlm_var  = NLM[var][:]
        tlm_var  = TLM[var][:]

        if ( nrows == 3 ): tlm2_var  = TLM2[var][:]

        if ( nrows == 3 ): cmax = np.max([np.max(np.max(np.abs(nlm_var))),np.max(np.max(np.abs(tlm_var))),np.max(np.max(np.abs(tlm2_var)))])
        else:              cmax = np.max([np.max(np.max(np.abs(nlm_var))),np.max(np.max(np.abs(tlm_var)))])

        cmax = np.ceil(cmax) if ( cmax > 1.0 ) else np.ceil(cmax*10.0**(np.abs(np.round(np.log10(cmax)))+1))*10**np.round(np.log10(cmax)-1)
        print 'cmax = %f' % cmax
        cmax = 5.0

        cntrs = np.arange(-cmax,cmax+cmax/10.0,cmax/10.0)

        fig = pyplot.figure()
        fig.clf()
        titlestr = 'fhr = %d, %s @ %d hPa' % (args.fhr, var.upper(), np.int(grid['plev'][args.level]))
        pyplot.suptitle(titlestr,fontsize=20)

        pyplot.subplot(nrows,1,1)
        map = drawMap(proj)
        map.drawcoastlines(color='0.5',linewidth=0.4)
        x,y = map(grid['glon'],grid['glat'])
        cx1 = map.contourf(x,y,nlm_var,cntrs,cmap=cmap,extend='both')
        pyplot.title('NLM',fontsize=16)

        pyplot.subplot(nrows,1,2)
        map = drawMap(proj)
        map.drawcoastlines(color='0.5',linewidth=0.4)
        x,y = map(grid['glon'],grid['glat'])
        cx2 = map.contourf(x,y,tlm_var,cntrs,cmap=cmap,extend='both')
        pyplot.title('CTL',fontsize=16)

        if ( nrows == 3 ):
            pyplot.subplot(nrows,1,3)
            map = drawMap(proj)
            map.drawcoastlines(color='0.5',linewidth=0.4)
            x,y = map(grid['glon'],grid['glat'])
            cx3 = map.contourf(x,y,tlm2_var,cntrs,cmap=cmap,extend='both')
            pyplot.title('GQ1',fontsize=16)

        pyplot.subplots_adjust(left=0.1,right=0.8,bottom=0.1,top=0.9)

        cbax = pyplot.axes([0.8, 0.1, 0.03, 0.8])
        cb = pyplot.colorbar(cax=cbax, ticks=cntrs[::2])
        cb.set_clim(vmin=-cmax,vmax=cmax)

        if ( args.savefigure ): savecomparisonmapfigure(args,grid,var)

    return

def show_NLMTLM_correlation_map(args,grid,gq1,gq2,ssig=False):

    cmap = cm.get_cmap('jet')

    for i, var in enumerate(gq1):

        if ( var not in ['q'] ): continue

        proj = setProj('mill',llcrnrlat=-85.0,urcrnrlat=85.0,resolution='l')

        cmax = 1.0
        cntrs = np.arange(-cmax,cmax+cmax/10,cmax/10)

        fig = pyplot.figure()
        fig.clf()
#        titlestr = 'fhr = %d, %s @ %d hPa' % (args.fhr, var.upper(), np.int(grid['plev'][args.level]))
#        pyplot.suptitle(titlestr,fontsize=20)

        pyplot.subplot(2,1,1)
        map = drawMap(proj)
        map.drawcoastlines(color='0.5',linewidth=0.4)
        x,y = map(grid['glon'],grid['glat'])
        cx1 = map.contourf(x,y,gq1[var],cntrs,cmap=cmap)
        pyplot.title('CTL',fontsize=16)

        pyplot.subplot(2,1,2)
        map = drawMap(proj)
        map.drawcoastlines(color='0.5',linewidth=0.4)
        x,y = map(grid['glon'],grid['glat'])
        cx2 = map.contourf(x,y,gq2[var],cntrs,cmap=cmap)
        pyplot.title('GQ1',fontsize=16)

        pyplot.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)

#        cbax = pyplot.axes([0.8, 0.1, 0.03, 0.8])
#        cb = pyplot.colorbar(cax=cbax, ticks=cntrs[::2])
#        cb.set_clim(vmin=-cmax,vmax=cmax)

        if ( args.savefigure ): savecorrelationmapfigure(args,grid,var,ssig=ssig)

    return

def show_NLMTLM_correlation_profile(args,grid,gq1,gq2,sgq1=None,sgq2=None):

    [levs,lats,lons] = get_region_ij(args,grid)

    if ( args.verbose ):
        for var in gq1:
            print '=============== %s ================' % var.upper()
            print 'Level         GQ = 1         GQ = 2'
            for k, lev in enumerate(levs):
                print '%5d %14.6f %14.6f' % (lev,gq1[var][k],gq2[var][k])
            mgq1 = np.mean(np.ma.masked_array(gq1[var][:],np.isnan(gq1[var][:])))
            mgq2 = np.mean(np.ma.masked_array(gq2[var][:],np.isnan(gq2[var][:])))
            print ' Mean %14.6f %14.6f' % (mgq1, mgq2)
            print

    xvec = np.arange(0,110,20) / 100.0
    zvec = grid['plev'][levs+1]

    if ( args.begin_date == args.end_date ):
        date_str = '%s' % ( args.begin_date.format(format='%Y/%m/%d %Hz') )
    else:
        date_str = '%s - %s' % ( args.begin_date.format(format='%Y/%m/%d %Hz'), args.end_date.format(format='%Y/%m/%d %Hz') )
    suptitle_str = '%s, fhr = %02d\nregion = %s, norm = %s' % (date_str,args.fhr,args.region,args.norm)

    for var in gq1:
        pyplot.figure()
        pyplot.clf()
        pyplot.subplot(1,1,1)
        mysubplot(var,xvec,zvec,gq1[var][:],gq2[var][:],sgq1=sgq1[var][:],sgq2=sgq2[var][:])
        title_str = '%s, var = %s' %(suptitle_str,var.upper())
        pyplot.suptitle(title_str)
        if ( args.savefigure ): savecorrelationprofilefigure(args,var)

    return

def mysubplot(var,xvec,zvec,gq1,gq2,sgq1=None,sgq2=None):
    l1 = pyplot.semilogy(gq1,zvec,'m',basey=10,linewidth=2.0)
    l2 = pyplot.semilogy(gq2,zvec,'c',basey=10,linewidth=2.0)
    l1 = pyplot.semilogy(gq1,zvec,'o',basey=10,mfc='m',mec='m')
    l2 = pyplot.semilogy(gq2,zvec,'o',basey=10,mfc='c',mec='c')
    if ( sgq1 != None ):
        pyplot.fill_betweenx(zvec,gq1-sgq1,x2=gq1+sgq1,alpha=0.2,edgecolor='m',facecolor='m',lw=0.0,antialiased=True)
    if ( sgq2 != None ):
        pyplot.fill_betweenx(zvec,gq2-sgq2,x2=gq2+sgq2,alpha=0.2,edgecolor='c',facecolor='c',lw=0.0,antialiased=True)
    mgq1 = np.mean(np.ma.masked_array(gq1,np.isnan(gq1)))
    mgq2 = np.mean(np.ma.masked_array(gq2,np.isnan(gq2)))
    pyplot.semilogy(mgq1,1000.0,'*',basey=10,mfc='m',mec='m',ms=12,clip_on=False)
    pyplot.semilogy(mgq2,1000.0,'*',basey=10,mfc='c',mec='c',ms=12,clip_on=False)
    FormatAxes_(var,xvec,zvec)
    pyplot.legend((l1,l2),('CTL','GQ1'),frameon=False,numpoints=1)
    return

def FormatAxes_(var,xvec,zvec):
    ax = pyplot.gca()
    ax.set_xticks(xvec)
    ax.set_xlim([0.0,1.0])
    ax.invert_yaxis()
    ax.set_ylim([1000,100])
    majorLocator = ticker.MultipleLocator(100.)
    majorFormatter = ticker.FormatStrFormatter('%4.0f')
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.set_xlabel('correlation')
    ax.set_ylabel('pressure (hPa)')
    return

def savecorrelationprofilefigure(args,var):
    if ( args.begin_date == args.end_date ):
        date_str = '%s' % ( args.begin_date.format(format='%Y%m%d%H') )
    else:
        date_str = '%s-%s' % ( args.begin_date.format(format='%Y%m%d%H'), args.end_date.format(format='%Y%m%d%H') )
    fname = '%s-f%02d_%s_%s_%s' % (date_str,args.fhr,args.region,args.norm,var)
    print 'saving correlation figure ... %s.[png,eps,pdf]' % fname
    pyplot.savefig(fname+'.png',dpi=120,orientation='landscape',format='png')
    # saving figure as pdf and eps since transparency does not work with epstopdf
    pyplot.savefig(fname+'.eps',dpi=300,orientation='landscape',format='eps')
    pyplot.savefig(fname+'.pdf',dpi=300,orientation='landscape',format='pdf')

    return

def savecorrelationmapfigure(args,grid,var,ssig=False):
    date_str = '%s-%s' % ( args.begin_date.format(format='%Y%m%d%H'), args.end_date.format(format='%Y%m%d%H') )
    if ( ssig ):
        fname = '%s-f%02d_%dhPa_%s-sig' % (date_str,args.fhr,np.int(grid['plev'][args.level]),var)
    else:
        fname = '%s-f%02d_%dhPa_%s' % (date_str,args.fhr,np.int(grid['plev'][args.level]),var)
    print 'saving correlation map figure ... %s.[png]' % fname
    pyplot.savefig(fname+'.png',dpi=120,orientation='landscape',format='png')
    # saving figure as pdf and eps since transparency does not work with epstopdf
#    pyplot.savefig(fname+'.pdf',dpi=300,orientation='landscape',format='pdf')

    return

def savecomparisonmapfigure(args,grid,var):
    fname = '%s-f%d_%dhPa_%s' % (args.begin_date.format(format='%Y%m%d%Hz'),args.fhr,np.int(grid['plev'][args.level]),var)
    print 'saving comparison map figure ... %s.[png]' % fname
    pyplot.savefig(fname+'.png',dpi=120,orientation='landscape',format='png')
    # saving figure as pdf and eps since transparency does not work with epstopdf
#    pyplot.savefig(fname+'.pdf',dpi=300,orientation='landscape',format='pdf')

def _exec(args):
    try:
        result = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        result = e.output

    return result

if __name__ == '__main__': main()
