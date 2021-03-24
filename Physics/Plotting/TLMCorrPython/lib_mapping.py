#!/usr/bin/env python

###############################################################
# < next few lines under version control, D O  N O T  E D I T >
# $Date: 2013-06-12 16:37:12 -0400 (Wed, 12 Jun 2013) $
# $Revision: 21 $
# $Author: rmahajan $
# $Id: lib_mapping.py 21 2013-06-12 20:37:12Z rmahajan $
###############################################################

'''
mapping.py contains mapping related functions:
setProj(projection)
drawMap(proj = setProj(projection))
'''

import sys
import numpy as np
from mpl_toolkits.basemap  import Basemap

__author__    = "Rahul Mahajan"
__email__     = "rahul.mahajan@nasa.gov"
__copyright__ = "Copyright 2011, NASA / GSFC / GMAO"
__license__   = "GPL"
__status__    = "Prototype"

def setProj(projection,              \
            resolution       = 'i',  \
            cenlat           = None, \
            cenlon           = None, \
            stdlat1          = None, \
            stdlat2          = None, \
            llcrnrlon        = None, \
            llcrnrlat        = None, \
            urcrnrlon        = None, \
            urcrnrlat        = None, \
            boundinglat      = None, \
            width            = None, \
            height           = None, \
            mlabel_int       = 45.0, \
            plabel_int       = 30.0, \
            meridians_labels = None, \
            parallels_labels = None, \
            **kwargs                 \
            ):

    proj = type('PROJECTION', (), {})
    proj.projection = projection
    proj.resolution = resolution

    if   proj.projection == 'stere':
        print '%s projection has been implemented but not supported yet' % (proj.projection)

    elif ( (proj.projection == 'npstere') or (proj.projection == 'spstere') ):
        proj.meridians        = np.arange(-180, 180, 30)
        proj.parallels        = np.arange(-90 ,90, 15)
        proj.lon_0            = 0.0       if ( cenlon           == None ) else cenlon
        proj.meridians_labels = [0,0,1,1] if ( meridians_labels == None ) else meridians_labels
        if (proj.projection == 'npstere'):
            proj.boundinglat      = 45.0      if ( boundinglat      == None ) else abs(boundinglat)
            proj.parallels_labels = [1,0,0,0] if ( parallels_labels == None ) else parallels_labels
        if (proj.projection == 'spstere'):
            proj.boundinglat      = -45.0     if ( boundinglat      == None ) else -abs(boundinglat)
            proj.parallels_labels = [0,1,0,0] if ( parallels_labels == None ) else parallels_labels

    elif ( (proj.projection == 'mill') or (proj.projection == 'merc') ):
        proj.llcrnrlon        = -180.0 if ( llcrnrlon == None ) else llcrnrlon
        proj.llcrnrlat        =  -90.0 if ( llcrnrlat == None ) else llcrnrlat
        proj.urcrnrlon        =  180.0 if ( urcrnrlon == None ) else urcrnrlon
        proj.urcrnrlat        =   90.0 if ( urcrnrlat == None ) else urcrnrlat
        proj.meridians        = np.arange(-180, 180, 45)
        proj.parallels        = np.arange(-90 ,90, 30)
        proj.meridians_labels = [0,0,0,1] if ( meridians_labels == None ) else meridians_labels
        proj.parallels_labels = [1,0,0,0] if ( parallels_labels == None ) else parallels_labels

    elif proj.projection == 'lcc':
        proj.cenlat    = 0.0  if ( cenlat  == None ) else cenlat
        proj.cenlon    = 0.0  if ( cenlon  == None ) else cenlon
        proj.stdlat1   = 0.0  if ( stdlat1 == None ) else stdlat1
        proj.stdlat2   = 0.0  if ( stdlat2 == None ) else stdlat2
        proj.width     = 8.0  if ( width   == None ) else width
        proj.height    = 11.0 if ( height  == None ) else height
        proj.meridians        = np.arange(-180, 180, 45)
        proj.parallels        = np.arange(-90 ,90, 30)
        proj.meridians_labels = [0,0,0,1] if ( meridians_labels == None ) else meridians_labels
        proj.parallels_labels = [0,1,0,0] if ( parallels_labels == None ) else parallels_labels

    elif proj.projection == 'ortho':
        print '{0} projection has been implemented but not supported yet'.format(proj.projection)

    elif proj.projection == 'robin':
        print '{0} projection has been implemented but not supported yet'.format(proj.projection)

    else:
        print 'Error message from : setProj(projection)'
        print '   projection "{0}" has not been implemented yet'.format(proj.projection)
        print '   valid options are: '
        print '   "stere" | "npstere" | "spstere" | "mill" | "merc" | "lcc" | "ortho" | "robin"'
        sys.exit(1)

    proj.meridians = np.arange(-180, 180, mlabel_int)
    proj.parallels = np.arange(-90 ,90,   plabel_int)

    return proj

def drawMap(proj, **kwargs):

    if proj.projection == 'stere':
        map = Basemap(projection = proj.projection, **kwargs)

    elif ( (proj.projection == 'npstere') or (proj.projection == 'spstere') ):
        map = Basemap(projection  = proj.projection, \
                      resolution  = proj.resolution, \
                      lon_0       = proj.cenlon,     \
                      boundinglat = proj.boundinglat,\
                      **kwargs)

    elif ( (proj.projection == 'mill') or (proj.projection == 'merc') ):
        map = Basemap(projection = proj.projection,  \
                      resolution = proj.resolution,  \
                      llcrnrlon  = proj.llcrnrlon,   \
                      llcrnrlat  = proj.llcrnrlat,   \
                      urcrnrlon  = proj.urcrnrlon,   \
                      urcrnrlat  = proj.urcrnrlat,   \
                      **kwargs)

    elif proj.projection == 'lcc':
        map = Basemap(projection = proj.projection,  \
                      resolution = proj.resolution,  \
                      width      = proj.width,       \
                      height     = proj.height,      \
                      lat_0      = proj.cenlat,      \
                      lon_0      = proj.cenlon,      \
                      lat_1      = proj.stdlat1,     \
                      lat_2      = proj.stdlat2,     \
                      **kwargs)

    elif proj.projection == 'ortho':
        map = Basemap(projection = proj.projection, \
                      resolution = proj.resolution, \
                      **kwargs)

    elif proj.projection == 'robin':
        map = Basemap(projection = proj.projection, \
                      resolution = proj.resolution, \
                      **kwargs)

    else:
        print 'Error message from : drawMap(proj)'
        print '   projection "%s" has not been implemented yet' % (proj.projection)
        print '   valid options fort projection are: '
        print '   "stere" | "npstere" | "spstere" | "mill" | "merc" | "lcc" | "ortho" | "robin"'
        sys.exit(1)

    return map
