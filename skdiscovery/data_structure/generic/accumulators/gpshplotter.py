# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Authors: Victor Pankratius, Justin Li, Cody Rude
# This software is part of the NSF DIBBS Project "An Infrastructure for
# Computer Aided Discovery in Geoscience" (PI: V. Pankratius) and 
# NASA AIST Project "Computer-Aided Discovery of Earth Surface 
# Deformation Phenomena" (PI: V. Pankratius)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import matplotlib.pyplot as plt
import numpy as np
import skdaccess.utilities.pbo_util as pbo_utils
import skdiscovery.data_structure.series.analysis.mogi as mogi
from mpl_toolkits.basemap import Basemap
from skdiscovery.data_structure.framework import PipelineItem

from skdiscovery.utilities.patterns import pbo_tools


class GPSHPlotter(PipelineItem):
    ''' 
    Plots results from General_Component_Analysis, for the GPS horizontal or vertical components
    '''
    
    def __init__(self, str_description, comp_name, mogi_name = None, pca_dir = 'H', pca_comp = 0, scaleFactor = 2.5, offset=.15, KF_tau = 0, errorEllipses = False, map_resolution='i'):
        '''
        Initialize GPHSHPlotter

        @param str_description: String describing accumulator
        @param comp_name: Name of the GPCA results for accessing the GPCA output
        @param mogi_name: Name of the Mogi results (optional)
        @param pca_dir: PCA direction to plot, horizontal (H) or vertical (V)
        @param pca_comp: The PCA component that will be plotted
        @param scaleFactor: Scale factor for arrows
        @param offset: Offset for plotting larger area on map
        @param KF_tau: Tau used in kalman filter
        @param errorEllipses: Boolean indicating whether or not to plot errorEllipses
        @param map_resolution: Resolution of map features (coastline) to use
        '''

        self.dir_sign = 0        
        self.pca_dir = pca_dir
        self.pca_comp = pca_comp
        self.scaleFactor = scaleFactor
        self.offset = offset
        self.errorE = errorEllipses
        self.KF_tau = KF_tau
        self._bmap_res = map_resolution
        self.comp_name = comp_name
        self.mogi_name = mogi_name
        super(GPSHPlotter, self).__init__(str_description)
        
    def process(self, obj_data):
        ''' 
        Plot the General Component Analysis results present stored in obj_data. 

        Saves the basemap in obj_data results.
        
        @param obj_data: Data Wrapper that holds component analysis HPCA
        '''
        HPCA_name = self.comp_name
        Mogi_name = self.mogi_name
        pca_comp  = self.pca_comp

        plt.figure()

        meta_data = obj_data.info()
        try:
            station_list = obj_data.get().minor_axis
        except AttributeError:
            station_list = list(obj_data.get().keys())

        lat_range, lon_range = pbo_utils.getLatLonRange(meta_data, station_list)
        coord_list = pbo_utils.getStationCoords(meta_data, station_list)

        # Create a map projection of area
        offset = self.offset
        bmap = Basemap(llcrnrlat=lat_range[0] - offset, urcrnrlat=lat_range[1] + offset, llcrnrlon=lon_range[0] - offset, urcrnrlon=lon_range[1] + offset,
                       projection='gnom', lon_0=np.mean(lon_range), lat_0=np.mean(lat_range), resolution=self._bmap_res)

        # bmap.fillcontinents(color='white')
        bmap.drawmapboundary(fill_color='white')

        # Draw just coastlines, no lakes
        for i,cp in enumerate(bmap.coastpolygons):
             if bmap.coastpolygontypes[i]<2:
                bmap.plot(cp[0],cp[1],'k-')

        parallels = np.arange(np.round(lat_range[0]-offset,decimals=1),np.round(lat_range[1]+offset,decimals=1),.1)
        meridians = np.arange(np.round(lon_range[0]-offset,decimals=1),np.round(lon_range[1]+offset,decimals=1),.1)

        bmap.drawmeridians(meridians, labels=[0,0,0,1],fontsize=14, rotation=90)
        bmap.drawparallels(parallels, labels=[1,0,0,0],fontsize=14)

        
        pca_results = obj_data.getResults()[HPCA_name]
        pca = pca_results['CA']

        lonscale = 1
        latscale = 1
        scaleFactor = self.scaleFactor

        if self.pca_dir == 'V':
            station_lat_list, station_lon_list, ev_lat_list, ev_lon_list, dir_sign = pbo_tools.dirEigenvectors(coord_list, pca.components_[pca_comp],pdir='V')
            self.dir_sign = dir_sign            
            pca_results['Projection'] *= dir_sign
            ev_lat_list *= latscale
            ev_lon_list *= lonscale
                
            # Plot station coords
            for coord in coord_list:
                bmap.plot(coord[1], coord[0], 'bo', markersize=8, latlon=True)
                x,y = bmap(coord[1], coord[0])
                plt.text(x-(1+np.sign(ev_lon_list[coord_list.index(coord)]))*900+250,
                         y-(1+np.sign(ev_lat_list[coord_list.index(coord)]))*100+450,
                         station_list[coord_list.index(coord)],fontsize=14)
    
            bmap.quiver(station_lon_list, station_lat_list, ev_lon_list, ev_lat_list, latlon=True, scale = scaleFactor)
            
            ax_x = plt.gca().get_xlim()
            ax_y = plt.gca().get_ylim()
            x,y = bmap(ax_x[0]+.1*(ax_x[1]-ax_x[0]), ax_y[0]+.1*(ax_y[1]-ax_y[0]),inverse=True)
            bmap.quiver(x, y, 0, .2, latlon=True, scale = scaleFactor, headwidth=3,headlength=3)
            plt.text(ax_x[0]+.1*(ax_x[1]-ax_x[0])-650, ax_y[0]+.1*(ax_y[1]-ax_y[0])-1000,'20%', fontsize=14)
            
        else:
            station_lat_list, station_lon_list, ev_lat_list, ev_lon_list, dir_sign = pbo_tools.dirEigenvectors(coord_list, pca.components_[pca_comp])
            self.dir_sign = dir_sign            
            pca_results['Projection'] *= dir_sign
            ev_lat_list *= latscale
            ev_lon_list *= lonscale
                
            # Plot station coords
            for coord in coord_list:
                bmap.plot(coord[1], coord[0], 'bo', markersize=8, latlon=True)
                x,y = bmap(coord[1], coord[0])
                plt.text(x-(1+np.sign(ev_lon_list[coord_list.index(coord)]))*900+250,
                         y-(1+np.sign(ev_lat_list[coord_list.index(coord)]))*800+450,
                         station_list[coord_list.index(coord)], fontsize=14)
    
            bmap.quiver(station_lon_list, station_lat_list, ev_lon_list, ev_lat_list, latlon=True, scale = scaleFactor)
            
            ax_x = plt.gca().get_xlim()
            ax_y = plt.gca().get_ylim()
            x,y = bmap(ax_x[0]+.1*(ax_x[1]-ax_x[0]), ax_y[0]+.1*(ax_y[1]-ax_y[0]),inverse=True)
            bmap.quiver(x, y, 0, .2, latlon=True, scale = scaleFactor, headwidth=3,headlength=3)
            plt.text(ax_x[0]+.1*(ax_x[1]-ax_x[0])-650, ax_y[0]+.1*(ax_y[1]-ax_y[0])-1000,'20%', fontsize=14)
    
            # Plotting Mogi source
            if Mogi_name != None:
                mogi_res = obj_data.getResults()[Mogi_name]
                bmap.plot(mogi_res['lon'], mogi_res['lat'], "g^", markersize = 10, latlon=True)
                mogi_x_disp, mogi_y_disp = pbo_tools.MogiVectors(mogi_res,station_lat_list,station_lon_list)
                bmap.quiver(station_lon_list, station_lat_list, mogi_x_disp*dir_sign, mogi_y_disp*dir_sign,
                            latlon=True, scale=scaleFactor,color='red')
                            
            # Plot error ellipses for the PCA
            if self.errorE:
                ax = plt.gca()
                yScale = (bmap.urcrnrlat - bmap.llcrnrlat)/scaleFactor
                xScale = (bmap.urcrnrlon - bmap.llcrnrlon)/scaleFactor
                midY = (bmap.urcrnrlat + bmap.llcrnrlat)/2
                midX = (bmap.urcrnrlon + bmap.llcrnrlon)/2
                from matplotlib.patches import Ellipse
 
                n=len(pca_results['Projection'][:,0])
                tau=self.KF_tau
                delta_t = np.arange(-(n-1),n+1)
                mseq = (1-np.abs(delta_t)/n)
                rdelt = np.exp(-np.abs(delta_t)/tau)
                neff = n/np.sum(mseq*rdelt)
                eigval = pca.explained_variance_
                aaTs = [np.outer(pca.components_[ii,:],pca.components_[ii,:].T) for ii in range(pca.components_.shape[0])]
                VVs = [eigval[ii]/neff*np.sum([eigval[k]/(eigval[k]-eigval[ii])**2*aaTs[k] for k in (j for j in range(pca.components_.shape[0]) if j != ii)],axis=0) for ii in range(pca.components_.shape[0])]
                sigmas = np.diag(VVs[0])**(1/2)                
                
                for kk in range(len(station_lon_list)):
                    vlon = ev_lon_list[kk]
                    vlat = ev_lat_list[kk]
                    slon = station_lon_list[kk]
                    slat = station_lat_list[kk]
                    Elat = sigmas[2*kk]
                    Elon = sigmas[2*kk+1]
                    cir_w, cir_h = np.array(bmap(midX+Elon/scaleFactor,midY+Elat/scaleFactor*.85))-np.array(bmap(midX,midY))
                
                    x,y = bmap(slon+vlon*xScale*.95,slat+vlat*yScale*.85)
                    # if need to rotate ellipse, np.arctan2(vlat,vlon)*180/np.pi
                    etest = Ellipse(xy=(x,y),width=cir_w,height=cir_h,angle=0,
                                    edgecolor='k',fc='w',lw=1,zorder=-1)
                    ax.add_artist(etest);
                    
                           
        obj_data.addResult(self.str_description, bmap)
