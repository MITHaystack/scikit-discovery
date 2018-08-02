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


import skdaccess.utilities.pbo_util as pbo_tools
import skdiscovery.data_structure.series.analysis.mogi as mogi

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt



def multiCaPlot(pipeline, mogiFlag=False, offset=.15, direction='H',pca_comp=0,scaleFactor=2.5,map_res='i'):
    '''
    The multiCaPlot function generates a geographic eigenvector plot of several pipeline runs
    
    This function plots multiple pipeline runs over perturbed pipeline
    parameters. The various perturbations are plotted more
    transparently (alpha=.5), while the median eigen_vector and Mogi
    inversion are plotted in solid blue and red

    @param pipeline: The pipeline object with multiple runs
    @param mogiFlag: Flag to indicate plotting the Mogi source as well as the PCA
    @param offset: Offset for padding the corners of the generated map
    @param direction: Indicates the eigenvectors to plot. Only Horizontal component is currently supported ('H')
    @param pca_comp: Choose the PCA component to use (integer)
    @param scaleFactor: Size of the arrow scaling factor
    @param map_res: Map data resolution for Basemap ('c', 'i', 'h', 'f', or None)
    '''
    
    # as this is a multi_ca_plot function, assumes GPCA
    plt.figure();
    
    meta_data = pipeline.data_generator.meta_data
    station_list = pipeline.data_generator.station_list
    
    lat_range, lon_range = pbo_tools.getLatLonRange(meta_data, station_list)
    coord_list = pbo_tools.getStationCoords(meta_data, station_list)
    
    # Create a map projection of area
    bmap = Basemap(llcrnrlat=lat_range[0] - offset, urcrnrlat=lat_range[1] + offset, llcrnrlon=lon_range[0] - offset, urcrnrlon=lon_range[1] + offset,
               projection='gnom', lon_0=np.mean(lon_range), lat_0=np.mean(lat_range), resolution=map_res)
    
    # bmap.fillcontinents(color='white')
    # bmap.drawmapboundary(fill_color='white')
    bmap.drawmapboundary(fill_color='#41BBEC');
    bmap.fillcontinents(color='white')
    
    # Draw just coastlines, no lakes
    for i,cp in enumerate(bmap.coastpolygons):
        if bmap.coastpolygontypes[i]<2:
            bmap.plot(cp[0],cp[1],'k-')
    
    parallels = np.arange(np.round(lat_range[0]-offset,decimals=1),np.round(lat_range[1]+offset,decimals=1),.1)
    meridians = np.arange(np.round(lon_range[0]-offset,decimals=1),np.round(lon_range[1]+offset,decimals=1),.1)
    
    bmap.drawmeridians(meridians, labels=[0,0,0,1])
    bmap.drawparallels(parallels, labels=[1,0,0,0])
    
    # Plot station coords
    for coord in coord_list:
        bmap.plot(coord[1], coord[0], 'ko', markersize=6, latlon=True,zorder=12)
        x,y = bmap(coord[1], coord[0])
        plt.text(x+250,y-450,station_list[coord_list.index(coord)],zorder=12)
    
    # loop over each pipeline run
    elatmean = np.zeros(len(station_list))
    elonmean = np.zeros_like(elatmean)
    # check if want to plot Mogi as well
    if mogiFlag:
        avg_mogi = np.array([0.,0.])
        mlatmean = np.zeros_like(elatmean)
        mlonmean = np.zeros_like(elatmean)
    
    for nrun in range(len(pipeline.RA_results)):
        pca = pipeline.RA_results[nrun]['GPCA']['CA']
        station_lat_list, station_lon_list, ev_lat_list, ev_lon_list, dir_sign = pbo_tools.dirEigenvectors(coord_list, pca.components_[pca_comp])
    
        elatmean += ev_lat_list
        elonmean += ev_lon_list
        # plot each run in light blue
        bmap.quiver(station_lon_list, station_lat_list, ev_lon_list, ev_lat_list, latlon=True,
                    scale = scaleFactor, alpha = .25, color = 'blue',zorder=11)
    
        if mogiFlag:
            mogi_res = pipeline.RA_results[nrun]['Mogi']
            avg_mogi += np.array([mogi_res['lon'], mogi_res['lat']])
            mogi_x_disp, mogi_y_disp = mogi.MogiVectors(mogi_res,station_lat_list,station_lon_list)
            mlatmean += mogi_y_disp
            mlonmean += mogi_x_disp
    
            bmap.plot(mogi_res['lon'], mogi_res['lat'], "g^", markersize = 10, latlon=True, alpha = .25,zorder=12)
            bmap.quiver(station_lon_list, station_lat_list, mogi_x_disp*dir_sign, mogi_y_disp*dir_sign,
                        latlon=True, scale=scaleFactor,color='red', alpha = .25,zorder=11)
    
    #plot the mean ev in blue
    elatmean = elatmean/len(pipeline.RA_results)
    elonmean = elonmean/len(pipeline.RA_results)
    bmap.quiver(station_lon_list, station_lat_list, elonmean, elatmean,
                latlon=True, scale = scaleFactor, color = 'blue', alpha = 1,zorder=11)
    if mogiFlag:
        # plot mean mogi results
        avg_mogi = avg_mogi/len(pipeline.RA_results)
        mlatmean = mlatmean/len(pipeline.RA_results)
        mlonmean = mlonmean/len(pipeline.RA_results)
        bmap.plot(avg_mogi[0], avg_mogi[1], "g^", markersize = 10, latlon=True, alpha = 1,zorder=12)
        bmap.quiver(station_lon_list, station_lat_list, mlonmean*dir_sign, mlatmean*dir_sign,
                    latlon=True, scale=scaleFactor,color='red', alpha = 1,zorder=11)
    
    
    ax_x = plt.gca().get_xlim()
    ax_y = plt.gca().get_ylim()
    x,y = bmap(ax_x[0]+.1*(ax_x[1]-ax_x[0]), ax_y[0]+.1*(ax_y[1]-ax_y[0]),inverse=True)
    bmap.quiver(x, y, 0, .2, latlon=True, scale = scaleFactor, headwidth=3,headlength=3,zorder=11)
    plt.text(ax_x[0]+.1*(ax_x[1]-ax_x[0])-650, ax_y[0]+.1*(ax_y[1]-ax_y[0])-1000,'20%',zorder=11)
