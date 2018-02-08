# The MIT License (MIT)
# Copyright (c) 2018 Massachusetts Institute of Technology
#
# Authors: Justin Li, Juha Vierinen, Bill Rideout, and Shun-Rong Zhang
#
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


import glob
import time
import warnings
import os
import traceback
import sys
from optparse import OptionParser


# third party imports
from mpl_toolkits.basemap import Basemap
import numpy
import matplotlib.pyplot as plt
import scipy.interpolate
import scipy.signal
import pandas as pd
import numpy as np
import tqdm
from scipy.stats import ks_2samp
from scipy.stats import anderson_ksamp

# Framework imports
from skdiscovery.utilities.planetary.map_util import wgs84_distance

# for calculating geographical distance
def geocalc(lat1, lon1, lat2, lon2):

    return wgs84_distance((lat1, lon1), (lat2, lon2))


def get_lp_tec(tvec, vtec_est, window_length=481, polyorder=3):
    """get_lp_tec returns a low pass version of the vertical tec at the same time spacing
    as vtec_est (that is, at the times given by tvec).  If problem, returns None.  Where data
    cannot be low pass filtered, returns numpy.nan values

    Inputs:
        tvec - input time array in float days
        vtec_est - input vertical tec arr, len = len(tvec)
        window_length - number of 15 second intervals to window over. Default is 481 (2 hours)
            Must be odd
        polyorder - order of polynomial fit to window. Default is 3.
    """
    gap_size = 1.0 / (24*4) # tvec is in float days, so gap_size is 15 minutes
    reg_step = 1.0 / (24*60*4) # step size = 15 seconds in float days
    if window_length % 2 == 0:
        raise ValueError('Window len must be odd, not %i' % (window_length))

    #ret_times = None
    ret_values = None # the numpy values to return

    # first task - break line into gap free sections
    indices = [0]

    for i in range(1, len(tvec)):
        if tvec[i] - tvec[i-1] > gap_size:
            indices.append(i)
    indices.append(len(tvec))

    # deal with each continuous segment separately
    for i in range(len(indices) - 1):
        sub_tvec = tvec[indices[i]:indices[i+1]]
        if len(sub_tvec) < 10:
            lp_values = numpy.zeros((len(sub_tvec),), dtype=numpy.float)
            lp_values[:] = numpy.nan
            if ret_values is None:
                ret_values = lp_values
            else:
                ret_values = numpy.concatenate((ret_values, lp_values))
            continue
        sub_vtec_est = vtec_est[indices[i]:indices[i+1]]
        # create regular time array covering the same span
        reg_tvec = numpy.arange(sub_tvec[0], sub_tvec[-1], reg_step)

        reg_func = scipy.interpolate.interp1d(sub_tvec, sub_vtec_est)
        reg_values = reg_func(reg_tvec)

        # now apply window filter
        if window_length > len(reg_values):
            this_window = len(reg_values)
            if this_window % 2 == 0:
                this_window -= 1
        else:
            this_window = window_length
        lp_reg_values = scipy.signal.savgol_filter(reg_values, this_window, polyorder)

        # now interp windowed data
        lp_func = scipy.interpolate.interp1d(reg_tvec, lp_reg_values)
        for i in range(10):
            try:
                lp_values = lp_func(sub_tvec[:-1 + (-1*i)])
                break
            except ValueError:
                # we may be beyond the interpolation area, drop last point
                continue

        # set values to nan if window was not long enough
        if window_length > len(reg_values):
            lp_values[:] = numpy.nan

        # pad with nan if needed
        if len(lp_values) < len(sub_tvec):
            pad_arr = numpy.zeros((len(sub_tvec)-len(lp_values),), dtype=numpy.float)
            pad_arr[:] = numpy.nan
            lp_values = numpy.concatenate((lp_values, pad_arr))

        if ret_values is None:
            ret_values = lp_values
        else:
            ret_values = numpy.concatenate((ret_values, lp_values))

    if ret_values is None:
        return(None)

    # make sure the length is right by padding begining and end with nan
    while(len(ret_values) < len(tvec)):
        ret_values = numpy.insert(ret_values, 0, numpy.nan)
        if len(ret_values) >= len(tvec):
            break
        ret_values = numpy.insert(ret_values, len(ret_values), numpy.nan)

    return(ret_values)





# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# a bunch of my functions
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

def getRawStitch(DOYs,llat,ulat,llon,rlon,year=2016):
    anasite = []
    anarray = []
    for adoy in DOYs:
        sites,tecs = getRawTEC(adoy,llat,ulat,llon,rlon,year)
        anasite.append(sites)
        anarray.append(tecs)

    catsite = anasite[-1].copy()
    for ii in np.arange(len(anasite)-1):
        appsites = []
        for asite in anasite[ii].index:
            appsites.append(asite)
        catsite.append(anasite[ii].loc[appsites])

    return catsite,pd.concat(anarray)


def fixTECoffset(siteprnTEC,doyN,dchk=3,dcut=.25,mjump=1):
    indices = list(np.hstack((0,np.where((np.abs(np.diff(siteprnTEC.vtec.values))>dcut)==1)[0]+1,len(siteprnTEC)-1)))
    for kk in np.arange(-dchk,dchk+1):
        atmp = np.argmin(np.abs(siteprnTEC.index.values-doyN-kk))
        if atmp not in indices:
            indices.append(atmp)
    indices.sort()

    itemp = indices.copy()
    for ii in np.arange(1,len(itemp[1:-1])):
        aidx = indices[ii]
        if (siteprnTEC.index.values[aidx]-siteprnTEC.index.values[aidx-1])<(mjump/60/24):
            siteprnTEC.vtec.values[indices[ii-1]:indices[ii]] += siteprnTEC.vtec.values[indices[ii]] - siteprnTEC.vtec.values[indices[ii]-1]
            indices[ii] = indices[ii-1]

# # need to check in between days to update index?
# dayindx = []
# for kk in np.arange(-dchk,dchk+1):
#     atmp = np.argmin(np.abs(siteprnTEC.index.values-4-kk))
#     if atmp not in indices:
#         dayindx.append(atmp)
#         indices.append(atmp)

#     elif aidx in dayindx:
#         print(siteprnTEC.vtec.values[indices[ii]] - siteprnTEC.vtec.values[indices[ii]-1])
#         indices[ii] = indices[ii-1]


def findTECevents(rawdata,dayNum,hrEvent,pwin=200,nstd=10,thrstd=.75,verbose=False, fixOffset=False):
    if len(rawdata)==2:
        raw_1 = rawdata[0]
        raw_2 = rawdata[1]
        raw_3 = None
    if len(rawdata)==3:
        # will code actual use of this later maybe
        raw_3 = rawdata[2]

    resbuf = []
    for asite in tqdm.tqdm(np.sort(raw_2[0].index.values)):
        if asite in raw_1[0].index.values:
            sitedata = raw_2[1].iloc[(raw_2[1].site.values==asite),:]
            p_sitedata = raw_1[1].iloc[(raw_1[1].site.values==asite),:]
            p_prns = np.unique(p_sitedata.prn.values)
            for aprn in np.unique(sitedata.prn.values):
                if aprn not in p_prns:
                    if verbose:
                        print('Previous day missing PRN'+str(aprn))
                    pass
                else:
                    pairdat = sitedata.iloc[sitedata.prn.values==aprn,:]
                    p_pairdat = p_sitedata.iloc[p_sitedata.prn.values==aprn,:]
                    if fixOffset:
                        fixTECoffset(pairdat,1)
                        fixTECoffset(p_pairdat,1)
                    if np.argmin(np.abs((pairdat.index.values-(dayNum-1))*24-hrEvent))<=0 or np.argmin(np.abs((pairdat.index.values-(dayNum-1))*24-hrEvent)) >= len(pairdat) or np.argmin(np.abs((p_pairdat.index.values-(dayNum-2))*24-hrEvent))<=0 or np.argmin(np.abs((p_pairdat.index.values-(dayNum-2))*24-hrEvent)) >= len(p_pairdat):
                        pass
                    else:
                        dtec = pairdat.vtec.values - get_lp_tec((pairdat.index.values-(dayNum-1))*24,pairdat.vtec.values)
                        p_dtec = p_pairdat.vtec.values - get_lp_tec((p_pairdat.index.values-(dayNum-2))*24,p_pairdat.vtec.values)

                        tidx = np.argmin(np.abs((pairdat.index.values-(dayNum-1))*24-hrEvent))-int(pwin/2)
                        p_tidx = np.argmin(np.abs((p_pairdat.index.values-(dayNum-2))*24-hrEvent))-int(pwin/2)
                        # zero out nan's and make the other half zero too
                        dtec[np.isnan(dtec)] = 0
                        p_dtec[np.isnan(p_dtec)] = 0
                        # okay to do by multiply because same fs and time aligned via tidx
                        dayLen = len(dtec[tidx:tidx+pwin]); p_dLen = len(p_dtec[p_tidx:p_tidx+pwin])
                        if dayLen > p_dLen:
                            p_dtec = np.hstack((p_dtec,np.zeros(dayLen-p_dLen)))
                        elif p_dLen > dayLen:
                            dtec = np.hstack((dtec,np.zeros(p_dLen-dayLen)))

                        try:
                            nanmask = (dtec[tidx:tidx+pwin]!=0)*(p_dtec[p_tidx:p_tidx+pwin]!=0)
                            dtec[tidx:tidx+pwin] *= nanmask
                            p_dtec[p_tidx:p_tidx+pwin] *= nanmask
                            siglvl = anderson_ksamp([dtec[tidx:tidx+pwin],p_dtec[p_tidx:p_tidx+pwin]]).significance_level
                            pval = ks_2samp(dtec[tidx:tidx+pwin],p_dtec[p_tidx:p_tidx+pwin]).pvalue
                        except:
                            siglvl = 1; pval = 1

                        if siglvl == 1 and pval == 1:
                            if verbose: print('No Matching Data')
                            pass
                        elif np.sum(pairdat.iloc[tidx:tidx+pwin].vtec.rolling(nstd).std().values>thrstd)>nstd or np.sum(p_pairdat.iloc[p_tidx:p_tidx+pwin].vtec.rolling(nstd).std().values>thrstd)>nstd:
                            if verbose: print('Too Noisy for PRN'+str(aprn))
                            pass
                        elif np.abs((pairdat.iloc[max([tidx,0])].name-(dayNum-1))*24-(p_pairdat.iloc[max([p_tidx,0])].name-(dayNum-2))*24)>.5:
                            # if start is off by more than 1/2 hour, pass
                            if verbose: print('Not Time Aligned')
                            pass
                        elif (pairdat.iloc[min([tidx+pwin,len(pairdat)-1])].name-(dayNum-1))*24<hrEvent:
                            # if day of event, data ends before the event, pass
                            if verbose: print('Not Within the Event Time')
                            pass
                        elif np.max(np.diff((pairdat.iloc[tidx:tidx+pwin].index.values-(dayNum-1))*24))>.25 or np.max(np.diff((p_pairdat.iloc[p_tidx:p_tidx+pwin].index.values-(dayNum-2))*24))>.25:
                            # if big (>15min) gap in middle of chunk (ie if because start after event)
                            if verbose:print('Missing Data')
                            pass
                        elif siglvl<.05 and pval<.05:
                            if verbose: print(asite,aprn)
                            resbuf.append([(pairdat.index.values-(dayNum-1))*24,(p_pairdat.index.values-(dayNum-2))*24,dtec,p_dtec,
                                          asite,aprn])
        else:
            if verbose: print('Previous day missing '+asite)
            pass

    print('Found Events for '+str(len(resbuf))+' Site-PRN Pairs')
    return resbuf


def plotTECres(pidx,resbuf,hrEvent,pwin=200):
    tidx = np.argmin(np.abs(resbuf[pidx][0]-hrEvent))-int(pwin/2);
    p_tidx = np.argmin(np.abs(resbuf[pidx][1]-hrEvent))-int(pwin/2)
    plt.plot(resbuf[pidx][0][tidx:tidx+pwin],resbuf[pidx][2][tidx:tidx+pwin]);
    plt.plot(resbuf[pidx][1][p_tidx:p_tidx+pwin],resbuf[pidx][3][p_tidx:p_tidx+pwin]+.15);
    plt.title(resbuf[pidx][4]+'-'+str(resbuf[pidx][5]));
    plt.xlabel('Time (Hrs)'); plt.ylabel('dTEC (TEC Units)');



def makeMap(lat_0,lon_0,dbuffer=5,projection='gnom',resolution='i'):
    # broader region of event
    llat = lat_0-dbuffer
    ulat = lat_0+dbuffer
    llon = lon_0-dbuffer
    rlon = lon_0+dbuffer

    bmap = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0, urcrnrlat=ulat, urcrnrlon=rlon,
                   llcrnrlat=llat,llcrnrlon=llon, resolution=resolution)
    latlons = [llat,ulat,llon,rlon]

    return bmap, latlons


def findPRNs(raw_tec,eventHr,doyN,lat_0,lon_0,latWin=5,lonWin=5,nThreshold=1000):
    potentialPRN = []
    for aprn in np.arange(33):
        site_day = raw_tec[1].iloc[(raw_tec[1].prn.values==aprn),:]
        site_day = site_day.iloc[((site_day.index.values-(doyN-1))*24>eventHr)*((site_day.index.values-(doyN-1))*24<(eventHr+1)),:]
        site_day = site_day.iloc[(np.abs(site_day.pplat.values-lat_0)<latWin)*(np.abs(site_day.pplon.values-lon_0)<lonWin),:]
        if len(site_day)>nThreshold:
            potentialPRN.append(aprn)
    return potentialPRN


def genDTecs(aprn,raw_tec,doyN):
    site_day = raw_tec[1].iloc[(raw_tec[1].prn.values==aprn),:]
    dtecDat = site_day.assign(dtec=0)

    for asite in np.unique(site_day.site.values):
        prdat = site_day.iloc[site_day.site.values==asite]
        dtec = prdat.vtec - get_lp_tec((prdat.index.values-(doyN-1))*24,prdat.vtec.values)
        dtecDat.loc[dtecDat.site.values==asite,'dtec'] = dtec
    return dtecDat


def plotPRNd(raw_tec,dtecDat,eventHr,doyN,lat_0,lon_0,m,fsize=(10,10),clim=.1,ms=5):
    plt.figure(figsize=fsize);m.drawcoastlines();m.plot(lon_0,lat_0,'r^',latlon=True,markersize=15);
    m.plot(raw_tec[0].rec_lon.values,raw_tec[0].rec_lat.values,'bs',latlon=True)
    site_day = dtecDat.iloc[((dtecDat.index.values-(doyN-1))*24>(eventHr-.5))*((dtecDat.index.values-(doyN-1))*24<(eventHr+1.5)),:]
    m.scatter(site_day.pplon.values,site_day.pplat.values,c=site_day.dtec.values,latlon=True,edgecolor='none',
              s=ms,vmin=-clim,vmax=clim);
    parallels = np.arange(-80,80,5.);meridians = np.arange(0.,360.,5.);
    m.drawparallels(parallels,labels=[False,True,True,False])
    m.drawmeridians(meridians,labels=[True,False,False,True])


def plotTracks(prns,asite,raw_tec,eventHr,doyN,lat_0,lon_0,m,fsize=(10,10),ms=[15,5,8]):
    #prns = array of prn's, e.g. [20,21,29]
    plt.figure(figsize=fsize);m.drawcoastlines();m.plot(lon_0,lat_0,'r*',latlon=True,markersize=ms[0]);
    m.plot(raw_tec[0].rec_lon.values,raw_tec[0].rec_lat.values,'w^',latlon=True,markersize=ms[1]);
    m.plot(D_Kaikoura2016[1][0].loc[asite].rec_lon,D_Kaikoura2016[1][0].loc[asite].rec_lat,'b^',latlon=True,markersize=ms[2]);
    for aprn in prns:
        site_day = raw_tec[1].iloc[(raw_tec[1].prn.values==aprn),:]
        dtecDat = site_day.assign(dtec=0)
        prdat = site_day.iloc[site_day.site.values==asite]
        dtec = prdat.vtec - utec.get_lp_tec((prdat.index.values-(doyN-1))*24,prdat.vtec.values)
        dtecDat.loc[dtecDat.site.values==asite,'dtec'] = dtec
        site_day = dtecDat.iloc[(dtecDat.prn.values==aprn)*(dtecDat.site.values==asite),:]
        site_day = site_day.iloc[((site_day.index.values-(doyN-1))*24>(eventHr-.5))*((site_day.index.values-(doyN-1))*24<(eventHr+1.5)),:]
        m.scatter(site_day.pplon.values,site_day.pplat.values,c=site_day.dtec.values,latlon=True,edgecolor='none',s=5,vmin=-.1,vmax=.1);
        t_xy = m(site_day.pplon.values[0],site_day.pplat.values[0])
        plt.text(t_xy[0],t_xy[1],'PRN '+str(aprn));
    parallels = np.arange(-80,80,5.);meridians = np.arange(0.,360.,5.);
    m.drawparallels(parallels,labels=[False,True,True,False])
    m.drawmeridians(meridians,labels=[True,False,False,True])


def genHodochron(raw_data,aprn,doyN,lat_0,lon_0):
    dists_m1 = dict(); times_m1 = dict(); vtecs_m1 = dict(); dtecs_m1 = dict()
    for asite in list(raw_data[0].index):
        gtemp = raw_data[1].iloc[(raw_data[1].prn.values==aprn)*(raw_data[1].site.values==asite),:]
        gtdis = []
        for ii in np.arange(len(gtemp)):
            gtdis.append(geocalc(gtemp.iloc[ii].pplat,gtemp.iloc[ii].pplon,lat_0,lon_0))
        dsigns = (gtemp.pplat.values>lat_0)*2-1
        dists_m1[asite] = gtdis*dsigns
        times_m1[asite] = (gtemp.index.values-(doyN-1))*24
        vtecs_m1[asite] = gtemp.vtec.values
        dtecs_m1[asite] = gtemp.vtec.values - get_lp_tec((gtemp.index.values-(doyN-1))*24,gtemp.vtec.values)
    return dists_m1,times_m1,vtecs_m1,dtecs_m1


def plotHodochron(genRes,eventTime,propTime=None,ylim=[-1500,1500],clim=.1,figsize=(12,5),ms=5,nDir=True,fntsize=10):
    dists_m1,times_m1,vtecs_m1,dtecs_m1 = genRes

    plt.figure(figsize=figsize);
    for asite in list(times_m1.keys()):
        if nDir:
            plt.scatter(times_m1[asite],np.abs(dists_m1[asite]),c=dtecs_m1[asite],edgecolor='none',s=ms);
        else:
            plt.scatter(times_m1[asite],dists_m1[asite],c=dtecs_m1[asite],edgecolor='none',s=ms);
        plt.clim(-clim,clim);
    plt.plot([eventTime,eventTime],[ylim[0]-10,ylim[1]+10],'r-');
    if propTime:
        plt.plot([eventTime+propTime/60,eventTime+propTime/60],[ylim[0]-10,ylim[1]+10],'m--');
    plt.xlim(np.round(eventTime,decimals=1)-.2,np.round(eventTime,decimals=1)+1);
    plt.ylim(ylim[0],ylim[1]);cbar=plt.colorbar();
    plt.xlabel('Time (Hr)',fontsize=fntsize); plt.ylabel('Distance (km)',fontsize=fntsize);
    cbar.ax.set_title('dTEC',fontsize=fntsize);
