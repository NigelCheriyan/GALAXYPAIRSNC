# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 20:12:04 2023

@author: nigel
"""

import csv 
import aplpy
import numpy as np
from astropy import units as u
from astropy.io.votable import parse
from astropy.io import fits
import pandas as pd
from astropy import coordinates as coord
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import binned_statistic
from scipy.stats import norm
import matplotlib.mlab as mlab
import seaborn as sns
sns.set(style="darkgrid")



#Open TOPCAT   matches
with open('MatchT', newline='') as csvfile:
     reader = csv.DictReader(csvfile)
     file = pd.DataFrame(reader)

#Get all files that are a part of a group  ( isolated galaxies are all noted as have '0')
GroupsOnly  = file[file['GroupID_2'] != '0']


#Open Characteristic file for all galaxy groups
file2 = fits.open('G3CFoFGroupv10.fits')
data2 = pd.DataFrame(file2[1].data)
data2 = data2.set_index('GroupID')
#Open Characteristics of each galaxy   
file8 = fits.open('gkvInputCatv02.fits')
data8 = pd.DataFrame(file8[1].data)
data8 = data8.set_index('CATAID')
#Open Galaxy Stellar Masses
file9 = fits.open('StellarMassesv19.fits')
data9 = pd.DataFrame(file9[1].data)
data9 = data9.set_index('CATAID')

Grouped  = file.groupby('GroupID_2').indices


#define variables for groups with different memberships(n_2: 2 galaxies in each group etc.)

n_2_group_flux = []
n_2_group_dist = []
n_2_flux_nsum =[]
n_2_fluxerr =[]
n_2 = pd.DataFrame(columns = file.columns )
n_1 =[]
n_1_flux=[]
n_1_fluxerr = []
n_3_flux = []
size  = []
n_2_fe_nsum =[]
n_3_dist =[]
sizediff = 0.25



for n in Grouped:
        if  n== '0':
            for i in Grouped[n]:
                distance1 = coord.Distance(z = float(file.iloc[i]['Z'])).kpc
                CATAID1 = file.iloc[i]['CATAID']
                uberID1= data8.iloc[int(CATAID1)]['uberID']
                flux = (float(file.iloc[i]['flux_int'])*4*np.pi*distance1**2)/10**data9.loc[int(CATAID1)]['logmstar']
                
                gal1w = np.pi*float(file.iloc[i]['maj_axis'])*float(file.iloc[i]['min_axis'])
                gal1g3 = np.pi*data8.iloc[int(CATAID1)]['R100']**2
                size1 = gal1g3 - gal1w
                if flux< 10**9:
                     n_1_flux.append(flux)
                     fluxerr1 = (float(file.iloc[Grouped[n][0]]['flux_int_err'])*4*np.pi*distance1**2)/10**(data9.loc[int(CATAID1)]['logmstar'])
                     n_1_fluxerr.append(fluxerr1)
            
        if len(Grouped[n]) == 2 and data2.loc[int(n)]['Nfof'] == 2 :  
            CATAID1 = file.iloc[Grouped[n][0]]['CATAID']
            CATAID2 = file.iloc[Grouped[n][1]]['CATAID']

            distance1 = coord.Distance(z = float(file.iloc[Grouped[n][0]]['Z']))
            distance2 = coord.Distance(z = float(file.iloc[Grouped[n][1]]['Z']))
            r_t = data2.loc[int(n)]['Rad100']*1000
            
            gal1g3 = np.pi*data8.iloc[int(CATAID1)]['R100']**2
            gal2g3 = np.pi*data8.iloc[int(CATAID2)]['R100']**2
            
            distance1 = distance1.kpc
            distance2 = distance2.kpc
            gal1w = np.pi*float(file.iloc[Grouped[n][0]]['maj_axis'])*float(file.iloc[Grouped[n][0]]['min_axis'])
            gal2w = np.pi*float(file.iloc[Grouped[n][1]]['maj_axis'])*float(file.iloc[Grouped[n][1]]['min_axis'])
            size1 = gal1g3/gal1w
            size2 = gal2g3/gal2w
            size.append(size1)
            size.append(size2)
            
            flux1 = (float(file.iloc[Grouped[n][0]]['flux_int'])*4*np.pi*distance1**2)/10**(data9.loc[int(CATAID1)]['logmstar'])
            flux2 = (float(file.iloc[Grouped[n][1]]['flux_int'])*4*np.pi*distance2**2)/10**(data9.loc[int(CATAID2)]['logmstar'])
            fluxerr1 = (float(file.iloc[Grouped[n][0]]['flux_int_err'])*4*np.pi*distance1**2)/10**(data9.loc[int(CATAID1)]['logmstar'])
            fluxerr2 =  (float(file.iloc[Grouped[n][1]]['flux_int_err'])*4*np.pi*distance2**2)/10**(data9.loc[int(CATAID2)]['logmstar'])
            if size1 < sizediff and  size2<sizediff :
                if   14 >= gal1w >= 12 or 14 >= gal2w >= 12 :
                    print('beam size noted')
                if (flux1+flux2)>10**10:
                    print('flux too big')
                else:
                    n_2_group_flux.append(flux1 + flux2)
                    n_2_fluxerr.append(np.sqrt(fluxerr1**2 + fluxerr2**2))
                    n_2_flux_nsum.append(flux1)
                    n_2_flux_nsum.append(flux2)
                    n_2_group_dist.append(r_t)
                
        if len(Grouped[n]) ==3 and data2.loc[int(n)]['Nfof'] == 3 :
             CATAID1 = file.iloc[Grouped[n][0]]['CATAID']
             CATAID2 = file.iloc[Grouped[n][1]]['CATAID']
             CATAID3 = file.iloc[Grouped[n][2]]['CATAID']
             distance1 = coord.Distance(z = float(file.iloc[Grouped[n][0]]['Z'])).kpc
             distance2 = coord.Distance(z = float(file.iloc[Grouped[n][1]]['Z'])).kpc
             distance3 = coord.Distance(z = float(file.iloc[Grouped[n][2]]['Z'])).kpc
             flux1 = (float(file.iloc[Grouped[n][0]]['flux_int'])*4*np.pi*distance1**2)/10**(data9.loc[int(CATAID1)]['logmstar'])
             flux2 = (float(file.iloc[Grouped[n][1]]['flux_int'])*4*np.pi*distance2**2)/10**(data9.loc[int(CATAID2)]['logmstar'])
             flux3 = (float(file.iloc[Grouped[n][2]]['flux_int'])*4*np.pi*distance3**2)/10**(data9.loc[int(CATAID3)]['logmstar'])
             n_3_dist.append(data2.loc[int(n)]['Rad100']*1000)
             n_3_dist.append(data2.loc[int(n)]['Rad100']*1000)
             n_3_dist.append(data2.loc[int(n)]['Rad100']*1000)
             n_3_flux.append(flux1)
             n_3_flux.append(flux2)
             n_3_flux.append(flux3)
                                        
plt.figure()   
bins= np.linspace(0, 4, 80)   
sns.histplot(pd.to_numeric(file['Separation']).values, bins = bins , color = 'skyblue')
plt.title('SwagX to GAMA Catalogue Separation')
plt.xlabel('Separation (")')


plt.figure()
n_2 = pd.DataFrame({'Dist':n_2_group_dist,'Flux':n_2_group_flux})

bin_means, bin_edges, binnumber =  binned_statistic( n_2['Dist'],n_2['Flux'], bins=np.logspace(0,4,num = 20))
plt.scatter(n_2_group_dist,n_2_group_flux,color= 'skyblue')

plt.errorbar(n_2_group_dist,n_2_group_flux, yerr = n_2_fluxerr, fmt = 'o',markersize=4, capsize = 0.1, color = 'skyblue')

bin_means_err , bin_null, binnumber1 =  binned_statistic( n_2['Dist'],n_2_fluxerr, bins=np.logspace(0,4,num = 20))

plt.hlines(y = bin_means, xmin = bin_edges[:-1], xmax = bin_edges[1:], colors='g', lw= bin_means_err ,label='error', alpha = 0.5)
plt.hlines(y = bin_means, xmin = bin_edges[:-1], xmax = bin_edges[1:], colors='g', lw=2,
           label='Mean of Paired Galaxies')
plt.ylim([0,1500])
plt.xlim([0,1500])
plt.title('Paired Galaxy Flux Denisty across Large Group Radius')
plt.xlabel('Distance(Kpc)')
plt.ylabel('Luminosity(MJy/solarmass)')



plt.title('Paired Galaxy Luminosity Denisty ')
plt.yscale('log')
plt.xlabel('Distance(Kpc)')
plt.ylabel('Flux(MJy/log(solarmass))')


plt.figure()
sns.histplot(size,color = 'skyblue')
plt.xlim([0,2])
plt.title('Size Ratio from Optical to Radio Source')
plt.xlabel('Optical source/Radio Source')


plt.figure()
sns.histplot(n_2_group_dist, bins= np.linspace(min(n_2_group_dist), max(n_2_group_dist) , 75))
plt.xlabel('Separation Distance (Kpc)')

 
plt.figure()
bins= np.linspace(0, 4, 100)
(mu, sigma) = norm.fit(n_2_flux_nsum)
y = norm.pdf( bins, mu, sigma)
plt.plot(bins, y,label = 'Galaxy Pair')
(mu, sigma) = norm.fit(n_1_flux)
y = norm.pdf( bins, mu, sigma)
plt.plot(bins, y, label = 'Isolated Galaxy' )
(mu, sigma) = norm.fit(n_3_flux)
y = norm.pdf( bins, mu, sigma)
plt.plot(bins, y, label = 'Galaxy Trio' )
plt.xlabel('Luminosity(MJy/(solarmass)')
plt.title('Luminosity for Isolated and Paired Galaxies')
plt.ylabel('Density')
plt.legend()

plt.figure()
n_2 = pd.DataFrame({'Dist':n_2_group_dist,'Flux':n_2_group_flux})
bin_means, bin_edges, binnumber =  binned_statistic( n_2['Dist'],n_2['Flux'], bins=np.logspace(0,3,num = 12), statistic = 'median')
bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2
plt.errorbar(n_2_group_dist,n_2_group_flux, yerr = n_2_fluxerr, fmt = 'o',markersize=4, capsize = 4, color = 'skyblue')
bin_means_err , bin_null, binnumber1 =  binned_statistic( n_2['Dist'],n_2['Flux'], bins=np.logspace(0,3,num = 12), statistic = lambda x:np.median(np.absolute(x - np.median(x)) ))
bin_centers =np.array( [sum(x) for x in zip(bin_edges[1:].astype('float'),np.logspace(0,3,num = 12))]).astype('float')/2
plt.errorbar(bin_centers,bin_means, yerr = bin_means_err, fmt = 'o', capsize = 4, color = 'green')
plt.hlines(y = bin_means, xmin = bin_edges[:-1], xmax = bin_edges[1:], colors='g', lw=2,label='Median of Paired Galaxies')
mad = np.median(np.absolute(n_1_flux - np.median(n_1_flux)) )
plt.errorbar(75 ,np.median(n_1_flux) , yerr = mad, fmt = 'o', capsize = 4, color = 'red')
plt.hlines(y = np.median(n_1_flux), xmin = 0 , xmax = 1500,lw = 1, colors = 'r', label = 'Median of Isolated Galaxies')
plt.ylim([0,1000])
plt.xlim([0,200])
plt.title('Paired Galaxy Luminosity')
plt.xlabel('Distance(kpc)')
plt.ylabel('Luminosity(MJy/solarmass)')
plt.legend()



plt.figure()
n_2 = pd.DataFrame({'Dist':n_2_group_dist,'Flux':n_2_group_flux})
bin_means, bin_edges, binnumber =  binned_statistic( n_2['Dist'],n_2['Flux'], bins=np.logspace(0,3,num = 12), statistic = 'median')
bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2
plt.errorbar(n_2_group_dist,n_2_group_flux, yerr = n_2_fluxerr, fmt = 'o',markersize=4, capsize = 4, color = 'skyblue')
bin_means_err , bin_null, binnumber1 =  binned_statistic( n_2['Dist'],n_2['Flux'], bins=np.logspace(0,3,num = 12), statistic = lambda x:np.median(np.absolute(x - np.median(x)) ))
bin_centers =np.array( [sum(x) for x in zip(bin_edges[1:].astype('float'),np.logspace(0,3,num = 12))]).astype('float')/2
plt.errorbar(bin_centers,bin_means, yerr = bin_means_err, fmt = 'o', capsize = 4, color = 'green')
plt.hlines(y = bin_means, xmin = bin_edges[:-1], xmax = bin_edges[1:], colors='g', lw=2,label='Median of Paired Galaxies')
mad = np.median(np.absolute(n_1_flux - np.median(n_1_flux)) )
plt.errorbar(75 ,np.median(n_1_flux) , yerr = mad, fmt = 'o', capsize = 4, color = 'red')
plt.hlines(y = np.median(n_1_flux), xmin = 0 , xmax = 1500,lw = 1, colors = 'r', label = 'Median of Isolated Galaxies')
plt.yscale('log')
plt.xscale('log')
plt.title('Paired Galaxy Luminosity')
plt.xlabel('Distance(kpc)')
plt.ylabel('Luminosity(MJy/solarmass)')
plt.legend()


plt.figure()
bins= np.linspace(0, 1e13, 100)
sns.histplot(n_1_flux,label = 'Isolated Galaxy',bins =bins ,color="skyblue",stat = 'density' )
sns.histplot(n_2_flux_nsum,label = 'Galaxy Pair',bins = bins ,color="red",stat= 'density')
sns.histplot(n_3_flux,label = 'Galaxy Trio',bins = bins ,color="orange",stat= 'density')
plt.legend()

plt.figure()
plt.scatter(n_3_dist,n_3_flux,color= 'skyblue')

