import numpy as np
from scipy import stats
import scipy.constants as c
import scipy.special as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors

M = 1000 #kg
T = 10 #yr

T2v = 1.926e21
T0vS = 1e3*c.N_A*M*T*np.log(2)*1/76
Qbb = 2039
alpha = 1 - sp.erf(3/np.sqrt(2))
beta = 0.5
eff_tot = 0.8

B = 4.7e-3*M*T*4.14

x = np.arange(0,max(10,B + 7*np.sqrt(B)), 1)
h0 = stats.poisson.pmf(x,B)

i = 0
integral = 0

while integral < 1 - alpha:
    integral += h0[i]
    i += 1

S = 0
integral2 = np.sum(stats.poisson.pmf(x,B+S)[:i])




while integral2 > beta:
    S += 0.1
    integral2 = np.sum(stats.poisson.pmf(x,B+S)[:i])
    
h1 = stats.poisson.pmf(x,B+S)

plt.figure(figsize = (15,5))
plt.plot(x,h0)
plt.plot(x,h1)

Thalf = T0vS/S

exp = np.arange(-8,-2, 0.1)
mu_barray = 10**exp 
sigarray = np.arange(0.1,10,0.1)


Thalfarray = np.zeros((len(sigarray), len(exp)))
#Thalfarray = np.zeros(len(exp))
Sarray = np.zeros(len(exp))
Barray = np.zeros(len(exp))
narray = np.zeros(len(exp))
iarray = np.zeros(len(exp))


# for isig, sigma in enumerate(sigarray):
#     for imu_b, mu_b in enumerate(mu_barray):
#         B2v = (T0vS/T2v)*(sigma/Qbb)**6
#         n = 0.5*(2.8+10**(-0.48 - 0.32*np.log10(mu_b*sigma*M*T) - 0.046*(np.log10(mu_b*sigma*M*T))**2))
#         aB = mu_b*M*T*2*n*sigma
#         Barray[imu_b] = aB + B2v 
#         narray[imu_b] = n
        
#         eff = sp.erf(n/np.sqrt(2))
        
#         x = np.arange(0, max(10,aB+B2v + 5*np.sqrt(aB+B2v)), 1)
#         h0 = stats.poisson.pmf(x,aB + B2v)
        
#         i = 0
#         integral = 0
        
#         while integral < 1 - alpha:
#             integral += h0[i]
#             i += 1
#         iarray[imu_b] = i
        
#         s = 0
#         integral2 = np.sum(stats.poisson.pmf(x,aB + B2v + s)[:i])
#         while integral2 > beta:
#             s += 0.1
#             integral2 = np.sum(stats.poisson.pmf(x,aB + B2v + s)[:i])

        
#         Thalfarray[isig,imu_b] = T0vS*eff_tot*eff/s
#         #Thalfarray[imu_b] = T0vS*eff/s
#         Sarray[imu_b] = s

#plt.figure(figsize = (15,5))
#plt.plot(sigarray,Thalfarray[:,1])
##plt.xscale('log')
#
#plt.figure(figsize = (15,5))
#plt.plot(Barray,Sarray)
#plt.xscale('log')
#plt.yscale('log')
#
#plt.figure(figsize = (15,5))
#plt.plot(Barray,iarray)
#plt.xscale('log')
#plt.yscale('log')

import gammapy.stats as gstats

exp = np.arange(-8,-4, 0.1)
mu_barray = 10**exp 
sigarray = np.arange(1,9.5,0.5)
#Thalfarray = np.zeros(len(exp))
Thalfarray = np.zeros((len(sigarray), len(exp)))
Sarray = np.zeros(len(exp))
Barray = np.zeros(len(exp))

for isig, sigma in enumerate(sigarray):
    for imu_b, mu_b in enumerate(mu_barray):
        B2v = (T0vS/T2v)*(sigma/Qbb)**6
        n = 0.5*(2.8+10**(-0.48 - 0.32*np.log10(mu_b*sigma*M*T) - 0.046*(np.log10(mu_b*sigma*M*T))**2))
        aB = mu_b*M*T*2*n*sigma
        Barray[imu_b] = aB + B2v 
        narray[imu_b] = n
        
        eff = sp.erf(n/np.sqrt(2))
        
        #Brand = np.random.poisson( aB + B2v, 1000)
        x_bins = np.arange(0, 100)
        mu_bins = np.linspace(0, 50, int(50 / 0.1) + 1, endpoint=True)
        matrix = [stats.poisson(mu + aB + B2v).pmf(x_bins) for mu in mu_bins]
        acceptance_intervals = gstats.fc_construct_acceptance_intervals_pdfs(matrix, 0.9)
        LowerLimitNum, UpperLimitNum, _ = gstats.fc_get_limits(mu_bins, x_bins, acceptance_intervals)
        sarr = np.zeros(1000)
        for ex in np.arange(1,1000,1):
            sarr[ex] = gstats.fc_find_limit(np.random.poisson(aB + B2v), UpperLimitNum, mu_bins)
        s = np.median(sarr)
        Thalfarray[isig,imu_b] = T0vS*eff_tot*eff/s
        #Thalfarray[imu_b] = T0vS*eff/s
        Sarray[imu_b] = s

#matplotlib.rcParams['font.sans-serif'] = "Helvetica"
font = {'weight' : 'normal', #'family' : 'sans-serif'
        'size'   : 15}
matplotlib.rc('font', **font)

# plt.figure(figsize = (15,5))
# plt.plot(mu_barray,Sarray)
# plt.xscale('log')

plt.figure(figsize = (10,5))
plt.pcolormesh(mu_barray,2.355*sigarray,Thalfarray, cmap = 'inferno', norm=colors.LogNorm(vmin=Thalfarray.min(), vmax=Thalfarray.max()))
cbar = plt.colorbar()
#cbar.set_label('$S_{1/2}^{0\\nu(disc)}$', labelpad=-10, y=1.00, rotation=0)
cbar.set_label('$S_{1/2}^{0\\nu}$', labelpad=-35, y=1.00, rotation=0)
plt.xscale('log')
plt.xlabel("Background Rate [c/(keV-kg-yr)]")
plt.ylabel("Resolution FWHM [keV]") 
plt.savefig('sensitivity_small.png', dpi=200, bbox_inches='tight')

#
##plt.figure(figsize = (15,5))
##plt.plot(mu_barray,Barray)
##plt.xscale('log')
##plt.yscale('log')
##
##plt.figure(figsize = (15,5))
##plt.plot(mu_barray,narray)
##plt.xscale('log')
#
