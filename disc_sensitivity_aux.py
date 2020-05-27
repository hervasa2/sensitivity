import numpy as np
from scipy import stats
import scipy.constants as c
import scipy.special as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import gammapy.stats as gstats

M = 1000 #kg
T = 10 #yr

T2v = 1.926e21
T0vS = 1e3*c.N_A*M*T*np.log(2)*1/76
Qbb = 2039
alpha = 1 - sp.erf(3/np.sqrt(2))
beta = 0.5
eff_tot = 0.8

# plt.figure(figsize = (10,5))
# plt.pcolormesh(mu_barray,2.355*sigarray,Thalfarray, cmap = 'inferno', norm=colors.LogNorm(vmin=Thalfarray.min(), vmax=Thalfarray.max()))
# cbar = plt.colorbar()
# #cbar.set_label('$S_{1/2}^{0\\nu(disc)}$', labelpad=-10, y=1.00, rotation=0)
# cbar.set_label('$S_{1/2}^{0\\nu}$', labelpad=-35, y=1.02, rotation=0)
# plt.xscale('log')
# plt.xlabel("Background Rate [c/(keV-kg-yr)]")
# plt.ylabel("Resolution FWHM [keV]") 
# plt.savefig('sensitivity_small.png', dpi=200, bbox_inches='tight')

mu_b = 8e-6
sigarray = np.arange(0.1,10,0.1)
Thalfarray = np.zeros(len(sigarray))
Shalfarray = np.zeros(len(sigarray))
Sarray = np.zeros(len(sigarray))

for isig, sigma in enumerate(sigarray):
    B2v = (T0vS/T2v)*(sigma/Qbb)**6
    n = 0.5*(2.8+10**(-0.48 - 0.32*np.log10(mu_b*sigma*M*T) - 0.046*(np.log10(mu_b*sigma*M*T))**2))
    aB = mu_b*M*T*2*n*sigma
    
    eff = sp.erf(n/np.sqrt(2))
    
    x = np.arange(0, max(10,aB+B2v + 5*np.sqrt(aB+B2v)), 1)
    h0 = stats.poisson.pmf(x,aB + B2v)
    
    i = 0
    integral = 0
    
    while integral < 1 - alpha:
        integral += h0[i]
        i += 1
    
    s = 0
    integral2 = np.sum(stats.poisson.pmf(x,aB + B2v + s)[:i])
    while integral2 > beta:
        s += 0.01
        integral2 = np.sum(stats.poisson.pmf(x,aB + B2v + s)[:i])
    
    Thalfarray[isig] = T0vS*eff_tot*eff/s
    
    x_bins = np.arange(0, 50)
    mu_bins = np.linspace(0, 20, int(20 / 0.01) + 1, endpoint=True)
    matrix = [stats.poisson(mu + aB + B2v).pmf(x_bins) for mu in mu_bins]
    acceptance_intervals = gstats.fc_construct_acceptance_intervals_pdfs(matrix, 0.9)
    LowerLimitNum, UpperLimitNum, _ = gstats.fc_get_limits(mu_bins, x_bins, acceptance_intervals)
    sarr = np.zeros(10000)
    for ex in np.arange(1,10000,1):
        sarr[ex] = gstats.fc_find_limit(np.random.poisson(aB + B2v), UpperLimitNum, mu_bins)
    s = np.median(sarr)
    Shalfarray[isig] = T0vS*eff_tot*eff/s
        
plt.figure(figsize = (5,5))
plt.plot(2.355*sigarray,Shalfarray/Shalfarray[0], label = '$S_{1/2}^{0\\nu}$')
plt.plot(2.355*sigarray,Thalfarray/Thalfarray[0], label = '$S_{1/2}^{0\\nu(disc)}$')
plt.legend()
plt.xlabel("Resolution FWHM [keV]") 
plt.savefig('discvslim.png', dpi=200, bbox_inches='tight')
#plt.xscale('log')
