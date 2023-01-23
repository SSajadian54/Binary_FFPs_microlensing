import numpy as np 
import pylab as py 
import matplotlib.pyplot as plt 
import pandas as pd 
import seaborn as sns
from matplotlib import rcParams
import time
import matplotlib
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.ticker import StrMethodFormatter
cmap=plt.get_cmap('viridis')
import collections
import scipy.stats as stats
from scipy.optimize import curve_fit
import scipy.optimize as optimize
import matplotlib.dates as mdates
from matplotlib.transforms import Transform
from matplotlib.ticker import (AutoLocator, AutoMinorLocator)
#######################################################

def func(x, a1, a2):
    return( a1*x + a2)
#######################################################    


nm=24

array=np.zeros((nm,3))
array=np.loadtxt("./files/MONT/Moon2.txt")
array[:,1]=np.log10(array[:,1])
array[:,2]=np.log10(array[:,2])


ini1 = np.array([0.25,0.1])
fitt, pcov = curve_fit(func, array[:,1], array[:,2], ini1, maxfev = 200000)
print ("*****************************************************")
print ("Slope:   ",   fitt[0] ,  "+/-:  ",  np.sqrt(pcov[0,0]))
print ("Offset:  ",   fitt[1] ,  "+/-:  ",  np.sqrt(pcov[1,1]))
print ("*****************************************************")    

#up= np.array([ fitt[0]+np.sqrt(pcov[0,0]), fitt[1]+np.sqrt(pcov[1,1]) ])

up2= np.array([ fitt[0]-np.sqrt(pcov[0,0]), fitt[1]+np.sqrt(pcov[1,1]) ])
do= np.array([ fitt[0]+np.sqrt(pcov[0,0]), fitt[1]-np.sqrt(pcov[1,1]) ])
#do2= np.array([ fitt[0]-np.sqrt(pcov[0,0]), fitt[1]-np.sqrt(pcov[1,1]) ])

par=np.array([fitt[0],  np.sqrt(pcov[0,0]) , fitt[1], np.sqrt(pcov[1,1])  ])
qq=np.arange( np.min(array[:,1]) ,   np.max(array[:,1])  ,  0.0005, dtype=float)

plt.clf()
fig= plt.figure(figsize=(8,6)) 
plt.scatter(array[:,1],array[:,2],s=70 ,c="g" , label=r"$\rm{Moon}-\rm{Planet}~\rm{systems}$")

plt.plot(qq,func(qq,*fitt),"m--", lw=2.4, label=r"$\rm{Fit}:\rm{Slope}=$"+str(round(fitt[0],2))+r"$\pm$"+str(round(np.sqrt(pcov[0,0]),2))+r"$,~\rm{Offset}=$"+str(round(fitt[1],2))+r"$\pm$"+str(round(np.sqrt(pcov[1,1]),2)) )

# {0:%3.1f}\pm{1:%3.1f},~Offset={2:%3.1f}\pm{3:%3.1f}$".format(par[0], par[1], par[2], par[3]) )
#plt.plot(qq,func(qq,*up), 'm-', lw=1.0, alpha=0.5)
plt.plot(qq,func(qq,*do), 'm-', lw=1.0, alpha=0.5)
plt.plot(qq,func(qq,*up2), 'm-', lw=1.0, alpha=0.5)
#plt.plot(qq,func(qq,*do2), 'k-', lw=1.0, alpha=0.5)
plt.fill_between(qq, func(qq,*up2), func(qq,*do), color='magenta', alpha=0.3)

plt.xlabel(r"$\log_{10}[q]$", fontsize=18)
plt.ylabel(r"$\log_{10}[d\big/R_{\rm{p}}]$", fontsize=18)
plt.xlim([-9.02,-0.8])
plt.ylim([0.25,2.0])
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.grid(True)
plt.grid(linestyle='dashed')
plt.legend()
plt.legend(loc=4,fancybox=True, shadow=False, framealpha=0.6,frameon=True, prop={"size":17}, labelspacing=0.1)
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./scatter1.jpg")
print(">>>>>>> scatter plot was made <<<<<<<<<<<<<<<")

##########################################################
qq=np.arange( np.min(array[:,1]) ,   np.max(array[:,1])  ,  0.0005, dtype=float)
ideal, isd=np.mean(qq), np.std(qq)
real,  rsd=np.mean(array[:,1]), np.std(array[:,1])
print ("Log10[q]    real average :  ",  real, rsd)
print ("*******     ideal average :  ",  ideal, isd,   np.min(qq),   np.max(qq)  )
print("*********************************************")

plt.clf()
fig= plt.figure(figsize=(8,6)) 
ax= plt.gca()              
plt.hist(array[:,1],10, histtype='bar',facecolor='green', alpha=0.7, rwidth=0.95 )
y_vals = ax.get_yticks()
ax.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/nm)) for x in y_vals]) 
#y_vals = ax.get_yticks()
#plt.ylim([np.min(y_vals), np.max(y_vals)])
plt.xlim(np.min(array[:,1]) , np.max(array[:,1]))
plt.axvline(x=ideal,  color="k", ls='--', lw=2.0)
plt.axvline(x=real,    color="darkgreen", ls='--', lw=2.0)
plt.xlabel(r"$\log_{10}[q]$", fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=18,labelpad=0.1)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.text(-4.0, 6.0, r"$\left<\log_{10}[q]\right>$="+str(round(real,2))+r"$\pm$"+str(round(rsd, 2)),color="darkgreen", fontsize=18)
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./HistQ.jpg")
print(">>>>>>> Histogram 1 was made <<<<<<<<<<<<<<<")




##########################################################
dis=np.arange( np.min(array[:,2]) ,   np.max(array[:,2]) , 0.0005, dtype=float)
ideal, isd=np.mean(dis), np.std(dis)
real,  rsd=np.mean(array[:,2]), np.std(array[:,2])
print ("log10[dis/Rp]  real average :  ",  real, rsd)
print ("*****          ideal average :  ",  ideal, isd, np.min(dis), np.max(dis)     )
print("*********************************************")

plt.clf()
fig= plt.figure(figsize=(8,6)) 
ax= plt.gca()              
plt.hist(array[:,2],8, histtype='bar',facecolor='green', alpha=0.7, rwidth=0.95 )
y_vals = ax.get_yticks()
ax.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/nm)) for x in y_vals]) 
#y_vals = ax.get_yticks()
#plt.ylim([np.min(y_vals), np.max(y_vals)])
plt.xlim(np.min(array[:,2]) , np.max(array[:,2]))

plt.xlabel(r"$\log_{10}[d\big/R_{\rm{p}}]$", fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=18,labelpad=0.1)
plt.axvline(x=ideal,  color="k", ls='--', lw=2.0)
plt.axvline(x=real,    color="darkgreen", ls='--', lw=2.0)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.text(1.26, 3.8, r"$\left<\log_{10}[d/R_{\rm p}]\right>$="+str(round(real,2))+r"$\pm$"+str(round(rsd, 2)),color="darkgreen", fontsize=18)
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./HistD.jpg")
print(">>>>>>> Histogram 2 was made <<<<<<<<<<<<<<<")








