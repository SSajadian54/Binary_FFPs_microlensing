import numpy as np 
import pylab as py 
import matplotlib.pyplot as plt 
from matplotlib import rcParams
from numpy import ma
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
from matplotlib import gridspec
import matplotlib.dates as mdates
from matplotlib.transforms import Transform
from matplotlib.ticker import (AutoLocator, AutoMinorLocator)
import sys

#######################################################################3
def func1(q): 
    dw=0.0;  dclo=0.0
    dw=np.power(1.0 + pow(q,float(1.0/3.0)),1.5)/np.sqrt(1.0+q)
    min1=1000000000.0;  
    for w in range(500): 
        dc= float(0.6+w*(0.4/500.0));
        tt=abs(pow(dc,8.0)*27.0*q - (1.0+q)*(1.0+q)*pow(1.0 - dc*dc*dc*dc , 3.0))
        if(tt<min1):  
            min1=tt;  
            dclo=dc
    return(dw , dclo)        
#######################################################################3
f1=open("./Moon2.txt","r")
nm= sum(1 for line in f1) 
res= np.zeros((nm , 4))
res=np.loadtxt("./Moon2.txt")
res[:,2]=np.power(10.0,res[:,2])
res[:,3]=np.power(10.0,res[:,3])
res[:,1]=np.power(10.0,res[:,1])
#######################################################################

logq=np.arange(-13.0,0.0, 0.05)
q=np.power(10.0,logq)
nq=len(q)
dwide=np.zeros(nq)
dclos=np.zeros(nq)
for i in range(nq): 
    dwide[i], dclos[i]= func1(q[i])
#######################################################################
lab=[r"$\rm{Earth}-\rm{Moon}$", r"$\rm{Mars}-\rm{Phobos}$" , r"$\rm{Mars}-\rm{Deimos}$",  
r"$\rm{Jupiter}-\rm{Europa}$", r"$\rm{Jupiter}-\rm{Ganymede}$",r"$\rm{Jupiter}-\rm{Callisto}$", r"$\rm{Jupiter}-\rm{IO}$", 
#r"$\rm{Jupiter}-\rm{Carpo}$",
#r"$\rm{Jupiter}-\rm{Adrastea}$", 
#r"$\rm{Jupiter}-\rm{Amalthea}$",   
r"$\rm{Saturn}$-\rm{Titan}",r"$\rm{Saturn}$-\rm{Rhea}", r"$\rm{Saturn}$-\rm{Iapetus}", r"$\rm{Saturn}$-\rm{Pione}",
r"$\rm{Saturn}$-\rm{Tethys}",r"$\rm{Saturn}$-\rm{Enceladus}", r"$\rm{Saturn}$-\rm{Mimas}",
r"$\rm{Uranus}-\rm{Titania}$",r"$\rm{Uranus}-\rm{Oberon}$",r"$\rm{Uranus}-\rm{Umbriel}$", r"$\rm{Uranus}-\rm{Ariel}$",r"$\rm{Uranus}-\rm{Miranda}$", 
r"$\rm{Neptune}-\rm{Triton}$", r"$\rm{Neptune}-\rm{Proteus}$", 
r"$\rm{Pluto}-\rm{Charon}$",r"$\rm{Pluto}-\rm{Hydra}$",
#r"$\rm{Saturn}$-\rm{Janus}", 
#r"$\rm{Saturn}$-\rm{Hyperion}",#r"$\rm{Saturn}$-\rm{Epimetheus}", r"$\rm{Saturn}$-\rm{Pheobe}",
r"$\rm{Neptune}$-\rm{Nereid}"]


plt.clf()
fig= plt.figure(figsize=(8,6))
ax = fig.add_subplot()
plt.plot(dwide, q, "r--", lw=2.0, label= r"$s_{\rm w}(q)$")
plt.plot(dclos, q, "b--", lw=2.0, label= r"$s_{\rm c}(q)$")
plt.ylim(2.0e-9, 0.19)
plt.xlim(0.014 , 1.80 )
plt.xscale('log')
plt.yscale('log')
for i in range(nm):
    plt.errorbar(res[i,2],res[i,1], xerr=abs(res[i,3]),fmt='o',color='darkgreen',ecolor='m',capsize=0.0, ms=7.0, lw=1.2)
    if(i==15 or i==17):     
        ax.text(res[i,2]*0.6,res[i,1]*0.6,lab[i],fontsize=12)
    elif(i==7): 
        ax.text(res[i,2]*0.9,res[i,1]*1.15,lab[i],fontsize=12)
    else:
        ax.text(res[i,2]*0.6,res[i,1]*1.1,lab[i],fontsize=12)        
plt.xlabel(r"$s(R_{\rm E})$", fontsize=18)
plt.ylabel(r"$q$",fontsize=18)
fig=plt.gcf()
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True, fontsize=19)
#plt.grid("True")
#plt.grid(linestyle='dashed')
fig.tight_layout()
fig.savefig("fig1b.jpg",dpi=200)
py.clf()






















