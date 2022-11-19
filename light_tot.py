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
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib import gridspec
import matplotlib.dates as mdates
from matplotlib.transforms import Transform
from matplotlib.ticker import (AutoLocator, AutoMinorLocator)
import sys
############################################################################################
f1=open("./files/Moon.txt","r")
nm= sum(1 for line in f1) 
arr= np.zeros((nm,54))
arr=np.loadtxt("./files/Moon.txt")
fil1=open('./info.txt',"w");   
fil1.close(); 
#######################################################################3
def func1(q): 
    dw=0.0;  dclo=0.0
    dw=np.power(1.0 + pow(q,float(1.0/3.0)),1.5)/np.sqrt(1.0+q)
    min1=1000000000.0;  
    for w in range(200): 
        dc= float(0.6+w*(0.4/200.0));
        tt=abs(pow(dc,8.0)*27.0*q - (1.0+q)*(1.0+q)*pow(1.0 - dc*dc*dc*dc , 3.0))
        if(tt<min1):  
            min1=tt;  
            dclo=dc
    return(dw , dclo)        
#######################################################################3
ns=int(11)

lab=[r"$\rm{Earth}-\rm{Moon}$",  r"$\rm{Mars}-\rm{Phobos}$" , r"$\rm{Mars}-\rm{Deimos}$",  r"$\rm{Jupiter}-\rm{Ganymede}$", 
     r"$\rm{Jupiter}-\rm{Callisto}$", r"$\rm{Jupiter}-\rm{Europa}$", r"$\rm{Saturn}$-\rm{Titan}",r"$\rm{Uranus}-\rm{Titania}$", 
     r"$\rm{Uranus}-\rm{Oberon}$", r"$\rm{Neptune}-\rm{Triton}$", r"$\rm{Pluto}-\rm{Charon}$"]##11

res=np.zeros((ns, 10))
for j in range(ns):
    

    clos=0.0; inte=0.0; wide=0.0; count=0.000003434512; 
    qave=0.0; dave=0.0; tave=0.0;  pave=0.0;  rave=0.0; roave=0.0;  siz1a=0.0;  siz2a=0.0
    dq= np.zeros((1000000,2));  k=0;
    
    for i in range(nm): 
        icon, lat, lon= int(arr[i,0]), arr[i,1], arr[i,2]
        strucl, logMl, Dl, vl= arr[i,3], arr[i,4], arr[i,5], arr[i,6]
        strucs, cl, Ms, Ds, Tstar, Rstar = arr[i,7], arr[i,8], arr[i,9], arr[i,10], arr[i,11], arr[i,12]
        logl, types, col, vs, MI, MW149, mI, mW149= arr[i,13], arr[i,14], arr[i,15], arr[i,16], arr[i,17], arr[i,18], arr[i,19], arr[i,20]
        magbI, magbW, fbI, fbW, NbI, NbW, ExI, ExW=arr[i,21], arr[i,22], arr[i,23], arr[i,24], arr[i,25], arr[i,26], arr[i,27], arr[i,28]
        tE,RE,t0, mul, Vt, u0= arr[i,29], arr[i,30], arr[i,31], arr[i,32], arr[i,33], arr[i,34]
        opd,ros,tetE=  arr[i,35], arr[i,36], arr[i,37]
        li,mus1, mus2, xi, mul1,mul2,piE=arr[i,38], arr[i,39], arr[i,40], arr[i,41], arr[i,42], arr[i,43], arr[i,44]
        q, dis, peri, siz1, siz2= arr[i,45], arr[i,46], arr[i,47], arr[i,48], arr[i,49]
        con,delta, Asym, case= arr[i,50], arr[i,51], arr[i,52], arr[i,53]

        if(case==j): 
            qave+=q;   
            dave+=np.log10(dis);     
            tave+=tE;   
            pave+=peri;  
            rave+=RE*10.0
            roave+=ros; 
            siz1a+=q*pow(dis-1.0/dis,-2.0)
            siz2a+=np.sqrt(q)*pow(dis,3.0)
         
            dq[k,0]= np.log10(arr[i,46]);### dis 
            dq[k,1]= arr[i,45];### q
            k+=1

            dwid, dclo= func1( arr[i,45] )
            if(  dis<=dclo):   clos+=1.0;
            elif( dis<dwid):   inte+=1.0;
            elif(dis>=dwid):   wide+=1.0;
            count+=1.0; 

    clos=float(clos*100.0/count)    
    inte=float(inte*100.0/count)    
    wide=float(wide*100.0/count)            
    ave=np.array([ np.log10(qave/count), float(dave/count), float(tave/count), float(pave/count), float(rave/count),
                 float(roave/count), np.log10(float(siz1a/count)), np.log10(float(siz2a/count)) ]) 
    res[j,:8]= ave 
    fil1=open('./info.txt',"a+");   
    np.savetxt(fil1,ave.reshape((1,8)),"$m{}$ & $%.3f$ &  $%.3f$ &  $%.2f$ &  $%.2f$ &  $%.3f$ &  $%.3f$ &  $%.4f$ &  $%.4f$\n")
    fil1.close()
    res[j,8]=np.std(dq[:k,0])## dis
    res[j,9]=np.std(dq[:k,1])## log q
    print(res[j,:])
    print("***********************************************")
#######################################################################33
logq=np.arange(-6.0,0.0, 0.1)
q=np.power(10.0,logq)
nq=len(q)
dwide=np.zeros(nq)
dclos=np.zeros(nq)
for i in range(nq): 
    dwide[i], dclos[i]= func1(q[i])
#######################################################################33

plt.clf()
fig= plt.figure(figsize=(8,6))
plt.plot(np.log10(dwide), logq, "r--", lw=2.0, label= r"$\rm{d}_{\rm w}(q)$")
plt.plot(np.log10(dclos), logq, "b--", lw=2.0, label= r"$\rm{d}_{\rm c}(q)$")


for i in range(int(ns)):
    if(i!=1 and i!=2):
        plt.errorbar(res[i,1],res[i,0], xerr=res[i,8],fmt='o',color='darkgreen',ecolor='yellowgreen',capsize=0.0, ms=7.0, lw=2.2)
        if(i==8):     
            plt.text(res[i,1]-0.7*res[i,8],res[i,0]-0.15,lab[i],fontsize=12)
        elif(i==7):     
            plt.text(res[i,1]+0.7*res[i,8],res[i,0],     lab[i],fontsize=12)    
        elif(i==9 or i==10):     
            plt.text(res[i,1]-0.9*res[i,8],res[i,0]-0.1, lab[i],fontsize=12)         
        else: 
            plt.text(res[i,1]-0.7*res[i,8],res[i,0]+0.07,lab[i],fontsize=12)
plt.xlabel(r"$\log_{10}[\rm{d}(R_{\rm E})]$", fontsize=18)
plt.ylabel(r"$\log_{10}[q]$",fontsize=18)
plt.ylim(-5.0,-0.8)
#plt.xscale('log')
plt.xlim(-1.7,0.24)
fig=plt.gcf()
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True, fontsize=18)
plt.grid("True")
plt.grid(linestyle='dashed')
fig.tight_layout()
fig.savefig("fig1.jpg",dpi=200)
py.clf()


