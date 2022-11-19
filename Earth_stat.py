import math
import numpy as np 
from numpy import conj
import matplotlib
import matplotlib.pyplot as plt
import pylab as py 
from matplotlib import rcParams
from matplotlib import gridspec
import matplotlib.gridspec as gridspec
import VBBinaryLensingLibrary
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import FixedLocator, FixedFormatter
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
matplotlib.rcParams['text.usetex']=True
cmap=plt.get_cmap('viridis')
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import matplotlib.font_manager
matplotlib.get_cachedir()
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ( AutoLocator, AutoMinorLocator)
import pandas as pd
import math


epsi=0.00000001;
ns=12
nw=int(123)
sig=np.zeros((123,2))
sig=np.loadtxt("./files/sigma_WFIRST.txt")
######################################################################
nam1=[r"$counter$", r"$line$",  r"$lat$",  r"$long$",  r"$strucl$",  r"$Ml$" , r"$Rl$",  r"$Dl$",  r"$vl$",  r"$mul1$", r"$mul2$", r"$strucs$", r"$CL$", r"$Ms$", r"$Ds$",  r"$T*$", r"$R*$",  r"$logl$", r"$type*$", r"$col$", r"$vs$",  r"$mus1$", r"$mus2$",  r"$q$", r"$dis$",  r"$period$",  r"$semi(Rl)$",  r"$Config$", r"$tE$", r"$RE$", r"$mul$", r"$Vt$", r"$u0$", r"$ksi$", r"$piE$", r"$tetE$", r"$opd$", r"$ros$", r"$MabI$", r"$MabW$", r"$mapI$", r"$mapW$", r"$magbI$", r"$magbW$", r"$fbI$", r"$fbW$", r"$NbI$", r"$NbW$", r"$ExI$", r"$ExW$", r"$\Delta \chi^{2}$", r"$Asym$"]    

nam2=[r"$\log_{10}[t_{\rm E}(\rm{days})]$", r"$u_{0}$", r"$\xi(degree)$",  r"$\log_{10}[q]$", r"$\log_{10}[d(R_{\rm E})]$", r"$\log_{10}[\rho_{\star}]$", r"$m_{\rm{W149}}(mag)$", r"$f_{\rm{W149}}$", r"$s(R_{\rm h})$", r"$\log_{10}[\Delta \chi^{2}]$" , r"$\log_{10}[\Delta mathcal{A}]$", r"$Configuration$"] 
 

f1=open("./files/MONT/MonteEarth.txt","r")
nm= sum(1 for line in f1) 
arr= np.zeros((nm,52))
arr=np.loadtxt("./files/MONT/MonteEarth.txt")
fil1=open('./EarthPap.txt',"w"); 
fil2=open('./files/MONT/Earth.txt',"w"); 
fil1.close();  fil2.close()


vbb=VBBinaryLensingLibrary.VBBinaryLensing()
vbb.Tol=1.e-5;
vbb.a1=0.0;
vbb.LoadESPLTable("./files/ESPL.tbl");
######################################################################  
def histo(array, i,  flag): 
    plt.clf()
    plt.cla()
    fig, ax1= plt.subplots(figsize=(8 , 6))
    plt.hist(array,30, density=True, histtype='bar',ec='black',facecolor='green', alpha=0.8, rwidth=0.95)
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    plt.ylabel(r"Histogram",fontsize=18,labelpad=0.1)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig3= plt.gcf()
    if(flag==0):  
        plt.xlabel(str(nam1[i]),fontsize=18,labelpad=0.1)
        fig3.savefig("./lights/Histo/Earth/Total_{0:d}.jpg".format(i) ,dpi=200)
    if(flag==1):  
        plt.xlabel(str(nam2[i]),fontsize=18,labelpad=0.1)
        fig3.savefig("./lights/Histo/Earth/Detect_{0:d}.jpg".format(i) , dpi=200)
    return(0)    
########################################################################      
def histo2(array1, array0, i):
    plt.clf()
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    ax1= plt.gca()
    plt.hist(array0, 27, histtype='bar',ec='darkgreen',facecolor='green',alpha=0.45, rwidth=1.5)
    plt.hist(array1,27, histtype='step',color='k', alpha=0.9, lw=2.5)
    ax1.set_ylabel(r"$\rm{Normalized}~\rm{distribution}$",fontsize=19,labelpad=0.1)
    ax1.set_xlabel(str(nam2[i]),fontsize=19, labelpad=0.1)
    y_vals = ax1.get_yticks()
    ax1.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/len(array0))) for x in y_vals]) 
    y_vals = ax1.get_yticks()
    plt.ylim([np.min(y_vals), np.max(y_vals)])
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    plt.axvline(x=np.mean(array0) , color='darkgreen', linestyle='--', lw=1.5)
    plt.axvline(x=np.mean(array1) , color='k',         linestyle='--', lw=1.5)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig3= plt.gcf()
    fig3.savefig("./lights/Histo/Earth/Both_{0:d}.jpg".format(i) , dpi=200)
    return(0)        
########################################################################  
'''
def Eroman(mag):
    w=-1;
    if(mag<sig[0,0] or mag==sig[0,0]):   w=0;
    elif(mag>sig[nw-1,0] or mag==sig[nw-1,0]): w=nw-1;
    else: 
        for ia in range(nw-1): 
            if(((mag-sig[ia,0])*(mag-sig[ia+1,0]))<0.0 or mag==sig[ia,0]):
                w=ia; 
                break
    if(w==-1):
        print("Error maq_wfirst: ", mag)
        input("Enter a number ") 
    return(sig[w,1])     
#######################################################################
dt= float(1.16/60.0/24.0)
def light(tE, ksi, u0, dis, q, ros, fbW, magbW):   
    if(ros>1.0): 
        ntim=int(5.0*tE*ros/dt)+1
        pt1=float(-2.5*tE*ros)
    else: 
        ntim=int(5.0*tE/dt)+1
        pt1=float(-2.5*tE)
    chi1=0.0; chi2=0.0; chi0=0.0;  imax=0;  maxm=0.0
    dat=np.zeros((ntim,4))
    for j in range(ntim):
        tim=float(pt1+j*dt)  
        xcent= float(tim/tE *np.cos(ksi) -u0*np.sin(ksi) )
        ycent= float(tim/tE *np.sin(ksi) +u0*np.cos(ksi) )
        xcent2=float(xcent + abs(dis*q/(1.0+q)) )
        
        uu=np.sqrt(xcent2*xcent2 + ycent*ycent);
        Astar0= vbb.ESPLMag2(uu,ros) 
        a0=-2.5*np.log10(fbW*Astar0+1.0-fbW) + magbW
        
        uu=np.sqrt(xcent*xcent + ycent*ycent);
        Bstar0= vbb.ESPLMag2(uu,ros) 
        b0=-2.5*np.log10(fbW*Bstar0+1.0-fbW) + magbW

        Astar1= fbW*vbb.BinaryMag2(dis,q,xcent,ycent,ros) + 1.0-fbW
        a1=-2.5*np.log10(Astar1) + magbW
        sigma=Eroman(a1)
        flac=np.random.normal(0.0,sigma)
         
        chi0+= (a1+flac-b0)**2.0/sigma/sigma
        chi1+= (a1+flac-a0)**2.0/sigma/sigma
        chi2+= (a1+flac-a1)**2.0/sigma/sigma
        
        sigA=float(Astar1*abs(pow(10.0,-0.4*sigma)-1.0))
        flac=np.random.normal(0.0,sigA) 
        dat[j, :]=np.array([ tim , Astar1+flac , sigA, fbW*Astar0+1.0-fbW ])
        print( dat[j,:] )
        if(Astar1>maxm): 
            maxm=Astar1; 
            imax=int(j)
    plt.clf()
    plt.plot(dat[:,0], dat[:,1], "r--")        
    plt.plot(dat[:,0], dat[:,3], "b--") 
    plt.show()
    asym=0.0
    for k in range(imax): 
        if(int(k+imax+1)<ntim and int(imax-k-1)>0): 
            err=   0.5*(dat[k+imax+1,2] +dat[imax-k-1,2] )
            asym+=float(dat[k+imax+1,1] -dat[imax-k-1,1] )**2.0/err/err
    dchi1=abs(chi1-chi2)
    dchi2=abs(chi0-chi2)
    
    print (dchi1, dchi2, asym)
    input("Enter a number")              
    if(dchi1<dchi2): 
        return(dchi1 , asym)
    else: 
        return(dchi2 , asym) 
'''
####################################################################### 
for i in range(52): 
    if(i==23 or i==24 or i==28 or i==37): 
        histo(np.log10(arr[:,i]), i, 0)
    else:  
        histo(arr[:,i], i, 0)
    print("histogram is plotted ", i)    
det= np.zeros((nm, ns)) 
######################################################################
k1=0;  ntot=int(0)
for i in range(nm): 
    icon, lin, lat, lon, strucl=       arr[i,0], arr[i,1], arr[i,2], arr[i,3], arr[i,4]
    Ml, Rl, Dl, vl,  mul1, mul2=       arr[i,5], arr[i,6], arr[i,7], arr[i,8], arr[i,9], arr[i,10]
    strucs, cl, Ms, Ds, Tstar, Rstar=  arr[i,11],arr[i,12],arr[i,13],arr[i,14],arr[i,15],arr[i,16]
    logl, types, col, vs, mus1, mus2=  arr[i,17],arr[i,18],arr[i,19],arr[i,20],arr[i,21],arr[i,22]
    q, dis, peri, semi, config, tE=    arr[i,23],arr[i,24],arr[i,25],arr[i,26],arr[i,27],arr[i,28]
    RE, mul, Vt, u0, ksi, piE=         arr[i,29],arr[i,30],arr[i,31],arr[i,32],arr[i,33],arr[i,34]
    tetE,opd,ros,MabI, MabW, mapI=     arr[i,35],arr[i,36],arr[i,37],arr[i,38],arr[i,39],arr[i,40]
    mapW, magbI, magbW, fbI, fbW, NbI= arr[i,41],arr[i,42],arr[i,43],arr[i,44],arr[i,45],arr[i,46]
    NbW, ExI, ExW, dchi, asym=         arr[i,47],arr[i,48],arr[i,49],arr[i,50],arr[i,51]
    if(ros>0.0):
        ntot+=1
        test=np.array([tE, u0, ksi, q, dis, ros, magbW, fbW, semi, dchi, asym, config])
        fil2=open('./files/MONT/Earth.txt',"a")
        np.savetxt(fil2 , test.reshape(-1,12), fmt="%.5f  %.5f   %.3f   %.7f    %.5f   %.6f    %.5f   %.5f   %.4f    %.3f    %.3f    %.1f")
        fil2.close() 
        if(dchi>800.0 and asym>1200.0):
            det[k1,:]= test
            k1+=1    
            
if(1>0): 
    print ("tot number:    ",  nm,    ntot,    k1 )
    print ("Efficiency:  ",  k1*100.0/ntot/1.0)        
    print("**************************************")           
    ave=np.zeros((ns,2))
    clos=0.0;  wide=0.0;  inte=0.0
    for j in range(k1): 
        if(det[j,11]==1):   clos+=1
        elif(det[j,11]==3): wide+=1
        else:               inte+=1
    for j in range(ns):
        ave[j,0]=np.mean(det[:k1,j])
        ave[j,1]= np.std(det[:k1,j])
        if(j==0 or j==3 or j==5 or j==9  or j==10): 
            histo(np.log10(det[:k1,j]), j, 1)
        else:  
            histo(det[:k1,j], j, 1)         
    histo2( np.log10(det[:k1,0]) ,  np.log10(arr[:ntot,28]) , 0 )##tE
    histo2( np.log10(det[:k1,3]) ,  np.log10(arr[:ntot,23]) , 3 )##q
    histo2( np.log10(det[:k1,5]) ,  np.log10(arr[:ntot,37]) , 5 )##ros
    histo2( np.log10(det[:k1,4]) ,  np.log10(arr[:ntot,24]) , 4 )##dis     
    histo2(det[:k1,2]*180.0/np.pi ,  arr[:ntot,33]*180.0/np.pi,2)##ksi    
    histo2(det[:k1,8] ,  arr[:ntot,26] , 8 )##semi     
    histo2(det[:k1,1] ,  arr[:ntot,32] , 1 )##u0       
    clos=float(clos*100.0/(k1+epsi)/1.0)
    wide=float(wide*100.0/(k1+epsi)/1.0)
    inte=float(inte*100.0/(k1+epsi)/1.0)
    effi=float(k1*100.0/ntot/1.0)
########################################################################
    sav=np.array([1.0, ave[0,0],ave[0,1],  ave[1,0],ave[1,1], np.log10(ave[3,0]),np.log10(ave[3,1]),  
                           np.log10(ave[4,0]),np.log10(ave[4,1]),  np.log10(ave[5,0]),np.log10(ave[5,1]),  
                           ave[6,0],ave[6,1], ave[8,0],ave[8,1], clos , wide, inte, effi ])##19     
    fil1=open('./EarthPap.txt',"a")
    np.savetxt(fil1, sav.reshape(-1,19), fmt="$%.1f$ & $%.2f \pm %.2f$ &  $%.2f \pm %.2f$  &  $%.2f \pm %.2f$ &  $%.2f \pm %.2f$ &  $%.2f\pm %.2f$ &  $%.1f \pm %.1f$ & $%.1f \pm %.1f$  & $%.1f$:$%.1lf$:$%.1f$ &  $%.3f$")  
    fil1.close()   
## tE,   u0,   log q,  log dis,   log ros,   magbw,  semi, close:wide:intermediate  
########################################################################
   
   
   
   
   
   
   
   
