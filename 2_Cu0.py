import numpy as np 
import math
from numpy import conj
import matplotlib 
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as py 
from matplotlib import rcParams
from matplotlib import gridspec
import matplotlib.gridspec as gridspec
import VBBinaryLensingLibrary
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
matplotlib.rcParams['text.usetex']=True
cmap=plt.get_cmap('viridis')
############################################################
parcs= 30856775814913673.0
G=6.67384*pow(10.,-11.0)#;// in [m^3/s^2*kg].
velocity=3.0*pow(10.0,8.0)#;//velosity of light
Msun=1.98892*pow(10.,30)#; //in [kg].
AU=  1.495978707*(10.0**11.0)
const= float(np.sqrt(1000.0*parcs*G*Msun)/velocity/AU)
Ds=8.0#*1000.0*parsec
Dl=4.0#*1000.0*parsec
Dls=Ds-Dl

nw=int(123)
sig=np.zeros((123,2))
sig=np.loadtxt("./files/sigma_WFIRST.txt")
nam=[r"$u_{0}$", r"$,~~~\Delta \rm{FWHM}(t_{\rm E})$",  r"$,~~~\Delta \chi^{2}$"]

vbb=VBBinaryLensingLibrary.VBBinaryLensing()
vbb.Tol=1.e-5;
vbb.a1=0.0;
vbb.LoadESPLTable("./files/ESPL.tbl");
############################################################ 
### Earth-moon systems 
cade=float(15.16/60.0/24.0)##days
epsi=float(0.006)
q=0.03
dis=0.7
tE=0.1
ros=0.4
u0=0.0
teta=0.0
dt= float(tE/120.0)
ntim=int(5.0*tE/dt)+1
angle = np.linspace(0, 2*np.pi , 150 ) 
xsou = ros * np.cos( angle ) 
ysou = ros * np.sin( angle )
tim=np.zeros((ntim))
Astar=np.zeros((ntim,2))
traj= np.zeros((ntim,2))
data= np.zeros((ntim,4)) 
############################################################ 
xmin=-2.1
xmax= 2.1
ymin=-2.1
ymax= 2.1
nx=6000
ny=6000 
dx=float(xmax-xmin)/nx
dy=float(ymax-ymin)/ny
xl1= -dis*q/(1.0+q)
xl2=    dis/(1.0+q)
yl1=0.0 
yl2=0.0
mu1=1.0/(1.0+q)
mu2=  q/(1.0+q)
xsmin= -0.8
xsmax=  0.2
ysmin= -0.5
ysmax=  0.5
dsx=(xsmax-xsmin)/nx
dsy=(ysmax-ysmin)/ny
############################################################ 
def LensEq(xl, yl, xl1, yl1, xl2, yl2, mu1, mu2):
    x1=float(xl-xl1)
    y1=float(yl-yl1)
    x2=float(xl-xl2)
    y2=float(yl-yl2)
    d1=float(x1*x1+ y1*y1+1.0e-52)
    d2=float(x2*x2+ y2*y2+1.0e-52)
    xs=float(xl-mu1*x1/d1 - mu2*x2/d2)
    ys=float(yl-mu1*y1/d1 - mu2*y2/d2)
    z=complex(xl,yl);  z1=complex(xl1,yl1);    z2=complex(xl2,yl2)
    zb=conj(z);        z1b=conj(z1);           z2b=conj(z2);
    f= mu1/((zb-z1b)**2.0) + mu2/((zb-z2b)**2.0)
    ep= abs(np.sqrt((f.real)**2.0 + (f.imag)**2.0)-1.0)
    if(ep<epsi):  
        flag=1
    else:  
        flag=0    
    return(xs, ys, flag)
############################################################
def fwhm(mag,tim):  
    ffl=-1;    
    FWHM=0.0
    magm= abs(np.max(mag)-np.min(mag) )
    for i in range(len(mag)):  
        if(abs(mag[i]-np.max(mag))>abs(0.5*magm) and ffl<0): 
            FWHM=tim[i]
            ffl=0
        if(abs(mag[i]-np.max(mag))<abs(0.5*magm) and ffl==0): 
            FWHM=abs(FWHM-tim[i])
            ffl=1  
    return(FWHM)
############################################################    
def Eroman(mag):
    w=-1;
    if(mag<sig[0,0] or mag==sig[0,0]):   w=0;
    elif(mag>sig[nw-1,0] or mag==sig[nw-1,0]): w=nw-1;
    else: 
        for i in range(nw-1): 
            if(float((mag-sig[i,0])*(mag-sig[i+1,0]))<0.0 or mag==sig[i,0]):
                w=i; 
                break
    if(w==-1):
        print("Error maq_wfirst: ", mag)
        input("Enter a number ") 
    return(sig[w,1])     
############################################################
cau1=[]
cau2=[]
crt1=[]
crt2=[]
cau=np.zeros((ny, nx))+0.01
for i1 in range(nx):
    xl=float(xmin+i1*dx)
    for j1 in range(ny):
        yl=float(ymin+j1*dy)
        xs,ys, flag=LensEq(xl , yl, xl1, yl1, xl2, yl2, mu1, mu2)   
        rx=int((xs-xsmin)/dsx)           
        ry=int((ys-ysmin)/dsy)
        if(rx>=0 and rx<nx and ry>0 and ry<ny):
            cau[ny-ry, rx] +=1.0 
        if(flag>0 and xs>xsmin and xs<xsmax and ys>ysmin and ys<ysmax):  
            cau1.append(xs);  cau2.append(ys)
            crt1.append(xl);  crt2.append(yl) 
cau=np.log10(cau)
vbb.PrintCau(dis, q, 0)
f30=open("./files/MONT/cau_{0:d}.dat".format(0),"r")
nc=  sum(1 for line in f30) 
caus= np.zeros((nc,2))
caus= np.loadtxt("./files/MONT/cau_{0:d}.dat".format(0) )
#############################################################
plt.clf()
plt.figure(figsize=(12,6))
for j in range(24): 
    u0=0.0+j*0.024##float(j*7.31456456741345*np.pi/180.0)
    chi1=0.0; chi2=0.0;  time=0.0; k=0;     
    for i in range (ntim):
        tim[i]=float(-2.5*tE+ i*dt)/tE
        traj[i,0]= tim[i]*np.cos(teta) -u0*np.sin(teta)
        traj[i,1]= tim[i]*np.sin(teta) +u0*np.cos(teta)  
        xcent2= traj[i,0]+abs(dis*q/(1.0+q))
        uu=np.sqrt(xcent2**2.0 + traj[i,1]**2.0)
        Astar[i,0]= -2.5*np.log10(vbb.ESPLMag2(uu,ros))+21.0 
        Astar[i,1]= -2.5*np.log10(vbb.BinaryMag2(dis , q, traj[i,0], traj[i,1] , ros))+21.0
        sigma=Eroman(Astar[i,1])
        flac= np.random.normal(0.0,sigma) 
        time+=dt
        if(time>cade):
            time-=cade;  
            data[k,:]=np.array([tim[i],  Astar[i,1]+flac, sigma,  Astar[i,1]+flac-Astar[i,0] ])  
            k+=1
            chi1+= float(Astar[i,1]+flac-Astar[i,0])**2.0/sigma/sigma
            chi2+= flac**2.0/sigma/sigma
##############################################################   
    dchi=abs(chi1-chi2)
    FWHM0= fwhm(Astar[:,0], tim)
    FWHM1= fwhm(Astar[:,1], tim)
    ymin= np.min(Astar);   ymax= np.max(Astar);   dy0= float(ymax-ymin)/6.0
    ymin2=np.min(Astar[:,1]-Astar[:,0]);  ymax2=max(Astar[:,1]-Astar[:,0]);  dy2= float(ymax2-ymin2)/4.0
    
    
    plt.clf()
    plt.cla()
    plt.gcf()
    spec3 = gridspec.GridSpec(nrows=3, ncols=2)
    ax1=plt.subplot(spec3[:2, 0])
    ####################################################################################
    ax1.errorbar(data[:k,0], data[:k,1], yerr=data[:k,2], fmt='mo',capsize=0, ms=1.5, lw=1.0) 
    ax1.plot(tim, Astar[:,0], "k--", lw=1.4, label=r"\rm{ESPL}")
    ax1.plot(tim, Astar[:,1], "b-",  lw=1.4, label=r"$\rm{Binary}~\rm{Lens}$")
   
    ax1.set_ylim(ymax+dy0/4.0,ymin-dy0/4.0)
    ax1.set_xlim( np.min(tim) , np.max(tim) )
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    
    ax1.set_ylabel(r"$\rm{W149}-band$",fontsize=17) 
    ax1.set_title(str(nam[0])+ '={0:.1f}'.format(u0) + str(nam[1])+ '={0:.2f}'.format(float(FWHM1-FWHM0)) +str(nam[2])+ '={0:.1f}'.format(dchi),fontsize=15,color='k')
    plt.grid(True)
    plt.grid(linestyle='dashed')
    plt.legend()
    plt.legend(loc='best',fancybox=True, shadow=True, fontsize=15)
    ###################################################################
    
    ax2=plt.subplot(spec3[2, 0], sharex=ax1)
    ax2.errorbar(data[:k,0], data[:k,3], yerr=data[:k,2], fmt='mo',capsize=0, ms=1.5, lw=1.0) 
    ax2.plot(tim, Astar[:,1]-Astar[:,0] ,"k-", lw=2.0)
    
    ax2.set_ylim(ymax2+dy2*0.2 , ymin2-dy2*0.2)
    ax2.set_xlim(xmin,xmax)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    ax2.set_xlabel(r"$time(t_{\rm E})$",fontsize=17, labelpad=0.1)
    ax2.set_ylabel(r"$\rm{\Delta m}(mag)$",  fontsize=17,labelpad=0.1)
    plt.grid(True)
    plt.grid(linestyle='dashed')
    plt.setp(ax1.get_xticklabels(), visible=False)
    yticks = ax1.yaxis.get_major_ticks()
    
    ####################################################################################       
    ax3=plt.subplot(spec3[:, 1])
    ax3.imshow(cau[:,:],extent=(xsmin,xsmax,ysmin,ysmax),cmap=cmap)
    ax3.scatter(caus[800:1600,0], caus[800:1600,1], color="y", s=0.4)
    ax3.plot(traj[:,0], traj[:,1], "k--", lw=2.0)
    ax3.plot(xsou, ysou, "k-", lw=1.2)

    ax3.set_ylim(-0.5,0.5)
    ax3.set_xlim(-0.8,0.2)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    ax3.set_xlabel(r"$\rm{x}(R_{E})$",fontsize=17, labelpad=0.0)
    ax3.set_ylabel(r"$\rm{y}(R_{E})$",fontsize=17, labelpad=0.0)
    mpl.rcParams['axes.formatter.use_mathtext'] = 'True'
    mpl.rcParams['axes.formatter.useoffset'] = 'False'
    mpl.rcParams['figure.subplot.left'] = 1.0
    mpl.rcParams['figure.subplot.right'] = .95
    mpl.rcParams['figure.subplot.bottom'] = .20
    mpl.rcParams['figure.subplot.top'] = .90
    fig=plt.gcf()
    fig=plt.savefig("./lights/central_u0/LIc_{0:d}.jpg".format(j),dpi=200)
    plt.clf()
    plt.cla()
    py.clf()
    print("First light curve was plotted ***** ")
    print("dchi:  ",    dchi)        
    print("FWHM:    ",  FWHM0, FWHM1)
    print("*****************************")
############################################################################
    






