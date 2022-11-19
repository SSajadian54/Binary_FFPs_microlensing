#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
#include "VBBinaryLensingLibrary.h"
using namespace std;

const int Num=2000;
const double MaxD=20.0;///kpc
const double RA=180.0/M_PI;
const double pi= M_PI;
const double step=MaxD/(double)Num/1.0;///step in kpc
const double Hp=6.62607004*pow(10.0,-34);
const double KP=3.08568025*pow(10.,19); // in meter.
const double G=6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity= 3.0*pow(10.0,8.0);//velosity of light
const double Msun=1.98892*pow(10.,30); //in [kg].
const double Rsun= 6.957*pow(10.0,8.0); ///solar radius [meter]
const double Mjupiter=1.898*pow(10.0,27.0); 
const double Rjupiter= 71.492*pow(10.,6);//radius of jupiter in [m]
const double Mearth=  5.9722*pow(10.0,24.0); 
const double Rearth= 6.3781*pow(10.0,6.0);//meter
const double vro_sun=226.0;
const double AU=1.4960*pow(10.0,11.0);
const double year=365.2421875;
const double binary_fraction=double(2.0/3.0);
const double Avks=double(8.20922);
const double VSunR =11.1;
const double VSunT =vro_sun*(1.00762+0.00712)+ 12.24;
const double VSunZ =7.25;
///============================ Besancon constant ==============///
const double R_sun=8.5;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
///============================ WFIRST & OGLE &  KMTNet ===================
const int M=5;///number of filters  V,I,K,H,W149
const double satu[M]={12.0, 12.0, 13.0, 13.0, 14.8}; //it should be changed
const double thre[M]={19.5, 20.5, 21.0, 21.0, 26.0};
const double FWHM[M]={0.22, 0.22, 0.22, 0.22 , 0.22};//3*pixel_size (0.11") of WFIRST, VIKH W149
const double AlAv[M]={1.009,0.600,0.118,0.184,0.225};///From besancon model[VIKH W149]
const double sigma[M]={0.022,0.022,0.02,0.025,0.025};//MOAاستفاده از مقاله کاردلی
////=======================================================
const int Nw=int(123);///sigma_WFIRST.dat
const int YZ=3615;///number of rows in file yzma
const int N1=36224, N2=25000, N3=3818, N4=3500;///CMD_BESANCON, thinD, bulge, thickD, halo
///======================================================


struct source{
    int nums,struc, cl;
    double Ds,TET,FI, lat, lon;
    double od_disk,od_ThD,od_bulge,od_halo,opd;
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs,nstart;
    double Nblend[M], blend[M], Fluxb[M], magb[M], Ai[M], Mab[M], Map[M];
    double type, Tstar, logl, col, Rstar, mass;
    double Romax,ros, deltao;
    double SV_n1, LV_n1, VSun_n1;
    double SV_n2, LV_n2, VSun_n2;
    double mus1, mus2, mul1, mul2, xi;
};
struct lens{
    int numl,struc, con;
    double Ml, Mp, Rl, semi,Dl, vl , vs, Vt, xls;
    double rhomaxl,tE,RE;
    double mul,u0;
    double tetE;
    double piE; 
    double q, dis, ksi; 
    double  pt1, pt2; 
    double peri;
        double magw[Nw], errw[Nw];
};
struct CMD{
    double Teff_d[N1],logl_d[N1],Mab_d[M][N1],Rs_d[N1],mass_d[N1],type_d[N1]; int cl_d[N1];  ///thin disk
    double Teff_b[N2],logl_b[N2],Mab_b[M][N2],Rs_b[N2],mass_b[N2],type_b[N2]; int cl_b[N2];  /// bulge
    double Teff_t[N3],logl_t[N3],Mab_t[M][N3],Rs_t[N3],mass_t[N3],type_t[N3]; int cl_t[N3];  ///thick disk
    double Teff_h[N4],logl_h[N4],Mab_h[M][N4],Rs_h[N4],mass_h[N4],type_h[N4]; int cl_h[N4];  /// halo
};
struct extinc{
   double dis[100];///distance
   double Extks[100];///ks-band extinction
   double Aks;
};
///===================== FUNCTION =================================
int Extinction(extinc & ex,source & s);
void read_cmd(CMD & cm);
void func_source(source & s, CMD & cm , extinc & ex);
void func_lens(lens & l, source & s, CMD &  cm, extinc & ex);
void vrel(source & s , lens & l);
void Disk_model(source & s);
void optical_depth(source & s);
double Interpol(double ds, extinc & ex);
double RandN(double sigma, double Nn);
double RandR(double down, double up);
double ErrorRoman(lens & l, double mag);
//////////////////////////////////////////////////////
    time_t _timeNow;
      unsigned int _randVal;
    unsigned int _dummyVal;
    FILE * _randStream;
///////////////////////////////////////////////////////
///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//
int main()
{
///****************************************************************************
//gettimeofday(&newTime,NULL);
//ftime(&tp);
	time(&_timeNow);
	_randStream = fopen("/dev/urandom", "r");
	_dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
	srand(_randVal);
    _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
	srand(_randVal);
///******************************************************************************     
    source s;
    lens l;
    CMD cm;
    extinc ex;
    read_cmd(cm);
    
     VBBinaryLensing vbb;
     vbb.Tol=1.e-5; 
     vbb.a1=0.0; 
     vbb.LoadESPLTable("./files/ESPL.tbl");
    

    char filnam[40]; 
    FILE* fil1;  FILE* data; 
    fil1=fopen("./files/MONT/MonteJupiter.txt","a+");
    fclose(fil1);
    
///===========================================================
    FILE * wfirstf;
    wfirstf=fopen("./files/sigma_WFIRST.txt","r");
    for(int i=0; i<Nw;++i){
    fscanf(wfirstf,"%lf  %lf\n",&l.magw[i],&l.errw[i]);}
    fclose(wfirstf);
    cout<<"*****  sigma_wfirst was read *****  "<<endl;
///===========================================================
    
    
    int counter=0;
    double lonn0, lonn, lat0, test;
    double chi0, chi1, chi2,  maxm, tmax;  
    double tim, xcent, ycent, xcent2, uu,  Astar0, Astar1, Bstar0; 
    double a0, b0, a1, sigma, flac, asym0, asym1;    
    double err, dchi1, dchi2, nnm; 
    double dt=double(15.16/60.0/24.0);
    int  imax, nd; 
///*********************************************************************************             
   for(int icon=1; icon<20000; ++icon){
    for (int li=1; li<=7; li++){
    if(li==1) {lonn0=1.3;  lat0=-0.875; }
    if(li==2) {lonn0=0.9;  lat0=-0.875; }
    if(li==3) {lonn0=1.3;  lat0=-1.625; }
    if(li==4) {lonn0=0.9;  lat0=-1.675; }
    if(li==5) {lonn0=0.5;  lat0=-1.675; }
    if(li==6) {lonn0=0.1;  lat0=-1.675; }
    if(li==7) {lonn0=-0.3; lat0=-1.675; }
    
    
    s.lat= lat0+ double(RandR(-0.375,0.375));
    lonn =lonn0+ double(RandR(-0.200,0.200));
    if(lonn<=0.0)   s.lon=360.0+lonn;
    else            s.lon=lonn;
    s.TET=(360.0-s.lon)/RA;///radian
    s.FI=s.lat/RA;///radian
    Disk_model(s);
    if(int(Extinction(ex,s))==1){


    do{
    func_source(s, cm, ex);
    func_lens(l, s,cm, ex);
    test=RandR(0.0 ,  1.0);
    }while(s.magb[4]>thre[4] or  test>s.blend[4] or l.tE<0.0001 or l.tE>1000.0 );
    optical_depth(s);
    counter+=1;

    
///######################################################################################
    //sprintf(filnam,"./files/MONT/%c%c%c%c%d.dat",'d','a','t','_', counter);
    //data=fopen(filnam,"w");
    chi0=0.0; chi1=0.0; chi2=0.0; maxm=0.0; tmax=0.0; imax=0;     asym0=0.0; asym1=0.0;
    l.pt1=double(-4.5*l.tE);  
    l.pt2=double(+4.5*l.tE);  
    nd= int((l.pt2-l.pt1)/dt)+1 ; 
    double astar[nd]={0.0}; 
    
    for(int i=0; i<nd;  ++i){
    tim=double(l.pt1+i*dt);  
    xcent= double(tim/l.tE * cos(l.ksi) -l.u0*sin(l.ksi) );
    ycent= double(tim/l.tE * sin(l.ksi) +l.u0*cos(l.ksi) );
    xcent2=double(xcent + fabs(l.dis*l.q/(1.0+l.q)) ); 
    uu=sqrt(xcent2*xcent2 + ycent*ycent);
    Astar0= vbb.ESPLMag2(uu, s.ros) ;
    uu=sqrt(xcent*xcent + ycent*ycent);
    Bstar0= vbb.ESPLMag2(uu, s.ros); 
    a0=-2.5*log10(s.blend[4]*Astar0+1.0-s.blend[4]) + s.magb[4];
    b0=-2.5*log10(s.blend[4]*Bstar0+1.0-s.blend[4]) + s.magb[4];
    Astar1= s.blend[4]*vbb.BinaryMag2(l.dis,l.q,xcent,ycent, s.ros) + 1.0-s.blend[4];
    a1=-2.5*log10(Astar1) + s.magb[4]; 
    sigma= ErrorRoman(l, a1);
    flac=RandN(sigma, 2.0);
    astar[i]=a1;
    chi0+= double(a1+flac-b0)*(a1+flac-b0)/sigma/sigma;
    chi1+= double(a1+flac-a0)*(a1+flac-a0)/sigma/sigma;
    chi2+= double(a1+flac-a1)*(a1+flac-a1)/sigma/sigma;
    if(Astar1>maxm){maxm=Astar1;  tmax=tim;    imax=i;}
    //fprintf(data,"%.4lf    %.5lf   %.5lf    %.5lf    %.5lf     %.5lf\n",tim, a1, a0, b0, xcent, ycent);
    }
    //fclose(data);
    
    dchi1=fabs(chi1-chi2);
    dchi2=fabs(chi0-chi2);
    if(dchi1>dchi2) dchi1=dchi2; 
    nnm=0.0;
    for(int i=0; i<imax; ++i){  
    if(int(i+imax+1)<nd and int(imax-i-1)>=0){
    tim=double(tmax+i*dt);  
    xcent= double(tim/l.tE * cos(l.ksi) -l.u0*sin(l.ksi) );
    ycent= double(tim/l.tE * sin(l.ksi) +l.u0*cos(l.ksi) );
    Astar1=-2.5*log10(s.blend[4]*vbb.BinaryMag2(l.dis,l.q,xcent,ycent, s.ros) +1.0-s.blend[4] ) +s.magb[4];
    tim=double(tmax-i*dt);  
    xcent= double(tim/l.tE * cos(l.ksi) -l.u0*sin(l.ksi) );
    ycent= double(tim/l.tE * sin(l.ksi) +l.u0*cos(l.ksi) );
    Astar0=-2.5*log10(s.blend[4]*vbb.BinaryMag2(l.dis,l.q,xcent,ycent, s.ros) +1.0-s.blend[4] ) +s.magb[4];
    err=ErrorRoman(l , Astar1); 
    flac=RandN(err , 2.5); 
    asym1+=double(Astar1+flac-Astar0)*(Astar1+flac- Astar0)/err/err;
    nnm+=1.0; }}
    //asym0+=double(Astar1+flac-Astar1)*(Astar1+flac- Astar1)/err/err;}}
    asym1=fabs(asym1-nnm); 
///######################################################################################   
        
    fil1=fopen("./files/MONT/MonteJupiter.txt","a+");
    fprintf(fil1,
    "%d  %d  %.5lf  %.5lf  "//4
    "%d  %.5lf   %.5lf   %.5lf  %.5lf  %.5lf  %.5lf  "//11
    "%d  %d  %.5lf   %.5lf  %.4lf  %.4lf  %.4lf  %.2lf  %.3lf  %.4lf  %.4lf   %.4lf  "//23
    "%.6lf  %.7lf  %.5lf   %.5lf  %d "  //28
    "%.8lf  %.7lf   %.5lf    %.5lf  %.5lf  %.5lf  %.7lf  %.7lf  %.6lf  %.6lf  "//38
    "%.5lf   %.5lf  %.5lf   %.5lf   %.5lf  %.5lf  %.5lf  %.5lf  %.2lf   %.2lf   %.4lf   %.4lf  %.5lf  %.5lf\n", //52
    counter, li, s.lat, s.lon, //49
    l.struc, l.Ml/Mearth, l.Rl/Rearth, l.Dl, l.vl, s.mul1, s.mul2, //11
    s.struc, s.cl, s.mass, s.Ds, s.Tstar, s.Rstar, s.logl, s.type, s.col, l.vs,s.mus1, s.mus2,//23
    l.q, l.dis, l.peri, l.semi/l.Rl, l.con, //28
    l.tE, l.RE/AU, l.mul*30.0, l.Vt, l.u0, l.ksi, l.piE, l.tetE,s.opd*1.0e6, s.ros,///38
    s.Mab[1], s.Mab[4], s.Map[1], s.Map[4], s.magb[1], s.magb[4],s.blend[1], s.blend[4],s.Nblend[1],s.Nblend[4],s.Ai[1],s.Ai[4],  
    dchi1, asym1);//52
    fclose(fil1);
   
    if(int(counter)%1==0){
    cout<<"dchi:  "<<dchi1<<"\t asym:   "<<asym1<<endl;
    cout<<"chi1:  "<<chi1<<"\t  chi2:   "<<chi2<<"\t chi0:  "<<chi0<<endl;
    cout<<"============================================================="<<endl;
    cout<<"counter:     "<<counter<<"\t icon:  "<<icon<<endl;
    cout<<"lat:  "<<s.lat<<"\t lon:  "<<s.lon<<endl;
    cout<<"********************** SOURCE **************************"<<endl;
    cout<<"Ds:  "<<s.Ds<<"\t nums:  "<<s.nums<<"\t strucs:  "<<s.struc<<endl;
    cout<<"cl:  "<<s.cl<<"\t mass:  "<<s.mass<<"\t Tstar:  "<<s.Tstar<<endl;
    cout<<"Rstar:  "<<s.Rstar<<"\t logl:  "<<s.logl<<"\t type:  "<<s.type<<endl;
    cout<<"col:  "<<s.col<<"\t vs:  "<<l.vs<<"\t rho*:  "<<s.ros<<endl;
    cout<<"********************** LENS **************************"<<endl;
    cout<<"Configurations:  "<<l.con<<endl;
    cout<<"Dl:    "<<l.Dl<<"\t numl:  "<<l.numl<<"\t strucl:  "<<l.struc<<endl;
    cout<<"mass:  "<<l.Ml<<"vl:        "<<l.vl<<"\t l.RE:     "<<l.RE/AU<<endl;
    cout<<"Period(days):  "<<l.peri<<"\t q:  "<<l.q<<"\t dis:  "<<l.dis<<endl;
    cout<<"*********************** LENSING ***********************"<<endl;
    cout<<"tE:  "<<l.tE<<"\t RE(AU):  "<<l.RE/AU<<endl;
    cout<<"Vt:  "<<l.Vt<<"\t mul(mas/days):  "<<l.mul<<"\t u0:  "<<l.u0<<endl;
    cout<<"****************************************"<<endl;
    for(int i=0; i<M; ++i){
    cout<<"i:  "<<i<<"\t\t blend[i]:  "<<s.blend[i]<<"\t s.magb[i]:  "<<s.magb[i]<<"\t Nblend:  "<<s.Nblend[i]<<endl;}
    cout<<"==============================================================="<<endl;}
    
       // int eel;  cin>>eel; 
    
    }//end of EXtinction
    }//loop icon
    }///end of loop J
    return(0);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double RandR(double down,double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return(p*(up-down)+down);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double ErrorRoman(lens & l, double mag){
     int w=-1;
     if(mag<l.magw[0] or mag==l.magw[0]) w=0;
     else if(mag>l.magw[Nw-1] or mag==l.magw[Nw-1]) w=Nw-1;
     else{
     for(int i=1; i<Nw; ++i){
     if(double((mag-l.magw[i])*(mag-l.magw[i-1]))<0.0 or  mag==l.magw[i-1]){w=i-1; break;}}}
     if(w==-1){cout<<"mapt[0]: "<<mag<<"\t w: "<<w<<endl;  int ww; cin>>ww;}
     return(l.errw[w]); 
}
///==============================================================//6
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s){ 
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*Msun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opd=0.0;
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC;}
    s.opd= s.od_disk+s.od_ThD+s.od_bulge+s.od_halo;///total
    //cout<<"total_opticalD: "<<s.opd<<"\t od_disk: "<<s.od_disk<<endl;
   // cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
}
///==============================================================//
///                                                              //
///                  func_source   Initial amounts               //
///                                                              //
///==============================================================//
void func_source(source & s, CMD & cm, extinc & ex)
{
    int num,struc,nums,yye;
    double rho,rf;
    double Ds,Ai[M],Av;
    double Mab[M],Map[M];
    double maxnb=0.0;

    for(int i=0; i<M; ++i){
    s.Fluxb[i]=s.Nblend[i]=0.0;
    s.Nblend[i]=s.Nstart*pow(FWHM[i]*0.5,2)*M_PI/(3600.0*3600.0);
    if(s.Nblend[i]<=1.0) s.Nblend[i]=1.0;
    if(s.Nblend[i]>maxnb) maxnb=s.Nblend[i];}
    for(int i=0; i<M; ++i){s.magb[i]=s.Ai[i]=s.Map[i]=s.Mab[i]=0.0;}


    for(int k=1; k<=int(maxnb); ++k){

    do{
    num=int(RandR(10.0 , Num-5.0));
    rho=RandR(0.0,s.Romaxs);
    Ds=(double)(num*step);
    }while(rho>s.Rostari[num] or Ds<0.1 or Ds>MaxD);
    if(Ds>MaxD or Ds<0.1){cout<<"ERROR (1): Ds: "<<Ds<<"\t MaxD: "<<MaxD<<"\t step: "<<step<<"\t num: "<<num<<endl; cin>>yye;}
    nums=num;
    if(k==1){s.Ds=Ds;  s.nums=nums;}



     rf=RandR(0.0,s.Rostar0[nums]);
     if (rf<= s.rho_disk[nums])                    struc=0;///thin disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums])) struc=3;///halo
    if(k==1) s.struc=struc;


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(struc==0){//thin disk
    num=int(RandR(0.0,N1-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_d[i][num];}
    if(k==1){
    s.Rstar=cm.Rs_d[num];
    s.type= cm.type_d[num];
    s.mass= cm.mass_d[num];
    s.Tstar=cm.Teff_d[num];
    s.logl= cm.logl_d[num];
    s.col=  Mab[0]-Mab[1];
    s.cl=   cm.cl_d[num];}
    if(cm.mass_d[num]<0.0 or int(cm.cl_d[num])>5 or cm.Teff_d[num]<0.0 or cm.type_d[num]>8.0 or cm.mass_d[num]>1000.0){
    cout<<"Error(thin disk) temp: "<<cm.Teff_d[num]<<"\t mass: "<<cm.mass_d[num]<<"\t counter: "<<num<<endl; cin>>yye; }}


    if(struc==1){///bulge
    num=int(RandR(0.0,N2-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_b[i][num]; }
    if(k==1){
    s.Rstar=cm.Rs_b[num];
    s.type=cm.type_b[num];
    s.mass=cm.mass_b[num];
    s.Tstar=cm.Teff_b[num];
    s.col= Mab[0]-Mab[1];
    s.logl=cm.logl_b[num];
    s.cl= cm.cl_b[num];}
    if(cm.mass_b[num]<0.0 or int(cm.cl_b[num])>5 or cm.Teff_b[num]<0.0 or cm.type_b[num]>8.0 or cm.mass_b[num]>10000.0){
    cout<<"Error(bulge) temp: "<<cm.Teff_b[num]<<"\t mass: "<<cm.mass_b[num]<<"\t counter: "<<num<<endl;   cin>>yye; }}




    if(struc==2){///thick disk
    num=int(RandR(0.0,N3-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_t[i][num]; }
    if(k==1){
    s.Rstar= cm.Rs_t[num];
    s.type=cm.type_t[num];
    s.mass=cm.mass_t[num];
    s.Tstar=cm.Teff_t[num];
    s.col=Mab[0]-Mab[1];
    s.logl=cm.logl_t[num];
    s.cl= cm.cl_t[num];}
    if( cm.mass_t[num]<0.0 or int(cm.cl_t[num])>5 or cm.Teff_t[num]<0.0 or cm.type_t[num]>8.0 or cm.mass_t[num]>100000.0){
    cout<<"Error(thick disk) temp: "<<cm.Teff_t[num]<<"\t mass: "<<cm.mass_t[num]<<"\t counter: "<<num<<endl;  cin>>yye;}}



    if(struc==3){///stellar halo
    num=int(RandR(0.0,N4-2.0));
    for(int i=0; i<M; ++i){ Mab[i]=cm.Mab_h[i][num]; }
    if(k==1){
    s.Rstar=cm.Rs_h[num];
    s.type=cm.type_h[num];
    s.mass=cm.mass_h[num];
    s.Tstar=cm.Teff_h[num];
    s.col=Mab[0]-Mab[1];
    s.logl=cm.logl_h[num];
    s.cl= cm.cl_h[num];}
    if(cm.mass_h[num]<0.0 or int(cm.cl_h[num])>5 or cm.Teff_h[num]<0.0 or cm.type_h[num]>8.0 or cm.mass_h[num]>10000.0){
    cout<<"Error(Galactic halo) temp: "<<cm.Teff_h[num]<<"\t mass: "<<cm.mass_h[num]<<"\t counter: "<<num<<endl;  cin>>yye;}}

    Mab[4]=(Mab[2]+Mab[3]+Mab[4])/3.0;//absolute magnitude in W149: (K+H+J)/3
    ex.Aks=Interpol(Ds,ex);///extinction in Ks-band
    Av=ex.Aks*Avks;
    if(Av<0.0) Av=0.0;
    
    
    
    for(int i=0;  i<M; ++i){
    Ai[i]=fabs(Av*AlAv[i])+RandN(sigma[i], 1.5); //extinction in other bands
    if(Ai[i]<0.0) Ai[i]=0.0;
    Map[i]=Mab[i]+5.0*log10(Ds*100.0)+Ai[i];
    if(s.Nblend[i]>=k){s.Fluxb[i]+=fabs(pow(10.0,-0.4*Map[i]));}}
    if(k==1){
    for(int i=0; i<M;  ++i){s.Ai[i]=Ai[i]; s.Map[i]=Map[i]; s.Mab[i]= Mab[i];}
    s.col=s.col+s.Ai[0]-s.Ai[1];}
    }///loop over the stars

    for(int i=0; i<M; ++i){
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    s.magb[i]= -2.5*log10(fabs(s.Fluxb[i])); 
    
    if(int(s.Nblend[i])<1.0 or (s.Nblend[i]==1.0 and s.blend[i]<1.0) or s.blend[i]>1.0 or s.blend[i]<0.0 or
    s.Fluxb[i]<-0.009543 or s.Ai[i]<0.0){
    cout<<"BIGG ERRROR nsbl: "<<s.Nblend[i]<<"\t Nblend[i]: "<<s.Nblend[i]<<"\t belnd: "<<s.blend[i]<<endl;
    cout<<"BIG ERROR Flux is negative: "<<s.Fluxb[i]<<"\t No. filter:  "<<i<<"\t ext:  "<<s.Ai[i]<<endl; cin>>yye;} }
    
    if(s.type>8.0 or s.type<1.0 or s.mass<=0.0 or s.Tstar<0.0 or s.Rstar<0.0 or s.mass>10000.0 or
    s.nums>Num or s.nums<=0 or Av<0.0 or s.cl<0 or s.cl>=6){
    cout<<"ERROR(source):  type: "<<s.type<<"\t struc: "<<struc<<"\t num: "<<num<<"\t cl:  "<<s.cl<<endl;   cin>>yye;}
    //cout<<"************** End of func_source  ****************"<<endl;
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts               //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s, CMD & cm, extinc & ex)
{ 
    double test, fs;
    double rholens[s.nums+2];
    l.rhomaxl=0.0;
    int yye;
    double cosi; 
    double dwide, dclos, min1, dc, tt;
    double MlE, MlJ, MlS; 

    for(int k=1;k<=s.nums;++k){
    rholens[k]=0.0;
    l.Dl=k*step;
    l.xls=l.Dl/s.Ds;
    if(l.Dl>s.Ds){cout<<"ERROR (Dl>Ds) Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;  cin>>yye;}
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}
    

    do{
    l.numl=int(RandR(1.0,s.nums-2.0));
    test= RandR(0.0,l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){
    cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;  int ue; cin>>ue;}
    }while(test>rholens[l.numl]);
    l.Dl=double(l.numl*step);


   double randflag=RandR(0.0,s.Rostar0[l.numl]);
       if (randflag<=s.rho_disk[l.numl]) l.struc=0;///thin disk
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1; // bulge structure
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2; //thick disk
  else if (randflag<= s.Rostar0[l.numl]) l.struc=3;//halo
  else {cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}
  
 
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH   
    l.Ml=Mjupiter;
    l.q=RandR(-5.0,-2.0); 
    l.q=pow(10.0,l.q); 
    l.Mp=double(l.Ml*l.q);//kg
    MlE= l.Ml/Mearth; 
    MlJ= l.Ml/Mjupiter; 
    MlS= l.Ml/Msun; 
    if(MlE<2.0)      {l.Rl=pow(MlE,0.28)*Rearth;       }
    else if(MlJ<0.4) {l.Rl=0.8*pow(MlE,0.59)*Rearth;   } 
    else if(MlS<0.08){l.Rl=16.9*pow(MlE,-0.04)*Rearth; }
    else             {l.Rl=0.0014*pow(MlE,0.88)*Rearth;}    
    

    l.semi=RandR(3.0,100.0)*l.Rl; 
    l.u0=  RandR(0.0, 1.0);
    l.ksi= RandR(0.0, 2.0*M_PI);
    l.xls=double(l.Dl/s.Ds);
    l.RE=sqrt(4.0*G*(l.Ml+l.Mp)*s.Ds*KP)/velocity;//meter
    l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
    vrel(s,l);
    l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0);///in day
    s.ros=double(s.Rstar*Rsun*l.xls/l.RE);
    l.mul=l.Vt*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days
    if(l.u0==0.0) l.u0=1.0e-5; 
    l.tetE= double(l.RE/AU/l.Dl);///marcs 
    s.mus1= double(s.SV_n1- s.VSun_n1)*1000.0*3600.0*24.0/(s.Ds*AU);///mas/days
    s.mus2= double(s.SV_n2- s.VSun_n2)*1000.0*3600.0*24.0/(s.Ds*AU);///mas/days
    s.mul1= double(s.LV_n1- s.VSun_n1)*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days
    s.mul2= double(s.LV_n2- s.VSun_n2)*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days 
    l.piE=double(1.0/l.Dl -1.0/s.Ds)/l.tetE; 
    cosi=RandR(0.0,0.9999999);
    l.dis=l.semi*sqrt(1.0-cosi*cosi)/l.RE;
    l.peri= 2.0*M_PI*pow(l.semi,1.5)/sqrt(G*(l.Ml+l.Mp))/(3600.0*24.0);//days 
    
    //##################################################################  
    dwide=dclos=-1.0;   
    dwide=pow(1.0+pow(l.q,1.0/3.0),1.5)/sqrt(1.0+l.q); 
    min1=1000000000.0;
    if(l.q<0.0001)   dclos=double(1.0-3.0*pow(l.q , 1.0/3.0)/4.0);   
    else{ 
    for(int w=0; w<200; ++w){
    dc= double(0.6+w*(0.4/200.0)); 
    tt=fabs( pow(dc,8.0)*27.0*l.q - (1.0+l.q)*(1.0+l.q)*pow(1.0-dc*dc*dc*dc,3.0)); 
    if(tt<min1){min1=tt;  dclos=dc;  }}} 
    if(l.dis<=dclos)     l.con=1;  
    else if(l.dis>=dwide)l.con=3; 
    else                 l.con=2;   
    
    if(dwide<dclos or dwide==dclos or l.peri<0.0 or  l.dis<0.0 or l.q<0.00001){
    cout<<"Error dwide:    "<<dwide<<"dclose:  "<<dclos<<endl;
    cout<<"peri:  "<<l.peri<<"\t dis:  "<<l.dis<<"\t q:  "<<l.q<<endl;
    int uue ; cin>>uue; }
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm)
{
    //mass, Teff, Age, logL,  log(g),  Z,  Rs,  MB, MV, MI, MK, Cl, type (13)
    int yye, uui, k1, k2, h, g;
    double metal, age,gravity, MB;
    char filename[40];
    FILE *fp2;


    double Age1[YZ]={0.0}; double B1[YZ]={0.0};  double M1[YZ]={0.0};   double mm1[YZ]={0.0};
    double Age2[YZ]={0.0}; double B2[YZ]={0.0};  double M2[YZ]={0.0};   double mm2[YZ]={0.0};
    int number[70]={0};   int count[70]={0};   double Metal[70]={0.0};
    FILE *meta;
    meta=fopen("./files/CMD_WFIRST/metal.txt","r");
    for(int i=0; i<70; ++i){
    fscanf(meta,"%lf   %d  %d\n",&Metal[i],&count[i],&number[i]);
    if((Metal[i]<Metal[i-1] and i>0) or float(Metal[i])<-0.001 or number[i]==0 or count[i]>YZ or (abs(count[i-1]+number[i-1]-count[i])>2 and i>0)){
    cout<<"ERROR Metal[i]: "<<Metal[i]<<"\t count[i]: "<<count[i]<<"\t number[i]: "<<number[i]<<"\t i: "<<i<<endl; cin>>uui;} }
    fclose(meta);
    FILE *hks;
    hks=fopen("./files/CMD_WFIRST/HKS.txt", "r");
    for(int i=0; i<YZ; ++i){
    fscanf(hks,"%lf  %lf  %lf  %lf\n",&Age1[i],&mm1[i],&B1[i],&M1[i]);
    if(Age1[i]<0.0 or mm1[i]<0.0 or fabs(B1[i])>0.3 or M1[i]<0.5 or Age1[i]>18.0){
    cout<<"ERROR Age(HKS): "<<Age1[i]<<"\t metal: "<<mm1[i]<<"\t B[i]"<<B1[i]<<"\t M[i]: "<<M1[i]<<"\t i: "<<i<<endl; cin>>uui; }}
    fclose(hks);
    FILE *ji;
    ji=fopen("./files/CMD_WFIRST/JI.txt", "r");
    for(int i=0; i<YZ; ++i){
    fscanf(ji,"%lf   %lf   %lf  %lf\n",&Age2[i],&mm2[i],&B2[i],&M2[i]);
    if(Age2[i]<0.0 or mm2[i]<0.0 or fabs(B2[i])>1.7 or M2[i]<0.5 or Age2[i]>18.0  or Age1[i]!=Age2[i] or mm1[i]!=mm2[i]){
    cout<<"ERROR Age(JI): "<<Age2[i]<<"\t metal: "<<mm2[i]<<"\t B[i]"<<B2[i]<<"\t M[i]: "<<M2[i]<<"\t i: "<<i<<endl;
    cout<<"Age1[i]:  "<<Age1[i]<<"\t mm1[i]:  "<<mm1[i]<<endl;  cin>>uui;}}
    fclose(ji);

////=================================== THIN DISK ==============================
    int j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','T','i','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTiW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_d[j],&cm.Teff_d[j],&age,&cm.logl_d[j],&gravity,&metal,&cm.Rs_d[j],&MB,
    &cm.Mab_d[0][j],&cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.cl_d[j],&cm.type_d[j]);
    ///*******************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0])         h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl; cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_d[3][j]= double(B1[g]+M1[g]*cm.Mab_d[2][j]*1.05263157894737); ///H-band  Mks/Mk=1.05263157894737
    cm.Mab_d[4][j]= double(B2[g]+M2[g]*cm.Mab_d[1][j]);   ///J-band
    if(fabs(cm.Mab_d[3][j]-cm.Mab_d[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_d[4][j]-cm.Mab_d[1][j])>1.5){
    cout<<"ERROR: Mab_d(y-band): "<<cm.Mab_d[4][j]<<"\t Mab_d(z-band): "<<cm.Mab_d[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
    ///********************************************************
    if(cm.mass_d[j]<0.0 or cm.mass_d[j]==0.0 or cm.Teff_d[j]<0.0 or metal>0.1 or age>10 or
    int(cm.cl_d[j])==6 or float(cm.type_d[j])>8.0 or cm.type_d[j]<1.0){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_d: "<<cm.type_d[j]<<"\t CL(thin disk): "<<cm.cl_d[j]<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N1){
    cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;






////=================================== BULGE ==================================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c.dat",'C','M','D','b','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDbW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_b[j],&cm.Teff_b[j],&age,&cm.logl_b[j],&gravity,&metal,&cm.Rs_b[j],&MB,
    &cm.Mab_b[0][j],&cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.cl_b[j],&cm.type_b[j]);
    ///*****************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0]) h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl; cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_b[3][j]= double(B1[g]+M1[g]*cm.Mab_b[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_b[4][j]= double(B2[g]+M2[g]*cm.Mab_b[1][j]);   ///J-band
    if(fabs(cm.Mab_b[3][j]-cm.Mab_b[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_b[4][j]-cm.Mab_b[1][j])>1.5){
    cout<<"ERROR: Mab_b(y-band): "<<cm.Mab_b[4][j]<<"\t Mab_b(z-band): "<<cm.Mab_b[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///*****************************************************
    if(cm.mass_b[j]<0.0 or cm.mass_b[j]==0.0 or cm.Teff_b[j]<0.0 or age>10 or metal>0.9 or cm.cl_b[j]==6 or
    cm.type_b[j]>=8.0 or (cm.cl_b[j]==5 and int(cm.type_b[j])>8.0) or (cm.cl_b[j]==6 and int(cm.type_b[j])!=9) ){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_b: "<<cm.type_b[j]<<"\t CL(bulge): "<<cm.cl_b[j]<<endl;  cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
    cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;


////=================================== THICK DISK =============================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','T','k','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTkW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_t[j],&cm.Teff_t[j],&age,&cm.logl_t[j],&gravity,&metal,&cm.Rs_t[j],&MB,
    &cm.Mab_t[0][j],&cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.cl_t[j],&cm.type_t[j]);
    ///*********************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0]) h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl;  cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_t[3][j]= double(B1[g]+M1[g]*cm.Mab_t[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_t[4][j]= double(B2[g]+M2[g]*cm.Mab_t[1][j]);   ///J-band
    if(fabs(cm.Mab_t[3][j]-cm.Mab_t[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_t[4][j]-cm.Mab_t[1][j])>1.5){
    cout<<"ERROR: Mab_t(y-band): "<<cm.Mab_t[4][j]<<"\t Mab_t(z-band): "<<cm.Mab_t[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///********************************************************
    if(cm.mass_t[j]<0.0||  cm.mass_t[j]==0.0 or cm.Teff_t[j]<0.0 or metal>0.2||cm.cl_t[j]==6|| cm.type_t[j]>=8.0 or
    (cm.cl_t[j]==5 and int(cm.type_t[j])>8) or (cm.cl_t[j]==6 and int(cm.type_t[j])!=9)){
    cout<<"type_thick: "<<cm.type_t[j]<<"\t CL(thick): "<<cm.cl_t[j]<<endl;
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;


////=================================== STELLAR HALO ===========================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c.dat",'C','M','D','h','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDhW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.Teff_h[j],&age,&cm.logl_h[j],&gravity,&metal,&cm.Rs_h[j],&MB,
    &cm.Mab_h[0][j],&cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.cl_h[j],&cm.type_h[j]);
    ///************************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0]) h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl;  cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_h[3][j]= double(B1[g]+M1[g]*cm.Mab_h[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_h[4][j]= double(B2[g]+M2[g]*cm.Mab_h[1][j]);   ///J-band
    if(fabs(cm.Mab_h[3][j]-cm.Mab_h[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_h[4][j]-cm.Mab_h[1][j])>1.5){
    cout<<"ERROR: Mab_h(y-band): "<<cm.Mab_h[4][j]<<"\t Mab_h(z-band): "<<cm.Mab_h[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///***********************************************************
    if(cm.mass_h[j]<0.0 || cm.mass_h[j]==0.0 || age<0 or cm.cl_h[j]<0  or cm.cl_h[j]==6  or  cm.Teff_h[j]<0.0 or
    metal>0.1 || cm.cl_h[j]>7|| cm.type_h[j]>9 or (cm.cl_h[j]==5 and int(cm.type_h[j])>8) or (cm.cl_h[j]==6 and int(cm.type_h[j])!=9) or
    (cm.cl_h[j]<5 and int(cm.type_h[j])==9)){
    cout<<"type_halo: "<<cm.type_h[j]<<"\t CL(halo): "<<cm.cl_h[j]<<endl;
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
   cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
   cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;


}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nn=0.4/0.8;
   double mBarre;///stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   double fd=1.0; ///see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;///0.657066;////just stars brighter than V=11.5, but we change to consider all stars
   double fh=1.0;///No limitation
   double Rdd=2.17;///2.53;///2.17;
   double Rhh=1.33;///1.32;//1.33;
/*
I assume that the up-limit of the mass is indicated by the simulation. Because it depends on the structre, .... all details can be found in mass_averaged.cpp  code. */

   /*
    char filename[40];
    FILE *fill;
    sprintf(filename,"./files/density/%c%.2lf%c%.2lf.dat",'D',s.lat,'_',s.lon);
    fill=fopen(filename,"w");
    if(!fill){cout<<"cannot open file longtitude : "<<s.lon<<"\t latitude: "<<s.lat<<endl;  exit(0);}
   */



for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;

   x=i*step;
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = R_sun-x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///M_sun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-R_sun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nn)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-R_sun)/2.5)*exp(nn)*exp(-fabs(zb)/0.8)/(1.0+0.5*nn);///M_sun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(0.5/R_sun,-2.44);
   else            s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(rdi/R_sun,-2.44);///M_sun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///M_sun/pc^3
///=================================================
///در اینجا اینکه تعداد ستاره ه
///ا به این بستگی دارد که ما  نمودار قدر رنگ مطلق ستاره ها را چگونه درست کرده باشیم.اگر هیچ گونه محدودیتی برای
///درست کردن آن  در نظر نگرفته ایم،  پس تعداد کل ستاره ها را نظر  میگیریم.
/// ولی بهتر است که ما رابطه بین قدر مطلق و جرم را تعیین کنیم. در این صورت می توانیم  ورودی قدر رنگ
///وارد شده به کد را خودمان محدود به ستاره های روشن کنیم تا سرعت اجرای برنامه بالارود.
///averaged mass are the same as the previous work!!! because we did not change the besancon model


s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]);///[M_sun/pc^3]
s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[M_sun/deg^2]
s.Nstari[i]=binary_fraction*(s.rho_disk[i]*fd/0.403445+s.rho_ThD[i]*fh/0.4542+s.rho_halo[i]*fh/0.4542+s.rho_bulge[i]*fb/0.308571);////[Nt/pc^3]

s.Nstari[i]= s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]

s.Nstart  +=s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]
if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection
//fprintf(fill,"%e   %e   %e   %e   %e  %e   %e\n",x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.Rostar0[i],s.Nstari[i]);
   }
 // fclose(fill);
}
///===========================================================================//
double RandN(double sigma, double nnd){
   double rr,f,frand;
   do{
   rr=RandR(-sigma*nnd , sigma*nnd);
   f= exp(-0.5*rr*rr/sigma/sigma);
   frand=RandR(0.0,1.0); 
   }while(frand>f);
   return rr;
}
///==============================================================//
///                                                              //
///                  Linear interpolarion                        //
///                                                              //
///==============================================================//
double Interpol(double ds, extinc & ex){
  double F=-1.0;
  if(ds<ex.dis[0])        F=ex.Extks[0];
  else if(ds>=ex.dis[99]) F=ex.Extks[99];
  else{
  for(int i=0; i<99; ++i){
  if(ex.dis[i]>=ex.dis[i+1]){
  cout<<"ERROR dis[i]: "<<ex.dis[i]<<"\t disi+1: "<<ex.dis[i+1]<<endl;  int yye; cin>>yye; }
  if(ds>=ex.dis[i] && ds<ex.dis[i+1]){
  F = ex.Extks[i]+(double)(ds-ex.dis[i])*(ex.Extks[i+1]-ex.Extks[i])/(ex.dis[i+1]-ex.dis[i]);
  break;}}}
  if(F==-1.0||F<0.0){cout<<"ERROR big Extinction(ds): "<<F<<"\t ds: "<<ds<<endl; int uut; cin>>uut;}
  return(F);
}
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
int Extinction(extinc & ex,source & s)
{
     double sig,Lon,Lat;
     int uue, flag=0;
     if(s.lon<0.0){sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
     else sig=1.0;
     double delt=fabs(s.lon)-floor(fabs(s.lon));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
     else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
     else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
     else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
     else               Lon=(floor(fabs(s.lon))+0.75)*sig;
     if(s.lon==0.0 or Lon<0.0 or Lon<0.1  or s.lon<0.1) Lon=360.00;

     if(s.lat<0.0) sig=-1.0;
     else sig=1.0;
     delt=fabs(s.lat)-floor(fabs(s.lat));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
     else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
     else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
     else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
     else                Lat=(floor(fabs(s.lat))+0.75)*sig;
     if(Lon>360.000 || Lon<0.25 || fabs(Lat)>10.0 || (Lon>100 && Lon<260)){
     cout<<"BIG error (stopped program) s.lon: "<<Lon<<"\t s.lat: "<<Lat<<endl;   cin>>uue;}


     char filename[40];
     FILE *fpd;
     sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',Lat,'_',Lon);
     fpd=fopen(filename,"r");
     //cout<<"Lat:  "<<Lat<<"\t Lon:    "<<Lon<<endl;

     double lonti,latit;
     if(!fpd){
     cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
     FILE *SD;
     SD=fopen("./files/Ext/saved_direction.txt","r");
     for(int i=0; i<64881; ++i) {
     fscanf(SD,"%lf %lf \n",&latit,&lonti);
     if(fabs(Lat-latit)<0.1 && fabs(Lon-lonti)<0.1){
     cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
     cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
     flag=-1;}
     else{
     flag=1;
     for(int i=0; i<100; ++i){
     fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
     //cout<<"distance: "<<ex.dis[i]<<"\t Extks: "<<ex.Extks[i]<<endl;
     if(ex.dis[i]<0.2  || ex.dis[i]>50.0 || ex.Extks[i]<0.0){
     cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
     cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0; }}
     fclose(fpd);}
     //cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<flag<<endl;
     return(flag);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l)
{
  if (l.Dl==0.0) l.Dl=0.00034735;
  double pi=M_PI;
  double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*l.Dl*cos(s.TET)*cos(s.FI));
  double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*s.Ds*cos(s.TET)*cos(s.FI));
  if(Rlc==0.0) Rlc=0.0000000000034346123;
  if(Rsc==0.0) Rsc=0.000000000004762654134; 
 
  double LVx, SVx;
  double SVT, SVR, SVZ, LVT, LVR, LVZ;
  double fv, testfv, test, age;
  double  VSunx, vls2, vls1;
  double betal, betas, deltal, deltas, tetd ;



  double NN=2.5;
  double sigma_R_Disk,   sigma_T_Disk,  sigma_Z_Disk;
  double sigma_R_DiskL,  sigma_T_DiskL, sigma_Z_DiskL;
  double sigma_R_DiskS,  sigma_T_DiskS, sigma_Z_DiskS;
  double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
  double sigma_R_halo= 131.0, sigma_T_halo=106.0, sigma_Z_halo=85.0;
  double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0, sigma_Z_Bulge=100.0;
  double Rho[8]={00.0}; double maxr=0.0;
  for(int i=0; i<8; ++i){Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}

 
for (int i=0;i<2; ++i){
 test=RandR(0.0,maxr);
     if(test<=Rho[0])                       {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0; age= 0.075;}
else if(test<=(Rho[0]+Rho[1]))              {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0; age=0.575; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]))       {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;age=1.5;  }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3])){sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2; age=2.5; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]))       {sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8; age=4.0; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5])){sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4; age=6.0; }
else if(test<=maxr)                                       {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5; age=8.5; }
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}
    if(i==0) {
    sigma_R_DiskS= sigma_R_Disk;
    sigma_T_DiskS= sigma_T_Disk;
    sigma_Z_DiskS= sigma_Z_Disk;}
    if(i==1){
    sigma_R_DiskL= sigma_R_Disk;
    sigma_T_DiskL= sigma_T_Disk;
    sigma_Z_DiskL= sigma_Z_Disk; }}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.struc==0){///Galactic disk
    SVR= RandN(sigma_R_DiskS, NN);
    SVT= RandN(sigma_T_DiskS, NN);
    SVZ= RandN(sigma_Z_DiskS, NN); }

    else if(s.struc==1){///Galactic bulge
    SVR= RandN(sigma_R_Bulge, NN);
    SVT= RandN(sigma_T_Bulge, NN);
    SVZ= RandN(sigma_Z_Bulge, NN); }

    else if(s.struc==2){///thick disk
    SVR= RandN(sigma_R_TDisk, NN);
    SVT= RandN(sigma_T_TDisk, NN);
    SVZ= RandN(sigma_Z_TDisk, NN); }

    else if(s.struc==3){///stellar halo
    SVR= RandN(sigma_R_halo, NN);
    SVT= RandN(sigma_T_halo, NN);
    SVZ= RandN(sigma_Z_halo, NN); }
    if(s.struc==0 or s.struc==2)  SVT =SVT+ vro_sun*(1.00762*pow(Rsc/R_sun,0.0394) + 0.00712);
    l.vs=sqrt( SVR*SVR + SVT*SVT + SVZ*SVZ );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   if(l.struc==0){///Galactic disk
   LVR= RandN(sigma_R_DiskL, NN) ;
   LVT= RandN(sigma_T_DiskL, NN) ;
   LVZ= RandN(sigma_Z_DiskL, NN) ; }

   else if(l.struc==1){///Galactic bulge
   LVR= RandN(sigma_R_Bulge, NN) ;
   LVT= RandN(sigma_T_Bulge, NN) ;
   LVZ= RandN(sigma_Z_Bulge, NN) ; }

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN) ;
   LVT= RandN(sigma_T_TDisk, NN) ;
   LVZ= RandN(sigma_Z_TDisk, NN) ; }

   else if(l.struc==3){///stellar halo
   LVR= RandN(sigma_R_halo, NN);
   LVT= RandN(sigma_T_halo, NN);
   LVZ= RandN(sigma_Z_halo, NN); }
   if(l.struc==0 or l.struc==2)  LVT = LVT+vro_sun *(1.00762*pow(Rlc/R_sun,0.0394) + 0.00712);
   l.vl=sqrt( LVT*LVT + LVZ*LVZ + LVR*LVR );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  BETA  HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    betal=betas=0.0;
    tetd= s.TET;
    test= double(l.Dl*cos(s.FI)*sin(tetd)/Rlc);
    if(fabs(test-1.0)<0.01)       betal= pi/2.0;
    else if(fabs(test+1.0)<0.01)  betal=-pi/2.0;
    else                          betal=asin(test);
    
    test= double(s.Ds*cos(s.FI)*sin(tetd)/Rsc); 
    if( fabs(test-1.0)<0.01)     betas=pi/2.0;
    else if(fabs(test+1.0)<0.01) betas=-pi/2.0;
    else                         betas=asin(test);
    
    if(R_sun < fabs(l.Dl*cos(s.FI)*cos(tetd))) betal= pi-betal; 
    if(R_sun < fabs(s.Ds*cos(s.FI)*cos(tetd))) betas= pi-betas; 

   if(fabs(l.Dl*cos(s.FI)*sin(tetd))>Rlc or fabs(test)>1.0){
   cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
   cout<<"FI: "<<s.FI<<"\t TET: "<<tetd<<"\t betal:  "<<betal<<endl;
   cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<"\t betas:   "<<betas<<endl;
   cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(tetd)/Rlc<<"\t sin(s): "<<test<<endl;
   int ew; cin>>ew;}
///HHHHHHHHHHHHHHHHHHHHHHHHHH  DELTA   HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.TET>pi)  tetd=s.TET-2.0*pi; 
    deltal= pi - fabs(tetd) -fabs(betal);
    deltas= pi - fabs(tetd) -fabs(betas);  
    if(betal<0.0)  deltal= -1.0*deltal;
    if(betas<0.0)  deltas= -1.0*deltas;
    s.deltao= double( pi-fabs(tetd) );
    if(tetd<0.0)  s.deltao=-1.0*s.deltao; 

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    s.SV_n1 =+SVR * sin(deltas)- SVT * cos(deltas);
    s.LV_n1 =+LVR * sin(deltal)- LVT * cos(deltal);
    s.VSun_n1=+VSunR*sin(s.deltao)-VSunT*cos(s.deltao);
    
    SVx= -SVR*cos(deltas)- SVT*sin(deltas);
    LVx= -LVR*cos(deltal)- LVT*sin(deltal);
    VSunx= -VSunR*cos(s.deltao) -VSunT*sin(s.deltao);
    
    s.SV_n2=-sin(s.FI)*(SVx) + cos(s.FI)*SVZ;
    s.LV_n2=-sin(s.FI)*(LVx) + cos(s.FI)*LVZ;
    s.VSun_n2=-sin(s.FI)*(VSunx)+cos(s.FI)*(VSunZ);
 
    
    vls1= l.xls*s.SV_n1 - s.LV_n1 +(1.0-l.xls)*s.VSun_n1;  ///Source - lens 
    vls2= l.xls*s.SV_n2 - s.LV_n2 +(1.0-l.xls)*s.VSun_n2;  /// Source -lens
    l.Vt=sqrt(fabs( vls1*vls1 + vls2*vls2 ) );
    

    if(vls1>0.0 and vls2>=0.0)         s.xi=atan(fabs(vls2)/fabs(vls1));
    else if(vls1<=0.0 and vls2>0.0)    s.xi=atan(fabs(vls1)/fabs(vls2)) + M_PI/2.0;
    else if(vls1<0.0 and vls2<=0.0)    s.xi=atan(fabs(vls2)/fabs(vls1)) + M_PI;
    else if(vls1>=0.0 and vls2<0.0)    s.xi=atan(fabs(vls1)/fabs(vls2)) + 3.0*M_PI/2.0;
    else if(s.xi>2.0*M_PI)             s.xi-= 2.0*M_PI; 
    
    if(vls1==0.0 and vls2==0.0){
    cout<<"Big error both zeros:  "<<vls1<<"\t vls_b:  "<<vls2<<endl;
    int iee;  cin>>iee;}

    if (l.Vt<0.0 or l.Vt>1.0e6){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;
    cout<<"source type:  "<<s.type<<"\t s.struc:  "<<s.struc<<endl;
    int yee; cin>>yee;}
    //cout<<"(vrel_func):   s.deltao:  "<<s.deltao<<"FI: "<<s.FI<<endl;
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
