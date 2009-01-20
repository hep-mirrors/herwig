/* This program generates top pairs at NLO accuracy from polarized e+e- annihilation 
at the ILC with C.O.M energy of 500 GeV (an update for arbitary energy will be presented soon).
 The events are produced in the LH xml format */
#include <cstdlib> 
#include <ctime> 
#include <iostream>
#include <fstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <iomanip>
#include "ILCTTBAR500_INPUTS.h"
using namespace std;

void filegen(double &rintgend1a, double &rintgend2a,double &rintgend3a, double &rintgend4a,double &rintgend1v, double &rintgend2v,double &rintgend3v, double &rintgend4v); // generates production interpolation files

 void filegendec (double &rintgend1, double &rintgend2, double &rintgend3, double &rintgend4);// generates decay interpolation files

double interp(int N,double XXX, double Xmax1, double Xmax2, double Xmax3, double Xmax4, bool axi); // interpolates for production process

 double interpd(int N,double XXX, double Xmax1, double Xmax2, double Xmax3, double Xmax4);// interpolates for decay process

void xmaxmin (double pt2, double d2, double &rxmax, double &rxmin5,double &rxmin6); // finds maximum and minimum values for variable x in production process

double ME (double d2, double pt2,double &rx, double &ry1, double &ry2,double xmin5, double &rM1, double &rM2, double &rw, bool axi); //calculates matrix element  for production process

 double MEd (double pt2,double &rx3, double &rx11, double &rx12, double &rMd1, double &rMd2);//calculates matrix element for decay process

 void PolarizedME (double xt, double xtb, double xg, double &rct, double &rst, double &rsphi, double &rcphi,double &rLLL,double &rLRR,double &rLRL,double &rLLR,double &rRLL,double &rRRR, double &rRRL,double &rRLR);//calculates polarized matrix element for production process

 double PolMEdec(double xw1, double xb1, double xg1,double &rctw, double &rctl,double &rbphi, double &rbgphi,int tm, double xd0, double zd0);//calculates polarized matrix element for decay process

double random (int &rseed); // random number generator

bool ax(int rseed1); // axial-vector or vector?

/***************************************************************************************************/
// PARAMETERS
   Input user;
   double pi = 3.14152654;
   double b = 23./(12.*pi);
   double bp = 58./(46.*pi);
   double mt =175.;
   double emcm = 500.;
   double ms = 2.*mt/emcm;
   double rho = pow(mt/emcm,2);
   double d2 = pow(ms,2);
   double beta=sqrt(1.-d2);
   double CF=4./3.;
// Max and min kt^2 for production. Note this is kappa. Hard wired for 500 GeV
   double pt2max = 0.03;
   double pt2min=0.00000016;
   double mw=80.;
   double a =pow(mw/mt,2);
// Max and min kt^2 for decays. Note this is kappa*mt^2
   double pt2maxd = pow((1.-sqrt(a))/2.,2)*pow(mt,2);
   double pt2mind = 0.04;
   double alphaem=0.007297352;
   double e2 = 4.*pi*alphaem;
   double s=pow(emcm,2);
   double Mz=91.2;
   double Yz=2.486;
   double GF=0.0000116639;
   double sinsq2thw=0.2024;
   double sinsqthw=0.2312;
   double costhw=0.8768; 
   double k = sqrt(2.)*GF*pow(Mz,2)/e2;
   double At=0.5;
   double Qt = 2./3.;
   double Vt=0.5-2.*Qt*sinsqthw;
   double Ae=-0.5;
   double Ve=-0.5+2.*sinsqthw;
   double Ref=s*(s-pow(Mz,2))/(4.*sinsq2thw*(pow((s-pow(Mz,2)),2)+pow(Mz*Yz,2)));
   double Imf=-Mz*Yz*s/(4.*sinsq2thw*(pow((s-pow(Mz,2)),2)+pow(Mz*Yz,2)));
   double ve = 2.*Ve;
   double ae = 2.*Ae;
   double vt = 2.*Vt;
   double at = 2.*At;
   int Pm =  user.Pem();
   int Pp =  user.Pep();
   double qel=(2.*sinsqthw-1.)/(2.*costhw);
   double qer=sinsqthw/costhw;
   double qtl=(3.-4.*sinsqthw)/(6.*costhw);
   double qtr= -2.*sinsqthw/(3.*costhw);
   double c1=3.5;
   double d1=2.25;
   
   double chi1 = k*s*(s-pow(Mz,2))/(pow((s-pow(Mz,2)),2)+pow(Yz*Mz,2));
   double chi2 = pow(k,2)*pow(s,2)/(pow((s-pow(Mz,2)),2)+pow(Yz*Mz,2));
     double alphaprod =1./(b*log(1./pt2min))*(1.-bp*log(log(1./pt2min))/(b*log(1./pt2min)));
   
   double sV=pow(Qt,2)-2.*Qt*Ve*Vt*chi1+(pow(Ae,2)+pow(Ve,2))*pow(Vt,2)*chi2;
   double sA=(pow(Ae,2)+pow(Ve,2))*pow(At,2)*chi2;
   double A=(pow(beta,3)+d1*alphaprod/pi)*sA*4.*pi*pow(alphaem,2)/s;
   double V=(beta*(1.+2.*pow(mt,2)/s)+c1*alphaprod/pi)*sV*4.*pi*pow(alphaem,2)/s;


//cross-section in picobarns
  double AV=(A+V)*389400000.;
    
// Parameters and integrals for Bbar production and decay born variables
double asmt = 0.107; // alphas at top mass
   double xg0=0.2*2./(pow(beta,2)*emcm);
    double F1Vg=2./3.;
double  F1Vz=(0.25-2.*sinsqthw/3.)/(sqrt(sinsqthw)*costhw);
double F1Az=(-0.25)/(sqrt(sinsqthw)*costhw);
double I1virt=-2.+(1.+pow(b,2))/(2.*b)*(-1.5+log(4.*pow(b,2)/(1.-pow(b,2)))*log((1.-b)/(1.+b))+2.*pow(pi,2)/3.+2*0.174285939+0.5*pow((log((1.-b)/(1.+b))),2));
double I2=(1.-pow(b,2))/(4.*b)*(log((1.-b)/(1.+b)));
double I1soft=-2.*log(xg0)*(1.+(1.+pow(b,2))/(2.*b)*log((1.-b)/(1.+b)))+log((1.-pow(b,2))/(4.*pow(b,2)))-1./b*log((1.-b)/(1.+b))+(1.+pow(b,2))/(2.*b)*(-log(pow(b,2))*log((1.-b)/(1.+b))-pow(pi,2)/3+2.*0.174285939+0.5*pow((log((1.-b)/(1.+b))),2));
//NB: Hard code ddilog((1-b)/(1+b))=0.174285939
double I11=I1virt+I1soft;
     double dF1Vg=asmt*4./(3.*2.*pi)*F1Vg*(I11+I2);
    double dF1Vz=asmt*4./(3.*2.*pi)*F1Vz*(I11+I2);
    double dF1Az=asmt*4./(3.*2.*pi)*F1Az*(I11-I2);

/***************************************************************************************************/

int main() {
 
  // Switches
  int nevg = user.nevgen();
  double emcm = user.cme();
  int seed = user.rseed();
  bool ptruncate = user.truncpro();
  bool dtruncate = user.truncdec();
  bool POWprod = user.POWHEGprod();
  bool POWdecay = user.POWHEGdecay();
 //=========================================================//
  //Parameters for production and decay Bbar born variables
   bool findmaxprod=false;
   bool findmaxdec=false;
   double int1=0.,int2=0.,int3=0.,int4=0.;
   double max1=0.,max2=0.,max3=0.,max4=0.;
   double maxd1=0.,maxd2=0.,maxd3=0.,maxd4=0.;
 //=========================================================//
  bool axi = ax(seed);
  bool emit; // true if truncated emission in production
  bool demitt, demittb; // true if truncated emission in t and tbar decay respectively
  double qhp = emcm; //Scale of hardest production emission for use in decay truncated shower
  int EVT = 0; //Final state polarization identifier
 
  double T1[6],T2[6],PUP[5][18],AQEDUP,AQCDUP,XWGTUP,SCALUP,mass[18];
  int ICOLUP[3][18],MOTHUP[3][18],ISTUP[18],ISPINUP[18],IDUP[18],NUP;
  
   int ii=0; //event count
 
   //Orientation sines and cosines for top pairs in production and decays

  double cost, sint, costz,costx,costbz,costbx,sintx,sintbx,sintbz,sintz;
  costz = 0.;costbz = 0.; costx = 0.;costbx = 0.;
  sintz = 0.;sintbz = 0.; sintx = 0.;sintbx = 0.;
  double xg1,xg2,xw11,xw12,pt2r1,pt2r2;
  xg1 = 0.; xg2 = 0.;xw11 = 0.; xw12 = 0.;pt2r1=0.; pt2r2 = 0.;
  double xgrt, xgrtb, xwrt, xwrtb,ptr2,ptr2b;
  xgrt = 0.; xgrtb = 0.;ptr2=0.; ptr2b = 0. ;
  // Integrals at end of pt2-integral tables:production
   double end1a,end2a,end3a,end4a,end1v,end2v,end3v,end4v ;

 /*********************************START PRODUCTION PROCESS*************************************/
   if (POWprod) {
  
  filegen(end1a,end2a, end3a, end4a,end1v,end2v,end3v,end4v);}

// Integrals at end of pt2-integral tables:decay
   double dend1,dend2,dend3,dend4;
   if (POWdecay) {
    
    filegendec (dend1, dend2, dend3, dend4);}
  
  
  int mm=0;
  ofstream outdata;
  bool xml = true ;
  /**********************************************************************************************/
    outdata.open("ILCTTBAR500.dat", ios::trunc); 
   outdata << "<LesHouchesEvents version =\"1.0\">" << endl; outdata << "<!--" << endl;
    if (Pm == 1) {outdata << "RH e-" << "\t" ; } else {outdata << "LH e-" << "\t" ; }
     if (Pp == 1) {outdata << "RH e+" << endl; } else {outdata << "LH e+" ; }
     outdata << "annihilation at the ILC" << endl;
     outdata << "File generated with ILCTTBAR500.cxx" << endl; outdata << "-->" << endl;
    outdata << "<init>" << endl;
     outdata << "\t11\t" << "-11\t" <<"\t" << user.cme()/2. << "\t" << user.cme()/2.  << "\t" << "0 \t 0 \t 7\t 7 \t 1 \t 1" << endl;
   outdata << "\t" << AV  << "\t" << 0.000000 << "\t1.00000 \t11" << endl; 
    outdata << "</init>" << endl;
    /**********************************************************************************************/
   
  for (int ixx=1; ixx < nevg+1; ixx++) {
    double x1, x2;
    if (POWprod) {
    bool iterate =false;
    emit = false;   
    double intg=0.;
   
    
    for (int ia=0; iterate==false; ia++){
      double logn=fabs(log(random(seed)));
      intg+=logn;
      double end1=end1v, end2=end2v, end3=end3v, end4=end4v;
      if(axi) {end1=end1a, end2=end2a, end3=end3a, end4=end4a;}
       double pt2=interp(1000,intg,end1,end2,end3,end4,axi);
      
      if (0.*intg !=0) {break;}
      if (intg > end4) {break;}
      if (pt2 < pt2min) {break;}
      double xmax, xmin5, xmin6;

       xmaxmin(pt2,d2,xmax, xmin5, xmin6);
       double x, y1, y2, M, M1, M2, w;
       double Mmax=0.;
       for (int ix=0; ix < 1000 ; ix++) {
       x=(xmax-xmin6)*random(seed)+xmin6;
       M= ME(d2,pt2,x,y1,y2,xmin5,M1,M2,w,axi);
     
       if (fabs(M) > fabs(Mmax)) {Mmax=M;}
       }
	do {
       x=(xmax-xmin6)*random(seed)+xmin6;
       M = ME(d2,pt2,x,y1,y2,xmin5,M1,M2,w,axi);} while (random(seed) > fabs(M/Mmax));
 
	int Jt; // Identifier for solution
	if (x > xmin5 && pt2 < 0.024){
	  if (random(seed) < (fabs(M1)/fabs(M))) {
	      //Kleiss trick
	       if (random(seed) <  pow(x,2)/(pow(x,2)+pow(y1,2))) {
		x1=x;
		x2=y1;
		Jt=1;}
	      else{
	        x2=x;
	        x1=y1;
		Jt=2;}
	        mm=mm+1;
	   
	      iterate=true;
	     
	    }
          else{
	   //Kleiss trick
	      if (random(seed) <  pow(x,2)/(pow(x,2)+pow(y2,2))) {
		x1=x;
		x2=y2;
		Jt=1;}
	      else{
	        x2=x;
	        x1=y2;
	        Jt=2;}
	       mm=mm+1;
	       iterate=true;}}
	else {
	   //Kleiss trick
	    if (random(seed) <  pow(x,2)/(pow(x,2)+pow(y2,2))) {
	         x1=x;
		 x2=y2;
		 Jt=1;}
           else{
	        x2=x;
	        x1=y2;
	        Jt=2;}
	        mm=mm+1;
	        iterate=true;  }
      if (iterate) {
	
	double z1, ptr2, ptsq ,z, q;
	z1 =0.; ptr2 =0.;
	double x3=2.-x1-x2;
	double Qg=0.75;
	double u=mt;
	double qc=2.8;
	double CF =4./3.;
	double alpham=1.;
	ptsq=pow(emcm,2)*pt2;
	double l=sqrt(1.-d2);
	double ac, ab, ag;
	if (Jt==2) {
	ac=(x2*(1.+l)-sqrt(pow(x2*(1.+l),2)-2.*d2*(1.+l-d2/2.)))/(2.*(1.+l-d2/2.));
	ab=0.5*(x1*(1.+l)+sqrt(pow((1.+l)*x1,2)-8.*(d2/4.+pt2)*(1.+l-d2/2.)))/((1.+l-d2/2.));
	ag=2./(1.+l)-ab-ac;
	z=ab/(ab+ag);
      }else {
	ac=(x1*(1.+l)-sqrt(pow(x1*(1.+l),2)-2.*d2*(1.+l-d2/2.)))/(2.*(1.+l-d2/2.));
	ab=0.5*(x2*(1.+l)+sqrt(pow((1.+l)*x2,2)-8.*(d2/4.+pt2)*(1.+l-d2/2.)))/((1.+l-d2/2.));
	ag=2./(1.+l)-ab-ac;
	z=ab/(ab+ag);
      }

	/**********************START PRODUCTION TRUNCATED SHOWER****************************/
	// Scale of hardest emission.
	qhp=sqrt(ptsq/pow(z*(1.-z),2)+pow(u/z,2)+pow(Qg/(1.-z),2)/z);
	// Truncated Shower for production
      if (ptruncate) {
	double qi=emcm;

       
	if (ptsq/pow(z*(1.-z),2)+pow(u/z,2)+pow(Qg/(1.-z),2)/z < 0.) {emit=false;}
	double qh=sqrt(ptsq/pow(z*(1.-z),2)+pow(u/z,2)+pow(Qg/(1.-z),2)/z);
        double C=CF*alpham*log(Qg/(qi-u))/pi;
	  double D=exp(2.*log(qh/qc)*C);
	  double E=exp(2.*log(qi/qc)*C);
	 if (random(seed) < 1.-E/D) {
	  do {
	  q=sqrt(exp(log(E/random(seed))/C)*pow(qc,2));} while (q < qh);
	  double Pzmax=(1.+pow((1.-Qg/emcm),2))/(Qg/emcm)-2.*pow(u,2)/(Qg/emcm*(1.-Qg/emcm)*pow(q,2));
	  double Pz;
	  do{
	    z1=random(seed)*(1.-Qg/emcm-u/emcm)+u/emcm;
	    Pz=(1.+pow(z1,2))/(1.-z1)-2.*pow(u,2)/(z1*(1.-z1)*pow(q,2));} while (random(seed) > Pz/Pzmax);

	  ptr2 = pow((1.-z1)*z1*q,2)-pow((1-z1)*u,2)-z1*pow(Qg,2);
	  emit = true;
	   double kapp=pow(sqrt(ptr2)/emcm,2);
	  //energy fractions
	  double ze=ab+2.*((d2/4.+kapp)/ab-d2/4.*ab)/(1.+l);
	  double ab1=z1*(ab+ag);
	  double z1e=ab1+2.*((d2/4.+kapp)/ab1-d2/4.*ab1)/(1.+l);  
	  if (Jt==1) {if ((1.-z1e)*(x2+x3)*emcm/2.< sqrt(ptr2)){emit=false;}
	  if (ze*z1e*(x2+x3)*emcm/2.< sqrt(ptsq))    {emit=false;}}
	  if (Jt==2) {if ((1.-z1e)*(x1+x3)*emcm/2.< sqrt(ptr2)){emit=false;}
	  if (ze*z1e*(x1+x3)*emcm/2.< sqrt(ptsq))    {emit=false;}}
	  double alphas = 1./(b*log(ptr2/(pt2min*pow(emcm,2))))*(1.-bp*log(log(ptr2/(pt2min*pow(emcm,2))))/(b*log(ptr2/(pt2min*pow(emcm,2)))));
	  alphas=alphas+pow(alphas,2)*(67./3.-pow(pi,2)-50./9.)/(4.*pi);
	  if (ptr2 < 0 || q > qi || ptr2 > ptsq || random(seed) > alphas/alpham || random(seed) > Pz*(1.-z1)/2.) {emit=false;}
	 
	 
          z1=z1e;
	 }
	
      }
      /**********************END PRODUCTION TRUNCATED SHOWER****************************/
       // Polarizer

      double LLL,LRR,LRL,LLR,RLL,RRR,RRL,RLR,ME;
      ME =0.; 
      double xt, xtb, xg, ct,st,sphi,cphi,phig,maxprod;
      if (!findmaxprod) {
       for (int im=0; im < 10000 ; im++) {
      do { 
      xt=random(seed);
      xtb=random(seed);
      xg=2.-xt-xtb;} while (xg < xg0 || (1.-xt)*(1.-xtb)*(xt+xtb-1.) < pow(mt/emcm*xg,2));
      phig=2.*random(seed)*pi;
      cphi=cos(phig);
      sphi=sin(phig);
      ct=1.-2.*random(seed);
      st=sin(acos(ct));
      PolarizedME (xt,xtb,xg,ct,st,sphi,cphi,LLL,LRR,LRL,LLR,RLL,RRR,RRL,RLR);
      if (Pm==-1) {if(fabs(LRR) > fabs(max1)) {max1=LRR;}
		   if(fabs(LLR) > fabs(max2)) {max2=LLR;}
		   if(fabs(LRL) > fabs(max3)) {max3=LRL;}
		   if(fabs(LLL) > fabs(max4)) {max4=LLL;}
		   int1+=LRR;
		   int2+=LLR;
		   int3+=LRL;
		   int4+=LLL;
      } else if (Pm==1) {if(fabs(RRR) > fabs(max1)) {max1=RRR;}
		   if(fabs(RLR) > fabs(max2)) {max2=RLR;}
		   if(fabs(RRL) > fabs(max3)) {max3=RRL;}
		   if(fabs(RLL) > fabs(max4)) {max4=RLL;}
		   int1+=RRR;
		   int2+=RLR;
		   int3+=RRL;
		   int4+=RLL;
      }
       }

      }
      findmaxprod=true;
  
      double intt=int1+int2+int3+int4;
      if(random(seed)< (int1+int2)/intt) { 
	if(random(seed) < int1/(int1+int2)) {
	  EVT=1;
	  maxprod=max1;
} else {
            EVT=2;
	    maxprod=max2;} }
      else {if(random(seed) < int3/(int3+int4)) {
	EVT=3;
	maxprod=max3;} else {
	  EVT=4;
	  maxprod=max4;	}}
    
      do{
      do {
      xt=random(seed);
      xtb=random(seed);
      xg=2.-xt-xtb;} while(xg < xg0 || (1.-xt)*(1.-xtb)*(xt+xtb-1.) < pow(mt/emcm*xg,2));
      phig=2.*random(seed)*pi;
      cphi=cos(phig);
      sphi=sin(phig);
      ct=1.-2.*random(seed);
      st=sin(acos(ct));
      PolarizedME (xt,xtb,xg,ct,st,sphi,cphi,LLL,LRR,LRL,LLR,RLL,RRR,RRL,RLR);
      if (Pm==-1) {if(EVT==1) {ME=LRR;}
		   if(EVT==2) {ME=LLR;}
		   if(EVT==3) {ME=LRL;}
		   if(EVT==4) {ME=LLL;}
      } else if (Pm==1) {if(EVT==1) {ME=RRR;}
                   if(EVT==2) {ME=RLR;}
		   if(EVT==3) {ME=RRL;}
		   if(EVT==4) {ME=RLL;}
      }

      }while (random(seed)> fabs(ME/maxprod));

      double maxphi=0.;
      double mg=0.75;
       xg=2.-x1-x2;
       cphi=1.;
       sphi=1.;
      PolarizedME (x1,x2,xg,ct,st,cphi,sphi,LLL,LRR,LRL,LLR,RLL,RRR,RRL,RLR);
      if (Pm==-1) {if(EVT==1) {maxphi=LRR;}
		   if(EVT==2) {maxphi=LLR;}
		   if(EVT==3) {maxphi=LRL;}
		   if(EVT==4) {maxphi=LLL;}
      } else if (Pm==1) {if(EVT==1) {maxphi=RRR;}
                   if(EVT==2) {maxphi=RLR;}
		   if(EVT==3) {maxphi=RRL;}
		   if(EVT==4) {maxphi=RLL;}
      }
      cphi=-1.;
      PolarizedME (x1,x2,xg,ct,st,cphi,sphi,LLL,LRR,LRL,LLR,RLL,RRR,RRL,RLR);
      if (Pm==-1) {if(EVT==1) {if(LRR > maxphi){maxphi=LRR;}}
		   if(EVT==2) {if(LLR > maxphi){maxphi=LLR;}}
		   if(EVT==3) {if(LRL > maxphi){maxphi=LRL;}}
		   if(EVT==4) {if(LLL > maxphi){maxphi=LLL;}}
		   } else if (Pm==1) {if(EVT==1) {if(RRR > maxphi){maxphi=RRR;}}
                   if(EVT==2) {if(RLR > maxphi){maxphi=RLR;}}
		   if(EVT==3) {if(RRL > maxphi){maxphi=RRL;}}
		   if(EVT==4) {if(RLL > maxphi){maxphi=RLL;}}

      }
      do {
      phig=2.*random(seed)*pi;
      cphi=cos(phig);
      sphi=sin(phig);
      PolarizedME (x1,x2,xg,ct,st,cphi,sphi,LLL,LRR,LRL,LLR,RLL,RRR,RRL,RLR);
      if (Pm==-1) {if(EVT==1) {ME=LRR;}
		   if(EVT==2) {ME=LLR;}
		   if(EVT==3) {ME=LRL;}
		   if(EVT==4) {ME=LLL;}
      } else if (Pm==1) {if(EVT==1) {ME=RRR;}
                   if(EVT==2) {ME=RLR;}
		   if(EVT==3) {ME=RRL;}
		   if(EVT==4) {ME=RLL;}
      }}while(random(seed)> fabs(ME/maxphi));

      // Assign momenta to production process

       double ppt,ctt,cf;
       double x=x2+x3;
       //energy fraction z
       z=x2/x;
       if(Jt==2) {x=x1+x3; z=x1/x;}
       double px1 = sqrt(pow(0.5*emcm*x1,2)-pow(mt,2));
       double px2 = sqrt(pow(0.5*emcm*x2,2)-pow(mt,2));
       double c12=(x1*x2-2.*(1.-x3)+d2)/sqrt((pow(x1,2)-d2)*(pow(x2,2)-d2));
       double c13=(x1*x3-2.*(1.-x2))/(sqrt(pow(x1,2)-d2)*x3);   
       double c23=(x2*x3-2.*(1.-x1))/(sqrt(pow(x2,2)-d2)*x3);
     
          if (Jt==1) {
    cf = -c13; ctt = -c12;
    }else
    {cf =-c23; ctt =-c12;}
	
       double phi=2.*pi*random(seed);
            if (Jt==1) {
      double x=x2+x3;
      PUP[3][4]=px1*ct;
      double pt = sqrt(pow(px1,2)-pow(PUP[3][4],2));
      PUP[1][4]=pt*cos(phi);
      PUP[2][4]=pt*sin(phi);
      double ppl,ppl2,ppe1,ppe2; 
      ppl2 =0.; ppe1 =0.;ppe2 = 0.;
      if (!emit){
	ppl = px2*ctt;
	if (random(seed) < 0.5) {
	 ppt = sqrt (pow(px2,2)-pow(ppl,2));}else {ppt = -sqrt (pow(px2,2)-pow(ppl,2));}
} else { 
	  ppe1=z1*x*emcm/2.;
	  ppe2=(1.-z1)*x*emcm/2.;
	  ppl2=sqrt(pow(ppe2,2)-ptr2);
	  ppl=sqrt(pow(ppe1,2)-ptr2);
	  if (ppe1 < sqrt(ptr2)) {
	    ppl=0.; ppe1=sqrt(ptr2); ppe2=x*emcm/2.-ppe1; ppl2=sqrt(pow(ppe2,2)-ptr2);}
	  if (ppe2 < sqrt(ptr2)) {
	    ppl2=0.; ppe2=sqrt(ptr2); ppe1=x*emcm/2.-ppe2; ppl=sqrt(pow(ppe1,2)-ptr2);}
	  if (random(seed) < 0.5) {
	 ppt = sqrt(ptr2);}else {ppt = -sqrt(ptr2) ;}
}
      double ppt1 = ppt*cos(phig);
      double ppt2 = ppt*sin(phig);
      PUP[3][5]=-(ppl*ct-ppt1*sin(acos(ct)));
      PUP[1][5]=-((ppl*sin(acos(ct))+ppt1*ct)*cos(phi)-ppt2*sin(phi));
      PUP[2][5]=-((ppl*sin(acos(ct))+ppt1*ct)*sin(phi)+ppt2*cos(phi));
       if (emit){
	PUP[4][4] = x1*emcm/2.;
	PUP[4][5]=ppe1;
	PUP[4][6]=ppe2;
      PUP[3][6]=-(ppl2*ct+ppt1*sin(acos(ct)));
      PUP[1][6]=-((ppl2*sin(acos(ct))-ppt1*ct)*cos(phi)+ppt2*sin(phi));
      PUP[2][6]=-((ppl2*sin(acos(ct))-ppt1*ct)*sin(phi)-ppt2*cos(phi));
     }
          }
    else {
         ct = -ct;
      double x=x1+x3;
      PUP[3][5]=px2*ct;
      double pt = sqrt(pow(px2,2)-pow(PUP[3][5],2));
      PUP[1][5]=pt*cos(phi);
      PUP[2][5]=pt*sin(phi);
      double ppl, ppl2, ppe1,ppe2;
      ppl2 =0.; ppe1 =0.;ppe2 = 0.;
      if (!emit){
       ppl = px1*ctt;
       if (random(seed) < 0.5) {
	ppt = sqrt (pow(px1,2)-pow(ppl,2));
	}else {ppt = -sqrt (pow(px1,2)-pow(ppl,2));}} else {
	  ppe1=z1*x*emcm/2.;
	  ppe2=(1.-z1)*x*emcm/2.;
	  ppl2=sqrt(pow(ppe2,2)-ptr2);
	  ppl=sqrt(pow(ppe1,2)-ptr2);
	  if (ppe1 < sqrt(ptr2)) {
	    ppl=0.; ppe1=sqrt(ptr2); ppe2=x*emcm/2.-ppe1; ppl2=sqrt(pow(ppe2,2)-ptr2);}
	  if (ppe2 < sqrt(ptr2)) {
	    ppl2=0.; ppe2=sqrt(ptr2); ppe1=x*emcm/2.-ppe2; ppl=sqrt(pow(ppe1,2)-ptr2);}
	  if (random(seed) < 0.5) {
	 ppt = sqrt(ptr2);}else {ppt = -sqrt(ptr2) ;}
	}
      
      double ppt1 = ppt*cos(phig);
      double ppt2 = ppt*sin(phig);
      PUP[3][4]=-(ppl*ct-ppt1*sin(acos(ct)));
      PUP[1][4]=-((ppl*sin(acos(ct))+ppt1*ct)*cos(phi)-ppt2*sin(phi));
      PUP[2][4]=-((ppl*sin(acos(ct))+ppt1*ct)*sin(phi)+ppt2*cos(phi));

      if (emit){
      PUP[4][4]=ppe1;
      PUP[4][5] = x2*emcm/2.;
      PUP[4][6]=ppe2;
      PUP[3][6]=-(ppl2*ct+ppt1*sin(acos(ct)));
      PUP[1][6]=-((ppl2*sin(acos(ct))-ppt1*ct)*cos(phi)+ppt2*sin(phi));
      PUP[2][6]=-((ppl2*sin(acos(ct))-ppt1*ct)*sin(phi)-ppt2*cos(phi));
      }}
    if (!emit){
      PUP[3][6]=-PUP[3][4]-PUP[3][5];
      PUP[1][6]=-PUP[1][4]-PUP[1][5];
      PUP[2][6]=-PUP[2][4]-PUP[2][5];
      PUP[4][4] = x1*emcm/2.;
      PUP[4][5] = x2*emcm/2.;
      PUP[4][6] = x3*emcm/2.;
     
        double alpha = 1.;
      int nit = 5;
     //boost masses if needed. Need to boost to nominal gluon mass always.
     double Eq,Eqb,Eg;
     double pq = sqrt(pow(PUP[1][4],2)+pow(PUP[2][4],2)+pow(PUP[3][4],2));
     double pqb = sqrt(pow(PUP[1][5],2)+pow(PUP[2][5],2)+pow(PUP[3][5],2));
     double pg = sqrt(pow(PUP[1][6],2)+pow(PUP[2][6],2)+pow(PUP[3][6],2));
     
         for (int ix=1; ix < nit+1; ix++) {
       Eq = sqrt(pow(mt,2)+alpha*pow(pq,2));
       Eqb = sqrt(pow(mt,2)+alpha*pow(pqb,2));
       Eg =  sqrt(pow(mg,2)+alpha*pow(pg,2));
       alpha = alpha+(2.*(emcm-Eq-Eqb-Eg))/(pow(pq,2)/Eq+pow(pqb,2)/Eqb+pow(pg,2)/Eg);
         }
     PUP[1][4] =sqrt(alpha)*PUP[1][4];
     PUP[2][4] =sqrt(alpha)*PUP[2][4];
     PUP[3][4] =sqrt(alpha)*PUP[3][4];
     PUP[4][4] =Eq;
     PUP[1][5] =sqrt(alpha)*PUP[1][5];
     PUP[2][5] =sqrt(alpha)*PUP[2][5];
     PUP[3][5] =sqrt(alpha)*PUP[3][5];
     PUP[4][5] =Eqb;
     PUP[1][6] =sqrt(alpha)*PUP[1][6];
     PUP[2][6] =sqrt(alpha)*PUP[2][6];
     PUP[3][6] =sqrt(alpha)*PUP[3][6];
     PUP[4][6] =Eg;
     
      } else {double phi3=2.*pi*random(seed);
	      double ph=2.*pi*random(seed);
	      double phh=2.*pi*random(seed);
	      if (Jt==2) {PUP[4][4]=PUP[4][5];}   
	      double ptgg;
	      if (random(seed) < 0.5) {ptgg=sqrt(ptr2);}else{ptgg=-sqrt(ptr2);}
              double btgg=ptr2/(4.*(1.-z1)*pow(emcm-PUP[4][4],2));
              double pegg=(1.-z1)*(emcm-PUP[4][4])+btgg*(emcm-PUP[4][4]);
              double pzgg=(1.-z1)*(emcm-PUP[4][4])-btgg*(emcm-PUP[4][4]);
              double sth=sqrt(ptr2)/(z1*(emcm-PUP[4][4]));
	      double ptq, ptg;
	      if (random(seed) < 0.5) {
	      ptq=sqrt(ptsq)*cos(asin(sth))+z*z1*(emcm-PUP[4][4])*sth;
	      ptg=-sqrt(ptsq)*cos(asin(sth))+(1-z)*z1*(emcm-PUP[4][4])*sth;
              } else {
	      ptq=-sqrt(ptsq)*cos(asin(sth))+z*z1*(emcm-PUP[4][4])*sth;
              ptg=sqrt(ptsq)*cos(asin(sth))+(1.-z)*z1*(emcm-PUP[4][4])*sth;}

	      double btq=pow(ptq,2)/(4.*z1*z*pow(emcm-PUP[4][4],2));
	      double btg=pow(ptg,2)/(4.*(1.-z)*z1*pow(emcm-PUP[4][4],2));
	      double pe=z*z1*(emcm-PUP[4][4])+btq*(emcm-PUP[4][4]);
	      double pzq=z*z1*(emcm-PUP[4][4])-btq*(emcm-PUP[4][4]);
	      double peg=z1*(1.-z)*(emcm-PUP[4][4])+btg*(emcm-PUP[4][4]);
	      double pzg=z1*(1.-z)*(emcm-PUP[4][4])-btg*(emcm-PUP[4][4]);
	      PUP[4][6]=pegg;
	      double p6t1=ptgg*cos(ph);
	      double p6t2=ptgg*sin(ph);
	      PUP[1][6]=-((pzgg*sin(acos(ct))+p6t1*ct)*cos(phh)-p6t2*sin(phh));
	      PUP[2][6]=-((pzgg*sin(acos(ct))+p6t1*ct)*sin(phh)+p6t2*cos(phh));
	      PUP[3][6]=-(pzgg*ct-p6t1*sin(acos(ct)));
	      PUP[4][5]=pe;
	      double p5t1=ptq*cos(phig);
	      double p5t2=ptq*sin(phig);
	      PUP[1][5]=-((pzq*sin(acos(ct))+p5t1*ct)*cos(phi3)-p5t2*sin(phi3));
	      PUP[2][5]=-((pzq*sin(acos(ct))+p5t1*ct)*sin(phi3)+p5t2*cos(phi3));
	      PUP[3][5]=-(pzq*ct-p5t1*sin(acos(ct)));
              PUP[4][7]=peg;
	      double p7t1=ptg*cos(phig);
	      double p7t2=ptg*sin(phig);
	      PUP[1][7]=-((pzg*sin(acos(ct))+p7t1*ct)*cos(phi3)-p7t2*sin(phi3));
	      PUP[2][7]=-((pzg*sin(acos(ct))+p7t1*ct)*sin(phi3)+p7t2*cos(phi3));
	      PUP[3][7]=-(pzg*ct-p7t1*sin(acos(ct)));
	      PUP[1][4]=-PUP[1][7]-PUP[1][6]-PUP[1][5];
	      PUP[2][4]=-PUP[2][7]-PUP[2][6]-PUP[2][5];
	      PUP[3][4]=-PUP[3][7]-PUP[3][6]-PUP[3][5];
	      if (Jt==2) {
         PUP[1][8]= PUP[1][4]; PUP[2][8]= PUP[2][4]; PUP[3][8]= PUP[3][4];
	 PUP[1][4]= PUP[1][5]; PUP[2][4]= PUP[2][5]; PUP[3][4]= PUP[3][5];
	 PUP[1][5]= PUP[1][8]; PUP[2][5]= PUP[2][8]; PUP[3][5]= PUP[3][8];}
	 double pq = sqrt(pow(PUP[1][4],2)+pow(PUP[2][4],2)+pow(PUP[3][4],2));
	 double pqb = sqrt(pow(PUP[1][5],2)+pow(PUP[2][5],2)+pow(PUP[3][5],2));
	 double pg = sqrt(pow(PUP[1][6],2)+pow(PUP[2][6],2)+pow(PUP[3][6],2));     
	 double pgg = sqrt(pow(PUP[1][7],2)+pow(PUP[2][7],2)+pow(PUP[3][7],2));     
	 double Eq,Eqb,Eg,Egg;
	 double alpha=1;
	 double nit=5;
	 for (int ixx=1; ixx < nit+1; ixx++) {
       Eq = sqrt(pow(mt,2)+alpha*pow(pq,2));
       Eqb = sqrt(pow(mt,2)+alpha*pow(pqb,2));
       Eg =  sqrt(pow(mg,2)+alpha*pow(pg,2));
       Egg= sqrt(pow(mg,2)+alpha*pow(pgg,2));
       alpha = alpha+(2.*(emcm-Eq-Eqb-Eg-Egg))/(pow(pq,2)/Eq+pow(pqb,2)/Eqb+pow(pg,2)/Eg+pow(pgg,2)/Egg);
         }
     PUP[1][4] =sqrt(alpha)*PUP[1][4];
     PUP[2][4] =sqrt(alpha)*PUP[2][4];
     PUP[3][4] =sqrt(alpha)*PUP[3][4];
     PUP[4][4] =Eq;
     PUP[1][5] =sqrt(alpha)*PUP[1][5];
     PUP[2][5] =sqrt(alpha)*PUP[2][5];
     PUP[3][5] =sqrt(alpha)*PUP[3][5];
     PUP[4][5] =Eqb;
     PUP[1][6] =sqrt(alpha)*PUP[1][6];
     PUP[2][6] =sqrt(alpha)*PUP[2][6];
     PUP[3][6] =sqrt(alpha)*PUP[3][6];
     PUP[4][6] =Eg;
     PUP[1][7] =sqrt(alpha)*PUP[1][7];
     PUP[2][7] =sqrt(alpha)*PUP[2][7];
     PUP[3][7] =sqrt(alpha)*PUP[3][7];
     PUP[4][7] =Egg;
    }
     //Momenta of incoming particles and Z boson
   
     PUP[4][1]=emcm/2.;
     PUP[1][1]=0.;
     PUP[2][1]=0.;
     PUP[3][1]=emcm/2;
     PUP[4][2]=emcm/2;
     PUP[1][2]=0.;
     PUP[2][2]=0.;
     PUP[3][2]=-emcm/2;
     PUP[4][3]=emcm;
     PUP[1][3]=0.;
     PUP[2][3]=0.;
     PUP[3][3]=0.;

      double sum1, sum2, sum3, sum4;
     sum1=0.;sum2=0.;sum3=0.;sum4=0.;
     ii+=1;
     // get cosines and sines for production process
           T1[1]=0.;
	   T1[2]=0.;
	   T1[3]=sqrt(pow(PUP[1][4],2)+pow(PUP[2][4],2)+pow(PUP[3][4],2));
	   T2[1]=0.;
	   T2[2]=0.;
	   T2[3]=sqrt(pow(PUP[1][5],2)+pow(PUP[2][5],2)+pow(PUP[3][5],2));
	   costz=PUP[3][4]/T1[3];
	   costx=PUP[1][4]/sqrt(pow(PUP[1][4],2)+pow(PUP[2][4],2));
           sintz=sin(acos(costz));
           sintx=sin(acos(costx));
	   costbz=PUP[3][5]/T2[3];
	   costbx=PUP[1][5]/sqrt(pow(PUP[1][5],2)+pow(PUP[2][5],2));
           sintbz=sin(acos(costbz));
           sintbx=sin(acos(costbx));
    }

      
    }}

    /*************************************END OF PRODUCTION*************************************/

    /*********************************START DECAY PROCESS*************************************/

      if (POWdecay) {

	int Jtd, Jtbd; // Identifier for solutions
	Jtd =0; Jtbd =0;
      bool diterate =false;
      double intg=0.;
      int ev = 0; // ev =1, 2 for top and tbar decay respectively
      for (int ia=0; diterate==false || ev != 2; ia++){
      double logn=fabs(log(random(seed)));
      intg+=logn;
      double pt2=interpd(1000,intg,dend1,dend2,dend3,dend4)*pow(mt,2);
      if (0.*intg !=0) {break;}
      if (intg > dend4) {break;}
      if (pt2 < pt2mind) {break;}
      double xdmax=(mt-a*mt-2.*sqrt(a*pt2))/mt; // maximum value of x
      double xdmin=2.*sqrt(pt2)/mt; // minimum value of x
       double Mdmax=0.;
       double x3, x11, x12, Md, Md1, Md2;
      for (int ix=0; ix < 1000 ; ix++) {
       x3=(xdmax-xdmin)*random(seed)+xdmin;
       Md= MEd (pt2, x3, x11, x12, Md1, Md2);
         if (fabs(Md) > fabs(Mdmax)) {Mdmax=Md;}
       }
      do {
       x3=(xdmax-xdmin)*random(seed)+xdmin;
       Md= MEd (pt2, x3, x11, x12, Md1, Md2);}while (random(seed) > fabs(Md/Mdmax));
   

    if (random(seed) < fabs(Md1/Md)) {
    	   x1=x11;
	    if (ev==0){ xg1 = x3; xw11=x1; pt2r1=pt2/pow(mt,2); Jtd = 1;}
	    if (ev==1){ xg2 = x3; xw12=x1; pt2r2=pt2/pow(mt,2); Jtbd =1;}
	    ev++;
	    diterate =true;
       }  else {
	 x1=x12;
         if (ev==0){xg1 = x3; xw11=x1; pt2r1=pt2/pow(mt,2); Jtd = 2;}
	 if (ev==1){ xg2 = x3; xw12=x1; pt2r2=pt2/pow(mt,2); Jtbd = 2;}
	 ev++;
	diterate =true;
       }
      }
      
      /******************************START DECAY TRUNCATED SHOWER*****************************************/   
	double z1, zt1, ztb1, zt2, ztb2, q;
	zt1 =0.; ztb1 =0.;zt2 =0.; ztb2 =0.;
	double qht2,qhtb2,qi, kapp;
	qht2 =0.; qhtb2 =0.;
	//	emitt =false;
	demitt=false;
        demittb=false;

      	//Decay truncated shower
      if (dtruncate) {
        double Qg=0.75;
	double u=mt;
	double ub=0.75; 
	double qc=sqrt(1.53*pow(mt,2));
	double qcc=2.8;
	double CF =4./3.;
	double alpham=1.;
	double mb = 0.; 
	double c = pow(mb/mt,2);

	
	if (Jtd ==1 ) {
	  zt1= -0.5*(2.*xg1-2.+xw11+a-sqrt(pow(xw11+a,2)-4.*a))/sqrt(pow(xw11+a,2)-4.*a);
	  zt2=-0.5*(2.*pow((xw11+a),2)+4*a*xg1-8*a-pow((xw11+a),2)*xg1-sqrt(-48.*(xw11+a)*a*xg1-16.*a-32.*a*c*xg1+8.*pow((xw11+a),2)*c*xg1+16.*(xw11+a)*a*c*xg1+32.*pow(a,2)*c-16.*a*pow(c,2)-8.*pow((xw11+a),2)*a-8.*pow((xw11+a),2)*c+32.*(xw11+a)*pow(a,2)-4.*a*pow(xg1,2)*pow((xw11+a),2)-32.*(xw11+a)*a*c-8.*pow((xw11+a),2)*a*c+16.*(xw11+a)*pow(xg1,2)*a+4.*pow((xw11+a),3)*a*xg1-4.*pow((xw11+a),3)*c*xg1+8.*pow((xw11+a),2)*a*xg1-16.*(xw11+a)*pow(a,2)*xg1-32.*pow(a,2)+4.*pow((xw11+a),2)+4.*pow((xw11+a),4)-16.*pow(a,3)-8.*pow((xw11+a),3)-8.*pow((xw11+a),2)*xg1+pow((xw11+a),4)*pow(xg1,2)+32.*a*xg1+32.*a*c+32.*(xw11+a)*a+32.*pow(a,2)*xg1-16.*a*pow(xg1,2)-8.*pow((xw11+a),3)*a+8.*pow((xw11+a),3)*c-4.*pow((xw11+a),4)*xg1-4.*pow((xw11+a),3)*pow(xg1,2)+12.*pow((xw11+a),3)*xg1+4.*pow((xw11+a),2)*pow(a,2)+4.*pow((xw11+a),2)*pow(c,2)+4.*pow((xw11+a),2)*pow(xg1,2)))/(4.*a-pow((xw11+a),2));}
	if (Jtd ==2){zt2=-0.5*(2.*pow((xw11+a),2)+4*a*xg1-8*a-pow((xw11+a),2)*xg1-sqrt(-48.*(xw11+a)*a*xg1-16.*a-32.*a*c*xg1+8.*pow((xw11+a),2)*c*xg1+16.*(xw11+a)*a*c*xg1+32.*pow(a,2)*c-16.*a*pow(c,2)-8.*pow((xw11+a),2)*a-8.*pow((xw11+a),2)*c+32.*(xw11+a)*pow(a,2)-4.*a*pow(xg1,2)*pow((xw11+a),2)-32.*(xw11+a)*a*c-8.*pow((xw11+a),2)*a*c+16.*(xw11+a)*pow(xg1,2)*a+4.*pow((xw11+a),3)*a*xg1-4.*pow((xw11+a),3)*c*xg1+8.*pow((xw11+a),2)*a*xg1-16.*(xw11+a)*pow(a,2)*xg1-32.*pow(a,2)+4.*pow((xw11+a),2)+4.*pow((xw11+a),4)-16.*pow(a,3)-8.*pow((xw11+a),3)-8.*pow((xw11+a),2)*xg1+pow((xw11+a),4)*pow(xg1,2)+32.*a*xg1+32.*a*c+32.*(xw11+a)*a+32.*pow(a,2)*xg1-16.*a*pow(xg1,2)-8.*pow((xw11+a),3)*a+8.*pow((xw11+a),3)*c-4.*pow((xw11+a),4)*xg1-4.*pow((xw11+a),3)*pow(xg1,2)+12.*pow((xw11+a),3)*xg1+4.*pow((xw11+a),2)*pow(a,2)+4.*pow((xw11+a),2)*pow(c,2)+4.*pow((xw11+a),2)*pow(xg1,2)))/(4.*a-pow((xw11+a),2));
}
	if (Jtbd==1){
	  ztb1= -0.5*(2.*xg2-2.+xw12+a-sqrt(pow(xw12+a,2)-4.*a))/sqrt(pow(xw12+a,2)-4.*a);
	  ztb2=-0.5*(2.*pow((xw12+a),2)+4*a*xg2-8*a-pow((xw12+a),2)*xg2-sqrt(-48.*(xw12+a)*a*xg2-16.*a-32.*a*c*xg2+8.*pow((xw12+a),2)*c*xg2+16.*(xw12+a)*a*c*xg2+32.*pow(a,2)*c-16.*a*pow(c,2)-8.*pow((xw12+a),2)*a-8.*pow((xw12+a),2)*c+32.*(xw12+a)*pow(a,2)-4.*a*pow(xg2,2)*pow((xw12+a),2)-32.*(xw12+a)*a*c-8.*pow((xw12+a),2)*a*c+16.*(xw12+a)*pow(xg2,2)*a+4.*pow((xw12+a),3)*a*xg2-4.*pow((xw12+a),3)*c*xg2+8.*pow((xw12+a),2)*a*xg2-16.*(xw12+a)*pow(a,2)*xg2-32.*pow(a,2)+4.*pow((xw12+a),2)+4.*pow((xw12+a),4)-16.*pow(a,3)-8.*pow((xw12+a),3)-8.*pow((xw12+a),2)*xg2+pow((xw12+a),4)*pow(xg2,2)+32.*a*xg2+32.*a*c+32.*(xw12+a)*a+32.*pow(a,2)*xg2-16.*a*pow(xg2,2)-8.*pow((xw12+a),3)*a+8.*pow((xw12+a),3)*c-4.*pow((xw12+a),4)*xg2-4.*pow((xw12+a),3)*pow(xg2,2)+12.*pow((xw12+a),3)*xg2+4.*pow((xw12+a),2)*pow(a,2)+4.*pow((xw12+a),2)*pow(c,2)+4.*pow((xw12+a),2)*pow(xg2,2)))/(4.*a-pow((xw12+a),2));}
	if (Jtbd == 2){
	  ztb2=-0.5*(2.*pow((xw12+a),2)+4*a*xg2-8*a-pow((xw12+a),2)*xg2-sqrt(-48.*(xw12+a)*a*xg2-16.*a-32.*a*c*xg2+8.*pow((xw12+a),2)*c*xg2+16.*(xw12+a)*a*c*xg2+32.*pow(a,2)*c-16.*a*pow(c,2)-8.*pow((xw12+a),2)*a-8.*pow((xw12+a),2)*c+32.*(xw12+a)*pow(a,2)-4.*a*pow(xg2,2)*pow((xw12+a),2)-32.*(xw12+a)*a*c-8.*pow((xw12+a),2)*a*c+16.*(xw12+a)*pow(xg2,2)*a+4.*pow((xw12+a),3)*a*xg2-4.*pow((xw12+a),3)*c*xg2+8.*pow((xw12+a),2)*a*xg2-16.*(xw12+a)*pow(a,2)*xg2-32.*pow(a,2)+4.*pow((xw12+a),2)+4.*pow((xw12+a),4)-16.*pow(a,3)-8.*pow((xw12+a),3)-8.*pow((xw12+a),2)*xg2+pow((xw12+a),4)*pow(xg2,2)+32.*a*xg2+32.*a*c+32.*(xw12+a)*a+32.*pow(a,2)*xg2-16.*a*pow(xg2,2)-8.*pow((xw12+a),3)*a+8.*pow((xw12+a),3)*c-4.*pow((xw12+a),4)*xg2-4.*pow((xw12+a),3)*pow(xg2,2)+12.*pow((xw12+a),3)*xg2+4.*pow((xw12+a),2)*pow(a,2)+4.*pow((xw12+a),2)*pow(c,2)+4.*pow((xw12+a),2)*pow(xg2,2)))/(4.*a-pow((xw12+a),2));
      }
	
         qi=mt;
	 double ptsqt,ptsqtb;
	 ptsqt = 0.; ptsqtb =0.;
	 if (Jtd == 2) {
	   ptsqt = pow(mt,2)*pt2r1;
	   qht2= sqrt(xg1/(1.-zt2)*pow(mt,2));}
	 if (Jtbd == 2) {
	   ptsqtb = pow(mt,2)*pt2r2;
	   qhtb2=sqrt(xg2/(1.-ztb2)*pow(mt,2));} 
	 if (Jtd == 1) {
	   qht2= sqrt(1.53*pow(mt,2));
	   ptsqt =sqrt(1.53*pow(mt,2));}
	   if(Jtbd == 1) {
	   qhtb2= sqrt(1.53*pow(mt,2));
	   ptsqtb =sqrt(1.53*pow(mt,2));}
	   double zmin=1.-(1.-a)/(qi/mt+2.*sqrt(a*(qi/mt-1)));
	  double zmax = 1.-Qg/qi;
         double C=CF*alpham*log((1-zmax)/(1.-zmin))/pi;
	  double D=exp(-2.*log(qht2/qc)*C);
	  double E=exp(-2.*log(qi/qc)*C);
	if (random(seed) < 1.-E/D) {
	 double Pz, Pzmax; 
	 do {
	    
	  q=sqrt(exp(log(E/random(seed))/C)*pow(qc,2));
	 Pzmax=(1.+pow((1.-Qg/mt),2))/(Qg/mt)+2.*pow(u,2)/(Qg/mt*(1.-Qg/mt)*pow(q,2));
	  z1=random(seed)*(+Qg/q-(1.-a)/(pow(q/mt,2)+2.*sqrt(a*(pow(q/mt,2)-1))))+1.-Qg/q;
	   kapp=pow(q/mt,2);
	   ptr2 = pow(mt,2)*(kapp-1.)*pow(1.-z1,2);
	   Pz=(1.+pow(z1,2))/(1.-z1)+2.*pow(u,2)/(z1*(1.-z1)*pow(q,2));} while (random(seed) > Pz/Pzmax || (kapp > 1.+a*pow(1.-sqrt(z1*(1.-a)/(a*(1.-z1))),2)) || kapp < 1 || q > qht2 || q < qi) ;
	 demitt = true;
	    
	  //energy fractions
	   xgrt=(1.-z1)*kapp;
	   double uu=1.+a-c-xgrt;
	   double wr=1.-(1.-z1)*(kapp-1.);
	   double vv=sqrt(pow(uu,2)-4.*a*wr*z1);
	   xwrt=(uu+vv)/(2.*wr)+(uu-vv)/(2.*z1)-a;
	   double alphas = 1./(b*log(ptr2/(pt2min*pow(mt,2))))*(1.-bp*log(log(ptr2/(pt2min*pow(mt,2))))/(b*log(ptr2/(pt2min*pow(mt,2)))));
	  alphas=alphas+pow(alphas,2)*(67./3.-pow(pi,2)-50./9.)/(4.*pi);
	 if (ptr2 < 0 || q < qi || ptr2 > ptsqt || random(seed) > alphas/alpham || random(seed) > Pz*(1.-z1)/2.|| xgrt*mt/2. < sqrt(ptr2)) {demitt=false;}
	 }
	 
	  C=CF*alpham*log((1-zmax)/(1.-zmin))/pi;
	  D=exp(-2.*log(qhtb2/qc)*C);
	  E=exp(-2.*log(qi/qc)*C);
	 if (random(seed) < 1.-E/D) {
	      double Pz, Pzmax;
	  do {
	  q=sqrt(exp(log(E/random(seed))/C)*pow(qc,2));
	  Pzmax=(1.+pow((1.-Qg/mt),2))/(Qg/mt)+2.*pow(u,2)/(Qg/mt*(1.-Qg/mt)*pow(q,2));
	    z1=random(seed)*(+Qg/q-(1.-a)/(pow(q/mt,2)+2.*sqrt(a*(pow(q/mt,2)-1))))+1.-Qg/q;
	     kapp=pow(q/mt,2);
	    ptr2b = pow(mt,2)*(kapp-1.)*pow(1.-z1,2);
	   Pz=(1.+pow(z1,2))/(1.-z1)+2.*pow(u,2)/(z1*(1.-z1)*pow(q,2));} while (random(seed) > Pz/Pzmax || (kapp > 1.+a*pow(1.-sqrt(z1*(1.-a)/(a*(1.-z1))),2)) || kapp < 1 || q > qhtb2 || q < qi);
	  demittb = true;
	  
	  //energy fractions
	   xgrtb=(1.-z1)*kapp;
	   double uu=1.+a-c-xgrtb;
	   double wr=1.-(1.-z1)*(kapp-1.);
	   double vv=sqrt(pow(uu,2)-4.*a*wr*z1);
	 
	   xwrtb=(uu+vv)/(2.*wr)+(uu-vv)/(2.*z1)-a;

	  double alphas = 1./(b*log(ptr2b/(pt2min*pow(mt,2))))*(1.-bp*log(log(ptr2b/(pt2min*pow(mt,2))))/(b*log(ptr2b/(pt2min*pow(mt,2)))));
	  alphas=alphas+pow(alphas,2)*(67./3.-pow(pi,2)-50./9.)/(4.*pi);
	 
	  if (ptr2b < 0 || q < qi || ptr2b > ptsqtb || random(seed) > alphas/alpham || random(seed) > Pz*(1.-z1)/2.|| xgrtb*mt/2. < sqrt(ptr2b)) {demittb=false;}
	 
	 }
	
	 if (!demitt && Jtd == 1) {
	   qi=sqrt(1.16*pow(mt,2));
	   ptsqt = pow(mt,2)*pt2r1;
	   qht2=sqrt(pow(mt,2)*(xw11-1.)/(-zt1*(1.-zt1)));
	   C=CF*alpham*log(Qg/(qi-ub))/pi;
	   D=exp(2.*log(qht2/qcc)*C);
	   E=exp(2.*log(qi/qcc)*C);
	   
	 if (random(seed) < 1.-E/D) {
	   double Pz, Pzmax;
	 
	  do {
	    
	  q=sqrt(exp(log(E/random(seed))/C)*pow(qcc,2));
	  Pzmax=(1.+pow((1.-Qg/qi),2))/(Qg/qi);
	  
	  
	    z1=random(seed)*(1.-Qg/qi-ub/qi)+ub/qi;
	    Pz=(1.+pow(z1,2))/(1.-z1);} while (random(seed) > Pz/Pzmax || q < qht2 || q >qi);
         
	     
	   demitt = true;
	   kapp=pow(q/mt,2);
	   ptr2 = pow(mt,2)*kapp*pow(z1*(1.-z1),2);
	  //energy fractions
	   xwrt=1.-c-z1*(1.-z1)*kapp;
	  
	   xgrt=0.5*(2.-xwrt-a)-(z1-0.5)*sqrt(pow(xwrt+a,2)-4.*a);

	  double alphas = 1./(b*log(ptr2/(pt2min*pow(mt,2))))*(1.-bp*log(log(ptr2/(pt2min*pow(mt,2))))/(b*log(ptr2/(pt2min*pow(mt,2)))));
	  alphas=alphas+pow(alphas,2)*(67./3.-pow(pi,2)-50./9.)/(4.*pi);
	  if (ptr2 < 0 || q > qi || ptr2 > ptsqt || random(seed) > alphas/alpham || random(seed) > Pz*(1.-z1)/2.|| xgrt*mt/2. < sqrt(ptr2)) {demitt=false;}
	 	
	 }	   
	 } 
	 if (!demittb && Jtbd == 1) {
	   qi=sqrt(1.16*pow(mt,2));
	   ptsqtb = pow(mt,2)*pt2r2;
	   qhtb2=sqrt(pow(mt,2)*(xw12-1.)/(-ztb1*(1.-ztb1)));
	   C=CF*alpham*log(Qg/(qi-ub))/pi;
	   D=exp(2.*log(qhtb2/qcc)*C);
	   E=exp(2.*log(qi/qcc)*C);
	 if (random(seed) < 1.-E/D) {
	   double Pz, Pzmax;
	  do {
	  q=sqrt(exp(log(E/random(seed))/C)*pow(qcc,2));
	  Pzmax=(1.+pow((1.-Qg/qi),2))/(Qg/qi);
	 
	    z1=random(seed)*(1.-Qg/qi-ub/qi)+ub/qi;
	    Pz=(1.+pow(z1,2))/(1.-z1);} while (random(seed) > Pz/Pzmax|| q < qhtb2 || q >qi);
	  
	   demittb = true;
	   kapp=pow(q/mt,2);
	    ptr2b = pow(mt,2)*kapp*pow(z1*(1.-z1),2);
	  //energy fractions
           xwrtb=1-c-z1*(1.-z1)*kapp;
	   xgrtb=0.5*(2.-xwrtb-a)-(z1-0.5)*sqrt(pow(xwrtb+a,2)-4.*a);
double alphas = 1./(b*log(ptr2b/(pt2min*pow(mt,2))))*(1.-bp*log(log(ptr2b/(pt2min*pow(mt,2))))/(b*log(ptr2b/(pt2min*pow(mt,2)))));
	  alphas=alphas+pow(alphas,2)*(67./3.-pow(pi,2)-50./9.)/(4.*pi);
	 if (ptr2b < 0 || q > qi || ptr2b > ptsqtb || random(seed) > alphas/alpham || random(seed) > Pz*(1.-z1)/2. || xgrtb*mt/2. < sqrt(ptr2b)) {demittb=false;}
	 
	  }	   
	 }

	 
      }
   /**********************END DECAY TRUNCATED SHOWER****************************/
    }
      // Polarizer 
      if (!POWprod) {
   double qelr=qel;
   double qtlr1=qtl; 
   double qtlr2=qtr;
   double MM11=1.22446;double MM22=1.22446;double MM21=9.1491; double MM12=2.1955;
   if (Pm==1 && Pp==-1)
  {qelr=qer;
   qtlr1=qtr; 
   qtlr2=qtl;
   MM11=0.88966; MM22=0.88966; MM21=1.6669; MM12=11.4006;
}
   double RA=-Qt+qelr*qtlr1*Ref*16.*pow(costhw,2);
   double RB=qelr*qtlr1*Imf*16.*pow(costhw,2);
   double LA=-Qt+qelr*qtlr2*Ref*16.*pow(costhw,2);
   double LB=qelr*qtlr2*Imf*16.*pow(costhw,2);
   double M11max = beta*(pow(RA+LA,2)+pow(RB+LB,2))*(1.-pow(beta,2));
   double M22max=beta*(pow(RA+LA,2)+pow(RB+LB,2))*(1.-pow(beta,2));
   double M21max = beta*(pow((RA*(1.+beta)+LA*(1.-beta)),2)+pow((RB*(1.+beta)+LB*(1.-beta)),2))*4.;
   double M12max = beta*(pow((RA*(1.-beta)+LA*(1.+beta)),2)+pow((RB*(1.-beta)+LB*(1.+beta)),2))*4.;
   if (Pm==1 && Pp==-1) {double Mdum=M21max; M21max=M12max; M12max=Mdum;}
	double MM, MMmax;
	double MMM=fabs(MM11)+fabs(MM12)+fabs(MM21)+fabs(MM22);
	if (random(seed) < (fabs(MM11)+fabs(MM12))/MMM){
	  if(random(seed) < fabs(MM11)/(fabs(MM11)+fabs(MM12))){
	    EVT=1; MM=MM11; MMmax=M11max; }else {
	      EVT=3; MM=MM12; MMmax=M12max; }}else{
           if(random(seed) < fabs(MM21)/(fabs(MM21)+fabs(MM22))){
	     EVT =2;MM=MM21; MMmax=M21max;}else {
	       EVT=4; MM=MM22; MMmax=M22max; }}

	double Mm;
	do {
	 cost = 1.-2.*random(seed);
	  if (EVT==1) {
	    Mm = beta*(pow(RA+LA,2)+pow(RB+LB,2))*(1.-pow(beta,2))*(1.-pow(cost,2));}
	  if (EVT==3) {
	    Mm = beta*(pow((RA*(1.+beta)+LA*(1.-beta)),2)+pow((RB*(1.+beta)+LB*(1.-beta)),2))*pow(1.+cost,2);}
	  if(EVT==2) {
	    Mm =  beta*(pow((RA*(1.-beta)+LA*(1.+beta)),2)+pow((RB*(1.-beta)+LB*(1.+beta)),2))*pow(1.-cost,2);}
	  if(EVT==4) {
	Mm = beta*(pow(RA+LA,2)+pow(RB+LB,2))*(1.-pow(beta,2))*(1.-pow(cost,2));}
	if (Pm==-1 && Pp==1) {
	  if(EVT==2) {
	    Mm = beta*(pow((RA*(1.+beta)+LA*(1.-beta)),2)+pow((RB*(1.+beta)+LB*(1.-beta)),2))*pow(1.+cost,2);}
	  if (EVT==3) {
	    Mm =  beta*(pow((RA*(1.-beta)+LA*(1.+beta)),2)+pow((RB*(1.-beta)+LB*(1.+beta)),2))*pow(1.-cost,2);}
	}} while (random(seed) > Mm/MMmax);
		T1[1]=0.5*emcm*beta*sin(acos(cost));
		T1[2]=0.;
	        T1[3]=0.5*emcm*beta*cost;
	        T2[1]=-0.5*emcm*beta*sin(acos(cost));
	        T2[2]=0.;
		T2[3]=-0.5*emcm*beta*cost;
			sint=sin(acos(cost));
}
      /**********************END DECAY PROCESS ****************************/ 

      /***********************START DISTRIBUTING BORN EVENTS ACCORDING TO BBAR FUNCTION AND ASSIGNING MMT TO DECAY PRODUCTS********************/
     
      double vt=sqrt(pow(T1[3],2)+pow(T1[1],2)+pow(T1[2],2))/sqrt(pow(T1[3],2)+pow(T1[1],2)+pow(T1[2],2)+pow(mt,2));
      double gammat=1./sqrt(1.-vt*vt);
      double vtb=sqrt(pow(T2[3],2)+pow(T2[1],2)+pow(T2[2],2))/sqrt(pow(T2[3],2)+pow(T2[1],2)+pow(T2[2],2)+pow(mt,2));
      double gammatb=1./sqrt(1.-vtb*vtb);
     double Ew,Eg,Eb,El,Ew2,Eg2,Eb2,xw2,xw1,xb1,xb2,bgphi,pw;
      int tm, tbm, st,stb;
      tm = 0; tbm = 0;
        if (EVT==1) {tm=1; tbm=-1; st=1;stb=-1;}
 	if (EVT==2) {tm=-1;tbm=-1; st=-1;stb=-1;}
 	if (EVT==3) {tm=1;tbm=1; st=1; stb=1;}
 	if (EVT==4){tm=-1;tbm=1; st=-1; stb=1;}
	if (POWdecay) {
           xw1=xw11+a;
	   xw2=xw12+a;
	   xb2=2.-xw2-xg2;
	   xb1= 2.-xw1-xg1;
	   Ew=xw1*mt/2.;
           Eg=xg1*mt/2.;
           Eb=xb1*mt/2.;
        if (POWdecay && demitt) {
           Ew = (2.-xb1-xg1-xgrt)*mt/2.; 
	     if (Ew < mw) {demitt = false; Ew=xw1*mt/2.;}
}
	   
           Ew2=xw2*mt/2.;
           Eg2=xg2*mt/2.;
           Eb2=xb2*mt/2.;
	   if (POWdecay  && demittb) {
	     Ew2 = (2.-xb2-xg2-xgrtb)*mt/2.;
	    if (Ew2 < mw) {demittb = false; Ew2=xw2*mt/2.;}
}

           El=0.5*mw;


} else {Ew = 0.5*(mt+mw*mw/mt);Eb=0.5*(mt-mw*mw/mt); Eg=0; El=0.5*mw;Ew2=Ew; Eb2=Eb; Eg2=Eg;}
	 double  xww1, xww2;
	 xww1 = xw1 ;
	 xww2 = xw2 ;
        pw=sqrt(pow(Ew,2)-pow(mw,2));

	double ctw,ctl,max,bphi,TL,TO,TI,ctw2,ctl2,bphi2,TL2,TO2,TI2,bgphi2;
	// Production process only with LO decays
	if (!POWdecay) {
        for(int ix=0;;ix++) {
        ctw=1.-2.*random(seed);
 	ctl=1.-2.*random(seed);
 	bphi=2*pi*random(seed);
	ctw2=1.-2.*random(seed);
 	ctl2=1.-2.*random(seed);
 	bphi2=2*pi*random(seed);
	TL=pow(mt,2)*(1.+tm*ctw)*(1.-pow(ctl,2));
	TO=pow(mw,2)*(1.-tm*ctw)*pow((1.-ctl),2);
	TI=tm*2.*mt*mw*(1.-ctl)*sqrt(1.-pow(ctw,2))*sqrt(1.-pow(ctl,2))*cos(bphi);
        TL2=pow(mt,2)*(1.+tbm*ctw2)*(1.-pow(ctl2,2));
        TO2=pow(mw,2)*(1.-tbm*ctw2)*pow((1.-ctl2),2);
        TI2=tbm*2.*mt*mw*(1.-ctl2)*sqrt(1.-pow(ctw2,2))*sqrt(1.-pow(ctl2,2))*cos(bphi2); 
	max=pow(pow(mt,2)*2.+pow(mw,2)*8.+2.*mt*mw*2.,2);
 	if (random(seed)< (TL+TO+TI)*(TL2+TO2+TI2)/(max)) {break;}}} else 
	  // NLO decay process hence Bbar function
	  if (POWdecay)
	    
       {double xd0=0.2/(mt/2.*(1.-a));
	double zd0=0.0028;
	 double cbgst1, cbgst2;
	if (!findmaxdec) {
	  for (int idec=1;idec <1000; idec++) {
	    do {
	xg1=random(seed)*(1.-a);
	xw1=a+random(seed)*(1.-(a*xg1/(1-xg1)+(1.-xg1)))+a*xg1/(1-xg1)+(1.-xg1);
	xb1=2.-xw1-xg1;
	double bw11=sqrt(1.-4.*a/pow(xw1,2));
	double wb=acos(1./(xw1*bw11*xb1)*(xg1-xw1-xb1+xw1*xb1+2.*a));
        double wg=acos(1./(xw1*bw11*xg1)*(xb1-xw1-xg1+xw1*xg1+2.*a));
      	double cbg=cos(2.*pi*-wb-wg);
	double sbg=sin(2.*pi*-wb-wg);
	double xtp=sqrt(pow((2.-xg1)/2.*mt,2)-pow(xg1*0.5*mt,2));
	double xte=(2.-xg1)/2.*mt;
	double vtep=xtp/xte;
	double gamt=sqrt(1./(1.-pow(vtep,2)));
	double pb =gamt*(xb1*mt*0.5*cbg-vtep*xb1*mt*0.5*cbg);
	double pt = xb1*mt*0.5*sbg;
        cbgst1=cos(atan2(pt,pb));
	xg2=random(seed)*(1.-a);
	xw2=a+random(seed)*(1.-(a*xg2/(1-xg2)+(1.-xg2)))+a*xg2/(1-xg2)+(1.-xg2);
	xb2=2.-xw2-xg2;	
	double bw22=sqrt(1.-4.*a/pow(xw2,2));		  
	wb=acos(1./(xw2*bw22*xb2)*(xg2-xw2-xb2+xw2*xb2+2.*a));
        wg=acos(1./(xw2*bw22*xg2)*(xg2-xw2-xg2+xw2*xg2+2.*a));
      	cbg=cos(2.*pi*-wb-wg);
	sbg=sin(2.*pi*-wb-wg);
        xtp=sqrt(pow((2.-xg2)/2.*mt,2)-pow(xg2*0.5*mt,2));
        xte=(2.-xg2)/2.*mt;
        vtep=xtp/xte;
        gamt=sqrt(1./(1.-pow(vtep,2)));
	pb =gamt*(xb2*mt*0.5*cbg-vtep*xb2*mt*0.5*cbg);
        pt = xb2*mt*0.5*sbg;
	cbgst2=cos(atan2(pt,pb));
}
	while (xg2 < xd0 || xg1 < xd0 || (1.-cbgst1)/2. < zd0|| (1.-cbgst2)/2. < zd0);
	  ctw=1.-2.*random(seed);
 	  ctl=1.-2.*random(seed);
 	  bphi=2*pi*random(seed);
	  ctw2=1.-2.*random(seed);
 	  ctl2=1.-2.*random(seed);
 	  bphi2=2*pi*random(seed);
	  bgphi=2*pi*random(seed);
          bgphi2=2*pi*random(seed);
	  double T1=PolMEdec(xw1,xb1,xg1,ctw,ctl,bphi,bgphi,1,xd0,zd0);
          double T2=PolMEdec(xw1,xb1,xg1,ctw,ctl,bphi,bgphi,-1,xd0,zd0);
	  double T3=PolMEdec(xw2,xb2,xg2,ctw2,ctl2,bphi2,bgphi2,1,xd0,zd0);
          double T4=PolMEdec(xw2,xb2,xg2,ctw2,ctl2,bphi2,bgphi2,-1,xd0,zd0);
	  if (fabs(T1*T4) > maxd1) {maxd1=T1*T4;}
	  if (fabs(T2*T4) > maxd2) {maxd2=T2*T4;}
	  if (fabs(T1*T3) > maxd3) {maxd3=T1*T3;}
	  if (fabs(T2*T3) > maxd4) {maxd4=T2*T3;}
	  }
	}
	findmaxdec=true;
	double maxd,TTB;
	maxd =0.;
	if (EVT==1) {maxd=maxd1;}
	if (EVT==2) {maxd=maxd2;}
	if (EVT==3) {maxd=maxd3;}
	if (EVT==4) {maxd=maxd4;}
	do{
	  do {
	xg1=random(seed)*(1.-a);
	xw1=a+random(seed)*(1.-(a*xg1/(1-xg1)+(1.-xg1)))+a*xg1/(1-xg1)+(1.-xg1);
	xb1=2.-xw1-xg1;
	double bw11=sqrt(1.-4.*a/pow(xw1,2));
	double wb=acos(1./(xw1*bw11*xb1)*(xg1-xw1-xb1+xw1*xb1+2.*a));
        double wg=acos(1./(xw1*bw11*xg1)*(xb1-xw1-xg1+xw1*xg1+2.*a));
      	double cbg=cos(2.*pi*-wb-wg);
	double sbg=sin(2.*pi*-wb-wg);
	double xtp=sqrt(pow((2.-xg1)/2.*mt,2)-pow(xg1*0.5*mt,2));
	double xte=(2.-xg1)/2.*mt;
	double vtep=xtp/xte;
	double gamt=sqrt(1./(1.-pow(vtep,2)));
	double pb =gamt*(xb1*mt*0.5*cbg-vtep*xb1*mt*0.5*cbg);
	double pt = xb1*mt*0.5*sbg;
        cbgst1=cos(atan2(pt,pb));
	xg2=random(seed)*(1.-a);
	xw2=a+random(seed)*(1.-(a*xg2/(1-xg2)+(1.-xg2)))+a*xg2/(1-xg2)+(1.-xg2);
	xb2=2.-xw2-xg2;	
	double bw22=sqrt(1.-4.*a/pow(xw2,2));		  
	wb=acos(1./(xw2*bw22*xb2)*(xg2-xw2-xb2+xw2*xb2+2.*a));
        wg=acos(1./(xw2*bw22*xg2)*(xg2-xw2-xg2+xw2*xg2+2.*a));
      	cbg=cos(2.*pi*-wb-wg);
	sbg=sin(2.*pi*-wb-wg);
        xtp=sqrt(pow((2.-xg2)/2.*mt,2)-pow(xg2*0.5*mt,2));
        xte=(2.-xg2)/2.*mt;
        vtep=xtp/xte;
        gamt=sqrt(1./(1.-pow(vtep,2)));
	pb =gamt*(xb2*mt*0.5*cbg-vtep*xb2*mt*0.5*cbg);
        pt = xb2*mt*0.5*sbg;
	cbgst2=cos(atan2(pt,pb));
}while (xg2 < xd0 || xg1 < xd0 || (1.-cbgst1)/2. < zd0|| (1.-cbgst2)/2. < zd0);
	  ctw=1.-2.*random(seed);
 	  ctl=1.-2.*random(seed);
 	  bphi=2*pi*random(seed);
	  ctw2=1.-2.*random(seed);
 	  ctl2=1.-2.*random(seed);
 	  bphi2=2.*pi*random(seed);
	  bgphi=2.*pi*random(seed);
          bgphi2=2.*pi*random(seed);
	  TTB=PolMEdec(xw1,xb1,xg1,ctw,ctl,bphi,bgphi,tm,xd0,zd0)*PolMEdec(xw2,xb2,xg2,ctw2,ctl2,bphi2,bgphi2,tbm,xd0,zd0);}
	   while (random(seed) > fabs(TTB/maxd));
       }
	// Assign momenta
	double W1[4],E[4],EE1[4],NN1[4],BB1[4],G1[4],B1[4];
       double W2[4],EE2[4],NN2[4],BB2[4],G2[4],B2[4];	
           double phi=2.*pi*random(seed);
	// top rest frame
 	  W1[1]=pw*sin(acos(ctw))*cos(phi);
 	  W1[2]=pw*sin(acos(ctw))*sin(phi);
	  W1[3]=pw*(ctw);
	//W rest frame
 	E[1]=El*sin(acos(ctl))*cos(bphi);
 	E[2]=El*sin(acos(ctl))*sin(bphi);
	E[3]=El*ctl;
	double vv=pw/Ew;
	double gammaw=1./sqrt(1.-pow(vv,2));
        //top rest frame 
	E[3]=gammaw*(E[3]+vv*El);
	
	  EE1[3]=E[3]*ctw+E[1]*sin(acos(ctw));
	  EE1[1]=E[3]*sin(acos(ctw))-E[1]*ctw;
      	  EE1[0]=EE1[1];
	  EE1[1]=EE1[0]*cos(phi)-E[2]*sin(phi);
	  EE1[2]=E[2]*cos(phi)+EE1[0]*sin(phi);
	  NN1[1]=W1[1]-EE1[1];
	  NN1[2]=W1[2]-EE1[2];
	  NN1[3]=W1[3]-EE1[3];
	  BB1[1]=-W1[1];
	  BB1[2]=-W1[2];
	  BB1[3]=-W1[3];
  	  
	  double KT,phb,GG1[5],G[4],KT1,phb2, GT[4],thetr,BT[4], GGT[4],TG1[4];
	  KT1 =0.; phb2 =0.;
	 if (POWdecay) {
	     GT[1]=0.;
	     GT[2]=0.;
	     GT[3]=0.;
	     TG1[1]=0.;
	     TG1[2]=0.;
	     TG1[3]=0.;
	       KT=sqrt(pt2r1*pow(mt,2));// transverse momenta for decay gluon
	  phb=2.*pi*random(seed);
	  B1[1]=-KT*cos(phb);
	  B1[2]=-KT*sin(phb);
	  B1[3]=-sqrt(pow(Eb,2)-pow(B1[1],2)-pow(B1[2],2));
	 G[1]=KT*cos(phb);
  	  G[2]=KT*sin(phb);
	  G[3]=-sqrt(pow(Eg,2)-pow(G[1],2)-pow(G[2],2));
	  // If truncated emission
	  if (!demitt) {
	    if (fabs(pw-sqrt(pow(B1[1]+G[1],2)+pow(B1[2]+G[2],2)+pow(B1[3]+G[3],2))) > 0.1) {G[3]=sqrt(pow(Eg,2)-pow(G[1],2)-pow(G[2],2));} } 
	  if (demitt) {
	    if (fabs(sqrt(pow(xww1*mt/2.,2)-pow(mw,2))-sqrt(pow(B1[1]+G[1],2)+pow(B1[2]+G[2],2)+pow(B1[3]+G[3],2))) > 0.1) {G[3]=sqrt(pow(Eg,2)-pow(G[1],2)-pow(G[2],2));}}
	  	  if (demitt) {
	  if (pow(B1[3]+G[3],2) < ptr2) {demitt = false; pw = sqrt (pow(xww1*mt/2.,2)-pow(mw,2));
	 

	// top rest frame
	  W1[1]=pw*sin(acos(ctw))*cos(phi);
	  W1[2]=pw*sin(acos(ctw))*sin(phi);
	  W1[3]=pw*(ctw);
	//W rest frame
	E[1]=El*sin(acos(ctl))*cos(bphi);
	E[2]=El*sin(acos(ctl))*sin(bphi);
	E[3]=El*ctl;
	 vv=pw/(xww1*mt/2.);
	 gammaw=1./sqrt(1.-pow(vv,2));
        //top rest frame 
	E[3]=gammaw*(E[3]+vv*El);

	 EE1[3]=E[3]*ctw+E[1]*sin(acos(ctw));
	 EE1[1]=E[3]*sin(acos(ctw))-E[1]*ctw;


      	  EE1[0]=EE1[1];

	  EE1[1]=EE1[0]*cos(phi)-E[2]*sin(phi);
	  EE1[2]=E[2]*cos(phi)+EE1[0]*sin(phi);
	  NN1[1]=W1[1]-EE1[1];
	  NN1[2]=W1[2]-EE1[2];
	  NN1[3]=W1[3]-EE1[3];
	 
	  B1[1]=-KT*cos(phb);
	  B1[2]=-KT*sin(phb);
	  B1[3]=-sqrt(pow(Eb,2)-pow(B1[1],2)-pow(B1[2],2));
	  G[1]=KT*cos(phb);
  	  G[2]=KT*sin(phb);
	  G[3]=-sqrt(pow(Eg,2)-pow(G[1],2)-pow(G[2],2));

    if (fabs(pw-sqrt(pow(B1[1]+G[1],2)+pow(B1[2]+G[2],2)+pow(B1[3]+G[3],2))) > 0.1) {G[3]=sqrt(pow(Eg,2)-pow(G[1],2)-pow(G[2],2));} 
  	  }
	  } 
	  
	   if(demitt) {
	   KT1 = sqrt(ptr2);// transverse momenta for decay truncated gluon
	    phb2 = 2.*pi*random(seed); 
	     BT[1]=KT1*cos(phb2);
	     BT[2]=KT1*sin(phb2);
	     BT[3]= -sqrt(pow(B1[3]+G[3],2)-pow(KT1,2));
	     GT[1]=-KT1*cos(phb2);
	     GT[2]=-KT1*sin(phb2);
	     GT[3]= -pw - BT[3];
}
	   thetr = asin(KT1/(B1[3]+G[3])); // angle of truncated emission
	   double ctwr = cos(thetr);
	  if (demitt) {
	  BT[3]=B1[3]*ctwr+B1[1]*sin(acos(ctwr));
	  BT[1]=-B1[3]*sin(acos(ctwr))+B1[1]*ctwr;
	  GGT[3]=G[3]*ctwr+G[1]*sin(acos(ctwr));
	  GGT[1]=-G[3]*sin(acos(ctwr))+G[1]*ctwr;
          BT[0]=BT[1];
	  BT[1]=BT[0]*cos(phb2)-B1[2]*sin(phb2);
	  BT[2]=B1[2]*cos(phb2)+BT[0]*sin(phb2);
	  GGT[0]=GGT[1];
	  GGT[1]=GGT[0]*cos(phb2)-G[2]*sin(phb2);
	  GGT[2]=G[2]*cos(phb2)+GGT[0]*sin(phb2);
	  B1[3] = BT[3];
	  B1[1] = BT[1];
	  B1[2] = BT[2];
	  G[3] = GGT[3];
	  G[2] = GGT[2];
	  G[1] = GGT[1];
}	
	  // Rotate through W angle in top rest frame 
	  BB1[3]=B1[3]*ctw+B1[1]*sin(acos(ctw));
	  BB1[1]=B1[3]*sin(acos(ctw))-B1[1]*ctw;
	  GG1[3]=G[3]*ctw+G[1]*sin(acos(ctw));
	  GG1[1]=G[3]*sin(acos(ctw))-G[1]*ctw;
	  TG1[3]=GT[3]*ctw+GT[1]*sin(acos(ctw));
	  TG1[1]=GT[3]*sin(acos(ctw))-GT[1]*ctw;
          BB1[0]=BB1[1];
	  BB1[1]=BB1[0]*cos(phi)-B1[2]*sin(phi);
	  BB1[2]=B1[2]*cos(phi)+BB1[0]*sin(phi);
	  GG1[0]=GG1[1];
	  GG1[1]=GG1[0]*cos(phi)-G[2]*sin(phi);
	  GG1[2]=G[2]*cos(phi)+GG1[0]*sin(phi);
	  TG1[0]=TG1[1];
	  TG1[1]=TG1[0]*cos(phi)-GT[2]*sin(phi);
	  TG1[2]=GT[2]*cos(phi)+TG1[0]*sin(phi);
	  if (fabs(W1[1]+GG1[1]+BB1[1]+TG1[1]) > 0.1 || fabs(W1[2]+GG1[2]+BB1[2]+TG1[2]) >0.1) {B1[3]=sqrt(pow(Eb,2)-pow(KT,2));
	   G[3]=-sqrt(pow(Eg,2)-pow(KT,2)); B1[1]=-KT*cos(phb) ; B1[2] = -KT*sin(phb); G[1]=KT*cos(phb); G[2] = KT*sin(phb); thetr = asin(KT1/(sqrt(pow(B1[1]+G[1],2)+pow(B1[2]+G[2],2)+pow(B1[3]+G[3],2)))); 
	   ctwr = cos(thetr);
	
	  if (demitt) {
	    // rotate through truncated emission angle
	  BT[3]=B1[3]*ctwr+B1[1]*sin(acos(ctwr));
	  BT[1]=-B1[3]*sin(acos(ctwr))+B1[1]*ctwr;
	  GGT[3]=G[3]*ctwr+G[1]*sin(acos(ctwr));
	  GGT[1]=-G[3]*sin(acos(ctwr))+G[1]*ctwr;
          BT[0]=BT[1];
	  BT[1]=BT[0]*cos(phb2)-B1[2]*sin(phb2);
	  BT[2]=B1[2]*cos(phb2)+BT[0]*sin(phb2);
	  GGT[0]=GGT[1];
	  GGT[1]=GGT[0]*cos(phb2)-G[2]*sin(phb2);
	  GGT[2]=G[2]*cos(phb2)+GGT[0]*sin(phb2);
	  B1[3] = BT[3];
	  B1[1] = BT[1];
	  B1[2] = BT[2];
	  G[3] = GGT[3];
	  G[2] = GGT[2];
	  G[1] = GGT[1];}
	  
	  BB1[3]=B1[3]*ctw+B1[1]*sin(acos(ctw));
	  BB1[1]=B1[3]*sin(acos(ctw))-B1[1]*ctw;
	  GG1[3]=G[3]*ctw+G[1]*sin(acos(ctw));
	  GG1[1]=G[3]*sin(acos(ctw))-G[1]*ctw;
          TG1[3]=GT[3]*ctw+GT[1]*sin(acos(ctw));
	  TG1[1]=GT[3]*sin(acos(ctw))-GT[1]*ctw;
          BB1[0]=BB1[1];
	  BB1[1]=BB1[0]*cos(phi)-B1[2]*sin(phi);
	  BB1[2]=B1[2]*cos(phi)+BB1[0]*sin(phi);
	  GG1[0]=GG1[1];
	  GG1[1]=GG1[0]*cos(phi)-G[2]*sin(phi);
	  GG1[2]=G[2]*cos(phi)+GG1[0]*sin(phi);
	  TG1[0]=TG1[1];
	  TG1[1]=TG1[0]*cos(phi)-GT[2]*sin(phi);
	  TG1[2]=GT[2]*cos(phi)+TG1[0]*sin(phi);
	   }
	 }
	 if (POWdecay && demitt){
	   // Boost momenta to preserve top mass.
	 double bp = sqrt(pow(BB1[3],2)+pow(BB1[2],2)+pow(BB1[1],2));
	 double gp = sqrt(pow(GG1[3],2)+pow(GG1[2],2)+pow(GG1[1],2));
	 double gtrp = sqrt(pow(TG1[3],2)+pow(TG1[2],2)+pow(TG1[1],2));	
	 double wp = sqrt(pow(W1[3],2)+pow(W1[2],2)+pow(W1[1],2));
	double mfac = (175.*gp + 175.*bp+175.*gtrp-sqrt(-pow(wp*mw,2)+30625.*pow(wp,2)+pow(bp*mw,2)+2.*bp*gp*pow(mw,2)+2.*bp*gtrp*pow(mw,2)+pow(gp*mw,2)+2.*gp*gtrp*pow(mw,2)+pow(gtrp*mw,2)))/(-pow(wp,2)+pow(bp,2)+2.*bp*gp+2.*bp*gtrp+pow(gp,2)+2.*gp*gtrp+pow(gtrp,2));		  
	BB1[3]=mfac*BB1[3]; BB1[2]=mfac*BB1[2]; BB1[1]=mfac*BB1[1];GG1[3]=mfac*GG1[3]; GG1[2]=mfac*GG1[2]; GG1[1]=mfac*GG1[1];
         TG1[3]=mfac*TG1[3]; TG1[2]=mfac*TG1[2]; TG1[1]=mfac*TG1[1];W1[3]=mfac*W1[3]; W1[2]=mfac*W1[2]; W1[1]=mfac*W1[1];
        E[1]=El*sin(acos(ctl))*cos(bphi);
 	E[2]=El*sin(acos(ctl))*sin(bphi);
	E[3]=El*ctl;
	double vv=sqrt((pow(W1[1],2)+pow(W1[2],2)+pow(W1[3],2))/(pow(W1[1],2)+pow(W1[2],2)+pow(W1[3],2)+pow(mw,2)));
	double gammaw=1./sqrt(1.-pow(vv,2));
        //top rest frame 
	E[3]=gammaw*(E[3]+vv*El);
	
	  EE1[3]=E[3]*ctw+E[1]*sin(acos(ctw));
	  EE1[1]=E[3]*sin(acos(ctw))-E[1]*ctw;
      	  EE1[0]=EE1[1];
	  EE1[1]=EE1[0]*cos(phi)-E[2]*sin(phi);
	  EE1[2]=E[2]*cos(phi)+EE1[0]*sin(phi);
	  NN1[1]=W1[1]-EE1[1];
	  NN1[2]=W1[2]-EE1[2];
	  NN1[3]=W1[3]-EE1[3];
}
	 // Boost to labframe
	  W1[3]=gammat*(W1[3]+vt*sqrt(pow(W1[3],2)+pow(W1[1],2)+pow(W1[2],2)+pow(mw,2)));
	  EE1[3]=gammat*(EE1[3]+vt*sqrt(pow(EE1[3],2)+pow(EE1[1],2)+pow(EE1[2],2)));
	  BB1[3]=gammat*(BB1[3]+vt*sqrt(pow(BB1[3],2)+pow(BB1[1],2)+pow(BB1[2],2)));

	  if(POWdecay && demitt) {TG1[3]=gammat*(TG1[3]+vt*sqrt(pow(TG1[3],2)+pow(TG1[1],2)+pow(TG1[2],2)));}
	  double WNE,WNNE,EEE,BBB,EEEE,BBBB,GGG,GGGG;
	  double phcm = 2.*pi*random(seed);
	   double cphcm = cos(phcm);
	   double sphcm = sin(phcm);
	  if(!POWprod){
	  T1[1]=250.*beta*sint*cphcm;
	  T1[2]=250.*beta*sint*sphcm;
	  T1[3]=250.*beta*cost;
	  PUP[1][3]=T1[1];
	  PUP[2][3]=T1[2];
	  PUP[3][3]=T1[3];
	  costx=T1[1]/sqrt(pow(T1[1],2)+pow(T1[2],2));
	  sintx=sin(acos(costx));
	  costz=T1[3]/sqrt(pow(T1[1],2)+pow(T1[2],2)+pow(T1[3],2));
	  sintz=sin(acos(costz));
	 WNE=W1[1];
	W1[1]=-WNE*costz+W1[3]*sintz;
	W1[3]=WNE*sintz+W1[3]*costz;
	WNNE=W1[1];
	W1[1]=WNNE*cos(phcm)-W1[2]*sin(phcm);
	W1[2]=WNNE*sin(phcm)+W1[2]*cos(phcm);
	EEE=EE1[1];
	EE1[1]=-EEE*costz+EE1[3]*sintz;
	EE1[3]=EEE*sintz+EE1[3]*costz;
	 EEEE=EE1[1];
	EE1[1]=EEEE*cos(phcm)-EE1[2]*sin(phcm);
	EE1[2]=EEEE*sin(phcm)+EE1[2]*cos(phcm);
	BBB=BB1[1];
	BB1[1]=-BBB*costz+BB1[3]*sintz;
	BB1[3]=BBB*sintz+BB1[3]*costz;
	BBBB=BB1[1];
	BB1[1]=BBBB*cos(phcm)-BB1[2]*sin(phcm);
	BB1[2]=BBBB*sin(phcm)+BB1[2]*cos(phcm);
	G1[1]=T1[1]-W1[1]-BB1[1];
        G1[2]=T1[2]-W1[2]-BB1[2];
	G1[3]=T1[3]-W1[3]-BB1[3];
	if (POWdecay && demitt) {
	  double GGG;
	GGG=TG1[1];
	TG1[1] = -GGG*costz+TG1[3]*sintz;
	TG1[3] = GGG*sintz+TG1[3]*costz;
	GGGG=TG1[1];
	TG1[1]=GGGG*cos(phcm)-TG1[2]*sin(phcm);
	TG1[2]=GGGG*sin(phcm)+TG1[2]*cos(phcm);
	G1[1]=T1[1]-W1[1]-BB1[1]-TG1[1];
	G1[2]=T1[2]-W1[2]-BB1[2]-TG1[2];
	G1[3]=T1[3]-W1[3]-BB1[3]-TG1[3];
}
}else {
	    T1[4]=T1[3];
	    T1[5]=T1[3]*sintz;
	    T1[3]=T1[3]*costz;
	    T1[1]=T1[5]*costx;
	   if(fabs(T1[1]-PUP[1][4]) > 0.001) {T1[1]=-T1[1];sintz=-sintz;}
	       T1[5]=T1[4]*sintz;
	       T1[2]=-T1[5]*sintx;
	     if(fabs(T1[2]-PUP[2][4]) > 0.001) {T1[2]=-T1[2];sintx=-sintx;} 
	     WNE=W1[1];
	     W1[1]=WNE*costz+W1[3]*sintz;
	W1[3]=-WNE*sintz+W1[3]*costz;
	WNNE=W1[1];
	W1[1]=WNNE*costx+W1[2]*sintx;
	W1[2]=-WNNE*sintx+W1[2]*costx;
	EEE=EE1[1];
	EE1[1]=EEE*costz+EE1[3]*sintz;
	EE1[3]=-EEE*sintz+EE1[3]*costz;
	 EEEE=EE1[1];
	EE1[1]=EEEE*costx+EE1[2]*sintx;
	EE1[2]=-EEEE*sintx+EE1[2]*costx;
	BBB=BB1[1];
	BB1[1]=BBB*costz+BB1[3]*sintz;
	BB1[3]=-BBB*sintz+BB1[3]*costz;
	BBBB=BB1[1];
	BB1[1]=BBBB*costx+BB1[2]*sintx;
	BB1[2]=-BBBB*sintx+BB1[2]*costx;
	G1[1]=T1[1]-W1[1]-BB1[1];
        G1[2]=T1[2]-W1[2]-BB1[2];
	G1[3]=T1[3]-W1[3]-BB1[3];
	if (POWdecay && demitt){
	GGG=TG1[1];
	TG1[1]=GGG*costz+TG1[3]*sintz;
	TG1[3]=-GGG*sintz+TG1[3]*costz;
	GGGG=TG1[1];
	TG1[1]=GGGG*costx+TG1[2]*sintx;
	TG1[2]=-GGGG*sintx+TG1[2]*costx;
	G1[1]=T1[1]-W1[1]-BB1[1]-TG1[1];
	G1[2]=T1[2]-W1[2]-BB1[2]-TG1[2];
	G1[3]=T1[3]-W1[3]-BB1[3]-TG1[3];
	 }
	  }
        NN1[1]=W1[1]-EE1[1];
	NN1[2]=W1[2]-EE1[2];
	NN1[3]=W1[3]-EE1[3];
 	 if(!POWdecay) {
	BB1[1]=T1[1]-W1[1];
	BB1[2]=T1[2]-W1[2];
       	BB1[3]=T1[3]-W1[3];} 
	 pw=sqrt(pow(Ew2,2)-pow(mw,2));
         phi=2.*pi*random(seed);

	 // Do the same as above for tbar

	// tbar rest frame

	  W2[1]=pw*sin(acos(ctw2))*cos(phi);
	  W2[2]=pw*sin(acos(ctw2))*sin(phi);
	  W2[3]=pw*(ctw2);
	 
	  //W rest frame

	E[1]=El*sin(acos(ctl2))*cos(bphi2);
	E[2]=El*sin(acos(ctl2))*sin(bphi2);
	E[3]=El*ctl2;
	 vv=pw/Ew2;
	 gammaw=1./sqrt(1.-pow(vv,2));
        //tbar rest frame 
	E[3]=gammaw*(E[3]+vv*El);

	 EE2[3]=E[3]*ctw2+E[1]*sin(acos(ctw2));
	 EE2[1]=E[3]*sin(acos(ctw2))-E[1]*ctw2;


      	  EE2[0]=EE2[1];

	  EE2[1]=EE2[0]*cos(phi)-E[2]*sin(phi);
	  EE2[2]=E[2]*cos(phi)+EE2[0]*sin(phi);
	  NN2[1]=W2[1]-EE2[1];
	  NN2[2]=W2[2]-EE2[2];
	  NN2[3]=W2[3]-EE2[3];

	  double GG2[5],TG2[4];;
	 if (POWdecay) {
	   GT[1]=0.;
	   GT[2]=0.;
	   GT[3]=0.;
	   KT=sqrt(pt2r2*pow(mt,2)); // transverse momenta for decay gluon
	   phb=2.*pi*random(seed);
          B2[1]=-KT*cos(phb);
	  B2[2]=-KT*sin(phb);
	  B2[3]=-sqrt(pow(Eb2,2)-pow(B2[1],2)-pow(B2[2],2));
	  G[1]=KT*cos(phb);
  	  G[2]=KT*sin(phb);
	  G[3]=-sqrt(pow(Eg2,2)-pow(G[1],2)-pow(G[2],2));

          if (!demittb) {
	    if (fabs(pw-sqrt(pow(B2[1]+G[1],2)+pow(B2[2]+G[2],2)+pow(B2[3]+G[3],2))) > 0.1) {G[3]=sqrt(pow(Eg2,2)-pow(G[1],2)-pow(G[2],2));} } 
	  if (demittb) {
	    if (fabs(sqrt(pow(xww2*mt/2.,2)-pow(mw,2))-sqrt(pow(B2[1]+G[1],2)+pow(B2[2]+G[2],2)+pow(B2[3]+G[3],2))) > 0.1) {G[3]=sqrt(pow(Eg2,2)-pow(G[1],2)-pow(G[2],2));}}

	  if (demittb) {
	  if (pow(B2[3]+G[3],2) < ptr2b) {demittb = false; pw = sqrt (pow(xww2*mt/2.,2)-pow(mw,2));
	 
	  
	// tbar rest frame
	  W2[1]=pw*sin(acos(ctw2))*cos(phi);
	  W2[2]=pw*sin(acos(ctw2))*sin(phi);
	  W2[3]=pw*(ctw2);
	//W rest frame
	E[1]=El*sin(acos(ctl2))*cos(bphi2);
	E[2]=El*sin(acos(ctl2))*sin(bphi2);
	E[3]=El*ctl2;
	 vv=pw/(xww2*mt/2.);
	 gammaw=1./sqrt(1.-pow(vv,2));
        //tbar rest frame 
	E[3]=gammaw*(E[3]+vv*El);

	 EE2[3]=E[3]*ctw2+E[1]*sin(acos(ctw2));
	 EE2[1]=E[3]*sin(acos(ctw2))-E[1]*ctw2;


      	  EE2[0]=EE2[1];

	  EE2[1]=EE2[0]*cos(phi)-E[2]*sin(phi);
	  EE2[2]=E[2]*cos(phi)+EE2[0]*sin(phi);
	  NN2[1]=W2[1]-EE2[1];
	  NN2[2]=W2[2]-EE2[2];
	  NN2[3]=W2[3]-EE2[3];
	 
	  B2[1]=-KT*cos(phb);
	  B2[2]=-KT*sin(phb);
	  B2[3]=-sqrt(pow(Eb2,2)-pow(B2[1],2)-pow(B2[2],2));
	G[1]=KT*cos(phb);
  	  G[2]=KT*sin(phb);
	  G[3]=-sqrt(pow(Eg2,2)-pow(G[1],2)-pow(G[2],2));

    if (fabs(pw-sqrt(pow(B2[1]+G[1],2)+pow(B2[2]+G[2],2)+pow(B2[3]+G[3],2))) > 0.1) {G[3]=sqrt(pow(Eg2,2)-pow(G[1],2)-pow(G[2],2));} 
  	  }
	  } 
	   if(demittb) {
	     KT1 = sqrt(ptr2b);// transverse momenta for decay truncated gluon

	     phb2 = 2.*pi*random(seed); 
	     BT[1]=KT1*cos(phb2);
	     BT[2]=KT1*sin(phb2);
	     BT[3]= -sqrt(pow(B2[3]+G[3],2)-pow(KT1,2));
	     GT[1]=-KT1*cos(phb2);
	     GT[2]=-KT1*sin(phb2);
	     GT[3]= -pw - BT[3];
}
	    thetr = asin(KT1/(B2[3]+G[3]));
	    double ctwr = cos(thetr);
	  if (POWdecay &&demittb) {
	  BT[3]=B2[3]*ctwr+B2[1]*sin(acos(ctwr));
	  BT[1]=-B2[3]*sin(acos(ctwr))+B2[1]*ctwr;
	  GGT[3]=G[3]*ctwr+G[1]*sin(acos(ctwr));
	  GGT[1]=-G[3]*sin(acos(ctwr))+G[1]*ctwr;
          BT[0]=BT[1];
	  BT[1]=BT[0]*cos(phb2)-B2[2]*sin(phb2);
	  BT[2]=B2[2]*cos(phb2)+BT[0]*sin(phb2);
	  GGT[0]=GGT[1];
	  GGT[1]=GGT[0]*cos(phb2)-G[2]*sin(phb2);
	  GGT[2]=G[2]*cos(phb2)+GGT[0]*sin(phb2);
	  B2[3] = BT[3];
	  B2[1] = BT[1];
	  B2[2] = BT[2];
	  G[3] = GGT[3];
	  G[2] = GGT[2];
	  G[1] = GGT[1];

}	
          //detour to W rest frame
	  B2[0]=gammaw*(B2[3]-vv*Eb2);
	 
	  BB2[3]=B2[3]*ctw2+B2[1]*sin(acos(ctw2));
	  BB2[1]=B2[3]*sin(acos(ctw2))-B2[1]*ctw2;
	  GG2[3]=G[3]*ctw2+G[1]*sin(acos(ctw2));
	  GG2[1]=G[3]*sin(acos(ctw2))-G[1]*ctw2;
	  TG2[3]=GT[3]*ctw2+GT[1]*sin(acos(ctw2));
	  TG2[1]=GT[3]*sin(acos(ctw2))-GT[1]*ctw2;
	  GG2[0]=GG2[1];
	  GG2[1]=GG2[0]*cos(phi)-G[2]*sin(phi);
	  GG2[2]=G[2]*cos(phi)+GG2[0]*sin(phi);
          BB2[0]=BB2[1];
	  BB2[1]=BB2[0]*cos(phi)-B2[2]*sin(phi);
	  BB2[2]=B2[2]*cos(phi)+BB2[0]*sin(phi);
	  TG2[0]=TG2[1];
	  TG2[1]=TG2[0]*cos(phi)-GT[2]*sin(phi);
	  TG2[2]=GT[2]*cos(phi)+TG2[0]*sin(phi);
	  
	  if (fabs(W2[1]+GG2[1]+BB2[1]+TG2[1]) > 0.1 || fabs(W2[2]+GG2[2]+BB2[2]+TG2[2]) >0.1) {B2[3]=sqrt(pow(Eb2,2)-pow(KT,2));
	  G[3]=-sqrt(pow(Eg2,2)-pow(KT,2)); B2[1]=-KT*cos(phb) ; B2[2] = -KT*sin(phb); G[1]=KT*cos(phb); G[2] = KT*sin(phb); thetr = asin(KT1/(sqrt(pow(B2[1]+G[1],2)+pow(B2[2]+G[2],2)+pow(B2[3]+G[3],2))));
	  
	  ctwr = cos(thetr);  
          if (POWdecay &&demittb) {
	  BT[3]=B2[3]*ctwr+B2[1]*sin(acos(ctwr));
	  BT[1]=-B2[3]*sin(acos(ctwr))+B2[1]*ctwr;
	  GGT[3]=G[3]*ctwr+G[1]*sin(acos(ctwr));
	  GGT[1]=-G[3]*sin(acos(ctwr))+G[1]*ctwr;
          BT[0]=BT[1];
	  BT[1]=BT[0]*cos(phb2)-B2[2]*sin(phb2);
	  BT[2]=B2[2]*cos(phb2)+BT[0]*sin(phb2);
	  GGT[0]=GGT[1];
	  GGT[1]=GGT[0]*cos(phb2)-G[2]*sin(phb2);
	  GGT[2]=G[2]*cos(phb2)+GGT[0]*sin(phb2);
	  B2[3] = BT[3];
	  B2[1] = BT[1];
	  B2[2] = BT[2];
	  G[3] = GGT[3];
	  G[2] = GGT[2];
	  G[1] = GGT[1];}

          BB2[3]=B2[3]*ctw2+B2[1]*sin(acos(ctw2));
	  BB2[1]=B2[3]*sin(acos(ctw2))-B2[1]*ctw2;
	  GG2[3]=G[3]*ctw2+G[1]*sin(acos(ctw2));
	  GG2[1]=G[3]*sin(acos(ctw2))-G[1]*ctw2;
	  TG2[3]=GT[3]*ctw2+GT[1]*sin(acos(ctw2));
	  TG2[1]=GT[3]*sin(acos(ctw2))-GT[1]*ctw2;
	  GG2[0]=GG2[1];
	  GG2[1]=GG2[0]*cos(phi)-G[2]*sin(phi);
	  GG2[2]=G[2]*cos(phi)+GG2[0]*sin(phi);
          BB2[0]=BB2[1];
	  BB2[1]=BB2[0]*cos(phi)-B2[2]*sin(phi);
	  BB2[2]=B2[2]*cos(phi)+BB2[0]*sin(phi);
	  TG2[0]=TG2[1];
	  TG2[1]=TG2[0]*cos(phi)-GT[2]*sin(phi);
	  TG2[2]=GT[2]*cos(phi)+TG2[0]*sin(phi);
	  }	
	 }
	 if (POWdecay && demittb){
         double bp = sqrt(pow(BB2[3],2)+pow(BB2[2],2)+pow(BB2[1],2));
	 double gp = sqrt(pow(GG2[3],2)+pow(GG2[2],2)+pow(GG2[1],2));
	 double gtrp = sqrt(pow(TG2[3],2)+pow(TG2[2],2)+pow(TG2[1],2));	
	 double wp = sqrt(pow(W2[3],2)+pow(W2[2],2)+pow(W2[1],2));
	 double mfac = (175.*gp + 175.*bp+175.*gtrp-sqrt(-pow(wp*mw,2)+30625.*pow(wp,2)+pow(bp*mw,2)+2.*bp*gp*pow(mw,2)+2.*bp*gtrp*pow(mw,2)+pow(gp*mw,2)+2.*gp*gtrp*pow(mw,2)+pow(gtrp*mw,2)))/(-pow(wp,2)+pow(bp,2)+2.*bp*gp+2.*bp*gtrp+pow(gp,2)+2.*gp*gtrp+pow(gtrp,2));		  
	 BB2[3]=mfac*BB2[3]; BB2[2]=mfac*BB2[2]; BB2[1]=mfac*BB2[1];GG2[3]=mfac*GG2[3]; GG2[2]=mfac*GG2[2]; GG2[1]=mfac*GG2[1];
         TG2[3]=mfac*TG2[3]; TG2[2]=mfac*TG2[2]; TG2[1]=mfac*TG2[1];W2[3]=mfac*W2[3]; W2[2]=mfac*W2[2]; W2[1]=mfac*W2[1];
	 E[1]=El*sin(acos(ctl2))*cos(bphi2);
 	E[2]=El*sin(acos(ctl2))*sin(bphi2);
	E[3]=El*ctl2;
	double vv=sqrt((pow(W2[1],2)+pow(W2[2],2)+pow(W2[3],2))/(pow(W2[1],2)+pow(W2[2],2)+pow(W2[3],2)+pow(mw,2)));
	double gammaw=1./sqrt(1.-pow(vv,2));
        //t bar rest frame 
	E[3]=gammaw*(E[3]+vv*El);
	
	  EE2[3]=E[3]*ctw2+E[1]*sin(acos(ctw2));
	  EE2[1]=E[3]*sin(acos(ctw2))-E[1]*ctw2;
      	  EE2[0]=EE2[1];
	  EE2[1]=EE2[0]*cos(phi)-E[2]*sin(phi);
	  EE2[2]=E[2]*cos(phi)+EE2[0]*sin(phi);
	  NN2[1]=W2[1]-EE2[1];
	  NN2[2]=W2[2]-EE2[2];
	  NN2[3]=W2[3]-EE2[3];
	 
 }

	   //  lab frame
	  W2[3]=gammatb*(W2[3]+vtb*sqrt(pow(W2[3],2)+pow(W2[1],2)+pow(W2[2],2)+pow(mw,2)));
	  EE2[3]=gammatb*(EE2[3]+vtb*sqrt(pow(EE2[3],2)+pow(EE2[1],2)+pow(EE2[2],2)));
	  BB2[3]=gammatb*(BB2[3]+vtb*sqrt(pow(BB2[3],2)+pow(BB2[1],2)+pow(BB2[2],2)));
	   if(POWdecay && demittb) {TG2[3]=gammatb*(TG2[3]+vtb*sqrt(pow(TG2[3],2)+pow(TG2[1],2)+pow(TG2[2],2)));}
	  
	  if(!POWprod){
	  
	  T2[1]=-250.*beta*sint*cphcm;
	  T2[2]=-250.*beta*sint*sphcm;
	  T2[3]=-250.*beta*cost;
	  PUP[1][4]=T2[1];
	  PUP[2][4]=T2[2];
	  PUP[3][4]=T2[3];
	  costbx=T2[1]/sqrt(pow(T2[1],2)+pow(T2[2],2));
	  sintbx = sin(acos(costbx));
	  costbz=-T2[3]/sqrt(pow(T2[1],2)+pow(T2[2],2)+pow(T2[3],2));
	  sintbz = sin(acos(costbz));

	WNE=W2[1];	  
	W2[1]=(-WNE*costbz+W2[3]*sintbz);
	W2[3]=(WNE*sintbz+W2[3]*costbz);
	WNNE=W2[1];
	W2[1]=(WNNE*cos(phcm)-W2[2]*sin(phcm));
	W2[2]=(WNNE*sin(phcm)+W2[2]*cos(phcm));
	EEE=EE2[1];
	EE2[1]=(-EEE*costbz+EE2[3]*sintbz);
	EE2[3]=(EEE*sintbz+EE2[3]*costbz);
	EEEE=EE2[1];
	EE2[1]=(EEEE*cos(phcm)-EE2[2]*sin(phcm));
	EE2[2]=(EEEE*sin(phcm)+EE2[2]*cos(phcm));
	BBB=BB2[1];
	BB2[1]=(-BBB*costbz+BB2[3]*sintbz);
	BB2[3]=(BBB*sintbz+BB2[3]*costbz);
	BBBB=BB2[1];
	BB2[1]=(BBBB*cos(phcm)-BB2[2]*sin(phcm));
	BB2[2]=(BBBB*sin(phcm)+BB2[2]*cos(phcm));
	W2[1]= -W2[1];
	W2[2]= -W2[2];
	W2[3]= -W2[3];
	EE2[1]= -EE2[1];
        EE2[2]= -EE2[2];
	EE2[3]= -EE2[3];
	BB2[1]= -BB2[1];
        BB2[2]= -BB2[2];
	BB2[3]= -BB2[3],
	G2[1]=T2[1]-W2[1]-BB2[1];
        G2[2]=T2[2]-W2[2]-BB2[2];
	G2[3]=T2[3]-W2[3]-BB2[3];
	if (POWdecay && demittb) {
	GGG=TG2[1];
	TG2[1] = -GGG*costbz+TG2[3]*sintbz;
        TG2[3] = GGG*sintbz+TG2[3]*costbz;
	GGGG=TG2[1];
	TG2[1]=(GGGG*cos(phcm)-TG2[2]*sin(phcm));
	TG2[2]=(GGGG*sin(phcm)+TG2[2]*cos(phcm));
	TG2[1]= -TG2[1];
        TG2[2]= -TG2[2];
	TG2[3]= -TG2[3],
	G2[1]=T2[1]-W2[1]-BB2[1]-TG2[1];
	G2[2]=T2[2]-W2[2]-BB2[2]-TG2[2];
	G2[3]=T2[3]-W2[3]-BB2[3]-TG2[3];
}
	}else 
	    {
	    T2[4]=T2[3];
	    T2[5]=T2[3]*sintbz;
	    T2[3]=T2[3]*costbz;
	    T2[1]=T2[5]*costbx;
	    if(fabs(T2[1]-PUP[1][5]) > 0.001) {T2[1]=-T2[1];sintbz=-sintbz;}
	       T2[5]=T2[4]*sintbz;
	       T2[2]=-T2[5]*sintbx;
	     if(fabs(T2[2]-PUP[2][5]) > 0.001) {T2[2]=-T2[2];sintbx=-sintbx;} 
	     WNE=W2[1];
	     W2[1]=WNE*costbz+W2[3]*sintbz;
	W2[3]=-WNE*sintbz+W2[3]*costbz;
	WNNE=W2[1];
	W2[1]=WNNE*costbx+W2[2]*sintbx;
	W2[2]=-WNNE*sintbx+W2[2]*costbx;
	EEE=EE2[1];
	EE2[1]=EEE*costbz+EE2[3]*sintbz;
	EE2[3]=-EEE*sintbz+EE2[3]*costbz;
	EEEE=EE2[1];
	EE2[1]=EEEE*costbx+EE2[2]*sintbx;
	EE2[2]=-EEEE*sintbx+EE2[2]*costbx;
	BBB=BB2[1];
	BB2[1]=BBB*costbz+BB2[3]*sintbz;
	BB2[3]=-BBB*sintbz+BB2[3]*costbz;
	BBBB=BB2[1];
	BB2[1]=BBBB*costbx+BB2[2]*sintbx;
	BB2[2]=-BBBB*sintbx+BB2[2]*costbx;
	G2[1]=T2[1]-W2[1]-BB2[1];
        G2[2]=T2[2]-W2[2]-BB2[2];
	G2[3]=T2[3]-W2[3]-BB2[3];
		if (POWdecay && demittb){
	GGG=TG2[1];
	TG2[1]=GGG*costbz+TG2[3]*sintbz;
	TG2[3]=-GGG*sintbz+TG2[3]*costbz;
	GGGG=TG2[1];
	TG2[1]=GGGG*costbx+TG2[2]*sintbx;
	TG2[2]=-GGGG*sintbx+TG2[2]*costbx;
	G2[1]=T2[1]-W2[1]-BB2[1]-TG2[1];
	G2[2]=T2[2]-W2[2]-BB2[2]-TG2[2];
	G2[3]=T2[3]-W2[3]-BB2[3]-TG2[3];
	}
	
}

	NN2[1]=W2[1]-EE2[1];
	NN2[2]=W2[2]-EE2[2];
	NN2[3]=W2[3]-EE2[3];
 	 if(!POWdecay) {
 	BB2[1]=T2[1]-W2[1];
 	BB2[2]=T2[2]-W2[2];
 	BB2[3]=T2[3]-W2[3];}

	
	double TE1=sqrt(pow(T1[1],2)+pow(T1[2],2)+pow(T1[3],2)+pow(mt,2));
	double TE2=sqrt(pow(T2[1],2)+pow(T2[2],2)+pow(T2[3],2)+pow(mt,2));
	double NE1=sqrt(pow(NN1[1],2)+pow(NN1[2],2)+pow(NN1[3],2));
	double NE2=sqrt(pow(NN2[1],2)+pow(NN2[2],2)+pow(NN2[3],2));
	double GE1=sqrt(pow(G1[1],2)+pow(G1[2],2)+pow(G1[3],2));
	double GE2=sqrt(pow(G2[1],2)+pow(G2[2],2)+pow(G2[3],2));
       	double BE1=sqrt(pow(BB1[1],2)+pow(BB1[2],2)+pow(BB1[3],2));
	double BE2=sqrt(pow(BB2[1],2)+pow(BB2[2],2)+pow(BB2[3],2));
	double LE1=sqrt(pow(EE1[1],2)+pow(EE1[2],2)+pow(EE1[3],2));
	double LE2=sqrt(pow(EE2[1],2)+pow(EE2[2],2)+pow(EE2[3],2));
	
	if (POWdecay &&demitt) {PUP[4][15]=sqrt(pow(TG1[1],2)+pow(TG1[2],2)+pow(TG1[3],2));PUP[1][15]=TG1[1]; PUP[2][15]=TG1[2]; PUP[3][15]=TG1[3];}
	if (POWdecay &&demittb){PUP[4][16]=sqrt(pow(TG2[1],2)+pow(TG2[2],2)+pow(TG2[3],2));PUP[1][16]=TG2[1]; PUP[2][16]=TG2[2]; PUP[3][16]=TG2[3];}

	// Fill in Les Houches common block
	mass[1]=0.;
	mass[2]=0.;
	mass[3]=mt;
	mass[4]=mt;
	mass[5]=0.;
	mass[6]=0.;
	mass[7]=0.;
	mass[8]=0.;
	mass[9]=0.;
	mass[10]=0.;
	mass[11]=0.;
	mass[12]=0.;
	mass[13]=0.;
	mass[14]=0.;
	mass[17]=500.;
	mass[15]=0.;
	mass[16]=0.;
	PUP[4][13]=PUP[4][6];
	PUP[3][13]=PUP[3][6];
	PUP[2][13]=PUP[2][6];
	PUP[1][13]=PUP[1][6];
	PUP[4][14]=PUP[4][7];
	PUP[3][14]=PUP[3][7];
	PUP[2][14]=PUP[2][7];
	PUP[1][14]=PUP[1][7];
	PUP[4][1]=250.;
	PUP[3][1]=250.;
	PUP[2][1]=0.;
	PUP[1][1]=0.;
	PUP[4][2]=250.;
	PUP[3][2]=-250.;
	PUP[2][2]=0.;
	PUP[1][2]=0.;
	PUP[4][3]=TE1;
	PUP[3][3]=T1[3];
	PUP[2][3]=T1[2];
	PUP[1][3]=T1[1];
	PUP[4][4]=TE2;
	PUP[3][4]=T2[3];
	PUP[2][4]=T2[2];
	PUP[1][4]=T2[1];
        PUP[4][17]=emcm;
 	PUP[3][17]=0.;
 	PUP[2][17]=0.;
 	PUP[1][17]=0.;
	PUP[4][11]=GE1;
	PUP[3][11]=G1[3];
	PUP[2][11]=G1[2];
	PUP[1][11]=G1[1];
	PUP[4][12]=GE2;
	PUP[3][12]=G2[3];
	PUP[2][12]=G2[2];
	PUP[1][12]=G2[1];
	PUP[4][6]=BE1;
	PUP[3][6]=BB1[3];
	PUP[2][6]=BB1[2];
	PUP[1][6]=BB1[1];
	PUP[4][7]=BE2;
	PUP[3][7]=BB2[3];
	PUP[2][7]=BB2[2];
	PUP[1][7]=BB2[1];
	 PUP[4][9]=LE1;
	 PUP[3][9]=EE1[3];
	 PUP[2][9]=EE1[2];
	 PUP[1][9]=EE1[1];
	 PUP[4][8]=LE2;
	 PUP[3][8]=EE2[3];
	 PUP[2][8]=EE2[2];
	 PUP[1][8]=EE2[1];
	 PUP[4][10]=NE1;
	 PUP[3][10]=NN1[3];
	 PUP[2][10]=NN1[2];
	 PUP[1][10]=NN1[1];
	 PUP[4][5]=NE2;
	 PUP[3][5]=NN2[3];
	 PUP[2][5]=NN2[2];
	 PUP[1][5]=NN2[1];
	IDUP[1]=11;
	IDUP[2]=-11;
	IDUP[3]=6;
	IDUP[4]=-6;
	IDUP[11]=21;
	IDUP[6]=5;
	IDUP[7]=-5;
	IDUP[8]=11;
	IDUP[9]=-11;
	IDUP[10]=12;
	IDUP[5]=-12;
	IDUP[12]=21;
	IDUP[13]=21;
	  IDUP[14]=21;
	  IDUP[17]=22;
	  IDUP[15]=21;
	  IDUP[16]=21;
	MOTHUP[1][1]=0;
	MOTHUP[2][1]=0;
		 MOTHUP[1][2]=0;
		 MOTHUP[2][2]=0;
		 MOTHUP[1][3]=3;
		 MOTHUP[2][3]=3;
		 MOTHUP[1][4]=3;
		 MOTHUP[2][4]=3;
		 MOTHUP[1][17]=1;
		 MOTHUP[2][17]=2;
		 MOTHUP[1][11]=4;
		 MOTHUP[2][11]=4;
		 MOTHUP[1][6]=4;
		 MOTHUP[2][6]=4;
		 MOTHUP[1][7]=5;
		 MOTHUP[2][7]=5;
		 MOTHUP[1][8]=5;
		 MOTHUP[2][8]=5;
		 MOTHUP[1][9]=4;
		 MOTHUP[2][9]=4;
		 MOTHUP[1][10]=4;
		 MOTHUP[2][10]=4;
		 MOTHUP[1][5]=5;
		 MOTHUP[2][5]=5;
		 MOTHUP[1][12]=5;
		 MOTHUP[2][12]=5;
		 MOTHUP[1][13]=3;
		 MOTHUP[2][13]=3;
		 MOTHUP[1][14]=3;
		 MOTHUP[2][14]=3;
		 MOTHUP[1][15]=4;
		 MOTHUP[2][15]=4;
		 MOTHUP[1][16]=5;
		 MOTHUP[2][16]=5;
		 ISTUP[1]=-1;
		 ISTUP[2]=-1;
		 ISTUP[3]=2;
		 ISTUP[4]=2;
		 ISTUP[17]=2;
		 ISTUP[5]=1;
		 ISTUP[6]=1;
		 ISTUP[7]=1;
		 ISTUP[8]=1;
		 ISTUP[9]=1;
		 ISTUP[10]=1;
		 ISTUP[11]=1;
		 ISTUP[12]=1;
		 ISTUP[13]=1;
		 ISTUP[14]=1;
		 ISTUP[15]=1;
		 ISTUP[16]=1;
		 ISPINUP[1]=1;
		 ISPINUP[2]=-1;
		 ISPINUP[3]=1;
		 ISPINUP[4]=-1;
		 ISPINUP[5]=1;
		 ISPINUP[6]=1;
		 ISPINUP[7]=-1;
		 ISPINUP[8]=1;
		 ISPINUP[9]=-1;
		 ISPINUP[10]=1;
		 ISPINUP[11]=-1;
		 ISPINUP[12]=-1;
		 ISPINUP[13]=-1;
		 ISPINUP[14]=1;
		 ISPINUP[17]=1;
		 ISPINUP[15]=1;
		 ISPINUP[16]=1;
		 ICOLUP[1][1]=0;
		 ICOLUP[2][1]=0;
		 ICOLUP[1][2]=0;
		 ICOLUP[2][2]=0;
		 ICOLUP[1][17]=0;
		 ICOLUP[2][17]=0;
		 if (!POWprod && !POWdecay) {
		   ICOLUP[1][3]=501;
		   ICOLUP[2][3]=0;
		   ICOLUP[1][4]=0;
		   ICOLUP[2][4]=501;
		   ICOLUP[1][6]=501;
		   ICOLUP[2][6]=0;
		   ICOLUP[1][7]=0;
		   ICOLUP[2][7]=501;
		 }
		 if (!POWprod && POWdecay){
		 ICOLUP[1][3]=501;
		 ICOLUP[2][3]=0;
		 ICOLUP[1][4]=0;
		 ICOLUP[2][4]=501;
		 ICOLUP[1][11]=501;
		 ICOLUP[2][11]=503;
		 ICOLUP[1][6]=503;
		 ICOLUP[2][6]=0;
		 ICOLUP[1][7]=0;
		 ICOLUP[2][7]=504;
		 ICOLUP[1][12]=504;
		 ICOLUP[2][12]=501;
		 if (demitt) {
		 ICOLUP[1][11]=501;
		 ICOLUP[2][11]=505;
		 ICOLUP[1][15]=505;
		 ICOLUP[2][15]=503;}
		 if (demittb) {
		 ICOLUP[1][12]=506;
		 ICOLUP[2][12]=501;
		 ICOLUP[1][16]=504;
		 ICOLUP[2][16]=506;
		 }
		 }
		 if (POWprod && POWdecay) {
		 ICOLUP[1][3]=501;
		 ICOLUP[2][3]=0;
		 ICOLUP[1][4]=0;
		 ICOLUP[2][4]=502;
		 if (!emit){
	         ICOLUP[1][13]=502;
		 ICOLUP[2][13]=501;
		 ICOLUP[1][11]=501;
		 ICOLUP[2][11]=503;
		 ICOLUP[1][6]=503;
		 ICOLUP[2][6]=0;
		 ICOLUP[1][7]=0;
		 ICOLUP[2][7]=504;
		 ICOLUP[1][12]=504;
		 ICOLUP[2][12]=502;
		 if (demitt) {
		 ICOLUP[1][11]=505;
		 ICOLUP[2][11]=503;
		 ICOLUP[1][15]=501;
		 ICOLUP[2][15]=505;}
		 if (demittb) {
		 ICOLUP[1][12]=504;
		 ICOLUP[2][12]=506;
		 ICOLUP[1][16]=506;
		 ICOLUP[2][16]=502;
		 }
		 }else {
		 ICOLUP[1][13]=505;
		 ICOLUP[2][13]=501;
		 ICOLUP[1][14]=502;
		 ICOLUP[2][14]=505;
		 ICOLUP[1][11]=501;
		 ICOLUP[2][11]=503;
		 ICOLUP[1][6]=503;
		 ICOLUP[2][6]=0;
		 ICOLUP[1][7]=0;
		 ICOLUP[2][7]=504;
		 ICOLUP[1][12]=504;
		 ICOLUP[2][12]=502;
		 if (demitt) {
		 ICOLUP[1][11]=507;
		 ICOLUP[2][11]=503;
		 ICOLUP[1][15]=501;
		 ICOLUP[2][15]=507;}
		 if (demittb) {
		 ICOLUP[1][12]=504;
		 ICOLUP[2][12]=508;
		 ICOLUP[1][16]=508;
		 ICOLUP[2][16]=502;
		 }
}}
		 if (POWprod & !POWdecay) {
		   if (!emit){
		     ICOLUP[1][3]=501;
		   ICOLUP[2][3]=0;
		   ICOLUP[1][4]=0;
		   ICOLUP[2][4]=502;
	         ICOLUP[1][13]=502;
		 ICOLUP[2][13]=501;
		 ICOLUP[1][6]=501;
		 ICOLUP[2][6]=0;
		 ICOLUP[1][7]=0;
		 ICOLUP[2][7]=502;
		 } else {
		   ICOLUP[1][3]=501;
                   ICOLUP[2][3]=0;
                   ICOLUP[1][4]=0;
                   ICOLUP[2][4]=502;
		   ICOLUP[1][13]=505;
		   ICOLUP[2][13]=501;
		   ICOLUP[1][14]=502;
		   ICOLUP[2][14]=505;
		 ICOLUP[1][6]=501;
		 ICOLUP[2][6]=0;
		 ICOLUP[1][7]=0;
		 ICOLUP[2][7]=502;
		 }
		 }
		 ICOLUP[1][8]=0;
		 ICOLUP[2][8]=0;
		 ICOLUP[1][9]=0;
		 ICOLUP[2][9]=0;
		 ICOLUP[1][10]=0;
		 ICOLUP[2][10]=0;
		 ICOLUP[1][5]=0;
		 ICOLUP[2][5]=0;
		
		 

		    //*********************************************************************/
		    //MADGRAPH FORMAT
		    //********************************************************************/
		     
// 		     if (EVT==1) {count1++;}
// 		     if (EVT==2) {count2++;}
// 		     if (EVT==3) {count3++;}
// 		     if (EVT==4) {count4++;}
//                   	 NUP=11;
// 			 if (POWprod) {NUP++;
// 			 if (emit) {NUP++;}}
// 			 if (POWdecay) {NUP++; NUP++;
// 			 if (demitt) {NUP++;}
// 			 if (demittb) {NUP++;}}
// 			 //			 cout << "NUP" << NUP << endl;
		
// 			//  if (POWprod && POWdecay) {
// // 			   if (!emit) {NUP=13;} else {NUP=14;}}
// // 			  if (!POWprod && POWdecay) {NUP =12;}
// // 			 if (POWprod && !POWdecay) {if (!emit){NUP=11;}else{NUP=12;}}
// 			 XWGTUP=1;
// 			 SCALUP=500.;
// 			 AQEDUP=0.0073;
// 			 AQCDUP=0.118;
// 			 double sum1,sum2,sum3,sum4;
// 			 	 int NUPP=15;
// 			 	 if (!emit){NUPP=14;}
// 			 sum1=0.;sum2=0.;sum3=0.;sum4=0.;
// 			 	   outdata << NUP <<"\t"<<ixx<<"\t"<<XWGTUP<<"\t"<<SCALUP<<"\t"<<AQEDUP<<"\t"<<AQCDUP <<endl;
// 			 for (int ja = 1; ja < 11; ja++) {
// 			 outdata << IDUP[ja] << "\t" ;
// 			 if (ja==2)  {outdata << IDUP[17] << "\t" ;}}
// 			 if (POWdecay && POWprod) {{for (int ja=11; ja < 14; ja++)
// 			   {outdata << IDUP[ja] << "\t" ;}}
			
// 			 if (emit) {outdata << IDUP[14] << "\t" ;}
// 			  if (demitt) {outdata << IDUP[15] << "\t" ;}
// 			  if (demittb) {outdata << IDUP[16] << "\t" ;}}

// 			 //}
// 			 if (POWdecay && !POWprod) {{for (int ja=11; ja < 13; ja++)
// 			   {outdata << IDUP[ja] << "\t" ;}}
// 			  if (demitt) {
// 			 for (int ja=15; ja < NUP+2; ja++)
// 			   {outdata << IDUP[ja] << "\t" ;}}
// 			   if (!demitt && demittb) 
// 			    {outdata << IDUP[16] << "\t" ;} 
// 			 }
// 			 if (POWprod && !POWdecay) {for (int ja=13; ja<NUPP;ja++)
// 			 {outdata << IDUP[ja] << "\t" ;}}
// 			 outdata << endl;

// 			 for (int jb = 1; jb < 3; jb++) {
// 			 for (int jc = 1; jc < 11; jc++) {
// 			 outdata << MOTHUP[jb][jc]<<"\t" ;
// 			 if (jc==2)  {outdata << MOTHUP[jb][17] << "\t" ;}}
// 			 if (POWdecay && POWprod) {{for (int jc=11; jc < 14; jc++)
// 			   {outdata << MOTHUP[jb][jc] << "\t" ;} }
		
// 			 if (emit) {outdata << MOTHUP[jb][14] << "\t" ;}
// 			  if (demitt) {outdata << MOTHUP[jb][15] << "\t" ;}
// 			  if (demittb) {outdata << MOTHUP[jb][16] << "\t" ;}}



// 			  if (POWdecay && !POWprod) {for (int jc=11; jc < 13; jc++)
// 			    {outdata << MOTHUP[jb][jc] << "\t" ;}
// 			  if (demitt) {
// 			 for (int jc=15; jc < NUP+2; jc++)
// 			   {outdata << MOTHUP[jb][jc] << "\t" ;} }
// 			  if (!demitt && demittb) 
// 			    {outdata << MOTHUP[jb][16] << "\t" ;} 
// 			 }
// 			 if (POWprod && !POWdecay) {for (int jc=13; jc<NUPP;jc++)
// 			 {outdata << MOTHUP[jb][jc] << "\t" ;}}
// 			 outdata << endl;}
// 			 for (int jd = 1; jd < 3; jd++) {
// 			    for (int je = 1; je < 11; je++) {
// 			  outdata << ICOLUP[jd][je]<<"\t" ;
// 			  if (je==2)  {outdata << ICOLUP[jd][17] << "\t" ;}}
// 			 if (POWdecay && POWprod) {{for (int je=11; je < 14; je++)
// 			   {outdata << ICOLUP[jd][je] << "\t" ;}}
// 			 if (emit) {outdata << ICOLUP[jd][14] << "\t" ;}
// 			  if (demitt) {outdata << ICOLUP[jd][15] << "\t" ;}
// 			  if (demittb) {outdata << ICOLUP[jd][16] << "\t" ;}
// }
// 			 if (POWdecay && !POWprod) {for (int je=11; je < 13; je++)
// 			    {outdata << ICOLUP[jd][je] << "\t" ;}
// 			   if (demitt) {
// 			 for (int je=15; je < NUP+2; je++)
// 			   {outdata << ICOLUP[jd][je] << "\t" ;} }
// 			    if (!demitt && demittb) 
// 			    {outdata << ICOLUP[jd][16] << "\t" ;}  
// 			 }
// 			 if (POWprod && !POWdecay) {for (int je=13; je<NUPP;je++)
// 			 {outdata << ICOLUP[jd][je] << "\t" ;}}
// 			  outdata << endl;}
// 			 for (int jf = 1; jf < 11; jf++) {
// 			 outdata << ISTUP[jf] << "\t" ;
// 			 if (jf==2) { outdata << ISTUP[17] << "\t" ;}}
// 			 if (POWdecay && POWprod) {{for (int jf=11; jf< 14; jf++)
// 			   {outdata << ISTUP[jf] << "\t" ;}}
// 			 if (emit) {outdata << ISTUP[14] << "\t" ;}
// 			  if (demitt) {outdata << ISTUP[15] << "\t" ;}
// 			  if (demittb) {outdata << ISTUP[16] << "\t" ;}
// }
// 			 if (POWdecay && !POWprod) {{for (int jf=11; jf < 13; jf++)
// 			   {outdata << ISTUP[jf] << "\t" ;}}
// 			   if (demitt) {
// 			 for (int jf=15; jf < NUP+2; jf++)
// 			   {outdata << ISTUP[jf] << "\t" ;} }
// 			   if (!demitt && demittb) 
// 			    {outdata << ISTUP[16] << "\t" ;}  
// 			 }
// 			  if (POWprod && !POWdecay) {for (int jf=13; jf<NUPP;jf++)
// 			 {outdata << ISTUP[jf] << "\t" ;}}
// 			  outdata << endl;
// 			 for (int jg = 1; jg < 11; jg++) {
// 			  outdata << ISPINUP[jg] << "\t" ;
// 			  if (jg==2) { outdata << ISPINUP[17] << "\t" ;}}
// 			 if (POWdecay && POWprod) {{for (int jg=11; jg< 14; jg++)
// 			   {outdata << ISPINUP[jg] << "\t" ;}}
// 			  if (emit) {outdata << ISPINUP[14] << "\t" ;}
// 			  if (demitt) {outdata << ISPINUP[15] << "\t" ;}
// 			  if (demittb) {outdata << ISPINUP[16] << "\t" ;}
// }
// 			 if (POWdecay && !POWprod) {for (int jg=11; jg < 13; jg++)
// 			    {outdata << ISPINUP[jg] << "\t" ;}
// 			  if (demitt) {
// 			 for (int jg=15; jg < NUP+2; jg++)
// 			   {outdata << ISPINUP[jg] << "\t" ;} }
// 			    if (!demitt && demittb) 
// 			    {outdata << ISPINUP[16] << "\t" ;}  
// 			 }
// 			  if (POWprod && !POWdecay) {for (int jg=13; jg<NUPP;jg++)
// 			  {outdata << ISPINUP[jg] << "\t" ;}}
			  
// 			   outdata << endl;
// 			   int jhh=0;
//  			 for (int jh = 1; jh < 11; jh++) {
// 			   jhh=jhh+1;
// 			    outdata << jhh <<"\t" << PUP[4][jh] << "\t" << PUP[1][jh]<<"\t"<<PUP[2][jh]<<"\t"<<PUP[3][jh]<<endl;
// 			  if(jh==2){outdata << 3 <<"\t" << PUP[4][17] << "\t" << PUP[1][17]<<"\t"<<PUP[2][17]<<"\t"<<PUP[3][17]<<endl;jhh=jhh+1;}
// 			  if (jh != 3 && jh != 4){
// 			    sum1+=PUP[4][jh];sum2+=PUP[1][jh];sum3+=PUP[2][jh];sum3+=PUP[3][jh];}}
// 			   if (POWdecay && POWprod) {{for (int jh=11; jh< 14; jh++)
			     
// 			   {jhh=jh+1;outdata << jhh <<"\t" << PUP[4][jh] << "\t" << PUP[1][jh]<<"\t"<<PUP[2][jh]<<"\t"<<PUP[3][jh]<<endl;
			    
// 			   sum1+=PUP[4][jh];sum2+=PUP[1][jh];sum3+=PUP[2][jh];sum3+=PUP[3][jh];}}
// 			   if (emit) {
// 			     outdata << jhh+1 <<"\t" << PUP[4][14] << "\t" << PUP[1][14]<<"\t"<<PUP[2][14]<<"\t"<<PUP[3][14]<<endl;
// 			     sum1+=PUP[4][14];sum2+=PUP[1][14];sum3+=PUP[2][14];sum3+=PUP[3][14];
// }
// 			   if (demitt) {
// 			     outdata << jhh+1 <<"\t" << PUP[4][15] << "\t" << PUP[1][15]<<"\t"<<PUP[2][15]<<"\t"<<PUP[3][15]<<endl;
// 			     sum1+=PUP[4][15];sum2+=PUP[1][15];sum3+=PUP[2][15];sum3+=PUP[3][15];
// }
// 			   	   if (demittb) {
// 			     outdata << jhh+1 <<"\t" << PUP[4][16] << "\t" << PUP[1][16]<<"\t"<<PUP[2][16]<<"\t"<<PUP[3][16]<<endl;
// 			     sum1+=PUP[4][16];sum2+=PUP[1][16];sum3+=PUP[2][16];sum3+=PUP[3][16];
// }
// }
// 			    if (POWdecay && !POWprod) {for (int jh=11; jh< 13; jh++)
			     
// 			   {jhh=jh+1;outdata << jhh <<"\t" << PUP[4][jh] << "\t" << PUP[1][jh]<<"\t"<<PUP[2][jh]<<"\t"<<PUP[3][jh]<<endl;
// 			 sum1+=PUP[4][jh];sum2+=PUP[1][jh];sum3+=PUP[2][jh];sum3+=PUP[3][jh];}
// 			    if (demitt) {
// 			    for (int jh=15; jh< NUP+2; jh++)
			     
// 			   {jhh=jh+2;outdata << jhh-3 <<"\t" << PUP[4][jh] << "\t" << PUP[1][jh]<<"\t"<<PUP[2][jh]<<"\t"<<PUP[3][jh]<<endl;
// 			   sum1+=PUP[4][jh];sum2+=PUP[1][jh];sum3+=PUP[2][jh];sum3+=PUP[3][jh];} }
// 			    if (!demitt && demittb) {
// 			   outdata << 14 <<"\t" << PUP[4][16] << "\t" << PUP[1][16]<<"\t"<<PUP[2][16]<<"\t"<<PUP[3][16]<<endl;
// 			   sum1+=PUP[4][16];sum2+=PUP[1][16];sum3+=PUP[2][16];sum3+=PUP[3][16];} 
// 			    }



// 			 	 if (POWprod && !POWdecay) {for (int jh=13; jh<NUPP;jh++)
// 			   {jhh=jh+1;outdata << jhh-2 <<"\t" << PUP[4][jh] << "\t" << PUP[1][jh]<<"\t"<<PUP[2][jh]<<"\t"<<PUP[3][jh]<<endl;
			   
// 			         sum1+=PUP[4][jh];sum2+=PUP[1][jh];sum3+=PUP[2][jh];sum3+=PUP[3][jh];}}
// 				 //	 cout << count1 << "\t" << count2 << "\t" << count3 << "\t" << count4 << endl;
// 				 cout <<ixx  << "\t" << sum1 << "\t" <<sum2 <<"\t" <<sum3 <<"\t" << sum4 <<  endl;
//   }



		    //*********************************************************************/
		    //XML FORMAT
		    //********************************************************************/

                         NUP=11;
			 if (POWprod) {NUP++;
			 if (emit) {NUP++;}}
			 if (POWdecay) {NUP++; NUP++;
			 if (demitt) {NUP++;}
			 if (demittb) {NUP++;}}		      
// 			 NUP=11;
// 			 if (POWprod && POWdecay) {
// 			   if (!emit) {NUP=14;} else {NUP=15;}}
// 			  if (!POWprod && POWdecay) {NUP =13;}
// 			 if (POWprod && !POWdecay) {if (!emit){NUP=12;}else{NUP=13;}}
			 XWGTUP=1;
			 SCALUP=500.;
			 AQEDUP=0.0073;
			 AQCDUP=0.118;
			 double sum1,sum2,sum3,sum4;
			 int NUPP=15;
			 if (!emit){NUPP=14;}
			 sum1=0.;sum2=0.;sum3=0.;sum4=0.;
			 /************************************************************************************************/
			 if (xml) {
			 outdata << "<event>" << endl;
			 outdata << NUP <<"\t" << "11" <<"\t"<<XWGTUP<<"\t" << SCALUP << "\t" << AQEDUP << "\t" << AQCDUP << endl;
		
			 for (int ja = 1; ja < 11; ja++) {
			   if (ja == 3) {
			   outdata << IDUP[17] << "\t" << ISTUP[17] << "\t" << MOTHUP[1][17] <<"\t" << MOTHUP[2][17] << "\t" << ICOLUP[1][17] << "\t" << ICOLUP[2][17] <<"\t" <<setprecision (9)<< PUP[1][17] << "\t" << PUP[2][17]<<"\t"<<PUP[3][17]<<"\t"<<PUP[4][17]<<"\t" << mass[17] <<"\t"<< "0" << "\t" <<"9" <<endl;   
			   }
			   outdata << IDUP[ja] << "\t" << ISTUP[ja] << "\t" << MOTHUP[1][ja] <<"\t" << MOTHUP[2][ja] << "\t" << ICOLUP[1][ja] << "\t" << ICOLUP[2][ja] <<"\t" <<setprecision (9)<< PUP[1][ja] << "\t" << PUP[2][ja]<<"\t"<<PUP[3][ja]<<"\t"<<PUP[4][ja]<<"\t" << mass[ja] <<"\t"<< "0" << "\t" <<"9" <<endl;
			   if (ja !=3 && ja != 4){
			     sum1+=PUP[4][ja];sum2+=PUP[1][ja];sum3+=PUP[2][ja];sum3+=PUP[3][ja];}
		      }
			 if (POWdecay && POWprod) {{for (int ja=11; ja < 14; ja++)
			   {outdata << IDUP[ja] << "\t" << ISTUP[ja] << "\t" << MOTHUP[1][ja] <<"\t" << MOTHUP[2][ja] << "\t" << ICOLUP[1][ja] << "\t" << ICOLUP[2][ja] <<"\t" <<setprecision (9)<< PUP[1][ja] << "\t" << PUP[2][ja]<<"\t"<<PUP[3][ja]<<"\t"<<PUP[4][ja]<<"\t" << mass[ja] <<"\t"<< "0" << "\t" <<"9" <<endl;
			   sum1+=PUP[4][ja];sum2+=PUP[1][ja];sum3+=PUP[2][ja];sum3+=PUP[3][ja];
}}
			 if (emit) {outdata << IDUP[14] << "\t" << ISTUP[14] << "\t" << MOTHUP[1][14] <<"\t" << MOTHUP[2][14] << "\t" << ICOLUP[1][14] << "\t" << ICOLUP[2][14] <<"\t" <<setprecision (9)<< PUP[1][14] << "\t" << PUP[2][14]<<"\t"<<PUP[3][14]<<"\t"<<PUP[4][14]<<"\t" << mass[14] <<"\t"<< "0" << "\t" <<"9" <<endl;
			 sum1+=PUP[4][14];sum2+=PUP[1][14];sum3+=PUP[2][14];sum3+=PUP[3][14];
}
			  if (demitt) {outdata << IDUP[15] << "\t" << ISTUP[15] << "\t" << MOTHUP[1][15] <<"\t" << MOTHUP[2][15] << "\t" << ICOLUP[1][15] << "\t" << ICOLUP[2][15] <<"\t" <<setprecision (9)<< PUP[1][15] << "\t" << PUP[2][15]<<"\t"<<PUP[3][15]<<"\t"<<PUP[4][15]<<"\t" << mass[15] <<"\t"<< "0" << "\t" <<"9" <<endl;
			  sum1+=PUP[4][15];sum2+=PUP[1][15];sum3+=PUP[2][15];sum3+=PUP[3][15];
}
			  if (demittb) {outdata << IDUP[16] << "\t" << ISTUP[16] << "\t" << MOTHUP[1][16] <<"\t" << MOTHUP[2][16] << "\t" << ICOLUP[1][16] << "\t" << ICOLUP[2][16] <<"\t" <<setprecision (9)<< PUP[1][16] << "\t" << PUP[2][16]<<"\t"<<PUP[3][16]<<"\t"<<PUP[4][16]<<"\t" << mass[16] <<"\t"<< "0" << "\t" <<"9" <<endl;
			  sum1+=PUP[4][16];sum2+=PUP[1][16];sum3+=PUP[2][16];sum3+=PUP[3][16];
}}

			 if (POWdecay && !POWprod) {{for (int ja=11; ja < 13; ja++)
			   {outdata << IDUP[ja] << "\t" << ISTUP[ja] << "\t" << MOTHUP[1][ja] <<"\t" << MOTHUP[2][ja] << "\t" << ICOLUP[1][ja] << "\t" << ICOLUP[2][ja] <<"\t" <<setprecision (9)<< PUP[1][ja] << "\t" << PUP[2][ja]<<"\t"<<PUP[3][ja]<<"\t"<<PUP[4][ja]<<"\t" << mass[ja] <<"\t"<< "0" << "\t" <<"9" <<endl;
			   sum1+=PUP[4][ja];sum2+=PUP[1][ja];sum3+=PUP[2][ja];sum3+=PUP[3][ja];
}}
			  if (demitt) {
			 for (int ja=15; ja < NUP+2; ja++)
			   {outdata << IDUP[ja] << "\t" << ISTUP[ja] << "\t" << MOTHUP[1][ja] <<"\t" << MOTHUP[2][ja] << "\t" << ICOLUP[1][ja] << "\t" << ICOLUP[2][ja] <<"\t" <<setprecision (9)<< PUP[1][ja] << "\t" << PUP[2][ja]<<"\t"<<PUP[3][ja]<<"\t"<<PUP[4][ja]<<"\t" << mass[ja] <<"\t"<< "0" << "\t" <<"9" <<endl;
			   sum1+=PUP[4][ja];sum2+=PUP[1][ja];sum3+=PUP[2][ja];sum3+=PUP[3][ja];
}}
			   if (!demitt && demittb) 
			    {outdata << IDUP[16] << "\t" << ISTUP[16] << "\t" << MOTHUP[1][16] <<"\t" << MOTHUP[2][16] << "\t" << ICOLUP[1][16] << "\t" << ICOLUP[2][16] <<"\t" <<setprecision (9)<< PUP[1][16] << "\t" << PUP[2][16]<<"\t"<<PUP[3][16]<<"\t"<<PUP[4][16]<<"\t" << mass[16] <<"\t"<< "0" << "\t" <<"9" <<endl;
			    sum1+=PUP[4][16];sum2+=PUP[1][16];sum3+=PUP[2][16];sum3+=PUP[3][16];
} 
			 }
			 if (POWprod && !POWdecay) {for (int ja=13; ja<NUPP;ja++)
			 {outdata << IDUP[ja] << "\t" << ISTUP[ja] << "\t" << MOTHUP[1][ja] <<"\t" << MOTHUP[2][ja] << "\t" << ICOLUP[1][ja] << "\t" << ICOLUP[2][ja] <<"\t" <<setprecision (9)<< PUP[1][ja] << "\t" << PUP[2][ja]<<"\t"<<PUP[3][ja]<<"\t"<<PUP[4][ja]<<"\t" << mass[ja] <<"\t"<< "0" << "\t" <<"9" <<endl;
			 sum1+=PUP[4][ja];sum2+=PUP[1][ja];sum3+=PUP[2][ja];sum3+=PUP[3][ja];
}}
			 //	 cout <<ixx  << "\t" << sum1 << "\t" <<sum2 <<"\t" <<sum3 <<"\t" << sum4 <<  endl;
			 cout  << "Generated event : " << "\t" << ixx << "\t" << "of" << "\t" << nevg << "\r" << flush;  
			 outdata << "</event>" << endl;}
                       }
    outdata << "</LesHouchesEvents>" << endl;
  cout << endl;
   return 0;
}
/***************************************************************************************************/
/* Axial or Vector axial current ? */
 bool ax(int seed1) {
   bool axial=false;
       if (random(seed1) < A/(A+V)){axial = true;}
       return axial;    
}

/***************************************************************************************************/
/* Generates interpolation files for production */

  void filegen (double &rintgend1a, double &rintgend2a,double &rintgend3a, double &rintgend4a,double &rintgend1v, double &rintgend2v,double &rintgend3v, double &rintgend4v){
    cout << "Generating production interpolation files..." << endl;
  ofstream outdata,outdata1a,outdata2a,outdata3a,outdata4a,outdata1v,outdata2v,outdata3v,outdata4v;
  outdata1a.open("ttbaraxpro1.dat", ios::trunc); 
  outdata2a.open("ttbaraxpro2.dat", ios::trunc); 
  outdata3a.open("ttbaraxpro3.dat", ios::trunc); 
  outdata4a.open("ttbaraxpro4.dat", ios::trunc); 
  outdata1v.open("ttbarvecpro1.dat", ios::trunc); 
  outdata2v.open("ttbarvecpro2.dat", ios::trunc); 
  outdata3v.open("ttbarvecpro3.dat", ios::trunc); 
  outdata4v.open("ttbarvecpro4.dat", ios::trunc); 
  double pt2t = 0.; double pt2b = 0.;
   double intga,intgv, pt2;
   int ngen;
   intga = 0.;
   intgv = 0.;
   ngen = 1000;
   for (int iy=1; iy < 5; iy++) {

     if (iy==1) {pt2t = pt2max; pt2b=0.001;}
     if (iy==2) {pt2t = 0.001; pt2b=0.0001;}
     if (iy==3) {pt2t = 0.0001; pt2b=0.00001;}
     if (iy==4) {pt2t = 0.00001; pt2b=pt2min;}
     double div = (pt2t-pt2b)/ngen;  
     pt2 = pt2t+div;
   for (int ix=0; ix<ngen; ix++) {
   pt2=pt2-div;
   if (pt2 < pt2b) {break;}
    double xmax=1.-2.*pt2-2.*sqrt(pow(pt2,2)+rho*pt2);
   double alphas=1./(b*log(pt2/pt2min))*(1.-bp*log(log(pt2/pt2min))/(b*log(pt2/pt2min)));
   alphas+=pow(alphas,2)*(67./3.-pow(pi,2)-50./9.)/(4.*pi);
   double xminr=39.*pt2-1.+3.*d2-42.*pt2*d2-3.*pow(d2,2)+pow(pt2,3)+15.*pow(pt2,2)+3.*pow(pt2,2)*d2+3.*pt2*pow(d2,2)+pow(d2,3);
   double xmini=6.*sqrt (abs(-3.*pt2+33.*pow(pt2,2)+3.*pow(pt2,3)-18.*pt2*pow(d2,2)+12.*pt2*pow(d2,3)+12.*pt2*d2-42.*pow(pt2,3)*d2+51.*pow(pt2,2)*pow(d2,2)-75.*pow(pt2,2)*d2-3.*d2*pow(pt2,4)-3.*pt2*pow(d2,4)-9.*pow(pt2,3)*pow(d2,2)-9.*pow(pt2,2)*pow(d2,3)));    
   double xminm=sqrt(pow(xminr,2)+pow(xmini,2));
   double xmint=atan2(xmini,xminr);
   double xmin1=-1./12.*pow(xminm,(1./3.))*cos(xmint/3.);
   double xmin2=(-1./12.+d2/6.-pow(pt2,2)/12.-5.*pt2/6.-pt2*d2/6.-pow(d2,2)/12.)*cos(xmint/3.)/(pow(xminm,(1./3.)));
   double xmin3=pt2/6.+5./6.+d2/6.;
   double xmin4=sqrt(3.)/12.*pow(xminm,(1./3.))*sin(xmint/3.)
               -3.*sqrt(3.)*(-1./36.+d2/18.-pow(pt2,2)/36.-5.*pt2/18.-pt2/18.*d2-pow(d2,2)/36.)
               /pow(xminm,(1./3.))*sin(xmint/3.);
   double xmin5=xmin1+xmin2+xmin3+xmin4;
   double xmin6=xmin1+xmin2+xmin3-xmin4;
   
   if (pt2 > 0.024) {xmax=xmin5;}
    double div1 = (xmax - xmin6)/ngen;
    double div2 = (xmax - xmin5)/ngen;
     double xx ;
    xx=xmax +div1/2.; 
    for (int ix =0; ix < ngen; ix++) {
      xx=xx-div1;
      if (xx < xmin6) {break;}
      double yy = 1./(2.*(1.+rho-xx))*(pow(xx,2)-3.*xx-2.*rho*xx+2.+4.*rho-sqrt((pow(xx,2)-4.*rho)*(4.*pt2*(xx-1.-rho)+pow((xx-1.),2))));
 double dydkap = sqrt(pow(xx,2)-4.*rho)/sqrt(4.*pt2*(xx-1.-rho)+pow(xx-1.,2));
     

      double Waa = ((pow((xx+0.5*d2),2)+pow((yy+0.5*d2),2)+0.5*d2*(pow((5.-xx-yy),2)-19.+d2))/((1.-d2)*(1.-xx)*(1.-yy))-d2/
		   (2.*pow((1.-xx),2))-d2/(2.*pow((1.-yy),2)));
      double Wvv = ((pow((xx+0.5*d2),2)+pow((yy+0.5*d2),2)-2.*d2*(1.+d2/2.))/((1.+d2/2.)*(1.-xx)*(1.-yy))-d2/(2.*pow((1.-
														      xx),2))-d2/(2.*pow((1.-yy),2)));
           if (yy*0. != 0.) {dydkap=0.; Waa=0.; Wvv=0.;} 
      intga +=  dydkap*Waa*2.*alphas*div1/(3.*pi*sqrt(1.-d2))*div;
      intgv +=  dydkap*Wvv*2.*alphas*div1/(3.*pi*sqrt(1.-d2))*div;
     


}
    xx =xmax + div2/2.; 
    if (pt2 < 0.024) {
    for (int iz =0; iz < ngen; iz++) {
      xx=xx-div2;
double yy = 1./(2.*(1.+rho-xx))*(pow(xx,2)-3.*xx-2.*rho*xx+2.+4.*rho+sqrt((pow(xx,2)-4.*rho)*(4.*pt2*(xx-1.-rho)+pow((xx-1.),2))));


 if (yy*0. != 0.) {continue;} 
      if (xx < xmin5) {break;}
   double  dydkap = sqrt(pow(xx,2)-4.*rho)/sqrt(4.*pt2*(xx-1.-rho)+pow(xx-1.,2));
   double  Waa = ((pow((xx+0.5*d2),2)+pow((yy+0.5*d2),2)+0.5*d2*(pow((5.-xx-yy),2)-19.+d2))/((1.-d2)*(1.-xx)*(1.-yy))-d2/(2.*pow((1.-xx),2))-d2/(2.*pow((1.-yy),2)));
  double Wvv = ((pow((xx+0.5*d2),2)+pow((yy+0.5*d2),2)-2.*d2*(1.+d2/2.))/((1.+d2/2.)*(1.-xx)*(1.-yy))-d2/(2.*pow((1.-xx),2))-d2/(2.*pow((1.-yy),2)));
  if (yy*0. != 0.) {dydkap=0.;Waa=0.; Wvv=0.;}
      intga +=  dydkap*Waa*2.*alphas*div2/(3.*pi*sqrt(1.-d2))*div;
      intgv +=  dydkap*Wvv*2.*alphas*div2/(3.*pi*sqrt(1.-d2))*div;
     
    }}
   
 
   if (iy==1) {outdata1a  << pt2 << "\t" << abs(intga)  <<endl;
               outdata1v  << pt2 << "\t" << abs(intgv)  <<endl;
}
   if (iy==2) {outdata2a  << pt2 << "\t" << abs(intga)  <<endl;
               outdata2v  << pt2 << "\t" << abs(intgv)  <<endl;
}
   if (iy==3) {outdata3a  << pt2 << "\t" << abs(intga)  <<endl;
               outdata3v  << pt2 << "\t" << abs(intgv)  <<endl;
}
   if (iy==4) {outdata4a  << pt2 << "\t" << abs(intga)  <<endl;
               outdata4v  << pt2 << "\t" << abs(intgv)  <<endl;
}
   }
   
   if (iy==1) {rintgend1a=abs(intga);rintgend1v=abs(intgv);}
   if (iy==2) {rintgend2a=abs(intga);rintgend2v=abs(intgv);}
   if (iy==3) {rintgend3a=abs(intga);rintgend3v=abs(intgv);}
   if (iy==4) {rintgend4a=abs(intga);rintgend4v=abs(intgv);}
}
    cout << "...Done" << endl;
}
/***************************************************************************************************/
/* Generates interpolation files for decays */

void filegendec (double &rintgend1, double &rintgend2, double &rintgend3, double &rintgend4){
   cout << "Generating decay interpolation files..." << endl;
  ofstream outdata1d,outdata2d,outdata3d,outdata4d;
  outdata1d.open("ttbardecay1.dat", ios::trunc); 
  outdata2d.open("ttbardecay2.dat", ios::trunc); 
  outdata3d.open("ttbardecay3.dat", ios::trunc); 
  outdata4d.open("ttbardecay4.dat", ios::trunc); 
   double pt2t= 0.; double pt2b= 0.;
   double intg, pt2;
   int ngen;
   intg = 0.;
   ngen = 1000;
   for (int iy=1; iy < 5; iy++) {
     if (iy==1) {pt2t = pt2maxd; pt2b=0.001*pow(mt,2);}
     if (iy==2) {pt2t = 0.001*pow(mt,2); pt2b=0.0001*pow(mt,2);}
     if (iy==3) {pt2t = 0.0001*pow(mt,2); pt2b=0.00001*pow(mt,2);}
     if (iy==4) {pt2t = 0.00001*pow(mt,2); pt2b=pt2mind;}
     double div = (pt2t-pt2b)/ngen;  
     pt2 = pt2t+div/2.;
 for (int ix=0;ix<ngen;ix++) {
 pt2=pt2-div;
 if (pt2 < pt2b) {break;}
  double alphas=1./(b*log(pt2/0.04))*(1.-bp*log(log(pt2/0.04))/(b*log(pt2/0.04)));
   alphas+=pow(alphas,2)*(67./3.-pow(pi,2)-50./9.)/(4.*pi);
  double xmax=(mt-a*mt-2.*sqrt(a*pt2))/mt;
  double xmin=2.*sqrt(pt2)/mt;
 
    double xx ;
      double div1 = (xmax - xmin)/ngen;
	 xx=xmax +div1/2.;
	   for (int ix =0; ix < ngen; ix++) {
      xx=xx-div1;
      if (xx < xmin) {break;} 
      double kapp=pt2/pow(mt,2);
              double yy1 = (pow(xx,2)+a*xx+2.-3.*xx-2.*a*kapp+sqrt((pow(xx,2)-4.*kapp*(1.-a))*pow(xx-1.,2)+4.*a*kapp*(4.*kapp+1.-a)+pow(xx,2)*(a+2.*xx-2.-8.*kapp)*a))/(2.*(kapp+1.-xx));
           double yy2 = (pow(xx,2)+a*xx+2.-3.*xx-2.*a*kapp-sqrt((pow(xx,2)-4.*kapp*(1.-a))*pow(xx-1.,2)+4.*a*kapp*(4.*kapp+1.-a)+pow(xx,2)*(a+2.*xx-2.-8.*kapp)*a))/(2.*(kapp+1.-xx));
	   double gamm1 =  alphas*4./(pi*3.*(1.-yy1)*pow(xx,2))*(xx-((1.-yy1)*(1.-xx)+pow(xx,2))/(1.-a)+xx*pow(yy1+xx-1.,2)/(2.*pow(1.-a,2))+2.*a*(1.-yy1)*pow(xx,2)/(pow(1.-a,2)*(1.+2.*a)));
	   double gamm2 =  alphas*4./(pi*3.*(1.-yy2)*pow(xx,2))*(xx-((1.-yy2)*(1.-xx)+pow(xx,2))/(1.-a)+xx*pow(yy2+xx-1.,2)/(2.*pow(1.-a,2))+2.*a*(1.-yy2)*pow(xx,2)/(pow(1.-a,2)*(1.+2.*a)));
	       double dydkapp1=(-2.*a+(0.5*(-4.*pow(xx,2)+8.*xx-4.-4.*pow(xx,2)*a-8.*a*xx+32.*a*kapp+8.*a-4.*pow(a,2)))/sqrt((pow(xx,2)-4.*kapp*(1.-a))*pow(xx-1.,2)+4.*a*kapp*(4.*kapp+1.-a)+pow(xx,2)*(a+2.*xx-2.-8.*kapp)*a)-2.*yy1)/(2.*kapp+2.-2.*xx);
      double dydkapp2=(-2.*a-(0.5*(-4.*pow(xx,2)+8.*xx-4.-4.*pow(xx,2)*a-8.*a*xx+32.*a*kapp+8.*a-4.*pow(a,2)))/sqrt((pow(xx,2)-4.*kapp*(1.-a))*pow(xx-1.,2)+4.*a*kapp*(4.*kapp+1.-a)+pow(xx,2)*(a+2.*xx-2.-8.*kapp)*a)-2.*yy2)/(2.*kapp+2.-2.*xx);
      intg+=(gamm1*fabs(dydkapp1)+gamm2*fabs(dydkapp2))*div1*div/pow(mt,2);}

   if (iy==1) {outdata1d  << pt2/pow(mt,2) << "\t" << abs(intg) << endl;}
   if (iy==2) {outdata2d  << pt2/pow(mt,2) << "\t" << abs(intg) << endl;}
   if (iy==3) {outdata3d  << pt2/pow(mt,2) << "\t" << abs(intg) << endl;}
   if (iy==4) {outdata4d  << pt2/pow(mt,2) << "\t" << abs(intg) << endl;}
   }
   if (iy==1) {rintgend1=abs(intg);}
   if (iy==2) {rintgend2=abs(intg);}
   if (iy==3) {rintgend3=abs(intg);}
   if (iy==4) {rintgend4=abs(intg);}
   }
   cout << "...Done" << endl;
}

/***************************************************************************************************/
// Interpolation for production process.

double interp(int N,double XXX, double Xmax1, double Xmax2, double Xmax3, double Xmax4, bool axi) {
       double FF;
       double AA;
  int IX=0;
  int IY=N+1;
  int MID;
  double aa,bb,cc,dd,ee,ii,gg,hh;
  ifstream indata;

      for (int iy=0;(IY-IX) > 1;iy++) {
	 if (axi) {
	if (XXX < Xmax1 ) {
    indata.open("ttbaraxpro1.dat", ios::in); 
  }else if (XXX > Xmax1  && XXX < Xmax2 ) {
    indata.open("ttbaraxpro2.dat", ios::in);
  }else if (XXX > Xmax2 && XXX < Xmax3 ) {
    indata.open("ttbaraxpro3.dat", ios::in);
  }else if (XXX > Xmax3) {
    indata.open("ttbaraxpro4.dat", ios::in);
  }}  else {
    if (XXX < Xmax1 ) {
    indata.open("ttbarvecpro1.dat", ios::in); 
  }else if (XXX > Xmax1  && XXX < Xmax2 ) {
    indata.open("ttbarvecpro2.dat", ios::in);
  }else if (XXX > Xmax2 && XXX < Xmax3 ) {
    indata.open("ttbarvecpro3.dat", ios::in);
  }else if (XXX > Xmax3) {
    indata.open("ttbarvecpro4.dat", ios::in);
  }} 

        MID=(IX+IY)/2;
        for(int ix=0;ix < MID ;ix++) {
	  indata >> FF >> AA;}
	
	 indata.close();
	
        if (XXX >= AA) {
        IX=MID;} else {IY=MID;}
	if ((IY-IX) < 1) break;}
          if (IX < 5) {IX=5;}
    if (axi) {
      if (XXX < Xmax1 ) {
     indata.open("ttbaraxpro1.dat", ios::in); 
  }else if (XXX > Xmax1  && XXX < Xmax2 ) {
    indata.open("ttbaraxpro2.dat", ios::in);
  }else if (XXX > Xmax2 && XXX < Xmax3 ) {
    indata.open("ttbaraxpro3.dat", ios::in);
  }else if (XXX > Xmax3) {
    indata.open("ttbaraxpro4.dat", ios::in);
  } } else {
       if (XXX < Xmax1 ) {
  indata.open("ttbarvecpro1.dat", ios::in); 
  }else if (XXX > Xmax1  && XXX < Xmax2 ) {
    indata.open("ttbarvecpro2.dat", ios::in);
  }else if (XXX > Xmax2 && XXX < Xmax3 ) {
    indata.open("ttbarvecpro3.dat", ios::in);
  }else if (XXX > Xmax3) {
    indata.open("ttbarvecpro4.dat", ios::in);
  } 
  }
     for (int ix=0;ix<IX-3;ix++) {
        indata >> FF >> AA;}
        indata >> aa >> bb ;
        indata >> cc >> dd ;
        indata >> ee >> ii ;
        indata >> gg >> hh ;

        indata.close();
       double d1=(cc-aa)/(dd-bb);
       double d2=(ee-cc)/(ii-dd);
       double d3=(gg-ee)/(hh-ii) ;
       double dydx=(gg-aa)/(hh-bb);
       double dydx1=(ee-aa)/(ii-bb);
	 double dydx2=(gg-cc)/(hh-dd);
	 double d2ydx2=(dydx2-dydx1)/((hh+dd-ii-bb)/2.);
	 double d2ydx21=2.*(d2-d1)/(ii-bb);
	 double d2ydx22=2.*(d3-d2)/(hh-dd);
	 double d3ydx3=2*(d2ydx22-d2ydx21)/(hh+dd-ii-bb);
	 double ddd=d3ydx3/6.;
       double ccc=(d2ydx2-6.*ddd*(XXX-bb))/2.;
       double bbb=(dydx-2.*ccc*(XXX-bb)-3.*ddd*pow((XXX-bb),2));
       return bbb*(XXX-bb)+ccc*pow((XXX-bb),2)+ddd*pow((XXX-bb),3)+aa;}
/***************************************************************************************************/
// Interpolation for decay process.

double interpd(int N,double XXX, double Xmax1, double Xmax2, double Xmax3, double Xmax4) {
       double FF;
       double AA;
  int IX=0;
  int IY=N+1;
  int MID;
  double aa,bb,cc,dd,ee,ii,gg,hh;
  ifstream indata;

      for (int iy=0;(IY-IX) > 1;iy++) {
	if (XXX < Xmax1 ) {
    indata.open("ttbardecay1.dat", ios::in); 
  }else if (XXX > Xmax1  && XXX < Xmax2 ) {
    indata.open("ttbardecay2.dat", ios::in);
  }else if (XXX > Xmax2 && XXX < Xmax3 ) {
    indata.open("ttbardecay3.dat", ios::in);
  }else if (XXX > Xmax3) {
    indata.open("ttbardecay4.dat", ios::in);
  }  
        MID=(IX+IY)/2;
        for(int ix=0;ix < MID ;ix++) {
	  indata >> FF >> AA;}
	 
	 indata.close();
	 
        if (XXX >= AA) {
        IX=MID;} else {IY=MID;}
	if ((IY-IX) < 1) break;}
          if (IX < 5) {IX=5;}
    if (XXX < Xmax1 ) {
    indata.open("ttbardecay1.dat", ios::in); 
  }else if (XXX > Xmax1  && XXX < Xmax2 ) {
    indata.open("ttbardecay2.dat", ios::in);
  }else if (XXX > Xmax2 && XXX < Xmax3 ) {
    indata.open("ttbardecay3.dat", ios::in);
  }else if (XXX > Xmax3) {
    indata.open("ttbardecay4.dat", ios::in);
  }  
     for (int ix=0;ix<IX-3;ix++) {
        indata >> FF >> AA;}
        indata >> aa >> bb ;
        indata >> cc >> dd ;
        indata >> ee >> ii ;
        indata >> gg >> hh ;
	
        indata.close();
       double d1=(cc-aa)/(dd-bb);
       double d2=(ee-cc)/(ii-dd);
       double d3=(gg-ee)/(hh-ii) ;
       double dydx=(gg-aa)/(hh-bb);
       double dydx1=(ee-aa)/(ii-bb);
	 double dydx2=(gg-cc)/(hh-dd);
	 double d2ydx2=(dydx2-dydx1)/((hh+dd-ii-bb)/2.);
	 double d2ydx21=2.*(d2-d1)/(ii-bb);
	 double d2ydx22=2.*(d3-d2)/(hh-dd);
	 double d3ydx3=2*(d2ydx22-d2ydx21)/(hh+dd-ii-bb);
	 double ddd=d3ydx3/6.;
       double ccc=(d2ydx2-6.*ddd*(XXX-bb))/2.;
       double bbb=(dydx-2.*ccc*(XXX-bb)-3.*ddd*pow((XXX-bb),2));
       return bbb*(XXX-bb)+ccc*pow((XXX-bb),2)+ddd*pow((XXX-bb),3)+aa;}
/***************************************************************************************************/
// Find maximum and minimum values of variable x for production process

     void xmaxmin (double pt2, double d2, double &rxmax, double &rxmin5,double &rxmin6) {
         double xminr=39.*pt2-1.+3.*d2-42.*pt2*d2-3.*pow(d2,2)+pow(pt2,3)+15.*pow(pt2,2)+3.*pow(pt2,2)*d2+3.*pt2*pow(d2,2)+pow(d2,3);
   double xmini=6.*sqrt (abs(-3.*pt2+33.*pow(pt2,2)+3.*pow(pt2,3)-18.*pt2*pow(d2,2)+12.*pt2*pow(d2,3)+12.*pt2*d2-42.*pow(pt2,3)*d2+51.*pow(pt2,2)*pow(d2,2)-75.*pow(pt2,2)*d2-3.*d2*pow(pt2,4)-3.*pt2*pow(d2,4)-9.*pow(pt2,3)*pow(d2,2)-9.*pow(pt2,2)*pow(d2,3)));    
   double xminm=sqrt(pow(xminr,2)+pow(xmini,2));
   double xmint=atan2(xmini,xminr);
   double xmin1=-1./12.*pow(xminm,(1./3.))*cos(xmint/3.);
   double xmin2=(-1./12.+d2/6.-pow(pt2,2)/12.-5.*pt2/6.-pt2*d2/6.-pow(d2,2)/12.)*cos(xmint/3.)/(pow(xminm,(1./3.)));
   double xmin3=pt2/6.+5./6.+d2/6.;
   double xmin4=sqrt(3.)/12.*pow(xminm,(1./3.))*sin(xmint/3.)
               -3.*sqrt(3.)*(-1./36.+d2/18.-pow(pt2,2)/36.-5.*pt2/18.-pt2/18.*d2-pow(d2,2)/36.)
               /pow(xminm,(1./3.))*sin(xmint/3.);
      rxmin5=xmin1+xmin2+xmin3+xmin4;
      rxmin6=xmin1+xmin2+xmin3-xmin4;
      if (pt2 < 0.08333) {
        rxmax=1.-2.*pt2-sqrt(pt2)*sqrt(d2+4.*pt2);} else {
        rxmax=rxmin5;
        }
      }
/***************************************************************************************************/
// Calculate Matrix Element for production process

  double ME (double d2, double pt2,double &rx, double &ry1, double &ry2,double xmin5, double &rM1, double &rM2, double &rw, bool axi) {
         ry1=-(2.*pow(rx,2)-6.*rx-d2*rx+2.*d2+4.+2.*sqrt(-d2+pow(rx,2)+2.*rx*d2-d2*pow(rx,2)+4.*pt2*d2-4.*rx*pt2*d2-2.*pow(rx,3)+pow(rx,4)-d2*pow(rx,2)*pt2+4.*pt2*pow(rx,3)+pt2*pow(d2,2)-4*pt2*pow(rx,2)))/(4.*rx-4.-d2);
       ry2=-(2.*pow(rx,2)-6.*rx-d2*rx+2.*d2+4.-2.*sqrt(-d2+pow(rx,2)+2.*rx*d2-d2*pow(rx,2)+4.*pt2*d2-4.*rx*pt2*d2-2.*pow(rx,3)+pow(rx,4)-d2*pow(rx,2)*pt2+4.*pt2*pow(rx,3)+pt2*pow(d2,2)-4*pt2*pow(rx,2)))/(4.*rx-4.-d2);
//       Using the pt^2 delta function to eliminate x. w1,w2=d(pt)/dy
       rw=(4.*rx-4.-d2)/4.*(ry1-ry2)/(2.*sqrt(pt2)*(pow(rx,2)-d2));
         if (axi) {
	 rM1=((pow((rx+0.5*d2),2)+pow((ry1+0.5*d2),2)+0.5*d2*(pow((5.-rx-ry1),2)-19.+d2))/((1.-d2)*(1.-rx)*(1.-ry1))-d2/(2.*pow((1.-rx),2))-d2/(2.*pow((1.-ry1),2)))/(rw*sqrt(1.-d2));
	 rM2=((pow((rx+0.5*d2),2)+pow((ry2+0.5*d2),2)+0.5*d2*(pow((5.-rx-ry2),2)-19.+d2))/((1.-d2)*(1.-rx)*(1.-ry2))-d2/(2.*pow((1.-rx),2))-d2/(2.*pow((1.-ry2),2)))/(rw*sqrt(1.-d2)); } else {
	 rM1=((pow((rx+0.5*d2),2)+pow((ry1+0.5*d2),2)-2.*d2*(1.+d2/2.))/((1.+d2/2.)*(1.-rx)*(1.-ry1))-d2/(2.*pow((1.-rx),2))-d2/(2.*pow((1.-ry1),2)))/(rw*sqrt(1.-d2));
         rM2=((pow((rx+0.5*d2),2)+pow((ry2+0.5*d2),2)-2.*d2*(1.+d2/2.))/((1.+d2/2.)*(1.-rx)*(1.-ry2))-d2/(2.*pow((1.-rx),2))-d2/(2.*pow((1.-ry2),2)))/(rw*sqrt(1.-d2));	   }
       
       if (rx < xmin5 || pt2 > 0.024){rM1=0;}
       double M=fabs(rM1)+fabs(rM2);
       return M;}
/***************************************************************************************************/
// Calculate Matrix Element for decay process

double MEd (double pt2,double &rx3, double &rx11, double &rx12, double &rMd1, double &rMd2) {
  rx11=-0.5*(2*pow(mt,2)-3*pow(mt,2)*rx3+pow(mt,2)*pow(rx3,2)+pow(mt,2)*rx3*a-2*pt2*a+sqrt(-4*pt2*pow(mt,2)*pow(rx3,2)-4*pow(mt,2)*pt2*pow(a,2)-8*pow(mt,2)*rx3*pt2*a-4*pow(mt,2)*pow(rx3,2)*pt2*a+pow(mt,4)*pow(rx3,2)-2*pow(mt,4)*pow(rx3,3)+pow(mt,4)*pow(rx3,4)-4*pt2*pow(mt,2)+16*pow(pt2,2)*a+8*pow(mt,2)*pt2*a-2*pow(mt,4)*pow(rx3,2)*a+2*pow(mt,4)*pow(rx3,3)*a+pow(mt,4)*pow(rx3,2)*pow(a,2)+8*pt2*pow(mt,2)*rx3))/(-pt2-pow(mt,2)+pow(mt,2)*rx3);
   rx12=-0.5*(2*pow(mt,2)-3*pow(mt,2)*rx3+pow(mt,2)*pow(rx3,2)+pow(mt,2)*rx3*a-2*pt2*a-sqrt(-4*pt2*pow(mt,2)*pow(rx3,2)-4*pow(mt,2)*pt2*pow(a,2)-8*pow(mt,2)*rx3*pt2*a-4*pow(mt,2)*pow(rx3,2)*pt2*a+pow(mt,4)*pow(rx3,2)-2*pow(mt,4)*pow(rx3,3)+pow(mt,4)*pow(rx3,4)-4*pt2*pow(mt,2)+16*pow(pt2,2)*a+8*pow(mt,2)*pt2*a-2*pow(mt,4)*pow(rx3,2)*a+2*pow(mt,4)*pow(rx3,3)*a+pow(mt,4)*pow(rx3,2)*pow(a,2)+8*pt2*pow(mt,2)*rx3))/(-pt2-pow(mt,2)+pow(mt,2)*rx3);
   double diff1=0.5*(-pow(mt,2)*(rx11+rx3*(2.-rx11-a)-pow(rx3,2)-1.)/(pow(rx11+a,2)-4*a)+pow(mt,2)*(1.-rx11)*(1.-rx3)/(pow(rx11+a,2)-4.*a)-pow(mt,2)*(1.-rx11)*(rx11+rx3*(2.-rx11-a)-pow(rx3,2)-1.)*(2.*rx11+2.*a)/(pow((pow(rx11+a,2)-4.*a),2)))/sqrt(pow(mt,2)*(1.-rx11)*(rx11+rx3*(2.-rx11-a)-pow(rx3,2)-1.)/(pow(rx11+a,2)-4.*a));
  double diff2=0.5*(-pow(mt,2)*(rx12+rx3*(2.-rx12-a)-pow(rx3,2)-1.)/(pow(rx12+a,2)-4.*a)+pow(mt,2)*(1.-rx12)*(1.-rx3)/(pow(rx12+a,2)-4.*a)-pow(mt,2)*(1.-rx12)*(rx12+rx3*(2.-rx12-a)-pow(rx3,2)-1.)*(2.*rx12+2.*a)/(pow((pow(rx12+a,2)-4.*a),2)))/sqrt(pow(mt,2)*(1.-rx12)*(rx12+rx3*(2.-rx12-a)-pow(rx3,2)-1.)/(pow(rx12+a,2)-4.*a));
   rMd1=4./(pi*3.*(1.-rx11)*pow(rx3,2))*(rx3-((1.-rx11)*(1.-rx3)+pow(rx3,2))/(1.-a)+rx3*pow(rx11+rx3-1.,2)/(2.*pow(1.-a,2))+2.*a*(1.-rx11)*pow(rx3,2)/(pow(1.-a,2)*(1.+2.*a)))/diff1;
   rMd2=4./(pi*3.*(1.-rx12)*pow(rx3,2))*(rx3-((1.-rx12)*(1.-rx3)+pow(rx3,2))/(1.-a)+rx3*pow(rx12+rx3-1.,2)/(2.*pow(1.-a,2))+2.*a*(1.-rx12)*pow(rx3,2)/(pow(1.-a,2)*(1.+2.*a)))/diff2;
   double M=fabs(rMd1)+fabs(rMd2);
       return M;
}
/***************************************************************************************************/
// Calculate Polarized Matrix Element for production process

  void PolarizedME (double xt, double xtb, double xg, double &rct, double &rst, double &rsphi, double &rcphi,double &rLLL,double &rLRR,double &rLRL,double &rLLR,double &rRLL,double &rRRR, double &rRRL,double &rRLR) {
    double F1VL=-F1Vg+(-0.5+sinsqthw)/(sqrt(sinsqthw)*costhw)*pow(emcm,2)/(pow(emcm,2)-pow(Mz,2))*F1Vz;
    double F1VR=-F1Vg+sinsqthw/(sqrt(sinsqthw)*costhw)*pow(emcm,2)/(pow(emcm,2)-pow(Mz,2))*F1Vz;
    double F1AL=(-0.5+sinsqthw)/(sqrt(sinsqthw)*costhw)*pow(emcm,2)/(pow(emcm,2)-pow(Mz,2))*F1Az;
    double F1AR=sinsqthw/(sqrt(sinsqthw)*costhw)*pow(emcm,2)/(pow(emcm,2)-pow(Mz,2))*F1Az;
    double F1RL=F1VL+F1AL;
    double F1RR=F1VR+F1AR;
    double F1LL=F1VL-F1AL;
    double F1LR=F1VR-F1AR;
    //redfine vertex factor to include soft and virtual contributions
    double vF1VL=-(F1Vg+dF1Vg)+(-0.5+sinsqthw)/(sqrt(sinsqthw)*costhw)*pow(emcm,2)/(pow(emcm,2)-pow(Mz,2))*(F1Vz+dF1Vz);
    double vF1VR=-(F1Vg+dF1Vg)+sinsqthw/(sqrt(sinsqthw)*costhw)*pow(emcm,2)/(pow(emcm,2)-pow(Mz,2))*(F1Vz+dF1Vz);
    double vF1AL=(-0.5+sinsqthw)/(sqrt(sinsqthw)*costhw)*pow(emcm,2)/(pow(emcm,2)-pow(Mz,2))*(F1Az+dF1Az);
    double vF1AR=sinsqthw/(sqrt(sinsqthw)*costhw)*pow(emcm,2)/(pow(emcm,2)-pow(Mz,2))*(F1Az+dF1Az);

    double areafac = 16.8374; //Inverse of the area of the 3-body phase space. Hard coded for mt=175GeV, emcm=500GeV, Egmin=0.2GeV. 
       double bt=sqrt(1.-4.*pow(mt,2)/(pow(xt*emcm,2)));
       double btb=sqrt(1.-4.*pow(mt,2)/(pow(xtb*emcm,2)));
       double tt=acos((xg-xt-xtb+xt*xtb+4.*pow(mt/emcm,2))/(xt*bt*xtb*btb));
       double tg=acos((xtb-xt-xg+xt*xg)/(xt*bt*xg));
       double gam=emcm/(2.*mt);
       double GLLR=pow(((vF1VL-b*vF1AL)*(1.+rct)),2);
       double GLRL=pow(((vF1VL+b*vF1AL)*(1.-rct)),2);
       double GLLL=pow(((vF1VL)*rst/gam),2);
       double GLRR=pow(((vF1VL)*rst/gam),2);
       double GRLR=pow((-(vF1VR-b*vF1AR)*(1.-rct)),2);
       double GRRL=pow((-(vF1VR+b*vF1AR)*(1.+rct)),2);
       double GRLL=pow(((vF1VR)*rst/gam),2);
       double GRRR=pow(((vF1VR)*rst/gam),2);
       double Apm=xg*sqrt(xt*xtb*(1.+bt)*(1.-btb))/(emcm*(1.-xt)*(1.-xtb));
       double App=xg*sqrt(xt*xtb*(1.+bt)*(1.+btb))/(emcm*(1.-xt)*(1.-xtb));
       double Amm=xg*sqrt(xt*xtb*(1.-bt)*(1.-btb))/(emcm*(1.-xt)*(1.-xtb));
       double Amp=xg*sqrt(xt*xtb*(1.-bt)*(1.+btb))/(emcm*(1.-xt)*(1.-xtb));
       double LLLLL=-Apm*sin(tg/2.)*cos(tg/2.)*cos(tt/2.)*(xt*bt+(1.-xt));
       double LRLLL=Apm*pow((sin(tg/2.)),2)*sin(tt/2.)*(1.-xt);
       double LZLLL=-Apm/sqrt(2.)*sin(tg/2.)*(xt*bt*cos(tg/2.)*sin(tt/2.)+(1.-xt)*sin((tt-tg)/2.));
       double LLLLR=Apm*sin(tg/2.)*(xt*bt*cos(tg/2.)*cos(tt/2.)+(1.-xtb)*cos((tt+tg)/2.));
       double LZLLR=Apm/sqrt(2.)*cos(tg/2.)*(xt*bt*sin(tg/2.)*sin(tt/2.)+(1.-xtb)*cos((tt+tg)/2.));
       double LLRLL=-Amm*pow((cos(tg/2.)),2)*cos(tt/2.)*(1.-xt);
       double LRRLL=-Amm*sin(tg/2.)*cos(tg/2.)*sin(tt/2.)*(xt*bt-(1.-xt));
       double LZRLL=-Amm/sqrt(2.)*cos(tg/2.)*(xt*bt*sin(tg/2.)*cos(tt/2.)+(1.-xt)*sin((tt-tg)/2.));
       double LRRLR=Amm*cos(tg/2.)*(xt*bt*sin(tg/2.)*sin(tt/2.)+(1.-xtb)*cos((tt+tg)/2.));
       double LZRLR=Amm/sqrt(2.)*sin(tg/2.)*(xt*bt*cos(tg/2.)*cos(tt/2.)+(1.-xtb)*cos((tt+tg)/2.));
       double LLLRL=-App*sin(tg/2.)*cos(tg/2.)*sin(tt/2.)*(xt*bt+(1.-xt));
       double LRLRL=-App*pow((sin(tg/2.)),2)*cos(tt/2.)*(1.-xt);
       double LZLRL=App/sqrt(2.)*sin(tg/2.)*(xt*bt*cos(tg/2.)*cos(tt/2.)+(1.-xt)*cos((tt-tg)/2.));
       double LLLRR=App*sin(tg/2.)*(xt*bt*cos(tg/2.)*sin(tt/2.)+(1.-xtb)*sin((tt+tg)/2.));
       double LZLRR=-App/sqrt(2.)*cos(tg/2.)*(xt*bt*sin(tg/2.)*cos(tt/2.)-(1.-xtb)*sin((tt+tg)/2.));
       double LLRRL=-Amp*pow((cos(tg/2.)),2)*sin(tt/2.)*(1.-xt);
       double LRRRL=Amp*sin(tg/2.)*cos(tg/2.)*cos(tt/2.)*(xt*bt-(1.-xt));
       double LZRRL=-Amp/sqrt(2.)*cos(tg/2.)*(xt*bt*sin(tg/2.)*sin(tt/2.)-(1.-xt)*cos((tt-tg)/2.));
       double LRRRR=-Amp*cos(tg/2.)*(xt*bt*sin(tg/2.)*cos(tt/2.)-(1.-xtb)*sin((tt+tg)/2.));
       double LZRRR=Amp/sqrt(2.)*sin(tg/2.)*(xt*bt*cos(tg/2.)*sin(tt/2.)+(1.-xtb)*sin((tt+tg)/2.));
	double LLLl=LLLLL+LLLLR;
      double LLLr=LRLLL;
      double LLLz=LZLLL+LZLLR;
      double LRRl=LLRRL;
      double LRRr=LRRRL+LRRRR;
      double LRRz=LZRRL+LZRRR;
      double LRLl=LLRLL;
      double LRLr=LRRLL+LRRLR;
      double LRLz=LZRLL+LZRLR;
      double LLRl=LLLRR+LLLRL;
      double LLRr=LRLRL;
      double LLRz=LZLRR+LZLRL;
      double RRRr=LLLl;
      double RRRl=LLLr;
      double RRRz=-LLLz;
      double RLLr=LRRl;
	double RLLl=LRRr;
	double RLLz=-LRRz;
	double RLRr=-LRLl;
	double RLRl=-LRLr;
	double RLRz=LRLz;
	double RRLr=-LLRl;
	double RRLl=-LLRr;
	double RRLz=LLRz;  
	double MLLL[4],MLRR[4],MLRL[4],MLLR[4],MRLL[4],MRRR[4],MRRL[4],MRLR[4];
	MLLL[1]=(F1LL*LLLl+F1RL*RLLl)*(1.+rct)/sqrt(2.);
	MLLL[2]=(F1LL*LLLr+F1RL*RLLr)*(1.-rct)/sqrt(2.);
	MLLL[3]=+(F1LL*LLLz+F1RL*RLLz)*rst;
	MLRR[1]=(F1LL*LRRl+F1RL*RRRl)*(1.+rct)/sqrt(2.);
	MLRR[2]=+(F1LL*LRRr+F1RL*RRRr)*(1.-rct)/sqrt(2.);
	MLRR[3]=+(F1LL*LRRz+F1RL*RRRz)*rst;
	MLRL[1]=(F1LL*LRLl+F1RL*RRLl)*(1.+rct)/sqrt(2.);
	MLRL[2]=+(F1LL*LRLr+F1RL*RRLr)*(1.-rct)/sqrt(2.);
	MLRL[3]=+(F1LL*LRLz+F1RL*RRLz)*rst;
	MLLR[1]=(F1LL*LLRl+F1RL*RLRl)*(1.+rct)/sqrt(2.);
	MLLR[2]=+(F1LL*LLRr+F1RL*RLRr)*(1.-rct)/sqrt(2.);
	MLLR[3]=+(F1LL*LLRz+F1RL*RLRz)*rst;
	MRLL[1]=-(F1LR*LLLl+F1RR*RLLl)*(1.-rct)/sqrt(2.);
	MRLL[2]=-(F1LR*LLLr+F1RR*RLLr)*(1.+rct)/sqrt(2.);
	MRLL[3]=+(F1LR*LLLz+F1RR*RLLz)*rst;
	MRRR[1]=-(F1LR*LRRl+F1RR*RRRl)*(1.-rct)/sqrt(2.);
	MRRR[2]=-(F1LR*LRRr+F1RR*RRRr)*(1.+rct)/sqrt(2.);
	MRRR[3]=+(F1LR*LRRz+F1RR*RRRz)*rst;
	MRRL[1]=-(F1LR*LRLl+F1RR*RRLl)*(1.-rct)/sqrt(2.);
	MRRL[2]=-(F1LR*LRLr+F1RL*RRRr)*(1.+rct)/sqrt(2.);
	MRRL[3]=+(F1LR*LRLz+F1RR*RRLz)*rst;
	MRLR[1]=-(F1LR*LLRl+F1RR*RLRl)*(1.-rct)/sqrt(2.);
	MRLR[2]=-(F1LR*LLRr+F1RR*RLRr)*(1.+rct)/sqrt(2.);
	MRLR[3]=+(F1LR*LLRz+F1RR*RLRz)*rst;
	rLLL=(pow(((MLLL[1]+MLLL[2])*rcphi+MLLL[3]),2)+pow(((MLLL[1]+MLLL[2])*rsphi),2))*16.*pi*asmt*pow(emcm,2)/(16.*pow(pi,2))+GLLL*areafac;
	rLRR=(pow(((MLRR[1]+MLRR[2])*rcphi+MLRR[3]),2)+pow(((MLRR[1]+MLRR[2])*rsphi),2))*16.*pi*asmt*pow(emcm,2)/(16.*pow(pi,2))+GLRR*areafac;
	rLRL=(pow(((MLRL[1]+MLRL[2])*rcphi+MLRL[3]),2)+pow(((MLRL[1]+MLRL[2])*rsphi),2))*16.*pi*asmt*pow(emcm,2)/(16.*pow(pi,2))+GLRL*areafac;
	rLLR=(pow(((MLLR[1]+MLLR[2])*rcphi+MLLR[3]),2)+pow(((MLLR[1]+MLLR[2])*rsphi),2))*16.*pi*asmt*pow(emcm,2)/(16.*pow(pi,2))+GLLR*areafac;
	rRLL=(pow(((MRLL[1]+MRLL[2])*rcphi+MRLL[3]),2)+pow(((MRLL[1]+MRLL[2])*rsphi),2))*16.*pi*asmt*pow(emcm,2)/(16.*pow(pi,2))+GRLL*areafac;
	rRRR=(pow(((MRRR[1]+MRRR[2])*rcphi+MRRR[3]),2)+pow(((MRRR[1]+MRRR[2])*rsphi),2))*16.*pi*asmt*pow(emcm,2)/(16.*pow(pi,2))+GRRR*areafac;
	rRRL=(pow(((MRRL[1]+MRRL[2])*rcphi+MRRL[3]),2)+pow(((MRRL[1]+MRRL[2])*rsphi),2))*16.*pi*asmt*pow(emcm,2)/(16.*pow(pi,2))+GRRL*areafac;
	rRLR=(pow(((MRLR[1]+MRLR[2])*rcphi+MRLR[3]),2)+pow(((MRLR[1]+MRLR[2])*rsphi),2))*16.*pi*asmt*pow(emcm,2)/(16.*pow(pi,2))+GRLR*areafac;
      }
/***************************************************************************************************/
// Calculate Polarized Matrix Element for decay process

double PolMEdec(double xw1, double xb1, double xg1,double &rctw, double &rctl,double &rbphi, double &rbgphi,int tm, double xd0, double zd0){
     
        double fac=pow(mt,2)*(1.-a)/2.*34.606; //Includes inverse of 3body phase space area (34.606). hardcoded for mt=175GeV, mw=80GeV, Egmin=0.2GeV, zd0=0.0028;
        	double F1WL=1.;
		double dF1WL=asmt*2./(3.*pi)*(-3.-pow(log(1.-a),2)+1.5*log(1-a)+4.+(5.-3.*a)/(8.*(1.-a))+pow(log(1.-a),2)-2.5*log(1.-a)-a*(2.-3.*a)/(4.*pow((1.-a),2))*log(a)-pow(pi,2)/2.-(1.+log(xd0))*(1.+log(zd0))+0.25*log(zd0)+1.05688-0.22105); //hardcoded for dilog (1-a)=1.05688 and dilog(a)=0.22105
		double dF2WL=asmt*2./(3.*pi)*1./a*log(1.-a);
        double TL=1./a*(1.+tm*rctw)*(1.-pow(rctl,2))*fac*pow((F1WL+dF1WL-0.5*a*dF2WL),2);
	double TO=(1.-tm*rctw)*pow((1.-rctl),2)*fac*pow((F1WL+dF1WL-0.5*dF2WL),2);
	double TI=tm*2./sqrt(a)*(1.-rctl)*sqrt(1.-pow(rctw,2))*sqrt(1.-pow(rctl,2))*cos(rbphi)*fac*(F1WL+dF1WL-0.5*a*dF2WL)*(F1WL+dF1WL-0.5*dF2WL);
	double T1sum=TL+TO+TI;
       double bw1=sqrt(1.-4.*a/pow(xw1,2));
       double cwb=1./(xw1*bw1*xb1)*(xg1-xw1-xb1+xw1*xb1+2.*a);
       double cwg=1./(xw1*bw1*xg1)*(xb1-xw1-xg1+xw1*xg1+2.*a);
       double twb=acos(cwb);
       double twg=acos(cwg);
       double swg=sin(twg);
       double cwg2=cos(twg/2.);
       double cwb2=cos(twb/2.); 
       double cbg2=cos((twb+twg)/2.);
       double swg2=sin(twg/2.);
       double swb2=sin(twb/2.);
       double zet1=1.+a-xw1;
       double MR1=-2./sqrt(zet1*xg1)*(xg1*cwg2+xb1*cwb2*cbg2);
       double ML1=2./sqrt(zet1*xg1)*(-xg1*swg2+xb1*swb2*cbg2);
       double MR2=xw1*(1.+bw1)/sqrt(a*2.*zet1*xg1)*(-xg1*swg2+xb1*swb2*cbg2);
       double ML2=-xw1*(1.-bw1)/sqrt(a*2.*zet1*xg1)*(xg1*cwg2+xb1*cwb2*cbg2);
       double MR3=sqrt(xb1)*cwb2*(2.*sqrt(xb1/(zet1*xg1))*cbg2-swg);
       double ML3=-sqrt(xb1)*cwb2*(1.-cwg);
       double MR4=-sqrt(xb1)*swb2*(1.+cwg);
       double ML4=-sqrt(xb1)*swb2*(2.*sqrt(xb1/(zet1*xg1))*cbg2+swg);
       double MR5=xw1*sqrt(xb1)/(2.*sqrt(2.*a))*((1.+bw1)*swb2*(-2.*sqrt(xb1/(zet1*xg1))*cbg2+swg)+(1.-bw1)*cwb2*(1.+cwg));
       double ML5=xw1*sqrt(xb1)/(2.*sqrt(2.*a))*((1.-bw1)*cwb2*(2.*sqrt(xb1/(zet1*xg1))*cbg2+swg)+(1.+bw1)*swb2*(1.-cwg));
       double A=(MR2+MR5)*sqrt(1.-pow(rctl,2)); double B=(MR1+MR3)*1./sqrt(2.)*(1.+rctl); double C=MR4*1./sqrt(2.)*(1.-rctl);
       double D=(ML2+ML5)*sqrt(1.-pow(rctl,2)); double E=ML3*1./sqrt(2.)*(1.+rctl); double F=(ML1+ML4)*1./sqrt(2.)*(1.-rctl);
       double MMT=((pow(A,2)+pow(B,2)+pow(C,2)+2.*A*B*cos(rbphi)+2.*B*C*cos(2.*rbphi)+2.*A*C*cos(rbphi))*0.5*(1.+tm*rctw)+(pow(D,2)+pow(E,2)+pow(F,2)+2.*D*E*cos(rbphi)+2.*E*F*cos(2.*rbphi)+2.*D*F*cos(rbphi))*0.5*(1.-tm*rctw)+tm*((A*D+B*E+C*F)*cos(rbgphi)+(A*E+C*D)*cos(rbphi+rbgphi)+(B*D+A*F)*cos(rbphi-rbgphi)+B*F*cos(2.*rbphi-rbgphi)+C*E*cos(2.*rbphi+rbgphi))*sqrt(1-pow(rctw,2)))*pow(mt,2)*pi*asmt*8./(pow(pi,2)*16.)+T1sum;
	return MMT;
}

/***************************************************************************************************/
// Random number generator

  double random (int &rseed) {
    int M = 2147483647;
    int A = 16807;
    int Q = 127773;
    int R = 2836;
    double MINV = 0.46566128752458e-09;
    int HI = rseed/Q;
    int LO = rseed % Q;
    rseed = A*LO - R*HI;
      if (rseed <= 0) {rseed = rseed + M;}
    return rseed*MINV;}  
/***************************************************************************************************/

// 			 if(POWprod && !POWdecay) {
// 			   for (int ja = 1; ja < 12; ja++) {
// 			     outdata << IDUP[ja] << "\t" << ISTUP[ja] << "\t" << MOTHUP[1][ja] <<"\t" << MOTHUP[2][ja] << "\t" << ICOLUP[1][ja] << "\t" << ICOLUP[2][ja] <<"\t" <<setprecision (9)<< PUP[1][ja] << "\t" << PUP[2][ja]<<"\t"<<PUP[3][ja]<<"\t"<<PUP[4][ja]<<"\t" << mass[ja] <<"\t"<< "0" << "\t" <<"9" <<endl;
// 			     sum1+=PUP[4][ja];sum2+=PUP[1][ja];sum3+=PUP[2][ja];sum3+=PUP[3][ja];
// }
// 			   if (emit) {
// 			   for (int ja=14; ja < 16; ja++) {
// 			     outdata << IDUP[ja] << "\t" << ISTUP[ja] << "\t" << MOTHUP[1][ja] <<"\t" << MOTHUP[2][ja] << "\t" << ICOLUP[1][ja] << "\t" << ICOLUP[2][ja] <<"\t" <<setprecision (9)<< PUP[1][ja] << "\t" << PUP[2][ja]<<"\t"<<PUP[3][ja]<<"\t"<<PUP[4][ja]<<"\t" << mass[ja] <<"\t"<< "0" << "\t" <<"9" <<endl;
// 			     sum1+=PUP[4][ja];sum2+=PUP[1][ja];sum3+=PUP[2][ja];sum3+=PUP[3][ja];}}else {
// 			       for (int ja=14; ja < 15; ja++) {
// 			     outdata << IDUP[ja] << "\t" << ISTUP[ja] << "\t" << MOTHUP[1][ja] <<"\t" << MOTHUP[2][ja] << "\t" << ICOLUP[1][ja] << "\t" << ICOLUP[2][ja] <<"\t" <<setprecision (9)<< PUP[1][ja] << "\t" << PUP[2][ja]<<"\t"<<PUP[3][ja]<<"\t"<<PUP[4][ja]<<"\t" << mass[ja] <<"\t"<< "0" << "\t" <<"9" <<endl;
//  			     sum1+=PUP[4][ja];sum2+=PUP[1][ja];sum3+=PUP[2][ja];sum3+=PUP[3][ja];
// 			     }
			   
// 			     }}
			 
//  			 //   if (decay) {for (int jh=11; jh< NUP+1; jh++)
			     
// //  			   {jhh=jh+1;outdata << jhh <<"\t" << PUP[4][jh] << "\t" << PUP[1][jh]<<"\t"<<PUP[2][jh]<<"\t"<<PUP[3][jh]<<endl;
// //  			 sum1+=PUP[4][jh];sum2+=PUP[1][jh];sum3+=PUP[2][jh];sum3+=PUP[3][jh];}}
// // 			 	 if (POWprod && !decay) {for (int jh=13; jh<NUPP;jh++)
// //  			   {jhh=jh+1;outdata << jhh-2 <<"\t" << PUP[4][jh] << "\t" << PUP[1][jh]<<"\t"<<PUP[2][jh]<<"\t"<<PUP[3][jh]<<endl;
// //  			   sum1+=2.*PUP[4][jh];sum2+=2.*PUP[1][jh];sum3+=2.*PUP[2][jh];sum3+=2.*PUP[3][jh];}}

//  				 cout <<ixx  << "\t" << sum1 << "\t" <<sum2 <<"\t" <<sum3 <<"\t" << sum4 << "\t" << EVT <<  endl;
// 			outdata << "</event>" << endl;
//  }
  //  outdata << "</LesHouchesEvents>" << endl;
//  return 0;
//  } 

			 //			 	   outdata << NUP+1 <<"\t"<<ixx+1<<"\t"<<XWGTUP<<"\t"<<SCALUP<<"\t"<<AQEDUP<<"\t"<<AQCDUP <<endl;
	// 		 for (int ja = 1; ja < 11; ja++) {
// 			 outdata << IDUP[ja] << "\t" ;
// 			 if (ja==2)  {outdata << IDUP[15] << "\t" ;}}
// 			 if (decay) {for (int ja=11; ja < NUP+1; ja++)
// 			 {outdata << IDUP[ja] << "\t" ;}}
// 			 if (POWprod && !decay) {for (int ja=13; ja<NUPP;ja++)
// 			 {outdata << IDUP[ja] << "\t" ;}}
// 			 outdata << endl;
// 			 for (int jb = 1; jb < 3; jb++) {
// 			 for (int jc = 1; jc < 11; jc++) {
// 			 outdata << MOTHUP[jb][jc]<<"\t" ;
// 			 if (jc==2)  {outdata << MOTHUP[jb][15] << "\t" ;}}
// 			 if (decay) {for (int jc=11; jc < NUP+1; jc++)
// 			 {outdata << MOTHUP[jb][jc] << "\t" ;}}
// 			 if (POWprod && !decay) {for (int jc=13; jc<NUPP;jc++)
// 			 {outdata << MOTHUP[jb][jc] << "\t" ;}}
// 			 outdata << endl;}
// 			 for (int jd = 1; jd < 3; jd++) {
// 			    for (int je = 1; je < 11; je++) {
// 			  outdata << ICOLUP[jd][je]<<"\t" ;
// 			  if (je==2)  {outdata << ICOLUP[jd][15] << "\t" ;}}
// 			 if (decay) {for (int je=11; je < NUP+1; je++)
// 			    {outdata << ICOLUP[jd][je] << "\t" ;}}
// 			 if (POWprod && !decay) {for (int je=13; je<NUPP;je++)
// 			 {outdata << ICOLUP[jd][je] << "\t" ;}}
// 			  outdata << endl;}
// 			 for (int jf = 1; jf < 11; jf++) {
// 			 outdata << ISTUP[jf] << "\t" ;
// 			 if (jf==2) { outdata << ISTUP[15] << "\t" ;}}
// 			 if (decay) {for (int jf=11; jf< NUP+1; jf++)
// 			 {outdata << ISTUP[jf] << "\t" ;}}
// 			  if (POWprod && !decay) {for (int jf=13; jf<NUPP;jf++)
// 			 {outdata << ISTUP[jf] << "\t" ;}}
// 			  outdata << endl;
// 			 for (int jg = 1; jg < 11; jg++) {
// 			  outdata << ISPINUP[jg] << "\t" ;
// 			  if (jg==2) { outdata << ISPINUP[15] << "\t" ;}}
// 			 if (decay) {for (int jg=11; jg< NUP+1; jg++)
// 			 {outdata << ISPINUP[jg] << "\t" ;}}
// 			  if (POWprod && !decay) {for (int jg=13; jg<NUPP;jg++)
// 			  {outdata << ISPINUP[jg] << "\t" ;}}
// 			   outdata << endl;
// 			   int jhh=0;
//  			 for (int jh = 1; jh < 11; jh++) {
// 			   jhh=jhh+1;
// 			    outdata << jhh <<"\t" << PUP[4][jh] << "\t" << PUP[1][jh]<<"\t"<<PUP[2][jh]<<"\t"<<PUP[3][jh]<<endl;
// 			  if(jh==2){outdata << 3 <<"\t" << PUP[4][15] << "\t" << PUP[1][15]<<"\t"<<PUP[2][15]<<"\t"<<PUP[3][15]<<endl;jhh=jhh+1;}
// 			 sum1+=PUP[4][jh];sum2+=PUP[1][jh];sum3+=PUP[2][jh];sum3+=PUP[3][jh];}
// 			   if (decay) {for (int jh=11; jh< NUP+1; jh++)
			     
// 			   {jhh=jh+1;outdata << jhh <<"\t" << PUP[4][jh] << "\t" << PUP[1][jh]<<"\t"<<PUP[2][jh]<<"\t"<<PUP[3][jh]<<endl;
// 			 sum1+=PUP[4][jh];sum2+=PUP[1][jh];sum3+=PUP[2][jh];sum3+=PUP[3][jh];}}
// 			 	 if (POWprod && !decay) {for (int jh=13; jh<NUPP;jh++)
// 			   {jhh=jh+1;outdata << jhh-2 <<"\t" << PUP[4][jh] << "\t" << PUP[1][jh]<<"\t"<<PUP[2][jh]<<"\t"<<PUP[3][jh]<<endl;
// 			   sum1+=2.*PUP[4][jh];sum2+=2.*PUP[1][jh];sum3+=2.*PUP[2][jh];sum3+=2.*PUP[3][jh];}}

// 				 cout <<ixx  << "\t" << sum1 << "\t" <<sum2 <<"\t" <<sum3 <<"\t" << sum4 << "\t" << EVT <<  endl;}

// return 0; }
