#include <cstdlib> 
#include <ctime> 
#include <iostream>
#include <fstream>
#include <string>
//#include "LHAPDFWrap.h"
#include <stdio.h>
#include "MCHG_INPUTS.h"
#include <cmath>
#include <iomanip>
using namespace std;
// Checks if phase space point is in soft and collinear region? 
void QQgsc(double r1,double r2,double r3,double r4, double &rx, double &ry, double param[9], double &rwgt, double &rwgt1,  bool resolve);
// Checks if phase space point is in hard region? 
void QQgh(double r1,double r2,double r3,double r4, double &rx, double &ry, double param[9], double &rwgt,  double &rwgt1,  bool resolve);
// Checks if phase space point is in gap between phase space for massive partons and full half-triangle? 
void QQgsp(double r1,double r2,double r3,double r4, double &rx, double &ry, double param[9], double &rwgt, double &rwgt1,  bool resolve);
// Dilogarithm Function
double ddilog(double x);
// Random number generator
double random (int &rseed);

// Main program 

int main() {
  Input user;
  int nev = user.nevint();
  int nevg = user.nevgen();
  int nit = user.it();
  bool iface = user.iface();
  bool rest = user.rest();
  bool ifile = user.ifile();
  char* filename;
  if (iface && ifile) {
    filename = user.filename(); 
    rest = false;
    nevg = user.nfile();
  }
  if (iface && !ifile && !rest) 
  {cout << "Need interface file!" << endl;
  return 0;}
  double emcm = user.cme();
  double mH = user.mH();
  int seed = user.rseed();
  bool resolve = user.resolve();
  int ii=0; //event count
  // Constants 
  double pi = 3.14152654;
  double GF = 0.0000116639;
  double m  = user.mb();
  double Mz = user.Mz();  
  double Mw = user.Mw();
  double Yz = 2.486;
  double alphaem = 0.007297352;
  double alphas = user.alphasmH();
  //Lepton charge & coupling constants
  double Che = -1;
  double sin2thw = 0.2312;
  double cos2thw = 1-sin2thw;
  double Ae = -0.5;
  double Ve = Ae-2.*Che*sin2thw;
  double CLe = (Ae+Ve)/2.;
  double CRe = (Ve-Ae)/2.;
  double gz2 = sqrt(32.)*GF*pow(Mz,2);
  double e2 = 4.*pi*alphaem;
  double k = sqrt(2.)*GF*pow(Mz,2)/e2;
  double chi1 =0.;
  double chi2 = pow(k,2)*pow(Mz,4)/(pow(Yz*Mz,2));
  //quark masses ,velocities and couplings. NB:gluon mass[10] =0.75 
 
  double Ch[6];
  double mg = 0.75;
  Ch[1]=-1./3.; Ch[2]=2./3.; Ch[3]=Ch[1]; Ch[4]=Ch[2]; Ch[5]=Ch[1];
 /* w, w1 are the integrands. p,p1 are the corresponding integrals calculated from sum,sum1 & sd, sd1  are the corresponding errors 
     calculated via ssq, ssq1. bigwt is the maximum value of the integrand over the phase space (will be used for event generation) */
    
  double wgt, wgt1;
  double w[3],w1[3],bigwt[3],sum[3],sum1[3],ssq[3],ssq1[3],sd[3],sd1[3],p[3],p1[3];

  sum1[0]=0.;sum1[1]=0.;sum[2]=0.;
  sum1[0]=0.;sum1[1]=0.;sum1[2]=0.;
  ssq[0]=0.;ssq[1]=0.;ssq[2]=0.;
  ssq1[0]=0.;ssq1[1]=0.;ssq1[2]=0.;
  bigwt[0]=0.;bigwt[1]=0.;bigwt[2]=0.;

  cout << "Integrating..." << endl;

  for (int d = 0; d < nev; d++) {
  //Calculate  Born cross-section

  // Parameters for NLO cross-section calculation.
  double param[9];
  param[0] = pow(m,2)/pow(mH,2);
  param[1] = 2.*param[0];
  param[2] = sqrt(1.-4.*param[0]); // beta
  param[3] = (1.-pow(param[2],2))/param[2]*log((1.+param[2])/(1.-param[2]))+ (1.+pow(param[2],2))/param[2]*(4.*ddilog((1.-param[2])/(1.+param[2]))+2.*log((1.+param[2])/(2.*param[2]))*log((1.+param[2])/(1.-param[2]))+2.*ddilog((param[2]-1.)/(1.+param[2]))-log(2./(1.+param[2]))*log((1.+param[2])/(1.-param[2])))-3.*log(4./(1.-pow(param[2],2)))-4.*log(param[2])+1./(16.*pow(param[2],3))*(3.+18.*pow(param[2],2)+3.*pow(param[2],4))*log((1.+param[2])/(1.-param[2]))-3./8.*pow(param[2],2)+21./8.;
  param[4] = 0.5*(1.+param[2]);
  param[5] = (param[0]*param[4]+0.25*param[2]*pow(1+param[2],2))/(param[4]-param[0]);
  param[6] = 4.*alphas/(3.*pi); // colour factor
  param[7] = 3.*GF*pow(m,2)/(4.*sqrt(2)*pi)*mH*pow(param[2],3);
  param[8] = 1./(2.*pow(param[2],3));

 // Random numbers and phase space points x and y.
  double r1=random(seed);
  double r2=random(seed);
  double r3=random(seed);
  double r4=random(seed);
  double x,y;
  x = r1;
  y = r2;
  QQgsc(r1, r2, r3, r4, x, y, param, wgt, wgt1, resolve);
  w[0] = wgt; w1[0]=wgt1;
  QQgh(r1, r2, r3, r4, x, y, param, wgt, wgt1, resolve);
  w[1] = wgt; w1[1]=wgt1;
  QQgsp(r1, r2, r3, r4, x, y, param, wgt, wgt1, resolve);
  w[2] = wgt; w1[2]=wgt1;
    
  for (int iz = 0; iz < 3; iz++) {
    if (abs(w[iz]) > abs(bigwt[iz])) {bigwt[iz] = w[iz];}
    sum[iz]+=w[iz];
    sum1[iz]+=w1[iz];
    ssq[iz]+=pow(w[iz],2);
    ssq1[iz]+=pow(w1[iz],2);
    
}
  cout << "Integrated event : " << "\t" << d+1 << "\t" << "of" << "\t" << nev << "\r" << flush;  

}
  cout << endl;
  double xsec=0;
  double cstt=0;
  double xsecabs = 0.;
  for (int ic=0; ic<3; ic++) {
         p[ic] = sum[ic]/nev;
	 p1[ic] = sum1[ic]/nev;
        sd[ic] = sqrt((ssq[ic]/nev-pow(p[ic],2))/nev);
         sd1[ic] = sqrt((ssq1[ic]/nev-pow(p1[ic],2))/nev);  
	 xsec+=p[ic]+p1[ic];
	 xsecabs+=fabs(p[ic])+fabs(p1[ic]);
	 //Output to screen the integrals, errors and maximum weights
	 // 	 cout << "PHASE SPACE INTEGRALS:" << "\t" << p[ic] << "/" << p1[ic] <<"/" << sd[ic] <<"/" << sd1[ic]
	 //   <<"/" << bigwt[ic] << "\t" << xsec << endl;
}


  //*********************************************************************************
  /*Les Houches Accord event parameters:: PUP are the momenta, 
   ICOLUP is the colour code, MOTHUP is the mother information, ISTUP is 
   the status code. ISPINUP signifies spins. IDUP is the particle ID.
    NUP is the number of particles in the event. SCALUP is the scale of 
   the event.AQEDUP and AQCDUP are alphaem and alphas resp.XWGTUP is 
   the weight of the event=1 or -1 for unweighted events.*/
  
  double PUP[6][11],AQEDUP,AQCDUP,XWGTUP,SCALUP;
  int ICOLUP[3][11],MOTHUP[3][11],ISTUP[11],ISPINUP[11],IDUP[11],NUP, NUPP;
  NUP=0;
  NUPP=0;
  SCALUP = mH;
  AQEDUP=0.007297352;
  AQCDUP=alphas;
  IDUP[1]=11;
  IDUP[2]=-11;
  IDUP[3]=23;
  IDUP[4]=25;
  IDUP[5]=23;
  IDUP[6]=5;
  IDUP[7]=-5;
  IDUP[8]=11;
  IDUP[9]=-11;
  IDUP[10]=21;
  MOTHUP[1][1]=0;
  MOTHUP[1][2]=0;
  MOTHUP[2][1]=0;
  MOTHUP[2][2]=0;
  MOTHUP[1][3]=1;
  MOTHUP[2][3]=2;
  MOTHUP[1][4]=3;
  MOTHUP[2][4]=3;
  MOTHUP[1][5]=3;
  MOTHUP[2][5]=3;
  MOTHUP[1][6]=4;
  MOTHUP[2][6]=4;
  MOTHUP[1][7]=4;
  MOTHUP[2][7]=4;
  MOTHUP[1][8]=5;
  MOTHUP[2][8]=5;
  MOTHUP[1][9]=5;
  MOTHUP[2][9]=5;
  MOTHUP[1][10]=4;
  MOTHUP[2][10]=4;
  ICOLUP[1][1]=0;
  ICOLUP[2][1]=0;
  ICOLUP[1][2]=0;
  ICOLUP[2][2]=0;
  ICOLUP[1][3]=0;
  ICOLUP[2][3]=0;
  ICOLUP[1][4]=0;
  ICOLUP[2][4]=0;
  ICOLUP[1][5]=0;
  ICOLUP[2][5]=0;
  ICOLUP[1][8]=0;
  ICOLUP[2][8]=0;
  ICOLUP[1][9]=0;
  ICOLUP[2][9]=0;
  ISTUP[1]=-1;
  ISTUP[2]=-1;
  ISTUP[3]=2;
  ISTUP[4]=2;
  ISTUP[5]=2;
  ISTUP[6]=1;
  ISTUP[7]=1;
  ISTUP[8]=1;
  ISTUP[9]=1;
  ISTUP[10]=1;
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
  
  //Allocate phase space region
  double cst[3],bigwgt;
  bigwgt=0.;
  cst[0]=fabs(p[0])+fabs(p1[0]);
  cst[1]=fabs(p[1])+fabs(p1[1]);
  cst[2]=fabs(p[2])+fabs(p1[2]);
  cstt = cst[0]+cst[1]+cst[2];
  
  // Calculating the associated cross-section in picobarns...
  double Betaz, tcst;
  if (!iface){
  Betaz = sqrt((pow(emcm,2)-pow(Mz+mH,2))*(pow(emcm,2)-pow(Mz-mH,2)))/(pow(emcm,2)-pow(mH,2)+pow(Mz,2));
  double gammaz = 1./sqrt(1.-pow(Betaz,2));
  tcst = 389400000.*pow(GF*Mw,2)*Betaz*gammaz*pow(Mz,3)/(cos2thw*32.*pi*pow(emcm, 3))*pow((pow(emcm,2)+pow(Mz,2)-pow(mH,2))/(pow(emcm,2)-pow(Mz,2)),2)*(1.-4.*sin2thw+8.*pow(sin2thw,2))*(4.*(1.-2./3.*pow(Betaz,2))) ;
  }
  //=====================================================================================================================================
   ifstream indata;
   ofstream outdata, outdata2, outdata3;
   if (iface && rest) {
   outdata2.open ("MCHGiface.dat", ios::trunc); // output for b, \bar{b} and g momenta in Higgs rest frame
   outdata2 << "MC@NLO Higgs boson decay in rest frame" <<endl ;}

   if (iface && ifile) {
   indata.open (filename); // input file containing Higgs 5-momenta
   outdata2.open ("MCHGiface.dat", ios::trunc); // output for b, \bar{b} and g momenta boosted to the lab frame
   outdata2 << "MC@NLO Higgs boson decay in lab frame" <<endl ;}

   if (!iface) {
   outdata.open("eeZH.dat", ios::trunc);//output for associated production events
   outdata << "<LesHouchesEvents version =\"1.0\">" << endl; outdata << "<!--" << endl;
   outdata << "IPROC = " << "\t" << 23 << endl;
   outdata << "NORMFACTOR = " << "\t" << xsecabs/xsec << endl;
   outdata << "MC@NLO Higgs boson decay" <<endl ;
     outdata << "File generated with MCHG.cxx" << endl; outdata << "-->" << endl;
     outdata << "<init>" << endl;
     outdata << "\t11\t" << "-11\t" ;
     outdata << emcm/2. << "\t" << emcm/2.  << "\t" << "0 \t 0 \t 7\t 7 \t 3 \t 1" << endl;
     outdata << "\t" << tcst  << "\t" << 0.000000 << "\t1.00000 \t11" << endl;     outdata << "</init>" << endl;
   }

  cout << "Generating events..." << endl;


 for (int d = 0; d < nevg; d++) {
  int Jt, Jtt; //phase space identifiers
  Jtt = 0;
  if (random(seed) < cst[0]/cstt) {
    if (random(seed) < fabs(p[0])/cst[0]) {
      Jt=1; bigwgt=bigwt[0];} else
	{Jt=2;}} else
	  if (random(seed) < cst[1]/(cst[1]+cst[2])) {
	    if (random(seed) < fabs(p[1])/cst[1]) {
              Jt=3; bigwgt=bigwt[1];} else
		{Jt=4;}} else
		  {Jt=5;}
  bool labrest = false;
  if (iface && ifile && !rest) {
    indata  >> PUP[1][4] >> PUP[2][4] >> PUP[3][4] >> PUP[4][4] >>  mH ;
     if (PUP[4][4] == mH) {labrest=true;} }

  double phi6;
  if (!iface) {
  // Distributing the Z and H momenta 

  double cth3, Mtht;
  do {
    cth3= 2.*random(seed)-1.;
    Mtht = 2.*(1.-pow(Betaz,2))+pow(Betaz,2)*(1.-pow(cth3,2));
  } while (random(seed) > Mtht/(2.*(1.-pow(Betaz,2))+pow(Betaz,2)));
  phi6 = 2.*random(seed)*pi;
  double pZ;

 if (Betaz == 0.) {pZ = 0.;} else {pZ = sqrt(pow(0.5*(emcm+(pow(mH,2)-pow(Mz,2))/emcm),2)-pow(mH,2));}
  PUP[4][4]=0.5*(emcm+(pow(mH,2)-pow(Mz,2))/emcm);
  PUP[1][4]=-pZ*sin(acos(cth3))*cos(phi6);
  PUP[2][4]=-pZ*sin(acos(cth3))*sin(phi6);
  PUP[3][4]=-pZ*cth3;
  PUP[4][5]=emcm-PUP[4][4];
  PUP[1][5]=pZ*sin(acos(cth3))*cos(phi6);
  PUP[2][5]=pZ*sin(acos(cth3))*sin(phi6);
  PUP[3][5]=pZ*cth3;
  }
  if (iface && !rest) {phi6 = pi + acos(PUP[1][4]/sqrt(pow(PUP[1][4],2) +pow(PUP[2][4],2)));}
  //=================================================================================================================================
  //select flavour of process
  double mass[10];
  int ID = 5;
  mass[1]=0.;
  mass[2]=0.;
  mass[3]=emcm;
  mass[4]=mH;
  mass[5]=Mz;
  mass[8]=0.;
  mass[9]=0.;
  //quark IDs
  IDUP[6] = ID;
  IDUP[7] = -ID;
 
// Parameters for NLO cross-section calculation.
  double param[9];
  param[0] = pow(m,2)/pow(mH,2);
  param[1] = 2.*param[0];
  param[2] = sqrt(1.-4.*param[0]); // beta
  param[3] = (1.-pow(param[2],2))/param[2]*log((1.+param[2])/(1.-param[2]))+ (1.+pow(param[2],2))/param[2]*(4.*ddilog((1.-param[2])/(1.+param[2]))+2.*log((1.+param[2])/(2.*param[2]))*log((1.+param[2])/(1.-param[2]))+2.*ddilog((param[2]-1.)/(1.+param[2]))-log(2./(1.+param[2]))*log((1.+param[2])/(1.-param[2])))-3.*log(4./(1.-pow(param[2],2)))-4.*log(param[2])+1./(16.*pow(param[2],3))*(3.+18.*pow(param[2],2)+3.*pow(param[2],4))*log((1.+param[2])/(1.-param[2]))-3./8.*pow(param[2],2)+21./8.;
  param[4] = 0.5*(1.+param[2]);
  param[5] = (param[0]*param[4]+0.25*param[2]*pow(1+param[2],2))/(param[4]-param[0]);
  param[6] = 4.*alphas/(3.*pi); // colour factor
  param[7] = 3.*GF*pow(m,2)/(4.*sqrt(2)*pi)*mH*pow(param[2],3);
  param[8] = 1./(2.*pow(param[2],3));

 // phase space point x1,x2
  double x1,x2;
  double wt; //event weight
  wt = 0.;
  bool qg = true;
  double cth,pcm;
  cth=0.;pcm=0.;

   if (Jt==1||Jt==3) {
    double c13, c23, x3;
    do {
      // Random numbers 
   double r1=random(seed);
   double r3=random(seed);
   double r4=random(seed);
   double r2 =random(seed); 
    x1 = r1;  
    x2 = r2;
        if (Jt==1) {
	QQgsc(r1, r2, r3, r4, x1, x2, param, wgt, wgt1, resolve);
	wt = wgt;}
     if (Jt==3) {
        QQgh(r1, r2, r3, r4, x1, x2, param, wgt, wgt1, resolve);
	wt  = wgt;}
     x3=2.-x1-x2;
    c13=(x1*x3-2.*(1.-x2))/(sqrt(pow(x1,2)-4.*param[0])*x3);
    c23=(x2*x3-2.*(1.-x1))/(sqrt(pow(x2,2)-4.*param[0])*x3);
    } while (random(seed) > fabs(wt)/fabs(bigwgt) || isinf(c13) || isinf(c23));
    
    if (Jt==1){wt=-1.;}
    if (Jt==3){wt=1.;}    
    double cf,ct, ppt;
    double px1 = sqrt(pow(0.5*mH*x1,2)-pow(m,2));
    double px2 = sqrt(pow(0.5*mH*x2,2)-pow(m,2));
    double px3 = sqrt(pow(0.5*mH*x3,2));
    do {cth = 2.*random(seed)-1.;
    if (random(seed) < pow(x1,2)/(pow(x1,2)+pow(x2,2))) {
    cf = (px2*c23+px3)/px1; ct = (px1-px3*cf)/px2;}else
    {cf = (px1*c13+px3)/px2; ct = (px2-px3*cf)/px1; qg=false;}
     } while (ct > 1. || cf > 1.);
    double phi = 2.*pi*random(seed);
    double phi2 = 2.*pi*random(seed);
    if (qg) {
      pcm = px1;
      PUP[3][6]=px1*cth;
      double pt = sqrt(pow(px1,2)-pow(PUP[3][6],2));
      PUP[1][6]=pt*cos(phi);
      PUP[2][6]=pt*sin(phi);
      MOTHUP[1][10] = 4;
      MOTHUP[2][10] = 4;
      double ppl = px2*ct;
       if (random(seed) < 0.5) {
	ppt = sqrt (pow(px2,2)-pow(ppl,2));}else {ppt = -sqrt (pow(px2,2)-pow(ppl,2));}
      double ppt1 = ppt*cos(phi2);
      double ppt2 = ppt*sin(phi2);
      PUP[3][7]=-(ppl*cth-ppt1*sin(acos(cth)));
      PUP[1][7]=-((ppl*sin(acos(cth))+ppt1*cth)*cos(phi)-ppt2*sin(phi));
      PUP[2][7]=-((ppl*sin(acos(cth))+ppt1*cth)*sin(phi)+ppt2*cos(phi));} 
    else {
      pcm = px2;
      PUP[3][7]=px2*cth;
      double pt = sqrt(pow(px2,2)-pow(PUP[3][7],2));
      PUP[1][7]=pt*cos(phi);
      PUP[2][7]=pt*sin(phi);
      MOTHUP[1][10] = 4;
      MOTHUP[2][10] = 4;
      double ppl = px1*ct;
      if (random(seed) < 0.5) {
	ppt = sqrt (pow(px1,2)-pow(ppl,2));}else {ppt = -sqrt (pow(px1,2)-pow(ppl,2));}
      double ppt1 = ppt*cos(phi2);
      double ppt2 = ppt*sin(phi2);
      PUP[3][6]=-(ppl*cth-ppt1*sin(acos(cth)));
      PUP[1][6]=-((ppl*sin(acos(cth))+ppt1*cth)*cos(phi)-ppt2*sin(phi));
      PUP[2][6]=-((ppl*sin(acos(cth))+ppt1*cth)*sin(phi)+ppt2*cos(phi));}

      PUP[3][10]= -PUP[3][6] - PUP[3][7];
      PUP[1][10]= -PUP[1][6] - PUP[1][7];
      PUP[2][10]= -PUP[2][6] - PUP[2][7];
      PUP[4][6] = x1*mH/2.;
      PUP[4][7] = x2*mH/2.;
      PUP[4][10] = x3*mH/2.;
      //colour info
      ICOLUP[1][6] = 501;
      ICOLUP[2][6] = 0;
      ICOLUP[1][7] = 0;
      ICOLUP[2][7] = 502;
      ICOLUP[1][10] = 502;
      ICOLUP[2][10] = 501;
      //Number of particles
      NUP=10;
      NUPP = 4; // number in Higgs rest frame
      if (Jt == 1) { Jtt = 1;} // from soft collinear region 
      if (Jt == 3) { Jtt = 2;} // from hard non-collinear region
}
     

     if (Jt==2||Jt==4||Jt==5) {
       Jtt = 1;
     x1 = 1.;x2 = 1.; wt = 1.;
     pcm = sqrt(pow(0.5*mH,2)-pow(m,2));
     cth = 2.*random(seed)-1.;
     double phi = 2.*random(seed)*pi;
     PUP[3][6] = pcm*cth;
     PUP[1][6] = sqrt(pow(pcm,2)-pow(PUP[3][6],2))*cos(phi);
     PUP[2][6] = sqrt(pow(pcm,2)-pow(PUP[3][6],2))*sin(phi);  
     PUP[1][7] = -PUP[1][6];
     PUP[2][7] = -PUP[2][6];
     PUP[3][7] = -PUP[3][6];
     PUP[4][6] = mH/2.;
     PUP[4][7] = mH/2.;
     PUP[3][10]=0.;
     PUP[1][10]=0.;
     PUP[2][10]=0.;
     //colour info
      ICOLUP[1][6] = 501;
      ICOLUP[2][6] = 0;
      ICOLUP[1][7] = 0;
      ICOLUP[2][7] = 501;
    //Number of particles
      NUP=9;
      NUPP = 3;  // number in Higgs rest frame
      Jtt = 3;   // Born event
     }
      double alpha = 1.;

     //boost masses if needed. Need to boost to nominal gluon mass always.
     double mb = m;
     double Eq,Eqb,Eg;
     double pq = sqrt(pow(PUP[1][6],2)+pow(PUP[2][6],2)+pow(PUP[3][6],2));
     double pqb = sqrt(pow(PUP[1][7],2)+pow(PUP[2][7],2)+pow(PUP[3][7],2));
     double pg = sqrt(pow(PUP[1][10],2)+pow(PUP[2][10],2)+pow(PUP[3][10],2));
       for (int ix=1; ix < nit+1; ix++) {
       Eq = sqrt(pow(mb,2)+alpha*pow(pq,2));
       Eqb = sqrt(pow(mb,2)+alpha*pow(pqb,2));
       if (NUP == 9) {mg =0.;} else {mg =0.75;}
       Eg =  sqrt(pow(mg,2)+alpha*pow(pg,2));
       if (NUP==10) {
       alpha = alpha+(2.*(mH-Eq-Eqb-Eg))/(pow(pq,2)/Eq+pow(pqb,2)/Eqb+pow(pg,2)/Eg);}
       else {alpha = alpha+(2.*(mH-Eq-Eqb))/(pow(pq,2)/Eq+pow(pqb,2)/Eqb);}
       }
     PUP[1][6] =sqrt(alpha)*PUP[1][6];
     PUP[2][6] =sqrt(alpha)*PUP[2][6];
     PUP[3][6] =sqrt(alpha)*PUP[3][6];
     PUP[4][6] =Eq;
     PUP[1][7] =sqrt(alpha)*PUP[1][7];
     PUP[2][7] =sqrt(alpha)*PUP[2][7];
     PUP[3][7] =sqrt(alpha)*PUP[3][7];
     PUP[4][7] =Eqb;
     PUP[1][10] =sqrt(alpha)*PUP[1][10];
     PUP[2][10] =sqrt(alpha)*PUP[2][10];
     PUP[3][10] =sqrt(alpha)*PUP[3][10];
     PUP[4][10] =Eg;
     mass[6] = m;
     mass[7] = m;
     mass[10] = mg;
     if (iface && rest) {
       ii += 1;
      outdata2 << "<event>" << endl;
      outdata2 <<  NUPP   << "\t" << Jtt <<  "\t" << wt << endl;
      outdata2 << IDUP[4] << "\t" << 0. << "\t" << 0. << "\t" << 0. << "\t" << mass[4] << "\t" << mass[4] << endl;
      outdata2 << IDUP[6] << "\t" << PUP[1][6] << "\t" << PUP[2][6] << "\t" << PUP[3][6] << "\t" << PUP[4][6] << "\t" << mass[6] << endl;
      outdata2 << IDUP[7] << "\t" << PUP[1][7] << "\t" << PUP[2][7] << "\t" << PUP[3][7] << "\t" << PUP[4][7] << "\t" << mass[7] << endl;
     if (Jt == 1 || Jt == 3) {
      outdata2 << IDUP[10] << "\t" << PUP[1][10] << "\t" << PUP[2][10] << "\t" << PUP[3][10] << "\t" << PUP[4][10] << "\t" << mass[10] << endl;}
      outdata2 << "</event>" << endl;
      cout  << "Generated event : " << "\t" << ii << "\t" << "of" << "\t" << nevg << "\r" << flush;
}
     if (labrest) {
       ii+=1;
       outdata2 << "<event>" << endl;
       outdata2 <<  NUPP   << "\t" << Jtt <<  "\t" << wt << endl;
       outdata2 << IDUP[4] << "\t" << 0. << "\t" << 0. << "\t" << 0. << "\t" << mass[4] << "\t" << mass[4] << endl;
       outdata2 << IDUP[6] << "\t" << PUP[1][6] << "\t" << PUP[2][6] << "\t" << PUP[3][6] << "\t" << PUP[4][6] << "\t" << mass[6] << endl;
       outdata2 << IDUP[7] << "\t" << PUP[1][7] << "\t" << PUP[2][7] << "\t" << PUP[3][7] << "\t" << PUP[4][7] << "\t" << mass[7] << endl;
       if (Jt == 1 || Jt == 3) {
       outdata2 << IDUP[10] << "\t" << PUP[1][10] << "\t" << PUP[2][10] << "\t" << PUP[3][10] << "\t" << PUP[4][10] << "\t" << mass[10] << endl;}
       outdata2 << "</event>" << endl;
       cout  << "Generated event : " << "\t" << ii << "\t" << "of" << "\t" << nevg << "\r" << flush;
     }

        if (!labrest && ((!iface && Betaz != 0.) || (iface && !rest))) {
       // boost to lab frame
       double Eqq = PUP[4][6] ; double Eqqb = PUP[4][7] ; double Eggg= PUP[4][10];
       double vv2 = sqrt(pow(PUP[4][4],2)-pow(mH,2))/PUP[4][4];
       PUP[4][6]=(Eqq - vv2*PUP[3][6])/sqrt(1.-pow(vv2,2));
       PUP[3][6]=(PUP[3][6]-Eqq*vv2)/sqrt(1.-pow(vv2,2));
       PUP[4][7]=(Eqqb-vv2*PUP[3][7])/sqrt(1.-pow(vv2,2));
       PUP[3][7]=(PUP[3][7]-Eqqb*vv2)/sqrt(1.-pow(vv2,2));
       PUP[4][10]=(Eggg-vv2*PUP[3][10])/sqrt(1.-pow(vv2,2));
       PUP[3][10]=(PUP[3][10]-Eggg*vv2)/sqrt(1.-pow(vv2,2));
       // Rotate to beam axis.
       double costz = PUP[3][4]/(sqrt(pow(PUP[4][4],2)-pow(mH,2)));
       double sintz = sin(acos(costz));
       double costx = PUP[1][4]/(sqrt(pow(PUP[1][4],2)+pow(PUP[2][4],2)));
       double sintx = sin(acos(costx));
       costx = -cos(phi6);
       sintx = -sin(phi6);
       double qx = PUP[1][6];
       PUP[1][6] = -qx*costz+PUP[3][6]*sintz;
       PUP[3][6] = qx*sintz + PUP[3][6]*costz;
       qx = PUP[1][6];
       PUP[1][6]=qx*costx-PUP[2][6]*sintx;
       PUP[2][6]=qx*sintx+PUP[2][6]*costx;
       double qbx = PUP[1][7];
       PUP[1][7] = -qbx*costz+PUP[3][7]*sintz;
       PUP[3][7] = qbx*sintz + PUP[3][7]*costz;
       qbx = PUP[1][7];
       PUP[1][7]=qbx*costx-PUP[2][7]*sintx;
       PUP[2][7]=qbx*sintx+PUP[2][7]*costx;
       double gx = PUP[1][10];
       PUP[1][10] = -gx*costz+PUP[3][10]*sintz;
       PUP[3][10] = gx*sintz + PUP[3][10]*costz;
       gx = PUP[1][10];
       PUP[1][10]=gx*costx-PUP[2][10]*sintx;
       PUP[2][10]=gx*sintx+PUP[2][10]*costx;
       PUP[1][6] = -PUP[1][6];
       PUP[2][6] = -PUP[2][6];
       PUP[3][6] = -PUP[3][6];
       PUP[1][7] = -PUP[1][7];
       PUP[2][7] = -PUP[2][7];
       PUP[3][7] = -PUP[3][7];
       PUP[1][10] = -PUP[1][10];
       PUP[2][10] = -PUP[2][10];
       PUP[3][10] = -PUP[3][10];
       if (iface && !rest) {
	 ii+=1;
	   outdata2 << "<event>" << endl;
	   outdata2 <<  NUPP   << "\t" << Jtt <<  "\t" << wt << endl;
	   outdata2 << IDUP[4] << "\t" << PUP[1][4] << "\t" << PUP[2][4] << "\t" << PUP[3][4] << "\t" << PUP[4][4] << "\t" << mass[4] << endl;
	  outdata2 << IDUP[6] << "\t" << PUP[1][6] << "\t" << PUP[2][6] << "\t" << PUP[3][6] << "\t" << PUP[4][6] << "\t" << mass[6] << endl;
	   outdata2 << IDUP[7] << "\t" << PUP[1][7] << "\t" << PUP[2][7] << "\t" << PUP[3][7] << "\t" << PUP[4][7] << "\t" << mass[7] << endl;
	  if (Jt == 1 || Jt == 3) {
	  outdata2 << IDUP[10] << "\t" << PUP[1][10] << "\t" << PUP[2][10] << "\t" << PUP[3][10] << "\t" << PUP[4][10] << "\t" << mass[10] << endl;}
	  outdata2 << "</event>" << endl;
	 // cout << "here" << endl;
	  cout  << "Generated event : " << "\t" << ii << "\t" << "of" << "\t" << nevg << "\r" << flush;
       }
	}

     double CLf, CRf;
     CLf = 0.; CRf = 0.;
      //spin information (not needed yet!)
     double C[9], M[9];
     C[1] = pow(CRe,2)*(chi2/(pow((mH*k),2)))*pow(gz2,2);
     C[2] = pow(CLe,2)*(chi2/(pow((mH*k),2)))*pow(gz2,2);
     C[3] = CLf*(mH/2.-pcm)+CRf*(mH/2.+pcm);
     C[4] = CLf*(mH/2.+pcm)+CLf*(mH/2.-pcm);
     C[5] = (CLf+CRf)*m;
     C[6] = 2*e2*Ch[ID]*CRe*chi1*gz2/(k*mH);
     C[7] = 2*e2*Ch[ID]*CLe*chi1*gz2/(k*mH);
     C[8] = e2*Ch[ID];
     M[1] = (C[1]*pow(C[3],2)-C[6]*C[3]+pow(C[8],2))*pow(1+cth,2);
     M[5] = (C[2]*pow(C[4],2)-C[7]*C[4]+pow(C[8],2))*pow(1+cth,2);
     M[2] = (C[1]*pow(C[4],2)-C[6]*C[4]+pow(C[8],2))*pow(1+cth,2);
     M[6] = (C[2]*pow(C[3],2)-C[7]*C[3]+pow(C[8],2))*pow(1+cth,2);
     M[3] = (C[1]*pow(C[5],2)-2.*m/mH*C[6]*C[5]+pow(C[8],2)*4.*param[0])*(1.-pow(cth,2));
     M[7] = (C[2]*pow(C[5],2)-2.*m/mH*C[7]*C[5]+pow(C[8],2)*4.*param[0])*(1.-pow(cth,2));
     M[4] = M[3]; 
     M[8] = M[7];
     double Mud, Mdu;
     Mud = 0.; Mdu =0.;
     for (int ix=1; ix < 5; ix++) {Mud+=M[ix];}
     for (int ix=5; ix < 9; ix++) {Mud+=M[ix];}
     if (random(seed) < Mud/(Mud+Mdu)) {
       ISPINUP[1]=+1;
       ISPINUP[2]=-1;
       ISPINUP[3]=+1;
     if (random(seed) < (M[1]+M[3])/Mud) {
       ISPINUP[4] = +1;
     if (random(seed) < M[1]/(M[1]+M[3])) {
       ISPINUP[5] = -1;}else {ISPINUP[5] = +1;}} else {
	 ISPINUP[4] = -1;
	 if (random(seed) < M[2]/(M[2]+M[4])) {
	   ISPINUP[5] = +1;} else {ISPINUP[5] = -1;} }} else {
       ISPINUP[1]=-1;
       ISPINUP[2]=+1;
       ISPINUP[3]=-1;
     if (random(seed) < (M[6]+M[7])/Mdu) {
       ISPINUP[4] = +1;
     if (random(seed) < M[6]/(M[6]+M[7])) {
       ISPINUP[5] = -1;}else {ISPINUP[5] = +1;}} else {
	 ISPINUP[4] = -1;
	 if (random(seed) < M[5]/(M[5]+M[8])) {
	   ISPINUP[5] = +1;} else {ISPINUP[5] = -1;} }}
     if (qg) {ISPINUP[6] = ISPINUP[4];} else {ISPINUP[6] = ISPINUP[5];}

     if (!iface) {
     // Z boson rest frame
     
     double cth2,cofact, maxcofact;
     IDUP[3]=23;
     IDUP[5]=23;
     IDUP[8]=11;
     IDUP[9]=-11;
     do {cth2=2.*random(seed)-1.;
       cofact = 1.+pow(cth2,2);
       maxcofact = 2.;} while(random(seed) > cofact/maxcofact);
     double phi5 = 2.*random(seed)*pi;
     PUP[3][9] = Mz*cth2/2.;
     PUP[1][9] = sqrt(pow(Mz/2.,2)-pow(PUP[3][4],2))*cos(phi5);
     PUP[2][9] = sqrt(pow(Mz/2.,2)-pow(PUP[3][4],2))*sin(phi5);
     PUP[1][8] = -PUP[1][9];
     PUP[2][8] = -PUP[2][9];
     PUP[3][8] = -PUP[3][9];
     PUP[4][8] = Mz/2.;
     PUP[4][9] = Mz/2.;
     if (Betaz != 0.) {
       // boost to lab
       double vv2 = sqrt(pow(PUP[4][5],2)-pow(Mz,2))/PUP[4][5];
       PUP[4][9]=(Mz/2.-vv2*PUP[3][9])/sqrt(1.-pow(vv2,2));
       PUP[3][9]=(PUP[3][9]-Mz*vv2/2.)/sqrt(1.-pow(vv2,2));
       PUP[4][8]=(Mz/2.-vv2*PUP[3][8])/sqrt(1.-pow(vv2,2));
       PUP[3][8]=(PUP[3][8]-Mz*vv2/2.)/sqrt(1.-pow(vv2,2));
       double ktt1=-sqrt(pow(PUP[4][9],2)-pow(PUP[3][9],2));
       double ktt2=sqrt(pow(PUP[4][8],2)-pow(PUP[3][8],2));
       double kt = sqrt(pow(PUP[1][5],2)+pow(PUP[2][5],2));
       double tht=atan2(kt,PUP[3][5]);
       double pe1=PUP[4][9];
       double pe2=PUP[4][8];
       double pz1=PUP[3][9];
       double pz2=PUP[3][8];
       PUP[4][9]=pe2;
       PUP[4][8]=pe1;
       PUP[3][9]=-(pz2*cos(tht)+ktt2*sin(tht));
       PUP[3][8]=-(pz1*cos(tht)+ktt1*sin(tht));
       PUP[1][9]=-(pz2*sin(tht)-ktt2*cos(tht))*cos(phi6);
       PUP[2][9]=-(pz2*sin(tht)-ktt2*cos(tht))*sin(phi6);
       PUP[1][8]=-(pz1*sin(tht)-ktt1*cos(tht))*cos(phi6);
       PUP[2][8]=-(pz1*sin(tht)-ktt1*cos(tht))*sin(phi6);
       if (random(seed) < 0.5) {
	 PUP[4][9]=pe1;
	 PUP[4][8]=pe2;
	 PUP[3][9]=-(pz1*cos(tht)+ktt1*sin(tht));
	 PUP[3][8]=-(pz2*cos(tht)+ktt2*sin(tht));
	 PUP[1][9]=-(pz1*sin(tht)-ktt1*cos(tht))*cos(phi6);
	 PUP[2][9]=-(pz1*sin(tht)-ktt1*cos(tht))*sin(phi6);
	 PUP[1][8]=-(pz2*sin(tht)-ktt2*cos(tht))*cos(phi6);
	 PUP[2][8]=-(pz2*sin(tht)-ktt2*cos(tht))*sin(phi6);
       } }
     //weight info
     XWGTUP = wt;
     ii+=1;
     double sum1, sum2, sum3, sum4;
     sum1=0.;sum2=0.;sum3=0.;sum4=0.;
     
     outdata << "<event>" << endl;
     outdata << NUP <<"\t" << "11" <<"\t"<<XWGTUP<<"\t" << SCALUP << "\t" << AQEDUP << "\t" << AQCDUP << endl;
     for (int ja = 1; ja < NUP+1; ja++) {
       outdata << IDUP[ja] << "\t" << ISTUP[ja] << "\t" << MOTHUP[1][ja] <<"\t" << MOTHUP[2][ja] << "\t" << ICOLUP[1][ja] << "\t" << ICOLUP[2][ja] <<"\t" <<setprecision (9)<< PUP[1][ja] << "\t" << PUP[2][ja]<<"\t"<<PUP[3][ja]<<"\t"<<PUP[4][ja]<<"\t" << mass[ja] <<"\t"<< "0" << "\t" <<"9" <<endl;
       sum1+=PUP[4][ja];sum2+=PUP[1][ja];sum3+=PUP[2][ja];sum4+=PUP[3][ja];
       cout  << "Generated event : " << "\t" << ii << "\t" << "of" << "\t" << nevg << "\r" << flush;  
}
     outdata << "</event>" << endl;}} 
 if (!iface) {
 outdata << "</LesHouchesEvents>" << endl;
 cout << endl;}
  return 0;
}
/********************************************************************************/
// In soft and collinear phase space region ?

void QQgsc(double r1, double r2, double r3,double r4, double &rx, double &ry, double param[9], double &rwgt, double &rwgt1, bool resolve) { 
  
  //Initialize weights
  rwgt = 0.; rwgt1 = 0.;

   /*Remove 2-jet events */
  double y1 = 1.-r1;
  double y2 = 1.-r2;
  if ((1.-y1)==1 && (1.-y2)==1.) {return;}
  bool soft = false;
  int region = 0; // 1 for jet region cont. x=1; 2 for jet region cont. y=1
  double fac = 1.; // weight factor for soft regions

     if (y1 < 0.25 && y2 <0.25) { 
     soft = true;
     if (y2 < y1) {
     y1 = 0.25*r3;
     y2 = r4*y1;
      fac = 2.*r3;
    } else {
   	y2 = 0.25*r4;
   	y1 = r3*y2;
   fac = 2.*r4;
         }
        }

  /* Check if outside Dalitz boundary */

  if (y1*y2*(1.-y1-y2) < param[0]*pow((y1+y2),2)) {return;}
      
  /*Inside jet 1? */
  double w1 = 1.-y1;
  double w2 = 1.-y2;
  double tk = param[4];
  double d1 = 0.25-(1.-w2)/tk;
  double z, d2 , xx1, xx2;
  if (d1 > 0.) {
    z = 0.5+sqrt(d1);
    d2 = pow(w2,2)-4.*param[0];
    if (d2 > 0.) {
      d2 = sqrt(d2);
      double rr = 1.+param[0]/(1.-w2+param[0]);
      xx1 = (rr*(2.-w2)+(2.*z-rr)*d2)/2.;
      xx2 = (rr*(2.-w2)+(2.*(1.-z)-rr)*d2)/2.;
      if (xx2 < w1 && w1 < xx1) {
	region = 1;}
    }
  }
  if (region == 0) {
      w1 = 1.-y2;
      w2 = 1.-y1;
      tk = param[5];
      d1 = 0.25-(1.-w2)/tk;
     if (d1 > 0.) {
    z = 0.5+sqrt(d1);
    d2 = pow(w2,2)-4.*param[0];
    if (d2 > 0.) {
      d2 = sqrt(d2);
      double rr = 1.+param[0]/(1.-w2+param[0]);
      xx1 = (rr*(2.-w2)+(2.*z-rr)*d2)/2.;
      xx2 = (rr*(2.-w2)+(2.*(1.-z)-rr)*d2)/2.;
      if (xx2 < w1 && w1 < xx1) {
	region = 2;} else {return;}
    } else {return;}
     } else  {return;} }	
  // fac = 1.;
     /* Calculate weights for jet regions */
     double r, zz;
     rx = 1.-y1;
     ry = 1.-y2;
   
    if (region == 1) {
       r = 0.5*(1.+param[0]/(1.+param[0]-ry));
       zz = r+(rx-(2.-ry)*r)/sqrt(pow(ry,2)-2.*param[1]);
       rwgt = ((pow(1-rx,2) + pow(1-ry,2))/((1-rx)*(1-ry))+2.*(1.-param[1])*(rx+ry-1.-2.*param[1])/((1.-rx)*(1.-ry))+2.*param[1]*(1./(1.-rx)+1./(1.-ry))-param[1]*(1.-2.*param[1])*(1./pow(1-rx,2)+1./pow(1.-ry,2))+2.)*param[8];
	rwgt = rwgt-((((1.+pow(zz,2))/(1.-zz))-param[1]/y2)/(y2*sqrt(pow(ry,2)-2.*param[1])))/2.;
	rwgt = rwgt*param[6]*fac;
	     } else if (region == 2) {
       r = 0.5*(1.+param[0]/(1.+param[0]-rx));
       zz = r+(ry-(2.-rx)*r)/sqrt(pow(rx,2)-2.*param[1]);
        rwgt = ((pow(1-rx,2) + pow(1-ry,2))/((1-rx)*(1-ry))+2.*(1.-param[1])*(rx+ry-1.-2.*param[1])/((1.-rx)*(1.-ry))+2.*param[1]*(1./(1.-rx)+1./(1.-ry))-param[1]*(1.-2.*param[1])*(1./pow(1-rx,2)+1./pow(1.-ry,2))+2.)*param[8];
	 rwgt = rwgt-((((1.+pow(zz,2))/(1.-zz))-param[1]/y1)/(y1*sqrt(pow(rx,2)-2.*param[1])))/2.;
	 rwgt = rwgt*param[6]*fac;
	     }
         rwgt1 = (2.+2.*param[3]*param[6])*fac-rwgt;
       /* Implement soft resolution cut */
     if (resolve) {
           if (rx+ry > 1.9999) {return;}}
}
 /*********************************************************************/
 // In hard region?

void QQgh(double r1, double r2, double r3,double r4, double &rx, double &ry, double param[9], double &rwgt,double &rwgt1, bool resolve) { 
 
  //Initialize weights and phase space points
  rwgt = 0.;rwgt1 = 0.;

   /*Remove 2-jet events */
  double y1 = 1.-r1;
  double y2 = 1.-r2;
  if ((1.-y1)==1. && (1.-y2)==1.) {return;}
  bool soft = false;
  double fac = 1.; // weight factor for soft regions 

// Mappings for soft and collinear regions 
  if (y1 < 0.25){
     if (y2 < 0.25) { 
  soft = true;
   if (y1 < y2) {
     y1 = 0.25-y1;
     y2 = y1*(1.5-2.*y2);
      
  } else {
  	y2 = 0.25-y2;
  	y1 = y2*(1.5-2.*y1);
  }} else if (y2 < y1+2.*pow(y1,2)) {return;}
   } else if (y2 < 0.25) { if (y1 < y2+2.*pow(y2,2)) {return;} }
  

   if (soft) {if (y1 < y2) {fac = 2.*y1;} else {fac = 2.*y2;} }
       if (!soft) {
  
    if (y1 > 0.375 && y2 < 0.25) {
    if (y1 > 1.-2.5*y2 && y1 < 1.-y2) {
    y2 = r4*0.25;
    y1 = 1.-(1.+1.5*r3)*y2;
    fac = 2.*r4;
       } else {return;} }
  
    if (y2 > 0.375 && y1 < 0.25) {
     if (y2 > 1.-2.5*y1 && y2 < 1.-y1) {
     y1 = r3*0.25;
     y2 = 1.-(1.+1.5*r4)*y1;
      fac = 2.*r3;
      } else {return;}}}
   
//    Check if outside half-box */
     if (2.-y1-y2 < 1.) {return;}

  

 /*Inside jet 1? */
  double w1 = 1.-y1;
  double w2 = 1.-y2;
  double tk = param[4];
  double d1 = 0.25-(1.-w2)/tk;
  double z,d2,xx;

  
  if (d1 > 0.) {
    z = 0.5+sqrt(d1);
    d2 = pow(w2,2)-4.*param[0];
    if (d2 > 0.) {
      d2 = sqrt(d2);
      double rr = 1.+param[0]/(1.-w2+param[0]);
      xx = (rr*(2.-w2)+(2.*z-rr)*d2)/2.;
      if (w1 < xx) {
	z = 1.-z;
	xx = (rr*(2.-w2)+(2.*z-rr)*d2)/2.;
	if (w1 > xx) {return;}}}}

  /*Inside jet 2? */ 
  
      w1 = 1.-y2;
      w2 = 1.-y1;
      tk = param[5];
      d1 = 0.25-(1.-w2)/tk;
     if (d1 > 0.) {
    z = 0.5+sqrt(d1);
    d2 = pow(w2,2)-4.*param[0];
    if (d2 > 0.) {
      d2 = sqrt(d2);
      double rr = 1.+param[0]/(1.-w2+param[0]);
      xx = (rr*(2.-w2)+(2.*z-rr)*d2)/2.;
      if (w1 < xx) {
	z = 1.-z;
	xx = (rr*(2.-w2)+(2.*z-rr)*d2)/2.;
        if (w1 > xx) {return;}}}}
      /* Outside Dalitz plot? */
 if (y1*y2*(1.-y1-y2) < param[0]*pow((y1+y2),2)) { return;}

     /* Calculate weights for dead region */
          rx = 1.-y1;
          ry = 1.-y2;
	  rwgt = (pow(1-rx,2) + pow(1-ry,2))/((1-rx)*(1-ry))+2.*(1.-param[1])*(rx+ry-1.-2.*param[1])/((1.-rx)*(1.-ry))+2.*param[1]*(1./(1.-rx)+1./(1.-ry))-param[1]*(1.-2.*param[1])*(1./pow(1-rx,2)+1./pow(1.-ry,2))+2.;
         rwgt = rwgt*param[6]*param[8]*fac;
 	 rwgt1 = (2.+2.*param[3]*param[6])*fac-rwgt;
	/* Implement soft resolution cut */
     if (resolve) {
           if (rx+ry > 1.9999) {return;}}
}
 /*********************************************************************/
 // In gap between phase space for massive partons and full half-triangle? 

void QQgsp(double r1, double r2, double r3,double r4, double &rx, double &ry, double param[8], double &rwgt,double &rwgt1, bool resolve) {

  //Initialize weights and phase space points
  rwgt = 0.;rwgt1 = 0.;

   /*Remove 2-jet events */
  double y1 = 1.-r1;
  double y2 = 1.-r2;
  if ((1.-y1)==1. && (1.-y2)==1.) {return;}
  /* Check if outside half-box */
  if (2.-y1-y2 <= 1.) {return;}

/* Inside Dalitz plot? */
 if (y1*y2*(1.-y1-y2) > param[0]*pow((y1+y2),2)) {return;}
 /* Calculate weights */
 rwgt1 = 2.+2.*param[3]*param[6];}
     
 /*********************************************************************/
// Random number generator.

  double random (int &rseed) {
    int M = 2147483647;
    int A = 16807;
    int Q = 127773;
    int R = 2836;
    double MINV = 0.46566128752458e-09;
    int HI = rseed/Q;
    int LO = rseed % Q;
    rseed = A*LO - R*HI;
      if (rseed < 0) {rseed = rseed + M;}
    return rseed*MINV;}
/*********************************************************************/
double ddilog(double x)
{
  // The DiLogarithm function
  // Code translated by R.Brun from CERNLIB DILOG function C332

  const double hf  = 0.5;
  const double pi  = 3.14152654;
  const double pi2 = pi*pi;
  const double pi3 = pi2/3;
  const double pi6 = pi2/6;
  const double pi12 = pi2/12;
  const double c[20] = {0.42996693560813697, 0.40975987533077105,
			  -0.01858843665014592, 0.00145751084062268,-0.00014304184442340,
			  0.00001588415541880,-0.00000190784959387, 0.00000024195180854,
			  -0.00000003193341274, 0.00000000434545063,-0.00000000060578480,
			  0.00000000008612098,-0.00000000001244332, 0.00000000000182256,
			  -0.00000000000027007, 0.00000000000004042,-0.00000000000000610,
			  0.00000000000000093,-0.00000000000000014, 0.00000000000000002};

  double t,h,y,s,a,alfa,b1,b2,b0;

  if (x == 1) {
    h = pi6;
  } else if (x == -1) {
    h = -pi12;
  } else {
    t = -x;
    if (t <= -2) {
      y = -1/(1+t);
      s = 1;
      b1= log(-t);
      b2= log(1+1/t);
      a = -pi3+hf*(b1*b1-b2*b2);
    } else if (t < -1) {
      y = -1-t;
      s = -1;
      a = log(-t);
      a = -pi6+a*(a+log(1+1/t));
    } else if (t <= -0.5) {
      y = -(1+t)/t;
      s = 1;
      a = log(-t);
      a = -pi6+a*(-hf*a+log(1+t));
    } else if (t < 0) {
      y = -t/(1+t);
      s = -1;
      b1= log(1+t);
      a = hf*b1*b1;
    } else if (t <= 1) {
      y = t;
      s = 1;
      a = 0;
    } else {
      y = 1/t;
      s = -1;
      b1= log(t);
      a = pi6+hf*b1*b1;
    }
    h    = y+y-1;
    alfa = h+h;
    b1   = 0;
    b2   = 0;
    for (int i=19;i>=0;i--){
      b0 = c[i] + alfa*b1-b2;
      b2 = b1;
      b1 = b0;
    }
    h = -(s*(b0-h*b2)+a);
  }
  return h;
}




 /*********************************************************************/
 /*********************************************************************/
 /*********************************************************************/  
