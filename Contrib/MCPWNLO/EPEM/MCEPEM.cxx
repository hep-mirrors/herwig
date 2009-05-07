#include <cstdlib> 
#include <ctime> 
#include <iostream>
#include <fstream>
#include <string>
//#include "LHAPDFWrap.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "MCEPEM_INPUTS.h"
#include <cmath>
#include <iomanip>
using namespace std;
// Checks if phase space point is in soft and collinear region? 
void QQgsc(double r3,double r4, double &rx, double &ry, double param[8], double &rwgt, double &rwgtax, double &rwgt1, double &rwgtax1, bool resolve);
// Checks if phase space point is in hard region? 
void QQgh(double r3,double r4, double &rx, double &ry, double param[8], double &rwgt, double &rwgtax, double &rwgt1, double &rwgtax1, bool resolve);
// Checks if phase space point is in gap between phase space for massive partons and full half-triangle? 
void QQgsp(double r3,double r4, double &rx, double &ry, double param[8], double &rwgt, double &rwgtax, double &rwgt1, double &rwgtax1, bool resolve);
// Random number generator
double random (int &rseed);

// Main program 

int main() {
  Input user;
  int nev = user.nevint();
  int nevg = user.nevgen();
  int nit = user.it();
  double emcm = user.cme();
  int seed = user.rseed();
  int nf = user.nf();
  bool massiveME = user.massiveME();
  bool resolve = user.resolve();
  bool boost = user.bst();
  if (massiveME) {boost=false;}
  int ii=0; //event count
  // Constants 
  double pi = 3.14152654;
  double GF = 0.0000116639;
  double Mz = user.Mz();  
  double Yz = 2.486;
  double alphaem = 0.007297352;
  double alphas = user.alphasmz();
  //Lepton charge & coupling constants
  double Che = -1;
  double sin2thw = 0.2312;
  double Ae = -0.5;
  double Ve = Ae-2.*Che*sin2thw;
  double AV2 = pow(Ae,2)+pow(Ve,2);
  double CLe = (Ae+Ve)/2.;
  double CRe = (Ve-Ae)/2.;
  double gz2 = sqrt(32.)*GF*pow(Mz,2);
  double e2 = 4.*pi*alphaem;
  double k = sqrt(2.)*GF*pow(Mz,2)/e2;
  double chi1 = k*pow(emcm,2)*(pow(emcm,2)-pow(Mz,2))/(pow((pow(emcm,2)-pow(Mz,2)),2)+pow(Yz*Mz,2));
  double chi2 = pow(k,2)*pow(emcm,4)/(pow((pow(emcm,2)-pow(Mz,2)),2)+pow(Yz*Mz,2));
  //quark masses ,velocities and couplings. NB:gluon mass[6] =0.75 
 
  double mq[7],beta[6],Ch[6],Aq[6],Vq[6],cs[6],sigV[6],sigA[6],sigVA[6],css;
  mq[1] = 0.325; mq[2] = 0.325; mq[3] = 0.5; mq[4] = 1.6; mq[5] = 5.; mq[6] = 0.75;  

    Ch[1]=-1./3.; Ch[2]=2./3.; Ch[3]=Ch[1]; Ch[4]=Ch[2]; Ch[5]=Ch[1];
    Aq[1]=-0.5; Aq[2]=0.5; Aq[3]=Aq[1]; Aq[4]=Aq[2]; Aq[5]=Aq[1];

/* w, w1 are the integrands. p,p1 are the corresponding integrals calculated from sum,sum1 & sd, sd1  are the corresponding errors 
     calculated via ssq, ssq1. bigwt is the maximum value of the integrand over the phase space (will be used for event generation) */
    
  double wgt, wgtax, wgt1,wgtax1;
  double w[3],w1[3],bigwt[3],sum[3],sum1[3],ssq[3],ssq1[3],sd[3],sd1[3],p[3],p1[3];

  sum1[0]=0.;sum1[1]=0.;sum[2]=0.;
  sum1[0]=0.;sum1[1]=0.;sum1[2]=0.;
  ssq[0]=0.;ssq[1]=0.;ssq[2]=0.;
  ssq1[0]=0.;ssq1[1]=0.;ssq1[2]=0.;
  bigwt[0]=0.;bigwt[1]=0.;bigwt[2]=0.;

  cout << "Integrating..." << endl;

  for (int d = 0; d < nev; d++) {
  //Calculate  Born cross-section

  //initialize cross-sections...
  cs[1]=0.; cs[2]=0.; cs[3]=0.; cs[4]=0.; cs[5]=0.; css=0.;
  for (int ix=1; ix < nf+1 ; ix++) {
    Vq[ix] = Aq[ix]-(2.*Ch[ix]*sin2thw);
    beta[ix] = sqrt(1.-4.*pow(mq[ix]/emcm,2));
    sigV[ix] = pow(Ch[ix],2)-2.*Ch[ix]*Ve*Vq[ix]*chi1+AV2*pow(Vq[ix],2)*chi2;
    sigA[ix] = AV2*pow(Aq[ix],2)*chi2;
    sigVA[ix]= -2.*Ch[ix]*Ae*Aq[ix]*chi1+4.*Ae*Ve*Aq[ix]*Vq[ix]*chi2;
    cs[ix] = beta[ix]+0.5*pow(beta[ix],2)+0.5*(1.-pow(beta[ix],2))*beta[ix]*sigV[ix]+pow(beta[ix],3)*sigA[ix];
    css+=cs[ix];
}


  //select flavour of process
  double m;
  int ID;
  if (random(seed) < cs[1]/css) {ID = 1; } else {
    if (random(seed) < cs[2]/(css-cs[1])) {ID = 2;} else {
      if (random(seed) < cs[3]/(css-cs[1]-cs[2])) {ID = 3; } else {
	if (random(seed) < cs[4]/(cs[4]+cs[5])) {ID = 4; } else {
	  ID = 5; m = mq[5];}}}}

   m=mq[ID];
  if (!massiveME) {m = 0.;}
 
// Parameters for NLO cross-section calculation.
  double param[8];
  param[0] = pow(m,2)/pow(emcm,2);
  param[1] = 1.+12.*param[0]; // R = 1+param[1]*alphas/pi for vector current
  param[2] = 1.-22.*param[0]; // R = 1+param[2]*alphas/pi for axial vector current
  param[3] = 2.*param[0];
  param[4] = sqrt(1.-4.*param[0]); // beta
  param[5] = 0.5*(1.+param[4]);
  param[6] = (param[0]*param[5]+0.25*param[4]*pow(1+param[4],2))/(param[5]-param[0]);
  param[7] = 2.*alphas/(3.*pi); // colour factor
 
  // Random numbers and phase space points x and y.
  double r1=random(seed);
  double r2=random(seed);
  double r3=random(seed);
  double r4=random(seed);
  double x,y;
  x = r1;
  y = r2;
  QQgsc(r3, r4, x, y, param, wgt, wgtax, wgt1, wgtax1, resolve);
  w[0] = (wgt*sigV[ID]+wgtax*sigA[ID])/2.; w1[0]=(wgt1*sigV[ID]+wgtax1*sigA[ID])/2.;
  QQgh(r3, r4, x, y, param, wgt, wgtax, wgt1, wgtax1, resolve);
  w[1] = (wgt*sigV[ID]+wgtax*sigA[ID])/2.; w1[1]=(wgt1*sigV[ID]+wgtax1*sigA[ID])/2.;
  QQgsp(r3, r4, x, y, param, wgt, wgtax, wgt1, wgtax1, resolve);
  w[2] = (wgt*sigV[ID]+wgtax*sigA[ID])/2.; w1[2]=(wgt1*sigV[ID]+wgtax1*sigA[ID])/2.;
  if (!massiveME) {w[2]=0; w1[2]=0.;}
  
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
//               <<"/" << bigwt[ic] << "\t" << xsec << endl;
}


  //*********************************************************************************
  /*Les Houches Accord event parameters:: PUP are the momenta, 
   ICOLUP is the colour code, MOTHUP is the mother information, ISTUP is 
   the status code. ISPINUP signifies spins. IDUP is the particle ID.
    NUP is the number of particles in the event. SCALUP is the scale of 
   the event.AQEDUP and AQCDUP are alphaem and alphas resp.XWGTUP is 
   the weight of the event=1 or -1 for unweighted events.*/
  
  double PUP[6][7],AQEDUP,AQCDUP,XWGTUP,SCALUP;
  int ICOLUP[3][7],MOTHUP[3][7],ISTUP[7],ISPINUP[7],IDUP[7],NUP;
  NUP=0;
  SCALUP = emcm;
  AQEDUP=0.007297352;
  AQCDUP=alphas;
  IDUP[1]=11;
  IDUP[2]=-11;
  IDUP[3]=23;
  IDUP[6]=21;
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
  MOTHUP[1][6]=3;
  MOTHUP[2][6]=3;
  ICOLUP[1][1]=0;
  ICOLUP[2][1]=0;
  ICOLUP[1][2]=0;
  ICOLUP[2][2]=0;
  ICOLUP[1][3]=0;
  ICOLUP[2][3]=0;
  ISTUP[1]=-1;
  ISTUP[2]=-1;
  ISTUP[3]=2;
  ISTUP[4]=1;
  ISTUP[5]=1;
  ISTUP[6]=1;
  //Allocate phase space region
  double cst[3],bigwgt;
  bigwgt=0.;
  cst[0]=fabs(p[0])+fabs(p1[0]);
  cst[1]=fabs(p[1])+fabs(p1[1]);
  cst[2]=fabs(p[2])+fabs(p1[2]);
  cstt = cst[0]+cst[1]+cst[2];
   
   ofstream outdata;
   outdata.open("MCEPEM.dat", ios::trunc);//output for events


  outdata << "<LesHouchesEvents version =\"1.0\">" << endl; outdata << "<!--" << endl;
   outdata << "NORMFACTOR = " << "\t" << xsecabs/xsec << endl;
   outdata << "MC@NLO Z/gamma boson production from e+e- annihilation at LEP" <<endl ;
     outdata << "File generated with MCEPEM.cxx" << endl; outdata << "-->" << endl;
     outdata << "<init>" << endl;
     outdata << "\t11\t" << "-11\t" ;
     outdata << user.cme()/2. << "\t" << user.cme()/2.  << "\t" << "0 \t 0 \t 7\t 7 \t 3 \t 1" << endl;
     outdata << "\t" << 30590.  << "\t" << 0.000000 << "\t1.00000 \t11" << endl;     outdata << "</init>" << endl;


  cout << "Generating events..." << endl;


 for (int d = 0; d < nevg; d++) {
  int Jt; //phase space identifier
  if (random(seed) < cst[0]/cstt) {
    if (random(seed) < fabs(p[0])/cst[0]) {
      Jt=1; bigwgt=bigwt[0];} else
	{Jt=2;}} else
	  if (random(seed) < cst[1]/(cst[1]+cst[2])) {
	    if (random(seed) < fabs(p[1])/cst[1]) {
              Jt=3; bigwgt=bigwt[1];} else
		{Jt=4;}} else
		  {Jt=5;}
 
  //select flavour of process
  double m,mg,mass[8];
  int ID;
 if (random(seed) < cs[1]/css) {ID = 1; } else {
    if (random(seed) < cs[2]/(css-cs[1])) {ID = 2;} else {
      if (random(seed) < cs[3]/(css-cs[1]-cs[2])) {ID = 3;} else {
	if (random(seed) < cs[4]/(cs[4]+cs[5])) {ID = 4;} else {
	  ID = 5;}}}}
 if (massiveME) { m = mq[ID];} else { m = 0.;}
  mass[1]=0.;
  mass[2]=0.;
  mass[3]=emcm;
  //quark IDs
  IDUP[4] = ID;
  IDUP[5] = -ID;

// Parameters for NLO cross-section calculation.
  double param[8];
  param[0] = pow(m,2)/pow(emcm,2);
  param[1] = 1.+12.*param[0]; // R = 1+param[1]*alphas/pi for vector current
  param[2] = 1.-22.*param[0]; // R = 1+param[2]*alphas/pi for axial vector current
  param[3] = 2.*param[0];
  param[4] = sqrt(1.-4.*param[0]); // beta
  param[5] = 0.5*(1.+param[4]);
  param[6] = (param[0]*param[5]+0.25*param[4]*pow(1+param[4],2))/(param[5]-param[0]);
  param[7] = 2.*alphas/(3.*pi); // colour factor

  // cross-section parameters
  double sigU = param[4]*sigV[ID]+pow(param[4],3)*sigA[ID];
  double sigL = 0.5*(1.-pow(param[4],2))*param[4]*sigV[ID];
  double sigF = pow(param[4],2)*sigVA[ID];
  double CLf = (Aq[ID]+Vq[ID])/2.;
  double CRf = (Vq[ID]-Aq[ID])/2.;
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
   double r2=random(seed);
   double r3=random(seed);
   double r4=random(seed);
    x1 = r1;  
    x2 = r2;
      if (Jt==1) {
	QQgsc(r3, r4, x1, x2, param, wgt, wgtax, wgt1, wgtax1, resolve);
      wt = (wgt*sigV[ID]+wgtax*sigA[ID])/2.;}
     if (Jt==3) {
        QQgh(r3, r4, x1, x2, param, wgt, wgtax, wgt1, wgtax1, resolve);
      wt  = (wgt*sigV[ID]+wgtax*sigA[ID])/2.;}
     x3=2.-x1-x2;
    c13=(x1*x3-2.*(1.-x2))/(sqrt(pow(x1,2)-4.*param[0])*x3);
    c23=(x2*x3-2.*(1.-x1))/(sqrt(pow(x2,2)-4.*param[0])*x3);
    } while (random(seed) > fabs(wt)/fabs(bigwgt) || isinf(c13) || isinf(c23));
    
    if (Jt==1){wt=-1.;}
    if (Jt==3){wt=1.;}    
    double cofac, maxcofac,cf,ct, ppt;
    double px1 = sqrt(pow(0.5*emcm*x1,2)-pow(m,2));
    double px2 = sqrt(pow(0.5*emcm*x2,2)-pow(m,2));
    double px3 = sqrt(pow(0.5*emcm*x3,2));
    do {cth = 2.*random(seed)-1.;
    cofac = (1.+pow(cth,2))*sigU+2.*sigL*(1.-pow(cth,2))+2.*cth*sigF;  
    maxcofac = 2.*(sigU+sigF);
    if (random(seed) < pow(x1,2)/(pow(x1,2)+pow(x2,2))) {
    cf = (px2*c23+px3)/px1; ct = (px1-px3*cf)/px2; }else
    {cf = (px1*c13+px3)/px2; ct = (px2-px3*cf)/px1; qg=false;}
    if (ct > 1. || cf > 1.) {cofac = 0.;}
    } while (random(seed) > cofac/maxcofac);
    double phi = 2.*pi*random(seed);
    double phi2 = 2.*pi*random(seed);
    if (qg) {
      pcm = px1;
      PUP[3][4]=px1*cth;
      double pt = sqrt(pow(px1,2)-pow(PUP[3][4],2));
      PUP[1][4]=pt*cos(phi);
      PUP[2][4]=pt*sin(phi);
      MOTHUP[1][6] = 3;
      MOTHUP[2][6] = 3;
      double ppl = px2*ct;
       if (random(seed) < 0.5) {
	ppt = sqrt (pow(px2,2)-pow(ppl,2));}else {ppt = -sqrt (pow(px2,2)-pow(ppl,2));}
      double ppt1 = ppt*cos(phi2);
      double ppt2 = ppt*sin(phi2);
      PUP[3][5]=-(ppl*cth-ppt1*sin(acos(cth)));
      PUP[1][5]=-((ppl*sin(acos(cth))+ppt1*cth)*cos(phi)-ppt2*sin(phi));
      PUP[2][5]=-((ppl*sin(acos(cth))+ppt1*cth)*sin(phi)+ppt2*cos(phi));} 
    else {
      pcm = px2;
      PUP[3][5]=px2*cth;
      double pt = sqrt(pow(px2,2)-pow(PUP[3][5],2));
      PUP[1][5]=pt*cos(phi);
      PUP[2][5]=pt*sin(phi);
      MOTHUP[1][6] = 3;
      MOTHUP[2][6] = 3;
      double ppl = px1*ct;
      if (random(seed) < 0.5) {
	ppt = sqrt (pow(px1,2)-pow(ppl,2));}else {ppt = -sqrt (pow(px1,2)-pow(ppl,2));}
      double ppt1 = ppt*cos(phi2);
      double ppt2 = ppt*sin(phi2);
      PUP[3][4]=-(ppl*cth-ppt1*sin(acos(cth)));
      PUP[1][4]=-((ppl*sin(acos(cth))+ppt1*cth)*cos(phi)-ppt2*sin(phi));
      PUP[2][4]=-((ppl*sin(acos(cth))+ppt1*cth)*sin(phi)+ppt2*cos(phi));}

    // double ppl2 = px3*cf;
    // double ppt3 = -ppt/fabs(ppt)*sqrt(pow(px3,2)-pow(ppl2,2));
     // double ppt4 = ppt3*cos(phi2);
     //double ppt5 = ppt3*sin(phi2);
     //  PUP[3][6]=-(ppl2*cth-ppt4*sin(acos(cth)));
     // PUP[1][6]=-((ppl2*sin(acos(cth))+ppt4*cth)*cos(phi)-ppt5*sin(phi));
     // PUP[2][6]=-((ppl2*sin(acos(cth))+ppt4*cth)*sin(phi)+ppt5*cos(phi));
      PUP[3][6]= -PUP[3][4] - PUP[3][5];
      PUP[1][6]= -PUP[1][4] - PUP[1][5];
      PUP[2][6]= -PUP[2][4] - PUP[2][5];
      PUP[4][4] = x1*emcm/2.;
      PUP[4][5] = x2*emcm/2.;
      PUP[4][6] = x3*emcm/2.;
      // cout << PUP[1][6] +PUP[1][4]+PUP[1][5] << "\t" << PUP[2][6] +PUP[2][4]+PUP[2][5] << "\t" << PUP[3][6] +PUP[3][4]+PUP[3][5] << endl;
      //colour info
      ICOLUP[1][4] = 501;
      ICOLUP[2][4] = 0;
      ICOLUP[1][5] = 0;
      ICOLUP[2][5] = 502;
      ICOLUP[1][6] = 502;
      ICOLUP[2][6] = 501;
      //Number of particles
      NUP=6;
      
  }
     

     if (Jt==2||Jt==4||Jt==5) {
     x1 = 1.;x2 = 1.; wt = 1.;
     pcm = sqrt(pow(0.5*emcm,2)-pow(m,2));
     double  cofac, maxcofac;
     do {cth = 2.*random(seed)-1.;
       cofac = (1.+pow(cth,2))*sigU+2.*sigL*(1.-pow(cth,2))+2.*cth*sigF;  
       maxcofac = 2.*(sigU+sigF);
     } while (random(seed) > cofac/maxcofac);
     double phi = 2.*random(seed)*pi;
     PUP[3][4] = pcm*cth;
     PUP[1][4] = sqrt(pow(pcm,2)-pow(PUP[3][4],2))*cos(phi);
     PUP[2][4] = sqrt(pow(pcm,2)-pow(PUP[3][4],2))*sin(phi);  
     PUP[1][5] = -PUP[1][4];
     PUP[2][5] = -PUP[2][4];
     PUP[3][5] = -PUP[3][4];
     PUP[4][4] = emcm/2.;
     PUP[4][5] = emcm/2.;
     PUP[3][6]=0.;
     PUP[1][6]=0.;
     PUP[2][6]=0.;
     //colour info
      ICOLUP[1][4] = 501;
      ICOLUP[2][4] = 0;
      ICOLUP[1][5] = 0;
      ICOLUP[2][5] = 501;
    //Number of particles
      NUP=5;
     }
      double alpha = 1.;

     //boost masses if needed. Need to boost to nominal gluon mass always.
      double mb;
      if (massiveME || (!massiveME && boost)) {mb = mq[ID]; mg = mq[6]; }
      if (!massiveME && !boost) {mb = 0.; mg = 0.;}
     double Eq,Eqb,Eg;
     double pq = sqrt(pow(PUP[1][4],2)+pow(PUP[2][4],2)+pow(PUP[3][4],2));
     double pqb = sqrt(pow(PUP[1][5],2)+pow(PUP[2][5],2)+pow(PUP[3][5],2));
     double pg = sqrt(pow(PUP[1][6],2)+pow(PUP[2][6],2)+pow(PUP[3][6],2));
       for (int ix=1; ix < nit+1; ix++) {
       Eq = sqrt(pow(mb,2)+alpha*pow(pq,2));
       Eqb = sqrt(pow(mb,2)+alpha*pow(pqb,2));
       Eg =  sqrt(pow(mg,2)+alpha*pow(pg,2));
       if (NUP==6) {
       alpha = alpha+(2.*(emcm-Eq-Eqb-Eg))/(pow(pq,2)/Eq+pow(pqb,2)/Eqb+pow(pg,2)/Eg);}
       else {alpha = alpha+(2.*(emcm-Eq-Eqb))/(pow(pq,2)/Eq+pow(pqb,2)/Eqb);}
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
     if (massiveME || (!massiveME && boost)) {
     mass[4] = mq[ID];
     mass[5] = mq[ID];
     mass[6] = mq[6];}
     if (!massiveME && !boost)
      {mass[4] = 0.;
       mass[5] = 0.;
       mass[6] = 0.;}
     

      //spin information (not needed yet!)
     double C[9], M[9];
     C[1] = pow(CRe,2)*(chi2/(pow((emcm*k),2)))*pow(gz2,2);
     C[2] = pow(CLe,2)*(chi2/(pow((emcm*k),2)))*pow(gz2,2);
     C[3] = CLf*(emcm/2.-pcm)+CRf*(emcm/2.+pcm);
     C[4] = CLf*(emcm/2.+pcm)+CLf*(emcm/2.-pcm);
     C[5] = (CLf+CRf)*m;
     C[6] = 2*e2*Ch[ID]*CRe*chi1*gz2/(k*emcm);
     C[7] = 2*e2*Ch[ID]*CLe*chi1*gz2/(k*emcm);
     C[8] = e2*Ch[ID];
     M[1] = (C[1]*pow(C[3],2)-C[6]*C[3]+pow(C[8],2))*pow(1+cth,2);
     M[5] = (C[2]*pow(C[4],2)-C[7]*C[4]+pow(C[8],2))*pow(1+cth,2);
     M[2] = (C[1]*pow(C[4],2)-C[6]*C[4]+pow(C[8],2))*pow(1+cth,2);
     M[6] = (C[2]*pow(C[3],2)-C[7]*C[3]+pow(C[8],2))*pow(1+cth,2);
     M[3] = (C[1]*pow(C[5],2)-2.*m/emcm*C[6]*C[5]+pow(C[8],2)*4.*param[0])*(1.-pow(cth,2));
     M[7] = (C[2]*pow(C[5],2)-2.*m/emcm*C[7]*C[5]+pow(C[8],2)*4.*param[0])*(1.-pow(cth,2));
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

     //weight info
     XWGTUP = wt;
     ii+=1;
     // cout << ii << endl;
     double sum1, sum2, sum3, sum4;
     sum1=0.;sum2=0.;sum3=0.;sum4=0.;
     outdata << "<event>" << endl;
     outdata << NUP <<"\t" << "11" <<"\t"<<XWGTUP<<"\t" << SCALUP << "\t" << AQEDUP << "\t" << AQCDUP << endl;
     for (int ja = 1; ja < NUP+1; ja++) {
       outdata << IDUP[ja] << "\t" << ISTUP[ja] << "\t" << MOTHUP[1][ja] <<"\t" << MOTHUP[2][ja] << "\t" << ICOLUP[1][ja] << "\t" << ICOLUP[2][ja] <<"\t" <<setprecision (9)<< PUP[1][ja] << "\t" << PUP[2][ja]<<"\t"<<PUP[3][ja]<<"\t"<<PUP[4][ja]<<"\t" << mass[ja] <<"\t"<< "0" << "\t" <<"9" <<endl;
       sum1+=PUP[4][ja];sum2+=PUP[1][ja];sum3+=PUP[2][ja];sum4+=PUP[3][ja];
       cout  << "Generated event : " << "\t" << ii << "\t" << "of" << "\t" << nevg << "\r" << flush;  
}
         outdata << "</event>" << endl;

	 
	 //   cout << sum1 << "\t" << sum2 << "\t" << sum3 << "\t" << sum4 << endl;
 }

 outdata << "</LesHouchesEvents>" << endl;
 cout << endl;
  return 0;
}
/********************************************************************************/
// In soft and collinear phase space region ?

void QQgsc(double r3,double r4, double &rx, double &ry, double param[8], double &rwgt, double &rwgtax, double &rwgt1, double &rwgtax1, bool resolve) { 
  
  //Initialize weights
  rwgt = 0.;rwgtax = 0.;rwgt1 = 0.;rwgtax1 = 0.;

   /*Remove 2-jet events */
  double y1 = 1.-rx;
  double y2 = 1.-ry;
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
  double tk = param[5];
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
      tk = param[6];
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

     /* Calculate weights for jet regions */
     double zeta = -8.*param[0]*(1.+param[3]);
     double zetax = param[3]*(pow((3.+y1+y2),2)-19.+4.*param[0]);
     double r, zz;
     rx = 1.-y1;
     ry = 1.-y2;
     
    if (region == 1) {
       r = 0.5*(1.+param[0]/(1.+param[0]-ry));
       zz = r+(rx-(2.-ry)*r)/sqrt(pow(ry,2)-2.*param[3]);
      
	 rwgt = ((pow((rx+param[3]),2)+pow((ry+param[3]),2)+zeta)/((1.+param[3])*y1*y2)-param[3]*(1./pow(y1,2)+1./pow(y2,2)))/param[4];
         rwgt = rwgt-(((1.+pow(zz,2))/(1.-zz))-param[3]/y2)/(y2*sqrt(pow(ry,2)-2.*param[3]));
	 
	 rwgt = rwgt*param[7]*fac;

	 rwgtax = ((pow((rx+param[3]),2)+pow((ry+param[3]),2)+zetax)/(pow(param[4],2)*y1*y2)-param[3]*(1./pow(y1,2)+1./pow(y2,2)))/param[4];
         rwgtax = rwgtax-(((1.+pow(zz,2))/(1.-zz))-param[3]/y2)/(y2*sqrt(pow(ry,2)-2.*param[3]));
	 rwgtax = rwgtax*param[7]*fac;

     } else if (region == 2) {
       r = 0.5*(1.+param[0]/(1.+param[0]-rx));
       zz = r+(ry-(2.-rx)*r)/sqrt(pow(rx,2)-2.*param[3]);
       
	 rwgt = ((pow((rx+param[3]),2)+pow((ry+param[3]),2)+zeta)/((1.+param[3])*y1*y2)-param[3]*(1./pow(y1,2)+1./pow(y2,2)))/param[4];
         rwgt = rwgt-(((1.+pow(zz,2))/(1.-zz))-param[3]/y1)/(y1*sqrt(pow(rx,2)-2.*param[3]));
	 rwgt = rwgt*param[7]*fac;
	
	 rwgtax = ((pow((rx+param[3]),2)+pow((ry+param[3]),2)+zetax)/(pow(param[4],2)*y1*y2)-param[3]*(1./pow(y1,2)+1./pow(y2,2)))/param[4];
         rwgtax = rwgtax-(((1.+pow(zz,2))/(1.-zz))-param[3]/y1)/(y1*sqrt(pow(rx,2)-2.*param[3]));
	 rwgtax = rwgtax*param[7]*fac;
	
     }
       rwgt1 = (2.+3.*param[7]*param[1])*fac-rwgt;
       rwgtax1 = (2.+3.*param[7]*param[2])*fac-rwgtax;
   /* Implement soft resolution cut */
     if (resolve) {
           if (rx+ry > 1.9999) {return;}}
}
 /*********************************************************************/
 // In hard region?

void QQgh(double r3,double r4, double &rx, double &ry, double param[8], double &rwgt, double &rwgtax, double &rwgt1, double &rwgtax1, bool resolve) { 
 
  //Initialize weights and phase space points
  rwgt = 0.;rwgtax = 0.;rwgt1 = 0.;rwgtax1 = 0.;

   /*Remove 2-jet events */
  double y1 = 1.-rx;
  double y2 = 1.-ry;
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
}} else if (y2 < y1+2.*pow(y1,2)) {rx = 1.; ry = 1.; return;}
  } else if (y2 < 0.25) { if (y1 < y2+2.*pow(y2,2)) {rx = 1.; ry = 1.; return;} }
  

  if (soft) {if (y1 < y2) {fac = 2.*y1;} else {fac = 2.*y2;} }
  
  if (!soft) {
  
  if (y1 > 0.375 && y2 < 0.25) {
    if (y1 > 1.-2.5*y2 && y1 < 1.-y2) {
      y2 = r4*0.25;
      y1 = 1.-(1.+1.5*r3)*y2;
      fac = 2.*r4;
    } else {rx = 1.; ry = 1.; return;} }
    
   if (y2 > 0.375 && y1 < 0.25) {
	  if (y2 > 1.-2.5*y1 && y2 < 1.-y1) {
	  y1 = r3*0.25;
	  y2 = 1.-(1.+1.5*r4)*y1;
          fac = 2.*r3;
	  } else {rx = 1.; ry = 1.; return;}}}
   
//    Check if outside half-box */
     if (2.-y1-y2 < 1.) {return;}

  

 /*Inside jet 1? */
  double w1 = 1.-y1;
  double w2 = 1.-y2;
  double tk = param[5];
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
      tk = param[6];
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
     double zeta = -8.*param[0]*(1.+param[3]);
     double zetax = param[3]*(pow((3.+y1+y2),2)-19.+4.*param[0]);
          rx = 1.-y1;
          ry = 1.-y2;

         rwgt = ((pow((rx+param[3]),2)+pow((ry+param[3]),2)+zeta)/((1.+param[3])*y1*y2)-param[3]*(1./pow(y1,2)+1./pow(y2,2)))/param[4];
         rwgt = rwgt*param[7]*fac;

	 rwgtax = ((pow((rx+param[3]),2)+pow((ry+param[3]),2)+zetax)/(pow(param[4],2)*y1*y2)-param[3]*(1./pow(y1,2)+1./pow(y2,2)))/param[4];
         rwgtax = rwgtax*param[7]*fac;

	 rwgt1 = (2.+3.*param[7]*param[1])*fac-rwgt;
        rwgtax1 = (2.+3.*param[7]*param[2])*fac-rwgtax;

	/* Implement soft resolution cut */
     if (resolve) {
           if (rx+ry > 1.9999) {return;}}
}
 /*********************************************************************/
 // In gap between phase space for massive partons and full half-triangle? 

void QQgsp(double r3,double r4, double &rx, double &ry, double param[8], double &rwgt, double &rwgtax, double &rwgt1, double &rwgtax1, bool resolve) {

  //Initialize weights and phase space points
  rwgt = 0.;rwgtax = 0.;rwgt1 = 0.;rwgtax1 = 0.;

   /*Remove 2-jet events */
  double y1 = 1.-rx;
  double y2 = 1.-ry;
  if ((1.-y1)==1. && (1.-y2)==1.) {return;}
  /* Check if outside half-box */
  if (2.-y1-y2 <= 1.) {return;}

/* Inside Dalitz plot? */
 if (y1*y2*(1.-y1-y2) > param[0]*pow((y1+y2),2)) {return;}
 /* Calculate weights */
 rwgt1 = 2.+3.*param[7]*param[1];
 rwgtax1 = 2.+3.*param[7]*param[2];}
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
 /*********************************************************************/
 /*********************************************************************/  
