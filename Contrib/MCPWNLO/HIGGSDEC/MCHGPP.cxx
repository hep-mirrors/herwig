#include <cstdlib> 
#include <ctime> 
#include <iostream>
#include <fstream>
#include <string>
//#include "LHAPDFWrap.h"
//#include "/usera/seyi/lhapdf/include/LHAPDF/LHAPDF.h"
#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "MCHGPP_INPUTS.h"
#include <cmath>
#include <iomanip>

using namespace LHAPDF;
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
  double emcm = user.cme();
  double mH = user.mH();
  int seed = user.rseed();
  int proc = user.VBID(); /*process ID: Z production Z/gamma, W */
  bool pp = user.pp(); /*type of accelerator */
  bool resolve = user.resolve();
  char* pdfset = user.PDFset(); /* PDF set */
  int ii=0; //event count
  // Constants 
  double pi = 3.14152654;
  double GF = 0.0000116639;
  double m  = user.mb();
  double Mw = user.Mw();  
  double Mz = user.Mz();
  double Yz = 2.486;
  double alphaem = 0.007297352;
  double alphas = user.alphasmH();
  //Lepton charge & coupling constants
  double Che = -1;
  double sin2thw = 0.2312;
  double cos2thw = 1.-0.2312;
  double Ae = -0.5;
  double Ve = Ae-2.*Che*sin2thw;
  double CLe = (Ae+Ve)/2.;
  double CRe = (Ve-Ae)/2.;
  double Mv = Mw;
  if (proc == 23) {Mv = Mz;}
  double gz2 = sqrt(32.)*GF*pow(Mv,2);
  double e2 = 4.*pi*alphaem;
  double k = sqrt(2.)*GF*pow(Mv,2)/e2;
  double chi1 = k*pow(mH,2)*(pow(mH,2)-pow(Mv,2))/(pow((pow(mH,2)-pow(Mv,2)),2)+pow(Yz*Mv,2));
  double chi2 = pow(k,2)*pow(Mv,4)/(pow((pow(Mv,2)-pow(Mv,2)),2)+pow(Yz*Mv,2));
  //quark masses ,velocities and couplings. NB:gluon mass[10] =0.75 
  double mg, Ch[6];
    mg = 0.75; // nominal gluon mass
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
 	 cout << "PHASE SPACE INTEGRALS:" << "\t" << p[ic] << "/" << p1[ic] <<"/" << sd[ic] <<"/" << sd1[ic]
	      <<"/" << bigwt[ic] << "\t" << xsec << endl;
}


  //*********************************************************************************
  /*Les Houches Accord event parameters:: PUP are the momenta, 
   ICOLUP is the colour code, MOTHUP is the mother information, ISTUP is 
   the status code. ISPINUP signifies spins. IDUP is the particle ID.
    NUP is the number of particles in the event. SCALUP is the scale of 
   the event.AQEDUP and AQCDUP are alphaem and alphas resp.XWGTUP is 
   the weight of the event=1 or -1 for unweighted events.*/
  
  double PUP[6][11],AQEDUP,AQCDUP,XWGTUP,SCALUP;
  int ICOLUP[3][11],MOTHUP[3][11],ISTUP[11],ISPINUP[11],IDUP[11],NUP;
  NUP=0;
  SCALUP = mH;
  AQEDUP=0.007297352;
  AQCDUP=alphas;
  IDUP[4]=25;
  IDUP[6]=5;
  IDUP[7]=-5;
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
  //Allocate phase space region
  double cst[3],bigwgt;
  bigwgt=0.;
  cst[0]=fabs(p[0])+fabs(p1[0]);
  cst[1]=fabs(p[1])+fabs(p1[1]);
  cst[2]=fabs(p[2])+fabs(p1[2]);
  cstt = cst[0]+cst[1]+cst[2];

  /*************************************************************************/

  /* pdf routine. Calls external LHAPDF which calculates the
     parton distribution functions*/
  const int SUBSET = 0;
  //    const string NAME = "cteq4d.LHgrid";
  setVerbosity(SILENT); // or LOWKEY

  initPDFByName(pdfset, SUBSET);
  //     int subset = 0;
  //   char name[160];
  //     string byname = user.PDFset();
  // //     if (DIS) {
  // //      byname = "/usera/seyi/lhapdf-5.3.0/PDFsets/cteq5d.LHgrid";} else
  // //     {byname = "/usera/seyi/lhapdf-5.3.0/PDFsets/cteq5m.LHgrid";}
  //     strcpy(name,byname.c_str());
  //     LHAPDFWrap pdf = LHAPDFWrap();
  //     pdf.initPDFSet(name);
  //     pdf.initPDF(subset);

  /*************************************************************************/
  //cross-section estimator
  /* Type of accelerator */
  int sgn=1;
  if (pp) {sgn=-1;}
  double t0 = pow(Mv + mH,2)/pow(emcm,2);
  double tcst = 0.;
  double fwmax = 0.;
  for (int ix =0; ix < 10000 ; ix++) {
    t0 = pow((Mv+mH)/emcm,2);
    double xq1 = random(seed)*(1.-t0)+t0;
    double xq2 = random(seed)*(1.-t0/xq1) + t0/xq1;
   double FW = 1./3.*(xfx(xq1,Mv+mH,2)*xfx(xq2,Mv+mH,sgn*1)+xfx(xq2,Mv+mH,-sgn*2)*xfx(xq1,Mv+mH,-1)+xfx(xq1,Mv+mH,-2)*xfx(xq2,Mv+mH,-sgn*1)+xfx(xq2,Mv+mH,sgn*2)*xfx(xq1,Mv+mH,1)); 
   if (proc ==23){2./3.*((xfx(xq1,Mv+mH,1)*xfx(xq2,Mv+mH,sgn*1)+xfx(xq2,Mv+mH,-sgn*1)*xfx(xq1,Mv+mH,-1))*(1./4.-2./3.*sin2thw+8./9.*pow(sin2thw,2))+(xfx(xq1,Mv+mH,2)*xfx(xq2,Mv+mH,sgn*2)+xfx(xq2,Mv+mH,-sgn*2)*xfx(xq1,Mv+mH,-2))*(1./4.-1./3.*sin2thw+4./9.*pow(sin2thw,2)));}
  double s2s = xq1*xq2*pow(emcm,2);
  double betav = sqrt((s2s-pow(Mv+mH,2))*(s2s-pow(Mv-mH,2)))/(s2s-pow(mH,2)+pow(Mv,2));
  double gammav = 1./sqrt(1.-pow(betav,2));
  //  double fw = 1./s2s*sqrt(1.-pow(Mv+mH,2)/s2s)*sqrt(1.-pow(Mv-mH,2)/s2s)*((1.-(Mv+mH,2)/s2s)*(1.-(Mv-mH,2)/s2s)+12.*pow(Mv,2)/s2s)*(1./pow(1.-Mv/s2s,2));
  double fw = 1./sqrt(pow(s2s,3))*betav*gammav*pow((s2s+pow(Mv,2) - pow(mH,2))/(s2s-pow(Mv,2)),2)*(1.-2.*pow(betav,2)/3.); 
  double fww =1./sqrt(pow(s2s,3))*betav*gammav*pow((s2s+pow(Mv,2) - pow(mH,2))/(s2s-pow(Mv,2)),2)*(2.*(1-pow(betav,2))+pow(betav,2));
  tcst += pow(GF*Mw,2)*pow(Mv,3)/(4.*pi*cos2thw)*FW*fw;
  if (fww > fwmax) {fwmax = fww;}
 }
  tcst = tcst/10000.*fabs(1.-1./t0)*389400000.;

  /*************************************************************************/
   ofstream outdata;
    outdata.open("PPVH.dat", ios::trunc);//output for events
   outdata << "<LesHouchesEvents version =\"1.0\">" << endl; outdata << "<!--" << endl;
   outdata << "IPROC = " <<"\t" << proc << endl;
   outdata << "NORMFACTOR = " << "\t" << xsecabs/xsec << endl;
   outdata << "MC@NLO W+Higgs boson production" <<endl ;
     outdata << "File generated with MCHGPP.cxx" << endl; outdata << "-->" << endl;
     outdata << "<init>" << endl;
     outdata << "\t2212\t" << "-2212\t" ;
     outdata << emcm/2. << "\t" << emcm/2.  << "\t" << "0 \t 0 \t 7\t 7 \t 3 \t 1" << endl;
     outdata << "\t" << tcst   << "\t" << 0.000000 << "\t1.00000 \t11" << endl;     outdata << "</init>" << endl;


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

  /*************************************************************************/
  int ID = 5;
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




  //Momenta of incoming particles and W, W* boson
  double ctheta, yy, Mb, xq, xqb, MVH, Ev, Eh, pv, S, S2;
  S = emcm;
  S2 = pow(S,2);
  double Qmin = Mv + mH;
  double xqc=sqrt(pow(Qmin,2)/S2);
  double fvmax=0.;

  /*Define the coupling constants. */
  /*Define the coupling constants. */
  double v[6],a[6],va[6];
  double sin2thw=0.2312;
  v[1]=-0.5+(2./3.)*sin2thw;
  v[2]=0.5-(4./3.)*sin2thw;
  v[3]=v[1];v[4]=v[2];v[5]=v[1];a[1]=0.5;a[2]=-0.5;a[3]=a[1];a[4]=a[2];a[5]=a[1];
  for (int ij=1; ij <6;ij++) {
    va[ij]=(pow(v[ij],2)+pow(a[ij],2));}
  double vq[7] ;
  vq[1]=0.952; vq[2]=0.049; vq[3]=0.952;
  vq[4]=0.00168; vq[5]=0.049; vq[6]=0.00001063;

  /*************************************************************************/


  /*PDF terms. [fq1,fq2] and [fqb1,fqb2] are pdfs for quarks from hadrons [1,2]
    and antiquarks from hadrons [1,2] resp. ffq1 and ffq2 are the
    2 combinations of the products of the Born process pdfs.*/

  double fv, fv1,fv2,fq1[7],fq2[7],fqb1[7],fqb2[7], Q;
  fv1 =0.; fv2 =0.; 
  if (proc == 24) {
  for(int id=1; id<5; id++) {
    fvmax+=2.*(xfx(xqc,Qmin,id)*xfx(xqc,Qmin,sgn*(id+1))+xfx(xqc,Qmin,-id)*xfx(xqc,Qmin,-sgn*(id+1)))*vq[id]/(xqc*xqc);}
  for(int id=1; id<3; id++) {
    fvmax+=2.*(xfx(xqc,Qmin,id)*xfx(xqc,Qmin,sgn*(id+3))+xfx(xqc,Qmin,-id)*xfx(xqc,Qmin,-sgn*(id+3)))*vq[id+4]/(xqc*xqc);}}
  if (proc==23) {
    for(int id=1; id<6; id++) {
      fvmax+=(xfx(xqc,Qmin,id)*xfx(xqc,Qmin,sgn*id)+xfx(xqc,Qmin,-id)*xfx(xqc,Qmin,-sgn*id))*(pow(v[id],2)+pow(a[id],2))/(xqc*xqc);}}
  //    Mbmax = (pow(-pow(mH,2)-2.*pow(Mv,2),2)+pow(Mv,2)*(pow(mH,2)+pow(Mv,2)))*sqrt(pow(pow(mH,2)+pow(Mv,2),2) + pow(mH,4)+pow(Mv,4))/(pow(mH,4)*pow(Mv,2)*pow(pow(mH,2)+pow(Mv,2),2));
  do {
        
    //  xx = random(seed);
    //  yy = 2.*random(seed)-1.;
    // ss = (pow(mH,2)+pow(Mv,2))/xx;
    // ll = ss*(ss+pow(Mv,4)/ss+pow(mH,4)/ss-(2.*pow(Mv,2)*pow(mH,2)/ss+2.*pow(mH,2)+2.*pow(Mv,2)));
    // tt = -ss/2.*(1.-xx - sqrt(ll)*yy/ss);
    //  uu = -ss/2.*(1.-xx + sqrt(ll)*yy/ss);
    //  Mb = (sqrt(ll)*(tt-pow(Mv,2))*(uu-pow(Mv,2))+pow(Mv,2)*ss)/(pow(ss-pow(Mv,2),2)*pow(Mv,2)*pow(ss,2));
    //  Ew = (ss+pow(Mv,2)-pow(mH,2))/(2.*sqrt(ss));
    // pw = sqrt(ll)/(2.*sqrt(ss));
    //  Eh = sqrt(pow(pw,2)+pow(mH,2));
    // MWH = Ew+Eh;
    // Q = MWH;
    // if (ll > 0. && MWH < S) {
    //  do {
    //  yz = 6.*random(seed) - 3.0 ; //rapidity range -3:3
     
    //xq=sqrt(pow(MWH,2)/S2)*exp(yz);
    //  xqb=sqrt(pow(MWH,2)/S2)*exp(-yz);
    //  } while (xq > 1. || xqb > 1.);}
    double tau = pow((Mv+mH)/emcm, 2);
    ctheta = 2.*random(seed)-1.;
       do {
    xq = random(seed)*(1.-tau) + tau;
    xqb = random(seed)*(1.-tau/xq) +tau/xq;
     } while(fabs(0.5*log(xq/xqb)) > 3.); 
    yy = xq*xqb*pow(emcm,2);
    double betav = sqrt((yy-pow(Mv+mH,2))*(yy-pow(Mv-mH,2)))/(yy-pow(mH,2)+pow(Mv,2)) ;
    double gammav = 1./sqrt(1.-pow(betav,2));
    Mb =  (betav*gammav)/sqrt(pow(yy,3))*pow((yy+pow(Mv,2)-pow(mH,2))/(yy-pow(Mv,2)),2)*(2.*(1.-pow(betav,2))+pow(betav,2)*(1.-pow(ctheta,2)));
    Ev = (yy+pow(Mv,2)-pow(mH,2))/(2.*sqrt(yy));
    pv = sqrt(pow(Ev,2) - pow(Mv,2));
    Eh = sqrt(pow(pv,2)+pow(mH,2));
    MVH = Ev+Eh;
    Q = MVH;
    if (proc==23) {
      for (int e=1; e<6; e++) {
	fq1[e] = xfx(xq,Q,e);
	fqb2[e]= xfx(xqb,Q,sgn*e);
	fqb1[e] =xfx(xq,Q,-e);
	fq2[e] = xfx(xqb,Q,-sgn*e);
	fv1 += va[e]*fq1[e]*fqb2[e];
	fv2 += va[e]*fqb1[e]*fq2[e];
      }} else if(proc==24){
     for(int e=1; e<5; e++) {
      fq1[e] = xfx(xq,Q,e)*xfx(xqb,Q,sgn*(e+1))+xfx(xq,Q,-e)*xfx(xqb,Q,-sgn*(e+1))+xfx(xqb,Q,sgn*e)*xfx(xq,Q,e+1)+xfx(xqb,Q,-sgn*e)*xfx(xq,Q,-(e+1));}
     int ee=1;
     for(int e=5; e<7; e++) {fq1[e] = xfx(xq,Q,ee)*xfx(xqb,Q,sgn*(ee+3))+xfx(xq,Q,-ee)*xfx(xqb,Q,-sgn*(ee+3))+xfx(xqb,Q,sgn*ee)*xfx(xq,Q,ee+3)+xfx(xqb,Q,-sgn*ee)*xfx(xq,Q,-(ee+3));
      ee+=1;}
     for (int e=1; e<7; e++) {
      fv1 += vq[e]*fq1[e];
      fv2 = 0;}}
     fv = (fv1 + fv2)/(xq*xqb);
      } while (random(seed) > fv*Mb/(fvmax*fwmax));

 // Distributing the Vector boson and H momenta
  
  PUP[4][1]=xq*S/2.;
  PUP[1][1]=0.;
  PUP[2][1]=0.;
  PUP[3][1]=xq*S/2.;
  PUP[4][2]=xqb*S/2.;
  PUP[1][2]=0.;
  PUP[2][2]=0.;
  PUP[3][2]=-xqb*S/2.;
  PUP[4][3]=PUP[4][1]+PUP[4][2];
  PUP[1][3]=0.;
  PUP[2][3]=0.;
  PUP[3][3]=PUP[3][1]+PUP[3][2];
  // Still in parton C.O.M
  double phi3 = 2.*random(seed)*pi;
  PUP[4][5]=Ev;
  PUP[1][5]=pv*sin(acos(ctheta))*cos(phi3);
  PUP[2][5]=pv*sin(acos(ctheta))*sin(phi3);
  PUP[3][5]=pv*ctheta;
  PUP[4][4]=Eh;
  PUP[1][4]=-PUP[1][5];
  PUP[2][4]=-PUP[2][5];
  PUP[3][4]=-PUP[3][5];
  // Boost to lab frame
  double vv = PUP[3][3]/PUP[4][3];
  PUP[4][5] = (Ev + vv*PUP[3][5])/sqrt(1.-pow(vv,2));
  PUP[3][5] = (PUP[3][5] + vv*Ev)/sqrt(1.-pow(vv,2));
  PUP[4][4] = (Eh + vv*PUP[3][4])/sqrt(1.-pow(vv,2));
  PUP[3][4] = (PUP[3][4] + vv*Eh)/sqrt(1.-pow(vv,2));
  

 /*************************************************************************/


  /*PDF terms. [fq1,fq2] and [fqb1,fqb2] are pdfs for quarks from hadrons [1,2]
    and antiquarks from hadrons [1,2] resp. ffq1 and ffq2 are the
    2 combinations of the products of the Born process pdfs.*/

  fv1 = 0;
  double ffq1[7], ffq2[7];
  if (proc==23) {
    for (int e=1; e<6; e++) {
      fq1[e] = xfx(xq,Q,e);
      fqb2[e]= xfx(xqb,Q,sgn*e);
      fqb1[e] =xfx(xq,Q,-e);
      fq2[e] = xfx(xqb,Q,-sgn*e);
      ffq1[e] = fq1[e]*fqb2[e];
      ffq2[e] = fqb1[e]*fq2[e];

      fv1 += va[e]*ffq1[e];
      fv2 += va[e]*ffq2[e];
    }} else if (proc==24) {
  for(int e=1; e<5; e++) {
    fq1[e] = xfx(xq,Q,e)*xfx(xqb,Q,sgn*(e+1))+xfx(xq,Q,-e)*xfx(xqb,Q,-sgn*(e+1))+xfx(xqb,Q,sgn*e)*xfx(xq,Q,e+1)+xfx(xqb,Q,-sgn*e)*xfx(xq,Q,-(e+1));
  }
  int ee=1;
  for(int e=5; e<7; e++) {
    fq1[e] = xfx(xq,Q,ee)*xfx(xqb,Q,sgn*(ee+3))+xfx(xq,Q,-ee)*xfx(xqb,Q,-sgn*(ee+3))+xfx(xqb,Q,sgn*ee)*xfx(xq,Q,(ee+3))+xfx(xqb,Q,-sgn*ee)*xfx(xq,Q,-(ee+3));
    ee+=1;}
  for (int e=1; e<7; e++) {
    fv1 += vq[e]*fq1[e];
    fv2 = 0;}}

  /* Select incoming parton IDs */
  int ID1, ID2, IDD[3];
  ID1 = 0; ID2 =0;
  double g1= 0 ;
  double g2 = 0;
  double g = 0;
  if (proc==23) {

    if (fv1 !=0 ) {
      if (random(seed)< va[1]*ffq1[1]/fv1) {g1 = va[1]*ffq1[1];
      ID1 = 1;}
      else if (random(seed) < va[2]*ffq1[2]/(fv1-va[1]*ffq1[1])) {g1 = va[2]*ffq1[2];
      ID1 = 2;}
      else if (random(seed) < va[3]*ffq1[3]/(fv1-va[1]*ffq1[1]-va[2]*ffq1[2])) {g1 = va[3]*ffq1[3];
      ID1 = 3;}
      else if (random(seed) < va[4]*ffq1[4]/(va[4]*ffq1[4]+va[5]*ffq1[5])) {g1 = va[4]*ffq1[4];
      ID1 = 4;}
      else {g1 = va[5]*ffq1[5];
      ID1 = 5;}}
    if (fv2 !=0 ) {
      if (random(seed) < va[1]*ffq2[1]/fv2) {g2 = va[1]*ffq2[1];
      ID2 = -1;}
      else if (random(seed) < va[2]*ffq2[2]/(fv2-va[1]*ffq2[1])) {g2 = va[2]*ffq2[2];
      ID2 = -2;}
      else if (random(seed) < va[3]*ffq2[3]/(fv2-va[1]*ffq2[1]-va[2]*ffq2[2])) {g2 = va[3]*ffq2[3];
      ID2 = -3;}
      else if (random(seed) < va[4]*ffq2[4]/(va[4]*ffq2[4]+va[5]*ffq2[5])) {g2 = va[4]*ffq2[4];
      ID2 = -4;}
      else {g2 = va[5]*ffq2[5];
      ID2 = -5;}}
    if (random(seed)< g1/(g1+g2)) {g = g1;
    IDD[0] = ID1; IDD[1] = -ID1;}
    else {g = g2;
    IDD[0] = ID2; IDD[1]= -ID2; }} else if(proc == 24) {
 
  if (random(seed) < vq[1]*fq1[1]/fv1) {g1=vq[1]*fq1[1];
  ID1=1;}
  else if (random(seed) < vq[2]*fq1[2]/(fv1-vq[1]*fq1[1])) {g1 = vq[2]*fq1[2];
  ID1 = 2;}
  else if (random(seed) < vq[3]*fq1[3]/(fv1-vq[1]*fq1[1]-vq[2]*fq1[2])) {g1 = vq[3]*fq1[3];
  ID1 = 3;}
  else if (random(seed) < vq[4]*fq1[4]/(vq[4]*fq1[4]+vq[5]*fq1[5]+vq[6]+fq1[6])) {g1 = vq[4]*fq1[4];
  ID1 = 4;}
  else if (random(seed) < vq[5]*fq1[5]/(vq[5]*fq1[5]+vq[6]*fq1[6])) {g1 = vq[5]*fq1[5];
  ID1 = 5;}
  else {g1 = vq[6]*fq1[6];
  ID1 = 6;}
  if(ID1 < 5) {
    if (random(seed) < xfx(xq,Q,ID1)*xfx(xqb,Q,sgn*(ID1+1))/fq1[ID1]) {
      IDD[0]=ID1;IDD[1]=-(ID1+1);}
   else if (random(seed) < xfx(xq,Q,-ID1)*xfx(xqb,Q,-sgn*(ID1+1))/(fq1[ID1]-xfx(xq,Q,ID1)*xfx(xqb,Q,sgn*(ID1+1)))) {IDD[0]=-ID1;IDD[1]=ID1+1;}
    else if (random(seed) < xfx(xq,Q,ID1+1)*xfx(xqb,Q,sgn*ID1)/(xfx(xq,Q,ID1+1)*xfx(xqb,Q,sgn*ID1)+xfx(xq,Q,-(ID1+1))*xfx(xqb,Q,-sgn*ID1))) {IDD[0]=ID1+1;IDD[1]=-ID1;}
    else {IDD[0]=-(ID1+1);IDD[1]=ID1;} }

  if (ID1 > 4) {
    if (random(seed) < xfx(xq,Q,ID1-1)*xfx(xqb,Q,sgn*(ID1-4))/fq1[ID1]) {IDD[0]=ID1-1;IDD[1]=-(ID1-4);}
    else if (random(seed) < xfx(xq,Q,-(ID1-1))*xfx(xqb,Q,-sgn*(ID1-4))/(fq1[ID1]-xfx(xq,Q,ID1-1)*xfx(xqb,Q,sgn*(ID1-4)))) {IDD[0]=-(ID1-1);IDD[1]=ID1-4;}
    else if (random(seed) < xfx(xq,Q,ID1-4)*xfx(xqb,Q,sgn*(ID1-1))/(xfx(xq,Q,ID1-4)*xfx(xqb,Q,sgn*(ID1-1))+xfx(xq,Q,-(ID1-4))*xfx(xqb,Q,-sgn*(ID1-1))))
      {IDD[0]=ID1-4;IDD[1]=-(ID1-1);}
    else {IDD[0]=-(ID1-4);IDD[1]=ID1-1;} }}

  //set flavour and masses of process
  double mass[10];
  mass[1]=0.;
  mass[2]=0.;
  mass[3]=MVH;
  mass[4]=mH;
  mass[5]=Mv;
  mass[8]=0.;
  mass[9]=0.;
  //quark IDs
  IDUP[6] = ID;
  IDUP[7] = -ID;
  IDUP[1] = IDD[0];
  IDUP[2] = IDD[1];
  ICOLUP[1][2] =503;
  ICOLUP[2][2] =0;
  ICOLUP[1][1] =0;
  ICOLUP[2][1] =503;
  if (IDUP[1] > 0) {
    ICOLUP[1][1] =503;
    ICOLUP[2][1] =0;
    ICOLUP[1][2] =0;
    ICOLUP[2][2] =503;
}

  double cth2, cofac, maxcofac;
  if (proc == 24) {
    do {cth2=2.*random(seed)-1.;
    cofac=pow(1.+cth2,2);
    maxcofac=4.;} while(random(seed) > cofac/maxcofac);

    IDUP[3] = -24;
    IDUP[5] = -24;
    IDUP[8] = 11;
    IDUP[9] = -12;
   if (IDD[0]==-1||IDD[0]==-3||IDD[0]==-5 ||IDD[0]==2||IDD[0]==4){
       IDUP[3]=24;
       IDUP[5] =24;
       IDUP[8] = -11;
       IDUP[9] = 12;

     cth2=-cth2;
}} else if (proc==23) {
     IDUP[3]=23;
     IDUP[5]=23;
     IDUP[8]=11;
     IDUP[9]=-11;
     do
       {cth2=2.*random(seed)-1.;
       cofac = 1.+pow(cth2,2);
       maxcofac = 2.;}while(random(seed) > cofac/maxcofac);
   } 
 
  double phi2=2.*pi*random(seed);
  //Assign 4-momenta in boson rest frame
  PUP[3][9] = Mv*cth2/2.;
  PUP[1][9] = sqrt(pow(Mv/2.,2)-pow(PUP[3][9],2))*cos(phi2);
  PUP[2][9] = sqrt(pow(Mv/2.,2)-pow(PUP[3][9],2))*sin(phi2);
  PUP[1][8] = -PUP[1][9];
  PUP[2][8] = -PUP[2][9];
  PUP[3][8] = -PUP[3][9];
  PUP[4][8] = Mv/2.;
  PUP[4][9] = Mv/2.;

  double vv2 = sqrt(pow(PUP[4][5],2)-pow(Mv,2))/PUP[4][5];
  PUP[4][9]=(Mv/2.-vv2*PUP[3][9])/sqrt(1.-pow(vv2,2));
  PUP[3][9]=(PUP[3][9]-Mv*vv2/2.)/sqrt(1.-pow(vv2,2));
  PUP[4][8]=(Mv/2.-vv2*PUP[3][8])/sqrt(1.-pow(vv2,2));
  PUP[3][8]=(PUP[3][8]-Mv*vv2/2.)/sqrt(1.-pow(vv2,2));
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
  PUP[1][9]=-(pz2*sin(tht)-ktt2*cos(tht))*cos(phi3);
  PUP[2][9]=-(pz2*sin(tht)-ktt2*cos(tht))*sin(phi3);
  PUP[1][8]=-(pz1*sin(tht)-ktt1*cos(tht))*cos(phi3);
  PUP[2][8]=-(pz1*sin(tht)-ktt1*cos(tht))*sin(phi3);
  if (random(seed) < 0.5) {
   PUP[4][9]=pe1;
   PUP[4][8]=pe2;
   PUP[3][9]=-(pz1*cos(tht)+ktt1*sin(tht));
   PUP[3][8]=-(pz2*cos(tht)+ktt2*sin(tht));
   PUP[1][9]=-(pz1*sin(tht)-ktt1*cos(tht))*cos(phi3);
   PUP[2][9]=-(pz1*sin(tht)-ktt1*cos(tht))*sin(phi3);
   PUP[1][8]=-(pz2*sin(tht)-ktt2*cos(tht))*cos(phi3);
   PUP[2][8]=-(pz2*sin(tht)-ktt2*cos(tht))*sin(phi3);
 }


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
   double r2=random(seed); 
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
    cf = (px2*c23+px3)/px1; ct = (px1-px3*cf)/px2; }else
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
      
  }

     if (Jt==2||Jt==4||Jt==5) {
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
     PUP[4][10]=0.;
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
       if (NUP == 9) {mg = 0.;} else {mg =0.75;}
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
      // boost to lab frame
     double Eqq = PUP[4][6] ; double Eqqb = PUP[4][7] ; double Eggg= PUP[4][10]; 
     vv2 = sqrt(pow(PUP[4][4],2)-pow(mH,2))/PUP[4][4];
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
     costx = -cos(phi3);
     sintx = -sin(phi3);
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
     mass[6] = m;
     mass[7] = m;
     mass[10] = mg;
     double CLf, CRf;
     CLf =0.; CRf = 0.;
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
      outdata << "</event>" << endl;
      //   cout << sum1 << "\t" << sum2 << "\t" << sum3 << "\t" << sum4 << endl;
 }

 outdata << "</LesHouchesEvents>" << endl;
 cout << endl;
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
  //  fac = 1.;
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
