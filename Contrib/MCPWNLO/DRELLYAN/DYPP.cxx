/* This program calculates for Drell Yan lepton pair production the integrals 
defined in equations (7.3) and (7.4) on page 28 of arXiv:07084390v1 and 
Appendices F and G. It also calculates the maximum values of the integrand 
in different regions of phase space which will be used for event generation 
in a sister program regZZ.cxx. */

#include <cstdlib> 
#include <ctime> 
#include <iostream>
#include <string>  
#include "LHAPDF/LHAPDF.h"
#include <fstream>
#include "DYPP_INPUTS.h"
#include <cmath>
#include <iomanip>

using namespace LHAPDF;
using namespace std;

//POWHEG functions
double G(double &rp2, double &rNq, double &rbb, double &rqq, double &rpm, double &rll);
double DG(double &rp2, double &rNq,double &rbb,double &rqq,double &rpm,double &rll);

//MC@NLO functions
double DYBC(double &rx,double &ry, double r3, double r4);
double DYD(double &rx, double &ry, double r3, double r4);
double DYGBC(double &rx, double &ry, double r3, double r4);
double DYGD(double &rx, double &ry, double r3, double r4);

//Random number generator
double random (int &rseed);

 /*********************************************************************/
 /*********************************************************************/
 /*********************************************************************/

// Main program for cross-section calculation.

  int main() {
  Input user;
  int seed = user.rseed(); /* initial seed for random number generator */
  char* scheme = user.rscheme(); /* DIS scheme or MSbar scheme */
  char* pdfset = user.PDFset(); /* PDF set */
  int proc = user.VBID(); /*process ID: Z production Z/gamma, W */
  bool pp = user.acc(); /* pp or ppbar */

  /* Define constants: emcm=center of mass energy, CF and TF are 
    colour factors respectively for gluon and quark emission. */

  double emcm = user.cme(); 
  double S2 = emcm*emcm;
  double CF = 4./3.;
  double TF =0.5;
  double pi = 3.14152654;

   bool MCNLO = user.MCNLO();
   bool truncsh = user.trunc() ; /* Implement truncated shower? */
   int nev = user.nevint(); /*number of events integrated over */
   int nevg = user.nevgen(); /*number of events to generate */
   int ii=0; //event count
 
    /*Define the coupling constants. */
     double v[6],a[6],va[6],vq[6];
    double sin2thw=0.2312; 
    v[1]=-0.5+(2./3.)*sin2thw;
    v[2]=0.5-(4./3.)*sin2thw;
    v[3]=v[1];v[4]=v[2];v[5]=v[1];a[1]=0.5;a[2]=-0.5;a[3]=a[1];a[4]=a[2];a[5]=a[1];
    for (int ij=1; ij <6;ij++) {
      va[ij]=(pow(v[ij],2)+pow(a[ij],2));}

    vq[1]=0.952; vq[2]=0.049; vq[3]=0.952;
    vq[4]=0.00168; vq[5]=0.049; vq[6]=0.00001063;
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
   
    
    /*Define variables for boson */
    double mp, Yzz;
    if (proc==23) {
    mp = user.Mz();   /* Z boson pole mass*/
    Yzz = user.widthz(); /* Z boson width */ }
    else if (proc==24) {
    mp = user.Mw();   /* W boson pole mass*/
    Yzz =user.widthw(); /* W boson width */ }
    double Qmax=mp+user.halfwidth()*Yzz;
    double Qmin=mp-user.halfwidth()*Yzz;
    double rhomin=atan((pow(Qmin,2)-pow(mp,2))/(Yzz*mp));
    double rhomax=atan((pow(Qmax,2)-pow(mp,2))/(Yzz*mp));
    double xqc=sqrt(pow(Qmin,2)/S2);
    double fvmax=0.;
    int sgn=1;
    if (pp) {sgn=-1;}
 /* Calculate maximum value for the  Born cross-section fvmax */
    int id;
     if (proc==23) {
       for(id=1; id<6; id++) {  
	 fvmax+=(xfx(xqc,Qmin,id)*xfx(xqc,Qmin,sgn*id)+xfx(xqc,Qmin,-id)*xfx(xqc,Qmin,-sgn*id))*(pow(v[id],2)+pow(a[id],2))/(xqc*xqc);}}
     if (proc==24) {
       for(id=1; id<5; id++) {
	 fvmax+=2.*(xfx(xqc,Qmin,id)*xfx(xqc,Qmin,sgn*(id+1))+xfx(xqc,Qmin,-id)*xfx(xqc,Qmin,-sgn*(id+1)))*vq[id]/(xqc*xqc);}
       for(id=1; id<3; id++) {
	 fvmax+=2.*(xfx(xqc,Qmin,id)*xfx(xqc,Qmin,sgn*(id+3))+xfx(xqc,Qmin,-id)*xfx(xqc,Qmin,-sgn*(id+3)))*vq[id+4]/(xqc*xqc);}
}

/* Generate Born variables Q, yz and hence xq and xqb which are Born process mmt fractions. Q = generated boson mass, yz= boson rapidity generated for the Born process*/
 
    double xq,xqb,x,y,Q,yz;

 /*PDF terms. [fq1,fq2] and [fqb1,fqb2] are pdfs for quarks from hadrons [1,2]
    and antiquarks from hadrons [1,2] resp. ffq1 and ffq2 are the 
    2 combinations of the products of the Born process pdfs.*/

    double ffq1[7],ffq2[7],fv, fv1,fv2,fq1[7],fq2[7],fqb1[7],fqb2[7];
    
  /* w1 and w2 are the integrands. p1 and p2 are the corresponding integrals 
     calcaulated from sum1 and sum2 & sd1 and sd2 are the corresponding errors 
     calculated via ssq1 and ssq2.  bigwt is the maximum value of w1 over 
     the phase space (will be used for event generation) */

   double w1[4],w2[4],p1[4],p2[4],p3[4],sd1[4],sd2[4],sd3[4];
   double bigwt[4],sum1[4],sum2[4],sum3[4],ssq1[4],ssq2[4],ssq3[4];
   double sumwgt,jgmax;
   sumwgt = 0.;
   //POWHEG stuff
  double rmwgt;
  double bigrmwgt=0.;
   /* initialize sums and variances */

    sum1[0]=0.;sum1[1]=0.;sum1[2]=0.;sum1[3]=0.;
    sum2[0]=0.;sum2[1]=0.;sum2[2]=0.;sum2[3]=0.;
    sum3[0]=0.;sum3[1]=0.;sum3[2]=0.;sum3[3]=0.;
    ssq1[0]=0.;ssq1[1]=0.;ssq1[2]=0.;ssq1[3]=0.;
    ssq2[0]=0.;ssq2[1]=0.;ssq2[2]=0.;ssq2[3]=0.;
    ssq3[0]=0.;ssq3[1]=0.;ssq3[2]=0.;ssq3[3]=0.;
    bigwt[0]=0.;bigwt[1]=0.;bigwt[2]=0.;bigwt[3]=0.;jgmax=0.;
    w1[0]=0.;w1[1]=0.;w1[2]=0.;w1[3]=0.;
    w2[0]=0.;w2[1]=0.;w2[2]=0.;w2[3]=0.;
    double Born  = 0.;

    cout << "Integrating..." << endl;

    for (int d = 0; d < nev;  d++) {
      //  cout << d << endl;
       do {
    fv = 0;   /*initialize cross-section for Born process*/	 
    fv1 = 0;
    fv2 = 0;
    if (!user.zerowidth()) {
      Q=sqrt(mp*Yzz*tan(rhomin+random(seed)*(rhomax-rhomin))+pow(mp,2));} else
	{Q=mp;}
      yz=log(pow(Q,2)/S2)*random(seed)-0.5*log(pow(Q,2)/S2);
      xq=sqrt(pow(Q,2)/S2)*exp(yz);
      xqb=sqrt(pow(Q,2)/S2)*exp(-yz);
      int e;
      if (proc==23) {
    for (e=1; e<6; e++) {
      fq1[e] = xfx(xq,Q,e);
      fqb2[e]= xfx(xqb,Q,sgn*e);
      fqb1[e] =xfx(xq,Q,-e);
      fq2[e] = xfx(xqb,Q,-sgn*e);
      ffq1[e] = fq1[e]*fqb2[e];
      ffq2[e] = fqb1[e]*fq2[e];
      
      fv1 += va[e]*ffq1[e];
      fv2 += va[e]*ffq2[e];
      }} else if(proc==24) {
	for(e=1; e<5; e++) {
	  fq1[e] = xfx(xq,Q,e)*xfx(xqb,Q,sgn*(e+1))+xfx(xq,Q,-e)*xfx(xqb,Q,-sgn*(e+1))+xfx(xqb,Q,sgn*e)*xfx(xq,Q,e+1)+xfx(xqb,Q,-sgn*e)*xfx(xq,Q,-(e+1));
	}
        int ee=1;
	for(e=5; e<7; e++) {
	  fq1[e] = xfx(xq,Q,ee)*xfx(xqb,Q,sgn*(ee+3))+xfx(xq,Q,-ee)*xfx(xqb,Q,-sgn*(ee+3))+xfx(xqb,Q,sgn*ee)*xfx(xq,Q,ee+3)+xfx(xqb,Q,-sgn*ee)*xfx(xq,Q,-(ee+3));
	ee+=1;}
	for (e=1; e<7; e++) {
	fv1 += vq[e]*fq1[e];
	fv2 = 0;}
      }
      fv = (fv1 + fv2)/(xq*xqb);
      
      } while (random(seed) > fv/fvmax);
    
    /* random numbers for phase space integrations */

  double r1=random(seed);
  double r2=random(seed);
  double r3=random(seed);
  double r4=random(seed);
     x = r1;
     y = 2.*r2 - 1.;
    w1[0] = DYBC(x,y,r3,r4);
    w1[1] = DYD(x,y,r3,r4);
    w1[2] = DYGBC(x,y,r3,r4);
    w1[3] = DYGD(x,y,r3,r4);
   
 /* Mandelstam variables s,t,u  and scale for NLO emission*/

  double s = pow(Q,2)/x;
  double t = -s*(1.-x)*(1.-y)/2.;
  double u = -s*(1.-x)*(1.+y)/2.;
  double scale = sqrt(u*t/s+pow(Q,2));
   if (scale < sqrt(getQ2min(0))) {scale = sqrt(getQ2min(0));}
  /*NLO momentum fractions x1,x2 */
  double x1 = xq*sqrt((2.-(1.-x)*(1.-y))/((2.-(1.-x)*(1.+y))*x));
  double x2 = xqb*sqrt((2.-(1.-x)*(1.+y))/((2.-(1.-x)*(1.-y))*x));

// Terms for calculating NLO pdf factors = combinations of sums of vv1..vv8.

  double fg1,fg2,fn1[6],fnb1[6],fn2[6],fnb2[6],fs1[6],fsb1[6];
  double fs2[6],fsb2[6],ffn1[6],ffn2[6],fgn1[6],fgn2[6],fgnb1[6],fgnb2[6];
  double ffs1[6],ffs2[6],v1[6],v2[6],v3[6],v4[6],v5[6],v6[6],v7[6],v8[6];
  double fqq[26],fq[26],fw[26];
  /* Calculate pdf factors for NLO emission. jqq,jg, jq are 
     the ratios of the NLO pdf factors to the Born level pdf factors. 
     jqq=[q(x1,scale)*q(x2,scale)+x1<->x2]/[q(x1,Q)*q(x2,Q)]. 
     jg=[g(x1,scale)*q(x2,scale)+x1<->x2]/[q(x1,Q)*q(x2,Q)].
    jq=[q(x1,Q)*q(x2,Q)+x1<->x2]/[q(x1,Q)*q(x2,Q)]*/
 // initialize sums

  double vv1=0;double vv2=0;double vv3=0;double vv4=0;
  double vv5=0;double vv6=0;double vv7=0;double vv8=0;
  double jqq=0;double jq=0;double jg=0;


  if(proc==23) {
  for (int ix=1; ix<6; ix++) {
    fg1 = xfx(x1,scale,0);
    fg2 = xfx(x2,scale,0);

    fq1[ix] = xfx(xq,Q,ix);
    fqb2[ix] = xfx(xqb,Q,sgn*ix);

    fs1[ix] =  xfx(xq,scale,ix);
    fsb2[ix] = xfx(xqb,scale,sgn*ix);
    fn1[ix] = xfx(x1,scale,ix);
    fnb2[ix] = xfx(x2,scale,sgn*ix);
    fqb1[ix] = xfx(xq,Q,-ix);
    fq2[ix] = xfx(xqb,Q,-sgn*ix);

    fsb1[ix] = xfx(xq,scale,-ix);
    fs2[ix] = xfx(xqb,scale,-sgn*ix);
    fnb1[ix] = xfx(x1,scale,-ix);
    fn2[ix] = xfx(x2,scale,-sgn*ix);

    ffq1[ix] = fq1[ix]*fqb2[ix];
    ffq2[ix] = fqb1[ix]*fq2[ix]; 
    ffn1[ix] = fn1[ix]*fnb2[ix];
    ffn2[ix] = fnb1[ix]*fn2[ix]; 
    ffs1[ix] = fs1[ix]*fsb2[ix];
    ffs2[ix] = fsb1[ix]*fs2[ix];
    fgn1[ix] = fg2*fn1[ix];
    fgn2[ix] = fg1*fn2[ix];
    fgnb1[ix] = fg2*fnb1[ix];
    fgnb2[ix] = fg1*fnb2[ix];
    v1[ix] = ffn1[ix]/ffq1[ix];
    v2[ix] = ffn2[ix]/ffq2[ix];
    v3[ix] = fgn1[ix]/ffq1[ix];
    v4[ix] = fgn2[ix]/ffq2[ix];
    v5[ix] = fgnb1[ix]/ffq2[ix];
    v6[ix] = fgnb2[ix]/ffq1[ix];
    v7[ix] = ffs1[ix]/ffq1[ix];
    v8[ix] = ffs2[ix]/ffq2[ix];
    if (ffq1[ix] == 0) {
    v1[ix] = 0;
    v3[ix] = 0;
    v6[ix] = 0;
    v7[ix] = 0;}
    if (ffq2[ix] == 0) {
    v2[ix] = 0;
    v4[ix] = 0; 
    v5[ix] = 0;
    v8[ix] = 0;}
    vv1 += v1[ix];
    vv2 += v2[ix];
    vv3 += v3[ix];
    vv4 += v4[ix];
    vv5 += v5[ix];
    vv6 += v6[ix];
    vv7 += v7[ix];
    vv8 += v8[ix];

  }

        jqq = vv1+vv2;
        jg = vv3+vv4+vv5+vv6;
        jq = vv7+vv8;
  }
    else if (proc==24) {
      int ix, iy;
    for (ix=1; ix<5; ix++) {
      fw[ix]=xfx(x1,scale,ix)*xfx(x2,scale,sgn*(ix+1));
      fq[ix]=xfx(xq,Q,ix)*xfx(xqb,Q,sgn*(ix+1));
      fqq[ix]=xfx(xq,scale,ix)*xfx(xqb,scale,sgn*(ix+1));
      if(fq[ix]!=0) {
      jqq+=fw[ix]/fq[ix];
      jq+=fqq[ix]/fq[ix];}
}
    int f=1;
    for (ix=5; ix<9; ix++) {
      fw[ix]=xfx(x1,scale,-f)*xfx(x2,scale,-sgn*(f+1));
      fq[ix]=xfx(xq,Q,-f)*xfx(xqb,Q,-sgn*(f+1));
      fqq[ix]=xfx(xq,scale,-f)*xfx(xqb,scale,-sgn*(f+1));
      if (fq[ix]!=0) {
      jqq+=fw[ix]/fq[ix];
      jq+=fqq[ix]/fq[ix];}
      f+=1;}
    f=1;
    for (ix=9; ix<13; ix++) {
      fw[ix]=xfx(x2,scale,sgn*f)*xfx(x1,scale,(f+1));
      fq[ix]=xfx(xqb,Q,sgn*f)*xfx(xq,Q,(f+1));
      fqq[ix]=xfx(xqb,scale,sgn*f)*xfx(xq,scale,(f+1));
      if (fq[ix]!=0) {
      jqq+=fw[ix]/fq[ix];
      jq+=fqq[ix]/fq[ix];}
      f+=1;}
    f=1;
    for (ix=13; ix<17; ix++) {
      fw[ix]=xfx(x2,scale,-sgn*f)*xfx(x1,scale,-(f+1));
      fq[ix]=xfx(xqb,Q,-sgn*f)*xfx(xq,Q,-(f+1));
      fqq[ix]=xfx(xqb,scale,-sgn*f)*xfx(xq,scale,-(f+1));
      if (fq[ix]!=0) {
      jqq+=fw[ix]/fq[ix]; 
      jq+=fqq[ix]/fq[ix];}
      f+=1;}
    f=1;
    for(ix=17;ix<19;ix++) {
      fw[ix]=xfx(x1,scale,f)*xfx(x2,scale,sgn*(f+3));
      fq[ix]=xfx(xq,Q,f)*xfx(xqb,Q,sgn*(f+3));
      fqq[ix]=xfx(xq,scale,f)*xfx(xqb,scale,sgn*(f+3));
      if (fq[ix]!=0) {
      jqq+=fw[ix]/fq[ix];
      jq+=fqq[ix]/fq[ix];}
      f+=1;}
    f=1;
    for(ix=19;ix<21;ix++) {
      fw[ix]=xfx(x1,scale,-f)*xfx(x2,scale,-sgn*(f+3));
      fq[ix]=xfx(xq,Q,-f)*xfx(xqb,Q,-sgn*(f+3));
      fqq[ix]=xfx(xq,scale,-f)*xfx(xqb,scale,-sgn*(f+3));
      if (fq[ix]!=0) {
      jqq+=fw[ix]/fq[ix];
      jq+=fqq[ix]/fq[ix];}
      f+=1;} 
    f=1;
    for(ix=21;ix<23;ix++) {
      fw[ix]=xfx(x2,scale,sgn*f)*xfx(x1,scale,(f+3));
      fq[ix]=xfx(xqb,Q,sgn*f)*xfx(xq,Q,(f+3));
      fqq[ix]=xfx(xqb,scale,sgn*f)*xfx(xq,scale,(f+3));
      if (fq[ix]!=0) {
      jqq+=fw[ix]/fq[ix];
      jq+=fqq[ix]/fq[ix];}
      f+=1;} 
    f=1;
    for(ix=23;ix<25;ix++) {
      fw[ix]=xfx(x2,scale,-sgn*f)*xfx(x1,scale,-(f+3));
      fq[ix]=xfx(xqb,Q,-sgn*f)*xfx(xq,Q,-(f+3));
      fqq[ix]=xfx(xqb,scale,-sgn*f)*xfx(xq,scale,-(f+3));
      if (fq[ix]!=0) {
      jqq+=fw[ix]/fq[ix];
      jq+=fqq[ix]/fq[ix];}
      f+=1;}
    
     

       double  Vud=0.952;
       double  Vus=0.049;
       double  Vub=0.00001063;
       double  Vcd=0.049;
       double  Vcs=0.952;
       double  Vcb=0.00168; 
      
      jg = 0.; jq = 0.,jqq=0.;
      double CKM[6],CKM2[25];
      CKM[1]=Vud+Vcd;
      CKM[2]=Vud+Vus+Vub;
      CKM[3]=Vus+Vcs;
      CKM[4]=Vcb+Vcs+Vcd;
      CKM[5]=Vub+Vcb;
      for (iy=1; iy< 6 ; iy++){
	jg +=  CKM[iy]*(xfx(x2,scale,0)*xfx(x1,scale,iy)+xfx(x2,scale,0)*xfx(x1,scale,-iy)+xfx(x1,scale,0)*xfx(x2,scale,sgn*iy)+xfx(x1,scale,0)*xfx(x2,scale,-sgn*iy));
      }
      //   if (jg > jgmax) {jgmax = jg; cout << jgmax << '\t' << scale << endl;}
    
       CKM2[1]=Vud;CKM2[2]=Vus;CKM2[3]=Vcs;CKM2[4]=Vcb;CKM2[5]=Vud;CKM2[6]=Vus;CKM2[7]=Vcs;CKM2[8]=Vcb;CKM2[9]=Vud;CKM2[10]=Vus;CKM2[11]=Vcs;CKM2[12]=Vcb;CKM2[13]=Vud;CKM2[14]=Vus;CKM2[15]=Vcs;CKM2[16]=Vcb;
       CKM2[17]=Vcd;CKM2[18]=Vub;CKM2[19]=Vcd;CKM2[20]=Vub;CKM2[21]=Vcd;CKM2[22]=Vub;CKM2[23]=Vcd;CKM2[24]=Vub;
       for (iy=1;iy <25 ; iy++){
	jq += CKM2[iy]*fqq[iy];
	}
       for (iy=1;iy <25 ; iy++){
	jqq += CKM2[iy]*fw[iy];
	}
       double Bfq=0.;
       for (int iyy=1; iyy < 25; iyy++) {
	Bfq += CKM2[iyy]*fq[iyy];
      }
       jq = jq/Bfq; jqq=jqq/Bfq; jg = jg/Bfq;
 
    }
      
    double scale2 = sqrt(pow(Q,2)+u*t/s);
    double as = alphasPDF(scale2);
     double wgt1, wgt2, wgt3, wgt4, wgt5;
      wgt3 = jq;
      jqq = jqq*as*CF/(2.*pi);
      jq = jq*as*CF/(2.*pi);
      jg = jg*as*TF/(2.*pi);
      if (scheme=="DIS") {
	wgt1 = -6.-4.*x;
       wgt2 = -5.*x+4.5*pow(x,2)+1.5+(pow(x,2)+pow((1.-x),2))*log((1.-x)*pow(Q/scale,2));
      wgt4 = 1.+4.*pow(pi,2)/3.+3.*log(pow(Q/scale,2));
       wgt5 = (3.+2.*(1.+pow(x,2))*log(pow(Q/scale,2)))*(x*jqq-jq)/(1.-x)+2.*(1.+pow(x,2))*log(1.-x)/(1.-x)*jqq*x-4.*log(1.-x)/(1.-x)*jq;
} else if (scheme=="MSbar") {
      wgt1 = -2.*(1.+pow(x,2))/(1.-x)*log(x);
      wgt2 = 0.5+3.*x-3.5*pow(x,2)+(pow(x,2)+pow(1.-x,2))*log(pow(1.-x,2)/x*pow(Q/scale,2));
      wgt4 = -8.+2.*pow(pi,2)/3.+3.*log(pow(Q/scale,2)) ;
      wgt5 = (2.*(1.+pow(x,2))*log(pow(Q/scale,2)))*(x*jqq-jq)/(1.-x)+4.*(1+pow(x,2))*log(1.-x)/(1.-x)*jqq*x-8.*log(1.-x)/(1.-x)*jq;}
      /* Calculate integrals p1,p2 and calculte maximum weights bigwt in
        phase space regions */
      
         for (int iz=0; iz<2; iz++) {
          w1[iz] = w1[iz]*jqq;
         
          if (w1[iz] != 0.) {
	   w2[iz] = wgt3-w1[iz]+(wgt1*jqq*x+wgt4*jq)
                    +wgt5;} else {
 	     w2[iz] = 0.;}
          if (abs(w1[iz]) > abs(bigwt[iz])) {bigwt[iz]=w1[iz];}}
	 
	 
       for (int ia=2; ia<4; ia++) {
       	 w1[ia] = w1[ia]*jg;
    	 if (w1[ia] != 0.) {
	 w2[ia] = -w1[ia]+wgt2*jg*x;} else {
	     w2[ia] = 0.;}
         if (fabs(w1[3]) > fabs(bigwt[3])) {bigwt[ia]=w1[ia];
	
}}
	 for (int ib=0; ib<4; ib++) {
         sum1[ib] += w1[ib];
	 if (w2[ib] > 0) {
	   sum2[ib] += w2[ib];
	   ssq1[ib] += pow(w1[ib],2);
	   ssq2[ib] += pow(w2[ib],2);
	 } else {sum3[ib]+=w2[ib];
          ssq3[ib] += pow(w2[ib],2);} }
	 Born += wgt3;
	 cout << "Integrated event : " << "\t" << d+1 << "\t" << "of" << "\t" << nev << "\r" << flush;  
    }
         cout << endl;
    double xsec = 0.;
    double xsecabs = 0.;
       for (int ic=0; ic<4; ic++) {
         p1[ic] = sum1[ic]/nev;
	 p2[ic] = sum2[ic]/nev;
	 p3[ic] = sum3[ic]/nev;
	 sd1[ic] = sqrt((ssq1[ic]/nev-pow(p1[ic],2))/nev);
         sd2[ic] = sqrt((ssq2[ic]/nev-pow(p2[ic],2))/nev);
	 sd3[ic] = sqrt((ssq3[ic]/nev-pow(p3[ic],2))/nev);
	 xsec+=p1[ic]+p2[ic]+p3[ic];
	 xsecabs+=fabs(p1[ic])+fabs(p2[ic])+fabs(p3[ic]);
	 //Output to screen the integrals, errors and maximum weights
	 	 cout << p1[ic] << "/" << p2[ic] <<"/"  << p3 [ic]<<"/" << bigwt[ic] << endl;
	 //	 	 cout << Born << endl;
         }
       cout << "k factor = " << xsec <<  endl;
       	 cout << "norm factor = " << xsecabs/xsec << endl;
	 
	 
       //*******************************************************************
       ofstream outdata, outdata1;
       if (MCNLO){
	   outdata.open("MCDYPP.dat", ios::trunc);//output for MC@NLO events
       } else {outdata.open("PWDYPP.dat", ios::trunc);} //output for POWHEG events

       int Jt,Jtt; /*phase space identifiers*/ 
       Jtt = 0;
       double cst[4];
       double cstt=0.;double ratio=0.;
       for(int ix=0; ix<4 ; ix++){
	 cst[ix]=fabs(p1[ix])+fabs(p2[ix])+fabs(p3[ix]);
       cstt += cst[ix];
       ratio += p1[ix]+p2[ix]+p3[ix];}
    

       
    outdata << "<LesHouchesEvents version =\"1.0\">" << endl; outdata << "<!--" << endl;
    if (MCNLO){outdata << "NORMFACTOR =" <<"\t"<< xsecabs/xsec << endl;}else {outdata << "NORMFACTOR =" <<"\t"<< 1.0 << endl;}
    if (MCNLO) {outdata << "MC@NLO" << "\t" ;} else {outdata << "POWHEG" << "\t" ;}
    if (proc==23) {outdata << "Z boson production from" << "\t" ;} else {outdata << "W boson production from" << "\t" ;}
    if (pp && user.cme() == 14000.) {outdata << "proton-proton collisions at LHC";} else if (!pp && (user.cme() > 1700. ||user.cme() < 2000.)) {outdata << "proton-antiproton collisions at the Tevatron" << endl;}
     outdata << "File generated with DYPP.cxx" << endl; outdata << "-->" << endl;
     outdata << "<init>" << endl;
     int WGTOPT = 3;
     if (!MCNLO) {WGTOPT = 1;}
     outdata << "\t2212\t" << "-2212\t" << user.cme()/2. << "\t" << user.cme()/2.  << "\t" << "0 \t 0 \t 7\t 7 \t" << WGTOPT << "\t 1" << endl;
     double csn;
     if (user.cme()==1800.) {
       if (proc==23) {csn=248.;} else if (proc==24) {csn=2500.;}}
     if (user.cme()==1960.) {
       if (proc==23) {csn=264.9;} else if (proc==24) {csn=2865.2;}}
     if (user.cme()==14000.) {
       if (proc==23) {csn=1860.;} else if (proc==24) {csn=20300.;}}
     
     outdata << "\t" << csn  << "\t" << 0.000000 << "\t1.00000 \t" << proc << endl;     outdata << "</init>" << endl;

/**********************************************************************************************/

     cout << "Generating events..." << endl;

   for (int d = 0; d < nevg; d++) {
    
     fv = 0.; /*cross section for Born process*/
     // POWHEG stuff
      bool emit = false;
      double wgt, bigwgt,kt;
      double s = 0,t = 0, u = 0,x1 = 0,x2 = 0;/* Mandelstam variables s,t,u  and NLO momentum fractions x1,x2*/
	// Assign particle IDs to quarks according to Born pdf factors
        int ID[2],ID1, ID2;//ID2;
      /************************************************/
      do {
	rmwgt=0.;
     do {
    fv1 = 0;
    fv2 = 0;
    if (!user.zerowidth()) {
      Q=sqrt(mp*Yzz*tan(rhomin+random(seed)*(rhomax-rhomin))+pow(mp,2));} else
	{Q=mp;}
      yz=log(pow(Q,2)/S2)*random(seed)-0.5*log(pow(Q,2)/S2);
      xq=sqrt(pow(Q,2)/S2)*exp(yz);
      xqb=sqrt(pow(Q,2)/S2)*exp(-yz);
      int e;
      if (proc==23) {
    for (e=1; e<6; e++) {
      fq1[e] = xfx(xq,Q,e);
      fqb2[e]= xfx(xqb,Q,sgn*e);
      fqb1[e] =xfx(xq,Q,-e);
      fq2[e] = xfx(xqb,Q,-sgn*e);
      ffq1[e] = fq1[e]*fqb2[e];
      ffq2[e] = fqb1[e]*fq2[e];
      
      fv1 += va[e]*ffq1[e];
      fv2 += va[e]*ffq2[e];
      }} else if(proc==24) {
	for(e=1; e<5; e++) {
	  fq1[e] = xfx(xq,Q,e)*xfx(xqb,Q,sgn*(e+1))+xfx(xq,Q,-e)*xfx(xqb,Q,-sgn*(e+1))+xfx(xqb,Q,sgn*e)*xfx(xq,Q,e+1)+xfx(xqb,Q,-sgn*e)*xfx(xq,Q,-(e+1));
	}
        int ee=1;
	for(e=5; e<7; e++) {
	  fq1[e] = xfx(xq,Q,ee)*xfx(xqb,Q,sgn*(ee+3))+xfx(xq,Q,-ee)*xfx(xqb,Q,-sgn*(ee+3))+xfx(xqb,Q,sgn*ee)*xfx(xq,Q,(ee+3))+xfx(xqb,Q,-sgn*ee)*xfx(xq,Q,-(ee+3));
	ee+=1;}
	for (e=1; e<7; e++) {
	fv1 += vq[e]*fq1[e];
	fv2 = 0;}
      }
      /*recall PDFs give xf(x,Q) so need to divide out xq*xqb in product */

      fv = (fv1 + fv2)/(xq*xqb);;
      
     } while (random(seed) > fv/fvmax);

// Assign particle IDs to quarks according to Born pdf factors
    
      
      double g1 = 0;
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
       ID[0] = ID1; ID[1] = -ID1;}
       else {g = g2;
       ID[0] = ID2; ID[1]= -ID2; }}

      else if (proc==24) {
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
	   ID[0]=ID1;ID[1]=-(ID1+1);}

	 else if (random(seed) < xfx(xq,Q,-ID1)*xfx(xqb,Q,-sgn*(ID1+1))/(fq1[ID1]-xfx(xq,Q,ID1)*xfx(xqb,Q,sgn*(ID1+1)))) {ID[0]=-ID1;ID[1]=ID1+1;}

         else if (random(seed) < xfx(xq,Q,ID1+1)*xfx(xqb,Q,sgn*ID1)/(xfx(xq,Q,ID1+1)*xfx(xqb,Q,sgn*ID1)+xfx(xq,Q,-(ID1+1))*xfx(xqb,Q,-sgn*ID1))) {ID[0]=ID1+1;ID[1]=-ID1;}

	 else {ID[0]=-(ID1+1);ID[1]=ID1;} } 

	 if (ID1 > 4) {
	 
	  if (random(seed) < xfx(xq,Q,ID1-1)*xfx(xqb,Q,sgn*(ID1-4))/fq1[ID1])   {   ID[0]=ID1-1;ID[1]=-(ID1-4);}

	 else if (random(seed) < xfx(xq,Q,-(ID1-1))*xfx(xqb,Q,-sgn*(ID1-4))/(fq1[ID1]-xfx(xq,Q,ID1-1)*xfx(xqb,Q,sgn*(ID1-4)))) {ID[0]=-(ID1-1);ID[1]=ID1-4;}

         else if (random(seed) < xfx(xq,Q,ID1-4)*xfx(xqb,Q,sgn*(ID1-1))/(xfx(xq,Q,ID1-4)*xfx(xqb,Q,sgn*(ID1-1))+xfx(xq,Q,-(ID1-4))*xfx(xqb,Q,-sgn*(ID1-1))))
          {ID[0]=ID1-4;ID[1]=-(ID1-1);}

	 else {ID[0]=-(ID1-4);ID[1]=ID1-1;} }
              
      }
      
 /* Select phase space according to contribution to cross-section. 
     Jt=odd generates NLO events. Jt=even generates Born events.*/
	bigwgt=0; 
   if (random(seed) < cst[0]/cstt) {
    if (random(seed) < fabs(p1[0])/cst[0]) {
       Jt = 1;
       bigwgt=bigwt[0];}
    else { Jt = 2;} }
  else if (random(seed) < cst[1]/(cstt-cst[0])) {
       if (random(seed) < fabs(p1[1])/cst[1]) {
       Jt = 3;
       bigwgt=bigwt[1];}
       else {Jt = 4;} }
   else if (random(seed) < cst[2]/(cst[3]+cst[2])) {
       if (random(seed) < fabs(p1[2])/cst[2]) {
       Jt = 5;
       bigwgt=bigwt[2];}
       else {Jt = 6;} }
   else if (random(seed) < fabs(p1[3])/cst[3]) {
       Jt = 7;
       bigwgt=bigwt[3];}
       else {Jt = 8;}
    
    wgt=0.;
    

     if (Jt==1||Jt==3||Jt==5||Jt==7) { 

        /* Hereafter wgt=integrand*/
      

/* random numbers for event generation */

  double r1=random(seed);
  double r2=random(seed);
  double r3=random(seed);
  double r4=random(seed);
  double r5=random(seed);
  double r6=random(seed);
    x = r1;
    y = 2.*r2 - 1.;
      if (Jt==1) {
	wgt=DYBC(x,y,r3,r4);}
      if (Jt==3) {
	wgt=DYD(x,y,r3,r4);}
      if (Jt==5) {
	wgt=DYGBC(x,y,r3,r4);}
      if (Jt==7) {
	wgt=DYGD(x,y,r3,r4);}
  
 /* Mandelstam variables s,t,u  and scale for NLO emission*/  
   double as,scale;
   s = pow(Q,2)/x;
   t = -s*(1.-x)*(1.-y)/2.;
   u = -s*(1.-x)*(1.+y)/2.;
  scale = sqrt(u*t/s+pow(Q,2));
  //  if (scale < sqrt(pdf.getQ2min(0))) {scale = sqrt(pdf.getQ2min(0));} // pdfs undefined for scale < Q2min.
   if (scale < sqrt(getQ2min(0))) {scale = sqrt(getQ2min(0));} 
 /*NLO momentum fractions x1,x2 */
  x1 = xq*sqrt((2.-(1.-x)*(1.-y))/((2.-(1.-x)*(1.+y))*x));
  x2 = xqb*sqrt((2.-(1.-x)*(1.+y))/((2.-(1.-x)*(1.-y))*x));
  r1=random(seed);
  r2=random(seed);
  r3=random(seed);
  r4=random(seed);
  r5=random(seed);
  r6=random(seed);

// Terms for calculating NLO pdf factors = combinations of sums of vv1..vv8.
  
  double fg1,fg2,fn1[6],fnb1[6],fn2[6],fnb2[6],fs1[6],fsb1[6];
  double fs2[6],fsb2[6],ffn1[6],ffn2[6],fgn1[6],fgn2[6],fgnb1[6];
  double fgnb2[6],ffs1[6],ffs2[6];
  double v1[6],v2[6],v3[6],v4[6],v5[6],v6[6],v7[6],v8[6];
  double fqq[26],fq[26],fw[26];
  
  /* Calculate pdf factors for NLO emission. 
   jqq,jg, jq are the ratios of the NLO pdf factors 
   to the Born level pdf factors. 
   jqq=[q(x1,scale)*q(x2,scale)+x1<->x2]/[q(x1,Q)*q(x2,Q)]. 
   jg=[g(x1,scale)*q(x2,scale)+x1<->x2]/[q(x1,Q)*q(x2,Q)].
   jq=[q(x1,Q)*q(x2,Q)+x1<->x2]/[q(x1,Q)*q(x2,Q)] */
  //  initialize sums.
  
  double vv1=0;double vv2=0;double vv3=0;double vv4=0;
  double vv5=0;double vv6=0;double vv7=0;double vv8=0;
  double jqq=0;double jq=0;double jg=0;
   if(proc==23) {
  for (int ix=1; ix<6; ix++) {
    int f = ix ;
    fg1 = xfx(x1,scale,0);
    fg2 = xfx(x2,scale,0);

    fq1[ix] = xfx(xq,Q,f);
    fqb2[ix] = xfx(xqb,Q,sgn*f);

    fs1[ix] =  xfx(xq,scale,f);
    fsb2[ix] = xfx(xqb,scale,sgn*f);
    fn1[ix] = xfx(x1,scale,f);
    fnb2[ix] = xfx(x2,scale,sgn*f);
    f = -ix;

    fqb1[ix] = xfx(xq,Q,f);
    fq2[ix] = xfx(xqb,Q,sgn*f);

    fsb1[ix] = xfx(xq,scale,f);
    fs2[ix] = xfx(xqb,scale,sgn*f);
    fnb1[ix] = xfx(x1,scale,f);
    fn2[ix] = xfx(x2,scale,sgn*f);
    ffq1[ix] = fq1[ix]*fqb2[ix];
    ffq2[ix] = fqb1[ix]*fq2[ix]; 
    ffn1[ix] = fn1[ix]*fnb2[ix];
    ffn2[ix] = fnb1[ix]*fn2[ix]; 
    ffs1[ix] = fs1[ix]*fsb2[ix];
    ffs2[ix] = fsb1[ix]*fs2[ix];
    fgn1[ix] = fg2*fn1[ix];
    fgn2[ix] = fg1*fn2[ix];
    fgnb1[ix] = fg2*fnb1[ix];
    fgnb2[ix] = fg1*fnb2[ix];
    v1[ix] = ffn1[ix]/ffq1[ix];
    v2[ix] = ffn2[ix]/ffq2[ix];
    v3[ix] = fgn1[ix]/ffq1[ix];
    v4[ix] = fgn2[ix]/ffq2[ix];
    v5[ix] = fgnb1[ix]/ffq2[ix];
    v6[ix] = fgnb2[ix]/ffq1[ix];
    v7[ix] = ffs1[ix]/ffq1[ix];
    v8[ix] = ffs2[ix]/ffq2[ix];
    if (ffq1[ix] == 0) {
    v1[ix] = 0;
    v3[ix] = 0;
    v6[ix] = 0;
    v7[ix] = 0;}
    if (ffq2[ix] == 0) {
    v2[ix] = 0;
    v4[ix] = 0; 
    v5[ix] = 0;
    v8[ix] = 0;}
    vv1 += v1[ix];
    vv2 += v2[ix];
    vv3 += v3[ix];
    vv4 += v4[ix];
    vv5 += v5[ix];
    vv6 += v6[ix];
    vv7 += v7[ix];
    vv8 += v8[ix];
  }
  double Bfq = 0.;
  jqq=0.;jq=0.; jg=0.;
      for (int ix=1 ; ix < 6 ; ix++) {
	Bfq += ffq1[ix]+ffq2[ix];
	jqq += ffn1[ix]+ffn2[ix];
	jq += ffs1[ix]+ffs2[ix];
	jg += fgn1[ix]+fgn2[ix]+fgnb1[ix]+fgnb2[ix];
      }
    jq = jq/Bfq;
    jqq = jqq/Bfq;
    jg = jg/Bfq;

   }
   
else if (proc==24) {
      int ix, iy;
      
    for (ix=1; ix<5; ix++) {
      fw[ix]=xfx(x1,scale,ix)*xfx(x2,scale,sgn*(ix+1));
      fq[ix]=xfx(xq,Q,ix)*xfx(xqb,Q,sgn*(ix+1));
      fqq[ix]=xfx(xq,scale,ix)*xfx(xqb,scale,sgn*(ix+1));
 
}
    int f=1;
    for (ix=5; ix<9; ix++) {
      fw[ix]=xfx(x1,scale,-f)*xfx(x2,scale,-sgn*(f+1));
      fq[ix]=xfx(xq,Q,-f)*xfx(xqb,Q,-sgn*(f+1));
      fqq[ix]=xfx(xq,scale,-f)*xfx(xqb,scale,-sgn*(f+1));
  
      f+=1;}
    f=1;
    for (ix=9; ix<13; ix++) {
      fw[ix]=xfx(x2,scale,sgn*f)*xfx(x1,scale,(f+1));
      fq[ix]=xfx(xqb,Q,sgn*f)*xfx(xq,Q,(f+1));
      fqq[ix]=xfx(xqb,scale,sgn*f)*xfx(xq,scale,(f+1));
  
      f+=1;}
    f=1;
    for (ix=13; ix<17; ix++) {
      fw[ix]=xfx(x2,scale,-sgn*f)*xfx(x1,scale,-(f+1));
      fq[ix]=xfx(xqb,Q,-sgn*f)*xfx(xq,Q,-(f+1));
      fqq[ix]=xfx(xqb,scale,-sgn*f)*xfx(xq,scale,-(f+1));
 
      f+=1;}
    
    f=1;
    for(ix=17;ix<19;ix++) {
      fw[ix]=xfx(x1,scale,f)*xfx(x2,scale,sgn*(f+3));
      fq[ix]=xfx(xq,Q,f)*xfx(xqb,Q,sgn*(f+3));
      fqq[ix]=xfx(xq,scale,f)*xfx(xqb,scale,sgn*(f+3));
 
      f+=1;}
    
    f=1;
    for(ix=19;ix<21;ix++) {
      fw[ix]=xfx(x1,scale,-f)*xfx(x2,scale,-sgn*(f+3));
      fq[ix]=xfx(xq,Q,-f)*xfx(xqb,Q,-sgn*(f+3));
      fqq[ix]=xfx(xq,scale,-f)*xfx(xqb,scale,-sgn*(f+3));
   
      f+=1;} 
    
    f=1;
    for(ix=21;ix<23;ix++) {
      fw[ix]=xfx(x2,scale,sgn*f)*xfx(x1,scale,(f+3));
      fq[ix]=xfx(xqb,Q,sgn*f)*xfx(xq,Q,(f+3));
      fqq[ix]=xfx(xqb,scale,sgn*f)*xfx(xq,scale,(f+3));
  
      f+=1;} 
    
    f=1;
    for(ix=23;ix<25;ix++) {
      fw[ix]=xfx(x2,scale,-sgn*f)*xfx(x1,scale,-(f+3));
      fq[ix]=xfx(xqb,Q,-sgn*f)*xfx(xq,Q,-(f+3));
      fqq[ix]=xfx(xqb,scale,-sgn*f)*xfx(xq,scale,-(f+3));
 
      f+=1;}

      double  Vud=0.952;
      double  Vus=0.049;
      double  Vub=0.00001063;
      double  Vcd=0.049;
      double  Vcs=0.952;
      double  Vcb=0.00168;   
       jg = 0.; jq = 0.,jqq=0.;
      double CKM[6],CKM2[25];
      CKM[1]=Vud+Vcd;
      CKM[2]=Vud+Vus+Vub;
      CKM[3]=Vus+Vcs;
      CKM[4]=Vcb+Vcs+Vcd;
      CKM[5]=Vub+Vcb;
      for (iy=1; iy< 6 ; iy++){
	jg +=  CKM[iy]*(xfx(x2,scale,0)*xfx(x1,scale,iy)+xfx(x2,scale,0)*xfx(x1,scale,-iy)+xfx(x1,scale,0)*xfx(x2,scale,sgn*iy)+xfx(x1,scale,0)*xfx(x2,scale,-sgn*iy));
      }
      //  cout << "jg1" << jg << endl;
      CKM2[1]=Vud;CKM2[2]=Vus;CKM2[3]=Vcs;CKM2[4]=Vcb;CKM2[5]=Vud;CKM2[6]=Vus;CKM2[7]=Vcs;CKM2[8]=Vcb;CKM2[9]=Vud;CKM2[10]=Vus;CKM2[11]=Vcs;CKM2[12]=Vcb;CKM2[13]=Vud;CKM2[14]=Vus;CKM2[15]=Vcs;CKM2[16]=Vcb;
      CKM2[17]=Vcd;CKM2[18]=Vub;CKM2[19]=Vcd;CKM2[20]=Vub;CKM2[21]=Vcd;CKM2[22]=Vub;CKM2[23]=Vcd;CKM2[24]=Vub;
       for (iy=1;iy <25 ; iy++){
	jq += CKM2[iy]*fqq[iy];
	}
       for (iy=1;iy <25 ; iy++){
	jqq += CKM2[iy]*fw[iy];
	}
       double Bfq=0.;
       for (int iyy=1; iyy < 25; iyy++) {
	Bfq += CKM2[iyy]*fq[iyy];
      }
       jq = jq/Bfq; jqq=jqq/Bfq;jg = jg/Bfq; 
      
      
      
}

   
    /*Select whether  quark-gluon (Jtt=1) or gluon-quark 
      (Jtt=2) according to pdf factor. */
      
       if (Jt==5||Jt==7){
	 
     double qg1 = xfx(x1,scale,ID[0]);
     double qg2 = xfx(x2,scale,0);
     double gq1 = xfx(x2,scale,-sgn*ID[1]);
     double gq2 = xfx(x1,scale,0);
     double jqg,jqbg;
     
     jqg=qg1*qg2;
     jqbg=gq1*gq2;
      Jtt=2;
      if (random(seed) < jqg/(jqg+jqbg)) {Jtt = 1;}
     
       }
       

       //   as = pdf.alphasPDF(scale); // NLO alphas evaluated at Q=scale
       as = alphasPDF(scale);
      // Multiply integrands by colour factors

      jq = jq*as*CF/(2.*pi);
      jg = jg*as*TF/(2.*pi);
      if (Jt==1||Jt==3) {
 	wgt=wgt*jqq;}
    
       if (Jt==5||Jt==7) {
 	wgt=wgt*jg;}
	//POWHEG stuff
	/***************************************************************/
	double wgt1, wgt2, wgt3, wgt4, wgt5;
      wgt3 = jq;
      jqq = jqq*as*CF/(2.*pi);
      jq = jq*as*CF/(2.*pi);
      jg = jg*as*TF/(2.*pi);
      if (scheme=="DIS") {
      wgt1 = -6.-4.*x;
       wgt2 = -5.*x+4.5*pow(x,2)+1.5+(pow(x,2)+pow((1.-x),2))*log((1.-x)*pow(Q/scale,2));
       wgt4 = 1.+4.*pow(pi,2)/3.+3.*log(pow(Q/scale,2));
       wgt5 = (3.+2.*(1.+pow(x,2))*log(pow(Q/scale,2)))*(x*jqq-jq)/(1.-x)+2.*(1.+pow(x,2))*log(1.-x)/(1.-x)*jqq*x-4.*log(1.-x)/(1.-x)*jq;
} else if (scheme=="MSbar") {
      wgt1 = -2.*(1.+pow(x,2))/(1.-x)*log(x);
      wgt2 = 0.5+3.*x-3.5*pow(x,2)+(pow(x,2)+pow(1.-x,2))*log(pow(1.-x,2)/x*pow(Q/scale,2));
      wgt4 = -8.+2.*pow(pi,2)/3.+3.*log(pow(Q/scale,2));
      wgt5 = (2.*(1.+pow(x,2))*log(pow(Q/scale,2)))*(x*jqq-jq)/(1.-x)+4.*(1+pow(x,2))*log(1.-x)/(1.-x)*jqq*x-8.*log(1.-x)/(1.-x)*jq;}
      /* Calculate integrals p1,p2 and calculte maximum weights bigwt in
        phase space regions */
      
         for (int iz=0; iz<2; iz++) {
          w1[iz] = w1[iz]*jqq;
         
          if (w1[iz] != 0.) {
	   w2[iz] = wgt3-w1[iz]+(wgt1*jqq*x+wgt4*jq)
                    +wgt5;} else {
 	     w2[iz] = 0.;}
        }
	 
	 // if (jg > jgmax) {jgmax = jg; cout << jgmax << endl;}
       for (int ia=2; ia<4; ia++) {
       	 w1[ia] = w1[ia]*jg;
    	 if (w1[ia] != 0.) {
	 w2[ia] = -w1[ia]+wgt2*jg*x;} else {
	     w2[ia] = 0.;}
        }
       for (int ibb=0; ibb < 4; ibb++) {
         rmwgt += w1[ibb]+w2[ibb];
	 
       } 
         if (rmwgt < 0.) {rmwgt = 0.;}

       if (!MCNLO) {wgt = rmwgt; bigwgt = bigrmwgt;}
       
       

       /***************************************************************/

       /* UNWEIGHTING. For Jt=1,3,5,7: generate event according to 
        integrand by comparing generated integrand to maximum value. 
       If it fails veto, go back and generate new event else proceed...*/
 
     } else {wgt = 1.; bigwgt=1.;}} while (random(seed) > fabs(wgt)/fabs(bigwgt));
       //POWHEG stuff
       if (!MCNLO) {Jt = 3;}

    // Assign weights for unweighted events
       
      if (Jt==1||Jt==5) {wgt=-1.;}
      if (Jt==3||Jt==7) {wgt=1.;}
     
     if (Jt==2) {if (random(seed) < fabs(p2[0])/(cst[0]-fabs(p1[0]))) {wgt=1.;} else {wgt=-1.;}}
     if (Jt==4) {if (random(seed) < fabs(p2[1])/(cst[1]-fabs(p1[1]))) {wgt=1.;} else {wgt=-1.;}}
     if (Jt==6) {if (random(seed) < fabs(p2[2])/(cst[2]-fabs(p1[2]))) {wgt=1.;} else {wgt=-1.;}}
     if (Jt==8) {if (random(seed) < fabs(p2[3])/(cst[3]-fabs(p1[3]))) {wgt=1.;} else {wgt=-1.;}}
     double z0, z1, pt1, xp, yp, qs, xp1, xp2;
     double xx, yy, ww, M, MT, V, Vp;
     int pro;
     pt1=0.; qs=0.; xp1=0.;xp2=0.;pro=0; 
    //POWHEG stuff
    if (!MCNLO) {
     double ll,csfac,bb,bp,Nq,eps,Q2,r,ptm,qq,p2m,n,pt2,pt;
     int Its, MaxIts, converged,pro, veto;
     ll=0.2; csfac=2.*1./(3.*pi);eps=1.0e-6;Nq=40./(2.*pi);bb=23./(12.*pi);bp=58./(46.*pi); Q2 = pow (Q,2);r =Q2/S2;
     ptm = pow(1.-r, 2)*S2/4.;
     qq = ptm+Q2;
     pt2 = 0.;
     do {
     if (pt2 < pow(ll,2)+0.001) {p2m = ptm;} else {p2m = pt2;}
     n = random(seed);
     pt2 = pow(ll,2)+0.001;
     converged = 0; Its=0; MaxIts=1000000;
  
     do {
       pt2 = pt2 - (G(pt2,Nq,bb,qq,ptm,ll) - G(p2m,Nq,bb,qq,ptm,ll)+log(n))/DG(pt2,Nq,bb,qq,ptm,ll);
       Its=Its+1;
       if (fabs(G(pt2,Nq,bb,qq,ptm,ll)-G(p2m,Nq,bb,qq,ptm,ll)+log(n)) < eps) {converged=1;}
     } while (converged ==0 );    
   
     kt = sqrt(pt2);
     double xp=pow((sqrt(1.+pt2/Q2)+kt/Q),2);
     double xm=pow((sqrt(1.+pt2/Q2)-kt/Q),2);
     double alphas=alphasPDF(kt) + pow(alphasPDF(kt),2)*(67./3.-pow(pi,2)-50./9.)/(4.*pi);
     V=Nq*alphas*log((sqrt(xp-r)+sqrt(xm-r))/(sqrt(xp-r)-sqrt(xm-r)));
     Vp = Nq/(2.*bb)*(log(qq/pt2))/log(pt2/pow(ll,2));
     if (random(seed) <  V/Vp && pt2 > pow(ll,2)+0.001) { 
       veto = 0;
      double Mmax=0.;
     for (int ix=0; ix < 1000 ; ix++) {
       yy = 2.*random(seed)-1.;
       xx = (-Q2-2.*pt2+pow(yy,2)*Q2+2.*sqrt(Q2*pt2+pow(pt2,2)-pt2*pow(yy,2)*Q2))/(Q2*(pow(yy,2)-1.));
       ww = -pt2/xx-2.*pt2;
       M = Nq*2.*pi/((1.-xx)*(1.-pow(yy,2))*ww);
       if (fabs(M) > fabs(Mmax)) {Mmax=M;}
     }
     do {
      yy = 2.*random(seed)-1.;
      xx = (-Q2-2.*pt2+pow(yy,2)*Q2+2.*sqrt(Q2*pt2+pow(pt2,2)-pt2*pow(yy,2)*Q2))/(Q2*(pow(yy,2)-1.));
      ww = -pt2/xx-2.*pt2;
      M = Nq*2.*pi/((1.-xx)*(1.-pow(yy,2))*ww);} while (random(seed) > fabs(M/Mmax));
  
      x1=xq*sqrt((2.-(1.-xx)*(1.-yy))/((2.-(1.-xx)*(1.+yy))*xx));
      x2=xqb*sqrt((2.-(1.-xx)*(1.+yy))/((2.-(1.-xx)*(1.-yy))*xx));
      if (x1 > 1.|| x2 > 1.) {yy=-yy;}
      x1=xq*sqrt((2.-(1.-xx)*(1.-yy))/((2.-(1.-xx)*(1.+yy))*xx));
      x2=xqb*sqrt((2.-(1.-xx)*(1.+yy))/((2.-(1.-xx)*(1.-yy))*xx));
      if (x1 < 1. && x2 <1.){
		pt = kt;
	if (pt < sqrt(getQ2min(0))) {pt = sqrt(getQ2min(0));}
	double rjqq = xfx(x1,pt,ID[0])*xfx(x2,pt,-sgn*ID[1])/(xfx(xq,Q,ID[0])*xfx(xqb,Q,-sgn*ID[1]));
        double rjqg = xfx(x1,pt,ID[0])*xfx(x2,pt,0)/(xfx(xq,Q,ID[0])*xfx(xqb,Q,-sgn*ID[1]));
	double rjqbg = xfx(x1,pt,0)*xfx(x2,pt,-sgn*ID[1])/(xfx(xq,Q,ID[0])*xfx(xqb,Q,-sgn*ID[1]));
       double Mqq=(4./3.*(pow(yy*(1.-xx),2)+pow(1.+xx,2))/((1.-xx)*(1.-pow(yy,2))))*rjqq; 
       double Mqg=(((3.+pow(yy,2))*pow(1.-xx,2)-2.*yy*(1.-pow(xx,2))+2.*(1.+pow(xx,2)))/(8.*(1.-yy)))*rjqg;
       double Mqbg=(((3.+pow(yy,2))*pow(1.-xx,2)+2.*yy*(1.-pow(xx,2))+2.*(1.+pow(xx,2)))/(8.*(1.+yy)))*rjqbg;
       MT = fabs(Mqq)+fabs(Mqg)+fabs(Mqbg);
      double MTT;
      if (random(seed) < fabs(Mqq)/fabs(MT)) {pro=1; MTT = Mqq;}else 
      if (random(seed) < fabs(Mqg)/(fabs(Mqg)+fabs(Mqbg))) {pro=2; MTT = Mqg;} else {pro=3; MTT = Mqbg;};
} else {MT =0.; ww = 1.;}
      
} else {veto = 1;} }while (veto == 1 || x1 > 1. || x2 > 1. || random(seed) > fabs(MT/ww)/fabs(M));

     
      bool vetoed = false;

     if (truncsh) {
     double mu = 0.;
     double alpham = 1.;  
     double Qg = 0.75;
     double qc = 0.631;
       z0=1.-0.5*(1.-xx)*(1.+yy);
     if (yy < 0.) {z0=1.-0.5*(1.-xx)*(1.-yy);} 
     int Jt = 1;
     if (random(seed) < pow(xq,2)/(pow(xq,2)+pow(xqb,2))) {Jt = 2;} 
       kt = sqrt(pt2);
       s = pow(Q,2)/xx;
       t = -s*(1.-xx)*(1.-yy)/2.;
       u = -s*(1.-xx)*(1.+yy)/2.;
       double Pz, Pzm, P, g, pt12, pdfr, alphatr, qi;
       P=0.;g=0.;pt12=0.;alphatr=0.;pdfr=0.;
        do {
	  emit = false;
	   qi = Q;
       if (vetoed == true) {qi = qs;}
       do {
       z1 = random(seed)*((2.*pow(qi,2)+pow(Qg,2)-Qg*sqrt(4.*pow(qi,2)+pow(Qg,2)))/(2.*pow(qi,2))-Qg/qi)+Qg/qi;
       Pz = (1.+pow(z1,2))/(z1*(1.-z1));
       Pzm = (1.+pow(((2.*pow(qi,2)+pow(Qg,2)-Qg*sqrt(4.*pow(qi,2)+pow(Qg,2)))/(2.*pow(qi,2))),2))/((2.*pow(qi,2)+pow(Qg,2)-Qg*sqrt(4.*pow(qi,2)+pow(Qg,2)))/(2.*pow(qi,2))*(1.-(2.*pow(qi,2)+pow(Qg,2)-Qg*sqrt(4.*pow(qi,2)+pow(Qg,2)))/(2.*pow(qi,2))));
       } while (random(seed) > Pz/Pzm);
       double qh = sqrt((pt2+z0*pow(Qg,2))/pow(1.-z0,2));
       double C=CF*alpham*log(1.-(2.*pow(qi,2)+pow(Qg,2)-Qg*sqrt(4.*pow(qi,2)+pow(Qg,2)))/(2.*pow(qi,2)))/pi;
       double D=exp(log(pow(qh/qc,2))*C);
       double E=exp(log(pow(qi/qc,2))*C);
       qs = 0.;
       if (random(seed) < 1. -E/D) {
	 do {qs = sqrt(exp(log(E/random(seed))/C)*pow(qc,2));} while (qs < qh);
	if(qs < qi) {
	  emit = true;}
       }
	if (emit) {
	  P=CF*((1.+pow(z1,2))/(1.-z1)-2.*pow(mu,2)/(z1*(1.-z1)*pow(qs,2)));
          g=2.*CF/(1.-z1);
          pt12=(pow((1.-z1)*qs,2)-z1*pow(Qg,2));
	 
	  if (pt12 > 0.) {
	    pt1= sqrt(pt12);
	    double kap=pt12/pow(qi*(1.-z1),2);
	    xp=z1/(1.+(1.-z1)*kap);
	    yp=1.-2.*(1.-z1)/(1.-xp);
	    if (yy < 0.) {yp=-1.+2.*(1.-z1)/(1.-xp);}
            xp1=xq*sqrt((2.-(1.-xp)*(1.-yp))/((2.-(1.-xp)*(1.+yp))*xp));
            xp2=xqb*sqrt((2.-(1.-xp)*(1.+yp))/((2.-(1.-xp)*(1.-yp))*xp)); 
	    pdfr = xp1*xp2*xfx(xp1,qs,ID[0])*xfx(xp2,qs,ID[1])/(xq*xqb*xfx(xq,qi,ID[0])*xfx(xqb,qi,ID[1]));
	    alphatr = alphasPDF(z1*(1-z1)*qs) + pow(alphasPDF(z1*(1-z1)*qs),2)*(67./3.-pow(pi,2)-50./9.)/(4.*pi);
	  }
	}
	if (emit == true && (pt12 > pt2 || pt12 < 0 || random(seed) > pdfr || random(seed) > P/g || random(seed) > alphatr/alpham)){vetoed = true;} else {vetoed = false;}
	} while (vetoed == true && qs < qi); }}

 /**********************************************************************/    

   
  /*Les Houches Accord event parameters:: PUP are the momenta, 
   ICOLUP is the colour code, MOTHUP is the mother information, ISTUP is 
   the status code. ISPINUP signifies spins. IDUP is the particle ID.
    NUP is the number of particles in the event. SCALUP is the scale of 
   the event.AQEDUP and AQCDUP are alphaem and alphas resp.XWGTUP is 
   the weight of the event=1 or -1 for unweighted events.*/

	double PUP[5][9],AQEDUP,AQCDUP,XWGTUP,SCALUP,mass[8];
	int ICOLUP[3][7],MOTHUP[3][7],ISTUP[7],ISPINUP[7],IDUP[7],NUP;

    XWGTUP = wgt;
  //Default parameters

  PUP[1][1]=0.;PUP[2][1]=0.;PUP[1][2]=0.;PUP[2][2]=0.;
  MOTHUP[1][1]=0;
  MOTHUP[2][1]=0;
  MOTHUP[1][2]=0;
  MOTHUP[2][2]=0;
  MOTHUP[1][3]=1;
  MOTHUP[2][3]=2;
  ICOLUP[1][3]=0;
  ICOLUP[2][3]=0;
  ISTUP[1]=-1;
  ISTUP[2]=-1;
  ISTUP[3]=2;
  mass[1]=0.;mass[2]=0.;mass[4]=0.;mass[5]=0.,mass[6]=0.; mass[7]=0.;
  AQEDUP=0.007297352;
  AQCDUP=0.118;
  
	if (Jt==1||Jt==3||Jt==5||Jt==7) { 
	 
       //transverse mmt of emission
	  if (MCNLO) {	  
	    kt = sqrt(pow(Q,2)*(1.-pow(y,2))*pow((1.-x),2)/(4.*x));}
       if (emit) { 
      x1=xp1*sqrt((2.-(1.-xx)*(1.-yy))/((2.-(1.-xx)*(1.+yy))*xx));
      x2=xp2*sqrt((2.-(1.-xx)*(1.+yy))/((2.-(1.-xx)*(1.-yy))*xx));}
      if (!MCNLO) {
      s=x1*x2*S2;
      t=-s/2.*(1.-xx)*(1.+yy);
      u=-s/2.*(1.-xx)*(1.-yy);
} 
     
     
      double xge = -(t/(x1*S2)+u/(x2*S2));
      double xgz = t/(x1*S2)-u/(x2*S2);
      double xze = x1+x2-xge;
      double xzz = x1-x2-xgz;
       if (emit) {kt = sqrt(pow(xge,2)-pow(xgz,2))*emcm/2.;}
        SCALUP = kt;
      //4-MOMENTA

      PUP[3][1]=x1*emcm/2.;
      PUP[4][1]=PUP[3][1];
      PUP[3][2]=-x2*emcm/2.;
      PUP[4][2]=-PUP[3][2];
      double phi=2.*pi*random(seed); //uniform azimuth
      PUP[1][3]=kt*cos(phi);
      PUP[2][3]=kt*sin(phi);
      PUP[4][3]=xze*emcm/2.;
      PUP[3][3]=xzz*emcm/2.;
      PUP[1][4]=-PUP[1][3];
      PUP[2][4]=-PUP[2][3];
      PUP[4][4]=xge*emcm/2.;
      PUP[3][4]=xgz*emcm/2.;
       if (emit) {
      double ZZ=sqrt(pow(PUP[1][3],2)+pow(PUP[2][3],2)+pow(PUP[3][3],2));
    double cost=PUP[3][3]/sqrt(pow(PUP[1][3],2)+pow(PUP[2][3],2)+pow(PUP[3][3],2));
    double sint=kt/sqrt(pow(PUP[1][3],2)+pow(PUP[2][3],2)+pow(PUP[3][3],2));
    PUP[1][3]=pt1*cos(phi);
    PUP[2][3]=pt1*sin(phi);
    PUP[1][7]=-pt1*cos(phi);
    PUP[2][7]=-pt1*sin(phi);
    PUP[3][3]=(pow(Q,2)*ZZ+ZZ*pow(PUP[4][3],2)-pow(ZZ,3)+sqrt(2.*pow(Q*ZZ*PUP[4][3],2)-2.*pow(ZZ,2)*pow(PUP[4][3],4)+pow(ZZ,4)*pow(PUP[4][3],2)-4.*pow(PUP[4][3],4)*pow(pt1,2)+pow(PUP[4][3],2)*pow(Q,4)-2.*pow(PUP[4][3],4)*pow(Q,2)+pow(PUP[4][3],6)+4.*pow(ZZ*pt1*PUP[4][3],2)))/(2.*(pow(PUP[4][3],2)-pow(ZZ,2)));
      if ( random(seed) < 0.5){
	PUP[3][3]=(pow(Q,2)*ZZ+ZZ*pow(PUP[4][3],2)-pow(ZZ,3)-sqrt(2.*pow(Q*ZZ*PUP[4][3],2)-2.*pow(ZZ,2)*pow(PUP[4][3],4)+pow(ZZ,4)*pow(PUP[4][3],2)-4.*pow(PUP[4][3],4)*pow(pt1,2)+pow(PUP[4][3],2)*pow(Q,4)-2.*pow(PUP[4][3],4)*pow(Q,2)+pow(PUP[4][3],6)+4.*pow(ZZ*pt1*PUP[4][3],2)))/(2.*(pow(PUP[4][3],2)-pow(ZZ,2)));}
      PUP[3][7]=ZZ-PUP[3][3];
      PUP[1][3]=(pt1*cost+PUP[3][3]*sint)*cos(phi);
      PUP[2][3]=(pt1*cost+PUP[3][3]*sint)*sin(phi);
      PUP[3][3]=-pt1*sint+PUP[3][3]*cost;
      PUP[4][3]=sqrt(pow(PUP[1][3],2)+pow(PUP[2][3],2)+pow(PUP[3][3],2)+pow(Q,2));
      PUP[1][7]=(-pt1*cost+PUP[3][7]*sint)*cos(phi);
      PUP[2][7]=(-pt1*cost+PUP[3][7]*sint)*sin(phi);
      PUP[3][7]=pt1*sint+PUP[3][7]*cost;
      PUP[4][7]=sqrt(pow(PUP[1][7],2)+pow(PUP[2][7],2)+pow(PUP[3][7],2));}

      double cth, cofac, maxcofac;
	  if (proc==23) {do 
          {cth=2.*random(seed)-1.;
	  double GF = 0.0000116639;
	  double alphaem = 0.007297352;
	  double e2 = 4.*pi*alphaem;
	  double k = sqrt(2.)*GF*pow(mp,2)/e2;
	  double sin2thw = 0.2312;
          double chi1 = k*pow(Q,2)*(pow(Q,2)-pow(mp,2))/(pow((pow(emcm,2)-pow(mp,2)),2)+pow(Yzz*mp,2));
          double chi2 = pow(k,2)*pow(Q,4)/(pow((pow(Q,2)-pow(mp,2)),2)+pow(Yzz*mp,2));
	  double Aq,Vq,Ae,Ve,Qq;
	  Ae = -0.5;
	  Ve = Ae+2.*sin2thw;
	  if (fabs(ID1)==1||fabs(ID1)==3||fabs(ID1)==5) {Qq=-1./3.; Aq=-0.5; Vq=Aq-(2.*Qq*sin2thw);}
	  if (fabs(ID1==2)||fabs(ID1)==4) {Qq=2./3.; Aq=0.5; Vq=Aq-(2.*Qq*sin2thw);}
	  cofac = (1.+pow(cth,2))*(1.+2.*Vq*Ve*chi1+(pow(Aq,2)+pow(Vq,2))*(pow(Ae,2)+pow(Ve,2))*chi2)+cth*(4.*Qq*Aq*Ae*chi1+8.*Aq*Vq*Ae*Ve*chi2);
          maxcofac = 2.*(1.+2.*Vq*Ve*chi1+(pow(Aq,2)+pow(Vq,2))*(pow(Ae,2)+pow(Ve,2))*chi2)+(4.*Qq*Aq*Ae*chi1+8.*Aq*Vq*Ae*Ve*chi2);}while(random(seed) > cofac/maxcofac);}
          if (proc==24) {do
	  {cth=2.*random(seed)-1.;
	  cofac=pow(1.+cth,2);
	  maxcofac=4.;} while(random(seed) > cofac/maxcofac);} 
 	  if (proc==24) {
          if (ID[0]==-1||ID[0]==-3||ID[0]==-5 ||ID[0]==2||ID[0]==4){
	    cth=-cth;} }
	  //	  cout << "cth" << cth << endl;
	  double phi2=2.*pi*random(seed);
	  //Assign 4-momenta in boson rest frame
	  PUP[3][5]=Q*cth/2.;
	  PUP[1][5]=sqrt(pow(Q,2)/4.-pow(PUP[3][5],2))*cos(phi2);
	  PUP[2][5]=sqrt(pow(Q,2)/4.-pow(PUP[3][5],2))*sin(phi2);
	  PUP[4][5]=Q/2.;
	  PUP[1][6]=-sqrt(pow(Q,2)/4.-pow(PUP[3][5],2))*cos(phi2);
	  PUP[2][6]=-sqrt(pow(Q,2)/4.-pow(PUP[3][5],2))*sin(phi2);
	  PUP[3][6]=-PUP[3][5];
	  PUP[4][6]=Q/2.;  

      mass[3]=Q;
      // Boost to Lab frame
      double v=sqrt(pow(PUP[4][3],2)-pow(Q,2))/PUP[4][3];
      PUP[4][5]=(Q/2.-v*PUP[3][5])/sqrt(1.-pow(v,2));
      PUP[3][5]=(PUP[3][5]-Q*v/2.)/sqrt(1.-pow(v,2));
      PUP[4][6]=(Q/2.-v*PUP[3][6])/sqrt(1.-pow(v,2));
      PUP[3][6]=(PUP[3][6]-Q*v/2.)/sqrt(1.-pow(v,2));
      double ktt1=-sqrt(pow(PUP[4][5],2)-pow(PUP[3][5],2));
      double ktt2=sqrt(pow(PUP[4][6],2)-pow(PUP[3][6],2));
      if (emit) {kt = sqrt(pow(PUP[1][3],2)+pow(PUP[2][3],2));}
      double tht=atan2(kt,PUP[3][3]);
      double pe1=PUP[4][5];
      double pe2=PUP[4][6];
      double pz1=PUP[3][5];
      double pz2=PUP[3][6];
      PUP[4][5]=pe2;
      PUP[4][6]=pe1;
      PUP[3][5]=-(pz2*cos(tht)+ktt2*sin(tht));
      PUP[3][6]=-(pz1*cos(tht)+ktt1*sin(tht));
      PUP[1][5]=-(pz2*sin(tht)-ktt2*cos(tht))*cos(phi);
      PUP[2][5]=-(pz2*sin(tht)-ktt2*cos(tht))*sin(phi);
      PUP[1][6]=-(pz1*sin(tht)-ktt1*cos(tht))*cos(phi);
      PUP[2][6]=-(pz1*sin(tht)-ktt1*cos(tht))*sin(phi);
     
      //selecting between 2 solutions
        if (random(seed) < 0.5) {
	  PUP[4][5]=pe1;
	  PUP[4][6]=pe2;
	  PUP[3][5]=-(pz1*cos(tht)+ktt1*sin(tht));
	  PUP[3][6]=-(pz2*cos(tht)+ktt2*sin(tht));
	  PUP[1][5]=-(pz1*sin(tht)-ktt1*cos(tht))*cos(phi);
	  PUP[2][5]=-(pz1*sin(tht)-ktt1*cos(tht))*sin(phi);
	  PUP[1][6]=-(pz2*sin(tht)-ktt2*cos(tht))*cos(phi);
	  PUP[2][6]=-(pz2*sin(tht)-ktt2*cos(tht))*sin(phi);
	}
 	if (!MCNLO){ 
        if (pro == 2) {Jt = 5; Jtt = 1;}
        if (pro == 3) {Jt = 5; Jtt = 2;}}
        //Assign particle IDs
	if (Jt==1||Jt==3) {
	  IDUP[1]=ID[0];
	  IDUP[2]=ID[1];
	  IDUP[4]=21;
	  if (emit) {IDUP[4] = 21; IDUP[7] = 21;}
}
        if (Jt==5||Jt==7) {
	  if (Jtt==1) {
	    IDUP[1]=ID[0];
	    IDUP[2]=21;
	    IDUP[4]=IDUP[1];
	    if (emit) {IDUP[4] = 21 ; IDUP[7] = IDUP[1];}
}
	  else {
	    IDUP[1]=21;
	    IDUP[2]=ID[1];
	    IDUP[4]=IDUP[2];
	    if (emit) {IDUP[4] = 21 ; IDUP[7] = IDUP[2];}
}}

	if (proc==23) {
	  IDUP[3]=23;
	  IDUP[5]=11;
	  IDUP[6]=-11;
	  } else if (proc==24) {
	 if (ID[0]==1||ID[0]==3||ID[0]==5 ||ID[0]==-2||ID[0]==-4){
	   IDUP[3]=-24;
          IDUP[5]=11;
	  IDUP[6]=-12;} else {
          IDUP[3]=24;
          IDUP[5]=-11;
	  IDUP[6]=12;}}

	  ISTUP[4]=1;
	  ISTUP[5]=1;
	  ISTUP[6]=1;
	  ICOLUP[1][5]=0;
	  ICOLUP[2][5]=0;
	  ICOLUP[1][6]=0;
	  ICOLUP[2][6]=0;
	  MOTHUP[1][5]=3;
	  MOTHUP[2][5]=3;
	  MOTHUP[1][6]=3;
	  MOTHUP[2][6]=3;
	  ISPINUP[1]=-1;
	  ISPINUP[2]=-1;
	  ISPINUP[4]=-1;
	  ISPINUP[3]=+1;
	  ISPINUP[5]=+1;
	  ISPINUP[6]=-1;
	  NUP=6; 
	  if (emit) {NUP=7;
	  PUP[4][8] = PUP[4][4];
          PUP[3][8] = PUP[3][4];
          PUP[2][8] = PUP[2][4];
          PUP[1][8] = PUP[1][4];
          PUP[4][4] = PUP[4][7];
          PUP[3][4] = PUP[3][7];
          PUP[2][4] = PUP[2][7];
          PUP[1][4] = PUP[1][7];
          PUP[4][7] = PUP[4][8];
          PUP[3][7] = PUP[3][8];
          PUP[2][7] = PUP[2][8];
          PUP[1][7] = PUP[1][8];
	  ISPINUP[7]=-1;
	  ISTUP[7]=1;}
      // Assign mother information and colour using Kleiss trick.
	  if (Jt==5||Jt==7) {
	    if (IDUP[1]==21) {
	      MOTHUP[1][4]=1;
	      MOTHUP[2][4]=1;
	      if (emit) {
              MOTHUP[1][7]=2;
	      MOTHUP[2][7]=2;
	  }
	      if (IDUP[2]> 0.) {
		ICOLUP[1][2]=501;
		ICOLUP[2][2]=0;
		ICOLUP[1][1]=502;
		ICOLUP[2][1]=501;
		ICOLUP[1][4]=502;
		ICOLUP[2][4]=0;
		  if (emit) {
	    ICOLUP[1][2]=503;
	    ICOLUP[2][2]=0;
	    ICOLUP[1][1]=502;
	    ICOLUP[2][1]=501;
	    ICOLUP[1][7]=502;
	    ICOLUP[2][7]=0;
	    ICOLUP[1][4]=503;
	    ICOLUP[2][4]=501;
	  }
	      }        else {
		ICOLUP[1][2]=0;
		ICOLUP[2][2]=501;
		ICOLUP[1][1]=501;
		ICOLUP[2][1]=502;
		ICOLUP[1][4]=0;
		ICOLUP[2][4]=502;
	if (emit) {
	    ICOLUP[1][2]=0;
	    ICOLUP[2][2]=503;
	    ICOLUP[1][1]=501;
	    ICOLUP[2][1]=502;
	    ICOLUP[1][7]=0;
	    ICOLUP[2][7]=502;
	    ICOLUP[1][4]=501;
	    ICOLUP[2][4]=503;  
	}
	      }
	    }
	    if (IDUP[2]==21) {
	      MOTHUP[1][4]=2;
	      MOTHUP[2][4]=2;
	      if (IDUP[1] > 0) {
		ICOLUP[1][1]=501;
		ICOLUP[2][1]=0;
		ICOLUP[1][2]=502;
		ICOLUP[2][2]=501;
		ICOLUP[1][4]=502;
		ICOLUP[2][4]=0;
		if (emit) {
	    ICOLUP[1][1]=503;
	    ICOLUP[2][1]=0;
	    ICOLUP[1][2]=502;
	    ICOLUP[2][2]=501;
	    ICOLUP[1][7]=502;
	    ICOLUP[2][7]=0;
	    ICOLUP[1][4]=503;
	    ICOLUP[2][4]=501;  
	  }
}
	      else {
		ICOLUP[1][1]=0;
		ICOLUP[2][1]=501;
		ICOLUP[1][2]=501;
		ICOLUP[2][2]=502;
		ICOLUP[1][4]=0;
		ICOLUP[2][4]=502;
		if (emit) {
	    ICOLUP[1][1]=0;
	    ICOLUP[2][1]=503;
	    ICOLUP[1][2]=501;
	    ICOLUP[2][2]=502;
	    ICOLUP[1][7]=0;
	    ICOLUP[2][7]=502;
	    ICOLUP[1][4]=501;
	    ICOLUP[2][4]=503;  
	    }
	      }
	    } }

       if (Jt==1||Jt==3) {
       if (random(seed) < pow(PUP[4][1],2)/(pow(PUP[4][1],2)
           +pow(PUP[4][2],2))){
	  MOTHUP[1][4]=2;
	  MOTHUP[2][4]=2;
	  if (emit){
	  MOTHUP[1][7]=2;
	  MOTHUP[2][7]=2;  
	  }
          if (IDUP[1] > 0.) {
	  ICOLUP[1][1]=501;
	  ICOLUP[2][1]=0;
	  ICOLUP[1][2]=0;
	  ICOLUP[2][2]=502;
	  ICOLUP[1][4]=501;
	  ICOLUP[2][4]=502;
	  if (emit) {
	  ICOLUP[1][1]=501;
	  ICOLUP[2][1]=0;
	  ICOLUP[1][2]=0;
	  ICOLUP[2][2]=502;
	  ICOLUP[1][7]=503;
	  ICOLUP[2][7]=502;
	  ICOLUP[1][4]=501;
	  ICOLUP[2][4]=503;
	  }
}
	else{
	  ICOLUP[1][1]=0;
	  ICOLUP[2][1]=501;
	  ICOLUP[1][2]=502;
	  ICOLUP[2][2]=0;
	  ICOLUP[1][4]=502;
	  ICOLUP[2][4]=501;
	  if (emit) {
	  ICOLUP[1][1]=0;
	  ICOLUP[2][1]=501;
	  ICOLUP[1][2]=502;
	  ICOLUP[2][2]=0;
	  ICOLUP[1][7]=502;
          ICOLUP[2][7]=503;
	  ICOLUP[1][4]=503;
	  ICOLUP[2][4]=501;}
}	  
	}
	else {
	  MOTHUP[1][4]=1;
	  MOTHUP[2][4]=1;
	  if (emit) {
	  MOTHUP[1][7]=1;
	  MOTHUP[2][7]=1;}
	  	if (IDUP[1] > 0.) {
	  ICOLUP[1][1]=501;
	  ICOLUP[2][1]=0;
	  ICOLUP[1][2]=0;
	  ICOLUP[2][2]=502;
	  ICOLUP[1][4]=501;
	  ICOLUP[2][4]=502;
	  if (emit) {
	  ICOLUP[1][1]=503;
	  ICOLUP[2][1]=0;
	  ICOLUP[1][2]=0;
	  ICOLUP[2][2]=502;
	  ICOLUP[1][7]=503;
	  ICOLUP[2][7]=501;
	  ICOLUP[1][4]=501;
	  ICOLUP[2][4]=502;
	  }
}
	else{
	  ICOLUP[1][1]=0;
	  ICOLUP[2][1]=501;
	  ICOLUP[1][2]=502;
	  ICOLUP[2][2]=0;
	  ICOLUP[1][4]=502;
	  ICOLUP[2][4]=501;
	  if (emit) {
	  ICOLUP[1][1]=0;
	  ICOLUP[2][1]=503;
	  ICOLUP[1][2]=502;
	  ICOLUP[2][2]=0;
	  ICOLUP[1][7]=501;
          ICOLUP[2][7]=503;
	  ICOLUP[1][4]=502;
	  ICOLUP[2][4]=501;}
	}
	}
        }}else {

        //Generate Born events
	  
	  SCALUP=Q;
	  ISPINUP[1]=-1;
	  ISPINUP[2]=-1;
	  ISPINUP[4]=-1;
	  ISPINUP[3]=+1;
	  ISPINUP[5]=+1;
          PUP[3][1]=xq*emcm/2.;
	  PUP[4][1]=xq*emcm/2.;
	  PUP[3][2]=-xqb*emcm/2.;
	  PUP[4][2]=xqb*emcm/2.;
	  PUP[1][3]=0.;
	  PUP[2][3]=0.;
	  PUP[3][3]=PUP[3][1]+PUP[3][2];
	  PUP[4][3]=PUP[4][1]+PUP[4][2];
	  double cth, cofac, maxcofac;
	  if (proc==23) {do 
          {cth=2.*random(seed)-1.;
	  double GF = 0.0000116639;
	  double alphaem = 0.007297352;
	  double e2 = 4.*pi*alphaem;
	  double k = sqrt(2.)*GF*pow(mp,2)/e2;
	  double sin2thw = 0.2312;
          double chi1 = k*pow(Q,2)*(pow(Q,2)-pow(mp,2))/(pow((pow(emcm,2)-pow(mp,2)),2)+pow(Yzz*mp,2));
          double chi2 = pow(k,2)*pow(Q,4)/(pow((pow(Q,2)-pow(mp,2)),2)+pow(Yzz*mp,2));
	  double Aq,Vq,Ae,Ve,Qq;
	  Ae = -0.5;
	  Ve = Ae+2.*sin2thw;
	  if (ID[0]==1||ID[0]==3||ID[0]==5) {Qq=-1./3.; Aq=-0.5; Vq=Aq-(2.*Qq*sin2thw);}
	  if (ID[0]==2||ID[0]==4) {Qq=2./3.; Aq=0.5; Vq=Aq-(2.*Qq*sin2thw);}
	  cofac = (1.+pow(cth,2))*(1.+2.*Vq*Ve*chi1+(pow(Aq,2)+pow(Vq,2))*(pow(Ae,2)+pow(Ve,2))*chi2)+cth*(4.*Qq*Aq*Ae*chi1+8.*Aq*Vq*Ae*Ve*chi2);
          maxcofac = 2.*(1.+2.*Vq*Ve*chi1+(pow(Aq,2)+pow(Vq,2))*(pow(Ae,2)+pow(Ve,2))*chi2)+(4.*Qq*Aq*Ae*chi1+8.*Aq*Vq*Ae*Ve*chi2);}while(random(seed) > cofac/maxcofac);}
          if (proc==24) {do
	  {cth=2.*random(seed)-1.;
	  cofac=pow(1.+cth,2);
	  maxcofac=4.;} while(random(seed) > cofac/maxcofac);}
	  if (proc==24) {
          if (ID[0]==-1||ID[0]==-3||ID[0]==-5 ||ID[0]==2||ID[0]==4){
	    cth=-cth;} }
	  double phi=2.*pi*random(seed);
	  //Assign 4-momenta in boson rest frame
	  PUP[3][4]=Q*cth/2.;
	  PUP[1][4]=sqrt(pow(Q,2)/4.-pow(PUP[3][4],2))*cos(phi);
	  PUP[2][4]=sqrt(pow(Q,2)/4.-pow(PUP[3][4],2))*sin(phi);
	  PUP[4][4]=Q/2.;
	  PUP[1][5]=-sqrt(pow(Q,2)/4.-pow(PUP[3][4],2))*cos(phi);
	  PUP[2][5]=-sqrt(pow(Q,2)/4.-pow(PUP[3][4],2))*sin(phi);
	  PUP[3][5]=-PUP[3][4];
	  PUP[4][5]=Q/2.;
	  mass[3]=Q;
	  //Boost to lab frame
	 
	  double v = (xq - xqb) / (xq +xqb); 
	  PUP[3][4]=-(PUP[3][4]-PUP[4][4]*v)/sqrt(1.-pow(v,2));
	  PUP[4][4]=(PUP[4][4]-Q*cth*v/2.)/sqrt(1.-pow(v,2));
	  PUP[3][5]=-(PUP[3][5]-PUP[4][5]*v)/sqrt(1.-pow(v,2));
	  PUP[4][5]=(PUP[4][5]+Q*cth*v/2.)/sqrt(1.-pow(v,2));
	  if (fabs(PUP[3][1]+PUP[3][2]-PUP[3][4]-PUP[3][5]) > 0.001) {
	    PUP[3][4]=-PUP[3][4];
	    PUP[3][5]=-PUP[3][5];}
	  //Assign particle IDs
	  IDUP[1]=ID[0];
	  IDUP[2]=ID[1];
	  if (proc==23) {
	    IDUP[3]=23;
	  IDUP[4]=11;
	  IDUP[5]=-11;} 
         else if (proc==24) {
	 if (ID[0]==1||ID[0]==3||ID[0]==5 ||ID[0]==-2||ID[0]==-4){
	   IDUP[3]=-24;
          IDUP[4]=11;
	  IDUP[5]=-12;} else {
	    IDUP[3]=24;
          IDUP[4]=-11;
	  IDUP[5]=12;}}

	  ISTUP[4]=1;
	  ISTUP[5]=1;
	  NUP=5;
	  MOTHUP[1][4]=3;
	  MOTHUP[2][4]=3;
	  MOTHUP[1][5]=3;
	  MOTHUP[2][5]=3;
	  ICOLUP[1][4]=0;
	  ICOLUP[2][4]=0;
	  ICOLUP[1][5]=0;
	  ICOLUP[2][5]=0;
	  //Assing colours
	  if (IDUP[1]> 0.) {
	    ICOLUP[1][1]=501;
	    ICOLUP[2][1]=0;
	    ICOLUP[1][2]=0;
	    ICOLUP[2][2]=501;}
	  else {
	    ICOLUP[1][1]=0;
	    ICOLUP[2][1]=501;
	    ICOLUP[1][2]=501;
	    ICOLUP[2][2]=0;}
		}
	 ii+=1;
	
	 double sum1, sum2, sum3, sum4;
	 sum1=0.;sum2=0.;sum3=0.;sum4=0.;
	 

	 sumwgt += XWGTUP;
	 outdata << "<event>" << endl;
     outdata << NUP <<"\t" << proc <<"\t"<<XWGTUP<<"\t" << SCALUP << "\t" << AQEDUP << "\t" << AQCDUP << endl;
     for (int ja = 1; ja < NUP+1; ja++) {
       outdata << IDUP[ja] << "\t" << ISTUP[ja] << "\t" << MOTHUP[1][ja] <<"\t" << MOTHUP[2][ja] << "\t" << ICOLUP[1][ja] << "\t" << ICOLUP[2][ja] <<"\t" <<setprecision (9)<< PUP[1][ja] << "\t" << PUP[2][ja]<<"\t"<<PUP[3][ja]<<"\t"<<PUP[4][ja]<<"\t" << mass[ja] <<"\t"<< "0" << "\t" <<"9" <<endl;
      }
 outdata << "</event>" << endl;
  cout  << "Generated event : " << "\t" << ii << "\t" << "of" << "\t" << nevg << "\r" << flush;  
}
  
  outdata << "</LesHouchesEvents>" << endl;
    cout << endl;
  return 0;
   }
 /*********************************************************************/
   double G(double &rp2, double &rNq, double &rbb, double &rqq, double &rpm, double &rll) {
     double rG=rNq/(2.*rbb)*(log(rqq/pow(rll,2))*log(log(rpm/pow(rll,2))/log(rp2/pow(rll,2)))-log(rpm/rp2));
     return rG;}
   double DG(double &rp2, double &rNq,double &rbb,double &rqq,double &rpm,double &rll) {
     double rDG=rNq/(2.*rbb)*(-log(rqq/pow(rll,2))*(1./(log(rp2/pow(rll,2))*rp2))+1./rp2);
   return rDG;}
 /*********************************************************************/
/* rx and ry are the radiation variables. rx is given by Q^2/S where Q is the 
 boson mass and S is the centre of mass energy squared. ry is the cosine of 
 the angle of emission of the parton at NLO in the centre of mass frame.*/

/**************************************************************************/

/* Phase space region DYBC includes the soft and collinear regions for
   gluon emission. Returns the integrand.*/
 
       double DYBC(double &rx,double &ry, double r3, double r4) { 
       double wgt = 0; 
       bool soft = false;
       double xr=rx;
       double yr=ry;
// Change of variables to regulate divergences for soft = true

    if (xr < 0.125 && yr > 0.75) {
      soft = true;
      xr=0.125*r3;
      yr=1.-2.*xr*r4;

       } else {
      if (xr < 0.125 && yr < -0.75) {
        soft = true;
	 xr=0.125*r3;
       yr=-(1.-2.*xr*r4);

      }
    }

      /* Check if phase space point is within the region defined by k < 1 , 
       calculate integrand (wgt) and multiply by Jacobian if soft = true*/

      double k=(1.-yr)/(xr*(1.+yr));
     
      if (k < 1.) {
	wgt=-(pow(yr,2)-2.*xr*pow(yr,2)+pow(xr*yr,2)+2.*xr*yr-2.*yr+8.*xr
            +3.*pow(xr,2)+5.)/(4.*xr*(1.+yr));
	if (soft) {
	 wgt=wgt*16.*xr ;
          }
      rx=xr; ry=yr;} else {
	k=(1.+yr)/(xr*(1.-yr));
	if (k < 1.) {
	 wgt=-(pow(yr,2)-2.*xr*pow(yr,2)+pow(xr*yr,2)-2.*xr*yr+2.*yr+8.*xr
             +3.*pow(xr,2)+5.)/(4.*xr*(1.-yr));
	 if (soft) {
	 wgt=wgt*16.*xr ;
	 }
	rx=xr; ry=yr;}else{
	  wgt= 0.;
	 	}
      }
     
      return wgt;
       }

 /**************************************************************************/

/* Phase space region DYD includes the hard regions for gluon emission.
   fac=Jacobian from change of variables. Returns the integrand */
 
     double DYD(double &rx, double &ry, double r3, double r4) {
       double xr=rx;
       double yr=ry;
       double wgt = 0; 
       bool soft=false;
       double fac = 0;

 // Change of variables to regulate divergences for soft = true

     if (xr < 0.25 && yr > 0.6) {
      soft = true;
      yr=1.-0.4*r3;
      xr=0.625*r4*(1.-yr);
      fac=5.*(1.-yr);
    }else{
      if (xr > 0.6 && yr < 0.25 && yr > 0.) { 
	soft = true;
	xr=(1.-0.4*r3);
	yr=0.625*r4*(1.-xr);
        fac=5.*(1.-xr);
      }else{
	if (xr < 0.25 && yr < -0.6) {
	  soft = true;
	  yr=-(1.-0.4*r3);
	  xr=0.625*r4*(1.+yr);
	  fac=5.*(1.+yr);
	}else{
	  if (xr > 0.6 && yr > -0.25 && yr < 0.) {
	    soft = true;
	    xr=1.-0.4*r3;
	    yr=-0.625*r4*(1.-xr);
	    fac=5.*(1.-xr);
    }
    }
    }
    }
   
     /* Check if phase space point is within the region defined by k > 1 ,
	calculate integrand (wgt) and multiply by Jacobian if soft = true*/

    double k=(1.-yr)/(xr*(1.+yr));
      if (k > 1.) {
	k=(1.+yr)/(xr*(1.-yr));
	wgt=0.;
	if (k > 1.) {
       wgt=(pow(yr,2)*pow((1.-xr),2)+pow((1.+xr),2))
           /((1.-xr)*(1.-pow(yr,2)));
	  
	  if (soft) {
	  wgt=wgt*fac;
          }
      rx=xr; ry=yr;} else {
	 wgt=0.;
	  
	}
      }
      return wgt;
  }

 /**********************************************************************/

/* Phase space region DYGBC includes the soft and collinear regions for 
   quark emission. Returns the integrand */

      double DYGBC(double &rx, double &ry, double r3, double r4) { 
      double xr=rx;
      double yr=ry;
      double wgt = 0; 
      bool soft = false;

 // Change of variables to regulate divergences for soft = true

     if(xr < 0.125 && yr > 0.75) {
       soft = true;
       xr = 0.125*r3;
       yr = 1.-2.*xr*r4;
      }

 /* Check if phase space point is within the region defined by k < 1 , 
    calculate integrand (wgt) and multiply by Jacobian if soft = true */
    
    double k=(1.-yr)/(xr*(1.+yr));
    if (k < 1.) {
      wgt=(-3.*pow(xr,2)*pow(yr,2)+3.*xr*pow(yr,2)-pow(yr,2)
         +pow(xr,3)*pow(yr,2)+3.*pow(xr,3)*yr-6.*pow(xr,2)*yr
         +3.*xr*yr-1.+4.*xr-7.*pow(xr,2)+4.*pow(xr,3))/(4.*xr);  
      if (soft){
    wgt=wgt*16.*xr ;
      }
    rx=xr; ry=yr;}else{
      wgt = 0.;
    }
    return wgt;
  }

 /*********************************************************************/

/* Phase space region DYGD includes the hard region for quark emission.
   fac=Jacobian from change of variables. Returns the integrand */

    double DYGD(double &rx, double &ry, double r3, double r4) {
     double xr=rx;
     double yr=ry;
     double fac = 0;
     double wgt = 0; 
     bool soft = false;
    if(xr < 0.25 && yr > 0.6) {
      soft = true;
      yr = 1.-0.4*r3;
      xr = 0.625*r4*(1.-yr);
      fac = 5.*(1.-yr);
    }

    /* Check if phase space point is within the region defined by k > 1 ,
    calculate integrand (wgt) and multiply by Jacobian if soft = true */
   
     double k=(1.-yr)/(xr*(1.+yr));
    if (k > 1.) {
      wgt=((3.+pow(yr,2))*(1.-xr)*(1.-xr)-2.*yr*(1.-pow(xr,2))
         +2.*(1.+pow(xr,2)))/(4.*(1.-yr));
      if (soft) {
	 wgt = fac*wgt;
	  }
    rx=xr; ry=yr;}else { wgt=0.;}
    return wgt;
  }

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
