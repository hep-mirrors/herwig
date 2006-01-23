#include "MRST.h"
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <istream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace ThePEG;
using namespace Herwig;

MRST::MRST() {}

MRST::MRST(const MRST &x) : PDFBase(x) {
   for(int i=0; i<=np; i++) for(int j=0; j<=nx; j++) for(int k=0; k<=nq; k++)   
      data[i][j][k] = x.data[i][j][k];
   initialize(false);
}//, dataPtr(x.dataPtr) {}

MRST::~MRST() {}

ClassDescription<MRST> MRST::initMRST;

bool MRST::canHandleParticle(tcPDPtr particle) const {
  // Return true if this PDF can handle the extraction of parton from the
  // given particle ie. if the particle is a proton or neutron.
  return ( abs(particle->id()) == abs(long(ParticleID::pplus)) ||
           abs(particle->id()) == abs(long(ParticleID::n0)) );
}

cPDVector MRST::partons(tcPDPtr p) const {
  // Return the parton types which are described by these parton
  // densities.
  cPDVector ret;
  if ( canHandleParticle(p) ) {
    ret.push_back(getParticleData(ParticleID::g));
    for ( int i = 1; i <= 5; ++i ) {
      ret.push_back(getParticleData(i));
      ret.push_back(getParticleData(-i));
    }
  }
  return ret;
}

double MRST::xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		 double l, Energy2 particleScale) const {
  //cout << "Calling MRST::xfl = ";
  double x = exp(-l);
  double a = xfx(particle,parton,partonScale,x,0.0,particleScale);
  //cout << a << endl;
  return a;
}

double MRST::xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                 double x, double eps, Energy2 particleScale) const {
  // Return the value of the density of parton at the given a scale
  // and fractional momentum l (the optional virtuality of the
  // incoming particle is not used).
  // Must de-constify this. Don't know why this function is const, seems 
  // logical to change internal structure...
  update(x, partonScale/GeV2);
  //if ( S() < 0.0 ) return 0.0;
  bool anti = particle->id() < 0;
  bool neutron = abs(particle->id()) == ParticleID::n0;
  switch(parton->id()) {
  case ParticleID::b:
  case ParticleID::bbar:
    return table[bot];
  case ParticleID::c:
  case ParticleID::cbar:
    return table[chm];
  case ParticleID::s:
  case ParticleID::sbar:
    return table[str];
  case ParticleID::u:
    return (neutron? (table[dnSea] + (anti? 0.0: table[dnValence])) :
	    (table[upSea] + (anti? 0.0: table[upValence])));
  case ParticleID::ubar:
    return (neutron? (table[dnSea] + (anti? table[dnValence]: 0.0)) :
	       (table[upSea] + (anti? table[upValence]: 0.0)));
  case ParticleID::d:
    return (neutron? (table[upSea] + (anti? 0.0: table[upValence])) :
	       (table[dnSea] + (anti? 0.0: table[dnValence])));
  case ParticleID::dbar:
    return (neutron? (table[upSea] + (anti? table[upValence]: 0.0)) :
	       (table[dnSea] + (anti? table[dnValence]: 0.0)));
  case ParticleID::g:
    return table[glu];
  }
  return 0.0;
}

double MRST::xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		  double l, Energy2 particleScale) const {
  double x = exp(-l);
  return xfvx(particle,parton,partonScale,x,0.0,particleScale);
}

double MRST::xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                  double x, double eps, Energy2 particleScale) const {
  // Return the valaens part of the density of parton at the given a
  // scale and fractional momentum x (the optional virtuality of
  // the incoming particle is not used).
  update(x, partonScale/GeV2);
  //if ( S() < 0.0 ) return 0.0;
  bool anti = particle->id() < 0;
  bool neutron = abs(particle->id()) == ParticleID::n0;
  switch(parton->id()) {
  case ParticleID::u:
    return (neutron? (anti? 0.0: table[dnValence]): 
	    (anti? 0.0: table[upValence]));
  case ParticleID::ubar:
    return (neutron? (anti? table[dnValence]: 0.0): 
	    (anti? table[upValence]: 0.0));
  case ParticleID::d:
    return (neutron? (anti? 0.0: table[upValence]): 
	    (anti? 0.0: table[dnValence]));
  case ParticleID::dbar:
    return (neutron? (anti? table[upValence]: 0.0): 
	    (anti? table[dnValence]: 0.0));
  }
  return 0.0;
}

int MRST::locate(double xx[],int n,double x) const {
  // returns an integer j such that x lies inbetween xx[j] and xx[j+1]. unit
  // offset of increasing ordered array xx assumed. n is the length of the
  // array (xx[n] highest element)
  int ju,jm,jl(0),j;
  ju=n+1;

  while (ju-jl>1) {
    jm=(ju+jl)/2; // compute a mid point.
    if(x >= xx[jm]) jl=jm;
    else ju=jm;
  }
  if(x==xx[1]) j=1;
  else if(x==xx[n]) j=n-1;
  else j=jl;

  return j;
}

double MRST::polderivative(double x1, double x2, double x3,
                           double y1, double y2, double y3) {
  // returns the estimate of the derivative at x2 obtained by a polynomial 
  // interpolation using the three points (x_i,y_i)
  return (x3*x3*(y1-y2) - 2.0*x2*(x3*(y1-y2) + x1*(y2-y3)) +
          x2*x2*(y1-y3) + x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3));
}

void MRST::persistentOutput(PersistentOStream &out) const {
  //out << dataPtr;
  out << _file;
  for(int i=0; i<=np; i++) for(int j=0; j<=nx; j++) for(int k=0; k<=nq; k++)
     out << data[i][j][k];
}

void MRST::persistentInput(PersistentIStream & in, int version) {
  //in >> dataPtr;
  in >> _file;
  for(int i=0; i<=np; i++) for(int j=0; j<=nx; j++) for(int k=0; k<=nq; k++)
     in >> data[i][j][k];
  initialize(false);
}

void MRST::Init() {
  static ClassDocumentation<MRST> docMRST("MRST");

  /*static Reference<MRST,MRSTData> interfaceData
    ("Data",
     "This is the MRSTData derived class with the desired data loaded.",
     &MRST::dataPtr, false, false, true, false);
  */
}

void MRST::doupdate() throw(UpdateException) { PDFBase::doupdate(); }
void MRST::doinit() throw(InitException) { 
  PDFBase::doinit();
  //  initialize();
}
//void MRST::doinitrun() { PDFBase::doinitrun(); initialize(); }

void MRST::dofinish() { PDFBase::dofinish(); }

void MRST::rebind(const TranslationMap & trans) throw(RebindException) {
  PDFBase::rebind(trans);
}

IBPtr MRST::clone() const { return new_ptr(*this); }
IBPtr MRST::fullclone() const { return clone(); }

IVector MRST::getReferences() { 
  IVector rval = PDFBase::getReferences(); 
  //rval.push_back(dataPtr);
  return rval;
}

void MRST::readSetup(istream &is) throw(SetupException) {
  _file = dynamic_cast<istringstream*>(&is)->str();
  initialize();
}

void MRST::initialize(bool reread) {
  //  if(reread) cout << "Opening file " << _file << endl;
  int i,n,m,k,l,j; // counters
  double dx,dq,dtemp;
  int wt[][16] = {{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		  { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
		  {-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0},
		  { 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
		  { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
		  { 0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1},
		  { 0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1},
		  {-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		  { 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
		  { 9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2},
		  {-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2},
		  { 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		  { 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
		  {-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1},
		  { 4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1}};
  // Variables used for initialising c_ij array: 
  double f1[np+1][nx+1][nq+1];//derivative wrt.x
  double f2[np+1][nx+1][nq+1];//derivative wrt.qq
  double f12[np+1][nx+1][nq+1];//cross derivative

  double xxd,d1d2,cl[16],x[16],d1,d2,y[5],y1[5],y2[5],y12[5];

  // Set xx and qq to the logs of their values. Makes entering below easier
  if(xx[1] == 1E-5) {
    for(n=1; n<=nx; n++) xx[n] = log10(xx[n]);
    for(n=1; n<=nq; n++) qq[n] = log10(qq[n]);
  }

  if(reread) {
     ifstream datafile;
     datafile.open(_file.c_str());

     if(datafile.bad()) { 
       cerr << "Could not open file " << _file << "\n";
       return;
     }
     for(int nn=1; nn<nx; nn++) {
       for(int mm=1; mm<=nq; mm++) {
         datafile >> data[1][nn][mm];
         datafile >> data[2][nn][mm];
         datafile >> data[3][nn][mm];
         datafile >> data[4][nn][mm];
         datafile >> data[5][nn][mm];
         datafile >> data[7][nn][mm];
         datafile >> data[6][nn][mm];
         datafile >> data[8][nn][mm];
         if(datafile.eof()) {
	   cerr << "Error while reading " << _file << endl;
	   return;
         }
       }
     }
  
     datafile >> dtemp;
     if(!datafile.eof()) {
       cerr << "Error reading end of " << _file << endl;
       return;
     }
     datafile.close();

     //     cout << "File read!" << endl;
  }

  // Now calculate the derivatives used for bicubic interpolation
  for (i=1;i<=np;i++) {
    // Start by calculating the first x derivatives
    // along the first x value
    dx=xx[2]-xx[1];
    for (m=1;m<=nq;m++)
      f1[i][1][m]=(data[i][2][m]-data[i][1][m])/dx;
    // The along the rest (up to the last)
    for (k=2;k<nx;k++) {
      for (m=1;m<=nq;m++) {
	f1[i][k][m]=polderivative(xx[k-1],xx[k],xx[k+1],
				  data[i][k-1][m],
				  data[i][k][m],
				  data[i][k+1][m]);
      }
    }
    // Then for the last column
    dx=xx[nx]-xx[nx-1];
    for (m=1;m<=nq;m++)
      f1[i][nx][m]=(data[i][nx][m]-data[i][nx-1][m])/dx;
    

    if ((i!=5)&&(i!=7)) {
      // then calculate the qq derivatives
      // Along the first qq value
      dq=qq[2]-qq[1];
      for (k=1;k<=nx;k++)
	f2[i][k][1]=(data[i][k][2]-data[i][k][1])/dq;
      // The rest up to the last qq value
      for (m=2;m<nq;m++) {
	for (k=1;k<=nx;k++)
	  f2[i][k][m]=polderivative(qq[m-1],qq[m],qq[m+1],
				    data[i][k][m-1],
				    data[i][k][m],
				    data[i][k][m+1]);
      }
      // then for the last row
      dq=qq[nq]-qq[nq-1];
      for (k=1;k<=nx;k++)
	f2[i][k][nq]=(data[i][k][nq]-data[i][k][nq-1])/dq;
      
      // Now, calculate the cross derivatives.
      // Calculate these as x-derivatives of the y-derivatives
      // ?? Could be improved by taking the average between dxdy and dydx ??
      
      // Start by calculating the first x derivatives
      // along the first x value
      dx=xx[2]-xx[1];
      for (m=1;m<=nq;m++)
	f12[i][1][m]=(f2[i][2][m]-f2[i][1][m])/dx;
      // The along the rest (up to the last)
      for (k=2;k<nx;k++) {
	for (m=1;m<=nq;m++)
	  f12[i][k][m]=polderivative(xx[k-1],xx[k],xx[k+1],
				     f2[i][k-1][m],f2[i][k][m],f2[i][k+1][m]);
      }
      // Then for the last column
      dx=xx[nx]-xx[nx-1];
      for (m=1;m<=nq;m++)
	f12[i][nx][m]=(f2[i][nx][m]-f2[i][nx-1][m])/dx;
    }
    
    if (i==5) {
      // zero all elements below the charm threshold
      for (m=1;m<=nqc0;m++) 
	for (k=1;k<=nx;k++)
	  f2[i][k][m]=0.0; 

      // then calculate the qq derivatives
      // Along the first qq value above the threshold
      dq=qq[nqc0+1]-qq[nqc0];
      for (k=1;k<=nx;k++)
	f2[i][k][m]=(data[i][k][m+1]-data[i][k][m])/dq;

      // The rest up to the last qq value
      for (m=nqc0+1;m<nq;m++) {
	for (k=1;k<=nx;k++)
	  f2[i][k][m]=polderivative(qq[m-1],qq[m],qq[m+1],
				    data[i][k][m-1],
				    data[i][k][m],
				    data[i][k][m+1]);
      }
      // then for the last row
      dq=qq[nq]-qq[nq-1];
      for (k=1;k<=nx;k++)
	f2[i][k][nq]=(data[i][k][nq]-data[i][k][nq-1])/dq;
      
      // Now, calculate the cross derivatives.
      // Calculate these as x-derivatives of the y-derivatives
      // ?? Could be improved by taking the average between dxdy and dydx ??

      dx=xx[2]-xx[1];
      for (m=1;m<=nq;m++)
	f12[i][1][m]=(f2[i][2][m]-f2[i][1][m])/dx;
      // The along the rest (up to the last)
      for (k=2;k<nx;k++) {
	for (m=1;m<=nq;m++)
	  f12[i][k][m]=polderivative(xx[k-1],xx[k],xx[k+1],
				     f2[i][k-1][m],f2[i][k][m],f2[i][k+1][m]);
      }
      // Then for the last column
      dx=xx[nx]-xx[nx-1];
      for (m=1;m<=nq;m++)
	f12[i][nx][m]=(f2[i][nx][m]-f2[i][nx-1][m])/dx;
    }

    if (i==7) {
      // zero all elements below the bottom threshold
      for (m=1;m<=nqb0;m++) 
	for (k=1;k<=nx;k++)
	  f2[i][k][m]=0.0; 

      // then calculate the qq derivatives
      // Along the first qq value above the threshold
      dq=qq[nqb0+1]-qq[nqb0];
      for (k=1;k<=nx;k++)
	f2[i][k][m]=(data[i][k][m+1]-data[i][k][m])/dq;

      // The rest up to the last qq value
      for (m=nqb0+1;m<nq;m++) {
	for (k=1;k<=nx;k++)
	  f2[i][k][m]=polderivative(qq[m-1],qq[m],qq[m+1],
				    data[i][k][m-1],
				    data[i][k][m],
				    data[i][k][m+1]);
      }
      // then for the last row
      dq=qq[nq]-qq[nq-1];
      for (k=1;k<=nx;k++)
	f2[i][k][nq]=(data[i][k][nq]-data[i][k][nq-1])/dq;
      
      // Now, calculate the cross derivatives.
      // Calculate these as x-derivatives of the y-derivatives
      // ?? Could be improved by taking the average between dxdy and dydx ??

      dx=xx[2]-xx[1];
      for (m=1;m<=nq;m++)
	f12[i][1][m]=(f2[i][2][m]-f2[i][1][m])/dx;
      // The along the rest (up to the last)
      for (k=2;k<nx;k++) {
	for (m=1;m<=nq;m++)
	  f12[i][k][m]=polderivative(xx[k-1],xx[k],xx[k+1],
				     f2[i][k-1][m],f2[i][k][m],f2[i][k+1][m]);
      }
      // Then for the last column
      dx=xx[nx]-xx[nx-1];
      for (m=1;m<=nq;m++)
	f12[i][nx][m]=(f2[i][nx][m]-f2[i][nx-1][m])/dx;
    }
	
	
    // Now calculate the coefficients c_ij
    for(n=1;n<=nx-1;n++) {
      for(m=1;m<=nq-1;m++) {
	d1=xx[n+1]-xx[n];
	d2=qq[m+1]-qq[m];
	d1d2=d1*d2;
	
	// Iterate around the grid and store the values of f, f_x, f_y and f_xy
	y[1]=data[i][n][m];
	y[2]=data[i][n+1][m];
	y[3]=data[i][n+1][m+1];
	y[4]=data[i][n][m+1];
	
	y1[1]=f1[i][n][m];
	y1[2]=f1[i][n+1][m];
	y1[3]=f1[i][n+1][m+1];
	y1[4]=f1[i][n][m+1];
	
	y2[1]=f2[i][n][m];
	y2[2]=f2[i][n+1][m];
	y2[3]=f2[i][n+1][m+1];
	y2[4]=f2[i][n][m+1];
	
	y12[1]=f12[i][n][m];
	y12[2]=f12[i][n+1][m];
	y12[3]=f12[i][n+1][m+1];
	y12[4]=f12[i][n][m+1];
	
	for (k=1;k<=4;k++) {
	  x[k-1]=y[k];
	  x[k+3]=y1[k]*d1;
	  x[k+7]=y2[k]*d2;
	  x[k+11]=y12[k]*d1d2;
	}
	
	for (l=0;l<=15;l++) {
	  xxd=0.0;
	  for (k=0;k<=15;k++) xxd+= wt[l][k]*x[k];
	  cl[l]=xxd;
	}
	
	l=0;
	for (k=1;k<=4;k++) 
	  for (j=1;j<=4;j++) c[i][n][m][k][j]=cl[l++];
      } //m
    } //n
  } // i
}

void MRST::update(double x, double q2) const {
  // Updates the parton content 
  double qsq;
  double xxx;
  int i,n,m,l;
  double t,u;
  
  //qsq=q*q;
  if(x<xmin) {
    x=xmin;
    //generator()->log() << "x   VALUE IS OUT OF RANGE (TOO LOW)";
  } else if(x>xmax) {
    x=xmax;
    //generator()->log() << "x  VALUE IS OUT OF RANGE (TOO HIGH)";  
  }
  
  if(q2<qsqmin) { 
    q2=qsqmin;
    //generator()->log() << "Q^2 VALUE IS OUT OF RANGE (TOO LOW) = " 
    //	       << q2 << endl;
  } else if(q2>qsqmax) {
    q2=qsqmax;
    //generator()->log() << "Q^2 VALUE IS OUT OF RANGE (TOO HIGH) = "
    //		       << q2 << endl;
  }
  
  xxx=x;
  
  // interpolation in logx, log qsq:
  xxx=log10(xxx);
  qsq=log10(q2);

  // NEW BIT STARTS HERE
  n=locate(xx,nx,xxx);
  m=locate(qq,nq,qsq);

  for(i=1;i<=np;i++) {
    t=(xxx-xx[n])/(xx[n+1]-xx[n]);
    u=(qsq-qq[m])/(qq[m+1]-qq[m]);
    table[i]=0.0;
    for(l=4;l>=1;l--) {
      table[i]=t*table[i]+((c[i][n][m][l][4]*u+c[i][n][m][l][3])*u +
			   c[i][n][m][l][2])*u+c[i][n][m][l][1];
    }
    //cout << "Table " << i << " = " << table[i] << endl;
  }
}

double MRST::xx[] =
  { 0.0, 1E-5, 2E-5, 4E-5, 6E-5, 8E-5, 1E-4, 2E-4, 4E-4, 6E-4, 8E-4,
    1E-3, 2E-3, 4E-3, 6E-3, 8E-3, 1E-2, 1.4E-2, 2E-2, 3E-2, 4E-2, 6E-2, 8E-2,
    .1, .125, 0.15, .175, .2, .225, 0.25, .275, .3, .325, 0.35, .375,
    .4, .425, 0.45, .475, .5, .525, 0.55, .575, .6, .65, .7, .75, 
    .8, .9, 1. };

double MRST::qq[] = 
  { 0.0, 1.25, 1.5, 2., 2.5, 3.2, 4., 5., 6.4, 8., 10., 12., 18., 26., 40., 
    64., 1E2, 1.6E2, 2.4E2, 4E2, 6.4E2, 1E3, 1.8E3, 3.2E3, 5.6E3,
    1E4, 1.8E4, 3.2E4, 5.6E4, 1E5, 1.8E5, 3.2E5, 5.6E5, 1E6, 1.8E6,
    3.2E6, 5.6E6, 1E7 };
