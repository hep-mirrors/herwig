// -*- C++ -*-
//
// MRST.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "MRST.h"
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Switch.h>
#include <istream>
#include <iostream>
#include <sstream>
#include <string>

using namespace ThePEG;
using namespace Herwig;

/**
 *  Minimum value of \f$x\f$
 */
const double MRST::xmin=1E-5;

/**
 *  Maximum value of \f$x\f$
 */
const double MRST::xmax=1.0;

/**
 *  Minimum value of \f$q^2\f$.
 */
const Energy2 MRST::qsqmin = 1.25 * GeV2;

/**
 *  Maximum value of \f$q^2\f$.
 */
const Energy2 MRST::qsqmax = 1E7 * GeV2;

/**
 *  Mass squared of the charm quark
 */
const Energy2 MRST::mc2 = 2.045 * GeV2;

/**
 *  Mass squared of the bottom quark
 */
const Energy2 MRST::mb2 = 18.5 * GeV2;

ClassDescription<MRST> MRST::initMRST;

MRST::MRST() : _inter(2), _xswitch(0.9),
	       data(np+1,vector<vector<double> >
		    (nx+1,vector<double>
		     (nq+1,0.0))),
	       fdata(np+1,vector<vector<double> >
		     (nx+1,vector<double>
		      (nq+1,0.0))) {
  if ( ! initialized ) {
    for ( int jj=1; jj < ntenth; ++jj ) {
      lxxb[jj] = log10(xx[jj]/xx[ntenth]) + xx[ntenth];
    }
    lxxb[ntenth] = xx[ntenth];
    for ( int n=1; n<=nx; n++ ) lxx[n] = log10(xx[n]);
    for ( int n=1; n<=nq; n++ ) lqq[n] = log10(qq[n]);
    initialized = true;
  }
}

bool MRST::canHandleParticle(tcPDPtr particle) const {
  // Return true if this PDF can handle the extraction of parton from the
  // given particle ie. if the particle is a proton or neutron.
  return ( abs(particle->id()) == ParticleID::pplus ||
           abs(particle->id()) == ParticleID::n0 );
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

double MRST::xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                 double x, double, Energy2) const {
  return pdfValue(x, partonScale, particle, parton,Total);
}

double MRST::xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                  double x, double, Energy2) const {
  return pdfValue(x, partonScale, particle, parton,Valence);
}

double MRST::xfsx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                  double x, double, Energy2) const {
  return pdfValue(x, partonScale, particle, parton,Sea);
}

double MRST::pdfValue(double x, Energy2 q2, 
		      tcPDPtr particle, tcPDPtr parton, PDFType type) const {
  assert(!isnan(x) && !isinf(x));
  // reset x  to min or max if outside range
  if(x<xmin)      x=xmin;
  else if(x>xmax) x=xmax;
  // reset q2 to min or max if outside range
  if(q2<qsqmin)      q2=qsqmin;
  else if(q2>qsqmax) q2=qsqmax;
  // c++ interpolation
  double output(0.);
  if(_inter==0||(_inter==2&&x<_xswitch)) {
    // interpolation is in logx, log qsq:
    double xxx=log10(x);
    double qsq=log10(q2/GeV2);
    
    // bin position
    int n=locate(lxx,nx,xxx);
    int m=locate(lqq,nq,qsq);
    
    // fraction along the bin
    double t=(xxx-lxx[n])/(lxx[n+1]-lxx[n]);
    double u=(qsq-lqq[m])/(lqq[m+1]-lqq[m]);
    
    bool anti = particle->id() < 0;
    bool neutron = abs(particle->id()) == ParticleID::n0;
    
    if (type==Valence) {
      switch(parton->id()) {
      case ParticleID::u:
	output= (neutron? 
		 (anti? 0.0: lookup(dnValence,n,m,u,t)): 
		 (anti? 0.0: lookup(upValence,n,m,u,t)));
	break;
      case ParticleID::ubar:
	output= (neutron? 
		 (anti? lookup(dnValence,n,m,u,t): 0.0): 
		 (anti? lookup(upValence,n,m,u,t): 0.0));
	break;
      case ParticleID::d:
	output= (neutron? 
		 (anti? 0.0: lookup(upValence,n,m,u,t)): 
		 (anti? 0.0: lookup(dnValence,n,m,u,t)));
	break;
      case ParticleID::dbar:
	output= (neutron? 
		 (anti? lookup(upValence,n,m,u,t): 0.0): 
		 (anti? lookup(dnValence,n,m,u,t): 0.0));
	break;
      }
    } 
    else if(type==Sea) {
      switch(parton->id()) {
      case ParticleID::b:
      case ParticleID::bbar:
	output= lookup(bot,n,m,u,t);
	break;
      case ParticleID::c:
      case ParticleID::cbar:
	output= lookup(chm,n,m,u,t);
	break;
      case ParticleID::s:
      case ParticleID::sbar:
	output= lookup(str,n,m,u,t);
	break;
      case ParticleID::u:
      case ParticleID::ubar:
	output= (neutron? lookup(dnSea,n,m,u,t) : lookup(upSea,n,m,u,t));
	break;
      case ParticleID::d:
      case ParticleID::dbar:
	output= (neutron? lookup(upSea,n,m,u,t) : lookup(dnSea,n,m,u,t));
	break;
      case ParticleID::g:
	output= lookup(glu,n,m,u,t);
	break;
      }
    }
    else if(type==Total) {
      switch(parton->id()) {
      case ParticleID::b:
      case ParticleID::bbar:
	output= lookup(bot,n,m,u,t);
	break;
      case ParticleID::c:
      case ParticleID::cbar:
	output= lookup(chm,n,m,u,t);
	break;
      case ParticleID::s:
      case ParticleID::sbar:
	output= lookup(str,n,m,u,t);
	break;
      case ParticleID::u:
	output= (neutron? 
		 (lookup(dnSea,n,m,u,t) + (anti? 0.0: lookup(dnValence,n,m,u,t))) :
		 (lookup(upSea,n,m,u,t) + (anti? 0.0: lookup(upValence,n,m,u,t))));
	break;
      case ParticleID::ubar:
	output= (neutron? 
		 (lookup(dnSea,n,m,u,t) + (anti? lookup(dnValence,n,m,u,t): 0.0)) :
		 (lookup(upSea,n,m,u,t) + (anti? lookup(upValence,n,m,u,t): 0.0)));
	break;
      case ParticleID::d:
	output= (neutron? 
		 (lookup(upSea,n,m,u,t) + (anti? 0.0: lookup(upValence,n,m,u,t))) :
		 (lookup(dnSea,n,m,u,t) + (anti? 0.0: lookup(dnValence,n,m,u,t))));
	break;
      case ParticleID::dbar:
	output= (neutron? 
		 (lookup(upSea,n,m,u,t) + (anti? lookup(upValence,n,m,u,t): 0.0)) :
		 (lookup(dnSea,n,m,u,t) + (anti? lookup(dnValence,n,m,u,t): 0.0)));
	break;
      case ParticleID::g:
	output= lookup(glu,n,m,u,t);
	break;
      }
    }
  }
  else {
    double xxx=x;
    if(x<lxxb[ntenth]) xxx = log10(x/lxxb[ntenth])+lxxb[ntenth];
    int nn=0;
    do ++nn;
    while(xxx>lxxb[nn+1]);
    double a=(xxx-lxxb[nn])/(lxxb[nn+1]-lxxb[nn]);
    double qsq=q2/GeV2;
    int mm=0;
    do ++mm;
    while(qsq>qq[mm+1]);
    double b=(qsq-qq[mm])/(qq[mm+1]-qq[mm]);
    double g[np+1];
    for(int ii=1;ii<=np;++ii) {
      g[ii]= (1.-a)*(1.-b)*fdata[ii][nn  ][mm] + (1.-a)*b*fdata[ii][nn  ][mm+1]
	+         a*(1.-b)*fdata[ii][nn+1][mm] +      a*b*fdata[ii][nn+1][mm+1];
      if(nn<ntenth&&!(ii==5||ii==7)) {
	double fac=(1.-b)*fdata[ii][ntenth][mm]+b*fdata[ii][ntenth][mm+1];
	g[ii] = fac*pow(10.,g[ii]-fac);
      }
      g[ii] *= pow(1.-x,n0[ii]);
    }
    bool anti = particle->id() < 0;
    bool neutron = abs(particle->id()) == ParticleID::n0;
    if (type==Valence) {
      switch(parton->id()) {
      case ParticleID::u:
	output= (neutron? 
		 (anti? 0.0: g[2]): 
		 (anti? 0.0: g[1]));
	break;
      case ParticleID::ubar:
	output= (neutron? 
		 (anti? g[2]: 0.0): 
		 (anti? g[1]: 0.0));
	break;
      case ParticleID::d:
	output= (neutron? 
		 (anti? 0.0: g[1]): 
		 (anti? 0.0: g[2]));
	break;
      case ParticleID::dbar:
	output= (neutron? 
		 (anti? g[1]: 0.0): 
		 (anti? g[2]: 0.0));
	break;
      }
    } 
    else if(type==Sea) {
      switch(parton->id()) {
      case ParticleID::b:
      case ParticleID::bbar:
	output= g[7];
	break;
      case ParticleID::c:
      case ParticleID::cbar:
	output= g[5];
	break;
      case ParticleID::s:
      case ParticleID::sbar:
	output= g[6];
	break;
      case ParticleID::u:
      case ParticleID::ubar:
	output= (neutron ? g[8] : g[4] );
	break;
      case ParticleID::d:
      case ParticleID::dbar:
	output= (neutron?  g[4] : g[8] );
	break;
      case ParticleID::g:
	output= g[3];
	break;
      }
    }
    else if(type==Total) {
      switch(parton->id()) {
      case ParticleID::b:
      case ParticleID::bbar:
	output= g[7];
	break;
      case ParticleID::c:
      case ParticleID::cbar:
	output= g[5];
	break;
      case ParticleID::s:
      case ParticleID::sbar:
	output= g[6];
	break;
      case ParticleID::u:
	output= (neutron? 
		 (g[8] + (anti? 0.0: g[2])) :
		 (g[4] + (anti? 0.0: g[1])));
	break;
      case ParticleID::ubar:
	output= (neutron? 
		 (g[8] + (anti? g[2]: 0.0)) :
		 (g[4] + (anti? g[1]: 0.0)));
	break;
      case ParticleID::d:
	output= (neutron? 
		 (g[4] + (anti? 0.0: g[1])) :
		 (g[8] + (anti? 0.0: g[2])));
	break;
      case ParticleID::dbar:
	output= (neutron? 
		 (g[4] + (anti? g[1]: 0.0)) :
		 (g[8] + (anti? g[2]: 0.0)));
	break;
      case ParticleID::g:
	output= g[3];
	break;
      }
    }
  }
  output = max(output,0.);
  assert(!isnan(output));
  return output;
}

void MRST::persistentOutput(PersistentOStream &out) const {
  out << _file << data << fdata << _inter << _xswitch;
}

void MRST::persistentInput(PersistentIStream & in, int) {
  in >> _file >> data >> fdata >> _inter >> _xswitch;
  initialize(false);
}

void MRST::Init() {

  static ClassDocumentation<MRST> documentation
    ("Implementation of the MRST PDFs",
     "Implementation of the MRST LO* / LO** PDFs \\cite{Sherstnev:2007nd}.",
     "  %\\cite{Sherstnev:2007nd}\n"
     "\\bibitem{Sherstnev:2007nd}\n"
     "  A.~Sherstnev and R.~S.~Thorne,\n"
     "  ``Parton Distributions for LO Generators,''\n"
     "  Eur.\\ Phys.\\ J.\\  C {\\bf 55} (2008) 553\n"
     "  [arXiv:0711.2473 [hep-ph]].\n"
     "  %%CITATION = EPHJA,C55,553;%%\n"
     );

  static Switch<MRST,unsigned int> interfaceInterpolation
    ("Interpolation",
     "Whether to use cubic or linear (C++ or FORTRAN) interpolation",
     &MRST::_inter, 2, false, false);
  static SwitchOption interfaceInterpolationCubic
    (interfaceInterpolation,
     "Cubic",
     "Use cubic interpolation",
     0);
  static SwitchOption interfaceInterpolationLinear
    (interfaceInterpolation,
     "Linear",
     "Use Linear Interpolation",
     1);;
  static SwitchOption interfaceInterpolationMixed
    (interfaceInterpolation,
     "Mixed",
     "Use cubic below xswitch and linear interpolation above",
     2);

  static Parameter<MRST,double> interfaceXSwitch
    ("XSwitch",
     "Value of x to switch from cubic to linear interpolation",
     &MRST::_xswitch, 0.9, 0.0, 1.0,
     false, false, Interface::limited);

}

void MRST::doinitrun() {

  cerr << "Warning: Herwig::MRST is obsolete and will be removed in Herwig 7.1.\n"
       << "         Please switch to using a PDF set provided by LHAPDF.\n";

  PDFBase::doinitrun();
#ifdef MRST_TESTING
  unsigned int intersave=_inter;
  tPDPtr proton=getParticleData(ParticleID::pplus);
  for(unsigned int itype=0;itype<8;++itype) {
    tPDPtr parton;
    string name;
    if(itype==0) {
      name="u.top";
      parton=getParticleData(ParticleID::u);
    }
    else if(itype==1) {
      name="d.top";
      parton=getParticleData(ParticleID::d);
    }
    else if(itype==2) {
      name="ubar.top";
      parton=getParticleData(ParticleID::ubar);
    }
    else if(itype==3) {
      name="dbar.top";
      parton=getParticleData(ParticleID::dbar);
    }
    else if(itype==4) {
      name="s.top";
      parton=getParticleData(ParticleID::s);
    }
    else if(itype==5) {
      name="c.top";
      parton=getParticleData(ParticleID::c);
    }
    else if(itype==6) {
      name="b.top";
      parton=getParticleData(ParticleID::b);
    }
    else if(itype==7) {
      name="g.top";
      parton=getParticleData(ParticleID::g);
    }
    ofstream output(name.c_str());
    Energy qmin=2.0*GeV,qmax=3000.0*GeV;
    int nq=10;
    Energy qstep=(qmax-qmin)/nq;
    for(Energy q=qmin+qstep;q<=qmax;q+=qstep) {
      double nx=500;
      double xmin=1e-5,xmax=1.;
      double xstep=(log(xmax)-log(xmin))/nx;
      output << "NEW FRAME"  << endl;
      output << "SET FONT DUPLEX\n";
      output << "SET SCALE Y LOG\n";
      output << "SET LIMITS X " << xmin << " " << xmax << endl;
      if(itype==0)
	output << "TITLE TOP \" up      distribution for q=" 
	       <<  q/GeV << "\"\n";
      else if(itype==1)
	output << "TITLE TOP \" down    distribution for q=" 
	       <<  q/GeV << "\"\n";
      else if(itype==2)
	output << "TITLE TOP \" ubar    distribution for q=" 
	       <<  q/GeV << "\"\n";
      else if(itype==3)
	output << "TITLE TOP \" dbar    distribution for q=" 
	       <<  q/GeV << "\"\n";
      else if(itype==4)
	output << "TITLE TOP \" strange distribution for q=" 
	       <<  q/GeV << "\"\n";
      else if(itype==5)
	output << "TITLE TOP \" charm   distribution for q=" 
	       <<  q/GeV << "\"\n";
      else if(itype==6)
	output << "TITLE TOP \" bottom  distribution for q=" 
	       <<  q/GeV << "\"\n";
      else if(itype==7)
	output << "TITLE TOP \" gluon   distribution for q=" 
	       <<  q/GeV << "\"\n";
      _inter=0;
      for(double xl=log(xmin)+xstep;xl<=log(xmax);xl+=xstep) {
	double x=exp(xl);
	double val=xfl(proton,parton,q*q,-xl);
	if(val>1e5) val=1e5;
	output << x << '\t' << val << '\n';
      }
      output << "JOIN RED" << endl;
      _inter=1;
      for(double xl=log(xmin)+xstep;xl<=log(xmax);xl+=xstep) {
	double x=exp(xl);
	double val=xfl(proton,parton,q*q,-xl);
	if(val>1e5) val=1e5;
	output << x << '\t' << val << '\n';
      }
      output << "JOIN BLUE" << endl;
      _inter=2;
      for(double xl=log(xmin)+xstep;xl<=log(xmax);xl+=xstep) {
	double x=exp(xl);
	double val=xfl(proton,parton,q*q,-xl);
	if(val>1e5) val=1e5;
	output << x << '\t' << val << '\n';
      }
      output << "JOIN GREEN" << endl;
    }
  }
  _inter=intersave;
#endif
}

void MRST::readSetup(istream &is) {
  _file = dynamic_cast<istringstream*>(&is)->str();
  initialize();
}

void MRST::initialize(bool reread) {
  useMe();
  //  int i,n,m,k,l,j; // counters
  double dx,dq;
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

  if(reread) {
    ifstream datafile(_file.c_str());
    if(!datafile) throw Exception() << "Could not open file '" << _file 
				    << "' in MRST::initialize()"
				    << Exception::runerror;

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
         if(datafile.eof()) throw Exception() << "Error while reading " << _file 
					      << " too few data points in file" 
					      << "in MRST::initialize()"
					      << Exception::runerror;
	 for(int ii=1;ii<=np;++ii) {
	   fdata[ii][nn][mm] = _inter==0 ? 0. :
	     data[ii][nn][mm]/pow(1.-xx[nn],n0[ii]);
	 }
       }
     }
     
     for (int n=1; n<=8; ++n) {
       for(int mm=1; mm<=nq; ++mm) {
	 data[n][nx][mm]=0.0;
       }
     }
     double dtemp;
     datafile >> dtemp;
     if(!datafile.eof()) throw Exception() << "Error reading end of " << _file 
					   << " too many data points in file" 
					   << "in MRST::initialize()"
					   << Exception::runerror;
     datafile.close();
     // calculate the FORTRAN interpolation
     for(int jj=1;jj<=ntenth-1;++jj) {
       for(int ii=1;ii<=np;++ii) {
	 if(ii==5||ii==7) continue;
	 for(int kk=1;kk<=nq;++kk) {
	   fdata[ii][jj][kk] = _inter==0 ? 0. : 
	     log10( fdata[ii][jj][kk] / fdata[ii][ntenth][kk] ) + 
	     fdata[ii][ntenth][kk];
	 }
       }
     }
     for (int n=1; n<=np; ++n) {
       for(int mm=1; mm<=nq; ++mm) {
	 fdata[n][nx][mm]=0.0;
       }
     }
  }

  // Now calculate the derivatives used for bicubic interpolation
  for (int i=1;i<=np;i++) {
    // Start by calculating the first x derivatives
    // along the first x value
    dx=lxx[2]-lxx[1];
    for (int m=1;m<=nq;m++)
      f1[i][1][m]=(data[i][2][m]-data[i][1][m])/dx;
    // The along the rest (up to the last)
    for (int k=2;k<nx;k++) {
      for (int m=1;m<=nq;m++) {
	f1[i][k][m]=polderivative(lxx[k-1],lxx[k],lxx[k+1],
				  data[i][k-1][m],
				  data[i][k][m],
				  data[i][k+1][m]);
      }
    }
    // Then for the last column
    dx=lxx[nx]-lxx[nx-1];
    for (int m=1;m<=nq;m++)
      f1[i][nx][m]=(data[i][nx][m]-data[i][nx-1][m])/dx;
    

    if ((i!=5)&&(i!=7)) {
      // then calculate the qq derivatives
      // Along the first qq value
      dq=lqq[2]-lqq[1];
      for (int k=1;k<=nx;k++)
	f2[i][k][1]=(data[i][k][2]-data[i][k][1])/dq;
      // The rest up to the last qq value
      for (int m=2;m<nq;m++) {
	for (int k=1;k<=nx;k++)
	  f2[i][k][m]=polderivative(lqq[m-1],lqq[m],lqq[m+1],
				    data[i][k][m-1],
				    data[i][k][m],
				    data[i][k][m+1]);
      }
      // then for the last row
      dq=lqq[nq]-lqq[nq-1];
      for (int k=1;k<=nx;k++)
	f2[i][k][nq]=(data[i][k][nq]-data[i][k][nq-1])/dq;
      
      // Now, calculate the cross derivatives.
      // Calculate these as x-derivatives of the y-derivatives
      // ?? Could be improved by taking the average between dxdy and dydx ??
      
      // Start by calculating the first x derivatives
      // along the first x value
      dx=lxx[2]-lxx[1];
      for (int m=1;m<=nq;m++)
	f12[i][1][m]=(f2[i][2][m]-f2[i][1][m])/dx;
      // The along the rest (up to the last)
      for (int k=2;k<nx;k++) {
	for (int m=1;m<=nq;m++)
	  f12[i][k][m]=polderivative(lxx[k-1],lxx[k],lxx[k+1],
				     f2[i][k-1][m],f2[i][k][m],f2[i][k+1][m]);
      }
      // Then for the last column
      dx=lxx[nx]-lxx[nx-1];
      for (int m=1;m<=nq;m++)
	f12[i][nx][m]=(f2[i][nx][m]-f2[i][nx-1][m])/dx;
    }
    
    if (i==5) {
      // zero all elements below the charm threshold
      for (int m=1;m<nqc0;m++) 
	for (int k=1;k<=nx;k++)
	  f2[i][k][m]=0.0; 

      // then calculate the qq derivatives 
      // Along the first qq value above the threshold (m=ncq0)
      dq=lqq[nqc0+1]-lqq[nqc0];
      for (int k=1;k<=nx;k++)
	f2[i][k][nqc0]=(data[i][k][nqc0+1]-data[i][k][nqc0])/dq;

      // The rest up to the last qq value
      for (int m=nqc0+1;m<nq;m++) {
	for (int k=1;k<=nx;k++)
	  f2[i][k][m]=polderivative(lqq[m-1],lqq[m],lqq[m+1],
				    data[i][k][m-1],
				    data[i][k][m],
				    data[i][k][m+1]);
      }
      // then for the last row
      dq=lqq[nq]-lqq[nq-1];
      for (int k=1;k<=nx;k++)
	f2[i][k][nq]=(data[i][k][nq]-data[i][k][nq-1])/dq;
      
      // Now, calculate the cross derivatives.
      // Calculate these as x-derivatives of the y-derivatives
      // ?? Could be improved by taking the average between dxdy and dydx ??

      dx=lxx[2]-lxx[1];
      for (int m=1;m<=nq;m++)
	f12[i][1][m]=(f2[i][2][m]-f2[i][1][m])/dx;
      // The along the rest (up to the last)
      for (int k=2;k<nx;k++) {
	for (int m=1;m<=nq;m++)
	  f12[i][k][m]=polderivative(lxx[k-1],lxx[k],lxx[k+1],
				     f2[i][k-1][m],f2[i][k][m],f2[i][k+1][m]);
      }
      // Then for the last column
      dx=lxx[nx]-lxx[nx-1];
      for (int m=1;m<=nq;m++)
	f12[i][nx][m]=(f2[i][nx][m]-f2[i][nx-1][m])/dx;
    }

    if (i==7) {
      // zero all elements below the bottom threshold
      for (int m=1;m<nqb0;m++) 
	for (int k=1;k<=nx;k++)
	  f2[i][k][m]=0.0; 

      // then calculate the qq derivatives
      // Along the first qq value above the threshold (m=nqb0)
      dq=lqq[nqb0+1]-lqq[nqb0];
      for (int k=1;k<=nx;k++)
	f2[i][k][nqb0]=(data[i][k][nqb0+1]-data[i][k][nqb0])/dq;

      // The rest up to the last qq value
      for (int m=nqb0+1;m<nq;m++) {
	for (int k=1;k<=nx;k++)
	  f2[i][k][m]=polderivative(lqq[m-1],lqq[m],lqq[m+1],
				    data[i][k][m-1],
				    data[i][k][m],
				    data[i][k][m+1]);
      }
      // then for the last row
      dq=lqq[nq]-lqq[nq-1];
      for (int k=1;k<=nx;k++)
	f2[i][k][nq]=(data[i][k][nq]-data[i][k][nq-1])/dq;
      
      // Now, calculate the cross derivatives.
      // Calculate these as x-derivatives of the y-derivatives
      // ?? Could be improved by taking the average between dxdy and dydx ??

      dx=lxx[2]-lxx[1];
      for (int m=1;m<=nq;m++)
	f12[i][1][m]=(f2[i][2][m]-f2[i][1][m])/dx;
      // The along the rest (up to the last)
      for (int k=2;k<nx;k++) {
	for (int m=1;m<=nq;m++)
	  f12[i][k][m]=polderivative(lxx[k-1],lxx[k],lxx[k+1],
				     f2[i][k-1][m],f2[i][k][m],f2[i][k+1][m]);
      }
      // Then for the last column
      dx=lxx[nx]-lxx[nx-1];
      for (int m=1;m<=nq;m++)
	f12[i][nx][m]=(f2[i][nx][m]-f2[i][nx-1][m])/dx;
    }
	
	
    // Now calculate the coefficients c_ij
    for(int n=1;n<=nx-1;n++) {
      for(int m=1;m<=nq-1;m++) {
	d1=lxx[n+1]-lxx[n];
	d2=lqq[m+1]-lqq[m];
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
	
	for (int k=1;k<=4;k++) {
	  x[k-1]=y[k];
	  x[k+3]=y1[k]*d1;
	  x[k+7]=y2[k]*d2;
	  x[k+11]=y12[k]*d1d2;
	}
	
	for (int l=0;l<=15;l++) {
	  xxd=0.0;
	  for (int k=0;k<=15;k++) xxd+= wt[l][k]*x[k];
	  cl[l]=xxd;
	}
	
	int l=0;
	for (int k=1;k<=4;k++) 
	  for (int j=1;j<=4;j++) c[i][n][m][k][j]=cl[l++];
      } //m
    } //n
  } // i
}

double MRST::xx[] =
  { 0.0, 1E-5, 2E-5, 4E-5, 6E-5, 8E-5, 1E-4, 2E-4, 4E-4, 6E-4, 8E-4,
    1E-3, 2E-3, 4E-3, 6E-3, 8E-3, 1E-2, 1.4E-2, 2E-2, 3E-2, 4E-2, 6E-2, 8E-2,
    .1, .125, 0.15, .175, .2, .225, 0.25, .275, .3, .325, 0.35, .375,
    .4, .425, 0.45, .475, .5, .525, 0.55, .575, .6, .65, .7, .75, 
    .8, .9, 1. };

double MRST::lxx[] =
  { 0.0, 1E-5, 2E-5, 4E-5, 6E-5, 8E-5, 1E-4, 2E-4, 4E-4, 6E-4, 8E-4,
    1E-3, 2E-3, 4E-3, 6E-3, 8E-3, 1E-2, 1.4E-2, 2E-2, 3E-2, 4E-2, 6E-2, 8E-2,
    .1, .125, 0.15, .175, .2, .225, 0.25, .275, .3, .325, 0.35, .375,
    .4, .425, 0.45, .475, .5, .525, 0.55, .575, .6, .65, .7, .75, 
    .8, .9, 1. };

double MRST::lxxb[] =
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

double MRST::lqq[] = 
  { 0.0, 1.25, 1.5, 2., 2.5, 3.2, 4., 5., 6.4, 8., 10., 12., 18., 26., 40., 
    64., 1E2, 1.6E2, 2.4E2, 4E2, 6.4E2, 1E3, 1.8E3, 3.2E3, 5.6E3,
    1E4, 1.8E4, 3.2E4, 5.6E4, 1E5, 1.8E5, 3.2E5, 5.6E5, 1E6, 1.8E6,
    3.2E6, 5.6E6, 1E7 };

double MRST::n0[] =
  {0,3,4,5,9,9,9,9,9};

bool MRST::initialized = false;
