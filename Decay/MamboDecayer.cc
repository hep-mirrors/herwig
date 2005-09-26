// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MamboDecayer class.
//

#include "MamboDecayer.h"
#include <ThePEG/Interface/ClassDocumentation.h>

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MamboDecayer.tcc"
#endif

#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/StandardRandom.h>
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"

using namespace Herwig;
using namespace ThePEG;

MamboDecayer::~MamboDecayer() {}

bool MamboDecayer::accept(const DecayMode & dm) const {
  return true;
}

void MamboDecayer::persistentOutput(PersistentOStream & os) const {
  os << _maxweight << global;
}

void MamboDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _maxweight >> global;
}

ClassDescription<MamboDecayer> MamboDecayer::initMamboDecayer;

void MamboDecayer::Init() {
  
  static ClassDocumentation<MamboDecayer> documentation
    ("This is the implementation of the MAMBO algorithm of Kleiss & "
     "Stirling for multi-particle phase-space decays.");

  
  static Parameter<MamboDecayer,double> interfaceMaximumWeight
    ("MaxWeight",
     "Maximum phase-space weight",
     &MamboDecayer::_maxweight, 10.0, 1.0, 50.,
     false, false, true);

  
  static Reference<MamboDecayer,GlobalParameters> interfaceGlobalParameters
    ("GlobalParameters",
     "The class that has effectiveGluonMass",
     &MamboDecayer::global, false, false, true, false);
 
}

ParticleVector MamboDecayer::decay(const DecayMode & dm,
				   const Particle & parent) const {
  PDVector children = dm.orderedProducts();
  const int N = children.size();
  ParticleVector out(N);
  if(N == 1) {
      out[0] = children[0]->produceParticle(parent.momentum());
      return out;
    }
  double totalMass(0.0);
  vector<Lorentz5Momentum> productMomentum(N);
  Energy gluMass = global->effectiveGluonMass();
  for(int i = 0; i < N; ++i) {
      if (children[i]->id() == 21) {
	productMomentum[i].setMass(gluMass);
	}
      else {
	productMomentum[i].setMass(children[i]->constituentMass());
      }
	totalMass += children[i]->constituentMass();
  }

  if(totalMass > parent.mass()) {
      generator()->log() << "MamboDecayer: The Decay mode " 
			 << dm.tag() << " cannot "
			 << "proceed, not enough phase space\n";
      out.clear();
      return out;
    }

  double wgt(0.);
  do {
    wgt = calculateMomentum(productMomentum,parent.mass());
    wgt *= UseRandom::rnd();
  }
  while(wgt > _maxweight);
  
  //set output momenta
  int iter  =  0;
  for(PDVector::const_iterator it=children.begin();iter<N;
      ++iter,++it) {
    out[iter] = (*it)->produceParticle(productMomentum[iter]);
  }

  // fix 3-gluon colour lines / general colour lines
  if(N == 3 && out[0]->id()==ParticleID::g && out[1]->id()==ParticleID::g
     && out[2]->id()== ParticleID::g ) 
    {
      out[0]->antiColourNeighbour(out[2]);
      out[0]->colourNeighbour(out[1]);
      out[1]->colourNeighbour(out[2]);
    }
  // fix top-quark decay colour-lines
  else if(abs(parent.id())==ParticleID::t) 
    {          
      for(int i=0;i < N;++i){ 
	  if(abs(out[i]->id())==ParticleID::b) {
	      parent.colourNeighbour(out[i]);	      
	    }
	  else {
	      if(out[i]->hasColour())
		out[i]->antiColourNeighbour(out[i + 1]);
	      if ( out[i]->hasAntiColour() )
		out[i]->colourNeighbour(out[i + 1]);
	      ++i;
	    }
      }
    }
  //everything else
  else  {
      for ( int i = 0; i < N; ++i ) { 	
	  if ( !out[i]->coloured() ) 
	    continue;
	  else if(i + 1 >= N) {
	      cerr << "Problem with colour lines in MamboDecayer for" 
		   << " decay mode\n"  << dm.tag() << endl;
	    }
	  if ( out[i]->hasColour() )
	    out[i]->antiColourNeighbour(out[i + 1]);
	  if ( out[i]->hasAntiColour() )
	    out[i]->colourNeighbour(out[i + 1]);
	  ++i; // skip the one that's linked up already
	}
    }
  finalBoost(parent, out); 
  setScales(parent, out); // in Decayer.cc

  return out;
}


double MamboDecayer::calculateMomentum(vector<Lorentz5Momentum> & mom,
				       const Energy & comEn) const
{
  const int N = mom.size();
  double rmtot(0.0);
  double rm2tot(0.0);
  for(int i = 0;i < N;++i) 
    {
      rmtot += mom[i].mass();
      rm2tot += mom[i].mass2();
    }

  long double wb = N*N*comEn*comEn-rmtot*rmtot;
  wb = wb/(N*N*(N-1.0));
  wb = sqrt(wb) - rmtot/N;
  long double wmin = 0.5*wb;
  long double wmax = (2.0/3.0)*wb;
  const long double tol(1e-12);
  long double sf,sf1,sff1,sm2f2;
  long double wold(wmax),r,f(0.),f1(0.),err,u0(0.),u1(0.),w(0.);
  do 
    {
      sf = 0.;sf1 = 0.;sff1 = 0.;sm2f2 = 0.;        
      for(int i = 0;i < N;++i) 
	{
	  r=abs(mom[i].mass()/wold);
	  if (r == 0.0) {
	    f=2.*wold;
	    f1=2.;
	  }
	  else {
	    f = wold*(2.0+r*(BessK0(r)/BessK1(r)));
	    f1 = 2.- (r*r*BessKPrime(r));
	  }
	  sf += f; sf1 += f1; sff1 += f*f1; sm2f2 -= f*f;
	}
      u0 = sf*sf+sm2f2+rm2tot;
      u1 = 2.*(sf*sf1-sff1);
      w = wold - (u0-comEn*comEn)/u1;
      err = abs(w-wold);
      wold = w;
    }
  while(err>tol && w > wmin);
  long double xu,xv;
  vector<long double> alpha(N),um(N),vm(N);
  for(int i = 0;i < N;++i) 
    { 
      alpha[i] = 2.*(mom[i].mass()/w);
      xu = (1.-alpha[i]+sqrt(1.+alpha[i]*alpha[i]))/2.;
      xv = (3.-alpha[i]+sqrt(9.+ 4.*alpha[i]+alpha[i]*alpha[i]))/2.;
      um[i] = exp(-xu/2.)*pow((xu*(xu+alpha[i])),0.25);
      vm[i] = xv*exp(-xv/2.)*pow((xv*(xv+alpha[i])),0.25);
    }
  
  //start k-momenta generation
  long double u(0.),v(0.),x(0.);
  vector<Lorentz5Momentum> qk(N);
  Lorentz5Momentum qktot;
  do 
    {
      qktot=LorentzMomentum();
      for(int i=0;i<N;++i) 
	{
	  long double usq(0.),bound(0.);
	  do 
	    {
	      u  =  UseRandom::rnd()*um[i];
	      v  =  UseRandom::rnd()*vm[i];
	      x  =  v/u;
	      usq  =  u*u;
	      bound  =  exp(-x)*sqrt(x*(x+alpha[i]));
	    }
	  while(usq>bound);
	  double ck,phi;
	  Kinematics::generateAngles(ck,phi);
	  double sk  =  sqrt(1.0-ck*ck);
	  double qkv  =  w*sqrt(x*(x+alpha[i]));
	  Lorentz5Momentum temp(qkv*sk*sin(phi),qkv*sk*cos(phi),qkv*ck,
				mom[i].mass()+w*x);
	  temp.rescaleMass();
	  qk[i]  =  temp;
	  qktot += qk[i];
	}
      qktot.rescaleMass();
      x = sqrt(comEn*comEn/qktot.mass2());
    }
  while(x>1.0);

      //Perform lorentz boost from k to q (use function from ThePEG????)
    vector<Lorentz5Momentum> q(N);
    long double q0=0.,q1=0.,q2=0.,q3=0.,t=0.;
    vector<long double> qsq(N);
   for(int i = 0;i<N;++i){
     q3 = qk[i]*qktot/qktot.mass(); t = (q3+qk[i](3))/(qktot(3)+qktot.mass());
     q2 = qk[i](2)-qktot(2)*t; q1 = qk[i](1)-qktot(1)*t;
     q0 = qk[i](0)-qktot(0)*t; 
     Lorentz5Momentum temp (x*q0,x*q1,x*q2,x*q3);
     temp.rescaleMass();
     q[i] = temp;
     qsq[i] = sqr(q[i][3])-x*x*mom[i].mass2();
   }
    
   long double xiold(1.),xi(0.); vector<long double> en(N);
   do {
     f = -comEn;
     f1 = 0.0;
     for(int i = 0;i<N;++i)	    {
       en[i] = sqrt((xiold*xiold*qsq[i]) + mom[i].mass2());
       f += en[i];
       f1 += qsq[i]/en[i];
     }
     xi = xiold - f/(xiold*f1);
     err = abs(xi-xiold);
     xiold = xi;
   }
   while(err>tol);
   //Now have desired momenta
   for(int i = 0;i < N;++i){
     Lorentz5Momentum temp(xi*q[i](0),xi*q[i](1),xi*q[i](2),en[i]);
     temp.rescaleMass();
     temp.rescaleMass();
     mom[i] = temp;
   }
   
   //Calculate weight of distribution
    double s1(1.),s2(0.),s3(0.),wxi(0.);
    for(int i=0;i<N;++i) {
      s1 *= q[i](3)/mom[i](3);
      s2 += mom[i].mass2()/q[i](3);
      s3 += mom[i].mass2()/mom[i](3);
    }
    wxi = pow(xi,(3*N-3))*s1*(comEn-x*x*s2)/(comEn-s3);
    
    return wxi;
}

  //  //Apply accept/reject 
//    if(reject(q,mom,comEn,xi,x)) 
//      {
//        MamboDecayer::calculateMomentum(mom,comEn);
//      }
//    else return;
// }

 // bool MamboDecayer::reject(const vector<Lorentz5Momentum> q,
// 			   const vector<Lorentz5Momentum> mom,
// 			   const Energy & comEn, const long double & xi,
// 			   const long double & r) const
// {
//   int N(mom.size());
//   double s1(1.),s2(0.),s3(0.),wxi(0.);
//   static StandardRandom gen;
//   long double maxweight(1.000001);
//   for(int i=0;i<N;++i) {
//       s1 *= q[i](3)/mom[i](3);
//       s2 += mom[i].mass2()/q[i](3);
//       s3 += mom[i].mass2()/mom[i](3);
//     }
//   wxi = pow(xi,(3*N-3))*s1*(comEn-r*r*s2)/(comEn-s3);
//   if(wxi>maxweight) {
//     generator()->log() << "MamxboDecayer::reject() - "
// 			<< "Weight max violation "
// 			<< wxi  << "   " << maxweight
// 			<< "!"  << endl;
//     }
//   double wgt = UseRandom::rnd()*maxweight;
//   if(wxi<wgt)
//     return true;
//   else
//     return false;
// }
