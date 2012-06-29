// -*- C++ -*-
//
// MamboDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MamboDecayer class.
//

#include "MamboDecayer.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"

using namespace Herwig;
using namespace ThePEG;

bool MamboDecayer::accept(tcPDPtr, const tPDVector & ) const {
  return true;
}

void MamboDecayer::persistentOutput(PersistentOStream & os) const {
  os << _maxweight;
}

void MamboDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _maxweight;
}

ClassDescription<MamboDecayer> MamboDecayer::initMamboDecayer;

void MamboDecayer::Init() {
  
  static ClassDocumentation<MamboDecayer> documentation
    ("Decayer class that implements MAMBO algorithm of Kleiss-"
     "Stirling.",
     "The MAMBO algorithm of \\cite{Kleiss:1991rn} was used for high"
     " multiplicity decays",
     "\\bibitem{Kleiss:1991rn} R.~Kleiss and W.~J.~Stirling,\n"
     "Nucl.\\ Phys.\\  B {\\bf 385} (1992) 413.\n"
     "%%CITATION = NUPHA,B385,413;%%\n");
  
  static Parameter<MamboDecayer,double> interfaceMaximumWeight
    ("MaxWeight",
     "Maximum phase-space weight",
     &MamboDecayer::_maxweight, 10.0, 1.0, 50.,
     false, false, true);

}

ParticleVector MamboDecayer::decay(const Particle & parent,
				   const tPDVector & children) const {
  useMe();
  const int N = children.size();
  ParticleVector out(N);
  if(N == 1) {
    out[0] = children[0]->produceParticle(parent.momentum());
    return out;
  }
  Energy totalMass(ZERO);
  for(int i = 0; i < N; ++i) {
    totalMass += children[i]->mass();
  }
  if(totalMass > parent.mass()) {
    generator()->log() << "MamboDecayer::decay - The Decay mode " 
		       << parent.PDGName() << "->";
    for(int i = 0; i < N; ++i) 
      generator()->log() << children[i]->PDGName() << " ";
    generator()->log() << " cannot proceed as there is not "
		       << "enough phase space.\n";

    out.clear();
    return out;
  }
  
  double wgt(0.);
  vector<Lorentz5Momentum> productMomentum(N);
  do {
    for(int i = 0; i < N; ++i)
      productMomentum[i].setMass(children[i]->mass());
    wgt = calculateMomentum(productMomentum,parent.mass());
  }
  while(wgt < _maxweight*UseRandom::rnd());
  
  //set output momenta
  int iter  =  0;
  for(tPDVector::const_iterator it=children.begin();iter<N;
      ++iter,++it) {
    out[iter] = (*it)->produceParticle(productMomentum[iter]);
  }

  colourConnections(parent, out);
  finalBoost(parent, out); 
  setScales(parent, out);
  return out;
}


double MamboDecayer::calculateMomentum(vector<Lorentz5Momentum> & mom,
				       Energy comEn) const {
  const int N = mom.size();
  Energy rmtot(ZERO);
  Energy2 rm2tot(ZERO);
  for(int i = 0;i < N;++i) {
    rmtot += mom[i].mass();
    rm2tot += mom[i].mass2();
  }
  Energy wb = sqrt( (N*sqr(comEn) - sqr(rmtot))/N/N/(N - 1.) ) - rmtot/N;
  Energy wmax = (2.0/3.0)*wb;
  const Energy tol(1e-12*MeV);
  long double r(0.), sf1(0.);
  Energy2 sm2f2(ZERO);
  Energy sf(ZERO), sff1(ZERO), w(ZERO), 
    wold(wmax), err(ZERO);
  unsigned int iter(0), maxiter(50);
  do {
    sf = ZERO; sf1 = 0.; sff1 = ZERO; sm2f2 = ZERO;        
    for(int i = 0;i < N;++i) {
      r = abs(mom[i].mass()/wold);
      Energy f(ZERO);
      long double f1(0.);
      if (r == 0.0) {
	f=2.*wold;
	f1=2.;
      }
      else {
	long double fk0(0.), fkp(0.);
	BesselFns(r, fk0, fkp);
	f = wold*(2.0 + r*fk0);
	f1 = 2.- r*r*fkp;
      }
      sf += f; 
      sf1 += f1; 
      sff1 += f*f1; 
      sm2f2 -= f*f;
    }
    Energy u1 = 2.*(sf*sf1 - sff1);
    Energy2 u0 = sf*sf + sm2f2 + rm2tot - sqr(comEn);
    w = wold - u0/u1;
    err = abs(w - wold);
    wold  = w;
    ++iter;
  }
  while(err > tol && iter < maxiter);
  long double xu,xv;
  vector<long double> alpha(N),um(N),vm(N);
  for(int i = 0;i < N;++i) { 
    alpha[i] = 2.*(mom[i].mass()/w);
    xu = (1.-alpha[i]+sqrt(1.+alpha[i]*alpha[i]))/2.;
    xv = (3.-alpha[i]+sqrt(9.+ 4.*alpha[i]+alpha[i]*alpha[i]))/2.;
    um[i] = exp(-xu/2.)*pow((xu*(xu+alpha[i])),0.25l);
    vm[i] = xv*exp(-xv/2.)*pow((xv*(xv+alpha[i])),0.25l);
  }
  
  //start k-momenta generation
  long double u(0.),v(0.),x(0.);
  vector<Lorentz5Momentum> qk(N);
  Lorentz5Momentum qktot;
  do {
    qktot=LorentzMomentum();
    for(int i=0;i<N;++i) {
      long double usq(0.),bound(0.);
      do {
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
      Energy qkv  =  w*sqrt(x*(x+alpha[i]));
      qk[i] = Lorentz5Momentum(qkv*sk*sin(phi),qkv*sk*cos(phi),qkv*ck,
			       mom[i].mass()+w*x);
      qktot += qk[i];
    }
    qktot.rescaleMass();
    x = sqrt(comEn*comEn/qktot.mass2());
  }
  while(x>1.0);

  //Perform lorentz boost from k to q
  vector<Lorentz5Momentum> q(N);
  Energy q0=ZERO, q1=ZERO, q2=ZERO, q3=ZERO;
  long double t=0.;
  vector<Energy2> qsq(N);
  for(int i = 0;i<N;++i){
    q3 = (qk[i]*qktot)/qktot.mass(); 
    t = (q3+qk[i].e())/(qktot.e()+qktot.mass());
    q2 = qk[i].z()-qktot.z()*t; 
    q1 = qk[i].y()-qktot.y()*t;
    q0 = qk[i].x()-qktot.x()*t; 
    q[i] = Lorentz5Momentum(x*q0,x*q1,x*q2,x*q3);
    qsq[i] = sqr(q[i].e())-x*x*mom[i].mass2();
  }
  
  long double xiold(1.),xi(0.); 
  vector<Energy> en(N);
  iter = 0;
  do {
    Energy f = -comEn;
    Energy f1 = ZERO;
    for(int i = 0; i < N; ++i)	    {
      en[i] = sqrt((xiold*xiold*qsq[i]) + mom[i].mass2());
      f += en[i];
      f1 += qsq[i]/en[i];
    }
    xi = xiold - f/(xiold*f1);
    err = abs(xi-xiold)*MeV;
    xiold = xi;
    ++iter;
  }
  while(err > tol && iter < maxiter);
  //Now have desired momenta
  for(int i = 0;i < N;++i)
    mom[i] = Lorentz5Momentum(xi*q[i].x(),xi*q[i].y(),xi*q[i].z(),en[i]);
  
  //Calculate weight of distribution
  double s1(1.);
  Energy s2(ZERO),s3(ZERO);
  double wxi(0.);
  for(int i=0;i<N;++i) {
    s1 *= q[i].e()/mom[i].e();
    s2 += mom[i].mass2()/q[i].e();
    s3 += mom[i].mass2()/mom[i].e();
  }
  wxi = pow(xi,(3*N-3))*s1*(comEn-x*x*s2)/(comEn-s3);
  return wxi;
}

void MamboDecayer::colourConnections(const Particle & parent, 
				     ParticleVector & out) const {
  const int N = out.size();
  // fix 3-gluon colour lines / general colour lines
  if(N == 3 && out[0]->id()==ParticleID::g && out[1]->id()==ParticleID::g
     && out[2]->id()== ParticleID::g ) {
    out[0]->antiColourNeighbour(out[2]);
    out[0]->colourNeighbour(out[1]);
    out[1]->colourNeighbour(out[2]);
    }
  // incoming colour3/3bar state
  else if(parent.data().iColour() == PDT::Colour3 || 
	  parent.data().iColour() == PDT::Colour3bar) {          
    PPtr pparent = const_ptr_cast<PPtr>(&parent);
    if(N==2) {
      if(out[0]->data().iColour() == parent.data().iColour()) {
	out[0]->incomingColour(pparent,out[0]->id() < 0);
      }
      else if(out[1]->data().iColour() == parent.data().iColour()) {
	out[1]->incomingColour(pparent,out[1]->id() < 0);
      }
      else if(parent.data().iColour() == PDT::Colour3 &&
	      out[0]->data().iColour() == PDT::Colour3bar &&
	      out[1]->data().iColour() == PDT::Colour3bar) {
	tColinePtr col[2] = {ColourLine::create(out[0],true),
			     ColourLine::create(out[1],true)};
	parent.colourLine()->setSinkNeighbours(col[0],col[1]);
      }
      else if(parent.data().iColour() == PDT::Colour3bar &&
	      out[0]->data().iColour() == PDT::Colour3 &&
	      out[1]->data().iColour() == PDT::Colour3) {
	tColinePtr col[2] = {ColourLine::create(out[0]),
			     ColourLine::create(out[1])};
	parent.antiColourLine()->setSourceNeighbours(col[0],col[1]);
      }
      else {
	throw Exception() << "MamboDecayer::decay() can't make colour connections for "
			  << "the two-body decay of the coloured particle " 
			  << parent << Exception::eventerror;
      }
    }
    else {
      for(int i=0;i < N;++i){ 
	if(out[i]->data().iColour() == PDT::Colour3 ||
	   out[i]->data().iColour() == PDT::Colour3bar) {
	  out[i]->incomingColour(pparent,out[i]->id() < 0);
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
  }
  //incoming octet
  else if(parent.data().iColour() == PDT::Colour8) {
    PPtr pparent = const_ptr_cast<PPtr>(&parent);
    for(int i=0; i < N;++i) {
      if(out[i]->data().iColour() == PDT::Colour8) {
	out[i]->incomingColour(pparent);
	out[i]->incomingAntiColour(pparent);
      }
      else if(out[i]->data().iColour() == PDT::Colour3 ||
	      out[i]->data().iColour() == PDT::Colour3bar) {
	out[i]->incomingColour(pparent,out[i]->id() < 0);
      }
    }
  }
  else if(N==3 && parent.data().iColour() == PDT::Colour0&&
	  out[0]->data().iColour() == PDT::Colour3 &&
	  out[1]->data().iColour() == PDT::Colour3 &&
	  out[2]->data().iColour() == PDT::Colour3) {
    tColinePtr col[3] = {ColourLine::create(out[0]),
			 ColourLine::create(out[1]),
			 ColourLine::create(out[2])};
    col[0]->setSourceNeighbours(col[1],col[2]);
  }
  else if(N==3 && parent.data().iColour() == PDT::Colour0&&
	  out[0]->data().iColour() == PDT::Colour3bar &&
	  out[1]->data().iColour() == PDT::Colour3bar &&
	  out[2]->data().iColour() == PDT::Colour3bar) {
    tColinePtr col[3] = {ColourLine::create(out[0],true),
			 ColourLine::create(out[1],true),
			 ColourLine::create(out[2],true)};
    col[0]->setSinkNeighbours(col[1],col[2]);
  }
  //everything else
  else  {
    for ( int i = 0; i < N; ++i ) { 	
      if ( !out[i]->coloured() ) 
	continue;
      else if(i + 1 >= N) {
	throw Exception() 
	  << "MamboDecayer::colourConnections() - "
	  << "Cannot find appropriate configuration."
	  << Exception::warning;
      }
      if ( out[i]->hasColour() )
	out[i]->antiColourNeighbour(out[i + 1]);
      if ( out[i]->hasAntiColour() )
	out[i]->colourNeighbour(out[i + 1]);
      ++i; // skip the one that's linked up already
    }
  }
  
}

void MamboDecayer::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  output << "set " << name() << ":MaxWeight "  << _maxweight << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void MamboDecayer::doinitrun() {
  HwDecayerBase::doinitrun();
  _a0[0] = 0.5;
  _a0[1] = 0.375;
  _a0[2] = 0.375;
  _a0[3] = 0.4921875;
  _a0[4] = 0.84375;
  _a0[5] = 1.854492188,
  _a0[6] = 5.0625;
  _a0[7] = 16.58578491;
  _a0[8] = 63.33398438; 
  _a0[9] = 275.6161079;

  _a1[0] = 0.5;
  _a1[1] = 0.75;
  _a1[2] = 1.125;
  _a1[3] = 1.96875;
  _a1[4] = 4.21875;
  _a1[5] = 11.12695313;
  _a1[6] = 35.4375;
  _a1[7] = 132.6862793;
  _a1[8] = 570.0058594;
  _a1[9] = 2756.161079;
}
