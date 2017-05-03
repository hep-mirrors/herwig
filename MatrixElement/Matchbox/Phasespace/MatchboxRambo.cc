// -*- C++ -*-
//
// MatchboxRambo.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxRambo class.
//

#include "MatchboxRambo.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/GSLBisection.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxRambo::MatchboxRambo() 
  : needToReshuffle(false), theMakeReferenceSample(false),
    referenceSample(0) {}

MatchboxRambo::~MatchboxRambo() {}

IBPtr MatchboxRambo::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxRambo::fullclone() const {
  return new_ptr(*this);
}

static double weights[7] = {

  -1.,-1.,
  0.039788735772973833942,
  0.00012598255637968550463,
  1.3296564302788840628E-7,
  7.0167897579949011130E-11,
  2.2217170114046130768E-14

};

void MatchboxRambo::setXComb(tStdXCombPtr xc) {
  MatchboxPhasespace::setXComb(xc);
  needToReshuffle = false;
  if ( xc ) {
    for ( cPDVector::const_iterator d = mePartonData().begin();
	  d != mePartonData().end(); ++d ) {
      if ( (**d).hardProcessMass() != ZERO ) {
	needToReshuffle = true;
	break;
      }
    }
  }
}

void MatchboxRambo::dumpReference(const vector<Lorentz5Momentum>& momenta, double weight) const {
  *referenceSample << lastX1() << " " << lastX2() << " ";
  Boost toLab = (lastPartons().first->momentum() + 
		 lastPartons().second->momentum()).boostVector();
  for ( vector<Lorentz5Momentum>::const_iterator p = momenta.begin();
	p != momenta.end(); ++p ) {
    Lorentz5Momentum pl = *p;
    if ( toLab.mag2() > Constants::epsilon )
      pl.boost(toLab);
    *referenceSample 
      << (pl.x()/GeV) << " "
      << (pl.y()/GeV) << " "
      << (pl.z()/GeV) << " "
      << (pl.t()/GeV) << " "
      << (pl.mass()/GeV) << " ";
  }
  double ymax = lastCuts().yHatMax();
  double ymin = lastCuts().yHatMin();
  double km = log(lastCuts().sHatMax()/lastCuts().sHatMin());
  ymax = min(ymax, log(lastCuts().x1Max()*sqrt(lastS()/lastSHat())));
  ymin = max(ymin, -log(lastCuts().x2Max()*sqrt(lastS()/lastSHat())));
  *referenceSample << weight*km*(ymax-ymin)/(lastX1()*lastX2()) << "\n" << flush;
}

double MatchboxRambo::generateTwoToNKinematics(const double* r,
					       vector<Lorentz5Momentum>& momenta) {

  if ( theMakeReferenceSample ) {
    map<cPDVector,ofstream*>::iterator ref =
      referenceSamples.find(mePartonData());
    if ( ref == referenceSamples.end() ) {
      ostringstream refname;
      for ( cPDVector::const_iterator p = mePartonData().begin();
	    p != mePartonData().end(); ++p ) {
	refname << (**p).PDGName();
      }
      refname << ".rambo";
      referenceSamples[mePartonData()] = new ofstream(refname.str().c_str(),std::ios_base::app);
      ref = referenceSamples.find(mePartonData());
      *(ref->second) << setprecision(26);
    }
    assert(ref != referenceSamples.end());
    referenceSample = ref->second;
  }

  size_t offset = dynamic_cast<const Tree2toNDiagram&>(*lastXComb().diagrams().front()).nSpace() > 0 ? 2 : 1;

  Energy w = sqrt(lastSHat());
  size_t count = 0;
  Lorentz5Momentum Q;
  for ( vector<Lorentz5Momentum>::iterator k = momenta.begin() + offset;
	k != momenta.end(); ++k ) {
    Energy q = -w*log(r[count]*r[count+1]);
    double ct = 2.*r[count+2]-1.;
    double st = sqrt(1.-sqr(ct));
    double phi = 2.*Constants::pi*r[count+3];
    double cphi = cos(phi);
    double sphi = sqrt(1.-sqr(cphi));
    if ( phi > Constants::pi )
      sphi = -sphi;
    (*k).setMass(ZERO);
    (*k).setT(q);
    (*k).setX(q*cphi*st);
    (*k).setY(q*sphi*st);
    (*k).setZ(q*ct);
    count += 4;
    Q += *k;
  }

  Energy M = sqrt(Q.m2());
  double x = w/M;
  Boost beta = -(Q.vect() * (1./M));
  double gamma = Q.t()/M;
  double a = 1./(1.+gamma);

  for ( vector<Lorentz5Momentum>::iterator k = momenta.begin() + offset;
	k != momenta.end(); ++k ) {
    Energy q = (*k).t();
    Energy bq = beta*(*k).vect();
    (*k).setT(x*(gamma*q+bq));
    (*k).setVect(x*((*k).vect()+(q+a*bq)*beta));
  }

  size_t n = momenta.size()-offset;
  double weight = weights[n];

  if ( !needToReshuffle ) {
    if ( !matchConstraints(momenta) )
      return 0.;
    fillDiagramWeights();
    if ( theMakeReferenceSample )
      dumpReference(momenta, weight);
    return weight;
  }

  double xi;

  ReshuffleEquation solve(w,mePartonData().begin()+offset,mePartonData().end(),
			  momenta.begin()+2,momenta.end());

  GSLBisection solver(1e-10,1e-8,10000);

  try {
    xi = solver.value(solve,0.0,1.1);
  } catch (GSLBisection::GSLerror) {
    return 0.;
  } catch (GSLBisection::IntervalError) {
    return 0.;
  }

  weight *= pow(xi,3.*(n-1.));

  Energy num = ZERO;
  Energy den = ZERO;

  cPDVector::const_iterator d = mePartonData().begin()+offset;
  for ( vector<Lorentz5Momentum>::iterator k = momenta.begin()+offset;
	k != momenta.end(); ++k, ++d ) {
    num += (*k).vect().mag2()/(*k).t();
    Energy q = (*k).t();
    (*k).setT(sqrt(sqr((**d).hardProcessMass())+xi*xi*sqr((*k).t())));
    (*k).setVect(xi*(*k).vect());
    weight *= q/(*k).t();
    den += (*k).vect().mag2()/(*k).t();
    (*k).setMass((**d).hardProcessMass());
  }

  if ( !matchConstraints(momenta) )
    return 0.;

  weight *= num/den;

  fillDiagramWeights();

  if ( theMakeReferenceSample )
    dumpReference(momenta, weight);
  
  return weight;

}

Energy MatchboxRambo::ReshuffleEquation::operator() (double xi) const {
  cPDVector::const_iterator d = dataBegin;
  vector<Lorentz5Momentum>::const_iterator p = momentaBegin;
  Energy res = -w;
  for ( ; d != dataEnd; ++d, ++p ) {
    res += sqrt(sqr((**d).hardProcessMass()) +
		xi*xi*sqr(p->t()));
  }
  return res;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxRambo::persistentOutput(PersistentOStream & os) const {
  os << needToReshuffle << theMakeReferenceSample;
}

void MatchboxRambo::persistentInput(PersistentIStream & is, int) {
  is >> needToReshuffle >> theMakeReferenceSample;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxRambo,MatchboxPhasespace>
  describeHerwigMatchboxRambo("Herwig::MatchboxRambo", "Herwig.so");

void MatchboxRambo::Init() {

  static ClassDocumentation<MatchboxRambo> documentation
    ("MatchboxRambo implements RAMBO phase space generation.");


  static Switch<MatchboxRambo,bool> interfaceMakeReferenceSample
    ("MakeReferenceSample",
     "Switch on generation of a reference sample of phase space points.",
     &MatchboxRambo::theMakeReferenceSample, false, false, false);
  static SwitchOption interfaceMakeReferenceSampleYes
    (interfaceMakeReferenceSample,
     "Yes",
     "Generate a reference sample.",
     true);
  static SwitchOption interfaceMakeReferenceSampleNo
    (interfaceMakeReferenceSample,
     "No",
     "Do not generate a reference sample.",
     false);
  interfaceMakeReferenceSample.rank(-1);

}

