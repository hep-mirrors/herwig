// -*- C++ -*-
//
// MatchboxPhasespace.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxPhasespace class.
//

#include "MatchboxPhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/ProcessData.h"
#include "Herwig++/MatrixElement/Matchbox/MatchboxFactory.h"

using namespace Herwig;

MatchboxPhasespace::MatchboxPhasespace() 
  : singularCutoff(10*GeV), theUseMassGenerators(false) {}

MatchboxPhasespace::~MatchboxPhasespace() {}

void MatchboxPhasespace::cloneDependencies(const std::string&) {}

Ptr<MatchboxFactory>::tcptr MatchboxPhasespace::factory() const {
  return lastMatchboxXComb()->factory();
}

Ptr<ProcessData>::tptr MatchboxPhasespace::processData() const {
  return factory()->processData();
}

double MatchboxPhasespace::generateKinematics(const double* r,
					      vector<Lorentz5Momentum>& momenta) {

  cPDVector::const_iterator pd = mePartonData().begin() + 2;
  vector<Lorentz5Momentum>::iterator p = momenta.begin() + 2;

  double massJacobian = 1.;

  if ( useMassGenerators() ) {
    Energy gmass = ZERO;
    tGenericMassGeneratorPtr mgen = 
      processData()->massGenerator(*pd);
    Energy maxMass = 
      (!haveX1X2() && momenta.size() > 3) ? 
      sqrt(lastSHat()) : sqrt(lastS());
    for ( ; pd != mePartonData().end(); ++pd, ++p ) {
      if ( mgen && !isInvertible() ) {
	Energy massMax = min((**pd).massMax(),maxMass);
	Energy massMin = (**pd).massMin();
	if ( massMin > massMax )
	  return 0.0;
	gmass = mgen->mass(massJacobian,**pd,massMin,massMax,r[0]);
	++r;
      } else if ( (**pd).width() != ZERO ) {
	Energy massMax = min((**pd).massMax(),maxMass);
	Energy massMin = (**pd).massMin();
	if ( massMin > massMax )
	  return 0.0;
	// use a standard Breit Wigner here which we can invert
	// see invertKinematics as well
	double bwILow = 
	  atan((sqr(massMin)-sqr((**pd).mass()))/((**pd).mass() * (**pd).width()));
	double bwIUp = 
	  atan((sqr(massMax)-sqr((**pd).mass()))/((**pd).mass() * (**pd).width()));
	gmass = sqrt(sqr((**pd).mass()) + 
		     (**pd).mass()*(**pd).width()*tan(bwILow+r[0]*(bwIUp-bwILow)));
	++r;
      } else {
	gmass = (**pd).mass();
      }
      maxMass -= gmass;
      p->setMass(gmass);
    }
  } else {
    for ( ; pd != mePartonData().end(); ++pd, ++p )
      p->setMass((**pd).mass());
  }

  double weight = momenta.size() > 3 ? 
    generateTwoToNKinematics(r,momenta) : 
    generateTwoToOneKinematics(r,momenta);

  return weight*massJacobian;

}

double MatchboxPhasespace::generateTwoToOneKinematics(const double* r,
						      vector<Lorentz5Momentum>& momenta) {

  double tau = momenta[2].mass2()/lastXCombPtr()->lastS();
  double ltau = log(tau)/2.;
  double y = ltau - 2.*r[0]*ltau;
  double x1 = sqrt(tau)*exp(y);
  double x2 = sqrt(tau)*exp(-y);

  ThreeVector<Energy> p1 =
    x1*(lastXCombPtr()->lastParticles().first->momentum().vect());

  ThreeVector<Energy> p2 =
    x2*(lastXCombPtr()->lastParticles().second->momentum().vect());

  ThreeVector<Energy> q = p1 + p2;

  momenta[0] = Lorentz5Momentum(momenta[0].mass(),p1);
  momenta[1] = Lorentz5Momentum(momenta[1].mass(),p2);
  momenta[2] = Lorentz5Momentum(momenta[2].mass(),q);

  lastXCombPtr()->lastX1X2(make_pair(x1,x2));
  lastXCombPtr()->lastSHat((momenta[0]+momenta[1]).m2());

  fillDiagramWeights();

  return -4.*Constants::pi*ltau;

}

double MatchboxPhasespace::invertKinematics(const vector<Lorentz5Momentum>& momenta,
					    double* r) const {

  if ( useMassGenerators() ) {

    Energy gmass = ZERO;
    Energy maxMass = 
      (!haveX1X2() && momenta.size() > 3) ? 
      sqrt((momenta[0]+momenta[1]).m2()) : sqrt(lastS());

    cPDVector::const_iterator pd = mePartonData().begin() + 2;
    vector<Lorentz5Momentum>::const_iterator p = momenta.begin() + 2;

    for ( ; pd != mePartonData().end(); ++pd, ++p ) {

      if ( (**pd).width() != ZERO ) {
	Energy massMax = min((**pd).massMax(),maxMass);
	Energy massMin = (**pd).massMin();
	if ( massMin > massMax )
	  return 0.0;
	double bwILow = 
	  atan((sqr(massMin)-sqr((**pd).mass()))/((**pd).mass() * (**pd).width()));
	double bwIUp = 
	  atan((sqr(massMax)-sqr((**pd).mass()))/((**pd).mass() * (**pd).width()));
	gmass = p->mass();
	double bw =
	  atan((sqr(gmass)-sqr((**pd).mass()))/((**pd).mass() * (**pd).width()));
	r[0] = (bw-bwILow)/(bwIUp-bwILow);
	++r;
      } else {
	gmass = (**pd).mass();
      }
      maxMass -= gmass;

    }

  }

  return momenta.size() > 3 ? 
    invertTwoToNKinematics(momenta,r) : 
    invertTwoToOneKinematics(momenta,r);
}

double MatchboxPhasespace::invertTwoToOneKinematics(const vector<Lorentz5Momentum>& momenta,
						    double* r) const {

  double tau = momenta[2].mass2()/lastXCombPtr()->lastS();
  double ltau = log(tau)/2.;

  r[0] = (ltau - (momenta[0]+momenta[1]).rapidity())/(2.*ltau);

  return -4.*Constants::pi*ltau;

}

pair<double,Lorentz5Momentum> 
MatchboxPhasespace::timeLikeWeight(const Tree2toNDiagram& diag,
				   int branch, double flatCut) const {

  pair<int,int> children = diag.children(branch);

  if ( children.first == -1 ) {
    return make_pair(1.,meMomenta()[diag.externalId(branch)]);
  }

  pair<double,Lorentz5Momentum> res
    = timeLikeWeight(diag,children.first,flatCut);

  pair<double,Lorentz5Momentum> other
    = timeLikeWeight(diag,children.second,flatCut);

  res.first *= other.first;
  res.second += other.second;

  Energy2 mass2 = sqr(diag.allPartons()[branch]->mass());
  Energy2 width2 = sqr(diag.allPartons()[branch]->width());

  if ( width2 == ZERO ) {
    if ( abs((res.second.m2()-mass2)/lastSHat()) > flatCut )
      res.first /=
	abs((res.second.m2()-mass2)/lastSHat());
  } else {
    res.first /=
      sqr((res.second.m2()-mass2)/lastSHat()) +
      mass2*width2/sqr(lastSHat());
  }

  return res;

}

double MatchboxPhasespace::spaceLikeWeight(const Tree2toNDiagram& diag,
					   const Lorentz5Momentum& incoming,
					   int branch, double flatCut) const {

  if ( branch == -1 )
    return 1.;

  pair<int,int> children = diag.children(branch);

  pair<double,Lorentz5Momentum> res =
    timeLikeWeight(diag,children.second,flatCut);

  if ( children.first == diag.nSpace() - 1 ) {
    return res.first;
  }

  res.second = incoming - res.second;

  Energy2 mass2 = sqr(diag.allPartons()[branch]->mass());
  Energy2 width2 = sqr(diag.allPartons()[branch]->width());

  if ( width2 == ZERO ) {
    if ( abs((res.second.m2()-mass2)/lastSHat()) > flatCut )
      res.first /=
	abs((res.second.m2()-mass2)/lastSHat());
  } else {
    res.first /=
      sqr((res.second.m2()-mass2)/lastSHat()) +
      mass2*width2/sqr(lastSHat());
  }

  return
    res.first * spaceLikeWeight(diag,res.second,children.first,flatCut);

}

void MatchboxPhasespace::fillDiagramWeights(double flatCut) {

  diagramWeights().clear();

  for ( StandardXComb::DiagramVector::const_iterator d =
	  lastXComb().diagrams().begin(); d != lastXComb().diagrams().end(); ++d ) {
    diagramWeights()[(**d).id()] = 
      spaceLikeWeight(dynamic_cast<const Tree2toNDiagram&>(**d),meMomenta()[0],0,flatCut);
  }

}

Selector<MEBase::DiagramIndex> 
MatchboxPhasespace::selectDiagrams(const MEBase::DiagramVector& diags) const {
  Selector<MEBase::DiagramIndex> ret;
  for ( MEBase::DiagramIndex d = 0; d < diags.size(); ++d ) {
    ret.insert(diagramWeight(dynamic_cast<const Tree2toNDiagram&>(*diags[d])),d);
  }
  return ret;
}

bool MatchboxPhasespace::matchConstraints(const vector<Lorentz5Momentum>& momenta) {

  if ( singularLimits().empty() )
    return true;

  lastSingularLimit() = singularLimits().begin();

  for ( ; lastSingularLimit() != singularLimits().end(); ++lastSingularLimit() ) {
    if ( lastSingularLimit()->first == lastSingularLimit()->second &&
	 momenta[lastSingularLimit()->first].t() < singularCutoff )
      break;
    if ( lastSingularLimit()->first != lastSingularLimit()->second &&
	 sqrt(momenta[lastSingularLimit()->first]*
	      momenta[lastSingularLimit()->second]) < singularCutoff ) {
      bool match = true;
      for ( set<pair<size_t,size_t> >::const_iterator other =
	      singularLimits().begin(); other != singularLimits().end(); ++other ) {
	if ( other == lastSingularLimit() )
	  continue;
	if ( other->first == other->second &&
	     momenta[other->first].t() < singularCutoff ) {
	  match = false;
	  break;
	}
	if ( other->first != other->second &&
	     sqrt(momenta[other->first]*
		  momenta[other->second]) < singularCutoff ) {
	  match = false;
	  break;
	}
      }
      if ( match )
	break;
    }
  }

  return lastSingularLimit() != singularLimits().end();

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void MatchboxPhasespace::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb
     << ounit(singularCutoff,GeV) << theUseMassGenerators;
}

void MatchboxPhasespace::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb
     >> iunit(singularCutoff,GeV) >> theUseMassGenerators;
  lastMatchboxXComb(theLastXComb);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<MatchboxPhasespace,HandlerBase>
  describeMatchboxPhasespace("Herwig::MatchboxPhasespace", "HwMatchbox.so");

void MatchboxPhasespace::Init() {

  static ClassDocumentation<MatchboxPhasespace> documentation
    ("MatchboxPhasespace defines an abstract interface to a phase "
     "space generator.");


  static Parameter<MatchboxPhasespace,Energy> interfaceSingularCutoff
    ("SingularCutoff",
     "[debug] Cutoff below which a region is considered singular.",
     &MatchboxPhasespace::singularCutoff, GeV, 10.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


  static Switch<MatchboxPhasespace,bool> interfaceUseMassGenerators
    ("UseMassGenerators",
     "Use mass generators instead of fixed masses.",
     &MatchboxPhasespace::theUseMassGenerators, false, false, false);
  static SwitchOption interfaceUseMassGeneratorsYes
    (interfaceUseMassGenerators,
     "Yes",
     "Use mass generators.",
     true);
  static SwitchOption interfaceUseMassGeneratorsNo
    (interfaceUseMassGenerators,
     "No",
     "Do not use mass generators.",
     false);

}

