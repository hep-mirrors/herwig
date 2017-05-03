// -*- C++ -*-
//
// MatchboxPhasespace.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/Utility/ProcessData.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

using namespace Herwig;

MatchboxPhasespace::MatchboxPhasespace() 
  : singularCutoff(10*GeV), theUseMassGenerators(false),
    theLoopParticleIdMin(200001), theLoopParticleIdMax(200100) {}

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
  Energy summ = ZERO;

  if ( useMassGenerators() ) {
    Energy gmass = ZERO;
    tGenericMassGeneratorPtr mgen;
    Energy maxMass = 
      (!haveX1X2() && momenta.size() > 3) ? 
      sqrt(lastSHat()) : sqrt(lastS());
    for ( ; pd != mePartonData().end(); ++pd, ++p ) {
      mgen = processData()->massGenerator(*pd);
      if ( mgen && !isInvertible() ) {
	Energy massMax = min((**pd).massMax(),maxMass);
	Energy massMin = (**pd).massMin();
	if ( massMin > massMax )
	  return 0.0;
	gmass = mgen->mass(massJacobian,**pd,massMin,massMax,r[0]);
	++r;
      } else if ( (**pd).hardProcessWidth() != ZERO ) {
	Energy massMax = min((**pd).massMax(),maxMass);
	Energy massMin = (**pd).massMin();
	if ( massMin > massMax )
	  return 0.0;
	// use a standard Breit Wigner here which we can invert
	// see invertKinematics as well
	double bwILow = 
	  atan((sqr(massMin)-sqr((**pd).hardProcessMass()))/((**pd).hardProcessMass() * (**pd).hardProcessWidth()));
	double bwIUp = 
	  atan((sqr(massMax)-sqr((**pd).hardProcessMass()))/((**pd).hardProcessMass() * (**pd).hardProcessWidth()));
	gmass = sqrt(sqr((**pd).hardProcessMass()) + 
		     (**pd).hardProcessMass()*(**pd).hardProcessWidth()*tan(bwILow+r[0]*(bwIUp-bwILow)));
	++r;
      } else {
	gmass = (**pd).hardProcessMass();
      }
      maxMass -= gmass;
      p->setMass(gmass);
      summ += gmass;
    }
  } else {
    for ( ; pd != mePartonData().end(); ++pd, ++p ) {
      summ += (**pd).hardProcessMass();
      p->setMass((**pd).hardProcessMass());
    }
  }

  if ( momenta.size() > 3 && !haveX1X2() ) {
    if ( summ > (momenta[0]+momenta[1]).m() )
      return 0.0;
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
  //old:  y = ltau - 2.*r[0]*ltau; x1 = sqrt(tau)*exp(y); x2 = sqrt(tau)*exp(-y);
  double x1=pow(tau,1.-r[0]);
  double x2=pow(tau,r[0]);

  // Due to the proton mass and P1.e() + P2.e() == lastS() we multiply here 
  // with the correction factor abs(P1.e()/P1.z()) to produce incoming 
  // p1/2 = (e1/2,0,0,+/- e1/2)  
  Lorentz5Momentum P1 = lastXCombPtr()->lastParticles().first->momentum();
  ThreeVector<Energy> p1 =  x1 * (P1.vect()) * abs(P1.e()/P1.z());

  Lorentz5Momentum P2 = lastXCombPtr()->lastParticles().second->momentum();
  ThreeVector<Energy> p2 =  x2 * (P2.vect()) * abs(P2.e()/P2.z());

  ThreeVector<Energy> q = p1 + p2;

  momenta[0] = Lorentz5Momentum(momenta[0].mass(),p1);
  momenta[1] = Lorentz5Momentum(momenta[1].mass(),p2);
  momenta[2] = Lorentz5Momentum(momenta[2].mass(),q);
  
  // check for energy conservation:
  if ((momenta[0]+momenta[1]-momenta[2]).e()>pow(10,-9)*GeV)
    generator()->log() 
    << "Warning: Momentum conservation in generateTwoToOneKinematics not precise.\n"
    << flush;
  

  lastXCombPtr()->lastX1X2({x1,x2});
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

      if ( (**pd).hardProcessWidth() != ZERO ) {
	Energy massMax = min((**pd).massMax(),maxMass);
	Energy massMin = (**pd).massMin();
	if ( massMin > massMax )
	  return 0.0;
	double bwILow = 
	  atan((sqr(massMin)-sqr((**pd).hardProcessMass()))/((**pd).hardProcessMass() * (**pd).hardProcessWidth()));
	double bwIUp = 
	  atan((sqr(massMax)-sqr((**pd).hardProcessMass()))/((**pd).hardProcessMass() * (**pd).hardProcessWidth()));
	gmass = p->mass();
	double bw =
	  atan((sqr(gmass)-sqr((**pd).hardProcessMass()))/((**pd).hardProcessMass() * (**pd).hardProcessWidth()));
	r[0] = (bw-bwILow)/(bwIUp-bwILow);
	++r;
      } else {
	gmass = (**pd).hardProcessMass();
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

void MatchboxPhasespace::setCoupling(long a, long b, long c, 
				     double coupling, bool includeCrossings) {
  cPDPtr A = getParticleData(a);
  cPDPtr B = getParticleData(b);
  cPDPtr C = getParticleData(c);
  if ( !A || !B || !C ) {
    generator()->log() << "Warning: could not determine particle data for ids "
		       << a << " " << b << " " << c << " when setting coupling in MatchboxPhasespace.\n"
		       << flush;
    return;
  }
  if ( !includeCrossings ) {
    theCouplings->couplings()[LTriple(a,b,c)] = coupling;
    return;
  }
  if ( A->CC() ) {
    theCouplings->couplings()[LTriple(-a,b,c)] = coupling;
    theCouplings->couplings()[LTriple(-a,c,b)] = coupling;
  } else {
    theCouplings->couplings()[LTriple(a,b,c)] = coupling;
    theCouplings->couplings()[LTriple(a,c,b)] = coupling;
  }
  if ( B->CC() ) {
    theCouplings->couplings()[LTriple(-b,a,c)] = coupling;
    theCouplings->couplings()[LTriple(-b,c,a)] = coupling;
  } else {
    theCouplings->couplings()[LTriple(b,a,c)] = coupling;
    theCouplings->couplings()[LTriple(b,c,a)] = coupling;
  }
  if ( C->CC() ) {
    theCouplings->couplings()[LTriple(-c,a,b)] = coupling;
    theCouplings->couplings()[LTriple(-c,b,a)] = coupling;
  } else {
    theCouplings->couplings()[LTriple(c,a,b)] = coupling;
    theCouplings->couplings()[LTriple(c,b,a)] = coupling;
  }
}

string MatchboxPhasespace::doSetCoupling(string in) {
  istringstream is(in);
  long a,b,c; double coupling;
  is >> a >> b >> c >> coupling;
  if ( !is )
    return "MatchboxPhasespace: error in setting coupling.";
  setCoupling(a,b,c,coupling,true);
  return "";
}

string MatchboxPhasespace::doSetPhysicalCoupling(string in) {
  istringstream is(in);
  long a,b,c; double coupling;
  is >> a >> b >> c >> coupling;
  if ( !is )
    return "MatchboxPhasespace: error in setting coupling.";
  setCoupling(a,b,c,coupling,false);
  return "";
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

  LTriple vertexKey(diag.allPartons()[branch]->id(),
		    diag.allPartons()[children.first]->id(),
		    diag.allPartons()[children.second]->id());
  map<LTriple,double>::const_iterator cit = theCouplings->couplings().find(vertexKey);
  if ( cit != theCouplings->couplings().end() ){
    res.first *= cit->second;
  }

  Energy2 mass2 = sqr(diag.allPartons()[branch]->hardProcessMass());
  Energy2 width2 = sqr(diag.allPartons()[branch]->hardProcessWidth());

  if ( abs(diag.allPartons()[branch]->id()) >= theLoopParticleIdMin
       && abs(diag.allPartons()[branch]->id()) <= theLoopParticleIdMax ) { // "loop particle"

    if ( abs((res.second.m2()-mass2)/lastSHat()) > flatCut ) {
      res.first /=
	abs((res.second.m2()-mass2)/GeV2);
      res.first *=
	log(abs((res.second.m2()-mass2)/GeV2)); // normal. of the argument in the log?
    }

  } else {

    if ( width2 == ZERO ) {
      if ( abs((res.second.m2()-mass2)/lastSHat()) > flatCut )
	res.first /=
	  abs((res.second.m2()-mass2)/GeV2);
    } else {
      res.first /=
	(sqr((res.second.m2()-mass2)/GeV2) +
	 mass2*width2/sqr(GeV2))/(abs(res.second.m2()/GeV2));
    }

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

  
  LTriple vertexKey(diag.allPartons()[branch]->id(),
		    diag.allPartons()[children.first]->id(),
		    diag.allPartons()[children.second]->id());
  if ( children.first == diag.nSpace() - 1 ) {
    if ( diag.allPartons()[children.first]->CC() )
      vertexKey = LTriple(diag.allPartons()[branch]->id(),
			  diag.allPartons()[children.second]->id(),
			  diag.allPartons()[children.first]->CC()->id());
    else
      vertexKey = LTriple(diag.allPartons()[branch]->id(),
			  diag.allPartons()[children.second]->id(),
			  diag.allPartons()[children.first]->id());
  }
  map<LTriple,double>::const_iterator cit = theCouplings->couplings().find(vertexKey);
  if ( cit != theCouplings->couplings().end() ){
    res.first *= cit->second;
  }
  if ( children.first == diag.nSpace() - 1 ) {
    return res.first;
  }

  res.second = incoming - res.second;

  Energy2 mass2 = sqr(diag.allPartons()[children.first]->hardProcessMass());
  Energy2 width2 = sqr(diag.allPartons()[children.first]->hardProcessWidth());

  if ( abs(diag.allPartons()[children.first]->id()) >= theLoopParticleIdMin
       && (diag.allPartons()[children.first]->id()) <= theLoopParticleIdMax ) { // "loop particle"

    if ( abs((res.second.m2()-mass2)/lastSHat()) > flatCut ) {
      res.first /=
	abs((res.second.m2()-mass2)/GeV2);
      res.first *=
	log(abs((res.second.m2()-mass2)/GeV2)); // normal. of the argument in the log?
    }

  } else {

    if ( width2 == ZERO ) {
      if ( abs((res.second.m2()-mass2)/lastSHat()) > flatCut )
	res.first /=
	  abs((res.second.m2()-mass2)/GeV2);
    } else {
      res.first /=
	(sqr((res.second.m2()-mass2)/GeV2) +
	 mass2*width2/sqr(GeV2))/(abs(res.second.m2()/GeV2));
    }

  }

  return
    res.first * spaceLikeWeight(diag,res.second,children.first,flatCut);

}

void MatchboxPhasespace::fillDiagramWeights(double flatCut) {

  diagramWeights().clear();

  for ( auto & d : lastXComb().diagrams() ) {
    diagramWeights()[d->id()] =
      spaceLikeWeight(dynamic_cast<const Tree2toNDiagram&>(*d),meMomenta()[0],0,flatCut);
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

int MatchboxPhasespace::nDim(const cPDVector& data) const {

  int ndimps = nDimPhasespace(data.size()-2);

  if ( useMassGenerators() ) {
    for ( cPDVector::const_iterator pd = data.begin();
	  pd != data.end(); ++pd ) {
      if ( (**pd).massGenerator() ||
	   (**pd).hardProcessWidth() != ZERO ) {
	++ndimps;
      }
    }
  }

  return ndimps;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void MatchboxPhasespace::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb
     << ounit(singularCutoff,GeV) << theUseMassGenerators 
     << theLoopParticleIdMin << theLoopParticleIdMax
     << theCouplings;
}

void MatchboxPhasespace::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb
     >> iunit(singularCutoff,GeV) >> theUseMassGenerators
     >> theLoopParticleIdMin >> theLoopParticleIdMax
     >> theCouplings;
  lastMatchboxXComb(theLastXComb);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<MatchboxPhasespace,HandlerBase>
  describeMatchboxPhasespace("Herwig::MatchboxPhasespace", "Herwig.so");

void MatchboxPhasespace::Init() {

  static ClassDocumentation<MatchboxPhasespace> documentation
    ("MatchboxPhasespace defines an abstract interface to a phase "
     "space generator.");


  static Parameter<MatchboxPhasespace,Energy> interfaceSingularCutoff
    ("SingularCutoff",
     "[debug] Cutoff below which a region is considered singular.",
     &MatchboxPhasespace::singularCutoff, GeV, 10.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  interfaceSingularCutoff.rank(-1);

  /*
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
  */

  static Command<MatchboxPhasespace> interfaceSetCoupling
    ("SetCoupling",
     "",
     &MatchboxPhasespace::doSetCoupling, false);

  static Command<MatchboxPhasespace> interfaceSetPhysicalCoupling
    ("SetPhysicalCoupling",
     "",
     &MatchboxPhasespace::doSetPhysicalCoupling, false);

  static Parameter<MatchboxPhasespace,int> interfaceLoopParticleIdMin
    ("LoopParticleIdMin",
     "First id in a range of id's meant to denote fictitious "
     "'ghost' particles to be used by the diagram generator "
     "in loop induced processes.",
     &MatchboxPhasespace::theLoopParticleIdMin, 200001, 0, 0,
     false, false, Interface::lowerlim);
  interfaceLoopParticleIdMin.rank(-1);

  static Parameter<MatchboxPhasespace,int> interfaceLoopParticleIdMax
    ("LoopParticleIdMax",
     "Last id in a range of id's meant to denote fictitious "
     "'ghost' particles to be used by the diagram generator "
     "in loop induced processes.",
     &MatchboxPhasespace::theLoopParticleIdMax, 200100, 0, 0,
     false, false, Interface::lowerlim);
  interfaceLoopParticleIdMax.rank(-1);


  static Reference<MatchboxPhasespace,PhasespaceCouplings> interfaceCouplingData
    ("CouplingData",
     "Set the storage for the couplings.",
     &MatchboxPhasespace::theCouplings, false, false, true, false, false);
  interfaceCouplingData.rank(-1);

}

