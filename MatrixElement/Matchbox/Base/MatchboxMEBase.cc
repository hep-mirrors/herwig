// -*- C++ -*-
//
// MatchboxMEBase.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMEBase class.
//

#include "MatchboxMEBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDF/PDF.h"
#include "ThePEG/PDT/PDT.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Handlers/StdXCombGroup.h"
#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/DiagramDrawer.h"
#include "Herwig++/MatrixElement/Matchbox/MatchboxFactory.h"

#include <iterator>
using std::ostream_iterator;

using namespace Herwig;

MatchboxMEBase::MatchboxMEBase() 
  : MEBase(), 
    theOneLoop(false),
    theOneLoopNoBorn(false) {}

MatchboxMEBase::~MatchboxMEBase() {}

Ptr<MatchboxFactory>::tcptr MatchboxMEBase::factory() const { return theFactory; }

void MatchboxMEBase::factory(Ptr<MatchboxFactory>::tcptr f) { theFactory = f; }

Ptr<Tree2toNGenerator>::tptr MatchboxMEBase::diagramGenerator() const { return factory()->diagramGenerator(); }

Ptr<ProcessData>::tptr MatchboxMEBase::processData() const { return factory()->processData(); }

unsigned int MatchboxMEBase::getNLight() const { return factory()->nLight(); }

double MatchboxMEBase::factorizationScaleFactor() const { return factory()->factorizationScaleFactor(); }

double MatchboxMEBase::renormalizationScaleFactor() const { return factory()->renormalizationScaleFactor(); }

bool MatchboxMEBase::fixedCouplings() const { return factory()->fixedCouplings(); }

bool MatchboxMEBase::fixedQEDCouplings() const { return factory()->fixedQEDCouplings(); }

bool MatchboxMEBase::checkPoles() const { return factory()->checkPoles(); }

bool MatchboxMEBase::verbose() const { return factory()->verbose(); }

bool MatchboxMEBase::initVerbose() const { return factory()->initVerbose(); }

void MatchboxMEBase::getDiagrams() const {

  if ( diagramGenerator() && processData() ) {

    vector<Ptr<Tree2toNDiagram>::ptr> diags;

    vector<Ptr<Tree2toNDiagram>::ptr>& res =
      processData()->diagramMap()[subProcess().legs];
    if ( res.empty() ) {
      res = diagramGenerator()->generate(subProcess().legs,orderInAlphaS(),orderInAlphaEW());
    }
    copy(res.begin(),res.end(),back_inserter(diags));
    processData()->fillMassGenerators(subProcess().legs);

    if ( diags.empty() )
      return;

    for ( vector<Ptr<Tree2toNDiagram>::ptr>::iterator d = diags.begin();
	  d != diags.end(); ++d ) {
      add(*d);
    }

    return;

  }

  throw Exception()
    << "MatchboxMEBase::getDiagrams() expects a Tree2toNGenerator and ProcessData object.\n"
    << "Please check your setup." << Exception::abortnow;

}

Selector<MEBase::DiagramIndex> 
MatchboxMEBase::diagrams(const DiagramVector & diags) const {

  if ( phasespace() ) {
    return phasespace()->selectDiagrams(diags);
  }

  throw Exception()
    << "MatchboxMEBase::diagrams() expects a MatchboxPhasespace object.\n"
    << "Please check your setup." << Exception::abortnow;
  return Selector<MEBase::DiagramIndex>();

}

Selector<const ColourLines *>
MatchboxMEBase::colourGeometries(tcDiagPtr diag) const {

  if ( matchboxAmplitude() ) {
    if ( matchboxAmplitude()->haveColourFlows() ) {
      if ( matchboxAmplitude()->treeAmplitudes() )
	matchboxAmplitude()->prepareAmplitudes(this);
      return matchboxAmplitude()->colourGeometries(diag);
    }
  }

  Ptr<Tree2toNDiagram>::tcptr tdiag =
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::tcptr>(diag);
  assert(diag && processData());

  vector<ColourLines*>& flows = processData()->colourFlowMap()[tdiag];

  if ( flows.empty() ) {

    list<list<list<pair<int,bool> > > > cflows =
      ColourBasis::colourFlows(tdiag);

    for ( list<list<list<pair<int,bool> > > >::const_iterator fit =
	    cflows.begin(); fit != cflows.end(); ++fit ) {
      flows.push_back(new ColourLines(ColourBasis::cfstring(*fit)));
    }

  }

  Selector<const ColourLines *> res;
  for ( vector<ColourLines*>::const_iterator f = flows.begin();
	f != flows.end(); ++f )
    res.insert(1.0,*f);

  return res;

}

unsigned int MatchboxMEBase::orderInAlphaS() const {
  return subProcess().orderInAlphaS;
}

unsigned int MatchboxMEBase::orderInAlphaEW() const {
  return subProcess().orderInAlphaEW;
}

void MatchboxMEBase::setXComb(tStdXCombPtr xc) {

  MEBase::setXComb(xc);
  lastMatchboxXComb(xc);
  if ( phasespace() )
    phasespace()->setXComb(xc);
  if ( scaleChoice() )
    scaleChoice()->setXComb(xc);
  if ( matchboxAmplitude() )
    matchboxAmplitude()->setXComb(xc);

}

double MatchboxMEBase::generateIncomingPartons(const double* r1, const double* r2) {
  
  // shamelessly stolen from PartonExtractor.cc

  Energy2 shmax = lastCuts().sHatMax();
  Energy2 shmin = lastCuts().sHatMin();
  Energy2 sh = shmin*pow(shmax/shmin, *r1);
  double ymax = lastCuts().yHatMax();
  double ymin = lastCuts().yHatMin();
  double km = log(shmax/shmin);
  ymax = min(ymax, log(lastCuts().x1Max()*sqrt(lastS()/sh)));
  ymin = max(ymin, -log(lastCuts().x2Max()*sqrt(lastS()/sh)));

  double y = ymin + (*r2)*(ymax - ymin);
  double x1 = exp(-0.5*log(lastS()/sh) + y);
  double x2 = exp(-0.5*log(lastS()/sh) - y);

  Lorentz5Momentum P1 = lastParticles().first->momentum();
  LorentzMomentum p1 = lightCone((P1.rho() + P1.e())*x1, Energy());
  p1.rotateY(P1.theta());
  p1.rotateZ(P1.phi());
  meMomenta()[0] = p1;

  Lorentz5Momentum P2 = lastParticles().second->momentum();
  LorentzMomentum p2 = lightCone((P2.rho() + P2.e())*x2, Energy());
  p2.rotateY(P2.theta());
  p2.rotateZ(P2.phi());
  meMomenta()[1] = p2;

  lastXCombPtr()->lastX1X2(make_pair(x1,x2));
  lastXCombPtr()->lastSHat((meMomenta()[0]+meMomenta()[1]).m2());

  return km*(ymax - ymin);

}

bool MatchboxMEBase::generateKinematics(const double * r) {

  if ( phasespace() ) {

    jacobian(phasespace()->generateKinematics(r,meMomenta()));
    if ( jacobian() == 0.0 )
      return false;

    setScale();
    logGenerateKinematics(r);

    assert(lastMatchboxXComb());

    if ( nDimAmplitude() > 0 ) {
      amplitudeRandomNumbers().resize(nDimAmplitude());
      copy(r + nDimPhasespace(), 
	   r + nDimPhasespace() + nDimAmplitude(),
	   amplitudeRandomNumbers().begin());
    }

    if ( nDimInsertions() > 0 ) {
      insertionRandomNumbers().resize(nDimInsertions());
      copy(r + nDimPhasespace() + nDimAmplitude(), 
	   r + nDimPhasespace() + nDimAmplitude() + nDimInsertions(),
	   insertionRandomNumbers().begin());
    }

    return true;

  }

  throw Exception()
    << "MatchboxMEBase::generateKinematics() expects a MatchboxPhasespace object.\n"
    << "Please check your setup." << Exception::abortnow;

  return false;

}

int MatchboxMEBase::nDim() const { 

  if ( lastMatchboxXComb() )
    return nDimPhasespace() + nDimAmplitude() + nDimInsertions();

  int ampAdd = 0;
  if ( matchboxAmplitude() ) {
    ampAdd = matchboxAmplitude()->nDimAdditional();
  }

  int insertionAdd = 0;
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v ) {
    insertionAdd = max(insertionAdd,(**v).nDimAdditional());
  }

  return nDimBorn() + ampAdd + insertionAdd;

}

int MatchboxMEBase::nDimBorn() const { 

  if ( lastMatchboxXComb() )
    return nDimPhasespace();

  if ( phasespace() ) {
    size_t nout = diagrams().front()->partons().size()-2;
    int n = phasespace()->nDim(nout);
    if ( phasespace()->useMassGenerators() ) {
      for ( cPDVector::const_iterator pd = 
	      diagrams().front()->partons().begin();
	    pd != diagrams().front()->partons().end(); ++pd ) {
	if ( processData()->massGenerator(*pd) ||
	     (**pd).width() != ZERO ) {
	  ++n;
	}
      }
    }
    return n;
  }

  throw Exception()
    << "MatchboxMEBase::nDim() expects a MatchboxPhasespace object.\n"
    << "Please check your setup." << Exception::abortnow;

  return 0;

}

void MatchboxMEBase::setScale() const {
  if ( haveX1X2() ) {
    lastXCombPtr()->lastSHat((meMomenta()[0]+meMomenta()[1]).m2());
  }
  Energy2 fscale = factorizationScale()*sqr(factorizationScaleFactor());
  Energy2 rscale = renormalizationScale()*sqr(renormalizationScaleFactor());
  Energy2 ewrscale = renormalizationScaleQED();
  lastXCombPtr()->lastScale(fscale);
  if ( !fixedCouplings() ) {
    if ( rscale > lastCuts().scaleMin() )
      lastXCombPtr()->lastAlphaS(SM().alphaS(rscale));
    else
      lastXCombPtr()->lastAlphaS(SM().alphaS(lastCuts().scaleMin()));
  } else {
    lastXCombPtr()->lastAlphaS(SM().alphaS());
  }
  if ( !fixedQEDCouplings() ) {
    lastXCombPtr()->lastAlphaEM(SM().alphaEM(ewrscale));
  } else {
    lastXCombPtr()->lastAlphaEM(SM().alphaEMMZ());
  }
  logSetScale();
}

Energy2 MatchboxMEBase::factorizationScale() const {
  if ( scaleChoice() ) {
    return scaleChoice()->factorizationScale();
  }

  throw Exception()
    << "MatchboxMEBase::factorizationScale() expects a MatchboxScaleChoice object.\n"
    << "Please check your setup." << Exception::abortnow;

  return ZERO;

}

Energy2 MatchboxMEBase::renormalizationScale() const {
  if ( scaleChoice() ) {
    return scaleChoice()->renormalizationScale();
  }

  throw Exception()
    << "MatchboxMEBase::renormalizationScale() expects a MatchboxScaleChoice object.\n"
    << "Please check your setup." << Exception::abortnow;

  return ZERO;

}

Energy2 MatchboxMEBase::renormalizationScaleQED() const {
  if ( scaleChoice() ) {
    return scaleChoice()->renormalizationScaleQED();
  }
  return renormalizationScale();
}

void MatchboxMEBase::setVetoScales(tSubProPtr) const {}

void MatchboxMEBase::getPDFWeight(Energy2 factorizationScale) const {

  if ( !mePartonData()[0]->coloured() &&
       !mePartonData()[1]->coloured() ) {
    lastMEPDFWeight(1.0);
    logPDFWeight();
    return;
  }

  double w = 1.;

  if ( mePartonData()[0]->coloured() && havePDFWeight1() )
    w *= pdf1(factorizationScale);

  if ( mePartonData()[1]->coloured() && havePDFWeight2() )
    w *= pdf2(factorizationScale);

  lastMEPDFWeight(w);

  logPDFWeight();

}

double MatchboxMEBase::pdf1(Energy2 fscale, double xEx) const {

  assert(lastXCombPtr()->partonBins().first->pdf());

  if ( xEx < 1. && lastX1() >= xEx ) {
    return
      ( ( 1. - lastX1() ) / ( 1. - xEx ) ) *
      lastXCombPtr()->partonBins().first->pdf()->xfx(lastParticles().first->dataPtr(),
						     lastPartons().first->dataPtr(),
						     fscale == ZERO ? lastScale() : fscale,
						     xEx)/xEx;
  }

  return lastXCombPtr()->partonBins().first->pdf()->xfx(lastParticles().first->dataPtr(),
							lastPartons().first->dataPtr(),
							fscale == ZERO ? lastScale() : fscale,
							lastX1())/lastX1();
}

double MatchboxMEBase::pdf2(Energy2 fscale, double xEx) const {

  assert(lastXCombPtr()->partonBins().second->pdf());

  if ( xEx < 1. && lastX2() >= xEx ) {
    return
      ( ( 1. - lastX2() ) / ( 1. - xEx ) ) *
      lastXCombPtr()->partonBins().second->pdf()->xfx(lastParticles().second->dataPtr(),
						      lastPartons().second->dataPtr(),
						      fscale == ZERO ? lastScale() : fscale,
						      xEx)/xEx;
  }

  return lastXCombPtr()->partonBins().second->pdf()->xfx(lastParticles().second->dataPtr(),
							 lastPartons().second->dataPtr(),
							 fscale == ZERO ? lastScale() : fscale,
							 lastX2())/lastX2();

}

double MatchboxMEBase::me2() const {

  if ( matchboxAmplitude() ) {

    if ( matchboxAmplitude()->treeAmplitudes() )
      matchboxAmplitude()->prepareAmplitudes(this);

    double res = 
      matchboxAmplitude()->me2()*
      crossingSign()*
      me2Norm();

    return res;

  }

  throw Exception()
    << "MatchboxMEBase::me2() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return 0.;

}

double MatchboxMEBase::finalStateSymmetry() const {

  if ( symmetryFactor() > 0.0 )
    return symmetryFactor();

  double sFactor = 1.;

  map<long,int> counts;

  cPDVector checkData;
  copy(mePartonData().begin()+2,mePartonData().end(),back_inserter(checkData));

  cPDVector::iterator p = checkData.begin();
  while ( !checkData.empty() ) {
    if ( counts.find((**p).id()) != counts.end() ) {
      counts[(**p).id()] += 1;
    } else {
      counts[(**p).id()] = 1;
    }
    checkData.erase(p);
    p = checkData.begin();
    continue;
  }

  for ( map<long,int>::const_iterator c = counts.begin();
	c != counts.end(); ++c ) {
    if ( c->second == 1 )
      continue;
    if ( c->second == 2 )
      sFactor /= 2.;
    else if ( c->second == 3 )
      sFactor /= 6.;
    else if ( c->second == 4 )
      sFactor /= 24.;
  }

  symmetryFactor(sFactor);

  return symmetryFactor();

}

double MatchboxMEBase::me2Norm(unsigned int addAlphaS) const {

  // assume that we always have incoming
  // spin-1/2 or massless spin-1 particles
  double fac = 1./4.;

  if ( hasInitialAverage() )
    fac = 1.;

  double couplings = 1.0;
  if ( orderInAlphaS() > 0 || addAlphaS != 0 ) {
    fac *= pow(lastAlphaS()/SM().alphaS(),double(orderInAlphaS()+addAlphaS));
    couplings *= pow(lastAlphaS(),double(orderInAlphaS()+addAlphaS));
  }
  if ( orderInAlphaEW() > 0 ) {
    fac *= pow(lastAlphaEM()/SM().alphaEM(),double(orderInAlphaEW()));
    couplings *= pow(lastAlphaEM(),double(orderInAlphaEW()));
  }
  lastMECouplings(couplings);

  if ( !hasInitialAverage() ) {
    if ( mePartonData()[0]->iColour() == PDT::Colour3 || 
	 mePartonData()[0]->iColour() == PDT::Colour3bar )
      fac /= SM().Nc();
    else if ( mePartonData()[0]->iColour() == PDT::Colour8 )
      fac /= (SM().Nc()*SM().Nc()-1.);

    if ( mePartonData()[1]->iColour() == PDT::Colour3 || 
	 mePartonData()[1]->iColour() == PDT::Colour3bar )
      fac /= SM().Nc();
    else if ( mePartonData()[1]->iColour() == PDT::Colour8 )
      fac /= (SM().Nc()*SM().Nc()-1.);
  }

  return !hasFinalStateSymmetry() ? finalStateSymmetry()*fac : fac;

}

CrossSection MatchboxMEBase::dSigHatDR() const {

  getPDFWeight();

  if ( !lastXCombPtr()->willPassCuts() ) {
    lastME2(0.0);
    lastMECrossSection(ZERO);
    return lastMECrossSection();
  }

  double xme2 = me2();
  lastME2(xme2);

  if ( xme2 == 0. && !oneLoopNoBorn() ) {
    lastMECrossSection(ZERO);
    return lastMECrossSection();
  }

  double vme2 = 0.;
  if ( oneLoop() )
    vme2 = oneLoopInterference();

  CrossSection res = ZERO;

  if ( !oneLoopNoBorn() )
    res += 
      (sqr(hbarc)/(2.*lastSHat())) *
      jacobian()* lastMEPDFWeight() * xme2;

  if ( oneLoop() )
    res += 
      (sqr(hbarc)/(2.*lastSHat())) *
      jacobian()* lastMEPDFWeight() * vme2;

  if ( !onlyOneLoop() ) {
    for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	    virtuals().begin(); v != virtuals().end(); ++v ) {
      (**v).setXComb(lastXCombPtr());
      res += (**v).dSigHatDR();
    }
    if ( checkPoles() )
      logPoles();
  }

  double weight = 0.0;
  bool applied = false;
  for ( vector<Ptr<MatchboxReweightBase>::ptr>::const_iterator rw =
	  theReweights.begin(); rw != theReweights.end(); ++rw ) {
    (**rw).setXComb(lastXCombPtr());
    if ( !(**rw).apply() )
      continue;
    weight += (**rw).evaluate();
    applied = true;
  }
  if ( applied )
    res *= weight;
  lastMECrossSection(res);

  return lastMECrossSection();

}

double MatchboxMEBase::oneLoopInterference() const {

  if ( matchboxAmplitude() ) {

    if ( matchboxAmplitude()->oneLoopAmplitudes() )
      matchboxAmplitude()->prepareOneLoopAmplitudes(this);

    double res = 
      matchboxAmplitude()->oneLoopInterference()*
      crossingSign()*
      me2Norm(1);

    return res;

  }

  throw Exception()
    << "MatchboxMEBase::oneLoopInterference() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return 0.;

}

MatchboxMEBase::AccuracyHistogram::AccuracyHistogram(double low,
						     double up,
						     unsigned int nbins) 
  : lower(low), upper(up), 
    sameSign(0), oppositeSign(0), nans(0),
    overflow(0), underflow(0) {

  double step = (up-low)/nbins;

  for ( unsigned int k = 1; k <= nbins; ++k )
    bins[lower + k*step] = 0.0;

}

void MatchboxMEBase::AccuracyHistogram::book(double a, double b) {
  if ( isnan(a) || isnan(b) ||
       isinf(a) || isinf(b) ) {
    ++nans;
    return;
  }
  if ( a*b >= 0. )
    ++sameSign;
  if ( a*b < 0. )
    ++oppositeSign;
  double r = 1.;
  if ( abs(a) != 0.0 )
    r = abs(1.-abs(b/a));
  else if ( abs(b) != 0.0 )
    r = abs(b);
  if ( log(r) < lower || r == 0.0 ) {
    ++underflow;
    return;
  }
  if ( log(r) > upper ) {
    ++overflow;
    return;
  }
  map<double,double>::iterator bin =
    bins.upper_bound(log(r));
  if ( bin == bins.end() )
    return;
  bin->second += 1.;
}

void MatchboxMEBase::AccuracyHistogram::dump(const std::string& prefix, 
					     const cPDVector& proc) const {
  ostringstream fname("");
  for ( cPDVector::const_iterator p = proc.begin();
	p != proc.end(); ++p )
    fname << (**p).PDGName();
  ofstream out((prefix+fname.str()+".dat").c_str());
  out << "# same sign : " << sameSign << " opposite sign : "
      << oppositeSign << " nans : " << nans 
      << " overflow : " << overflow
      << " underflow : " << underflow << "\n";
  for ( map<double,double>::const_iterator b = bins.begin();
	b != bins.end(); ++b ) {
    map<double,double>::const_iterator bp = b; --bp;
    if ( b->second != 0. ) {
      if ( b != bins.begin() )
	out << bp->first;
      else
	out << lower;
      out << " " << b->first
	  << " " << b->second
	  << "\n" << flush;
    }
  }
}

void MatchboxMEBase::AccuracyHistogram::persistentOutput(PersistentOStream& os) const {
  os << lower << upper << bins
     << sameSign << oppositeSign << nans
     << overflow << underflow;
}

void MatchboxMEBase::AccuracyHistogram::persistentInput(PersistentIStream& is) {
  is >> lower >> upper >> bins
     >> sameSign >> oppositeSign >> nans
     >> overflow >> underflow;
}

void MatchboxMEBase::logPoles() const {
  double res2me = oneLoopDoublePole();
  double res1me = oneLoopSinglePole();
  double res2i = 0.;
  double res1i = 0.;
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v ) {
    res2i += (**v).oneLoopDoublePole();
    res1i += (**v).oneLoopSinglePole();
  }
  epsilonSquarePoleHistograms[mePartonData()].book(res2me,res2i);
  epsilonPoleHistograms[mePartonData()].book(res1me,res1i);
}

bool MatchboxMEBase::haveOneLoop() const {
  if ( matchboxAmplitude() )
    return matchboxAmplitude()->haveOneLoop();
  return false;
}

bool MatchboxMEBase::onlyOneLoop() const {
  if ( matchboxAmplitude() )
    return matchboxAmplitude()->onlyOneLoop();
  return false;
}

bool MatchboxMEBase::isDR() const {
  if ( matchboxAmplitude() )
    return matchboxAmplitude()->isDR();
  return false;
}

bool MatchboxMEBase::isCS() const {
  if ( matchboxAmplitude() )
    return matchboxAmplitude()->isCS();
  return false;
}

bool MatchboxMEBase::isBDK() const {
  if ( matchboxAmplitude() )
    return matchboxAmplitude()->isBDK();
  return false;
}

bool MatchboxMEBase::isExpanded() const {
  if ( matchboxAmplitude() )
    return matchboxAmplitude()->isExpanded();
  return false;
}

Energy2 MatchboxMEBase::mu2() const {
  if ( matchboxAmplitude() )
    return matchboxAmplitude()->mu2();
  return 0*GeV2;
}

double MatchboxMEBase::oneLoopDoublePole() const {

  if ( matchboxAmplitude() ) {

    return
      matchboxAmplitude()->oneLoopDoublePole()*
      crossingSign()*
      me2Norm(1);

  }

  return 0.;

}

double MatchboxMEBase::oneLoopSinglePole() const {

  if ( matchboxAmplitude() ) {

    return 
      matchboxAmplitude()->oneLoopSinglePole()*
      crossingSign()*
      me2Norm(1);

  }

  return 0.;

}

vector<Ptr<SubtractionDipole>::ptr> 
MatchboxMEBase::getDipoles(const vector<Ptr<SubtractionDipole>::ptr>& dipoles,
			   const vector<Ptr<MatchboxMEBase>::ptr>& borns) const {

  vector<Ptr<SubtractionDipole>::ptr> res;

  // keep track of the dipoles we already did set up
  set<pair<pair<pair<int,int>,int>,pair<Ptr<MatchboxMEBase>::tptr,Ptr<SubtractionDipole>::tptr> > > done;

  cPDVector rep = diagrams().front()->partons();
  int nreal = rep.size();

  // now loop over configs
  for ( int emitter = 0; emitter < nreal; ++emitter ) {
    for ( int spectator = 0; spectator < nreal; ++spectator ) {
      if ( emitter == spectator )
	continue;
      for ( int emission = 2; emission < nreal; ++emission ) {
	if ( emission == emitter || emission == spectator )
	  continue;
	for ( vector<Ptr<MatchboxMEBase>::ptr>::const_iterator b =
		borns.begin(); b != borns.end(); ++b ) {
	  if ( (**b).onlyOneLoop() )
	    continue;
	  for ( vector<Ptr<SubtractionDipole>::ptr>::const_iterator d =
		  dipoles.begin(); d != dipoles.end(); ++d ) {
	    if ( !rep[emitter]->coloured() ||
		 !rep[emission]->coloured() ||
		 !rep[spectator]->coloured() ) {
	      continue;
	    }
	    if ( noDipole(emitter,emission,spectator) ) {
	      continue;
	    }
	    if ( done.find(make_pair(make_pair(make_pair(emitter,emission),spectator),make_pair(*b,*d))) 
		 != done.end() ) {
	      continue;
	    }
	    if ( !(**d).canHandle(rep,emitter,emission,spectator) ) {
	      continue;
	    }
	    // now get to work
	    (**d).clearBookkeeping();
	    (**d).realEmitter(emitter);
	    (**d).realEmission(emission);
	    (**d).realSpectator(spectator);
	    (**d).realEmissionME(const_cast<MatchboxMEBase*>(this));
	    (**d).underlyingBornME(*b);
	    (**d).setupBookkeeping();
	    if ( !((**d).empty()) ) {
	      Ptr<SubtractionDipole>::ptr nDipole = (**d).cloneMe();
	      res.push_back(nDipole);
	      done.insert(make_pair(make_pair(make_pair(emitter,emission),spectator),make_pair(*b,*d)));
	      if ( nDipole->isSymmetric() )
		done.insert(make_pair(make_pair(make_pair(emission,emitter),spectator),make_pair(*b,*d)));
	      ostringstream dname;
	      dname << fullName() << "." << (**b).name() << "."
		    << (**d).name() << ".[(" 
		    << emitter << "," << emission << ")," << spectator << "]";
	      if ( ! (generator()->preinitRegister(nDipole,dname.str()) ) )
		throw InitException() << "Dipole " << dname.str() << " already existing.";
	      if ( !factory()->reweighters().empty() ) {
		for ( vector<ReweightPtr>::const_iterator rw = factory()->reweighters().begin();
		      rw != factory()->reweighters().end(); ++rw )
		  nDipole->addReweighter(*rw);
	      }
	      if ( !factory()->preweighters().empty() ) {
		for ( vector<ReweightPtr>::const_iterator rw = factory()->preweighters().begin();
		      rw != factory()->preweighters().end(); ++rw )
		  nDipole->addPreweighter(*rw);
	      }
	      nDipole->cloneDependencies(dname.str());
	    }
	  }
	}
      }
    }
  }

  for ( vector<Ptr<SubtractionDipole>::ptr>::iterator d = res.begin();
	d != res.end(); ++d )
    (**d).partnerDipoles(res);

  return res;

}

double MatchboxMEBase::colourCorrelatedME2(pair<int,int> ij) const {

  if ( matchboxAmplitude() ) {

    if ( matchboxAmplitude()->treeAmplitudes() )
      matchboxAmplitude()->prepareAmplitudes(this);

    double res = 
      matchboxAmplitude()->colourCorrelatedME2(ij)*
      crossingSign()*
      me2Norm();

    return res;

  }

  throw Exception()
    << "MatchboxMEBase::colourCorrelatedME2() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return 0.;

}

double MatchboxMEBase::largeNColourCorrelatedME2(pair<int,int> ij,
						 Ptr<ColourBasis>::tptr largeNBasis) const {

  if ( matchboxAmplitude() ) {

    if ( matchboxAmplitude()->treeAmplitudes() )
      matchboxAmplitude()->prepareAmplitudes(this);
    largeNBasis->prepare(mePartonData(),false);

    double res = 
      matchboxAmplitude()->largeNColourCorrelatedME2(ij,largeNBasis)*
      crossingSign()*
      me2Norm();

    return res;

  }

  throw Exception()
    << "MatchboxMEBase::largeNColourCorrelatedME2() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return 0.;

}

double MatchboxMEBase::spinColourCorrelatedME2(pair<int,int> ij,
					       const SpinCorrelationTensor& c) const {

  if ( matchboxAmplitude() ) {

    if ( matchboxAmplitude()->treeAmplitudes() )
      matchboxAmplitude()->prepareAmplitudes(this);

    double res = 
      matchboxAmplitude()->spinColourCorrelatedME2(ij,c)*
      crossingSign()*
      me2Norm();

    return res;

  }

  throw Exception()
    << "MatchboxMEBase::spinColourCorrelatedME2() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return 0.;

}

void MatchboxMEBase::flushCaches() { 
  MEBase::flushCaches();
  if ( matchboxAmplitude() )
    matchboxAmplitude()->flushCaches();
  for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator r =
	  reweights().begin(); r != reweights().end(); ++r ) {
    (**r).flushCaches();
  }
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v ) {
    (**v).flushCaches();
  }
}

void MatchboxMEBase::print(ostream& os) const {

  os << "--- MatchboxMEBase setup -------------------------------------------------------\n";

  os << " '" << name() << "' for subprocess:\n";

  os << "  ";
  for ( PDVector::const_iterator pp = subProcess().legs.begin();
	pp != subProcess().legs.end(); ++pp ) {
    os << (**pp).PDGName() << " ";
    if ( pp == subProcess().legs.begin() + 1 )
      os << "-> ";
  }
  os << "\n";

  os << " including " << (oneLoop() ? "" : "no ") << "virtual corrections";
  if ( oneLoopNoBorn() )
    os << " without Born contributions";
  os << "\n";

  if ( oneLoop() && !onlyOneLoop() ) {
    os << " using insertion operators\n";
    for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	    virtuals().begin(); v != virtuals().end(); ++v ) {
      os << " '" << (**v).name() << "' with " 
	 << ((**v).isDR() ? "" : "C") << "DR/";
      if ( (**v).isCS() )
	os << "CS";
      if ( (**v).isBDK() )
	os << "BDK";
      if ( (**v).isExpanded() )
	os << "expanded";
      os << " conventions\n";
    }
  }

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void MatchboxMEBase::printLastEvent(ostream& os) const {

  os << "--- MatchboxMEBase last event information --------------------------------------\n";

  os << " for matrix element '" << name() << "'\n";

  os << " process considered:\n ";

  int in = 0;
  for ( cPDVector::const_iterator p = mePartonData().begin();
	p != mePartonData().end(); ++p ) {
    os << (**p).PDGName() << " ";
    if ( ++in == 2 )
      os << " -> ";
  }

  os << " kinematic environment as set by the XComb " << lastXCombPtr() << ":\n"
     << " sqrt(shat)/GeV = " << sqrt(lastSHat()/GeV2)
     << " x1 = " << lastX1() << " x2 = " << lastX2() 
     << " alphaS = " << lastAlphaS() << "\n";

  os << " momenta/GeV generated from random numbers\n ";
  copy(lastXComb().lastRandomNumbers().begin(),
       lastXComb().lastRandomNumbers().end(),ostream_iterator<double>(os," "));
  os << ":\n ";

  for ( vector<Lorentz5Momentum>::const_iterator p = meMomenta().begin();
	p != meMomenta().end(); ++p ) {
    os << (*p/GeV) << "\n ";
  }

  os << "last cross section/nb calculated was:\n "
     << (lastMECrossSection()/nanobarn) << " (pdf weight " << lastMEPDFWeight() << ")\n";

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void MatchboxMEBase::logGenerateKinematics(const double * r) const {

  if ( !verbose() )
    return;

  generator()->log() << "'" << name() << "' generated kinematics\nfrom "
		     << nDim() << " random numbers:\n";
  copy(r,r+nDim(),ostream_iterator<double>(generator()->log()," "));
  generator()->log() << "\n";

  generator()->log() << "storing phase space information in XComb "
		     << lastXCombPtr() << "\n";

  generator()->log() << "generated phase space point (in GeV):\n";

  vector<Lorentz5Momentum>::const_iterator pit = meMomenta().begin();
  cPDVector::const_iterator dit = mePartonData().begin();

  for ( ; pit != meMomenta().end() ; ++pit, ++dit )
    generator()->log() << (**dit).PDGName() << " : "
		       << (*pit/GeV) << "\n";

  generator()->log() << "with x1 = " << lastX1() << " x2 = " << lastX2() << "\n"
		     << "and Jacobian = " << jacobian() << " sHat/GeV2 = "
		     << (lastSHat()/GeV2) << "\n" << flush;

}

void MatchboxMEBase::logSetScale() const {

  if ( !verbose() )
    return;

  generator()->log() << "'" << name() << "' set scales using XComb " << lastXCombPtr() << ":\n"
		     << "scale/GeV2 = " << (scale()/GeV2) << " xi_R = "
		     << renormalizationScaleFactor() << " xi_F = "
		     << factorizationScaleFactor() << "\n"
		     << "alpha_s = " << lastAlphaS() << "\n" << flush;

}

void MatchboxMEBase::logPDFWeight() const {

  if ( !verbose() )
    return;

  generator()->log() << "'" << name() << "' calculated pdf weight = "
		     << lastMEPDFWeight() << " from XComb "
		     << lastXCombPtr() << "\n"
		     << "x1 = " << lastX1() << " (" << (mePartonData()[0]->coloured() ? "" : "not ") << "used) "
		     << "x2 = " << lastX2() << " (" << (mePartonData()[1]->coloured() ? "" : "not ") << "used)\n"
		     << flush;

}

void MatchboxMEBase::logME2() const {

  if ( !verbose() )
    return;

  generator()->log() << "'" << name() << "' evaluated me2 using XComb "
		     << lastXCombPtr() << "\n"
		     << "and phase space point (in GeV):\n";

  vector<Lorentz5Momentum>::const_iterator pit = meMomenta().begin();
  cPDVector::const_iterator dit = mePartonData().begin();

  for ( ; pit != meMomenta().end() ; ++pit, ++dit )
    generator()->log() << (**dit).PDGName() << " : "
		       << (*pit/GeV) << "\n";

  generator()->log() << "with x1 = " << lastX1() << " x2 = " << lastX2() << "\n"
		     << "sHat/GeV2 = " << (lastSHat()/GeV2) 
		     << " me2 = " << lastME2() << "\n" << flush;

}

void MatchboxMEBase::logDSigHatDR() const {

  if ( !verbose() )
    return;

  generator()->log() << "'" << name() << "' evaluated cross section using XComb "
		     << lastXCombPtr() << "\n"
		     << "Jacobian = " << jacobian() << " sHat/GeV2 = "
		     << (lastSHat()/GeV2) << " dsig/nb = "
		     << (lastMECrossSection()/nanobarn) << "\n" << flush;

}

void MatchboxMEBase::cloneDependencies(const std::string& prefix) {

  if ( phasespace() ) {
    Ptr<MatchboxPhasespace>::ptr myPhasespace = phasespace()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myPhasespace->name();
    if ( ! (generator()->preinitRegister(myPhasespace,pname.str()) ) )
      throw InitException() << "Phasespace generator " << pname.str() << " already existing.";
    myPhasespace->cloneDependencies(pname.str());
    phasespace(myPhasespace);
  }

  theAmplitude = dynamic_ptr_cast<Ptr<MatchboxAmplitude>::ptr>(amplitude());

  if ( matchboxAmplitude() ) {
    Ptr<MatchboxAmplitude>::ptr myAmplitude = matchboxAmplitude()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myAmplitude->name();
    if ( ! (generator()->preinitRegister(myAmplitude,pname.str()) ) )
      throw InitException() << "Amplitude " << pname.str() << " already existing.";
    myAmplitude->cloneDependencies(pname.str());
    matchboxAmplitude(myAmplitude);
    amplitude(myAmplitude);
    matchboxAmplitude()->orderInGs(orderInAlphaS());
    matchboxAmplitude()->orderInGem(orderInAlphaEW());
  }

  if ( scaleChoice() ) {
    Ptr<MatchboxScaleChoice>::ptr myScaleChoice = scaleChoice()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myScaleChoice->name();
    if ( ! (generator()->preinitRegister(myScaleChoice,pname.str()) ) )
      throw InitException() << "Scale choice " << pname.str() << " already existing.";
    scaleChoice(myScaleChoice);
  }

  for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator rw =
	  theReweights.begin(); rw != theReweights.end(); ++rw ) {
    Ptr<MatchboxReweightBase>::ptr myReweight = (**rw).cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << (**rw).name();
    if ( ! (generator()->preinitRegister(myReweight,pname.str()) ) )
      throw InitException() << "Reweight " << pname.str() << " already existing.";
    myReweight->cloneDependencies(pname.str());
    *rw = myReweight;
  }

  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v ) {
    Ptr<MatchboxInsertionOperator>::ptr myIOP = (**v).cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << (**v).name();
    if ( ! (generator()->preinitRegister(myIOP,pname.str()) ) )
      throw InitException() << "Insertion operator " << pname.str() << " already existing.";
    *v = myIOP;
  }

}

void MatchboxMEBase::prepareXComb(MatchboxXCombData& xc) const {

  if ( phasespace() ) {
    size_t nout = diagrams().front()->partons().size()-2;
    xc.nDimPhasespace(phasespace()->nDim(nout));
  }

  if ( matchboxAmplitude() ) {
    xc.nDimAmplitude(matchboxAmplitude()->nDimAdditional());
    if ( matchboxAmplitude()->colourBasis() ) {
      size_t cdim = 
 	matchboxAmplitude()->colourBasis()->prepare(diagrams(),matchboxAmplitude()->noCorrelations());
      xc.colourBasisDim(cdim);
    }
  }

  int insertionAdd = 0;
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v ) {
    insertionAdd = max(insertionAdd,(**v).nDimAdditional());
  }
  xc.nDimInsertions(insertionAdd);

  xc.nLight(getNLight());

  xc.olpId(olpProcess());

  if ( initVerbose() ) {
    string fname = name() + ".diagrams";
    ifstream test(fname.c_str());
    if ( !test ) {
      test.close();
      ofstream out(fname.c_str());
      for ( vector<Ptr<DiagramBase>::ptr>::const_iterator d = diagrams().begin();
	    d != diagrams().end(); ++d ) {
	DiagramDrawer::drawDiag(out,dynamic_cast<const Tree2toNDiagram&>(**d));
	out << "\n";
      }
    }
  }

}

StdXCombPtr MatchboxMEBase::makeXComb(Energy newMaxEnergy, const cPDPair & inc,
				      tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
				      tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
				      const PBPair & newPartonBins, tCutsPtr newCuts,
				      const DiagramVector & newDiagrams, bool mir,
				      const PartonPairVec&,
				      tStdXCombPtr newHead,
				      tMEPtr newME) {

  if ( !newME )
    newME = this;

  Ptr<MatchboxXComb>::ptr xc =
    new_ptr(MatchboxXComb(newMaxEnergy, inc,
			  newEventHandler, newSubProcessHandler,
			  newExtractor, newCKKW,
			  newPartonBins, newCuts, newME,
			  newDiagrams, mir,
			  newHead));

  prepareXComb(*xc);

  return xc;

}

StdXCombPtr MatchboxMEBase::makeXComb(tStdXCombPtr newHead,
				      const PBPair & newPartonBins,
				      const DiagramVector & newDiagrams,
				      tMEPtr newME) {
  if ( !newME )
    newME = this;

  Ptr<MatchboxXComb>::ptr xc = 
    new_ptr(MatchboxXComb(newHead, newPartonBins, newME, newDiagrams));

  prepareXComb(*xc);

  return xc;
}

void MatchboxMEBase::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb << theFactory << thePhasespace 
     << theAmplitude << theScaleChoice << theVirtuals 
     << theReweights << theSubprocess << theOneLoop 
     << theOneLoopNoBorn
     << epsilonSquarePoleHistograms << epsilonPoleHistograms
     << theOLPProcess;
}

void MatchboxMEBase::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb >> theFactory >> thePhasespace 
     >> theAmplitude >> theScaleChoice >> theVirtuals 
     >> theReweights >> theSubprocess >> theOneLoop 
     >> theOneLoopNoBorn
     >> epsilonSquarePoleHistograms >> epsilonPoleHistograms
     >> theOLPProcess;
  lastMatchboxXComb(theLastXComb);
}

void MatchboxMEBase::Init() {

  static ClassDocumentation<MatchboxMEBase> documentation
    ("MatchboxMEBase is the base class for matrix elements "
     "in the context of the matchbox NLO interface.");

}

IBPtr MatchboxMEBase::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMEBase::fullclone() const {
  return new_ptr(*this);
}

void MatchboxMEBase::doinit() {
  MEBase::doinit();
  if ( !theAmplitude )
    theAmplitude = dynamic_ptr_cast<Ptr<MatchboxAmplitude>::ptr>(amplitude());
}

void MatchboxMEBase::dofinish() {
  MEBase::dofinish();
  for ( map<cPDVector,AccuracyHistogram>::const_iterator
	  b = epsilonSquarePoleHistograms.begin();
	b != epsilonSquarePoleHistograms.end(); ++b ) {
    b->second.dump(factory()->poleData() + "epsilonSquarePoles-",b->first);
  }
  for ( map<cPDVector,AccuracyHistogram>::const_iterator
	  b = epsilonPoleHistograms.begin();
	b != epsilonPoleHistograms.end(); ++b ) {
    b->second.dump(factory()->poleData() + "epsilonPoles-",b->first);
  }
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxMEBase,MEBase>
describeHerwigMatchboxMEBase("Herwig::MatchboxMEBase", "HwMatchbox.so");

