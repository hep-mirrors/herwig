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

#include <iterator>
using std::ostream_iterator;

using namespace Herwig;

MatchboxMEBase::MatchboxMEBase() 
  : MEBase(), 
    theFactorizationScaleFactor(1.0),
    theRenormalizationScaleFactor(1.0),
    theVerbose(false),
    theFixedCouplings(false),
    theNLight(0), theGetColourCorrelatedMEs(false) {}

MatchboxMEBase::~MatchboxMEBase() {}

void MatchboxMEBase::getDiagrams() const {

  if ( diagramGenerator() ) {

    vector<Ptr<Tree2toNDiagram>::ptr> diags;
    for ( vector<PDVector>::const_iterator p = subProcesses().begin();
	  p != subProcesses().end(); ++p ) {
      vector<Ptr<Tree2toNDiagram>::ptr> res =
	diagramGenerator()->generate(*p,orderInAlphaS(),orderInAlphaEW());
      copy(res.begin(),res.end(),back_inserter(diags));
    }

    if ( diags.empty() )
      return;

    for ( vector<Ptr<Tree2toNDiagram>::ptr>::iterator d = diags.begin();
	  d != diags.end(); ++d ) {
      add(*d);
    }

    if ( theVerbose ) {
      string fname = name() + ".diagrams";
      ifstream test(fname.c_str());
      if ( !test ) {
	test.close();
	ofstream out(fname.c_str());
	for ( vector<Ptr<Tree2toNDiagram>::ptr>::const_iterator d = diags.begin();
	      d != diags.end(); ++d ) {
	  DiagramDrawer::drawDiag(out,**d);
	  out << "\n";
	}
      }
    }

    return;

  }

  throw Exception()
    << "MatchboxMEBase::getDiagrams() expects a Tree2toNGenerator object.\n"
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
    if ( !matchboxAmplitude()->haveColourFlows() )
      throw Exception() << "A colour flow implementation is not present."
			<< Exception::abortnow;
    matchboxAmplitude()->prepareAmplitudes();
    return matchboxAmplitude()->colourGeometries(diag);
  }

  throw Exception()
    << "MatchboxMEBase::colourGeometries() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return Selector<const ColourLines *>();

}

unsigned int MatchboxMEBase::orderInAlphaS() const {
  if ( matchboxAmplitude() ) {
    return matchboxAmplitude()->orderInGs();
  }
  throw Exception()
    << "MatchboxMEBase::orderInAlphaS() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return 0;
}

unsigned int MatchboxMEBase::orderInAlphaEW() const {
  if ( matchboxAmplitude() ) {
    return matchboxAmplitude()->orderInGem();
  }
  throw Exception()
    << "MatchboxMEBase::orderInAlphaEW() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return 0;
}

void MatchboxMEBase::setXComb(tStdXCombPtr xc) {

  MEBase::setXComb(xc);
  for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator rw =
	  theReweights.begin(); rw != theReweights.end(); ++rw )
    (**rw).setXComb(xc);
  if ( phasespace() && !xc->head() )
    phasespace()->prepare(xc,theVerbose);
  if ( scaleChoice() )
    scaleChoice()->setXComb(xc);
  if ( cache() )
    cache()->setXComb(xc);
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

    return true;

  }

  throw Exception()
    << "MatchboxMEBase::generateKinematics() expects a MatchboxPhasespace object.\n"
    << "Please check your setup." << Exception::abortnow;

  return false;

}

int MatchboxMEBase::nDim() const { 
  if ( phasespace() ) {
    size_t nout = diagrams().front()->partons().size()-2;
    int n = phasespace()->nDim(nout);
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
  Energy2 fscale = factorizationScale()*factorizationScaleFactor();
  Energy2 rscale = renormalizationScale()*renormalizationScaleFactor();
  lastXCombPtr()->lastScale(fscale);
  if ( !fixedCouplings() ) {
    if ( rscale > lastCuts().scaleMin() )
      lastXCombPtr()->lastAlphaS(SM().alphaS(rscale));
    else
      lastXCombPtr()->lastAlphaS(SM().alphaS(lastCuts().scaleMin()));
    lastXCombPtr()->lastAlphaEM(SM().alphaEM(rscale));
  } else {
    lastXCombPtr()->lastAlphaS(SM().alphaS());
    lastXCombPtr()->lastAlphaEM(SM().alphaEM());
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

void MatchboxMEBase::setVetoScales(tSubProPtr subpro) const {
  for ( vector<Ptr<MatchboxReweightBase>::ptr>::const_iterator rw =
	  theReweights.begin(); rw != theReweights.end(); ++rw ) {
    if ( !(**rw).apply() )
      continue;
    (**rw).setVetoScales(subpro);
  }
}

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

double MatchboxMEBase::pdf1(Energy2 fscale) const {

  assert(lastXCombPtr()->partonBins().first->pdf());

  return lastXCombPtr()->partonBins().first->pdf()->xfx(lastParticles().first->dataPtr(),
							lastPartons().first->dataPtr(),
							fscale == ZERO ? lastScale() : fscale,
							lastX1())/lastX1();
}

double MatchboxMEBase::pdf2(Energy2 fscale) const {

  assert(lastXCombPtr()->partonBins().second->pdf());

  return lastXCombPtr()->partonBins().second->pdf()->xfx(lastParticles().second->dataPtr(),
							 lastPartons().second->dataPtr(),
							 fscale == ZERO ? lastScale() : fscale,
							 lastX2())/lastX2();

}

double MatchboxMEBase::me2() const {

  if ( matchboxAmplitude() ) {

    double res;
    if ( !calculateME2(res) )
      return res;

    matchboxAmplitude()->prepareAmplitudes();

    lastME2(matchboxAmplitude()->me2()*
	    matchboxAmplitude()->lastCrossingSign()*
	    me2Norm());

    cacheME2(lastME2());

    logME2();
    
    return lastME2();

  }

  throw Exception()
    << "MatchboxMEBase::me2() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return 0.;

}

double MatchboxMEBase::finalStateSymmetry() const {

  map<tStdXCombPtr,double>::const_iterator s =
    symmetryFactors.find(lastXCombPtr());

  if ( s != symmetryFactors.end() )
    return s->second;

  double symmetryFactor = 1.;

  map<long,int> counts;

  cPDVector checkData;
  copy(mePartonData().begin()+2,mePartonData().end(),back_inserter(checkData));

  cPDVector::iterator p = checkData.begin();
  while ( !checkData.empty() ) {
    if ( (**p).iSpin() == PDT::Spin0 ||
	 (**p).iSpin() == PDT::Spin1 ) {
      if ( counts.find((**p).id()) != counts.end() ) {
	counts[(**p).id()] += 1;
      } else {
	counts[(**p).id()] = 1;
      }
      checkData.erase(p);
      p = checkData.begin();
      continue;
    } else if ( (**p).iSpin() == PDT::Spin1Half ) {
      assert((**p).CC());
      long fid = abs((**p).id());
      cPDPtr findCC = (**p).CC();
      checkData.erase(p);
      p = checkData.begin();
      for ( ; p != checkData.end(); ++p )
	if ( (*p) == findCC )
	  break;
      if ( p != checkData.end() ) {
	if ( counts.find(fid) != counts.end() ) {
	  counts[fid] += 1;
	} else {
	  counts[fid] = 1;
	}
	checkData.erase(p);
      }
      p = checkData.begin();
      continue;
    } else assert(false);
  }

  for ( map<long,int>::const_iterator c = counts.begin();
	c != counts.end(); ++c ) {
    if ( c->second == 1 )
      continue;
    if ( c->second == 2 )
      symmetryFactor /= 2.;
    else if ( c->second == 3 )
      symmetryFactor /= 6.;
    else if ( c->second == 4 )
      symmetryFactor /= 24.;
  }

  symmetryFactors[lastXCombPtr()] = symmetryFactor;

  return symmetryFactor;

}

double MatchboxMEBase::me2Norm(unsigned int addAlphaS) const {

  // assume that we always have incoming
  // spin-1/2 or massless spin-1 particles
  double fac = 1./4.;

  if ( orderInAlphaS() > 0 || addAlphaS != 0 )
    fac *= pow(lastAlphaS()/SM().alphaS(),double(orderInAlphaS()+addAlphaS));
  if ( orderInAlphaEW() > 0 )
    fac *= pow(lastAlphaEM()/SM().alphaEM(),double(orderInAlphaEW()));

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

  return finalStateSymmetry()*fac;

}

void MatchboxMEBase::storeColourCorrelatedMEs(double xme2) const {

  map<pair<int,int>,double>& ccs = 
    colourCorrelatedMEs[lastXCombPtr()];
  int n = meMomenta().size();
  for ( int i = 0; i < n; ++i ) {
    if ( !mePartonData()[i]->coloured() )
      continue;
    for ( int j = i+1; j < n; ++j ) {
      if ( !mePartonData()[j]->coloured() )
	continue;
      if ( noDipole(i,j) )
	continue;
      ccs[make_pair(i,j)] = colourCorrelatedME2(make_pair(i,j));
    }
  }
  lastXCombPtr()->meta(MatchboxMetaKeys::ColourCorrelatedMEs,ccs);
  double& bme2 = bornMEs[lastXCombPtr()];
  bme2 = xme2 >= 0. ? xme2 : me2();
  lastXCombPtr()->meta(MatchboxMetaKeys::BornME,bme2);

}

CrossSection MatchboxMEBase::dSigHatDR() const {

  getPDFWeight();

  if ( !lastXComb().willPassCuts() ) {
    lastMECrossSection(ZERO);
    return lastMECrossSection();
  }

  double xme2 = me2();

  if ( xme2 == 0. ) {
    lastMECrossSection(ZERO);
    return lastMECrossSection();
  }

  if ( getColourCorrelatedMEs() )
    storeColourCorrelatedMEs(xme2);

  CrossSection res = 
    (sqr(hbarc)/(2.*lastSHat())) *
    jacobian() * xme2;

  double weight = 0.0;
  bool applied = false;
  for ( vector<Ptr<MatchboxReweightBase>::ptr>::const_iterator rw =
	  theReweights.begin(); rw != theReweights.end(); ++rw ) {
    if ( !(**rw).apply() )
      continue;
    weight += (**rw).evaluate();
    applied = true;
  }
  if ( applied )
    res *= weight;
  lastMECrossSection(res * lastMEPDFWeight());

  return lastMECrossSection();

}

double MatchboxMEBase::oneLoopInterference() const {

  if ( matchboxAmplitude() ) {

    double res;
    if ( !calculateME2(res,make_pair<int,int>(-1,-1)) )
      return res;

    matchboxAmplitude()->prepareOneLoopAmplitudes();
    lastME2(matchboxAmplitude()->oneLoopInterference()*
	    matchboxAmplitude()->lastCrossingSign()*
	    me2Norm(1));

    cacheME2(lastME2(),make_pair<int,int>(-1,-1));

    logME2();
    
    return lastME2();

  }

  throw Exception()
    << "MatchboxMEBase::oneLoopInterference() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return 0.;

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
	    Ptr<SubtractionDipole>::ptr nDipole = (**d).cloneMe();
	    nDipole->realEmitter(emitter);
	    nDipole->realEmission(emission);
	    nDipole->realSpectator(spectator);
	    nDipole->realEmissionME(const_cast<MatchboxMEBase*>(this));
	    nDipole->underlyingBornME(*b);
	    nDipole->setupBookkeeping();
	    if ( !(nDipole->empty()) ) {
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
	      nDipole->cloneDependencies(dname.str());
	    }
	  }
	}
      }
    }
  }

  return res;

}

double MatchboxMEBase::colourCorrelatedME2(pair<int,int> ij) const {

  if ( matchboxAmplitude() ) {

    double res;
    if ( !calculateME2(res,ij) )
      return res;

    matchboxAmplitude()->prepareAmplitudes();
    lastME2(matchboxAmplitude()->colourCorrelatedME2(ij)*
	    matchboxAmplitude()->lastCrossingSign()*
	    me2Norm());

    cacheME2(lastME2(),ij);

    logME2();
    
    return lastME2();

  }

  throw Exception()
    << "MatchboxMEBase::colourCorrelatedME2() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return 0.;

}

double MatchboxMEBase::spinColourCorrelatedME2(pair<int,int> ij,
					       const SpinCorrelationTensor& c) const {

  if ( matchboxAmplitude() ) {

    double res;
    if ( !calculateME2(res,ij) )
      return res;

    matchboxAmplitude()->prepareAmplitudes();
    lastME2(matchboxAmplitude()->spinColourCorrelatedME2(ij,c)*
	    matchboxAmplitude()->lastCrossingSign()*
	    me2Norm());

    cacheME2(lastME2(),ij);

    logME2();
    
    return lastME2();

  }

  throw Exception()
    << "MatchboxMEBase::spinColourCorrelatedME2() expects a MatchboxAmplitude object.\n"
    << "Please check your setup." << Exception::abortnow;
  return 0.;

}

void MatchboxMEBase::flushCaches() { 
  if ( cache() )
    cache()->flush();
  if ( matchboxAmplitude() )
    matchboxAmplitude()->flushCaches();
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
  copy(meInfo().begin(),meInfo().end(),ostream_iterator<double>(os," "));
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

  if ( !theVerbose )
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

  if ( !theVerbose )
    return;

  generator()->log() << "'" << name() << "' set scales using XComb " << lastXCombPtr() << ":\n"
		     << "scale/GeV2 = " << (scale()/GeV2) << " xi_R = "
		     << renormalizationScaleFactor() << " xi_F = "
		     << factorizationScaleFactor() << "\n"
		     << "alpha_s = " << lastAlphaS() << "\n" << flush;

}

void MatchboxMEBase::logPDFWeight() const {

  if ( !theVerbose )
    return;

  generator()->log() << "'" << name() << "' calculated pdf weight = "
		     << lastMEPDFWeight() << " from XComb "
		     << lastXCombPtr() << "\n"
		     << "x1 = " << lastX1() << " (" << (mePartonData()[0]->coloured() ? "" : "not ") << "used) "
		     << "x2 = " << lastX2() << " (" << (mePartonData()[1]->coloured() ? "" : "not ") << "used)\n"
		     << flush;

}

void MatchboxMEBase::logME2() const {

  if ( !theVerbose )
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

  if ( !theVerbose )
    return;

  generator()->log() << "'" << name() << "' evaluated cross section using XComb "
		     << lastXCombPtr() << "\n"
		     << "Jacobian = " << jacobian() << " sHat/GeV2 = "
		     << (lastSHat()/GeV2) << " dsig/nb = "
		     << (lastMECrossSection()/nanobarn) << "\n" << flush;

}

void MatchboxMEBase::dumpInfo(const string& prefix) const {
  generator()->log() << prefix << fullName()
		     << " [" << this << "]\n";
  generator()->log() << prefix << "  | XComb " << lastXCombPtr()
		     << " for ";
  if ( lastXCombPtr() ) {
    for ( cPDVector::const_iterator p = lastXComb().mePartonData().begin();
	  p != lastXComb().mePartonData().end(); ++p ) {
      generator()->log() << (**p).PDGName() << " ";
    }
  }
  generator()->log() << "\n";
  if ( !reweights().empty() ) {
    generator()->log() << prefix << "  | Reweights\n";
    for ( vector<Ptr<MatchboxReweightBase>::ptr>::const_iterator r =
	    reweights().begin(); r != reweights().end(); ++r ) {
      (**r).dumpInfo(prefix + "  | ");
    }
  }
  if ( phasespace() ) {
    generator()->log() << prefix << "  | Phasespace\n";
    phasespace()->dumpInfo(prefix + "  | ");
  }
  if ( matchboxAmplitude() ) {
    generator()->log() << prefix << "  | Amplitude\n";
    matchboxAmplitude()->dumpInfo(prefix + "  | ");
  }
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

}

void MatchboxMEBase::persistentOutput(PersistentOStream & os) const {
  os << theReweights << thePhasespace << theAmplitude << theScaleChoice
     << theDiagramGenerator << theSubprocesses
     << theFactorizationScaleFactor << theRenormalizationScaleFactor
     << theVerbose << theCache << theFixedCouplings
     << theNLight << theGetColourCorrelatedMEs << symmetryFactors;
}

void MatchboxMEBase::persistentInput(PersistentIStream & is, int) {
  is >> theReweights >> thePhasespace >> theAmplitude >> theScaleChoice
     >> theDiagramGenerator >> theSubprocesses
     >> theFactorizationScaleFactor >> theRenormalizationScaleFactor
     >> theVerbose >> theCache >> theFixedCouplings
     >> theNLight >> theGetColourCorrelatedMEs >> symmetryFactors;
}

void MatchboxMEBase::Init() {

  static ClassDocumentation<MatchboxMEBase> documentation
    ("MatchboxMEBase is the base class for matrix elements "
     "in the context of the matchbox NLO interface.");

  static RefVector<MatchboxMEBase,MatchboxReweightBase> interfaceReweights
    ("Reweights",
     "Reweight objects to be applied to this matrix element.",
     &MatchboxMEBase::theReweights, -1, false, false, true, true, false);


  static Reference<MatchboxMEBase,MatchboxPhasespace> interfacePhasespace
    ("Phasespace",
     "Set the phasespace generator to be used.",
     &MatchboxMEBase::thePhasespace, false, false, true, true, false);

  static Reference<MatchboxMEBase,Tree2toNGenerator> interfaceDiagramGenerator
    ("DiagramGenerator",
     "Set the diagram generator to be used.",
     &MatchboxMEBase::theDiagramGenerator, false, false, true, true, false);

  static Reference<MatchboxMEBase,MatchboxScaleChoice> interfaceScaleChoice
    ("ScaleChoice",
     "Set the scale choice to be used.",
     &MatchboxMEBase::theScaleChoice, false, false, true, true, false);

  static Reference<MatchboxMEBase,MatchboxMECache> interfaceMECache
    ("Cache",
     "Set the cache object to be used.",
     &MatchboxMEBase::theCache, false, false, true, true, false);

  static Parameter<MatchboxMEBase,double> interfaceFactorizationScaleFactor
    ("FactorizationScaleFactor",
     "The factorization scale factor.",
     &MatchboxMEBase::theFactorizationScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxMEBase,double> interfaceRenormalizationScaleFactor
    ("RenormalizationScaleFactor",
     "The renormalization scale factor.",
     &MatchboxMEBase::theRenormalizationScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Switch<MatchboxMEBase,bool> interfaceVerbose
    ("Verbose",
     "Print full infomation on each evaluated phase space point.",
     &MatchboxMEBase::theVerbose, false, false, false);
  static SwitchOption interfaceVerboseOn
    (interfaceVerbose,
     "On",
     "On",
     true);
  static SwitchOption interfaceVerboseOff
    (interfaceVerbose,
     "Off",
     "Off",
     false);

  static Switch<MatchboxMEBase,bool> interfaceFixedCouplings
    ("FixedCouplings",
     "Indicate that no running couplings should be used.",
     &MatchboxMEBase::theFixedCouplings, false, false, false);
  static SwitchOption interfaceFixedCouplingsOn
    (interfaceFixedCouplings,
     "On",
     "On",
     true);
  static SwitchOption interfaceFixedCouplingsOff
    (interfaceFixedCouplings,
     "Off",
     "Off",
     false);

}

IBPtr MatchboxMEBase::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMEBase::fullclone() const {
  return new_ptr(*this);
}

void MatchboxMEBase::doinit() {
  MEBase::doinit();
  theAmplitude = dynamic_ptr_cast<Ptr<MatchboxAmplitude>::ptr>(amplitude());
  if ( matchboxAmplitude() ) {
    if ( matchboxAmplitude()->colourBasis() ) {
      size_t dim = 
	matchboxAmplitude()->colourBasis()->prepare(diagrams(),matchboxAmplitude()->noCorrelations());
      matchboxAmplitude()->colourBasisDim(dim);
    }
    matchboxAmplitude()->nLight(nLight());
  }
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxMEBase,MEBase>
describeHerwigMatchboxMEBase("Herwig::MatchboxMEBase", "HwMatchbox.so");

