// -*- C++ -*-
//
// MatchboxAmplitude.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitude class.
//

#include "MatchboxAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SpinorHelicity.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"
#include "MatchboxMEBase.h"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <iterator>
using std::ostream_iterator;

using namespace Herwig;

MatchboxAmplitude::MatchboxAmplitude() 
  : Amplitude() {}

MatchboxAmplitude::~MatchboxAmplitude() {}

void MatchboxAmplitude::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb << theColourBasis;
}

void MatchboxAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb >> theColourBasis;
  lastMatchboxXComb(theLastXComb);
}

void MatchboxAmplitude::doinit() {
  Amplitude::doinit();
  if ( colourBasis() )
    colourBasis()->init();
}

void MatchboxAmplitude::doinitrun() {
  Amplitude::doinitrun();
  if ( colourBasis() )
    colourBasis()->initrun();
}

void MatchboxAmplitude::cloneDependencies(const std::string&) {}

Ptr<MatchboxMEBase>::ptr MatchboxAmplitude::makeME(const PDVector&) const {
  return new_ptr(MatchboxMEBase());
}

Selector<const ColourLines *> MatchboxAmplitude::colourGeometries(tcDiagPtr d) const {
  if ( haveColourFlows() )
    return colourBasis()->colourGeometries(d,lastLargeNAmplitudes());
  return Selector<const ColourLines *>();
}

void MatchboxAmplitude::olpOrderFileHeader(ostream& os) const {

  os << "# OLP order file created by Herwig++/Matchbox\n\n";

  os << "InterfaceVersion          BLHA2\n\n";

  os << "Model                     SM\n"
     << "CorrectionType            QCD\n"
     << "IRregularisation          " << (isDR() ? "DRED" : "CDR") << "\n"
     << "Extra HelAvgInitial       no\n"
     << "Extra ColAvgInitial       no\n"
     << "Extra MCSymmetrizeFinal   no\n";

  os << "\n";

}

void MatchboxAmplitude::olpOrderFileProcesses(ostream& os,
					      const map<pair<Process,int>,int>& proc) const {

  map<int,pair<Process,int> > sorted;

  for ( map<pair<Process,int>,int>::const_iterator p = proc.begin();
	p != proc.end(); ++p ) {
    sorted[p->second] = p->first;
  }

  unsigned int currentOrderInAlphaS = sorted.begin()->second.first.orderInAlphaS;
  unsigned int currentOrderInAlphaEW = sorted.begin()->second.first.orderInAlphaEW;
  int currentType = sorted.begin()->second.second;

  os << "AlphasPower               " << currentOrderInAlphaS << "\n"
     << "AlphaPower                " << currentOrderInAlphaEW << "\n"
     << "AmplitudeType             ";
  if ( currentType == ProcessType::treeME2 ) {
    os << "tree\n";
  } else if ( currentType == ProcessType::oneLoopInterference ) {
    os << "loop\n";
  } else if ( currentType == ProcessType::colourCorrelatedME2 ) {
    os << "cctree\n";
  } else if ( currentType == ProcessType::spinColourCorrelatedME2 ) {
    os << "sctree\n";
  } else assert(false);


  for ( map<int,pair<Process,int> >::const_iterator p = sorted.begin();
	p != sorted.end(); ++p ) {

    if ( currentOrderInAlphaS != p->second.first.orderInAlphaS ) {
      currentOrderInAlphaS = p->second.first.orderInAlphaS;
      os << "AlphasPower               " << currentOrderInAlphaS << "\n";
    }

    if ( currentOrderInAlphaEW != p->second.first.orderInAlphaEW ) {
      currentOrderInAlphaEW = p->second.first.orderInAlphaEW;
      os << "AlphaPower                " << currentOrderInAlphaEW << "\n";
    }

    if ( currentType != p->second.second ) {
      currentType = p->second.second;
      os << "AmplitudeType             ";
      if ( currentType == ProcessType::treeME2 ) {
	os << "tree\n";
      } else if ( currentType == ProcessType::oneLoopInterference ) {
	os << "loop\n";
      } else if ( currentType == ProcessType::colourCorrelatedME2 ) {
	os << "cctree\n";
      } else if ( currentType == ProcessType::spinColourCorrelatedME2 ) {
	os << "sctree\n";
      } else assert(false);
    }

    os << p->second.first.legs[0]->id() << " "
       << p->second.first.legs[1]->id() << " -> ";
    for ( PDVector::const_iterator o = p->second.first.legs.begin() + 2;
	  o != p->second.first.legs.end(); ++o ) {
      os << (**o).id() << " ";
    }
    os << "\n";

  }

}

bool MatchboxAmplitude::startOLP(const map<pair<Process,int>,int>& procs) {

  string orderFileName = name() + ".OLPOrder.lh";
  ofstream orderFile(orderFileName.c_str());

  olpOrderFileHeader(orderFile);
  olpOrderFileProcesses(orderFile,procs);

  string contractFileName = name() + ".OLPContract.lh";

  signOLP(orderFileName, contractFileName);

  // TODO check the contract file

  int status = 0;
  startOLP(contractFileName, status);

  if ( status != 1 )
    return false;

  return true;

}


struct orderPartonData {

  bool operator()(const pair<tcPDPtr,int>& a, 
		  const pair<tcPDPtr,int>& b) const {

    if ( a.first == b.first )
      return a.second < b.second;

    int acolour = a.first->iColour();
    int bcolour = b.first->iColour();

    if ( abs(acolour) != abs(bcolour) )
      return abs(acolour) < abs(bcolour);

    if ( a.first->iSpin() != b.first->iSpin() )
      return a.first->iSpin() < b.first->iSpin();

    int acharge = a.first->iCharge();
    int bcharge = b.first->iCharge();

    if ( abs(acharge) != abs(bcharge) )
      return abs(acharge) < abs(bcharge);

    if ( abs(a.first->id()) != abs(b.first->id()) )
      return abs(a.first->id()) < abs(b.first->id());

    return a.first->id() > b.first->id();

  }

};

void MatchboxAmplitude::setXComb(tStdXCombPtr xc) {
  theLastXComb = xc;
  lastMatchboxXComb(xc);
  fillCrossingMap();
  if ( treeAmplitudes() || oneLoopAmplitudes() )
    for ( size_t k = 0 ; k < meMomenta().size(); ++k )
      amplitudeMomenta()[k] = amplitudeMomentum(k);
}

void MatchboxAmplitude::fillCrossingMap(size_t shift) {

  if ( !amplitudePartonData().empty() )
    return;

  double csign = 1.;
  set<pair<tcPDPtr,int>,orderPartonData > processLegs;
  for ( unsigned int l = 0; l < mePartonData().size(); ++l ) {
    if ( l > 1 )
      processLegs.insert(make_pair(mePartonData()[l],l));
    else {
      if ( mePartonData()[l]->CC() ) {
	processLegs.insert(make_pair(mePartonData()[l]->CC(),l));
	if ( mePartonData()[l]->iSpin() == PDT::Spin1Half )
	  csign *= -1.;
      } else {
	processLegs.insert(make_pair(mePartonData()[l],l));
      }
    }
  }

  crossingSign(csign);

  set<pair<tcPDPtr,int> > amplitudeLegs;
  crossingMap().resize(mePartonData().size());
  amplitudePartonData().resize(mePartonData().size());
  amplitudeMomenta().resize(mePartonData().size());

  int ampCount = 0;

  // process legs are already sorted, we only need to arrange for
  // adjacent particles and anti-particles
  while ( !processLegs.empty() ) {
    set<pair<tcPDPtr,int>,orderPartonData >::iterator next
      = processLegs.begin();
    while ( next->first->id() < 0 ) {
      if ( ++next == processLegs.end() )
	break;
    }
    assert(next != processLegs.end());
    crossingMap()[ampCount] = next->second - shift;
    amplitudeLegs.insert(make_pair(next->first,ampCount));
    tcPDPtr check = next->first;
    processLegs.erase(next);
    ++ampCount;
    if ( check->CC() ) {
      set<pair<tcPDPtr,int>,orderPartonData>::iterator checkcc
	= processLegs.end();
      for ( set<pair<tcPDPtr,int>,orderPartonData>::iterator c = processLegs.begin();
	    c != processLegs.end(); ++c ) {
	if ( c->first == check->CC() ) {
	  checkcc = c; break;
	}
      }
      if ( checkcc == processLegs.end() )
	for ( set<pair<tcPDPtr,int>,orderPartonData>::iterator c = processLegs.begin();
	      c != processLegs.end(); ++c ) {
	  if ( !SU2Helper::SU2CC(check) )
	    continue;
	  if ( c->first == SU2Helper::SU2CC(check)->CC() ) {
	    checkcc = c; break;
	  }
	}
      if ( checkcc == processLegs.end() ) {
	int f = SU2Helper::family(check);
	for ( int i = 1 - f; i < 5 - f; i++ ) {
	  bool gotone = false;
	  for ( set<pair<tcPDPtr,int>,orderPartonData>::iterator c = processLegs.begin();
		c != processLegs.end(); ++c ) {
	    if ( !SU2Helper::SU2CC(check,i) )
	      continue;
	    if ( c->first == SU2Helper::SU2CC(check,i)->CC() ) {
	      checkcc = c; gotone = true; break;
	    }
	  }
	  if ( gotone )
	    break;
	}
      }
      // default to just pick the next available anti-particle
      if ( checkcc == processLegs.end() ) {
	checkcc = processLegs.begin();
	while ( checkcc->first->id() > 0 )
	  if ( ++checkcc == processLegs.end() )
	    break;
      }
      assert(checkcc != processLegs.end());
      crossingMap()[ampCount] = checkcc->second - shift;
      amplitudeLegs.insert(make_pair(checkcc->first,ampCount));
      processLegs.erase(checkcc);
      ++ampCount;
    }
  }

  for ( set<pair<tcPDPtr,int> >::const_iterator l = amplitudeLegs.begin();
	l != amplitudeLegs.end(); ++l )
    amplitudePartonData()[l->second] = l->first;

  if ( colourBasis() ) {
    assert(colourBasis()->indexMap().find(mePartonData()) !=
	   colourBasis()->indexMap().end());
    const map<size_t,size_t> colourCross = 
      colourBasis()->indexMap().find(mePartonData())->second;
    for ( size_t k = 0; k < crossingMap().size(); ++k ) {
      if ( colourCross.find(crossingMap()[k]) !=
	   colourCross.end() ) {
	size_t ccross = colourCross.find(crossingMap()[k])->second;
	amplitudeToColourMap()[k] = ccross;
	colourToAmplitudeMap()[ccross] = k;
      }
    }
  }

}

const string& MatchboxAmplitude::colourOrderingString(size_t id) const {

  static string empty = "";
  if ( !colourBasis() ) {
    return empty;
  }

  return colourBasis()->orderingString(mePartonData(),colourToAmplitudeMap(),id);

}

const set<vector<size_t> >& MatchboxAmplitude::colourOrdering(size_t id) const {

  static set<vector<size_t> > empty;
  if ( !colourBasis() ) {
    return empty;
  }

  return colourBasis()->ordering(mePartonData(),colourToAmplitudeMap(),id);

}

Lorentz5Momentum MatchboxAmplitude::amplitudeMomentum(int i) const {
  int iCrossed = crossingMap()[i];
  Lorentz5Momentum res = meMomenta()[iCrossed];
  if ( iCrossed < 2 )
    res = -res;
  res.setMass(meMomenta()[iCrossed].mass());
  Energy2 rho = res.t()*res.t() - res.mass2();
  res.setRho(sqrt(abs(rho)));
  return res;
}

set<vector<int> > MatchboxAmplitude::generateHelicities() const {
  set<vector<int> > res;
  vector<int> current(amplitudePartonData().size());
  doGenerateHelicities(res,current,0);
  return res;
}

void MatchboxAmplitude::doGenerateHelicities(set<vector<int> >& res,
					     vector<int>& current,
					     size_t pos) const {

  if ( pos == amplitudePartonData().size() ) {
    res.insert(current);
    return;
  }

  if ( amplitudePartonData()[pos]->iSpin() == PDT::Spin0 ||
       ( amplitudePartonData()[pos]->iSpin() == PDT::Spin1 &&
	 amplitudePartonData()[pos]->mass() != ZERO ) ) {
    current[pos] = 0;
    doGenerateHelicities(res,current,pos+1);
  } else if ( amplitudePartonData()[pos]->iSpin() == PDT::Spin1Half ||
	      amplitudePartonData()[pos]->iSpin() == PDT::Spin1 ) {
    current[pos] = 1;
    doGenerateHelicities(res,current,pos+1);
    current[pos] = -1;
    doGenerateHelicities(res,current,pos+1);
  }

}

void MatchboxAmplitude::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr) {

  if ( !calculateTreeAmplitudes() )
    return;

  bool initialized = !lastAmplitudes().empty();

  if ( !initialized ) {
    set<vector<int> > helicities = generateHelicities();
    for ( set<vector<int> >::const_iterator h = helicities.begin();
          h != helicities.end(); ++h ) {
      lastAmplitudes().insert(make_pair(*h,CVector(colourBasisDim())));
      lastLargeNAmplitudes().insert(make_pair(*h,CVector(colourBasisDim())));
    }
  }

  AmplitudeIterator amp = lastAmplitudes().begin();
  AmplitudeIterator lamp = lastLargeNAmplitudes().begin();
  for ( ;amp != lastAmplitudes().end(); ++amp, ++lamp ) {
    for ( size_t k = 0; k < colourBasisDim(); ++k )
      amp->second(k) = evaluate(k,amp->first,lamp->second(k));
  }

  if ( !initialized ) {
    map<vector<int>,CVector> clean;
    for ( map<vector<int>,CVector>::const_iterator amp = lastAmplitudes().begin();
          amp != lastAmplitudes().end(); ++amp ) {
      bool nonZero = false;
      for ( size_t k = 0; k < colourBasisDim(); ++k ) {
        if ( amp->second(k) != Complex(0.0) ) {
          nonZero = true;
          break;
        }
      }
      if ( nonZero )
        clean.insert(*amp);
    }
    lastAmplitudes() = clean;
    clean.clear();
    for ( map<vector<int>,CVector>::const_iterator amp = lastLargeNAmplitudes().begin();
          amp != lastLargeNAmplitudes().end(); ++amp ) {
      bool nonZero = false;
      for ( size_t k = 0; k < colourBasisDim(); ++k ) {
        if ( amp->second(k) != Complex(0.0) ) {
          nonZero = true;
          break;
        }
      }
      if ( nonZero )
        clean.insert(*amp);
    }
    lastLargeNAmplitudes() = clean;
  }

  haveTreeAmplitudes();

}

void MatchboxAmplitude::prepareOneLoopAmplitudes(Ptr<MatchboxMEBase>::tcptr) {

  if ( !calculateOneLoopAmplitudes() )
    return;

  bool initialized = !lastOneLoopAmplitudes().empty();

  if ( !initialized ) {
    set<vector<int> > helicities = generateHelicities();
    for ( set<vector<int> >::const_iterator h = helicities.begin();
          h != helicities.end(); ++h ) {
      lastOneLoopAmplitudes().insert(make_pair(*h,CVector(colourBasisDim())));
    }
  }

  for ( AmplitudeIterator amp = lastOneLoopAmplitudes().begin();
	amp != lastOneLoopAmplitudes().end(); ++amp ) {
    for ( size_t k = 0; k < colourBasisDim(); ++k )
      amp->second(k) = evaluateOneLoop(k,amp->first);
  }

  if ( !initialized ) {
    map<vector<int>,CVector> clean;
    for ( map<vector<int>,CVector>::const_iterator amp = lastOneLoopAmplitudes().begin();
          amp != lastOneLoopAmplitudes().end(); ++amp ) {
      bool nonZero = false;
      for ( size_t k = 0; k < colourBasisDim(); ++k ) {
        if ( amp->second(k) != Complex(0.0) ) {
          nonZero = true;
          break;
        }
      }
      if ( nonZero )
        clean.insert(*amp);
    }
    lastOneLoopAmplitudes() = clean;
  }

  haveOneLoopAmplitudes();

}

Complex MatchboxAmplitude::value(const tcPDVector&,
				 const vector<Lorentz5Momentum>&, 
				 const vector<int>&) {
  assert(false && "ThePEG::Amplitude interface is not sufficient at the moment.");
  throw Exception() << "ThePEG::Amplitude interface is not sufficient at the moment."
		    << Exception::abortnow;
  return 0.;
}

double MatchboxAmplitude::me2() const {
  if ( !calculateTreeME2() )
    return lastTreeME2();
  lastTreeME2(crossingSign()*colourBasis()->me2(mePartonData(),lastAmplitudes()));
  return lastTreeME2();
}

double MatchboxAmplitude::oneLoopInterference() const {
  if ( !calculateOneLoopInterference() )
    return lastOneLoopInterference();
  lastOneLoopInterference(crossingSign()*
			  colourBasis()->interference(mePartonData(),
						      lastOneLoopAmplitudes(),lastAmplitudes()));
  return lastOneLoopInterference();
}

double MatchboxAmplitude::colourCorrelatedME2(pair<int,int> ij) const {
  double cfac = 1.;
  double Nc = generator()->standardModel()->Nc();
  if ( mePartonData()[ij.first]->iColour() == PDT::Colour8 ) {
    cfac = Nc;
  } else if ( mePartonData()[ij.first]->iColour() == PDT::Colour3 ||
	      mePartonData()[ij.first]->iColour() == PDT::Colour3bar ) {
    cfac = (sqr(Nc)-1.)/(2.*Nc);
  } else assert(false);
  if ( !calculateColourCorrelator(ij) )
    return lastColourCorrelator(ij)/cfac;
  double res =
    crossingSign()*colourBasis()->colourCorrelatedME2(ij,mePartonData(),lastAmplitudes());
  lastColourCorrelator(ij,res);
  return res/cfac;
}

double MatchboxAmplitude::largeNColourCorrelatedME2(pair<int,int> ij,
						    Ptr<ColourBasis>::tptr largeNBasis) const {
  double cfac = 1.;
  double Nc = generator()->standardModel()->Nc();
  if ( mePartonData()[ij.first]->iColour() == PDT::Colour8 ) {
    cfac = Nc;
  } else if ( mePartonData()[ij.first]->iColour() == PDT::Colour3 ||
	      mePartonData()[ij.first]->iColour() == PDT::Colour3bar ) {
    cfac = Nc/2.;
  } else assert(false);
  if ( !calculateLargeNColourCorrelator(ij) )
    return lastLargeNColourCorrelator(ij)/cfac;
  double res =
    crossingSign()*largeNBasis->colourCorrelatedME2(ij,mePartonData(),lastLargeNAmplitudes());
  lastLargeNColourCorrelator(ij,res);
  return res/cfac;
}

// compare int vectors modulo certain element
// which needs to differe between the two
bool equalsModulo(unsigned int i, const vector<int>& a, const vector<int>& b) {
  assert(a.size()==b.size());
  if ( a[i] == b[i] )
    return false;
  for ( unsigned int k = 0; k < a.size(); ++k ) {
    if ( k == i )
      continue;
    if ( a[k] != b[k] )
      return false;
  }
  return true;
}

LorentzVector<Complex> MatchboxAmplitude::plusPolarization(const Lorentz5Momentum& p,
							   const Lorentz5Momentum& n,
							   int) const {

  using namespace SpinorHelicity;

  LorentzVector<complex<Energy> > num =
    PlusSpinorCurrent(PlusConjugateSpinor(n),MinusSpinor(p)).eval();

  complex<Energy> den =
    sqrt(2.)*PlusSpinorProduct(PlusConjugateSpinor(n),PlusSpinor(p)).eval();

  LorentzVector<Complex> polarization(num.x()/den,num.y()/den,num.z()/den,num.t()/den);

  return polarization;

}

double MatchboxAmplitude::spinColourCorrelatedME2(pair<int,int> ij,
						  const SpinCorrelationTensor& c) const {

  Lorentz5Momentum p = meMomenta()[ij.first];
  Lorentz5Momentum n = meMomenta()[ij.second];

  LorentzVector<Complex> polarization = plusPolarization(p,n,ij.first);

  Complex pFactor = (polarization*c.momentum())/sqrt(abs(c.scale()));

  double avg =
    colourCorrelatedME2(ij)*(-c.diagonal()+ (c.scale() > ZERO ? 1. : -1.)*norm(pFactor));

  int iCrossed = -1;
  for ( unsigned int k = 0; k < crossingMap().size(); ++k )
    if ( crossingMap()[k] == ij.first ) {
      iCrossed = k;
      break;
    }      
  assert(iCrossed >= 0);

  Complex csCorr = 0.0;

  if ( calculateColourSpinCorrelator(ij) ) {
    set<const CVector*> done;
    for ( AmplitudeConstIterator a = lastAmplitudes().begin();
	  a != lastAmplitudes().end(); ++a ) {
      if ( done.find(&(a->second)) != done.end() )
	continue;
      AmplitudeConstIterator b = lastAmplitudes().begin();
      while ( !equalsModulo(iCrossed,a->first,b->first) )
	if ( ++b == lastAmplitudes().end() )
	  break;
      if ( b == lastAmplitudes().end() || done.find(&(b->second)) != done.end() )
	continue;
      done.insert(&(a->second)); done.insert(&(b->second));
      if ( a->first[iCrossed] == 1 )
	swap(a,b);
      csCorr += colourBasis()->colourCorrelatedInterference(ij,mePartonData(),a->second,b->second);
    }
    lastColourSpinCorrelator(ij,csCorr);
  } else {
    csCorr = lastColourSpinCorrelator(ij);
  }

  double corr = 
    2.*real(csCorr*sqr(pFactor));

  double Nc = generator()->standardModel()->Nc();
  double cfac = 1.;
  if ( mePartonData()[ij.first]->iColour() == PDT::Colour8 ) {
    cfac = Nc;
  } else if ( mePartonData()[ij.first]->iColour() == PDT::Colour3 ||
	      mePartonData()[ij.first]->iColour() == PDT::Colour3bar ) {
    cfac = (sqr(Nc)-1.)/(2.*Nc);
  } else assert(false);

  return 
    avg + crossingSign()*(c.scale() > ZERO ? 1. : -1.)*corr/cfac;

}

void MatchboxAmplitude::Init() {

  static ClassDocumentation<MatchboxAmplitude> documentation
    ("MatchboxAmplitude is the base class for amplitude "
     "implementations inside Matchbox.");

  static Reference<MatchboxAmplitude,ColourBasis> interfaceColourBasis
    ("ColourBasis",
     "Set the colour basis implementation.",
     &MatchboxAmplitude::theColourBasis, false, false, true, true, false);

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<MatchboxAmplitude,Amplitude>
describeMatchboxAmplitude("Herwig::MatchboxAmplitude", "Herwig.so");
