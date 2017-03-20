// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the JetsPlusAnalysis class.
//

#include "JetsPlusAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/SubProcessGroup.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "Herwig/Sampling/GeneralSampler.h"

#include "Herwig/Utilities/XML/ElementIO.h"

using namespace Herwig;

JetsPlusAnalysis::JetsPlusAnalysis() 
  : theIsShowered(false) {}

JetsPlusAnalysis::~JetsPlusAnalysis() {}



#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

struct SortPt {

  inline bool operator()(const LorentzMomentum& a,
			 const LorentzMomentum& b) const {
    return a.perp() > b.perp();
  }

};

void JetsPlusAnalysis::reconstructJets(const ParticleVector& parts) {

  tcPDVector outType;
  vector<LorentzMomentum> outMomenta;

  for ( ParticleVector::const_iterator p = parts.begin();
	p != parts.end(); ++p ) {
    outType.push_back((**p).dataPtr());
    outMomenta.push_back((**p).momentum());
  }

  jetFinder()->cluster(outType, outMomenta,
		       tcCutsPtr(), tcPDPtr(),tcPDPtr());

  sort(outMomenta.begin(),outMomenta.end(),SortPt());

  for ( vector<Ptr<JetRegion>::ptr>::iterator j =
	  theJetRegions.begin(); j != theJetRegions.end(); ++j )
    (**j).reset();

  for ( size_t k = 0; k < outMomenta.size(); ++k ) {
    for ( vector<Ptr<JetRegion>::ptr>::const_iterator r = jetRegions().begin();
  	  r != jetRegions().end(); ++r ) {
      if ( (**r).matches(tCutsPtr(),k+1,outMomenta[k])) {
	jetMomentum(k+1) = outMomenta[k];
  	break;	
      }
    }
  }

}

void JetsPlusAnalysis::analyze(ParticleVector& parts, long id, double weight) {

  clear();
  reconstructHardObjects(parts);
  reconstructJets(parts);

  for ( map<string,LorentzMomentum>::const_iterator h = theHardObjects.begin();
	h != theHardObjects.end(); ++h ) {
    hardObjectProperties(h->first).count(h->second,weight,id);
    map<string,LorentzMomentum>::const_iterator g = h; ++g;
    for ( ; g != theHardObjects.end(); ++g ) {
      hardPairProperties(h->first,g->first).count(h->second,g->second,weight,id);
    }
  }

  unsigned int njets = 0;
  Energy jetSummedPerp = ZERO;
  double jetSummedRapidity = 0.0;
  double jetSummedPhi = 0.0;
  Energy jetSummedM = ZERO;

  nJetsInclusive().count(Statistics::EventContribution(njets,weight,0.0),id);
  if ( njets == theJets.size() )
    nJetsExclusive().count(Statistics::EventContribution(njets,weight,0.0),id);

  for ( map<unsigned int,LorentzMomentum>::const_iterator h = theJets.begin();
	h != theJets.end(); ++h ) {
    njets += 1;
    jetProperties(h->first).count(h->second,weight,id);
    jetInclusiveProperties().count(h->second,weight,id);
    nJetsInclusive().count(Statistics::EventContribution(njets,weight,0.0),id);
    if ( njets == theJets.size() ) {
      exclusiveJetProperties(h->first).count(h->second,weight,id);
      nJetsExclusive().count(Statistics::EventContribution(njets,weight,0.0),id);
    }
    jetSummedPerp += h->second.perp();
    jetSummedRapidity += h->second.rapidity();
    jetSummedPhi += h->second.phi();
    jetSummedM += h->second.m();
    map<unsigned int,LorentzMomentum>::const_iterator g = h; ++g;
    for ( ; g != theJets.end(); ++g ) {
      jetPairProperties(h->first,g->first).count(h->second,g->second,weight,id);
      map<unsigned int,LorentzMomentum>::const_iterator g1 = g; ++g1;
      for ( ; g1 != theJets.end(); ++g1 ) {
	LorentzMomentum p123 =
	  h->second + g->second + g1->second;
	threeJetProperties(h->first,g->first,g1->first).count(p123,weight,id);
	map<unsigned int,LorentzMomentum>::const_iterator g2 = g1; ++g2;
	for ( ; g2 != theJets.end(); ++g2 ) {
	  LorentzMomentum p1234 =
	    h->second + g->second + g1->second + g2->second;
	  fourJetProperties(h->first,g->first,g1->first,g2->first).count(p1234,weight,id);
	}
      }
    }
  }

  if ( njets > 0 )
    jetSummedProperties().count(jetSummedPerp,jetSummedRapidity,
				jetSummedPhi,jetSummedM,
				weight,id);

  if ( njets > 0 )
    jetAverageProperties().count(jetSummedPerp/njets,jetSummedRapidity/njets,
				 jetSummedPhi/njets,jetSummedM/njets,
				 weight,id);

  for ( map<string,LorentzMomentum>::const_iterator h = theHardObjects.begin();
	h != theHardObjects.end(); ++h ) {
    for ( map<unsigned int,LorentzMomentum>::const_iterator g = theJets.begin();
	  g != theJets.end(); ++g ) {
      jetHardPairProperties(g->first,h->first).count(g->second,h->second,weight,id);
    }
  }

  analyzeSpecial(id,weight);

}

void JetsPlusAnalysis::analyze(tEventPtr event, long ieve, int, int) {

  // doing nothing
  // AnalysisHandler::analyze(event, ieve, loop, state);

  Ptr<StandardEventHandler>::tptr seh =
    dynamic_ptr_cast<Ptr<StandardEventHandler>::tptr>(generator()->eventHandler());
  Ptr<GeneralSampler>::tptr sampler =
    dynamic_ptr_cast<Ptr<GeneralSampler>::tptr>(seh->sampler());

  double norm = sampler->maxXSec()/picobarn;

  if ( !theIsShowered ) {

    tSubProPtr sub = event->primarySubProcess();
    Ptr<SubProcessGroup>::tptr grp = 
      dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(sub);

    ParticleVector hfs = sub->outgoing();

    analyze(hfs,ieve,norm*event->weight()*sub->groupWeight());

    if ( grp ) {
    
      for ( SubProcessVector::const_iterator s = grp->dependent().begin();
	    s != grp->dependent().end(); ++s ) {
	ParticleVector fs = (**s).outgoing();
	analyze(fs,ieve,norm*event->weight()*(**s).groupWeight());
      }

    }

  } else {

    ParticleVector fs;
    event->getFinalState(fs);
    analyze(fs,ieve,norm*event->weight());

  }

}

void JetsPlusAnalysis::dofinish() {
  AnalysisHandler::dofinish();

  Ptr<StandardEventHandler>::tptr seh =
    dynamic_ptr_cast<Ptr<StandardEventHandler>::tptr>(generator()->eventHandler());
  Ptr<GeneralSampler>::tptr sampler =
    dynamic_ptr_cast<Ptr<GeneralSampler>::tptr>(seh->sampler());

  unsigned long attemptedPoints = sampler->attempts();
  double sumOfWeights = sampler->sumWeights();
  double sumOfSquaredWeights = sampler->sumWeights2();
  CrossSection maxXSection = sampler->maxXSec();

  XML::Element elem(XML::ElementTypes::Element,"Run");

  elem.appendAttribute("name",generator()->runName());
  elem.appendAttribute("attemptedPoints",attemptedPoints);
  elem.appendAttribute("sumOfWeights",sumOfWeights*maxXSection/picobarn);
  elem.appendAttribute("sumOfSquaredWeights",sumOfSquaredWeights*sqr(maxXSection/picobarn));

  XML::Element xhistos(XML::ElementTypes::Element,"Histograms");

  for ( map<string,ObjectProperties>::iterator h = theHardObjectProperties.begin();
	h != theHardObjectProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<unsigned int,ObjectProperties>::iterator h = theJetProperties.begin();
	h != theJetProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<unsigned int,ObjectProperties>::iterator h = theExclusiveJetProperties.begin();
	h != theExclusiveJetProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  if ( !theJetInclusiveProperties.pt.bins().empty() ) {
    theJetInclusiveProperties.finalize(xhistos);
  }

  if ( !theJetSummedProperties.pt.bins().empty() ) {
    theJetSummedProperties.finalize(xhistos);
  }

  if ( !theJetAverageProperties.pt.bins().empty() ) {
    theJetAverageProperties.finalize(xhistos);
  }

  if ( !theNJetsInclusive.bins().empty() ) {
    theNJetsInclusive.finalize();
    xhistos.append(theNJetsInclusive.toXML());
  }

  if ( !theNJetsExclusive.bins().empty() ) {
    theNJetsExclusive.finalize();
    xhistos.append(theNJetsExclusive.toXML());
  }

  for ( map<pair<string,string>,PairProperties>::iterator h = theHardPairProperties.begin();
	h != theHardPairProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = theJetPairProperties.begin();
	h != theJetPairProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<pair<unsigned int,string>,PairProperties>::iterator h = theJetHardPairProperties.begin();
	h != theJetHardPairProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int>,ObjectProperties>::iterator h =
	  theThreeJetProperties.begin(); h != theThreeJetProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties>::iterator h =
	  theFourJetProperties.begin(); h != theFourJetProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  finalize(xhistos);

  elem.append(xhistos);

  string fname = generator()->filename() + string("-") + name() + string(".xml");
  ofstream runXML(fname.c_str());
  runXML << setprecision(16);
  XML::ElementIO::put(elem,runXML);

}

IBPtr JetsPlusAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr JetsPlusAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void JetsPlusAnalysis::persistentOutput(PersistentOStream & os) const {
  os << theIsShowered << theJetFinder << theJetRegions;
}

void JetsPlusAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> theIsShowered >> theJetFinder >> theJetRegions;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<JetsPlusAnalysis,AnalysisHandler>
  describeHerwigJetsPlusAnalysis("Herwig::JetsPlusAnalysis", "JetCuts.so HwJetsAnalysis.so");

void JetsPlusAnalysis::Init() {

  static ClassDocumentation<JetsPlusAnalysis> documentation
    ("There is no documentation for the JetsPlusAnalysis class");

  static Reference<JetsPlusAnalysis,JetFinder> interfaceJetFinder
    ("JetFinder",
     "",
     &JetsPlusAnalysis::theJetFinder, false, false, true, false, false);

  static RefVector<JetsPlusAnalysis,JetRegion> interfaceJetRegions
    ("JetRegions",
     "",
     &JetsPlusAnalysis::theJetRegions, -1, false, false, true, false, false);

  static Switch<JetsPlusAnalysis,bool> interfaceIsShowered
    ("IsShowered",
     "",
     &JetsPlusAnalysis::theIsShowered, false, false, false);
  static SwitchOption interfaceIsShoweredYes
    (interfaceIsShowered,
     "Yes",
     "",
     true);
  static SwitchOption interfaceIsShoweredNo
    (interfaceIsShowered,
     "No",
     "",
     false);

}

