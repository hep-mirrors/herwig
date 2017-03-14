// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptonsJetsAnalysis class.
//

#include "LeptonsJetsAnalysis.h"
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

LeptonsJetsAnalysis::LeptonsJetsAnalysis() 
  : theIsShowered(false), theApplyCuts(false) {}

LeptonsJetsAnalysis::~LeptonsJetsAnalysis() {}



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

struct SortId {

  inline bool operator()(const pair<PID,LorentzMomentum>& a,
			 const pair<PID,LorentzMomentum>& b) const {
    // sort by abs(pid); if equal, first particle, then antiparticle
    // this puts pairs forming bosons next to each other
    long p1 = a.first;
    long p2 = b.first;
    if (abs(p1)==abs(p2)) {
      return p1 > p2;
    } else {
      return abs(p1) < abs(p2);
    }
  }

};

void LeptonsJetsAnalysis::reconstructJets(const ParticleVector& parts) {

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

void LeptonsJetsAnalysis::reconstructEWParticles(ParticleVector& parts) {

  vector< pair<PID,LorentzMomentum> > partall;
  vector<LorentzMomentum> partl, partnu, parth;
  LorentzMomentum ptmiss = LorentzMomentum(ZERO,ZERO,ZERO,ZERO);

  ParticleVector::iterator p = parts.begin();
  while (p != parts.end()) {
    PID pid = (**p).id();
    if ( ( static_cast<long>(pid) == ParticleID::eminus    ) || 
         ( static_cast<long>(pid) == ParticleID::eplus     ) || 
         ( static_cast<long>(pid) == ParticleID::muminus   ) || 
         ( static_cast<long>(pid) == ParticleID::muplus    ) || 
         ( static_cast<long>(pid) == ParticleID::tauminus  ) || 
         ( static_cast<long>(pid) == ParticleID::tauplus   ) ) {
      partall.push_back(pair<PID,LorentzMomentum>(pid,(**p).momentum()));
      partl.push_back((**p).momentum());
      p = parts.erase(p);
    } else
    if ( ( static_cast<long>(pid) == ParticleID::nu_e      ) || 
         ( static_cast<long>(pid) == ParticleID::nu_ebar   ) || 
         ( static_cast<long>(pid) == ParticleID::nu_mu     ) || 
         ( static_cast<long>(pid) == ParticleID::nu_mubar  ) || 
         ( static_cast<long>(pid) == ParticleID::nu_tau    ) || 
         ( static_cast<long>(pid) == ParticleID::nu_taubar ) ) {
      partall.push_back(pair<PID,LorentzMomentum>(pid,(**p).momentum()));
      partnu.push_back((**p).momentum());
      ptmiss += (**p).momentum();
      p = parts.erase(p);
    } else
    if (   static_cast<long>(pid) == ParticleID::h0          ) {
      partall.push_back(pair<PID,LorentzMomentum>(pid,(**p).momentum()));
      parth.push_back((**p).momentum());
      p = parts.erase(p);
    } else 
      p++;

  }

  sort(partall.begin(),partall.end(),SortId());
  sort(partl.begin(),partl.end(),SortPt());
  sort(partnu.begin(),partnu.end(),SortPt());
  sort(parth.begin(),parth.end(),SortPt());

  // make missing transverse momentum transverse and also add as last entry in EWID
  ptmiss.setE(ptmiss.perp());
  ptmiss.setZ(0*GeV);
  partall.push_back(pair<PID,LorentzMomentum>(ParticleID::nu_e,ptmiss));

  for ( size_t k = 0; k < partall.size(); ++k ) 
    eWIDMomentum(k+1) = partall[k].second;
  for ( size_t k = 0; k < partl.size(); ++k ) 
    chargedLeptonMomentum(k+1) = partl[k];
  for ( size_t k = 0; k < partnu.size(); ++k ) 
    neutrinoMomentum(k+1) = partnu[k];
  for ( size_t k = 0; k < parth.size(); ++k ) 
    higgsMomentum(k+1) = parth[k];
  pTmissMomentum() = ptmiss;

}

void LeptonsJetsAnalysis::analyze(ParticleVector& parts, long id, double weight) {

  clear();
  reconstructEWParticles(parts);
  reconstructJets(parts);

  if ( theApplyCuts ) {
// VBF cuts
    if ( nJets()<2 ) return;
    if ( (jetMomentum(1)+jetMomentum(2)).m() < 600*GeV ) return;
    if ( abs(jetMomentum(1).rapidity()-jetMomentum(2).rapidity()) < 3.6 ) return;
    // if ( jetMomentum(1).rapidity()*jetMomentum(2).rapidity() > 0 ) return;
    for ( map<unsigned int,LorentzMomentum>::const_iterator h = theChargedLeptons.begin();
          h != theChargedLeptons.end(); ++h ) {
      if ( h->second.perp() < 20*GeV ) return;
      if ( abs(h->second.rapidity()) > 2.5 ) return;
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
	threeJetProperties(h->first,g->first,g1->first).count(h->second,g->second,g1->second,weight,id);
	map<unsigned int,LorentzMomentum>::const_iterator g2 = g1; ++g2;
	for ( ; g2 != theJets.end(); ++g2 ) {
	  LorentzMomentum p1234 =
	    h->second + g->second + g1->second + g2->second;
	  fourJetProperties(h->first,g->first,g1->first,g2->first).count(p1234,weight,id);
	}
      }
      for ( map<unsigned int,LorentzMomentum>::const_iterator g1 = theEWIDs.begin();
          g1 != theEWIDs.end(); ++g1 ) 
        jetPairEWIDTripleProperties(h->first,g->first,g1->first).count(h->second,g->second,g1->second,weight,id);
      for ( map<unsigned int,LorentzMomentum>::const_iterator g1 = theChargedLeptons.begin();
          g1 != theChargedLeptons.end(); ++g1 ) 
        jetPairChargedLeptonTripleProperties(h->first,g->first,g1->first).count(h->second,g->second,g1->second,weight,id);
      for ( map<unsigned int,LorentzMomentum>::const_iterator g1 = theNeutrinos.begin();
          g1 != theNeutrinos.end(); ++g1 ) 
        jetPairNeutrinoTripleProperties(h->first,g->first,g1->first).count(h->second,g->second,g1->second,weight,id);
      jetPairPTmissTripleProperties(h->first,g->first).count(h->second,g->second,pTmissMomentum(),weight,id);
      for ( map<unsigned int,LorentzMomentum>::const_iterator g1 = theHiggs.begin();
          g1 != theHiggs.end(); ++g1 ) 
        jetPairHiggsTripleProperties(h->first,g->first,g1->first).count(h->second,g->second,g1->second,weight,id);
    }
    for ( map<unsigned int,LorentzMomentum>::const_iterator g = theEWIDs.begin();
	g != theEWIDs.end(); ++g ) 
      jetEWIDPairProperties(h->first,g->first).count(h->second,g->second,weight,id);
    for ( map<unsigned int,LorentzMomentum>::const_iterator g = theChargedLeptons.begin();
	g != theChargedLeptons.end(); ++g ) 
      jetChargedLeptonPairProperties(h->first,g->first).count(h->second,g->second,weight,id);
    for ( map<unsigned int,LorentzMomentum>::const_iterator g = theNeutrinos.begin();
	g != theNeutrinos.end(); ++g ) 
      jetNeutrinoPairProperties(h->first,g->first).count(h->second,g->second,weight,id);
    jetPTmissPairProperties(h->first).count(h->second,pTmissMomentum(),weight,id);
    for ( map<unsigned int,LorentzMomentum>::const_iterator g = theHiggs.begin();
	g != theHiggs.end(); ++g ) 
      jetHiggsPairProperties(h->first,g->first).count(h->second,g->second,weight,id);
  }

  if ( njets > 0 )
    jetSummedProperties().count(jetSummedPerp,jetSummedRapidity,
				jetSummedPhi,jetSummedM,
				weight,id);

  if ( njets > 0 )
    jetAverageProperties().count(jetSummedPerp/njets,jetSummedRapidity/njets,
				 jetSummedPhi/njets,jetSummedM/njets,
				 weight,id);

  for ( map<unsigned int,LorentzMomentum>::const_iterator h = theEWIDs.begin();
	h != theEWIDs.end(); ++h ) {
    eWIDProperties(h->first).count(h->second,weight,id);
    map<unsigned int,LorentzMomentum>::const_iterator g = h; ++g;
    for ( ; g != theEWIDs.end(); ++g ) {
      eWIDPairProperties(h->first,g->first).count(h->second,g->second,weight,id);
      map<unsigned int,LorentzMomentum>::const_iterator g1 = g; ++g1;
      for ( ; g1 != theEWIDs.end(); ++g1 ) {
	threeEWIDProperties(h->first,g->first,g1->first).count(h->second,g->second,g1->second,weight,id);
	map<unsigned int,LorentzMomentum>::const_iterator g2 = g1; ++g2;
	for ( ; g2 != theEWIDs.end(); ++g2 ) {
	  LorentzMomentum p1234 =
	    h->second + g->second + g1->second + g2->second;
	  fourEWIDProperties(h->first,g->first,g1->first,g2->first).count(p1234,weight,id);
	}
      }
    }
  }

  for ( map<unsigned int,LorentzMomentum>::const_iterator h = theChargedLeptons.begin();
	h != theChargedLeptons.end(); ++h ) {
    chargedLeptonProperties(h->first).count(h->second,weight,id);
    map<unsigned int,LorentzMomentum>::const_iterator g = h; ++g;
    for ( ; g != theChargedLeptons.end(); ++g ) {
      chargedLeptonPairProperties(h->first,g->first).count(h->second,g->second,weight,id);
      map<unsigned int,LorentzMomentum>::const_iterator g1 = g; ++g1;
      for ( ; g1 != theChargedLeptons.end(); ++g1 ) {
	threeChargedLeptonProperties(h->first,g->first,g1->first).count(h->second,g->second,g1->second,weight,id);
	map<unsigned int,LorentzMomentum>::const_iterator g2 = g1; ++g2;
	for ( ; g2 != theChargedLeptons.end(); ++g2 ) {
	  LorentzMomentum p1234 =
	    h->second + g->second + g1->second + g2->second;
	  fourChargedLeptonProperties(h->first,g->first,g1->first,g2->first).count(p1234,weight,id);
	}
      }
    }
  }

  for ( map<unsigned int,LorentzMomentum>::const_iterator h = theNeutrinos.begin();
	h != theNeutrinos.end(); ++h ) {
    neutrinoProperties(h->first).count(h->second,weight,id);
  }
  pTmissProperties().count(pTmissMomentum(),weight,id);

  for ( map<unsigned int,LorentzMomentum>::const_iterator h = theHiggs.begin();
	h != theHiggs.end(); ++h ) {
    higgsProperties(h->first).count(h->second,weight,id);
  }

  analyzeSpecial(id,weight);

}

void LeptonsJetsAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);

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

void LeptonsJetsAnalysis::dofinish() {
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

  for ( map<unsigned int,ObjectProperties>::iterator h = theEWIDProperties.begin();
	h != theEWIDProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<unsigned int,ObjectProperties>::iterator h = theChargedLeptonProperties.begin();
	h != theChargedLeptonProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<unsigned int,ObjectProperties>::iterator h = theNeutrinoProperties.begin();
	h != theNeutrinoProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  thePTmissProperties.finalize(xhistos);
 

  for ( map<unsigned int,ObjectProperties>::iterator h = theHiggsProperties.begin();
	h != theHiggsProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = theJetPairProperties.begin();
	h != theJetPairProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = theJetEWIDPairProperties.begin();
	h != theJetEWIDPairProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = theJetChargedLeptonPairProperties.begin();
	h != theJetChargedLeptonPairProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = theJetNeutrinoPairProperties.begin();
	h != theJetNeutrinoPairProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<unsigned int,PairProperties>::iterator h = theJetPTmissPairProperties.begin();
	h != theJetPTmissPairProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = theJetHiggsPairProperties.begin();
	h != theJetHiggsPairProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = theEWIDPairProperties.begin();
	h != theEWIDPairProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = theChargedLeptonPairProperties.begin();
	h != theChargedLeptonPairProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator h =
	  theThreeJetProperties.begin(); h != theThreeJetProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator h =
	  theJetPairEWIDTripleProperties.begin(); h != theJetPairEWIDTripleProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator h =
	  theJetPairChargedLeptonTripleProperties.begin(); h != theJetPairChargedLeptonTripleProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator h =
	  theJetPairNeutrinoTripleProperties.begin(); h != theJetPairNeutrinoTripleProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<pair<unsigned int,unsigned int>,TripleProperties>::iterator h =
	  theJetPairPTmissTripleProperties.begin(); h != theJetPairPTmissTripleProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator h =
	  theJetPairHiggsTripleProperties.begin(); h != theJetPairHiggsTripleProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator h =
	  theThreeEWIDProperties.begin(); h != theThreeEWIDProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator h =
	  theThreeChargedLeptonProperties.begin(); h != theThreeChargedLeptonProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties>::iterator h =
	  theFourJetProperties.begin(); h != theFourJetProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties>::iterator h =
	  theFourEWIDProperties.begin(); h != theFourEWIDProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  for ( map<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties>::iterator h =
	  theFourChargedLeptonProperties.begin(); h != theFourChargedLeptonProperties.end(); ++h ) {
    h->second.finalize(xhistos);
  }

  finalize(xhistos);

  elem.append(xhistos);

  string fname = generator()->filename() + string("-") + name() + string(".xml");
  ofstream runXML(fname.c_str());
  runXML << setprecision(16);
  XML::ElementIO::put(elem,runXML);

}

IBPtr LeptonsJetsAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr LeptonsJetsAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void LeptonsJetsAnalysis::persistentOutput(PersistentOStream & os) const {
  os << theIsShowered << theApplyCuts << theJetFinder << theJetRegions;
}

void LeptonsJetsAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> theIsShowered >> theApplyCuts >> theJetFinder >> theJetRegions;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<LeptonsJetsAnalysis,AnalysisHandler>
  describeHerwigLeptonsJetsAnalysis("Herwig::LeptonsJetsAnalysis", "JetCuts.so HwJetsAnalysis.so");

void LeptonsJetsAnalysis::Init() {

  static ClassDocumentation<LeptonsJetsAnalysis> documentation
    ("General-purpose analysis for processes with jets and leptons");

  static Reference<LeptonsJetsAnalysis,JetFinder> interfaceJetFinder
    ("JetFinder",
     "",
     &LeptonsJetsAnalysis::theJetFinder, false, false, true, false, false);

  static RefVector<LeptonsJetsAnalysis,JetRegion> interfaceJetRegions
    ("JetRegions",
     "",
     &LeptonsJetsAnalysis::theJetRegions, -1, false, false, true, false, false);

  static Switch<LeptonsJetsAnalysis,bool> interfaceIsShowered
    ("IsShowered",
     "",
     &LeptonsJetsAnalysis::theIsShowered, false, false, false);
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

  static Switch<LeptonsJetsAnalysis,bool> interfaceApplyCuts
    ("ApplyCuts",
     "",
     &LeptonsJetsAnalysis::theApplyCuts, false, false, false);
  static SwitchOption interfaceApplyCutsYes
    (interfaceApplyCuts,
     "Yes",
     "",
     true);
  static SwitchOption interfaceApplyCutsNo
    (interfaceApplyCuts,
     "No",
     "",
     false);

}

