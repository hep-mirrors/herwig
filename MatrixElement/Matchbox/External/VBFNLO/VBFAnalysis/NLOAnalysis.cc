// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOAnalysis class.
//

#include "NLOAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "fastjet/ClusterSequence.hh"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/SubProcessGroup.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include<sys/stat.h>
#include <string.h>
#include <unistd.h>

using namespace Herwig;

inline double deltaphi(double phi1, double phi2){
  double diff=phi1-phi2;
  if(diff<-Constants::pi){
    diff+=(2.0*Constants::pi);
  }
  else if (diff>Constants::pi){
    diff-=(2.0*Constants::pi);
    }
  return diff;
}

NLOAnalysis::NLOAnalysis()
  : AnalysisHandler(), events(), higgs(),
    theEtaDetector(5.0), theRparam(0.7),
    theJetDefPt(20.0), plainNLO(false), jetCache(), jets(),
    theJetCacheSet(false), debuginfo(false), theRjj_min(0.7),
    theRapidityGap(0.0),themjj_min(0.0),
    theopposite_dir(false),thenjets_min(2),thenjets_ex(0),
    theDelYhj_min(0) {}

NLOAnalysis::~NLOAnalysis() {}

bool NLOAnalysis::allowedParticle (const Particle &p) {
  if (p.id()!=12 && p.id()!=14 && p.id()!=16 && p.id()!=18 &&p.id()!=25 && p.id()!=36 && p.id()!=82) return true;
  return false;
}

bool NLOAnalysis::passCuts(vector<fastjet::PseudoJet> jets){

  if ( thenjets_ex != 0 ){
    if ( jets.size() != thenjets_ex ) return false;
  }
  else {
    if ( jets.size() < thenjets_min ) return false; 
    if ( jets.size() == 1 && thenjets_min == 1) return true;
  }

  assert(jets.size() > 1);

  double Rjj=sqrt(pow(deltaphi(jets[0].phi_std(),jets[1].phi_std()),2)+pow(jets[0].pseudorapidity()-jets[1].pseudorapidity(),2));
  if ( Rjj < theRjj_min ) return false;

  if ( abs(jets[0].pseudorapidity()-jets[1].pseudorapidity()) < theRapidityGap ) return false;
  if ( (jets[0]+jets[1]).m() < themjj_min ) return false;

  if (theopposite_dir)
    if ( ( jets[0].rapidity()*jets[1].rapidity() > 0 ) ? true : false ) return false;
  return true;
}


void NLOAnalysis::analyze(tEventPtr event, long, int, int) {
  //AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  
  tSubProPtr subpro = event->primarySubProcess();
  Ptr<SubProcessGroup>::tptr subgrp = dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(subpro);

  if ( plainNLO ){
    if ( !subgrp ) { // born + virtual 
      //      cerr << "born +virtual\n" << flush;
      const ParticleVector& toAnalyze = subpro->outgoing();
      vector<fastjet::PseudoJet> input_particles = recombinables(toAnalyze);
      const double w = event->weight();
      doAnalyze(event, input_particles, input_particles, w);
    } else { 
      //      cerr << "real emission\n" << flush;
      // real emission
      const ParticleVector& toAnalyze = subgrp->outgoing();
      vector<fastjet::PseudoJet> input_particles = recombinables(toAnalyze);
      const double w = event->weight()*subgrp->groupWeight();
      doAnalyze(event, input_particles, input_particles, w);
      // dipoles
      for ( SubProcessVector::const_iterator s = subgrp->dependent().begin();
	    s != subgrp->dependent().end(); ++s ) {
	//	cerr << "dipole\n" << flush;
	const ParticleVector& toAnalyze = (**s).outgoing();
	vector<fastjet::PseudoJet> input_particles = recombinables(toAnalyze);
	const double w = event->weight()*(**s).groupWeight();
	doAnalyze(event, input_particles, input_particles, w);
      }
    }
  }
  else {
    //    cerr << "Shower\n" << flush;
    const tPVector& toAnalyze = (*event).getFinalState();
    vector<fastjet::PseudoJet> input_particles = recombinables(toAnalyze);
    // ParticleVector toAnalyze;
    // if ( jetCache.size() != 0 )
    //   event->getFinalState(toAnalyze);
    const ParticleVector& partonsToAnalyze = event->primarySubProcess()->outgoing();
    vector<fastjet::PseudoJet> hard_partons = recombinables(partonsToAnalyze); // HIGGS NOT INCLUDED THIS WAY
    const double w = event->weight();
    doAnalyze(event, input_particles, hard_partons, w);
  }
}


void NLOAnalysis::doAnalyze(tEventPtr event, const vector<fastjet::PseudoJet>& input_particles, const vector<fastjet::PseudoJet>& hard_partons, const double w){

  vector<fastjet::PseudoJet> partons = sorted_by_pt(getInRange(recombine(hard_partons))); //no higgs in here
  
  bool cutsOk = false;
  if ( !cacheIsSet() || plainNLO) {
    jets = sorted_by_pt(getInRange(recombine(input_particles))); 
    cutsOk = passCuts( jets );
  }
  else {
    jets = jetCache;
    clearJetCache();
    cutsOk = true;
  }

  const ParticleVector& MEParticles = event->primarySubProcess()->outgoing();
  fastjet::PseudoJet theHiggs = getHiggs(MEParticles);

  //now do the analysis with the jets
  if (jets.size() > 1 && cutsOk) {
    events.fill(jets,partons,theHiggs,w); 
  }

}

vector<fastjet::PseudoJet> NLOAnalysis::recombinables(const ParticleVector& p, bool sort_out){
  vector<fastjet::PseudoJet> recombinables;
  for (ParticleVector::const_iterator iter=p.begin(); iter!=p.end(); iter++){
    if (sort_out && !allowedParticle(**iter))
      continue;
    recombinables.push_back(fastjet::PseudoJet( (**iter).momentum().x()/GeV,(**iter).momentum().y()/GeV,(**iter).momentum().z()/GeV, (**iter).momentum().t()/GeV ));
  }
  return recombinables;
}

vector<fastjet::PseudoJet> NLOAnalysis::recombinables(const tPVector& p, bool sort_out){
  vector<fastjet::PseudoJet> recombinables;
  for (tPVector::const_iterator iter=p.begin(); iter!=p.end(); iter++){
    if (sort_out && !allowedParticle(**iter))
      continue;
    recombinables.push_back(fastjet::PseudoJet( (**iter).momentum().x()/GeV,(**iter).momentum().y()/GeV,(**iter).momentum().z()/GeV, (**iter).momentum().t()/GeV ));
  }
  return recombinables;
}

fastjet::PseudoJet NLOAnalysis::getHiggs(const ParticleVector& p){
  fastjet::PseudoJet theHiggs;
  for (ParticleVector::const_iterator iter=p.begin(); iter!=p.end(); iter++){
    if ((**iter).id()!=25)
      continue;
    theHiggs = (fastjet::PseudoJet( (**iter).momentum().x()/GeV,(**iter).momentum().y()/GeV,(**iter).momentum().z()/GeV, (**iter).momentum().t()/GeV ));
    break;
  }
  return theHiggs;
}

vector<fastjet::PseudoJet> NLOAnalysis::recombine(const vector<fastjet::PseudoJet>& p){
  //use fastjet to recombine jets
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
  fastjet::JetDefinition jet_def(fastjet::kt_algorithm, theRparam, recomb_scheme, strategy);
  fastjet::ClusterSequence clust_seq(p, jet_def);
  
  vector<fastjet::PseudoJet> recoJets = clust_seq.inclusive_jets(theJetDefPt);
  recoJets = sorted_by_pt(recoJets);

  return recoJets;
}

vector<fastjet::PseudoJet> NLOAnalysis::getInRange(const vector<fastjet::PseudoJet>& j){
  vector<fastjet::PseudoJet> inRange;
  for (vector<fastjet::PseudoJet>::const_iterator iter = j.begin();
       iter != j.end(); iter++){
    //if ( abs(iter->pseudorapidity()) < theEtaDetector ) 
    inRange.push_back(*iter);
  }
  return inRange;
}

LorentzRotation NLOAnalysis::transform(tcEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void NLOAnalysis::analyze(const tPVector & particles, double) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void NLOAnalysis::analyze(tPPtr, double) {}

IBPtr NLOAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr NLOAnalysis::fullclone() const {
  return new_ptr(*this);
}

void NLOAnalysis::doinitrun(){
  mkdir("runresults",0777);
}

void NLOAnalysis::dofinish(){
  
  if( chdir("runresults") != 0 ) cerr << "NLOAnalysis could not move into subdirectory.\n" << flush;
  const char * cAnaName = name().c_str();
  char AnaName [100];
  strcpy(AnaName,cAnaName);
  events.writetofolder(name().c_str());
  //higgs.writetofolder(strcat(AnaName,"_higgs"));
  // events.writetofolder("events");
  // higgs.writetofolder("higgs");
  if( chdir("..") != 0 ) cerr << "NLOAnalysis could not move into top directory.\n" << flush;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOAnalysis::persistentOutput(PersistentOStream & os) const {
  os << theEtaDetector << theRparam << theJetDefPt << plainNLO 
     << theRjj_min << theRapidityGap << theopposite_dir << themjj_min
     << theDelYhj_min << thenjets_min << thenjets_ex; 
}

void NLOAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> theEtaDetector >> theRparam >> theJetDefPt >> plainNLO
     >> theRjj_min >> theRapidityGap >> theopposite_dir >> themjj_min
     >> theDelYhj_min >> thenjets_min >> thenjets_ex; 

}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOAnalysis,AnalysisHandler>
  describeThePEGNLOAnalysis("Herwig::NLOAnalysis", "VBFAnalysis.so");

void NLOAnalysis::Init() {

  typedef bool (arnold::NLOAnalysis::*IGFN)() const;
  typedef void (arnold::NLOAnalysis::*ISFN)(bool);
  typedef Energy (arnold::NLOAnalysis::*IGFNK)() const;
  typedef void (arnold::NLOAnalysis::*ISFNK)(Energy);

  static ClassDocumentation<NLOAnalysis> documentation
    ("There is no documentation for the NLOAnalysis class");

  static Parameter<NLOAnalysis,double> interfaceRparam
    ("Rparam",
     "The R parameter of the jet definition algorithm",
     &NLOAnalysis::theRparam, 0.7, 0.0, 2.0, true, false, true);

  static Parameter<NLOAnalysis,double> interfaceEtaDetector
    ("etaDetector",
     "The detector range",
     &NLOAnalysis::theEtaDetector, 5, 0.0, 100.0, true, false, true);

  static Parameter<NLOAnalysis,double> interfaceJetDefPt
    ("JetDefPt",
     "The minimum allowed pt for a jet",
     &NLOAnalysis::theJetDefPt, 20, 0, 14000.0, true, false, true);

  static Switch<NLOAnalysis,bool> interfacePlainNLO
    ("Mode",
     "Change from plain NLO mode to shower mode",
     &NLOAnalysis::plainNLO, true, true, false);
  static SwitchOption interfaceplainNLOOn
    (interfacePlainNLO,
     "PlainNLO",
     "Analyze a plain NLO calculation",
     true);
  static SwitchOption interfaceplainNLOOff
    (interfacePlainNLO,
     "Shower",
     "Analyze a showered calculation",
     false);

  static Parameter<NLOAnalysis,unsigned int> interfacenjets_min
    ("njets_min",
     "The minimum number of jets",
     &NLOAnalysis::thenjets_min, 2, 1, 10, true, false, true);

  static Parameter<NLOAnalysis,unsigned int> interfacenjets_ex
    ("njets_ex",
     "The exact number of jets",
     &NLOAnalysis::thenjets_ex, 0, 0, 10, true, false, true);

  static Parameter<NLOAnalysis,double> interfaceRjj_min
    ("Rjj_min",
     "The minimum lego plot separation of any jets in the detector",
     &NLOAnalysis::theRjj_min, 0.6, 0.0, 10.0, true, false, true);

  static Parameter<NLOAnalysis,double> interfacedeltay_min
    ("RapidityGap",
     "The distance in pseudorapidity between the two tagging jets",
     &NLOAnalysis::theRapidityGap, 0.0, 0.0, 20.0, true, false, true);

  static Parameter<NLOAnalysis,bool> interfaceopposite_dir
    ("opposite_dir",
     "Should the tagging jets reside in opposite directions (1/0)? ",
     &NLOAnalysis::theopposite_dir, false, false, true, true, false, false,
     (ISFN)0, (IGFN)0,  (IGFN)0, (IGFN)0, (IGFN)0);

  static Parameter<NLOAnalysis,double> interfacemjj_min
    ("mjj_min",
     "The minimum dijet invariant mass of the two tagging jets",
     &NLOAnalysis::themjj_min, 0, 0, 14000, true, false, true);

  static Parameter<NLOAnalysis,double> interfaceDelYhj_min
    ("Yhj_min",
     "The minimum rapidity separation of Higgs to closest jet",
     &NLOAnalysis::theDelYhj_min, 0.0, 0.0, 10.0, true, false, true);


}

