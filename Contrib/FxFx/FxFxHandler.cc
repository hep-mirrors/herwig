#include "FxFxHandler.h"
#include "FxFxReader.h"
#include "FxFxReader.fh"
#include "FxFxEventHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Shower/QTilde/Base/PartnerFinder.h"
#include "Herwig/PDF/HwRemDecayer.h"
#include <queue>
#include "ThePEG/Utilities/Throw.h"
#include "Herwig/Shower/QTilde/Base/KinematicsReconstructor.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include <boost/algorithm/string.hpp>


using namespace Herwig;
using namespace ThePEG;

bool recordEntry(PPtr i,PPtr j) {
  return (i->number()<j->number());
}
bool pTsortFunction(PPtr i,PPtr j) {
  return (i->momentum().perp2()>j->momentum().perp2());
}
bool ETsortFunction(pair<Energy, Lorentz5Momentum> i,pair<Energy, Lorentz5Momentum> j) {
  return (i.first>j.first);
}
bool isMomLessThanEpsilon(Lorentz5Momentum p,Energy epsilon) {
  return (abs(p.x())<epsilon&&abs(p.y())<epsilon&&
	  abs(p.z())<epsilon&&abs(p.t())<epsilon);
}

FxFxHandler::FxFxHandler()
  : ncy_(100),ncphi_(60),ihvy_(-999),nph_(-999),nh_(-999),
    etclusmean_(20*GeV),rclus_(0.4),etaclmax_(5.0),rclusfactor_(1.5),
    ihrd_(-999),njets_(-999),drjmin_(-999), highestMultiplicity_(false),
    ycmax_(5.4),ycmin_(-5.4),jetAlgorithm_(1),vetoIsTurnedOff_(false),vetoSoftThanMatched_(false), etclusfixed_(true),epsetclus_(2.5*GeV), mergemode_(0), vetoHeavyQ_(true), hpdetect_(true), vetoHeavyFlavour_(false)
  {}

void FxFxHandler::doinitrun() {
  QTildeShowerHandler::doinitrun();
  // et_ holds the ET deposited in the (ncy_ x ncphi_) calorimeter cells.
  et_.resize(ncy_);
  for(unsigned int ixx=0; ixx<et_.size(); ixx++) et_[ixx].resize(ncphi_);
  // jetIdx_ for a given calorimeter cell this holds the index of the jet
  // that the cell was clustered into.
  jetIdx_.resize(ncy_);
  for(unsigned int ixx=0; ixx<jetIdx_.size(); ixx++) jetIdx_[ixx].resize(ncphi_);
  if(mergemode_ == 0) {
    cout << "Merging mode is FxFx." << endl;
  } else if (mergemode_ == 1) {
    cout << "Merging mode is Tree level." << endl;
  } else if (mergemode_ == 2) {
    cout << "Merging mode is Tree level, using MadGraph Les Houches information." << endl;
  }
  if(hpdetect_ && mergemode_ != 2) {
    cout << "Automatic detection of hard process enabled." << endl;
    ihrd_ = -999;
  }
}

IBPtr FxFxHandler::clone() const {
  return new_ptr(*this);
}

IBPtr FxFxHandler::fullclone() const {
  return new_ptr(*this);
}

void FxFxHandler::persistentOutput(PersistentOStream & os) const {
  os  << alphaS_
      << ncy_ << ncphi_ << ihvy_ << nph_ << nh_
      << ounit(etclusmean_,GeV) << rclus_ << etaclmax_ << rclusfactor_
      << ihrd_ << njets_ << drjmin_ << highestMultiplicity_
      << ycmax_ << ycmin_ << jetAlgorithm_ << vetoIsTurnedOff_ << vetoSoftThanMatched_ << etclusfixed_
      << cphcal_ << sphcal_ << cthcal_ << sthcal_ << ounit(epsetclus_,GeV) << vetoHeavyQ_ << mergemode_ << hpdetect_ << vetoHeavyFlavour_;
}

void FxFxHandler::persistentInput(PersistentIStream & is, int) {
  is  >> alphaS_
      >> ncy_ >> ncphi_ >> ihvy_ >> nph_ >> nh_
      >> iunit(etclusmean_,GeV) >> rclus_ >> etaclmax_ >> rclusfactor_
      >> ihrd_ >> njets_ >> drjmin_ >> highestMultiplicity_
      >> ycmax_ >> ycmin_ >> jetAlgorithm_ >> vetoIsTurnedOff_ >> vetoSoftThanMatched_ >> etclusfixed_
      >> cphcal_ >> sphcal_ >> cthcal_ >> sthcal_ >> iunit(epsetclus_,GeV) >> vetoHeavyQ_ >> mergemode_ >> hpdetect_ >> vetoHeavyFlavour_;
}

ClassDescription<FxFxHandler> FxFxHandler::initFxFxHandler;
// Definition of the static class description member.

void FxFxHandler::Init() {

  static ClassDocumentation<FxFxHandler> documentation
    ("The FxFxHandler class performs MEPS merging "
     "using the MLM procedure.");

  static Reference<FxFxHandler,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &FxFxHandler::alphaS_, false, false, true, false, false);

  static Parameter<FxFxHandler,int> interfaceihvy
    ("ihvy",
     "heavy flavour in WQQ,ZQQ,2Q etc (4=c, 5=b, 6=t)",
     &FxFxHandler::ihvy_, -999, -999, 7,
     false, false, Interface::limited);

  static Parameter<FxFxHandler,int> interfacenph
    ("nph",
     "Number of photons in the AlpGen process",
     &FxFxHandler::nph_, -999, -999, 7,
     false, false, Interface::limited);

  static Parameter<FxFxHandler,int> interfacenh
    ("nh",
     "Number of higgses in the AlpGen process",
     &FxFxHandler::nph_, -999, -999, 7,
     false, false, Interface::limited);

  static Parameter<FxFxHandler,Energy> interfaceETClus
    ("ETClus",
     "The ET threshold defining a jet in the merging procedure",
     &FxFxHandler::etclusmean_, GeV, 20*GeV, 0*GeV, 14000*GeV,
     false, false, Interface::limited);

  static Parameter<FxFxHandler,double> interfaceRClus
    ("RClus",
     "The cone size used to define a jet in the merging procedure",
     &FxFxHandler::rclus_, 0.4, 0.0, 4.0,
     false, false, Interface::limited);

  static Parameter<FxFxHandler,double> interfaceEtaClusMax
    ("EtaClusMax",
     "The maximum |eta| used to define a jet in the merging procedure",
     &FxFxHandler::etaclmax_, 5.0, 0.0, 15.0,
     false, false, Interface::limited);

  static Parameter<FxFxHandler,double> interfaceRClusFactor
    ("RClusFactor",
     "The prefactor for RClus used to define the jet-parton matching "
     "distance",
     &FxFxHandler::rclusfactor_, 1.5, 0.0, 4.0,
     false, false, Interface::limited);

  static Parameter<FxFxHandler,int> interfaceihrd
    ("ihrd",
     "The hard process code",
     &FxFxHandler::ihrd_, -999, 0, 10000,
     false, false, Interface::limited);

  static Parameter<FxFxHandler,int> interfacenjetsmax
    ("njetsmax",
     "The number of light jets in the maximum-multiplicity process",
     &FxFxHandler::njets_, -999, 0, 10000,
     false, false, Interface::limited);

  static Parameter<FxFxHandler,double> interfacedrjmin
    ("drjmin",
     "Mimimum parton-parton R-sep used for generation.",
     &FxFxHandler::drjmin_, 0.7, 0.0, 4.0,
     false, false, Interface::limited);

  static Parameter<FxFxHandler,bool> interfacehighestMultiplicity
    ("highestMultiplicity",
     "If true it indicates that this is the highest multiplicity input "
     "ME-level configuration to be processed.",
     &FxFxHandler::highestMultiplicity_, 0, 0, 1,
     false, false, Interface::limited);

  static Parameter<FxFxHandler,bool> interfaceETClusFixed
    ("ETClusFixed",
     "If false, indicates that the jet merging scale, etclus_ is allowed to vary"
     "according to epsetclus_",
     &FxFxHandler::etclusfixed_, 1, 0, 1,
     false, false, Interface::limited);

 static Parameter<FxFxHandler,Energy> interfaceEpsilonETClus
    ("EpsilonETClus",
     "The ET threshold defining a jet in the merging procedure",
     &FxFxHandler::epsetclus_, GeV, 2.5*GeV, 0*GeV, 100.0*GeV,
     false, false, Interface::limited);

  static Switch<FxFxHandler,int> interfaceJetAlgorithm
    ("JetAlgorithm",
     "Determines the jet algorithm for finding jets in parton-jet "
     "matching in the MLM procedure.",
     &FxFxHandler::jetAlgorithm_, 1, false, false);
  static SwitchOption AntiKt
    (interfaceJetAlgorithm,
     "AntiKt",
     "The anti-kt jet algorithm.",
     -1);
  static SwitchOption CambridgeAachen
    (interfaceJetAlgorithm,
     "CambridgeAachen",
     "The Cambridge-Aachen jet algorithm.",
     0);
  static SwitchOption Kt
    (interfaceJetAlgorithm,
     "Kt",
     "The Kt jet algorithm.",
     1);

   static Switch<FxFxHandler,int> interfaceMergeMode
    ("MergeMode",
     "The choice of merging mode",
     &FxFxHandler::mergemode_, 0, false, false);
  static SwitchOption FxFx
    (interfaceMergeMode,
     "FxFx",
     "FxFx merging.",
     0);
  static SwitchOption Tree
    (interfaceMergeMode,
     "Tree",
     "Tree-level merging.",
     1);
  static SwitchOption TreeMG
    (interfaceMergeMode,
     "TreeMG5",
     "Tree-level merging using the MadGraph pt clustering information.",
     2);


     static Switch<FxFxHandler,bool> interfaceHardProcessDetection
    ("HardProcessDetection",
     "The choice of merging mode",
     &FxFxHandler::hpdetect_, true, false, false);
  static SwitchOption Automatic
    (interfaceHardProcessDetection,
     "Automatic",
     "Automatically determine which particles to include in the merging.",
     true);
  static SwitchOption Manual
    (interfaceHardProcessDetection,
     "Manual",
     "Use the ihrd code to determine which particles to include in the merging.",
     false);

  static Switch<FxFxHandler,bool> interfaceVetoIsTurnedOff
    ("VetoIsTurnedOff",
     "Allows the vetoing mechanism to be switched off.",
     &FxFxHandler::vetoIsTurnedOff_, false, false, false);
  static SwitchOption VetoingIsOn
    (interfaceVetoIsTurnedOff,
     "VetoingIsOn",
     "The MLM merging veto mechanism is switched ON.",
     false);
  static SwitchOption VetoingIsOff
    (interfaceVetoIsTurnedOff,
     "VetoingIsOff",
     "The MLM merging veto mechanism is switched OFF.",
     true);


  static Switch<FxFxHandler,bool> interfaceVetoHeavyFlavour
    ("VetoHeavyFlavour",
     "Allows the heavy flavour vetoing mechanism to be switched off.",
     &FxFxHandler::vetoHeavyFlavour_, false, false, false);
  static SwitchOption HeavyFVetoingIsOn
    (interfaceVetoHeavyFlavour,
     "Yes",
     "The MLM merging veto mechanism for heavy flavour is switched ON.",
     true);
  static SwitchOption HeavyFVetoingIsOff
    (interfaceVetoHeavyFlavour,
     "No",
     "The MLM merging veto mechanism for heavy flavour is switched OFF.",
     false);

    static Switch<FxFxHandler,bool> interfaceHeavyQVeto
    ("HeavyQVeto",
     "Allows the vetoing mechanism on the heavy quark products to be switched off.",
     &FxFxHandler::vetoHeavyQ_, false, false, false);
  static SwitchOption HQVetoingIsOn
    (interfaceHeavyQVeto,
     "Yes",
     "The MLM merging veto on Heavy quark decay produts mechanism is switched ON.",
     true);
  static SwitchOption HQVetoingIsOff
    (interfaceHeavyQVeto,
     "No",
     "The MLM merging veto on Heavy quark decay products mechanism is switched OFF.",
     false);

    static Switch<FxFxHandler,bool> interfaceVetoSoftThanMatched
    ("VetoSoftThanMatched",
     "Allows the vetoing mechanism to be switched off.",
     &FxFxHandler::vetoSoftThanMatched_, false, false, false);
  static SwitchOption VetoSoftIsOn
    (interfaceVetoSoftThanMatched,
     "VetoSoftIsOn",
     "The vetoing of highest-mult. events with jets softer than matched ones is ON",
     true);
  static SwitchOption VetoSoftIsOff
    (interfaceVetoSoftThanMatched,
     "VetoSoftIsOff",
     "The vetoing of highest-mult. events with jets softer than matched ones is OFF.",
     false);


}

void FxFxHandler::dofinish() {
  QTildeShowerHandler::dofinish();
}

void FxFxHandler::doinit() {

  //print error if HardProcID is not set in input file
  if(ihrd_ == -999 && !hpdetect_) { cout << "Error: FxFxHandler:ihrd not set and FxFx:HardProcessDetection set to Manual!" << endl; exit(1); }
  QTildeShowerHandler::doinit();

  // Compute calorimeter edges in rapidity for GetJet algorithm.
  ycmax_=etaclmax_+rclus_;
  ycmin_=-ycmax_;
}

// Throws a veto according to MLM strategy ... when we finish writing it.
bool FxFxHandler::showerHardProcessVeto() const {
  int debug_mode = 0;
  if(vetoIsTurnedOff_) {
    //    cout << "Vetoing is turned OFF." << endl;
    return false;
  }

  //if(debug_mode) { cout << "debug_mode = " << 5 << endl; } 

  // Skip veto for processes in which merging is not implemented:
  if(ihrd_==7||ihrd_==8||ihrd_==13) {
      ostringstream wstring;
      wstring << "FxFxHandler::showerHardProcessVeto() - warning."
	      << "MLM merging not implemented "
	      << "processes 4Q (ihrd=7), QQh (ihrd=8), "
	      << "(single) top (ihrd=13) \n";
      generator()->logWarning( Exception(wstring.str(), 
                                         Exception::warning) );
      return false;
  }

  // Fill preshowerISPs_ pair and preshowerFSPs_ particle pointer vector. 
  getPreshowerParticles();

  // Fill showeredISHs_, showeredISPs and showeredRems pairs, as well as
  // showeredFSPs_ particle pointer vector. 
  getShoweredParticles();

  // Turn on some screen output debugging: 0 = none ---> 5 = very verbose.
  doSanityChecks(debug_mode);

   
  // Dimensions of each calorimter cell in y and phi.
  dely_ = (ycmax_-ycmin_)/double(ncy_);
  delphi_ = 2*M_PI/double(ncphi_);
  
  // Fill partonsToMatch_ with only those pre-shower partons intended to
  // used in jet-parton matching and fill particlesToCluster_ using only
  // those final state particles (post-shower) which are supposed to go
  // in the jet clustering used to do merging.
  partonsToMatch_     = preshowerFSPs_;
  particlesToCluster_ = showeredFSPs_ ; 

  // Filter out all but the 'extra' light-parton progenitors and their
  // associated final state particles.
  if(mergemode_ == 0 || mergemode_ == 1) { caldel_m(); }
  else if(mergemode_ == 2) { caldel_mg(); }  
 
  double prob(1);
  //if etclusfixed_ then set the etclus_ to the fixed chosen value
  if(etclusfixed_) { 
    etclus_ = etclusmean_;
  } else {
    //else, if we wish to vary etclus_, we use the probability distribution
    //choose a probability between 0 and 1
    prob = rnd();
    etclus_ = etclusran_(prob);
  }  
  
  // Cluster particlesToCluster_ into jets with FastJet.
  getFastJets(rclus_,etclus_,etaclmax_);

  if(mergemode_ == 0) { 
    
    // Get npLO_ and npNLO_ for FxFx matching
    getnpFxFx();

    // print the npXLO_ values obtained 
    //  cout << "HANDLER:\t\t\t\t" << npLO_ << "\t\t" << npNLO_ << endl;

    //FxFx modifications start here. 

    // Sort partonsToMatch_ from high to low pT.
    sort(partonsToMatch_.begin(),partonsToMatch_.end(),pTsortFunction);

    // Count the number of jets.
    int njets_found(pjet_.size());
  
    // If the number of jets found is not equal to the number of partons in the Born 
    // (i.e., the number of partons in the S-event, or one less than the number of partons in an H-event), 
    // the jets cannot be matched and the event has to be rejected. The number of partons in the Born is written in the event file with a name “npNLO” 
    // if there are no jets to match and no jets have been found, do not veto. 
    if(njets_found == 0 && npNLO_ == 0)  {  /*cout << "njets_found = " << njets_found << " and npNLO = " << npNLO_ << ", accepting" << endl;*/ return false; }

    //if the number of jets is smaller than npNLO -> reject the event.
    if(njets_found < npNLO_) { /*cout << "njets_found = " << njets_found << " and npNLO = " << npNLO_ << ", rejecting" << endl;*/ return true; }
    // For the maximum-multiplicity sample, the number of jets obtained does not have to be exactly equal to npNLO, it may also be larger; 
    if(njets_found > npNLO_ && npNLO_ != njets_) { /*cout << "njets_found = " << njets_found << " and npNLO = " << npNLO_ << ", rejecting" << endl;*/ return true; }

    // Create the matrix-element jets.
    // Cluster also the partons at the hard-matrix element level into jets with the same algorithm as above, 
    // but without the requirement of a minimal pT on the jets (or set it very small). 
    // By construction, for S-events you should find exactly npNLO jets, while for the H-events it is either npNLO or npNLO+1.
    // Cluster partonsToMatch_ into jets with FastJet.
    getFastJetsToMatch(rclus_,0*GeV,etaclmax_);
  
    int me_njets_found(pjetME_.size());
  
    // cout << "number of ME jets found = " << me_njets_found << "partons to match: " << partonsToMatch_.size() << endl;

    // Match light progenitors to jets.
    vector<int> jetToPartonMap(pjetME_.size(),-999);
    Energy etmin(777e100*GeV);

    // Match the jets.
    // Try to match the “npNLO” hardest jets created post-shower with any of the jets pre-shower. Two jets are matched if the distance between them is smaller than 1.5*DeltaR. 
    // If not all the npNLO hardest shower jets are matched the event has to be rejected.
    // Note that if the current event does not belong to the maximum multiplicity sample, this means that all the shower jets need to be matched, because the requirement above already rejects 
    // events that do not have npNLO shower jets.
    // For those events, at the level of the matrix elements there can either be npNLO or npNLO+1 matrix-element jets, depending on S- or H-events and the kinematics of those partons.
    // Still only the shower jets need to be matched, so an event should not be rejected if a matrix-element jet cannot be matched.
    // For each parton, starting with the hardest one ...
    for(unsigned int ixx=0; ixx<npNLO_; ixx++) {
      // ... loop over all jets not already matched.
      double DRmin(777e100);
      int    jetIndexForDRmin(-999);
      for(unsigned int jxx=0; jxx<pjetME_.size(); jxx++) {
        // ... and tag closest of the remaining ones
        double DRpartonJet(partonJetDeltaR(pjetME_[ixx],pjet_[jxx]));
        if(jetToPartonMap[jxx]<0&&DRpartonJet<DRmin) {
          DRmin=DRpartonJet;
          jetIndexForDRmin=jxx;
        }
      }
      // If the parton-jet distance is less than the matching
      // distance, the parton and jet match.
      if(DRmin<rclus_*rclusfactor_&&jetIndexForDRmin>=0) {
        jetToPartonMap[jetIndexForDRmin]=ixx;
        if(ixx==0||etjet_[jetIndexForDRmin]<etmin)
          etmin=etjet_[jetIndexForDRmin]; 
        // Otherwise this parton is not matched so veto the event.
      } else return true;
    }

    // Veto events where matched jets are softer than non-matched ones,
    // in the inclusive (highestMultiplicity_ = true) mode, unless we
    // are dealing with NLO input events.
    if(npNLO_ == njets_ && vetoSoftThanMatched_) {
      //cout << "highest mult. event being tested for softer jets than matched..." << endl;
      for(unsigned int iyy=0; iyy<pjet_.size(); iyy++) {
        //cout << "etjet_[iyy] = " << etjet_[iyy]/GeV << endl;
        if(jetToPartonMap[iyy]<0&&etmin<etjet_[iyy]) { /*cout << "VETO!" << endl;*/ return true; }
      }
    }

  
    if(!vetoHeavyQ_) {
      //cout << "no heavy quark decay product veto!" << endl;
      return false;
    }
  
    //end of FxFx part
    // **************************************************************** //
    // * Now look to the non-light partons for heavy quark processes. * //
    // **************************************************************** //

    if( (ihrd_<=2||ihrd_==6||ihrd_==10||ihrd_==15||ihrd_==16) || (hpdetect_==true && hvqfound==true) ) {
      // Extract heavy quark progenitors and the radiation they
      // produce and put it in the calorimeter.
      caldel_hvq();
    
      // Cluster particlesToCluster_ into jets with FastJet.
      getFastJets(rclus_,etclus_,etaclmax_);
    
      // If the radiation from the heavy quarks does not give rise
      // to any jets we accept event.
      if(pjet_.size() == 0) return false;

      // If extra jets emerge from the jet clustering we only
      // accept events where the jets formed by radiation from
      // b and c quarks lies within drjmin_ of the heavy quark
      // progenitor.
      int nmjet(pjet_.size());
      for(unsigned int ixx=0; ixx<pjet_.size(); ixx++) {
        for(unsigned int jxx=0; jxx<partonsToMatch_.size(); jxx++) {
          if(!(abs(partonsToMatch_[jxx]->id())==4||abs(partonsToMatch_[jxx]->id())==5)) continue;
          if(partonJetDeltaR(partonsToMatch_[jxx],pjet_[ixx])<drjmin_) {
            nmjet--;           // Decrease the number of unmatched jets.
            etjet_[ixx]=0*GeV; // Set jet ET to zero to indicate it is 'matched'.
          }
        }
      }

      // If every jet matched to _at_least_one_ progenitor accept the event.
      if(nmjet<=0) return false;

    }
  }

   /*****
   *****
   ***** END of mergemode_ == 0 (i.e. FxFx MERGING BLOCK)
   *****/


  
  if(mergemode_ == 1 || mergemode_ == 2) {
    
    // determine whether event is of highest multiplicity or not
    if(partonsToMatch_.size()==njets_) { highestMultiplicity_ = true; }

    //  cout << "jet finding gives pjet_.size() = " << pjet_.size() << endl;
    // If there are less jets than partons then parton-jet matching is
    // bound to fail: reject the event already. Also, if the input is
    // an NLO event file it will 99.5% of the time contain a number of
    // light partons in the F.S. equal to that in the real emission
    // process in the NLO calculation, moreover, it has already
    // effectively merged njets_-1 and njets jet events. So in that
    // case we do not reject events on the grounds they have jet
    // multiplicity less than partonsToMatch_.size() but rather less
    // jets than partonsToMatch.size()-1; such events are better
    // described by the lower-by-one-unit (partonsToMatch_.size()-1)
    // of multiplicity NLO event file, or the lower-by-two-units
    // (partonsToMatch_.size()-2) of multiplicity LO event file. 

    // If it is not jet production apply rejection criterion as above.
    if(ihrd_!=9) {
      if(pjet_.size() < partonsToMatch_.size()) return true;
    }
    // Otherwise, in the case of jet production allow the lowest
    // contributing multiplicity process (just at NLO), namely,
    // dijet production, to give rise to 1-jet and even 0-jet 
    // events, since these can contribute to, for example, the
    // inclusive jet cross section i.e. in this case the rejection
    // is only applied in the case of the next-to-lowest multiplicity
    // processes (>2 parton events at LO and >3 parton events at NLO).
    else {
      // KH - March 5th
      // Removed the following line giving special treatment
      // also to the LO events, to maintain consistency with
      // the fortran algorithm, at least for now. So now jet
      // production at LO is being treated the same as all 
      // other processes.
      //      if(partonsToMatch_.size()==2 && pjet_.size()<2) return false;
      if(pjet_.size() < partonsToMatch_.size()) return true;
    }
    
    // Sort partonsToMatch_ from high to low pT.
    sort(partonsToMatch_.begin(),partonsToMatch_.end(),pTsortFunction);
    
    // Match light progenitors to jets.
    vector<int> jetToPartonMap(pjet_.size(),-999);
    Energy etmin(777e100*GeV);
  
    // For each parton, starting with the hardest one ...
    for(unsigned int ixx=0; ixx<partonsToMatch_.size(); ixx++) {
      // ... loop over all jets not already matched.
      double DRmin(777e100);
      int    jetIndexForDRmin(-999);
      for(unsigned int jxx=0; jxx<pjet_.size(); jxx++) {
	// ... and tag closest of the remaining ones
	double DRpartonJet(partonJetDeltaR(partonsToMatch_[ixx],pjet_[jxx]));
	if(jetToPartonMap[jxx]<0&&DRpartonJet<DRmin) {
	  DRmin=DRpartonJet;
	  jetIndexForDRmin=jxx;
	}
      }
      // If the parton-jet distance is less than the matching
      // distance, the parton and jet match.
      if(DRmin<rclus_*rclusfactor_&&jetIndexForDRmin>=0) {
	jetToPartonMap[jetIndexForDRmin]=ixx;
	if(ixx==0||etjet_[jetIndexForDRmin]<etmin)
	  etmin=etjet_[jetIndexForDRmin]; 
	// Otherwise this parton is not matched so veto the event.
      } else return true;
    }


    // Veto events with larger jet multiplicity from exclusive sample.
    if(!highestMultiplicity_&&pjet_.size()>partonsToMatch_.size()) return true; 
 
    // Veto events where matched jets are softer than non-matched ones,
    // in the inclusive (highestMultiplicity_ = true) mode, unless we
    // are dealing with NLO input events.
    if(highestMultiplicity_) {
      for(unsigned int ixx=0; ixx<pjet_.size(); ixx++) 
	if(jetToPartonMap[ixx]<0&&etmin<etjet_[ixx]) return true;
    }

    if(!vetoHeavyQ_) {
      //cout << "no heavy quark decay product veto!" << endl;
      return false;
    }
  
    
    // **************************************************************** //
    // * Now look to the non-light partons for heavy quark processes. * //
    // **************************************************************** //
  
    if( (ihrd_<=2||ihrd_==6||ihrd_==10||ihrd_==15||ihrd_==16) || (hpdetect_==true && hvqfound==true))  {
  
      // Extract heavy quark progenitors and the radiation they
      // produce and put it in the calorimeter.
      caldel_hvq();
    
      // Cluster particlesToCluster_ into jets with FastJet.
      getFastJets(rclus_,etclus_,etaclmax_);
    
      // If the radiation from the heavy quarks does not give rise
      // to any jets we accept event.
      if(pjet_.size() == 0) return false;

      // If extra jets emerge from the jet clustering we only
      // accept events where the jets formed by radiation from
      // b and c quarks lies within drjmin_ of the heavy quark
      // progenitor.
      int nmjet(pjet_.size());
      for(unsigned int ixx=0; ixx<pjet_.size(); ixx++) {
	for(unsigned int jxx=0; jxx<partonsToMatch_.size(); jxx++) {
	  if(!(abs(partonsToMatch_[jxx]->id())==4||abs(partonsToMatch_[jxx]->id())==5)) continue;
	  if(partonJetDeltaR(partonsToMatch_[jxx],pjet_[ixx])<drjmin_) {
	    nmjet--;           // Decrease the number of unmatched jets.
	    etjet_[ixx]=0*GeV; // Set jet ET to zero to indicate it is 'matched'.
	  }
	}
      }

      // If every jet matched to _at_least_one_ progenitor accept the event.
      if(nmjet<=0) return false;
      else {
	// If unmatched jets remain, reject the event if highestMultiplicity_!=1
	if(!highestMultiplicity_) return true;
	else {
	  // If unmatched jets remain and highestMultiplicity is true then check
	  // that these are softer than all the matched ones (from the light-parton
	  // matching round).
	  Energy etmax(0.*GeV);
	  for(unsigned int ixx=0; ixx<pjet_.size(); ixx++) etmax=max(etjet_[ixx],etmax);
	  if(etmax>etmin) return true;
	}
      }
    }

  }

  
  // Otherwise we accept the event ...
  return false;

}

/* Function that returns the R distance
   between a particle and a jet. */
double FxFxHandler::partonJetDeltaR(ThePEG::tPPtr partonptr, LorentzMomentum jetmom) const { 
  LorentzMomentum partonmom(partonptr->momentum());
  // Calculate DY, DPhi and then DR
  double DY(partonmom.eta()-jetmom.eta());
  double DPhi(partonmom.phi()-jetmom.phi());
  if(DPhi>M_PI) DPhi=2*M_PI-DPhi;
  double DR(sqrt(sqr(DY)+sqr(DPhi)));
  return DR;
}

double FxFxHandler::partonJetDeltaR(LorentzMomentum jetmom1, LorentzMomentum jetmom2) const { 
  // Calculate DY, DPhi and then DR
  double DY(jetmom1.eta()-jetmom2.eta());
  double DPhi(jetmom1.phi()-jetmom2.phi());
  if(DPhi>M_PI) DPhi=2*M_PI-DPhi;
  double DR(sqrt(sqr(DY)+sqr(DPhi)));
  return DR;
}

// Get FastJets
void FxFxHandler::getFastJets(double rjet, Energy ejcut, double etajcut) const {

  vector<fastjet::PseudoJet> particlesToCluster;
  for(unsigned int ipar=0; ipar<particlesToCluster_.size(); ipar++) { 
    double y(particlesToCluster_[ipar]->momentum().eta());
    if(y>=ycmin_&&y<=ycmax_) { 
      int absId(abs(particlesToCluster_[ipar]->id()));
      // If it's not a lepton / top / photon it may go in the jet finder.
      if(!(absId>=11&&absId<=16) && absId!=6 && absId!=22) { 
	// input particles into fastjet pseudojet
	fastjet::PseudoJet p(particlesToCluster_[ipar]->momentum().x()/GeV,
			     particlesToCluster_[ipar]->momentum().y()/GeV,
			     particlesToCluster_[ipar]->momentum().z()/GeV,
			     particlesToCluster_[ipar]->momentum().e()/GeV);
	p.set_user_index(ipar);
	particlesToCluster.push_back(p);
      }
    }
  }

  fastjet::RecombinationScheme recombinationScheme = fastjet::E_scheme;
  fastjet::Strategy            strategy            = fastjet::Best;
  double R(rjet);
  fastjet::JetDefinition theJetDefinition;
  switch (jetAlgorithm_) {
  case  -1: theJetDefinition=fastjet::JetDefinition(fastjet::antikt_algorithm,
						    R,
						    recombinationScheme,
						    strategy); break;
  case   0: theJetDefinition=fastjet::JetDefinition(fastjet::cambridge_algorithm,
						    R,
						    recombinationScheme,
						    strategy); break;
  case   1: theJetDefinition=fastjet::JetDefinition(fastjet::kt_algorithm,
						    R,
						    recombinationScheme,
						    strategy); break;
  default:  theJetDefinition=fastjet::JetDefinition(fastjet::cambridge_algorithm,
						    R,
						    recombinationScheme,
						    strategy); break;
  }
  fastjet::ClusterSequence fastjetEvent(particlesToCluster,theJetDefinition);
  vector<fastjet::PseudoJet> inclusiveJets = fastjetEvent.inclusive_jets();
  inclusiveJets = fastjet::sorted_by_pt(inclusiveJets);

  // Fill the array of jet momenta for the rest of the veto procedure.
  pjet_.clear();
  pjet_.resize(inclusiveJets.size());
  etjet_.clear();
  etjet_.resize(inclusiveJets.size());
  for(unsigned int ffj=0; ffj<pjet_.size();ffj++) {
    pjet_[ffj]  = Lorentz5Momentum(inclusiveJets[ffj].px()*GeV,
				   inclusiveJets[ffj].py()*GeV,
				   inclusiveJets[ffj].pz()*GeV,
				   inclusiveJets[ffj].e()*GeV);
    pjet_[ffj].rescaleMass();
    etjet_[ffj] = pjet_[ffj].et();
  }

  // Throw the jet away if it's outside required eta region or
  // has transverse energy below ejcut.
  for(unsigned int fj=0; fj<pjet_.size(); fj++)
    if(etjet_[fj]<ejcut||fabs(pjet_[fj].eta())>etajcut) {
      pjet_.erase(pjet_.begin()+fj);
      etjet_.erase(etjet_.begin()+fj);
      fj--;
    }

  // Sort jets from high to low ET.
  vector<pair<Energy, Lorentz5Momentum> > etjet_pjet;
  for(unsigned int ixx=0; ixx<etjet_.size(); ixx++)
    etjet_pjet.push_back(make_pair(etjet_[ixx],pjet_[ixx]));
  sort(etjet_pjet.begin(),etjet_pjet.end(),ETsortFunction);
  for(unsigned int ixx=0; ixx<etjet_.size(); ixx++) {
    etjet_[ixx]=etjet_pjet[ixx].first;
    pjet_[ixx]=etjet_pjet[ixx].second;
  }

  return;
}

// Get FastJets from partonsToMatch_
void FxFxHandler::getFastJetsToMatch(double rjet, Energy ejcut, double etajcut) const {

  vector<fastjet::PseudoJet> particlesToCluster;
  for(unsigned int ipar=0; ipar<partonsToMatch_.size(); ipar++) { 
    double y(partonsToMatch_[ipar]->momentum().eta());
    if(y>=ycmin_&&y<=ycmax_) { 
      int absId(abs(partonsToMatch_[ipar]->id()));
      // If it's not a lepton / top / photon it may go in the jet finder.
      if(!(absId>=11&&absId<=16) && absId!=6 && absId!=22) { 
	// input particles into fastjet pseudojet
	fastjet::PseudoJet p(partonsToMatch_[ipar]->momentum().x()/GeV,
			     partonsToMatch_[ipar]->momentum().y()/GeV,
			     partonsToMatch_[ipar]->momentum().z()/GeV,
			     partonsToMatch_[ipar]->momentum().e()/GeV);
	p.set_user_index(ipar);
	particlesToCluster.push_back(p);
      }
    }
  }

  fastjet::RecombinationScheme recombinationScheme = fastjet::E_scheme;
  fastjet::Strategy            strategy            = fastjet::Best;
  double R(rjet);
  fastjet::JetDefinition theJetDefinition;
  switch (jetAlgorithm_) {
  case  -1: theJetDefinition=fastjet::JetDefinition(fastjet::antikt_algorithm,
						    R,
						    recombinationScheme,
						    strategy); break;
  case   0: theJetDefinition=fastjet::JetDefinition(fastjet::cambridge_algorithm,
						    R,
						    recombinationScheme,
						    strategy); break;
  case   1: theJetDefinition=fastjet::JetDefinition(fastjet::kt_algorithm,
						    R,
						    recombinationScheme,
						    strategy); break;
  default:  theJetDefinition=fastjet::JetDefinition(fastjet::cambridge_algorithm,
						    R,
						    recombinationScheme,
						    strategy); break;
  }
  fastjet::ClusterSequence fastjetEvent(particlesToCluster,theJetDefinition);
  vector<fastjet::PseudoJet> inclusiveJets = fastjetEvent.inclusive_jets();
  inclusiveJets = fastjet::sorted_by_pt(inclusiveJets);
  
  // Fill the array of jet momenta for the rest of the veto procedure.
  pjetME_.clear();
  pjetME_.resize(inclusiveJets.size());
  for(unsigned int ffj=0; ffj<pjetME_.size();ffj++) {
    pjetME_[ffj]  = Lorentz5Momentum(inclusiveJets[ffj].px()*GeV,
                                     inclusiveJets[ffj].py()*GeV,
                                     inclusiveJets[ffj].pz()*GeV,
                                     inclusiveJets[ffj].e()*GeV);
    pjetME_[ffj].rescaleMass();
  }

  return;
}

// Deletes particles from partonsToMatch_ and particlesToCluster_
// vectors so that these contain only the partons to match to the
// jets and the particles used to build jets respectively. By and
// large the candidates for deletion are: vector bosons and their
// decay products, Higgs bosons, photons as well as _primary_, i.e.
// present in the lowest multiplicity process, heavy quarks and
// any related decay products.
void FxFxHandler::getDescendents(PPtr theParticle) const {
  ParticleVector theChildren(theParticle->children());
  for (unsigned int ixx=0; ixx<theChildren.size(); ixx++)
    if(theChildren[ixx]->children().size()==0)
      tmpList_.push_back(theChildren[ixx]);
    else
      getDescendents(theChildren[ixx]);
  return;
}

void FxFxHandler::caldel_m() const {
  
  preshowerFSPsToDelete_.clear();
  showeredFSPsToDelete_.clear();

  hvqfound = false;
  
  if(hpdetect_ && mergemode_!=2) {
    for(unsigned int ixx=0; ixx<preshowerFSPs_.size(); ixx++) {
      tmpList_.clear();
      /* Exclude the top quarks and any children
         they may have produced from the jet parton matching. */

     
      if(abs(preshowerFSPs_[ixx]->parents()[0]->id())==6) {
        hvqfound = true;
        // cout << "preshowerFSPs_[ixx]->id() = " << preshowerFSPs_[ixx]->id() << " " << preshowerFSPs_[ixx]->momentum().perp2() << endl;
        preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
        getDescendents(preshowerFSPs_[ixx]->parents()[0]);
        for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++) {
          // cout << "tmpList_[jxx]->id() = " << tmpList_[jxx]->id() << " " <<  tmpList_[jxx]->momentum().perp2() << endl;
          showeredFSPsToDelete_.push_back(tmpList_[jxx]);
        }
        continue; 
      }
             
      /* Exclude the v.bosons and any children they may
         have produced from the jet parton matching. */
      if( (abs(preshowerFSPs_[ixx]->parents()[0]->id())==23||
         abs(preshowerFSPs_[ixx]->parents()[0]->id())==24||
         abs(preshowerFSPs_[ixx]->id())==22||
           abs(preshowerFSPs_[ixx]->id())==25|| (abs(preshowerFSPs_[ixx]->id()) < 17 && abs(preshowerFSPs_[ixx]->id()) > 10)) && (abs(preshowerFSPs_[ixx]->parents()[0]->id())!=6))  {
        preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
        getDescendents(preshowerFSPs_[ixx]);
        for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
          showeredFSPsToDelete_.push_back(tmpList_[jxx]);
        continue;
      }
      /* Exclude the bottom quarks and any children
         they may have produced from the jet parton matching. */
      if(abs(preshowerFSPs_[ixx]->id())==5&&ixx<2) {
        hvqfound = true;   
        preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
        getDescendents(preshowerFSPs_[ixx]);
        for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
          showeredFSPsToDelete_.push_back(tmpList_[jxx]);
        continue;
      }
    }
  }
  if(!hpdetect_) { 
    for(unsigned int ixx=0; ixx<preshowerFSPs_.size(); ixx++) {
      tmpList_.clear();
      if(ihrd_<=2) {
        /* wqq... , zqq... */
        /* Exclude the heavy quarks and any children they may
           have produced as well as the v.boson and any children
           it may have produced from the jet parton matching. */
        if(abs(preshowerFSPs_[ixx]->parents()[0]->id())==23||
           abs(preshowerFSPs_[ixx]->parents()[0]->id())==24) {
          preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
          getDescendents(preshowerFSPs_[ixx]);
          for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
            showeredFSPsToDelete_.push_back(tmpList_[jxx]);
        }
        if(abs(preshowerFSPs_[ixx]->id())==ihvy_&&ixx<2) {
          preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
          getDescendents(preshowerFSPs_[ixx]);
          for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
            showeredFSPsToDelete_.push_back(tmpList_[jxx]);
        }
        if(ixx==preshowerFSPs_.size()-1&&preshowerFSPsToDelete_.size()!=4) {
          throw Exception()
            << "FxFxHandler::caldel_m() - ERROR!\n"
            << "wqq / zqq process should have 4 particles to omit from"
            << "jet-parton matching for ihvy=" << ihvy_ << "." << Exception::eventerror;
        }
      } else if(ihrd_<=4)  {
        /* zjet... */
        /* Exclude the v.boson and any children
           it may have produced from the jet parton matching. */
        if(abs(preshowerFSPs_[ixx]->parents()[0]->id())==23||
           abs(preshowerFSPs_[ixx]->parents()[0]->id())==24) {
          preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
          getDescendents(preshowerFSPs_[ixx]);
          for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
            showeredFSPsToDelete_.push_back(tmpList_[jxx]);
        }
        /*if(ixx==preshowerFSPs_.size()-1&&preshowerFSPsToDelete_.size()!=2) {
          throw Exception()
          << "FxFxHandler::caldel_m() - ERROR!\n"
          << "zjet process should have 2 particles to omit from"
          << "jet-parton matching." << Exception::eventerror;
          }*/
      } else if(ihrd_==5)  {
        /* vbjet... */
        /* Exclude the v.bosons and any children they may
           have produced from the jet parton matching. */
        if(abs(preshowerFSPs_[ixx]->parents()[0]->id())==23||
           abs(preshowerFSPs_[ixx]->parents()[0]->id())==24||
           abs(preshowerFSPs_[ixx]->id())==22||
           abs(preshowerFSPs_[ixx]->id())==25) {
          preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
          getDescendents(preshowerFSPs_[ixx]);
          for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
            showeredFSPsToDelete_.push_back(tmpList_[jxx]);
        }
      } else if(ihrd_==6)  {
        /* 2Q... */
        /* Exclude the heavy quarks and any children
           they may have produced from the jet parton matching. */
        if(ihvy_==6) {
          if(abs(preshowerFSPs_[ixx]->parents()[0]->id())==6||
             abs(preshowerFSPs_[ixx]->parents()[0]->id())==24) {
         //   cout << "preshowerFSPs_[ixx]->id() = " << preshowerFSPs_[ixx]->id() << " " << preshowerFSPs_[ixx]->momentum().perp2() << endl;
            preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
          }
          if(abs(preshowerFSPs_[ixx]->parents()[0]->id())==6) {
            getDescendents(preshowerFSPs_[ixx]->parents()[0]);
            for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++) {
           //   cout << "tmpList_[jxx]->id() = " << tmpList_[jxx]->id() << " " <<  tmpList_[jxx]->momentum().perp2() << endl;
              showeredFSPsToDelete_.push_back(tmpList_[jxx]);
            }
          }
          if(ixx==preshowerFSPs_.size()-1&&preshowerFSPsToDelete_.size()!=6) {
            throw Exception()
              << "FxFxHandler::caldel_m() - ERROR!\n"
              << "2Q process should have 6 particles to omit from"
              << "jet-parton matching for ihvy=" << ihvy_ << "." << Exception::eventerror;
          }
        } else {
          if(abs(preshowerFSPs_[ixx]->id())==ihvy_&&ixx<2) {
            preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
            getDescendents(preshowerFSPs_[ixx]);
            for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
              showeredFSPsToDelete_.push_back(tmpList_[jxx]);
          }
          if(ixx==preshowerFSPs_.size()-1&&preshowerFSPsToDelete_.size()!=2) {
            throw Exception()
              << "FxFxHandler::caldel_m() - ERROR!\n"
              << "2Q process should have 2 particles to omit from"
              << "jet-parton matching for ihvy=" << ihvy_ << "." << Exception::eventerror;
          }
          }
      } else if(ihrd_==7)  {
          /* 4Q... */
          /* There are no light jets for this process, so nothing to match. */
        } else if(ihrd_==9)  {
          /* Njet... */
        } else if(ihrd_==10) {
          /* wcjet... */
          /* Exclude the charm quark and any children it may 
             have produced as well as the v.boson and any children
             it may have produced from the jet parton matching. */
          if((abs(preshowerFSPs_[ixx]->id())==4&&ixx<1)||
             abs(preshowerFSPs_[ixx]->parents()[0]->id())==24) {
            preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
            getDescendents(preshowerFSPs_[ixx]);
            for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
              showeredFSPsToDelete_.push_back(tmpList_[jxx]);
          }
          if(ixx==preshowerFSPs_.size()-1&&preshowerFSPsToDelete_.size()!=3) {
            throw Exception()
              << "FxFxHandler::caldel_m() - ERROR!\n"
              << "wcjet process should have 3 particles to omit from"
              << "jet-parton matching." << Exception::eventerror;
          }
        } else if(ihrd_==11) {
          /* phjet... */
          /* Exclude the hard photons from the jet parton matching. */
          if(abs(preshowerFSPs_[ixx]->id())==22) {
            preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
            getDescendents(preshowerFSPs_[ixx]);
            for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
              showeredFSPsToDelete_.push_back(tmpList_[jxx]);
          }
          unsigned int tmpUnsignedInt(nph_);
          if(ixx==preshowerFSPs_.size()-1&&preshowerFSPsToDelete_.size()!=tmpUnsignedInt) {
            throw Exception()
              << "FxFxHandler::caldel_m() - ERROR!\n"
              << "phjet process should have " << nph_ << " particles to omit from"
              << "jet-parton matching." << Exception::eventerror;
          }
        } else if(ihrd_==12) {
          /* hjet... */
          /* Exclude the higgs and any children it may have
             produced from the jet parton matching. */
          if(abs(preshowerFSPs_[ixx]->id())==25) {
            preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
            getDescendents(preshowerFSPs_[ixx]);
            for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
              showeredFSPsToDelete_.push_back(tmpList_[jxx]);
          }
          if(ixx==preshowerFSPs_.size()-1&&preshowerFSPsToDelete_.size()!=1) {
            throw Exception()
              << "FxFxHandler::caldel_m() - ERROR!\n"
              << "hjet process should have 1 particle to omit from"
              << "jet-parton matching." << Exception::eventerror;
          }
        } else if(ihrd_==14) {
          /* wphjet... */
          /* Exclude the v.boson and any children it may have
             produced from the jet parton matching. */
          // AND WHAT ABOUT THE PHOTON? <--- CHECK THIS WITH AlpGen GUYs.
          if(abs(preshowerFSPs_[ixx]->id())==22||
             abs(preshowerFSPs_[ixx]->parents()[0]->id())==24) {
            preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
            getDescendents(preshowerFSPs_[ixx]);
            for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
              showeredFSPsToDelete_.push_back(tmpList_[jxx]);
          }
          unsigned int tmpUnsignedInt(2+nph_);
          if(ixx==preshowerFSPs_.size()-1&&preshowerFSPsToDelete_.size()!=tmpUnsignedInt) {
            throw Exception()
              << "FxFxHandler::caldel_m() - ERROR!\n"
              << "wphjet process should have " << 2+nph_ << " particles to omit from"
              << "jet-parton matching." << Exception::eventerror;
          }
        } else if(ihrd_==15) {
          /* wphqq... <--- N.B. if q = top, it is not decayed. */
          /* Exclude the heavy quarks and any children they may
             have produced as well as the v.boson and any children
             it may have produced from the jet parton matching. */
          // AND WHAT ABOUT THE PHOTON? <--- CHECK THIS WITH AlpGen GUYs.
          if(abs(preshowerFSPs_[ixx]->id())==22||
             (abs(preshowerFSPs_[ixx]->id())==ihvy_&&ixx==(preshowerFSPs_.size()-(2+nph_+1)))||
             (abs(preshowerFSPs_[ixx]->id())==ihvy_&&ixx==(preshowerFSPs_.size()-(2+nph_+2)))||
             abs(preshowerFSPs_[ixx]->parents()[0]->id())==24) {
            preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
            getDescendents(preshowerFSPs_[ixx]);
            for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
              showeredFSPsToDelete_.push_back(tmpList_[jxx]);
          }
          unsigned int tmpUnsignedInt(4+nph_);
          if(ixx==preshowerFSPs_.size()-1&&preshowerFSPsToDelete_.size()!=tmpUnsignedInt) {
            throw Exception()
              << "FxFxHandler::caldel_m() - ERROR!\n"
              << "wphjet process should have " << 4+nph_ << " particles to omit from"
              << "jet-parton matching." << Exception::eventerror;
          }
        } else if(ihrd_==16) {
          /* 2Qph... <--- if Q is a top it will decay. */
          /* Exclude the hard photons and any children they
             may have produced from the jet parton matching
             as well as the heavy quarks and any children it
             may have produced. */
          // AND WHAT ABOUT THE PHOTON?
          if(ihvy_==6) {
            if(abs(preshowerFSPs_[ixx]->id())==22||
               abs(preshowerFSPs_[ixx]->parents()[0]->id())==6||
               abs(preshowerFSPs_[ixx]->parents()[0]->id())==24) {
              preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
              getDescendents(preshowerFSPs_[ixx]);
              for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
                showeredFSPsToDelete_.push_back(tmpList_[jxx]);
            }
            unsigned int tmpUnsignedInt(6+nph_);
            if(ixx==preshowerFSPs_.size()-1&&preshowerFSPsToDelete_.size()!=tmpUnsignedInt) {
              throw Exception()
                << "FxFxHandler::caldel_m() - ERROR!\n"
                << "wphjet process should have " << 6+nph_ << " particles to omit from"
                << "jet-parton matching for ihvy=" << ihvy_ << "." << Exception::eventerror;
            }
          } else {
            if(abs(preshowerFSPs_[ixx]->id())==22||
               (abs(preshowerFSPs_[ixx]->id())==ihvy_&&ixx<2)) {
              preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
              getDescendents(preshowerFSPs_[ixx]);
              for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
                showeredFSPsToDelete_.push_back(tmpList_[jxx]);
            }
            unsigned int tmpUnsignedInt(2+nph_);
            if(ixx==preshowerFSPs_.size()-1&&preshowerFSPsToDelete_.size()!=tmpUnsignedInt) {
              throw Exception()
                << "FxFxHandler::caldel_m() - ERROR!\n"
                << "wphjet process should have " << 2+nph_ << " particles to omit from"
                << "jet-parton matching for ihvy=" << ihvy_ << "." << Exception::eventerror;
            }
          }
        }
    }
  }
    //  cout << "partonsToMatch_.size()= " << partonsToMatch_.size() << " preshowerFSPsToDelete_.size() = " << preshowerFSPsToDelete_.size() << endl;
    // cout << "deleting preshowerFSPs" << endl;
    for(unsigned int ixx=0; ixx<preshowerFSPsToDelete_.size(); ixx++) {
      for(unsigned int jxx=0; jxx<partonsToMatch_.size(); jxx++) {
        if(preshowerFSPsToDelete_[ixx]==partonsToMatch_[jxx]) {
          //    cout << "deleting " << preshowerFSPsToDelete_[ixx]->id() << " with parent " << preshowerFSPsToDelete_[ixx]->parents()[0]->id() << endl; // " energy = " << preshowerFSPsToDelete_[ixx]->momentum().e()*GeV << endl;
          partonsToMatch_.erase(partonsToMatch_.begin()+jxx);
          break;
        }
      }
    }
    //cout << "partonsToMatch_.size() (AFTER) = " << partonsToMatch_.size() << endl;
    //cout << "deleting showeredFSPs" << endl;
    for(unsigned int ixx=0; ixx<showeredFSPsToDelete_.size(); ixx++) {
      for(unsigned int jxx=0; jxx<particlesToCluster_.size(); jxx++) {
        if(showeredFSPsToDelete_[ixx]==particlesToCluster_[jxx]) {
          //   cout << "deleting " << showeredFSPsToDelete_[ixx]->id() << " with parent " << showeredFSPsToDelete_[ixx]->parents()[0]->id() << endl; //" energy = " << preshowerFSPsToDelete_[ixx]->momentum().e()*GeV << endl;
          particlesToCluster_.erase(particlesToCluster_.begin()+jxx);
          break;
        }
      }
    }

    // Sanity check!
    if(partonsToMatch_.size()>particlesToCluster_.size()) {
      throw Exception()
        << "FxFxHandler::caldel_m() - ERROR!\n"
        << "No. of ME level partons to be matched to jets = "
        << partonsToMatch_.size() << "\n"
        << "No. of showered particles to build jets from  = "
        << particlesToCluster_.size() << "\n"
        << "There should be at least as many partons to\n"
        << "cluster as there are partons to match to.\n"
        << Exception::eventerror;
    }
    //  cout << "partonsToMatch_.size() (AFTER2) = " << partonsToMatch_.size() << endl;


  // Acid test.
  /* unsigned int tmpUnsignedInt(njets_);
  if(!inputIsNLO_&&partonsToMatch_.size()!=tmpUnsignedInt) {
    for(unsigned int ixx=0; ixx<partonsToMatch_.size(); ixx++) {
      if(abs(partonsToMatch_[ixx]->id())>=6&&
	 abs(partonsToMatch_[ixx]->id())!=21)
	throw Exception()
	  << "FxFxHandler::caldel_m() - ERROR!\n"
	  << "Found a parton to match to which is not a quark or gluon!"
	  << *partonsToMatch_[ixx] << "\n"
	  << Exception::eventerror;
    }
    throw Exception()
      << "FxFxHandler::caldel_m() - ERROR!\n"
      << "No. of ME level partons to be matched to jets = "
      << partonsToMatch_.size() << "\n"
      << "No. of light jets (njets) in AlpGen process = "
      << njets_ << "\n"
      << "These should be equal." << "\n"
      << Exception::eventerror;
      }    */
  //cout << "partonsToMatch_.size() (AFTER3) = " << partonsToMatch_.size() << endl;

  return;
}

// This looks for all descendents of a top up to but not including
// the W and b children.
void FxFxHandler::getTopRadiation(PPtr theParticle) const {
  ParticleVector theChildren(theParticle->children());
  for (unsigned int ixx=0; ixx<theChildren.size(); ixx++)
    if(theChildren[ixx]->children().size()==0)
      tmpList_.push_back(theChildren[ixx]);
    else if(abs(theChildren[ixx]->id())==5||abs(theChildren[ixx]->id())==24)
      return;
    else
      getTopRadiation(theChildren[ixx]);
  return;
}

void FxFxHandler::caldel_mg() const {
  
  preshowerFSPsToDelete_.clear();
  showeredFSPsToDelete_.clear();
  /* 
   * Get the MadGraph clustering information
   * from the Les Houches event tags
   */
  ptclust_.clear();
  getptclust();
  
  //for(unsigned int izz = 0; izz<ptclust_.size(); izz++) cout << "ptclust_[" << izz << "] = " << ptclust_[izz] << endl; 

  /* 
   * Get the COM energy of the colliding hadrons 
   *
   */
  getECOM();
  // cout << "ECOM_ = " << ECOM_ << endl;

  for(unsigned int ixx=0; ixx<partonsToMatch_.size(); ixx++) {
    if(ptclust_[ixx] == ECOM_){
      preshowerFSPsToDelete_.push_back(partonsToMatch_[ixx]);
      tmpList_.clear();
      getDescendents(partonsToMatch_[ixx]);
      for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
	showeredFSPsToDelete_.push_back(tmpList_[jxx]);
    }
  }

  for(unsigned int ixx=0; ixx<preshowerFSPsToDelete_.size(); ixx++) {
    for(unsigned int jxx=0; jxx<partonsToMatch_.size(); jxx++) {
      if(preshowerFSPsToDelete_[ixx]==partonsToMatch_[jxx]) {
        //cout << "deleting " << preshowerFSPsToDelete_[ixx]->id() << endl;
        partonsToMatch_.erase(partonsToMatch_.begin()+jxx);
        break;
      }
    }
  }

  for(unsigned int ixx=0; ixx<showeredFSPsToDelete_.size(); ixx++) {
    for(unsigned int jxx=0; jxx<particlesToCluster_.size(); jxx++) {
      if(showeredFSPsToDelete_[ixx]==particlesToCluster_[jxx]) {
        //   cout << "deleting " << showeredFSPsToDelete_[ixx]->id() << " with parent " << showeredFSPsToDelete_[ixx]->parents()[0]->id() << endl; //" energy = " << preshowerFSPsToDelete_[ixx]->momentum().e()*GeV << endl;
        particlesToCluster_.erase(particlesToCluster_.begin()+jxx);
        break;
      }
    }
  }

  // Sanity check!
  if(partonsToMatch_.size()>particlesToCluster_.size()) {
    throw Exception()
      << "FxFxHandler::caldel_m() - ERROR!\n"
      << "No. of ME level partons to be matched to jets = "
      << partonsToMatch_.size() << "\n"
      << "No. of showered particles to build jets from  = "
      << particlesToCluster_.size() << "\n"
      << "There should be at least as many partons to\n"
      << "cluster as there are partons to match to.\n"
      << Exception::eventerror;
  }



}


void FxFxHandler::caldel_hvq() const {

  // Fill partonsToMatch_ with only those pre-shower partons intended to
  // be used in heavy-quark-jet matching and fill particlesToCluster_ using
  // only those final state particles (post-shower) which are supposed
  // in the heavy-quark-jet clustering used to do merging. To begin with
  // these are made from the corresponding sets of particles that were
  // omitted from the initial jet-parton matching run.
  partonsToMatch_     = preshowerFSPsToDelete_;
  particlesToCluster_.resize(showeredFSPsToDelete_.size());
  for(unsigned int ixx=0; ixx<showeredFSPsToDelete_.size(); ixx++)
    particlesToCluster_[ixx] = showeredFSPsToDelete_[ixx];

  // Reset the arrays of particles to delete so that the partonsToMatch_
  // and particlesToCluster_ vectors can undergo further filtering.
  preshowerFSPsToDelete_.clear();
  showeredFSPsToDelete_.clear();

  // Determine further particles in partonsToMatch_ and particlesToCluster_
  // for deletion.
  for(unsigned int ixx=0; ixx<partonsToMatch_.size(); ixx++) {
    // If the progenitor particle is not a heavy 
    // quark we delete it and its descendents.
    if(abs(partonsToMatch_[ixx]->id())<4||abs(partonsToMatch_[ixx]->id())>6) {
      preshowerFSPsToDelete_.push_back(partonsToMatch_[ixx]);
      tmpList_.clear();
      getDescendents(partonsToMatch_[ixx]);
      for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
	showeredFSPsToDelete_.push_back(tmpList_[jxx]);
    // If the progenitor is a b quark from a top decay drop
    // it & it's descendents too!
    } else if(abs(partonsToMatch_[ixx]->id())==5&&
	      partonsToMatch_[ixx]->parents().size()>0&&
	      abs(partonsToMatch_[ixx]->parents()[0]->id())==6) {
	preshowerFSPsToDelete_.push_back(partonsToMatch_[ixx]);
	tmpList_.clear();
	getDescendents(partonsToMatch_[ixx]);
	for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
	  showeredFSPsToDelete_.push_back(tmpList_[jxx]);
    // If (it's a hvy quark not from a top decay and) it has a W/Z/H
    // as a parent [ditto].
    } else if(partonsToMatch_[ixx]->parents().size()>0&&
	      (abs(partonsToMatch_[ixx]->parents()[0]->id())==23||
	       abs(partonsToMatch_[ixx]->parents()[0]->id())==24||
	       abs(partonsToMatch_[ixx]->parents()[0]->id())==25)) {
	preshowerFSPsToDelete_.push_back(partonsToMatch_[ixx]);
	tmpList_.clear();
	getDescendents(partonsToMatch_[ixx]);
	for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
	  showeredFSPsToDelete_.push_back(tmpList_[jxx]);
    }
  }

  // Now do the necessary deleting from partonsToMatch_ and
  // particlesToCluster_.
  for(unsigned int ixx=0; ixx<preshowerFSPsToDelete_.size(); ixx++) {
    for(unsigned int jxx=0; jxx<partonsToMatch_.size(); jxx++) {
      if(preshowerFSPsToDelete_[ixx]==partonsToMatch_[jxx]) {
	partonsToMatch_.erase(partonsToMatch_.begin()+jxx);
	break;
      }
    }
  }
  for(unsigned int ixx=0; ixx<showeredFSPsToDelete_.size(); ixx++) {
    for(unsigned int jxx=0; jxx<particlesToCluster_.size(); jxx++) {
      if(showeredFSPsToDelete_[ixx]==particlesToCluster_[jxx]) {
	particlesToCluster_.erase(particlesToCluster_.begin()+jxx);
	break;
      }
    }
  }

  // Now we return to get the decaying top quarks and any
  // radiation they produced.
  ParticleVector intermediates(lastXCombPtr()->subProcess()->intermediates());
  for(unsigned int ixx=0; ixx<intermediates.size(); ixx++) {
    if(abs(intermediates[ixx]->id())==6) {
      partonsToMatch_.push_back(intermediates[ixx]);
      tmpList_.clear();
      getTopRadiation(partonsToMatch_.back());
      for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
	particlesToCluster_.push_back(tmpList_[jxx]);
    }
  }

  // If there are any heavy quark progenitors we have to remove
  // the final (showered) instance of them from particlesToCluster.
  ParticleVector evolvedHeavyQuarks;
  for(unsigned int ixx=0; ixx<partonsToMatch_.size(); ixx++) {
    if(abs(partonsToMatch_[ixx]->id())>=4&&abs(partonsToMatch_[ixx]->id())<=6) {
      theProgenitor = partonsToMatch_[ixx];
      // Follow the heavy quark line down to where it stops branching.
      while(theProgenitor->children().size()>0) {
	theLastProgenitor = theProgenitor;
	for(unsigned int jxx=0; jxx<theProgenitor->children().size(); jxx++) {
	  if(theProgenitor->children()[jxx]->id()==theProgenitor->id())
	    theProgenitor=theProgenitor->children()[jxx];
	}
	// If the progenitor had children but none of them had
	// the same particle id as it, then it must have undergone
	// a decay rather than a branching, i.e. it is the end of
	// the evolution line, so,
	if(theProgenitor==theLastProgenitor) break;
      }
      evolvedHeavyQuarks.push_back(theProgenitor);
    }
  }
  // Now delete the evolved heavy quark from the particlesToCluster.
  for(unsigned int ixx=0; ixx<evolvedHeavyQuarks.size(); ixx++) {
    for(unsigned int jxx=0; jxx<particlesToCluster_.size(); jxx++) {
      if(evolvedHeavyQuarks[ixx]==particlesToCluster_[jxx]) {
	particlesToCluster_.erase(particlesToCluster_.begin()+jxx);
	break;
      }
    }
  }
  
  return;
}

// get npLO_ and npNLO_ 
void FxFxHandler::getnpFxFx() const {
  split_vector_type SplitVec; 
  // pull the optional weights from the current event
  map<string,double> optionalEventWeights = eventHandler()->currentEvent()->optionalWeights();
  // loop over the optional weights and find np values
  for (map<string,double>::const_iterator it=optionalEventWeights.begin(); it!=optionalEventWeights.end(); ++it){
    // split the line 
    boost::split( SplitVec, it->first, boost::is_any_of(" ") );	
    // if np is found, store the information 
    if(SplitVec[0] == "np") {
      npLO_ = atof(SplitVec[1].c_str());
      npNLO_ = atof(SplitVec[2].c_str());
    } 
  }
  return;
}

// get hadron COM energy
void FxFxHandler::getECOM() const {
  split_vector_type SplitVec; 
  // pull the optional weights from the current event
  map<string,double> optionalEventWeights = eventHandler()->currentEvent()->optionalWeights();
  // loop over the optional weights and find np values
  for (map<string,double>::const_iterator it=optionalEventWeights.begin(); it!=optionalEventWeights.end(); ++it){
    // split the line 
    boost::split( SplitVec, it->first, boost::is_any_of(" ") );	
    // if np is found, store the information 
    if(SplitVec[0] == "ecom") {
      ECOM_ = it->second;
    } 
  }
  return;
}


// get pt_clust information 
void FxFxHandler::getptclust() const {
  split_vector_type SplitVec; 
  // pull the optional weights from the current event
  map<string,double> optionalEventWeights = eventHandler()->currentEvent()->optionalWeights();
  // loop over the optional weights and find np values
  string str_eq = "=";
  string str_quote = "\"";
  int countptc_(0);
  for (map<string,double>::const_iterator it=optionalEventWeights.begin(); it!=optionalEventWeights.end(); ++it){
    // split the line
    double wgtid = it->second;
    //cout << "wgtdid = " << wgtid << endl;
    if(wgtid == -333) {
      //   cout << "it->first for -333 = " << it->first << endl;
      boost::split( SplitVec, it->first, boost::is_any_of(" ") );	
      // if np is found, store the information
      double ptclust(0.);
      string stringtohandle("");
      for(unsigned int i_ = 0; i_ < SplitVec.size(); ++i_) {
          if(SplitVec[i_].find("pt_clust_") != std::string::npos) {
            //cout << "pt_clust_ found in " << SplitVec[i_] << endl;
            countptc_++;
            //  cout << "SplitVec[i_] = " << SplitVec[i_] << endl;
            stringtohandle = SplitVec[i_];
            stringtohandle.erase(0,10);
            erase_substr(stringtohandle,str_eq);
            erase_substr(stringtohandle,str_quote);
            //            cout << "stringtohandle = " << stringtohandle << endl;
            ptclust = atof(stringtohandle.c_str());
            ptclust_.push_back(ptclust);
          }
      }
    }
  }
  return;
}

void FxFxHandler::getPreshowerParticles() const {
  // LH file initial-state partons:
  preshowerISPs_ = lastXCombPtr()->subProcess()->incoming();
  // LH file final-state partICLEs:
  preshowerFSPs_ = lastXCombPtr()->subProcess()->outgoing();
  return;
}

void FxFxHandler::getShoweredParticles() const {
  // Post-shower initial-state hadrons:
  showeredISHs_ = eventHandler()->currentEvent()->incoming();

  // Post-shower initial-state partons:
  for(unsigned int ixx=0; ixx<(showeredISHs_.first)->children().size(); ixx++)
    if(((showeredISHs_.first)->children()[ixx]->id())<6||
       ((showeredISHs_.first)->children()[ixx]->id())==21)
      showeredISPs_.first=(showeredISHs_.first)->children()[ixx];
  for(unsigned int ixx=0; ixx<(showeredISHs_.second)->children().size(); ixx++)
    if(((showeredISHs_.second)->children()[ixx]->id())<6||
       ((showeredISHs_.second)->children()[ixx]->id())==21)
      showeredISPs_.second=(showeredISHs_.second)->children()[ixx];

  // Post-shower final-state partICLEs plus remnants (to be removed later):
  showeredFSPs_ = eventHandler()->currentEvent()->getFinalState();

  // Post-shower final-state remnants:
  for(unsigned int ixx=0; ixx<showeredFSPs_.size(); ixx++) {
    if(showeredFSPs_[ixx]->PDGName()=="Rem:p+"||
       showeredFSPs_[ixx]->PDGName()=="Rem:pbar-") {
      if(showeredFSPs_[ixx]->parents()[0]->parents()[0]==
	 showeredISHs_.first)
	showeredRems_.first=showeredFSPs_[ixx];
      else if(showeredFSPs_[ixx]->parents()[0]->parents()[0]==
	      showeredISHs_.second)
	showeredRems_.second=showeredFSPs_[ixx];
    }
  }

  // Now delete found remnants from the showeredFSPs vector for consistency.  
  for(unsigned int ixx=0; ixx<showeredFSPs_.size(); ixx++)
    if(showeredFSPs_[ixx]->PDGName()=="Rem:p+")
      showeredFSPs_.erase(showeredFSPs_.begin()+ixx);
  for(unsigned int ixx=0; ixx<showeredFSPs_.size(); ixx++)
    if(showeredFSPs_[ixx]->PDGName()=="Rem:pbar-")
      showeredFSPs_.erase(showeredFSPs_.begin()+ixx);
  sort(showeredFSPs_.begin(),showeredFSPs_.end(),recordEntry);

  return;
}

void FxFxHandler::doSanityChecks(int debugLevel) const {

  // When checking momentum conservation in the form 
  // p_in - p_out, any momentum component bigger / less
  // than + / - epsilon will result in the p_in - p_out
  // vector being flagged as "non-null", triggering a
  // warning that momentum conservation is violated.
  Energy epsilon(0.5*GeV);
  if(debugLevel>=5) epsilon=1e-9*GeV;

  // Print out what was found for the incoming and outgoing
  // partons in the lastXCombPtr regardless.
  if(debugLevel>=5) {
    cout << "\n\n\n\n";
    cout << "****************************************************" << endl;
    cout << " The following are the hard subprocess momenta from " << "\n"
	 << " lastXCombPtr and should be basically identical to  " << "\n"
	 << " the input LH file momenta." << "\n\n";
    cout << " Incoming particles:" << "\n"
	 << *(preshowerISPs_.first)      << "\n"
	 << *(preshowerISPs_.second)     << endl;
    cout << " Outgoing particles:" << endl;
    for(unsigned int ixx=0; ixx<preshowerFSPs_.size(); ixx++)
      cout << *(preshowerFSPs_[ixx]) << endl;
  }
  // Print out what was found for the incoming and outgoing
  // partons after the shower.
  if(debugLevel>=5) {
    cout << "\n\n";
    cout << "****************************************************" << endl;
    cout << " The following are the particles left at the end of" << "\n"
	 << " the showering step." << "\n\n";
    cout << " Incoming hadrons:"   << "\n"
	 << *(showeredISHs_.first)  << "\n"
	 << *(showeredISHs_.second) << endl;
    cout << " Incoming partons:"   << "\n"
	 << *(showeredISPs_.first)  << "\n"
	 << *(showeredISPs_.second) << endl;
    cout << " Outgoing partons:" << endl;
    for(unsigned int ixx=0; ixx<showeredFSPs_.size(); ixx++)
      cout << *(showeredFSPs_[ixx]) << endl;
    cout << " Outgoing remnants:"   << "\n"
	 << *(showeredRems_.first)  << "\n"
	 << *(showeredRems_.second) << endl;
  }

  // Check if we correctly identified all initial and final-state
  // particles by testing momentum is conserved.
  if(debugLevel>=4) {
    Lorentz5Momentum tmpMom;
    tmpMom += showeredISPs_.first->momentum();
    tmpMom += showeredISPs_.second->momentum();
    for(unsigned int ixx=0; ixx<showeredFSPs_.size(); ixx++)
      tmpMom -= showeredFSPs_[ixx]->momentum();
    if(!isMomLessThanEpsilon(tmpMom,epsilon))
      cout << "Total parton mom.in - total parton mom.out = "
	   <<  tmpMom/GeV << endl;
    tmpMom = showeredISHs_.first->momentum()
           - showeredRems_.first->momentum() -showeredISPs_.first->momentum();
    if(!isMomLessThanEpsilon(tmpMom,epsilon))
      cout << "First  p_hadron-p_remnant-p_incoming " << tmpMom/GeV << endl;
    tmpMom = showeredISHs_.second->momentum()
           - showeredRems_.second->momentum()-showeredISPs_.second->momentum();
    if(!isMomLessThanEpsilon(tmpMom,epsilon))
      cout << "Second p_hadron-p_remnant-p_incoming " << tmpMom/GeV << endl;
  }

  // Check if what we found to be the remnant is consistent with
  // what we identified as the parent incoming hadron i.e. p+
  // goes with Rem:p+ and pbar- goes with Rem:pbar-.
  if(debugLevel>=0) {
    string tmpString;
    tmpString=showeredRems_.first->PDGName();
    tmpString=tmpString.substr(tmpString.find_first_of(":")+1,
			       string::npos);
    if(showeredISHs_.first->PDGName()!=tmpString) {
      cout << "FxFxHandler::showerHardProcessVeto" << "\n"
	   << "Fatal error in pairing of remnant and parent hadron." << "\n"
	   << "Remnant = " << *(showeredRems_.first) << "\n"
	   << "Parent hadron = " << *(showeredISHs_.first)
	   << endl;
      cout << showeredISHs_.first->PDGName() << endl;
      cout << tmpString << endl;
    }
    tmpString=showeredRems_.second->PDGName();
    tmpString=tmpString.substr(tmpString.find_first_of(":")+1,
			       string::npos);
    if(showeredISHs_.second->PDGName()!=tmpString) {
      cout << "FxFxHandler::showerHardProcessVeto" << "\n"
	   << "Fatal error in pairing of remnant and parent hadron." << "\n"
	   << "Remnant = " << *(showeredRems_.second) << "\n"
	   << "Parent hadron = " << *(showeredISHs_.second)
	   << endl;
      cout << showeredISHs_.second->PDGName() << endl;
      cout << tmpString << endl;
    }
  }
  return;
}

void FxFxHandler::printMomVec(vector<Lorentz5Momentum> momVec) {
  cout << "\n\n";

  // Label columns.
  printf("%5s %9s %9s %9s %9s %9s %9s %9s %9s %9s\n",
	 "jet #",
	 "px","py","pz","E",
	 "eta","phi","pt","et","mass");

  // Print out the details for each jet
  for (unsigned int ixx=0; ixx<momVec.size(); ixx++) {
    printf("%5u %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n",
	   ixx,
	   double(momVec[ixx].x()/GeV),double(momVec[ixx].y()/GeV),
	   double(momVec[ixx].z()/GeV),double(momVec[ixx].t()/GeV),
	   double(momVec[ixx].eta())  ,double(momVec[ixx].phi()),
	   double(momVec[ixx].perp()/GeV),
	   double(momVec[ixx].et()/GeV),
	   double(momVec[ixx].m()/GeV));
  }

}

Energy FxFxHandler::etclusran_(double petc) const { 
  return (((2 * epsetclus_)/M_PI) * asin(2 * petc - 1) + etclusmean_);
}


void FxFxHandler::erase_substr(std::string& subject, const std::string& search) const {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
      subject.erase( pos, search.length() );
    }
}
