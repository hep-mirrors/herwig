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
#include "Herwig/Shower/Base/PartnerFinder.h"
#include "Herwig/PDF/HwRemDecayer.h"
#include <queue>
#include "ThePEG/Utilities/Throw.h"
#include "Herwig/Shower/Base/KinematicsReconstructor.h"
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
    ycmax_(5.4),ycmin_(-5.4),jetAlgorithm_(1),vetoIsTurnedOff_(false),vetoSoftThanMatched_(false), etclusfixed_(true),epsetclus_(2.5*GeV)
  {}

void FxFxHandler::doinitrun() {
  ShowerHandler::doinitrun();
  // et_ holds the ET deposited in the (ncy_ x ncphi_) calorimeter cells.
  et_.resize(ncy_);
  for(unsigned int ixx=0; ixx<et_.size(); ixx++) et_[ixx].resize(ncphi_);
  // jetIdx_ for a given calorimeter cell this holds the index of the jet
  // that the cell was clustered into.
  jetIdx_.resize(ncy_);
  for(unsigned int ixx=0; ixx<jetIdx_.size(); ixx++) jetIdx_[ixx].resize(ncphi_);

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
      << cphcal_ << sphcal_ << cthcal_ << sthcal_ << ounit(epsetclus_,GeV);
}

void FxFxHandler::persistentInput(PersistentIStream & is, int) {
  is  >> alphaS_
      >> ncy_ >> ncphi_ >> ihvy_ >> nph_ >> nh_
      >> iunit(etclusmean_,GeV) >> rclus_ >> etaclmax_ >> rclusfactor_
      >> ihrd_ >> njets_ >> drjmin_ >> highestMultiplicity_
      >> ycmax_ >> ycmin_ >> jetAlgorithm_ >> vetoIsTurnedOff_ >> vetoSoftThanMatched_ >> etclusfixed_
      >> cphcal_ >> sphcal_ >> cthcal_ >> sthcal_ >> iunit(epsetclus_,GeV);
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
     "The AlpGen hard process code",
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
  ShowerHandler::dofinish();
}

void FxFxHandler::doinit() {

  //print error if HardProcID is not set in input file
  if(ihrd_ == -999) { cout << "Error: FxFxHandler:ihrd not set!" << endl; exit(1); }
  ShowerHandler::doinit();

  // Compute calorimeter edges in rapidity for GetJet algorithm.
  ycmax_=etaclmax_+rclus_;
  ycmin_=-ycmax_;

  // Initialise calorimeter.
  calini_m();
}

// Throws a veto according to MLM strategy ... when we finish writing it.
bool FxFxHandler::showerHardProcessVeto() {

  if(vetoIsTurnedOff_) return false;

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

  // Get npLO_ and npNLO_ for FxFx matching
  getnpFxFx();

  // print the npXLO_ values obtained 
  //  cout << "HANDLER:\t\t\t\t" << npLO_ << "\t\t" << npNLO_ << endl; 

  // Turn on some screen output debugging: 0 = none ---> 5 = very verbose.
  int debug_mode = 0;
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
  caldel_m();
  
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

  //end of FxFx part
  // **************************************************************** //
  // * Now look to the non-light partons for heavy quark processes. * //
  // **************************************************************** //

  if(ihrd_<=2||ihrd_==6||ihrd_==10||ihrd_==15||ihrd_==16) {

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
    
    //else {
      // If unmatched jets remain, reject the event if highestMultiplicity_!=1
      //if(!highestMultiplicity_) return true;
      //else {
      // If unmatched jets remain and highestMultiplicity is true then check
      // that these are softer than all the matched ones (from the light-parton
      // matching round).
    /*	Energy etmax(0.*GeV);
	for(unsigned int ixx=0; ixx<pjet_.size(); ixx++) etmax=max(etjet_[ixx],etmax);
	if(etmax>etmin) return true;
      }
    */

  }
  // Otherwise we accept the event ...
  return false;

}

/* Function that returns the R distance
   between a particle and a jet. */
double FxFxHandler::partonJetDeltaR(ThePEG::tPPtr partonptr, LorentzMomentum jetmom) { 
  LorentzMomentum partonmom(partonptr->momentum());
  // Calculate DY, DPhi and then DR
  double DY(partonmom.eta()-jetmom.eta());
  double DPhi(partonmom.phi()-jetmom.phi());
  if(DPhi>M_PI) DPhi=2*M_PI-DPhi;
  double DR(sqrt(sqr(DY)+sqr(DPhi)));
  return DR;
}

double FxFxHandler::partonJetDeltaR(LorentzMomentum jetmom1, LorentzMomentum jetmom2) { 
  // Calculate DY, DPhi and then DR
  double DY(jetmom1.eta()-jetmom2.eta());
  double DPhi(jetmom1.phi()-jetmom2.phi());
  if(DPhi>M_PI) DPhi=2*M_PI-DPhi;
  double DR(sqrt(sqr(DY)+sqr(DPhi)));
  return DR;
}


// Initialize calorimeter for calsim_m and getjet_m. Note that
// because initialization is separte calsim_m can be called more
// than once to simulate pileup of several events.
void FxFxHandler::calini_m() {

  // Making sure arrays are clear before filling;
  cphcal_.clear();  sphcal_.clear();
  cthcal_.clear();  sthcal_.clear();

  // Fill array holding phi values of calorimeter cell centres. 
  double deltaPhi(2*M_PI/ncphi_);
  for(unsigned int iphi=1; iphi<=ncphi_; iphi++) { 
    double phi(deltaPhi*(iphi-0.5)); // Goes phi~=0 to phi~=2*pi        (iphi=0--->ncphi).
    cphcal_.push_back(cos(phi));     // ==> goes from +1 ---> +1        (iphi=0--->ncphi).
    sphcal_.push_back(sin(phi));     // ==> goes 0 -> 1 -> 0 -> -1 -> 0 (iphi=0--->ncphi).
  }

  // Fill array holding theta values of calorimeter cell centres in Y. 
  double deltaY((ycmax_-ycmin_)/double(ncy_));
  for(unsigned int iy=1; iy<=ncy_; iy++) { 
    double Y(deltaY*(iy-0.5)+ycmin_);
    double th(2*atan(exp(-Y))); // Goes bwds th~=pi to fwds th~=0  (iy=0--->ncy).
    cthcal_.push_back(cos(th)); // ==> goes from -1 ---> +1        (iy=0--->ncy).
    sthcal_.push_back(sin(th)); // ==> goes from  0 ---> +1 ---> 0 (iy=0--->ncy).
  }
  return;
}

// Get FastJets
void FxFxHandler::getFastJets(double rjet, Energy ejcut, double etajcut) {

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
void FxFxHandler::getFastJetsToMatch(double rjet, Energy ejcut, double etajcut) {

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


// Simple calorimeter simulation - assume uniform Y and phi bins.
void FxFxHandler::calsim_m() {
  // Reset transverse energies of all calorimter cells ready for new fill.
  for(unsigned int ixx=0; ixx<et_.size(); ixx++)
    for(unsigned int iyy=0; iyy<et_[ixx].size(); iyy++)
      et_[ixx][iyy]=0*GeV;

  // Assign ET to each calorimeter cell (mostly 0's).
  for(unsigned int ipar=0; ipar<particlesToCluster_.size(); ipar++) { 
    double y(particlesToCluster_[ipar]->momentum().eta());
    if(y>=ycmin_&&y<=ycmax_) { 
      int absId(abs(particlesToCluster_[ipar]->id()));
      // If it's not a lepton / top / photon it goes in the calorimeter.
      if(!(absId>=11&&absId<=16) && absId!=6 && absId!=22) { 
	double phi(atan2(particlesToCluster_[ipar]->momentum().y()/GeV,
			 particlesToCluster_[ipar]->momentum().x()/GeV));
	if(phi<0) phi+=2*M_PI;
	unsigned int iy(int((y-ycmin_)/dely_));
	unsigned int iphi(int(phi/delphi_));
	et_[iy][iphi]+=particlesToCluster_[ipar]->momentum().e()*sthcal_[iy]; 
      }
    }
  }
  return;
}

// Find highest remaining cell > etstop and sum surrounding cells
// with -- delta(y)^2+delta(phi)^2 < Rjet^2 , ET>eccut. Keep sets
// with ET>ejcut and abs(eta)<etacut. 
void FxFxHandler::getjet_m(double rjet, Energy ejcut, double etajcut) {

  // Minimum ET the calorimeter can "see".
  Energy eccut(0.1*GeV);
  // So long as the cell remaining with the highest ET has
  // ET < etstop we try to cluster the surrounding cells into
  // it to potentially form a jet.
  Energy etstop(1.5*GeV);
  
  // Reset the vector holding the jet-index each calo cell
  // was clustered into.
  for(unsigned int iy=0; iy<ncy_; iy++)
    for(unsigned int iphi=0; iphi<ncphi_; iphi++)
      jetIdx_[iy][iphi]=-777;

  // Reset the vector that will hold the jet momenta.
  pjet_.clear();

  // # cells spanned by cone radius in phi dir., _rounded_down_.
  unsigned int nphi1(rjet/delphi_);
  // # cells spanned by cone radius in eta dir., _rounded_down_.
  unsigned int ny1(rjet/dely_);

  // Vector to hold the "ET" of each jet, where here ET really means
  // the scalar sum of ETs in each calo cell clustered into the jet.
  // Note that this is _not_ the same as the ET you would compute from
  // the final momentum worked out for each jet.
  etjet_.clear();

  // The ET of the highest ET cell found.
  Energy etmax(777e100*GeV);

  // Counter for number of highest ET calo cells found.
  unsigned int ipass(0);

  // Start finding jets.
  while(etmax>=etstop) {

    // Find the cell with the highest ET from
    // those not already assigned to a jet.
    etmax=0*GeV;
    int iymx(0), iphimx(0);
    for(unsigned int iphi=0; iphi<ncphi_; iphi++)
      for(unsigned int iy=0; iy<ncy_; iy++)
	if(et_[iy][iphi]>etmax&&jetIdx_[iy][iphi]<0) {
	  etmax  = et_[iy][iphi];
	  iymx   = iy;
	  iphimx = iphi;
       	}

    // If the remaining cell with the highest ET has ET < etstop, stop.
    if(etmax<etstop) break;

    // You cannot have more cells with the highest ET
    // so far than you have cells in the calorimeter.
    ipass++;
    if(ipass>(ncy_*ncphi_)) {
      cout << "FxFxHandler::getjet_m() - Fatal error." << endl;
      cout << "We found " << ipass << " calo cells with the highest ET so"
	   << "far\nbut the calorimeter only has " << ncy_*ncphi_ << " "
	   << "cells in it!" << endl;
      exit(10);
    }

    // Add a jet vector (may get deleted if jet fails ET / eta cuts).
    etjet_.push_back(0*GeV);
    pjet_.push_back(Lorentz5Momentum(0.*GeV,0.*GeV,0.*GeV,0.*GeV,0.*GeV));

    // Loop over all calo cells in range iphimx +/- nphi1 (inclusive)
    // wrapping round in azimuth if required.
    for(unsigned int iphi1=0; iphi1<=2*nphi1; iphi1++) {
      int iphix(iphimx-nphi1+iphi1);
      if(iphix<0)            iphix += ncphi_;
      if(iphix>=int(ncphi_)) iphix -= ncphi_;
      // Loop over all calo cells in range iymx +/- ny1 (inclusive).
      for(unsigned int iy1=0; iy1<=2*ny1; iy1++) {
	int iyx(iymx-ny1+iy1);
	// If the cell is outside the calorimeter OR if it was already
	// associated to a jet then skip to the next loop.
	if(iyx>=0&&iyx<int(ncy_)&&jetIdx_[iyx][iphix]<0) {
	  // N.B. iyx-iymx = iy1-ny1 and iphix-iphimx = iphi1-nphi1
	  //      hence the following is the distance in R between the
	  //      centre of the cell we are looking at and the one
	  //      with the highest ET.
	  double r2(sqr(  dely_*(double(iy1)  -double(ny1)  ))
		   +sqr(delphi_*(double(iphi1)-double(nphi1))));
	  if(r2<sqr(rjet)&&et_[iyx][iphix]>=eccut) {
	    Energy ECell(et_[iyx][iphix]/sthcal_[iyx]);
	    pjet_.back()+=LorentzMomentum(ECell*sthcal_[iyx]*cphcal_[iphix], // px
					  ECell*sthcal_[iyx]*sphcal_[iphix], // py
					  ECell*cthcal_[iyx],ECell);         // pz, E.
	    // N.B. This is the same reln as in ThePEG between phi and x,y.
	    etjet_.back()+=et_[iyx][iphix];
	    jetIdx_[iyx][iphix] = pjet_.size()-1; // Identify cell with this jet.
	  }
	}
      }
    }

    // Compute the current jet's mass.
    pjet_.back().rescaleMass();

    // Throw the jet away if it's ET is less than ejcut.
    if(etjet_.back()<ejcut||fabs(pjet_.back().eta())>etajcut) {
      pjet_.pop_back();
      etjet_.pop_back();
    }

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

// Deletes particles from partonsToMatch_ and particlesToCluster_
// vectors so that these contain only the partons to match to the
// jets and the particles used to build jets respectively. By and
// large the candidates for deletion are: vector bosons and their
// decay products, Higgs bosons, photons as well as _primary_, i.e.
// present in the lowest multiplicity process, heavy quarks and
// any related decay products.
void FxFxHandler::getDescendents(PPtr theParticle) {
  ParticleVector theChildren(theParticle->children());
  for (unsigned int ixx=0; ixx<theChildren.size(); ixx++)
    if(theChildren[ixx]->children().size()==0)
      tmpList_.push_back(theChildren[ixx]);
    else
      getDescendents(theChildren[ixx]);
  return;
}
void FxFxHandler::caldel_m() {


  preshowerFSPsToDelete_.clear();
  showeredFSPsToDelete_.clear();

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
	preshowerFSPsToDelete_.push_back(preshowerFSPs_[ixx]);
      }
      if(abs(preshowerFSPs_[ixx]->parents()[0]->id())==6) {
	getDescendents(preshowerFSPs_[ixx]->parents()[0]);
	for(unsigned int jxx=0; jxx<tmpList_.size(); jxx++)
	  showeredFSPsToDelete_.push_back(tmpList_[jxx]);
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
  // cout << "partonsToMatch_.size()= " << partonsToMatch_.size() << " preshowerFSPsToDelete_.size() = " << preshowerFSPsToDelete_.size() << endl;
  for(unsigned int ixx=0; ixx<preshowerFSPsToDelete_.size(); ixx++) {
    for(unsigned int jxx=0; jxx<partonsToMatch_.size(); jxx++) {
      if(preshowerFSPsToDelete_[ixx]==partonsToMatch_[jxx]) {
	partonsToMatch_.erase(partonsToMatch_.begin()+jxx);
	break;
      }
    }
  }
  //cout << "partonsToMatch_.size() (AFTER) = " << partonsToMatch_.size() << endl;

  for(unsigned int ixx=0; ixx<showeredFSPsToDelete_.size(); ixx++) {
    for(unsigned int jxx=0; jxx<particlesToCluster_.size(); jxx++) {
      if(showeredFSPsToDelete_[ixx]==particlesToCluster_[jxx]) {
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
void FxFxHandler::getTopRadiation(PPtr theParticle) {
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
void FxFxHandler::caldel_hvq() {

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
void FxFxHandler::getnpFxFx() {
  
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

void FxFxHandler::getPreshowerParticles() {
  // LH file initial-state partons:
  preshowerISPs_ = lastXCombPtr()->subProcess()->incoming();

  // LH file final-state partICLEs:
  preshowerFSPs_ = lastXCombPtr()->subProcess()->outgoing();

  return;
}

void FxFxHandler::getShoweredParticles() {
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

void FxFxHandler::doSanityChecks(int debugLevel) {

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

Energy FxFxHandler::etclusran_(double petc) { 
  return (((2 * epsetclus_)/M_PI) * asin(2 * petc - 1) + etclusmean_);
}
