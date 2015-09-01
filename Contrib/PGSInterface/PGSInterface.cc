// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PGSInterface class.
//

#include "PGSInterface.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Config/HepMCHelper.h"
#include "HepMC/HEPEVT_Wrapper.h"

using namespace Herwig;

extern "C" {
  // initialize PGS
  void pgs_initialize_();
  // call PGS trigger code
  void pgs_trigger_();
  // call PGS reconstruction code
  void pgs_recon_();
  // common block for PGS input parameters
  extern struct {
    int  numarg;              // number of arguments supplied to program
    char pgs_args[10][40];    // list of arguments (char*80)
    int  nevpgs;              // number of events to generate/read
    int  target_lum;          // target luminosity (in pb-1)
    int  nprpgs;              // number of events to print out 
    int  pgs_iseed,pgs_jseed; // seeds for pgs_ranmar
    int  pgs_log_unit;        // log file unit number
    char optpgs[6];           // type of run: 'PYTHIA', 'ISAJET', 'FILE', 
    char evtlum[6];           // number of events ('events') or luminosity ('pb-1')
    char pgs_input_file[80];  // input file
    char pgs_output_file[80]; // output file
    char pgs_log_file[80];    // log file
    char pgs_param_file[80];  // detector parameter file
    char pgs_isajet_decay[80];// ISAJET decay table file name
    char pgs_isajet_cards[80];// ISAJET card file name
    char pgs_pythia_cards[80];// PYTHIA card file name
    int  pgs_herwig_proc;     // HERWIG process to generate
    char pgs_herwig_susy[80]; // HERWIG SUSY data file
    char pgs_alpgen_stem[80]; // ALPGEN unweighted events file stem
  } pgsevt_;
  // common block for PGS reconstructed objects
  const int pgsrec_nojmx = 500;
  extern struct {
    int numobj,dumobj;               // number of reconstructed objects
    int indobj[pgsrec_nojmx];        // index to HEPEVT particle (where relevant)
    int typobj[pgsrec_nojmx];        // reconstructed type
    double pobj[pgsrec_nojmx][4];    // four vector of reconstructed object
    double qobj[pgsrec_nojmx];       // charge of reconstructed object
    double vecobj[pgsrec_nojmx][10]; // interesting object quantities
    int   unique[pgsrec_nojmx];     // true for object if it is uniquely identified
                                     // and passes cuts in pgs_object_cuts
  } pgsrec_;

  const int pgsnetamax=600;
  const int pgsnphimax=600;
  extern struct {
    double ecal[pgsnetamax][pgsnphimax];  // electromagnetic energy in each tower
    double hcal[pgsnetamax][pgsnphimax];  // hadronic energy in each tower
    double met_cal;                 // calorimeter missing ET
    double phi_met_cal;             // calorimeter missing ET phi
    double met_cor;                 // missing ET corrected for muons
    double phi_met_cor;             // corrected missing ET phi
  } pgscal_; 
  
}

void PGSInterface::persistentOutput(PersistentOStream & os) const {
  os << _pgs_param_file;
}

void PGSInterface::persistentInput(PersistentIStream & is, int) {
  is >> _pgs_param_file;
}

ClassDescription<PGSInterface> PGSInterface::initPGSInterface;
// Definition of the static class description member.

void PGSInterface::Init() {

  static ClassDocumentation<PGSInterface> documentation
    ("The PGSInterface class is designed to allow the PGS detector"
     " simulation to be used as an analysis handler in Herwig");

}

void PGSInterface::dofinish() {
  AnalysisHandler::dofinish();
  // delete the converter from HepMC to HEPEVT
  delete _converter;
}

void PGSInterface::doinitrun() {
  AnalysisHandler::doinitrun();
  // create the converter from HepMC to HEPEVT
  _converter = new HepMC::IO_HEPEVT();
  // convert variables and pass to PGS
  // name of detector parameters file
  const char *temp;
  temp = _pgs_param_file.c_str();
  for(unsigned int ix=0;ix<80;++ix) {
    if(temp[ix]=='\0') break;
    pgsevt_.pgs_param_file[ix]=temp[ix];
  }
  // initialize PGS
  pgs_initialize_(); 
}

void PGSInterface::analyze(tEventPtr event, long ieve, int loop, int state) {
  // check energy unit
  if(HepMC::Units::default_momentum_unit()!=HepMC::Units::GEV)
    throw Exception() << "Must be using GeV in HepMC if using PGS"
		      << " in PGSInterface::analyze()"
		      << Exception::runerror;
  // convert the event to HepMC
  HepMC::GenEvent * hepmc = HepMCConverter<HepMC::GenEvent>::convert(*event);
  HepMC::HEPEVT_Wrapper::set_max_number_entries(4000);
  // convert to HEPEVT
  _converter->write_event(hepmc);
  // call PGS trigger code
  pgs_trigger_();     
  // call PGS reconstruction code
  pgs_recon_();
  // convert the unique reconstructed objects
  _objects.clear();
  for(int ix=0;ix<pgsrec_.numobj;++ix) {
    if(!pgsrec_.unique[ix]) continue;
    _objects.push_back(ReconstructedObject());
    // set the type
    _objects.back().type = ObjectType(pgsrec_.typobj[ix]);
    // set the momentum
    _objects.back().momentum = 
      Lorentz5Momentum(pgsrec_.pobj[ix][0]*GeV, pgsrec_.pobj[ix][1]*GeV,
		       pgsrec_.pobj[ix][2]*GeV, pgsrec_.pobj[ix][3]*GeV);
    // by default mass is calculated from momentum
    _objects.back().momentum.rescaleMass();
    // set the charge
    _objects.back().charge = pgsrec_.qobj[ix];
    // the electromagentic energy
    _objects.back().emenergy  = pgsrec_.vecobj[ix][0]*GeV;
    // the hadronic energy
    _objects.back().hadenergy = pgsrec_.vecobj[ix][1]*GeV; 
    // energy of the tracks
    _objects.back().trackenergy = pgsrec_.vecobj[ix][2]*GeV;
    // the number of tracks
    _objects.back().numtracks = int(pgsrec_.vecobj[ix][3]);
    // now stuff which differs depending on the type of object
    // default is no B tagging
    _objects.back().btagging = None;
    // ET from momentum is the default
    _objects.back().ET = _objects.back().momentum.et();
    // isolation energy (zero by default)
    _objects.back().isolationET = Energy(); 
    // Transverse momentum in the isolation cone (zero by default)
    _objects.back().isolationpT = Energy();
    // Ratio of hadronic to electromagentic energy (zero by default)
    _objects.back().hadronicem = 0.;
    // Ratio of electromagnetic energy to track momentum (zero by default)
    _objects.back().ep = 0.;
    // track isolation energy (zero by default)
    _objects.back().trkisoEnergy = Energy(); 
    // Number of pi0 in cone for tau (zero by default)
    _objects.back().npi0 = 0;
    //  Sum of pt of pi0 not in cone for tau (zero by default)
    _objects.back().taupTpi0 = Energy();
    // Sum of pt of tracks not in cone for tau (zero by default)
    _objects.back().taupTtracks = Energy();
    // cluster width (zero by default)
    _objects.back().clusterWidth = Energy();
    // pT of highest track (zero by default)
    _objects.back().ptHightestTrack = Energy();
    // for photons
    if(_objects.back().type==Photon) {
      // PDG code
      _objects.back().PDGcode=ParticleID::gamma;
      // ET
      _objects.back().ET = pgsrec_.vecobj[ix][5]*GeV;
      // energy in the isolation cone
      _objects.back().isolationET = pgsrec_.vecobj[ix][6]*GeV;
      // pt in isolation cone
      _objects.back().isolationpT = pgsrec_.vecobj[ix][7]*GeV;
      // Ratio of hadronic to electromagentic energy
      _objects.back().hadronicem = pgsrec_.vecobj[ix][8];
      // Ratio of electromagnetic energy to track momentum 
      _objects.back().ep =  pgsrec_.vecobj[ix][9];
    }
    // for electrons
    else if(_objects.back().type==Electron) {
      // PDG code
      _objects.back().PDGcode = _objects.back().charge<0. ? 
	ParticleID::eminus : ParticleID::eplus;
      // ET
      _objects.back().ET = pgsrec_.vecobj[ix][5]*GeV;
      // energy in the isolation cone
      _objects.back().isolationET = pgsrec_.vecobj[ix][6]*GeV;
      // pt in isolation cone
      _objects.back().isolationpT = pgsrec_.vecobj[ix][7]*GeV;
      // Ratio of hadronic to electromagentic energy
      _objects.back().hadronicem = pgsrec_.vecobj[ix][8];
      // Ratio of electromagnetic energy to track momentum 
      _objects.back().ep =  pgsrec_.vecobj[ix][9];
    }
    // for muons
    else if(_objects.back().type==Muon) {
      // PDG code
      _objects.back().PDGcode = _objects.back().charge<0. ? 
	ParticleID::muminus : ParticleID::muplus;
      // pt in isolation cone
      _objects.back().isolationpT = pgsrec_.vecobj[ix][5]*GeV;
      // track isolation energy
      _objects.back().trkisoEnergy = pgsrec_.vecobj[ix][4]*GeV;
    }
    // for taus
    else if(_objects.back().type==Tau) {
      // PDG code
      _objects.back().PDGcode = _objects.back().charge<0. ? 
	ParticleID::tauminus : ParticleID::tauplus;
      // set the mass
      _objects.back().momentum.setMass(pgsrec_.vecobj[ix][5]*GeV);
      // Number of pi0 in cone for tau
      _objects.back().npi0 =int(pgsrec_.vecobj[ix][8]);
      //  Sum of pt of pi0 not in cone for tau
      _objects.back().taupTpi0 = pgsrec_.vecobj[ix][9]*GeV;
      // Sum of pt of tracks not in cone for tau
      _objects.back().taupTtracks = pgsrec_.vecobj[ix][7]*GeV;
      // cluster width
      _objects.back().clusterWidth = pgsrec_.vecobj[ix][4]*GeV;
      // pT of highest track
      _objects.back().ptHightestTrack = pgsrec_.vecobj[ix][6]*GeV;
    }
    // for jets
    else if(_objects.back().type==Jet) {
      // PDG code
      _objects.back().PDGcode = int(pgsrec_.vecobj[ix][5]);
      // b-tagging
      if(pgsrec_.vecobj[ix][6]>0.001&&pgsrec_.vecobj[ix][7]>0.001) 
	_objects.back().btagging = Both;
      else if(pgsrec_.vecobj[ix][6]>0.001)                         
	_objects.back().btagging = Loose;
      else if(pgsrec_.vecobj[ix][7]>0.001)                         
	_objects.back().btagging = Tight;
      else                                                         
	_objects.back().btagging = None;
      // cluster width
      _objects.back().clusterWidth = pgsrec_.vecobj[ix][4]*GeV;
    }
    // heavy charged
    else {
      _objects.back().PDGcode = 0;
    }
  }  
  // missing ET
  _calorimeterMET = make_pair(pgscal_.met_cal*GeV,pgscal_.phi_met_cal);
  _muonMET = make_pair(pgscal_.met_cor*GeV,pgscal_.phi_met_cor);
  // delete the HepMC event
  delete hepmc;
}
