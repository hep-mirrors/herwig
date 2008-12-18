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
#include "Herwig++/Config/HepMCHelper.h"
#include "HepMC/HEPEVT_Wrapper.h"

using namespace Herwig;

extern "C" {
  // initialize PGS
  void pgs_initialize__();
  // call PGS trigger code
  void pgs_trigger__();
  // call PGS reconstruction code
  void pgs_recon__();
  // common block for PGS input parameters
  extern struct {
    int  numarg;              // number of arguments supplied to program
    char pgs_args[25][80];    // list of arguments (char*80)
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
    bool   unique[pgsrec_nojmx];     // true for object if it is uniquely identified
                                     // and passes cuts in pgs_object_cuts
  } pgsrec_;
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
     " simulation to be used as an analysis handler in Herwig++");

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
  pgs_initialize__(); 
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
  pgs_trigger__();     
  // call PGS reconstruction code
  pgs_recon__();
  // convert the reconstructed objects
  _objects.clear();
  _objects.resize(pgsrec_.numobj);
  for(int ix=0;ix<pgsrec_.numobj;++ix) {
    // set the type
    _objects[ix].type =  ObjectType(pgsrec_.typobj[ix]);
    // set the momentum
    _objects[ix].momentum = Lorentz5Momentum(pgsrec_.pobj[ix][0]*GeV,
					     pgsrec_.pobj[ix][1]*GeV,
					     pgsrec_.pobj[ix][2]*GeV,
					     pgsrec_.pobj[ix][3]*GeV);
    // by default mass is calculated from momentum
    _objects[ix].momentum.rescaleMass();
    // set the charge
    _objects[ix].charge = pgsrec_.qobj[ix];
    // whether or not unique
    _objects[ix].unique = pgsrec_.unique[ix];
    // the electromagentic energy
    _objects[ix].emenergy  = pgsrec_.vecobj[ix][0]*GeV;
    // the hadronic energy
    _objects[ix].hadenergy = pgsrec_.vecobj[ix][1]*GeV; 
    // energy of the tracks
    _objects[ix].trackenergy = pgsrec_.vecobj[ix][2]*GeV;
    // the number of tracks
    _objects[ix].numtracks = int(pgsrec_.vecobj[ix][3]);
    // now stuff which differs depending on the type of object
    // default is no B tagging
    _objects[ix].btagging = None;
    // ET from momentum is the default
    _objects[ix].ET = _objects[ix].momentum.et();
    // isolation energy (zero by default)
    _objects[ix].isolationET = Energy(); 
    // Transverse momentum in the isolation cone (zero by default)
    _objects[ix].isolationpT = Energy();
    // Ratio of hadronic to electromagentic energy (zero by default)
    _objects[ix].hadronicem = 0.;
    // Ratio of electromagnetic energy to track momentum (zero by default)
    _objects[ix].ep = 0.;
    // track isolation energy (zero by default)
    _objects[ix].trkisoEnergy = Energy(); 
    // Number of pi0 in cone for tau (zero by default)
    _objects[ix].npi0 = 0;
    //  Sum of pt of pi0 not in cone for tau (zero by default)
    _objects[ix].taupTpi0 = Energy();
    // Sum of pt of tracks not in cone for tau (zero by default)
    _objects[ix].taupTtracks = Energy();
    // cluster width (zero by default)
    _objects[ix].clusterWidth = Energy();
    // pT of highest track (zero by default)
    _objects[ix].ptHightestTrack = Energy();
    // for photons
    if(_objects[ix].type==Photon) {
      // PDG code
      _objects[ix].PDGcode=ParticleID::gamma;
      // ET
      _objects[ix].ET = pgsrec_.vecobj[ix][5]*GeV;
      // energy in the isolation cone
      _objects[ix].isolationET = pgsrec_.vecobj[ix][6]*GeV;
      // pt in isolation cone
      _objects[ix].isolationpT = pgsrec_.vecobj[ix][7]*GeV;
      // Ratio of hadronic to electromagentic energy
      _objects[ix].hadronicem = pgsrec_.vecobj[ix][8];
      // Ratio of electromagnetic energy to track momentum 
      _objects[ix].ep =  pgsrec_.vecobj[ix][9];
    }
    // for electrons
    else if(_objects[ix].type==Electron) {
      // PDG code
      _objects[ix].PDGcode = _objects[ix].charge<0. ? ParticleID::eminus : ParticleID::eplus;
      // ET
      _objects[ix].ET = pgsrec_.vecobj[ix][5]*GeV;
      // energy in the isolation cone
      _objects[ix].isolationET = pgsrec_.vecobj[ix][6]*GeV;
      // pt in isolation cone
      _objects[ix].isolationpT = pgsrec_.vecobj[ix][7]*GeV;
      // Ratio of hadronic to electromagentic energy
      _objects[ix].hadronicem = pgsrec_.vecobj[ix][8];
      // Ratio of electromagnetic energy to track momentum 
      _objects[ix].ep =  pgsrec_.vecobj[ix][9];
    }
    // for muons
    else if(_objects[ix].type==Muon) {
      // PDG code
      _objects[ix].PDGcode = _objects[ix].charge<0. ? ParticleID::muminus : ParticleID::muplus;
      // pt in isolation cone
      _objects[ix].isolationpT = pgsrec_.vecobj[ix][5]*GeV;
      // track isolation energy
      _objects[ix].trkisoEnergy = pgsrec_.vecobj[ix][4]*GeV;
    }
    // for taus
    else if(_objects[ix].type==Tau) {
      // PDG code
      _objects[ix].PDGcode = _objects[ix].charge<0. ? ParticleID::tauminus : ParticleID::tauplus;
      // set the mass
      _objects[ix].momentum.setMass(pgsrec_.vecobj[ix][5]*GeV);
      // Number of pi0 in cone for tau
      _objects[ix].npi0 =int(pgsrec_.vecobj[ix][8]);
      //  Sum of pt of pi0 not in cone for tau
      _objects[ix].taupTpi0 = pgsrec_.vecobj[ix][9]*GeV;
      // Sum of pt of tracks not in cone for tau
      _objects[ix].taupTtracks = pgsrec_.vecobj[ix][7]*GeV;
      // cluster width
      _objects[ix].clusterWidth = pgsrec_.vecobj[ix][4]*GeV;
      // pT of highest track
      _objects[ix].ptHightestTrack = pgsrec_.vecobj[ix][6]*GeV;
    }
    // for jets
    else if(_objects[ix].type==Jet) {
      // PDG code
      _objects[ix].PDGcode = int(pgsrec_.vecobj[ix][5]);
      // b-tagging
      if(pgsrec_.vecobj[ix][6]>0.001&&pgsrec_.vecobj[ix][7]>0.001) _objects[ix].btagging = Both;
      else if(pgsrec_.vecobj[ix][6]>0.001)                         _objects[ix].btagging = Loose;
      else if(pgsrec_.vecobj[ix][7]>0.001)                         _objects[ix].btagging = Tight;
      else                                                         _objects[ix].btagging = None;
      // cluster width
      _objects[ix].clusterWidth = pgsrec_.vecobj[ix][4]*GeV;
    }
    // heavy charged
    else {
      _objects[ix].PDGcode = 0;
    }
  }  
  // delete the HepMC event
  delete hepmc;
}
