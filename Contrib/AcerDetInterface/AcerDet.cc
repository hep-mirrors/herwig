// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AcerDet class.
//

#include "AcerDet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Config/HepMCHelper.h"

using namespace Herwig;

// maximum number of leptons in the common block
const int acdleptt_maxlept = 12;
// maximum number of photons in the common block
const int acdphot_maxphot  = 12;
// maximum number of jets in the common block
const int acdjets_maxjet   = 20;

extern "C" {
  // main AcerDet driving routine
  void acerdet_(int*);
  // convert from HEPEVT to internal format
  void acdherwig6_();
  // open the files for input/output
  void acerdet_files_();
  // common block for reconstructed leptons
  extern struct {
    int nlept;
    int   kflept[acdleptt_maxlept];
    float pxlept[acdleptt_maxlept];
    float pylept[acdleptt_maxlept];
    float pzlept[acdleptt_maxlept];
    float eelept[acdleptt_maxlept];
  } acdleptt_;
  // common block for reconstructed photons
  extern struct {
    int nphot;
    int   kfphot[acdphot_maxphot];
    float pxphot[acdphot_maxphot];
    float pyphot[acdphot_maxphot];
    float pzphot[acdphot_maxphot];
    float eephot[acdphot_maxphot];
  } acdphot_;
  // common block for reconstructed jets
  extern struct {
    int njets;
    int   kfjets[acdjets_maxjet];
    float pxjets[acdjets_maxjet];
    float pyjets[acdjets_maxjet];
    float pzjets[acdjets_maxjet];
    float eejets[acdjets_maxjet];
  } acdjets_;
  // common block for reconstructed etmiss
  extern struct {
    float pxmiss;
    float pymiss;
    float pxnues;
    float pynues;
    float pxcalo;
    float pycalo;
  } acdmiss_;
}

void AcerDet::analyze(tEventPtr event, long ieve, int loop, int state) {
  // energy unit
  Energy eunit(GeV);
  if(HepMC::Units::default_momentum_unit()==HepMC::Units::MEV) eunit = GeV;
  // clear storage
  _nphoton=0;
  _photonMomentum.clear();
  _nlepton=0;
  _leptonMomentum.clear();
  _leptonID;
  _njet=0;
  _jetMomentum.clear();
  _jetID;
  // convert the event to HepMC
  HepMC::GenEvent * hepmc = HepMCConverter<HepMC::GenEvent>::convert(*event);
  // convert the event from HepMC to HEPEVT
  HepMC::HEPEVT_Wrapper::set_max_number_entries(4000);
  _converter->write_event(hepmc);
  // convert the event to AcerDet internal format
  acdherwig6_();
  int imode=0;
  acerdet_(&imode);
  // number of photons and momenta
  _nphoton = acdphot_.nphot;
  for(int ix=0;ix<_nphoton;++ix) {
    _photonMomentum.push_back(LorentzMomentum(acdphot_.pxphot[ix]*eunit,
					      acdphot_.pyphot[ix]*eunit,
					      acdphot_.pzphot[ix]*eunit,
					      acdphot_.eephot[ix]*eunit));
  }
  // number of leptons, momenta and PDG codes
  _nlepton = acdleptt_.nlept;
  for(int ix=0;ix<_nlepton;++ix) {
    _leptonMomentum.push_back(LorentzMomentum(acdleptt_.pxlept[ix]*eunit,
					      acdleptt_.pylept[ix]*eunit,
					      acdleptt_.pzlept[ix]*eunit,
					      acdleptt_.eelept[ix]*eunit));
    _leptonID.push_back(acdleptt_.kflept[ix]);
  }
  // number of jets, momenta and PDG codes
  _njet = acdjets_.njets;
  for(int ix=0;ix<_njet;++ix) {
    _jetMomentum.push_back(LorentzMomentum(acdjets_.pxjets[ix]*eunit,
					   acdjets_.pyjets[ix]*eunit,
					   acdjets_.pzjets[ix]*eunit,
					   acdjets_.eejets[ix]*eunit));
    _jetID.push_back(acdjets_.kfjets[ix]);
  }
  // missing ET
  _etcalo     = make_pair(acdmiss_.pxcalo*eunit,acdmiss_.pycalo*eunit);
  _etneutrino = make_pair(acdmiss_.pxnues*eunit,acdmiss_.pynues*eunit);
  _etstable   = make_pair(acdmiss_.pxmiss*eunit,acdmiss_.pymiss*eunit);
  delete hepmc;
}

NoPIOClassDescription<AcerDet> AcerDet::initAcerDet;
// Definition of the static class description member.

void AcerDet::Init() {

  static ClassDocumentation<AcerDet> documentation
    ("The AcerDet class provides an interface to the AcerDet fast"
     " dectector simulation.");

}

void AcerDet::dofinish() {
  AnalysisHandler::dofinish();
  // delete the converter from HepMC to HEPEVT
  delete _converter;
  int imode=1;
  acerdet_(&imode);
}

void AcerDet::doinitrun() {
  AnalysisHandler::doinitrun();
  // create the converter from HepMC to HEPEVT
  _converter = new HepMC::IO_HEPEVT();
  // initialise AcerDet
  int mode=-1;
  generator()->log() << "testing in init run \n";
  acerdet_files_();
  generator()->log() << "testing read files\n";
  acerdet_(&mode);
}
