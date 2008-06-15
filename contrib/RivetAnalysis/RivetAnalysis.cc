// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RivetAnalysis class.
//

#include "RivetAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Vectors/HepMCConverter.h"
#include "Herwig++/Config/HepMCHelper.h"
#include "HepMC/GenEvent.h"
//#include "Rivet/Rivet.hh"
#include "Rivet/AnalysisHandler.hh"
//#include "Rivet/RivetAIDA.hh"
//#include "AIDA/IManagedObject.h"

using namespace Herwig;

RivetAnalysis::RivetAnalysis() : _rivet() 
{}

void RivetAnalysis::analyze(ThePEG::tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // convert to hepmc
  HepMC::GenEvent * hepmc = ThePEG::HepMCConverter<HepMC::GenEvent>::convert(*event);
  // analyse the event
  _rivet->analyze(*hepmc);
  // delete hepmc event
  delete hepmc;
}

ThePEG::IBPtr RivetAnalysis::clone() const {
  return new_ptr(*this);
}

ThePEG::IBPtr RivetAnalysis::fullclone() const {
  return new_ptr(*this);
}

void RivetAnalysis::persistentOutput(ThePEG::PersistentOStream & os) const {
  os << _analyses;
}

void RivetAnalysis::persistentInput(ThePEG::PersistentIStream & is, int) {
  is >> _analyses;
}

ThePEG::ClassDescription<RivetAnalysis> RivetAnalysis::initRivetAnalysis;
// Definition of the static class description member.

void RivetAnalysis::Init() {

  static ThePEG::ClassDocumentation<RivetAnalysis> documentation
    ("The RivetAnalysis class is a simple class to allow analyses"
     " from the Rivet library to be called from Herwig++");

  static ThePEG::ParVector<RivetAnalysis,string> interfaceAnalyses
    ("Analyses",
     "The names of the Rivet analyses to use",
     &RivetAnalysis::_analyses, -1, "", "","" "",
     false, false, ThePEG::Interface::nolimits);

}

void RivetAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  _rivet->finalize();
  // this is needed but won't work because of the LWH crap is still in ThePEG
  //_rivet->tree().commit();
}

void RivetAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  if(_analyses.empty()) 
    throw ThePEG::Exception() << "Must have at least one analysis loaded in "
			      << "RivetAnalysis::doinitrun()"
			      << ThePEG::Exception::runerror;
  // create Rivet analysis handler
  _rivet = new Rivet::AnalysisHandler("Rivet");
  // specify the analyses to be used
  for(unsigned int ix=0;ix<_analyses.size();++ix) {
    _rivet->addAnalysis(_analyses[ix]);
  }
  // initialize the rivet analysis handler
  _rivet->init();
}
