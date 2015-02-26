// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RivetAnalysis class.
//

#include "RivetAnalysis.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/SamplerBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Vectors/HepMCConverter.h"
#include "ThePEG/Config/HepMCHelper.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "HepMC/GenEvent.h"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Tools/Logging.hh"
#include <config.h>

using namespace ThePEG;

RivetAnalysis::RivetAnalysis() :  debug(false), _rivet(), _nevent(0) 
{}

void RivetAnalysis::analyze(ThePEG::tEventPtr event, long ieve, int loop, int state) {
  ++_nevent;
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // convert to hepmc
  HepMC::GenEvent * hepmc = ThePEG::HepMCConverter<HepMC::GenEvent>::convert(*event);
  // analyse the event
  CurrentGenerator::Redirect stdout(cout);
  if ( _rivet ){
    _rivet->analyze(*hepmc);
    _rivet->finishEvent();
  }
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
  os << _analyses << _paths << filename << debug;
}

void RivetAnalysis::persistentInput(ThePEG::PersistentIStream & is, int) {
  is >> _analyses >> _paths >> filename >> debug;
}

ThePEG::ClassDescription<RivetAnalysis> RivetAnalysis::initRivetAnalysis;
// Definition of the static class description member.

void RivetAnalysis::Init() {

  static ThePEG::ClassDocumentation<RivetAnalysis> documentation
    ("The RivetAnalysis class is a simple class to allow analyses"
     " from the Rivet library to be called from ThePEG");

  static ThePEG::ParVector<RivetAnalysis,string> interfaceAnalyses
    ("Analyses",
     "The names of the Rivet analyses to use",
     &RivetAnalysis::_analyses, -1, "", "","" "",
     false, false, ThePEG::Interface::nolimits);

  static ThePEG::ParVector<RivetAnalysis,string> interfacePaths
    ("Paths",
     "The directory paths where Rivet should look for analyses.",
     &RivetAnalysis::_paths, -1, "", "","" "",
     false, false, ThePEG::Interface::nolimits);

  static Parameter<RivetAnalysis,string> interfaceFilename
    ("Filename",
#if ThePEG_RIVET_VERSION == 1
     "The name of the file where the AIDA histograms are put. If empty, "
     "the run name will be used instead. '.aida' will in any case be "
     "appended to the file name.",
#elif ThePEG_RIVET_VERSION > 1
     "The name of the file where the YODA histograms are put. If empty, "
     "the run name will be used instead. '.yoda' will in any case be "
     "appended to the file name.",
#else
#error "Unknown ThePEG_RIVET_VERSION"
#endif
     &RivetAnalysis::filename, "", true, false);


  static Switch<RivetAnalysis,bool> interfaceDebug
    ("Debug",
     "Enable debug information from Rivet",
     &RivetAnalysis::debug, false, true, false);
  static SwitchOption interfaceDebugNo
    (interfaceDebug,
     "No",
     "Disable debug information.",
     false);
  static SwitchOption interfaceDebugYes
    (interfaceDebug,
     "Yes",
     "Enable debug information from Rivet.",
     true);


  interfaceAnalyses.rank(10);

}

void RivetAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  if( _nevent > 0 && _rivet ) {
    CurrentGenerator::Redirect stdout(cout);
    _rivet->setCrossSection(generator()->integratedXSec()/picobarn);
    StdEHPtr tmp = dynamic_ptr_cast<StdEHPtr>(generator()->currentEventHandler());
    _rivet->setAttempts(long(tmp->sampler()->attempts()));
    _rivet->finalize();

    string fname = filename;
#if ThePEG_RIVET_VERSION == 1
    if ( fname.empty() ) fname = generator()->runName() + ".aida";
#elif ThePEG_RIVET_VERSION > 1
    if ( fname.empty() ) fname = generator()->runName() + ".yoda";
#else
#error "Unknown ThePEG_RIVET_VERSION"
#endif
    _rivet->writeData(fname);
  }
  delete _rivet;
  _rivet = 0;
}

void RivetAnalysis::doinit() {
  AnalysisHandler::doinit();
  if(_analyses.empty()) 
    throw ThePEG::Exception() << "Must have at least one analysis loaded in "
			      << "RivetAnalysis::doinitrun()"
			      << ThePEG::Exception::runerror;

  // check that analysis list is available
  _rivet = new Rivet::AnalysisHandler; //(fname);
  for ( int i = 0, N = _paths.size(); i < N; ++i ) Rivet::addAnalysisLibPath(_paths[i]);
  _rivet->addAnalyses(_analyses);
  if ( _rivet->analysisNames().size() != _analyses.size() ) {
    throw ThePEG::Exception() 
      << "Rivet could not find all requested analyses.\n"
      << "Use 'rivet --list-analyses' to check availability.\n"
      << ThePEG::Exception::runerror;
  }
  delete _rivet;
  _rivet = 0;
}

void RivetAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  // create Rivet analysis handler
  CurrentGenerator::Redirect stdout(cout);
  _rivet = new Rivet::AnalysisHandler; //(fname);
  for ( int i = 0, N = _paths.size(); i < N; ++i ) Rivet::addAnalysisLibPath(_paths[i]);
  _rivet->addAnalyses(_analyses);
  // check that analysis list is still available
  if ( _rivet->analysisNames().size() != _analyses.size() ) {
    throw ThePEG::Exception() 
      << "Rivet could not find all requested analyses.\n"
      << "Use 'rivet --list-analyses' to check availability.\n"
      << ThePEG::Exception::runerror;
  }
  if ( debug )
    Rivet::Log::setLevel("Rivet",Rivet::Log::DEBUG);
}
