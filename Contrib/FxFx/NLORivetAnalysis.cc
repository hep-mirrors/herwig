// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLORivetAnalysis class.
//

#include "NLORivetAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Vectors/HepMCConverter.h"
#include "ThePEG/Config/HepMCHelper.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/SubProcessGroup.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "HepMC/GenEvent.h"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Logging.hh"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/IO_GenEvent.h"
#include <config.h>

using namespace ThePEG;

NLORivetAnalysis::NLORivetAnalysis() 
  :  _remnantId(82), _format(1),_unitchoice(),
   _geneventPrecision(16), debug(false), _rivet(), _nevent(0) {}

HepMC::GenEvent * NLORivetAnalysis::makeEvent(tEventPtr event, tSubProPtr sub, long no,
					  Energy eUnit, Length lUnit, 
					  CrossSection xsec, CrossSection xsecErr) const {
  
  // generate beam particles
  const PPair& beam = event->incoming();
  HepMC::GenParticle * b1 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.first->momentum(),beam.first->id(),
					      1,eUnit);
  HepMC::GenParticle * b2 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.second->momentum(),beam.second->id(),
					      1,eUnit);

  // generate remnants
  HepMC::GenParticle * r1 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.first->momentum() - 
					      sub->incoming().first->momentum(),
					      _remnantId,1,eUnit);
  HepMC::GenParticle * r2 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.second->momentum() - 
					      sub->incoming().second->momentum(),
					      _remnantId,1,eUnit);

  // generate outgoing particles
  vector<HepMC::GenParticle*> outgoing;
  for ( ParticleVector::const_iterator p = sub->outgoing().begin();
	p != sub->outgoing().end(); ++p ) {
    outgoing.push_back(HepMCTraits<HepMC::GenEvent>::newParticle((**p).momentum(),(**p).id(),
								 1,eUnit));
  }

  // generate one blob vertex
  HepMC::GenVertex * vertex = HepMCTraits<HepMC::GenEvent>::newVertex();

  HepMCTraits<HepMC::GenEvent>::addIncoming(*vertex,b1);
  HepMCTraits<HepMC::GenEvent>::addIncoming(*vertex,b2);

  HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,r1);
  HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,r2);

  for ( vector<HepMC::GenParticle*>::const_iterator p = outgoing.begin();
	p != outgoing.end(); ++p )
    HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,*p);

  HepMC::GenEvent * ev = 
    HepMCTraits<HepMC::GenEvent>::newEvent(no,event->weight()*sub->groupWeight(),
					   event->optionalWeights());

  //  cout << "event->weight()*sub->groupWeight()" << event->weight()*sub->groupWeight() << endl;
  /*print optional weights here
   for (map<string,double>::const_iterator it= event->optionalWeights().begin(); it!=event->optionalWeights().end(); ++it){
    std::cout << it->first << "  => " << it->second << '\n';
   }
   cout << endl;*/
  
  HepMCTraits<HepMC::GenEvent>::setUnits(*ev,eUnit,lUnit);
  HepMCTraits<HepMC::GenEvent>::setBeamParticles(*ev,b1,b2);

  HepMCTraits<HepMC::GenEvent>::addVertex(*ev,vertex);

  HepMCTraits<HepMC::GenEvent>::setCrossSection(*ev,xsec/picobarn,
						xsecErr/picobarn);

  return ev;

}

HepMC::GenEvent * NLORivetAnalysis::makeEventW(tEventPtr event, tSubProPtr sub, long no,
					  Energy eUnit, Length lUnit, 
					       CrossSection xsec, CrossSection xsecErr, double weighttest) const {
  
  // generate beam particles
  const PPair& beam = event->incoming();
  HepMC::GenParticle * b1 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.first->momentum(),beam.first->id(),
					      1,eUnit);
  HepMC::GenParticle * b2 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.second->momentum(),beam.second->id(),
					      1,eUnit);

  // generate remnants
  HepMC::GenParticle * r1 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.first->momentum() - 
					      sub->incoming().first->momentum(),
					      _remnantId,1,eUnit);
  HepMC::GenParticle * r2 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.second->momentum() - 
					      sub->incoming().second->momentum(),
					      _remnantId,1,eUnit);

  // generate outgoing particles
  vector<HepMC::GenParticle*> outgoing;
  for ( ParticleVector::const_iterator p = sub->outgoing().begin();
	p != sub->outgoing().end(); ++p ) {
    outgoing.push_back(HepMCTraits<HepMC::GenEvent>::newParticle((**p).momentum(),(**p).id(),
								 1,eUnit));
  }

  // generate one blob vertex
  HepMC::GenVertex * vertex = HepMCTraits<HepMC::GenEvent>::newVertex();

  HepMCTraits<HepMC::GenEvent>::addIncoming(*vertex,b1);
  HepMCTraits<HepMC::GenEvent>::addIncoming(*vertex,b2);

  HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,r1);
  HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,r2);

  for ( vector<HepMC::GenParticle*>::const_iterator p = outgoing.begin();
	p != outgoing.end(); ++p )
    HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,*p);

  HepMC::GenEvent * ev = 
    HepMCTraits<HepMC::GenEvent>::newEvent(no,event->weight()*sub->groupWeight()*weighttest,
					   event->optionalWeights());

// cout << "event->weight()*sub->groupWeight()" << event->weight()*sub->groupWeight() << endl;
  /*print optional weights here
   for (map<string,double>::const_iterator it= event->optionalWeights().begin(); it!=event->optionalWeights().end(); ++it){
    std::cout << it->first << "  => " << it->second << '\n';
   }
   cout << endl;*/
  
  HepMCTraits<HepMC::GenEvent>::setUnits(*ev,eUnit,lUnit);
  HepMCTraits<HepMC::GenEvent>::setBeamParticles(*ev,b1,b2);

  HepMCTraits<HepMC::GenEvent>::addVertex(*ev,vertex);

  HepMCTraits<HepMC::GenEvent>::setCrossSection(*ev,xsec/picobarn,
						xsecErr/picobarn);

  return ev;

}

void NLORivetAnalysis::analyze(ThePEG::tEventPtr event, long ieve, int loop, int state) {
  Energy eUnit;
  Length lUnit;
  switch (_unitchoice) {
  default: eUnit = GeV; lUnit = millimeter; break;
  case 1:  eUnit = MeV; lUnit = millimeter; break;
  case 2:  eUnit = GeV; lUnit = centimeter; break;
  case 3:  eUnit = MeV; lUnit = centimeter; break;
  }

  tcEHPtr eh = dynamic_ptr_cast<tcEHPtr>(event->primaryCollision()->handler());
  assert(eh);

  CrossSection xsec = eh->integratedXSec();
  CrossSection xsecErr = eh->integratedXSecErr();

  tSubProPtr sub = event->primarySubProcess();
  Ptr<SubProcessGroup>::tptr grp = 
    dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(sub);

  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // convert to hepmc


  HepMC::GenEvent * hepmc = 
    makeEvent(event,sub,_nevent,eUnit,lUnit,xsec,xsecErr);

  //count weights here
  std::vector<std::string> strs;
  if(_numweights == -999) {
    _numweights = 0;
    for (map<string,double>::const_iterator it= event->optionalWeights().begin(); it!=event->optionalWeights().end(); ++it){
      //std::cout << it->first << " => " << it->second << '\n';
      string first_piece = it->first;
      string word;
      istringstream iss(first_piece, istringstream::in);
      while( iss >> word ) strs.push_back(word);
      if(strs[0] == "PDF" || strs[0] == "SC") { _numweights++; } 
    }
    strs.clear();
  }

  //find the weights
  _i = 0;//counter
  for (map<string,double>::const_iterator it= event->optionalWeights().begin(); it!=event->optionalWeights().end(); ++it){
    //  std::cout << it->first << " => " << it->second << '\n';
    string first_piece = it->first;
    string word;
    istringstream iss(first_piece, istringstream::in);
    while( iss >> word ) {
      strs.push_back(word);
    }
    std::pair<int,double> OptWeightsTemp;
    if(strs[0] == "PDF") {
      OptWeightsTemp.first = atoi(strs[2].c_str());
      OptWeightsTemp.second= it->second;
      OptWeights.push_back(OptWeightsTemp);
      OptXS[_i] += it->second;
      _i++;
    }
    if(strs[0] == "SC") {
      OptWeightsTemp.first = atoi(strs[3].c_str());
      OptWeightsTemp.second = it->second;
      OptWeights.push_back(OptWeightsTemp);
      OptXS[_i] += it->second;
      _i++;
    }
    strs.clear();
  }
  
 
  /* multiple hepmcs for scale/pdf variations
   */
  vector<HepMC::GenEvent*>  hepmcMULTI;// = new HepMC::GenEvent[_numweights];
  HepMC::GenEvent * hepmcMULTIi;
  for(int rr = 0; rr < _numweights; rr++) {
    hepmcMULTIi = makeEventW(event,sub,_nevent,eUnit,lUnit,xsec,xsecErr,OptWeights[rr].second);
    hepmcMULTI.push_back(hepmcMULTIi); 
  }
  
  CurrentGenerator::Redirect stdout(cout);

  if ( _rivet ) _rivet->analyze(*hepmc);

  for(int rr = 0; rr < _numweights; rr++) {
    if ( _rivetMULTI[rr] ) _rivetMULTI[rr]->analyze(*hepmcMULTI[rr]);
  }
  // delete hepmc events
  delete hepmc;
  for(int rr = 0; rr < _numweights; rr++) {
    delete hepmcMULTI[rr];
  }
  
  if ( grp ) {
    std::cout << "grp is true" << endl;
    std::cout << endl;
    for ( SubProcessVector::const_iterator s = grp->dependent().begin();
	  s != grp->dependent().end(); ++s ) {

      hepmc = makeEvent(event,*s,_nevent,eUnit,lUnit,xsec,xsecErr);

      if ( _rivet ) _rivet->analyze(*hepmc);
      // delete hepmc event
      delete hepmc;
      

    }

  }

  ++_nevent;
}


ThePEG::IBPtr NLORivetAnalysis::clone() const {
  return new_ptr(*this);
}

ThePEG::IBPtr NLORivetAnalysis::fullclone() const {
  return new_ptr(*this);
}

void NLORivetAnalysis::persistentOutput(ThePEG::PersistentOStream & os) const {
  os << _analyses << filename << debug;
}

void NLORivetAnalysis::persistentInput(ThePEG::PersistentIStream & is, int) {
  is >> _analyses >> filename >> debug;
}

ThePEG::ClassDescription<NLORivetAnalysis> NLORivetAnalysis::initNLORivetAnalysis;
// Definition of the static class description member.

void NLORivetAnalysis::Init() {

  static ThePEG::ClassDocumentation<NLORivetAnalysis> documentation
    ("The NLORivetAnalysis class is a simple class to allow analyses"
     " from the Rivet library to be called from ThePEG");

  static ThePEG::ParVector<NLORivetAnalysis,string> interfaceAnalyses
    ("Analyses",
     "The names of the Rivet analyses to use",
     &NLORivetAnalysis::_analyses, -1, "", "","" "",
     false, false, ThePEG::Interface::nolimits);

  static Parameter<NLORivetAnalysis,long> interfaceRemnantId
    ("RemnantId",
     "Set the PDG id to be used for remnants.",
     &NLORivetAnalysis::_remnantId, 82, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<NLORivetAnalysis,string> interfaceFilename
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
     &NLORivetAnalysis::filename, "", true, false);


  static Switch<NLORivetAnalysis,bool> interfaceDebug
    ("Debug",
     "Enable debug information from Rivet",
     &NLORivetAnalysis::debug, false, true, false);
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

void NLORivetAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  if( _nevent > 0 && _rivet ) {
    CurrentGenerator::Redirect stdout(cout);
    _rivet->setCrossSection(generator()->integratedXSec()/picobarn);
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
  cout << "finalizing with numweights = " << _numweights << endl;

  for(int rr = 0; rr < _numweights; rr++) {
    OptXS[rr] = OptXS[rr]/_nevent;
  }

  for(int rr = 0; rr < _numweights; rr++) {
    cout << "cross section = " << OptXS[rr] << endl;
    if( _nevent > 0 && _rivetMULTI[rr] ) {
      _rivetMULTI[rr]->setCrossSection(OptXS[rr]);
      _rivetMULTI[rr]->finalize();
      
      string fname = filename;
#if ThePEG_RIVET_VERSION == 1
      if ( fname.empty() ) fname = generator()->runName() + "_" + boost::lexical_cast<string>(OptWeights[rr].first) + ".aida";
#elif ThePEG_RIVET_VERSION > 1
      if ( fname.empty() ) fname = generator()->runName() + "_" + boost::lexical_cast<string>(OptWeights[rr].first) + ".yoda";
#else
#error "Unknown ThePEG_RIVET_VERSION"
#endif
      _rivetMULTI[rr]->writeData(fname);
    }
    delete _rivetMULTI[rr];
    _rivetMULTI[rr] = 0;
  }
}

void NLORivetAnalysis::doinit() {
  _numweights = -999;
  
  for(int rr = 0; rr < 120; rr++) {
    OptXS.push_back(0.);
  }
  
  AnalysisHandler::doinit();
  if(_analyses.empty()) 
    throw ThePEG::Exception() << "Must have at least one analysis loaded in "
			      << "NLORivetAnalysis::doinitrun()"
			      << ThePEG::Exception::runerror;

  // check that analysis list is available
  _rivet = new Rivet::AnalysisHandler; //(fname);
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

void NLORivetAnalysis::doinitrun() {
  _numweights = -999;
  AnalysisHandler::doinitrun();
  // create NLORivet analysis handler
  CurrentGenerator::Redirect stdout(cout);
  _rivet = new Rivet::AnalysisHandler; //(fname);
  _rivet->addAnalyses(_analyses);

  for(int rr = 0; rr < 120; rr++) {
    OptXS.push_back(0.);
    _rivetMULTI[rr] = new Rivet::AnalysisHandler; //(fname);
    _rivetMULTI[rr]->addAnalyses(_analyses);
  }
  
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
