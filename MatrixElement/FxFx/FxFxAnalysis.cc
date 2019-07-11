// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FxFxAnalysis class.
//

#include "FxFxAnalysis.h"
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
#include "ThePEG/Utilities/XSecStat.h"

using namespace ThePEG;

FxFxAnalysis::FxFxAnalysis() 
  :  _remnantId(82), _format(1),_unitchoice(),
     _geneventPrecision(16), debug(false), _rivet(), _nevent(0), useoptweights(false), normoptweights(false) {}

HepMC::GenEvent * FxFxAnalysis::makeEvent(tEventPtr event, tSubProPtr sub, long no,
					  Energy eUnit, Length lUnit, 
					  CrossSection xsec, CrossSection xsecErr) const {
  //convert the event from the Herwig format to the HepMC format and write it to the common block
  HepMC::GenEvent * ev = HepMCConverter<HepMC::GenEvent>::convert(*event, false,eUnit, lUnit);
 
  //reset the event 
  HepMCTraits<HepMC::GenEvent>::resetEvent(ev, no, event->weight()*sub->groupWeight(), event->optionalWeights());

  //set the cross section
  HepMCTraits<HepMC::GenEvent>::setCrossSection(*ev,xsec/picobarn,
						xsecErr/picobarn);
  return ev;
}

HepMC::GenEvent * FxFxAnalysis::makeEventW(tEventPtr event, tSubProPtr sub, long no,
					  Energy eUnit, Length lUnit, 
                                           CrossSection xsec, CrossSection xsecErr, double evoptweight, double centralweight) const {

  //convert the event from the Herwig format to the HepMC format and write it to the common block
  HepMC::GenEvent * ev = HepMCConverter<HepMC::GenEvent>::convert(*event, false,eUnit, lUnit);

  if(normoptweights) { evoptweight /= centralweight; }
  
  //reset the event 
  HepMCTraits<HepMC::GenEvent>::resetEvent(ev, no, evoptweight, event->optionalWeights());
  //set the cross section
  HepMCTraits<HepMC::GenEvent>::setCrossSection(*ev,xsec/picobarn,
						xsecErr/picobarn);
  return ev;
}

void FxFxAnalysis::analyze(ThePEG::tEventPtr event, long ieve, int loop, int state) {
  Energy eUnit;
  Length lUnit;
  switch (_unitchoice) {
  default: eUnit = GeV; lUnit = millimeter; break;
  case 1:  eUnit = MeV; lUnit = millimeter; break;
  case 2:  eUnit = GeV; lUnit = centimeter; break;
  case 3:  eUnit = MeV; lUnit = centimeter; break;
  }

  tcFxFxEventHandlerPtr eh = dynamic_ptr_cast<tcFxFxEventHandlerPtr>(event->primaryCollision()->handler());
  assert(eh);

  CrossSection xsec = eh->integratedXSec();
  CrossSection xsecErr = eh->integratedXSecErr();
  
  optxsec = eh->optintegratedXSecMap();

  
  int ii = 0;
  if(useoptweights) { 
    for (map<string,CrossSection>::const_iterator it= optxsec.begin(); it!=optxsec.end(); ++it){
      OptXS[ii] =  it->second/picobarn;
      ii++;
    }
  }

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
      string first_piece = it->first;
      string word;
      istringstream iss(first_piece, istringstream::in);
      while( iss >> word ) strs.push_back(word);
      if(strs[0] != "np") { _numweights++; }
      strs.clear();
    }
  }
  double CentralWeight = 1.;
  if(useoptweights) { 
    //find the weights
    _i = 0;//counter
    OptWeights.clear();
    for (map<string,double>::const_iterator it= event->optionalWeights().begin(); it!=event->optionalWeights().end(); ++it){
      // std::cout << it->first << " => " << it->second << '\n';
      string first_piece = it->first;
      string word;
      istringstream iss(first_piece, istringstream::in);
      while( iss >> word ) {
	strs.push_back(word);
      }
      std::pair<int,double> OptWeightsTemp;
      /*      if(strs[0] == "PDF") {
	OptWeightsTemp.first = atoi(strs[2].c_str());
	OptWeightsTemp.second= it->second;
	OptWeights.push_back(OptWeightsTemp);
	_i++;
      }
      if(strs[0] == "SC") {
	OptWeightsTemp.first = atoi(strs[3].c_str());
	OptWeightsTemp.second = it->second;
        //        cout << "OptWeightsTemp.first = " << OptWeightsTemp.first << " OptWeightsTemp.second = " << OptWeightsTemp.second << endl;
        if(OptWeightsTemp.first == 1001) { cout << "OptWeightsTemp.second = " << OptWeightsTemp.second << endl; CentralWeight = OptWeightsTemp.second; } 
	OptWeights.push_back(OptWeightsTemp);
	_i++;
      }*/
      if(strs[0] != "np") {
        OptWeightsTemp.first = atoi((it->first).c_str());
        OptWeightsTemp.second= it->second;
        if(OptWeightsTemp.first == 1001) { /*cout << "OptWeightsTemp.second = " << OptWeightsTemp.second << endl;*/ CentralWeight = OptWeightsTemp.second; } 
        OptWeights.push_back(OptWeightsTemp);
        _i++;
      }
      strs.clear();
    }
  }
  
  /* multiple hepmcs for scale/pdf variations
   */
  vector<HepMC::GenEvent*>  hepmcMULTI;// = new HepMC::GenEvent[_numweights];
  HepMC::GenEvent * hepmcMULTIi;
  if(useoptweights) { 
    for(int rr = 0; rr < _numweights; rr++) {
      double xsrr = optxsec[std::to_string(OptWeights[rr].first)]/picobarn;
      /* cout << "xsec = " << xsec/picobarn << endl;
      cout << "OptWeights[rr].second = " << OptWeights[rr].second << endl;
      cout << "xsrr = " << xsrr << endl;*/
      hepmcMULTIi = makeEventW(event,sub,_nevent,eUnit,lUnit,xsrr*picobarn,xsecErr,OptWeights[rr].second, CentralWeight);
      hepmcMULTI.push_back(hepmcMULTIi); 
    }
  }
  
  CurrentGenerator::Redirect stdout(cout);

  if ( _rivet ) _rivet->analyze(*hepmc);

  if(useoptweights) { 
    for(int rr = 0; rr < _numweights; rr++) {
      if ( _rivetMULTI[rr] ) _rivetMULTI[rr]->analyze(*hepmcMULTI[rr]);
    }
  }
  // delete hepmc events
  delete hepmc;
  if(useoptweights) {  
    for(int rr = 0; rr < _numweights; rr++) {
      delete hepmcMULTI[rr];
    }
  }

  if ( grp ) {
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


ThePEG::IBPtr FxFxAnalysis::clone() const {
  return new_ptr(*this);
}

ThePEG::IBPtr FxFxAnalysis::fullclone() const {
  return new_ptr(*this);
}

void FxFxAnalysis::persistentOutput(ThePEG::PersistentOStream & os) const {
  os << _analyses << filename << debug << useoptweights << normoptweights;
}

void FxFxAnalysis::persistentInput(ThePEG::PersistentIStream & is, int) {
  is >> _analyses >> filename >> debug >> useoptweights >> normoptweights;
}

ThePEG::ClassDescription<FxFxAnalysis> FxFxAnalysis::initFxFxAnalysis;
// Definition of the static class description member.

void FxFxAnalysis::Init() {
  static ThePEG::ClassDocumentation<FxFxAnalysis> documentation
    ("The FxFxAnalysis class is a simple class to allow analyses"
     " from the Rivet library to be called from ThePEG");

  static ThePEG::ParVector<FxFxAnalysis,string> interfaceAnalyses
    ("Analyses",
     "The names of the Rivet analyses to use",
     &FxFxAnalysis::_analyses, -1, "", "","" "",
     false, false, ThePEG::Interface::nolimits);

  static Parameter<FxFxAnalysis,long> interfaceRemnantId
    ("RemnantId",
     "Set the PDG id to be used for remnants.",
     &FxFxAnalysis::_remnantId, 82, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<FxFxAnalysis,string> interfaceFilename
    ("Filename",
     "The name of the file where the YODA histograms are put. If empty, "
     "the run name will be used instead. '.yoda' will in any case be "
     "appended to the file name.",
     &FxFxAnalysis::filename, "", true, false);


  static Switch<FxFxAnalysis,bool> interfaceDebug
    ("Debug",
     "Enable debug information from Rivet",
     &FxFxAnalysis::debug, false, true, false);
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

  
  static Switch<FxFxAnalysis,bool> interfaceUseOptWeights
    ("UseOptWeights",
     "Enable debug information from Rivet",
     &FxFxAnalysis::useoptweights, false, true, false);
  static SwitchOption interfaceUseOptWeightsNo
    (interfaceUseOptWeights,
     "No",
     "Disable optional weights",
     false);
  static SwitchOption interfaceUseOptWeightsYes
    (interfaceUseOptWeights,
     "Yes",
     "Enable debug information from Rivet.",
     true);

    static Switch<FxFxAnalysis,bool> interfaceNormOptWeights
    ("NormOptWeights",
     "Enable debug information from Rivet",
     &FxFxAnalysis::normoptweights, false, true, false);
  static SwitchOption interfaceNormOptWeightsNo
    (interfaceNormOptWeights,
     "No",
     "Do not normalize optional weights",
     false);
  static SwitchOption interfacedNormOptWeightsYes
    (interfaceNormOptWeights,
     "Yes",
     "Normalize optional weights.",
     true);




  interfaceAnalyses.rank(10);

}

void FxFxAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  if( _nevent > 0 && _rivet ) {
    CurrentGenerator::Redirect stdout(cout);
    _rivet->setCrossSection(generator()->integratedXSec()/picobarn);
    _rivet->finalize();

    string fname = filename;
    fname = generator()->runName() + ".yoda";
    _rivet->writeData(fname);
  }
  
  delete _rivet;
  _rivet = 0;

  if(useoptweights) {
    for(int rr = 0; rr < _numweights; rr++) {
      cout << (OptWeights[rr].first)  << ", cross section = " << OptXS[rr] << endl;
      if( _nevent > 0 && _rivetMULTI[rr] ) {
	double xsrr = optxsec[std::to_string(OptWeights[rr].first)]/picobarn;
	//      _rivetMULTI[rr]->setCrossSection(OptXS[rr]);
	_rivetMULTI[rr]->setCrossSection(xsrr);
	_rivetMULTI[rr]->finalize();
	
	string fname = filename;
	fname = generator()->runName() + "_" + std::to_string(OptWeights[rr].first) + ".yoda";
	
	_rivetMULTI[rr]->writeData(fname);
      }
      delete _rivetMULTI[rr];
      _rivetMULTI[rr] = 0;
    }
  }
}

void FxFxAnalysis::doinit() {
  _numweights = -999;
  
  for(int rr = 0; rr < 120; rr++) {
    OptXS.push_back(0.);
  }
  AnalysisHandler::doinit();
  if(_analyses.empty()) 
    throw ThePEG::Exception() << "Must have at least one analysis loaded in "
			      << "FxFxAnalysis::doinitrun()"
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

void FxFxAnalysis::doinitrun() {
  _numweights = -999;
  AnalysisHandler::doinitrun();
  // create FxFx analysis handler
  if(useoptweights) { cout << "Warning: Using optional weights launches multiple rivet analyses. This may slow down your run substantially!" << endl; }
  CurrentGenerator::Redirect stdout(cout);
  _rivet = new Rivet::AnalysisHandler; //(fname);
  _rivet->addAnalyses(_analyses);
  if(useoptweights) { 
    for(int rr = 0; rr < 110; rr++) {
      OptXS.push_back(0.);
      _rivetMULTI[rr] = new Rivet::AnalysisHandler; //(fname);
      _rivetMULTI[rr]->addAnalyses(_analyses);
    }
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
