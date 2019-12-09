// -*- C++ -*-
//
// HiggsPair.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2009-2019 The Herwig Collaboration
//
// Herwig is licenaaaced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiggsPair class.
//

#include "HiggsPair.h"
#include "MEHiggsPairJet.h"
#include "MEHiggsPairOL.h"

#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/Utilities/Maths.h"
#include <fstream>
#include <cmath>


using namespace Herwig;

void HiggsPair::doinitrun()  {  
 _theModel = generator()->standardModel();
  tcHiggsPairPtr hwHiggsPair=dynamic_ptr_cast<tcHiggsPairPtr>(_theModel);
  if(hwHiggsPair){
    _selfcoupling=hwHiggsPair->selfcoupling();
    _process=hwHiggsPair->process();
    _includeBquark=hwHiggsPair->includeBquark();
    _includeWidths=hwHiggsPair->includeWidths();
    _implementation=hwHiggsPair->implementation();
    _fixedscale=hwHiggsPair->fixedscale();
    _alphascale=hwHiggsPair->alphascale();
    _fixedalphaS=hwHiggsPair->fixedalphaS();
    _alphasfixedvalue = hwHiggsPair->alphasfixedvalue();
    _basescale=hwHiggsPair->basescale();
    _scalemultiplier=hwHiggsPair->scalemultiplier();
    
  }
  //cout << _selfcoupling << " " << _process << " " << _includeBquark << " " << _includeWidths << " " << _implementation << " " << "- " << _fixedscale << " " << _alphascale << " " << _fixedalphaS << endl;
  //cout << _basescale << endl;
  // cout << _alphascale/GeV << endl;

  if(_implementation == 0 && _process != 0) { cerr << "The OpenLoops implementation currently only contains SM production" << endl; exit(1); }

  if(_process == 4 || _process == 5) { cerr << "HH and hH production not implemented yet, please choose an hh production subprocess" << endl; exit(1); } 

  _higgs = getParticleData(ParticleID::h0);
  PDPtr top = getParticleData(ParticleID::t);
  _topmass = top->mass();
  //  cout << "_topmass = " << _topmass << endl;
  PDPtr zboson = getParticleData(ParticleID::Z0);
  PDPtr wboson = getParticleData(24);

  _zmass = zboson->mass();
  PDPtr bottom = getParticleData(ParticleID::b);
  _bottommass = bottom->mass();
  PDPtr heavyH = getParticleData(35);
  _heavyHmass = heavyH->mass();
  _heavyHwidth = heavyH->width();
  
  Mass_E = 0.0;
  Mass_M = 0.0;
  Mass_L = 0.0;

  Mass_T = top->mass()/GeV;
  Mass_U = 0.0;
  Mass_C = 0.0;
  Mass_D = 0.0;
  Mass_S=  0.0;
  Mass_B = bottom->mass()/GeV;
  if(_includeBquark == 0) { Mass_B = 0; cout << "Not including B quark mass." << endl; }

  Mass_Z = zboson->mass()/GeV;
  Mass_W = wboson->mass()/GeV;
  Mass_H = _higgs->mass()/GeV;

  Width_C = 0.0;
  Width_B = 0.0;
  Width_T = 0.0;
  Width_W = 0.0;
  Width_Z = 0.0;
  Width_H = 0.0;
  if(_includeWidths == 1 && _implementation == 0) {
    Width_T = top->width()/GeV;
    Width_W = wboson->width()/GeV;
    Width_Z = zboson->width()/GeV;
    Width_H = _higgs->width()/GeV;
  } 
  if(_includeWidths == 1 && _implementation != 0) {
    cout << "Warning: top quark width effects are only included in the OpenLoops implementation" << endl;
  }

  Coupl_Alpha_QED= 1/137.; //change to correct scale.
  Coupl_Alpha_QCD = 0.12360931820752918;//SM().alphaS(sqr(2*Mass_H*GeV)); //change to correct scale
  
  InitOpenLoops();

  _mh = _higgs->mass();
  _m1 = _higgs->mass();
  _m2 = _higgs->mass();
  _wh = _higgs->width();
    
}

void HiggsPair::doinit()  {  
  _theModel = generator()->standardModel();
  tcHiggsPairPtr hwHiggsPair=dynamic_ptr_cast<tcHiggsPairPtr>(_theModel);
  if(hwHiggsPair){
    _selfcoupling=hwHiggsPair->selfcoupling();
    _process=hwHiggsPair->process();
    _includeBquark=hwHiggsPair->includeBquark();
    _includeWidths=hwHiggsPair->includeWidths();
    _implementation=hwHiggsPair->implementation();
    _fixedscale=hwHiggsPair->fixedscale();
    _alphascale=hwHiggsPair->alphascale();
    _fixedalphaS=hwHiggsPair->fixedalphaS();
    _alphasfixedvalue = hwHiggsPair->alphasfixedvalue();
    _basescale=hwHiggsPair->basescale();
    _scalemultiplier=hwHiggsPair->scalemultiplier();
    
  }
  //cout << _selfcoupling << " " << _process << " " << _includeBquark << " " << _includeWidths << " " << _implementation << " " << "- " << _fixedscale << " " << _alphascale << " " << _fixedalphaS << endl;
  //cout << _basescale << endl;
  // cout << _alphascale/GeV << endl;

  if(_implementation == 0 && _process != 0) { cerr << "The OpenLoops implementation currently only contains SM production" << endl; exit(1); }

  if(_process == 4 || _process == 5) { cerr << "HH and hH production not implemented yet, please choose an hh production subprocess" << endl; exit(1); } 

  _higgs = getParticleData(ParticleID::h0);
  PDPtr top = getParticleData(ParticleID::t);
  _topmass = top->mass();
  //  cout << "_topmass = " << _topmass << endl;
  PDPtr zboson = getParticleData(ParticleID::Z0);
  PDPtr wboson = getParticleData(24);

  _zmass = zboson->mass();
  PDPtr bottom = getParticleData(ParticleID::b);
  _bottommass = bottom->mass();
  PDPtr heavyH = getParticleData(35);
  _heavyHmass = heavyH->mass();
  _heavyHwidth = heavyH->width();
  
  Mass_E = 0.0;
  Mass_M = 0.0;
  Mass_L = 0.0;

  Mass_T = top->mass()/GeV;
  Mass_U = 0.0;
  Mass_C = 0.0;
  Mass_D = 0.0;
  Mass_S=  0.0;
  Mass_B = bottom->mass()/GeV;
  if(_includeBquark == 0) { Mass_B = 0; cout << "Not including B quark mass." << endl; }

  Mass_Z = zboson->mass()/GeV;
  Mass_W = wboson->mass()/GeV;
  Mass_H = _higgs->mass()/GeV;

  Width_C = 0.0;
  Width_B = 0.0;
  Width_T = 0.0;
  Width_W = 0.0;
  Width_Z = 0.0;
  Width_H = 0.0;
  if(_includeWidths == 1 && _implementation == 0) {
    Width_T = top->width()/GeV;
    Width_W = wboson->width()/GeV;
    Width_Z = zboson->width()/GeV;
    Width_H = _higgs->width()/GeV;
  } 
  if(_includeWidths == 1 && _implementation != 0) {
    cout << "Warning: top quark width effects are only included in the OpenLoops implementation" << endl;
  }

  Coupl_Alpha_QED= 1/137.; //change to correct scale.
  Coupl_Alpha_QCD = 0.12360931820752918;//SM().alphaS(sqr(2*Mass_H*GeV)); //change to correct scale
  // cout << "init init init" << endl;
  //  if(_implementation == 0 || _implementation == 2) { InitOpenLoops(); }
  //  cout << "init init init 2" << endl;

  InitOpenLoops();

  _mh = _higgs->mass();
  _m1 = _higgs->mass();
  _m2 = _higgs->mass();
  _wh = _higgs->width();

  if(_higgs->massGenerator()) {
    _hmass=dynamic_ptr_cast<GenericMassGeneratorPtr>(_higgs->massGenerator());
  }

  BSMModel::doinit();
}

HiggsPair::HiggsPair()
  : _selfcoupling(1.0), _process(0), _fixedalphaS(0), _fixedscale(0), _alphascale(100.*GeV),
    _implementation(0),
    _includeWidths(0), _includeBquark(1), _maxflavour(5), _basescale(125.*GeV), _subprocessjet(0), _alphasreweight(0), _scalemultiplier(1.0)
   
{
}



IBPtr HiggsPair::clone() const {
    return new_ptr(*this);
  }
  
IBPtr HiggsPair::fullclone() const {
    return new_ptr(*this);
  }

void HiggsPair::persistentOutput(PersistentOStream & os) const {
  os << _selfcoupling << _process << _implementation << _fixedalphaS << _includeWidths << ounit(_alphascale,GeV) << _fixedscale << _includeBquark << _maxflavour << ounit(_basescale,GeV) << _subprocessjet << _alphasreweight << _scalemultiplier;
}

void HiggsPair::persistentInput(PersistentIStream & is, int) {
  is >> _selfcoupling  >> _process >> _implementation >> _fixedalphaS >> _includeWidths >> iunit(_alphascale,GeV) >> _fixedscale >> _includeBquark >> _maxflavour >> iunit(_basescale,GeV) >> _subprocessjet >> _alphasreweight >> _scalemultiplier;

}


// Definition of the static class description member.
ClassDescription<HiggsPair> HiggsPair::initHiggsPair;

void HiggsPair::Init() {
  static ClassDocumentation<HiggsPair> documentation
    ("The HiggsPair class implements the hh process in hadron-hadron"
     " collisions");

 static Parameter<HiggsPair, double> interfaceSelfCoupling
    ("SelfCoupling",
     "Multiplier for the SM Higgs triple coupling",
     &HiggsPair::_selfcoupling, 1.0, -10.0, 10.0,
     false, false, Interface::limited);
 
 //choose whether to include the triangle, box or both
  static Switch<HiggsPair,unsigned int> interfaceProcess
    ("Process",
     "Which processes to include",
     &HiggsPair::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all SM pp->hh subprocesses",
     0);
  static SwitchOption interfaceProcessTriangle
    (interfaceProcess,
     "ggToTriangleTohh",
     "Include only SM gg->hh triangle subprocess",
     1);
  static SwitchOption interfaceProcessBox
    (interfaceProcess,
     "ggToBoxTohh",
     "Include only SM gg->hh box subprocess",
     2);

  //choose between OpenLoops or HPAIR implementations
  static Switch<HiggsPair,unsigned int> interfaceImplementation
    ("Implementation",
     "Which implementation to use",
     &HiggsPair::_implementation, 0, false, false);
  static SwitchOption interfaceImplementationOpenLoops
    (interfaceImplementation,
     "OpenLoops",
     "Use the OpenLoops implementation",
     0);
  static SwitchOption interfaceImplementationHPAIR
    (interfaceImplementation,
     "HPAIR",
     "Use the HPAIR implementation",
     1);
  static SwitchOption interfaceImplementationTesting
    (interfaceImplementation,
     "Testing",
     "Calculate both implementations, return the HPAIR one",
     2);

  static Switch<HiggsPair,unsigned int> interfaceFixedAlphaS
    ("FixedAlphaS",
     "Which implementation to use",
     &HiggsPair::_fixedalphaS, 0, false, false);
  static SwitchOption interfaceFixedAlphaSNo
    (interfaceFixedAlphaS,
     "No",
     "Don't fix the scale of alphaS.",
     0);
  static SwitchOption interfaceFixedAlphaSScale
    (interfaceFixedAlphaS,
     "Scale",
     "Fix the scale of alphaS.",
     1);
  static SwitchOption interfaceFixedAlphaSAlphaS
    (interfaceFixedAlphaS,
     "AlphaS",
     "Fix alphaS directly via the AlphaS interface.",
     2);


 static Parameter<HiggsPair, double> interfaceAlphaS
    ("AlphaS",
     "The value of AlphaS if FixedAlphaS is AlphaS, option 2",
     &HiggsPair::_alphasfixedvalue, 0.12361, 0., 100000.0,
     false, false, Interface::limited);

  static Parameter<HiggsPair, Energy> interfaceAlphaSScale
    ("AlphaSScale",
     "Scale for alphaS if fixed",
     &HiggsPair::_alphascale, GeV, 100.0*GeV, 0.0*GeV, 10000000.0*GeV,
     false, false, Interface::limited);

  static Parameter<HiggsPair, Energy> interfaceBaseScale
    ("BaseScale",
     "Base scale if fixed",
     &HiggsPair::_basescale, GeV, 125.0*GeV, 0.0*GeV, 10000000.0*GeV,
     false, false, Interface::limited);

  static Switch<HiggsPair,unsigned int> interfaceIncludeB
    ("IncludeBottomLoop",
     "Include bottom quark loops or not",
     &HiggsPair::_includeBquark, 1, false, false);
  static SwitchOption interfaceIncludeBNo
    (interfaceIncludeB,
     "No",
     "Don't include b quark loops.",
     0);
  static SwitchOption interfaceIncludeBYes
    (interfaceIncludeB,
     "Yes",
     "Include b quark loops.",
     1);


  static Switch<HiggsPair,unsigned int> interfaceAlphaSReweighting
    ("AlphaSReweighting",
     "Include AlphaS reweighting or not",
     &HiggsPair::_alphasreweight, 0, false, false);
  static SwitchOption interfaceAlphaSReweightNo
    (interfaceAlphaSReweighting,
     "No",
     "Don't do AlphaS reweighting",
     0);
  static SwitchOption interfaceeAlphaSReweightYes
    (interfaceAlphaSReweighting,
     "Yes",
     "AlphaS reweighting enabled.",
     1);


  static Switch<HiggsPair,unsigned int> interfaceFixedScale
    ("FixedScale",
     "Whether to use fixed base scale or not",
     &HiggsPair::_fixedscale, 0, false, false);
  static SwitchOption interfaceFixedScaleNo
    (interfaceFixedScale,
     "sHat",
     "Don't fix the base scale, use sHat().",
     0);
  static SwitchOption interfaceFixedScaleFixed
    (interfaceFixedScale,
     "Fixed",
     "Fix the scale of the process to chosen (base scale)^2.",
     2);
  static SwitchOption interfaceFixedScaleBaseQuad
    (interfaceFixedScale,
     "FixedQuadPt",
     "Fix the scale of the process to (base scale)^2 + pt^2.",
     1);
  static SwitchOption interfaceFixedScaleBaseLin
    (interfaceFixedScale,
     "FixedLinPt",
     "Fix the scale of the process to (base scale + pt)^2.",
     3);
  static SwitchOption interfaceFixedScaleMHH
    (interfaceFixedScale,
     "MHH",
     "Don't fix the base scale, use MHH.",
     5);

  static Parameter<HiggsPair, double> interfaceScaleMultiplier
    ("ScaleMultiplier",
     "Multiplier scale",
     &HiggsPair::_scalemultiplier, 1.0, 0.0, 100000.0,
     false, false, Interface::limited);


  static Switch<HiggsPair,unsigned int> interfaceIncludeWidths
    ("IncludeWidths",
     "Include top/W/Z widths or not",
     &HiggsPair::_includeWidths, 0, false, false);
  static SwitchOption interfaceIncludeWidthsNo
    (interfaceIncludeWidths,
     "No",
     "Don't include widths.",
     0);
  static SwitchOption interfaceIncludeWidthsYes
    (interfaceIncludeWidths,
     "Yes",
     "Include widths.",
     1);

  static Parameter<HiggsPair,int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &HiggsPair::_maxflavour, 5, 0, 8, false, false, true);

  static Switch<HiggsPair,unsigned int> interfaceSubProcessJet
    ("SubProcessJet",
     "Which subprocesses to include in the HHj calculation",
     &HiggsPair::_subprocessjet, 0, false, false);
  static SwitchOption interfaceSubProcessJetAll
    (interfaceSubProcessJet,
     "All",
     "Include all SM pp->hhj subprocesses",
     0);
  static SwitchOption interfaceSubProcessJetgghhj
    (interfaceSubProcessJet,
     "gghhj",
     "Include only SM gg->hhj subprocess",
     1);
  static SwitchOption interfaceSubProcessJetqxghhqx
    (interfaceSubProcessJet,
     "qxghhqx",
     "Include only SM qxg->hhqx HH subprocess",
     2);
  static SwitchOption interfaceSubProcessJetqghhq
    (interfaceSubProcessJet,
     "qghhq",
     "Include only SM qg->hhq HH subprocess",
     3);
  static SwitchOption interfaceSubProcessJetqqxhhg
    (interfaceSubProcessJet,
     "qqxhhg",
     "Include only SM qqx->hhg HH subprocess",
     4);
}


void HiggsPair::InitOpenLoops() {

  ol_setparameter_double("lambda_hhh" , _selfcoupling);

  //ol_setparameter_int("order_ew", 2);
  ol_setparameter_int("redlib1",1);
  ol_setparameter_int("redlib2",7);
  ol_setparameter_int("verbose",2);

  ol_setparameter_int("stability_mode",21);
  ol_setparameter_double("mass(23)", Mass_Z);
  ol_setparameter_double("mass(24)", Mass_W);
  ol_setparameter_double("mass(25)", Mass_H);
  ol_setparameter_double("mass(1)", Mass_D);
  ol_setparameter_double("mass(2)", Mass_U);
  ol_setparameter_double("mass(3)", Mass_S);
  ol_setparameter_double("mass(4)", Mass_C);
  ol_setparameter_double("mass(5)", Mass_B);
  ol_setparameter_double("mass(6)", Mass_T);
  
  ol_setparameter_double("width(23)", Width_Z);
  ol_setparameter_double("width(24)", Width_W);
  ol_setparameter_double("width(25)", Width_H);

  // Initialize LO process
  _id = ol_register_process("21 21 -> 25 25", 12);
  char only3h[] = "only3h";
  ol_setparameter_string("approx", only3h);
  _idtri = ol_register_process("21 21 -> 25 25", 12);
  char no3h[] = "no3h";
  ol_setparameter_string("approx", no3h);
  _idbox = ol_register_process("21 21 -> 25 25", 12);
  char interf3h[] = "interf3h";
  ol_setparameter_string("approx", interf3h);
  _idint = ol_register_process("21 21 -> 25 25", 12);

  // Initialize real emission process
  ol_setparameter_string("approx", "");
  _id1 = ol_register_process("2 21 -> 2 25 25", 12);
  _id2 = ol_register_process("-2 21 -> -2 25 25", 12);
  _id3 = ol_register_process("21 21 -> 21 25 25", 12);
  _id4 = ol_register_process("2 -2 -> 21 25 25", 12);
  ol_setparameter_string("approx", only3h);
  _id1tri = ol_register_process("2 21 -> 2 25 25", 12);
  _id2tri = ol_register_process("-2 21 -> -2 25 25", 12);
  _id3tri = ol_register_process("21 21 -> 21 25 25", 12);
  _id4tri = ol_register_process("2 -2 -> 21 25 25", 12);
  ol_setparameter_string("approx", no3h);
  _id1box = ol_register_process("2 21 -> 2 25 25", 12);
  _id2box = ol_register_process("-2 21 -> -2 25 25", 12);
  _id3box = ol_register_process("21 21 -> 21 25 25", 12);
  _id4box = ol_register_process("2 -2 -> 21 25 25", 12); 
  ol_setparameter_string("approx", interf3h);
  _id1int = ol_register_process("2 21 -> 2 25 25", 12);
  _id2int = ol_register_process("-2 21 -> -2 25 25", 12);
  _id3int = ol_register_process("21 21 -> 21 25 25", 12);
  _id4int = ol_register_process("2 -2 -> 21 25 25", 12); 

  // Initialize OpenLoops
  ol_start();
}
