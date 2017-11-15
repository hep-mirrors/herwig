// -*- C++ -*-
//
// MSSM.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MSSM class.
//

#include "MSSM.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void MSSM::persistentOutput(PersistentOStream & os) const {
  os << theStopMix << theSbotMix << theStauMix << theAlpha 
     << ounit(theAtop,GeV) << ounit(theAbottom,GeV) << ounit(theAtau,GeV) 
     << theHiggsMix << HiggsAMix_ << HiggsPMix_ << createDiagonalMixing_;
}

void MSSM::persistentInput(PersistentIStream & is, int) {
  is >> theStopMix >> theSbotMix >> theStauMix >> theAlpha 
     >> iunit(theAtop,GeV) >> iunit(theAbottom,GeV) >> iunit(theAtau,GeV) 
     >> theHiggsMix >> HiggsAMix_ >> HiggsPMix_ >>  createDiagonalMixing_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MSSM,SusyBase>
describeMSSM("Herwig::MSSM", "HwSusy.so");

void MSSM::Init() {

  static ClassDocumentation<MSSM> documentation
    ("The MSSM class is the base class for the MSSM model.",
     "MSSM Feynman rules were taken from \\cite{Haber:1984rc,Gunion:1984yn}.",
     " %\\cite{Haber:1984rc}\n"
     "\\bibitem{Haber:1984rc}\n"
     "  H.~E.~Haber and G.~L.~Kane,\n"
     "  %``The Search For Supersymmetry: Probing Physics Beyond The Standard Model,''\n"
     "  Phys.\\ Rept.\\  {\\bf 117}, 75 (1985).\n"
     "  %%CITATION = PRPLC,117,75;%%\n"
     "%\\cite{Gunion:1984yn}\n"
     "\\bibitem{Gunion:1984yn}\n"
     "  J.~F.~Gunion and H.~E.~Haber,\n"
     "  %``Higgs Bosons In Supersymmetric Models. 1,''\n"
     "  Nucl.\\ Phys.\\  B {\\bf 272}, 1 (1986)\n"
     "  [Erratum-ibid.\\  B {\\bf 402}, 567 (1993)].\n"
     "  %%CITATION = NUPHA,B272,1;%%\n"
    );

  static Switch<MSSM,bool> interfaceCreateDiagonalMixingMatrices
    ("CreateDiagonalMixingMatrices",
     "Create diagonal stop, sbottom and stau mixings if not present.",
     &MSSM::createDiagonalMixing_, false, false, false);
  static SwitchOption interfaceCreateDiagonalMixingMatricesNo
    (interfaceCreateDiagonalMixingMatrices,
     "No",
     "Don't create them",
     false);
  static SwitchOption interfaceCreateDiagonalMixingMatricesYes
    (interfaceCreateDiagonalMixingMatrices,
     "Yes",
     "Create them if needed",
     true);

}

void MSSM::createMixingMatrices() {
  useMe();
  map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    // create the stop, sbottom and stau mixing matrices
    if(name == "stopmix"  ){
      createMixingMatrix(theStopMix,name,it->second.second,it->second.first);
    }
    else if (name == "sbotmix" ) {
      createMixingMatrix(theSbotMix,name,it->second.second,it->second.first);
    }
    else if (name == "staumix") {
      createMixingMatrix(theStauMix,name,it->second.second,it->second.first);
    }
    // Higgs mixing matrix in extended models
    else if (name == "nmhmix" || name == "rvhmix") {
      createMixingMatrix(theHiggsMix,name,it->second.second,it->second.first);
    }
  }
  // create stop, sbottom and stau mixing if needed and absent
  if(createDiagonalMixing_) {
    // stop
    if(!theStopMix) {
      theStopMix = new_ptr(MixingMatrix(2,2));
      (*theStopMix)(0,0) = 1.;
      (*theStopMix)(1,1) = 1.;
    }
    // sbottom
    if(!theSbotMix) {
      theSbotMix = new_ptr(MixingMatrix(2,2));
      (*theSbotMix)(0,0) = 1.;
      (*theSbotMix)(1,1) = 1.;
    }
    // stau
    if(!theStauMix) {
      theStauMix = new_ptr(MixingMatrix(2,2));
      (*theStauMix)(0,0) = 1.;
      (*theStauMix)(1,1) = 1.;
    }
  }
  // neutral higgs mixing if not already set
  if(!theHiggsMix) {
    MixingVector hmix;
    hmix.push_back(MixingElement(2,1, cos(theAlpha)));
    hmix.push_back(MixingElement(2,2, sin(theAlpha)));
    hmix.push_back(MixingElement(1,1,-sin(theAlpha)));
    hmix.push_back(MixingElement(1,2, cos(theAlpha)));
    vector<long> ids(2);
    ids[0] = 25; ids[1] = 35;
    theHiggsMix = new_ptr(MixingMatrix(2,2));
    (*theHiggsMix).setIds(ids);
    for(unsigned int ix=0; ix < hmix.size(); ++ix)
      (*theHiggsMix)(hmix[ix].row-1,hmix[ix].col-1) = hmix[ix].value;
    hmix.clear();
    double beta = atan(tanBeta());
    hmix.push_back(MixingElement(1,1,sin(beta)));
    hmix.push_back(MixingElement(1,2,cos(beta)));
    ids.clear();
    ids.resize(1,37);
    HiggsPMix_ = new_ptr(MixingMatrix(1,2));
    (*HiggsPMix_).setIds(ids);
    for(unsigned int ix=0; ix < hmix.size(); ++ix)
      (*HiggsPMix_)(hmix[ix].row-1,hmix[ix].col-1) = hmix[ix].value;
    ids.clear();
    ids.resize(1,36);
    HiggsAMix_ = new_ptr(MixingMatrix(1,2));
    (*HiggsAMix_).setIds(ids);
    for(unsigned int ix=0; ix < hmix.size(); ++ix)
      (*HiggsAMix_)(hmix[ix].row-1,hmix[ix].col-1) = hmix[ix].value;
  }
  // base class for neutralinos and charginos
  SusyBase::createMixingMatrices();
}

void MSSM::adjustMixingMatrix(long id) {
  switch (id) {
  case 1000006 :
  case 2000006 :
    if(theStopMix)
      theStopMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The stop mixing matrix pointer is null!" 
			     << Exception::runerror;
    break;
  case 1000005 :
  case 2000005 :
    if(theSbotMix)
      theSbotMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The sbottom mixing matrix pointer is null!" 
			     << Exception::runerror;
    break;
  case 1000015 :
  case 2000015 :
    if(theStauMix)
      theStauMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The stau mixing matrix pointer is null!" 
			     << Exception::runerror;
    break;
  default :
    SusyBase::adjustMixingMatrix(id);
    break;
  }
}

void MSSM::extractParameters(bool checkmodel) {
  map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
  // trilinear couplings
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    MixingVector::const_iterator vit;
    if(name=="au") {
      theAtop=ZERO;
      for(vit=it->second.second.begin();vit!=it->second.second.end();++vit) {
	if(vit->row==3&&vit->col==3) theAtop=vit->value*GeV;
      }
    }
    else if(name=="ad") {
      theAbottom=ZERO;
      for(vit=it->second.second.begin();vit!=it->second.second.end();++vit) {
	if(vit->row==3&&vit->col==3) theAbottom=vit->value*GeV;
      }
    }
    else if(name=="ae") {
      theAtau=ZERO;
      for(vit=it->second.second.begin();vit!=it->second.second.end();++vit) {
	if(vit->row==3&&vit->col==3) theAtau=vit->value*GeV;
      }
    }
  }
  // neutralino and chargino paramters in the base class
  SusyBase::extractParameters(false);
  // check the model
  map<string,ParamMap>::const_iterator pit;
  pit=parameters().find("modsel");
  if(pit==parameters().end()) return;
  // nmssm or mssm
  ParamMap::const_iterator jt;
  jt = pit->second.find(3);
  int inmssm = jt!=pit->second.end() ? int(jt->second) : 0; 
  // RPV
  jt = pit->second.find(4);
  int irpv = jt!=pit->second.end() ? int(jt->second) : 0;
  // CPV
  jt = pit->second.find(5);
  int icpv = jt!=pit->second.end() ? int(jt->second) : 0;
  // flavour violation
  jt = pit->second.find(6);
  int ifv = jt!=pit->second.end() ? int(jt->second) : 0;
  // the higgs mixing angle 
  theAlpha=0.;
  bool readAlpha = false;
  pit=parameters().find("alpha");
  if(pit!=parameters().end()) {
    ParamMap::const_iterator it = pit->second.find(1);
    if(it!=pit->second.end()) {
      readAlpha = true;
      theAlpha=it->second;
    }
  }
  if(inmssm==0&&irpv==0&&!readAlpha) 
    throw Exception() << "In the MSSM model BLOCK ALPHA which must be"
		      << " present in the SLHA file is missing"
		      << Exception::runerror;
  if(checkmodel) {
    if(inmssm!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				    << " used but NMSSM read in " 
				    << Exception::runerror;
    if(irpv!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				  << " used but RPV read in " 
				  << Exception::runerror; 
    if(icpv!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				  << " used but CPV read in " 
				  << Exception::runerror; 
    if(ifv!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				 << " used but flavour violation read in " 
				 << Exception::runerror;
  }
}
