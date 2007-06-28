// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MSSM class.
//

#include "MSSM.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void MSSM::persistentOutput(PersistentOStream & os) const {
  os << theStopMix << theSbotMix << theStauMix << theAlpha 
     << theAtop << theAbottom << theAtau << theHiggsMix;
}

void MSSM::persistentInput(PersistentIStream & is, int) {
  is >> theStopMix >> theSbotMix >> theStauMix >> theAlpha 
     >> theAtop >> theAbottom >> theAtau >> theHiggsMix;
}

ClassDescription<MSSM> MSSM::initMSSM;
// Definition of the static class description member.

void MSSM::Init() {

  static ClassDocumentation<MSSM> documentation
    ("There is no documentation for the MSSM class");

}

void MSSM::createMixingMatrices() {
  map<string,pair<pair<unsigned int,unsigned int>, 
    vector<MixingElement> > >::const_iterator it;
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
    else if (name == "nmhmix") {
      createMixingMatrix(theHiggsMix,name,it->second.second,it->second.first);
    }
  }
  // neutral higgs mixing if not already set
  if(!theHiggsMix) {
    vector<MixingElement> hmix;
    hmix.push_back(MixingElement(1,1, cos(theAlpha)));
    hmix.push_back(MixingElement(1,2,-sin(theAlpha)));
    hmix.push_back(MixingElement(2,1, sin(theAlpha)));
    hmix.push_back(MixingElement(2,2, cos(theAlpha)));
    vector<long> ids(2);
    ids[0] = 25; ids[1] = 35;
    theHiggsMix = new_ptr(MixingMatrix(2,2));
    (*theHiggsMix).setIds(ids);
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
    theSbotMix->adjustPhase(id);
    if(theStopMix)
      theStopMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The stop mixing matrix pointer is null!" 
			     << Exception::runerror;
    break;
  case 1000015 :
  case 2000015 :
    if(theStopMix)
      theStopMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The stop mixing matrix pointer is null!" 
			     << Exception::runerror;
    break;
  default :
    SusyBase::adjustMixingMatrix(id);
    break;
  }
}

void MSSM::extractParameters(bool checkmodel) {
  map<string,pair<pair<unsigned int,unsigned int>, 
    vector<MixingElement> > >::const_iterator it;
  // trilinear couplings
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    vector<MixingElement>::const_iterator vit;
    if(name=="au") {
      theAtop=0.*GeV;
      for(vit=it->second.second.begin();vit!=it->second.second.end();++vit) {
	if(vit->row==3&&vit->col==3) theAtop=vit->value*GeV;
      }
    }
    else if(name=="ad") {
      theAbottom=0.*GeV;
      for(vit=it->second.second.begin();vit!=it->second.second.end();++vit) {
	if(vit->row==3&&vit->col==3) theAbottom=vit->value*GeV;
      }
    }
    else if(name=="ae") {
      theAtau=0.*GeV;
      for(vit=it->second.second.begin();vit!=it->second.second.end();++vit) {
	if(vit->row==3&&vit->col==3) theAtau=vit->value*GeV;
      }
    }
  }
  // the higgs mixing angle
  map<string,ParamMap>::const_iterator pit;
  theAlpha=0.;
  pit=parameters().find("alpha");
  if(pit!=parameters().end()) {
    ParamMap::const_iterator it = pit->second.find(1);
    if(it!=pit->second.end()) theAlpha=it->second;
  }
  // neutralino and chargino paramters in thew base class
  SusyBase::extractParameters(false);
  if(checkmodel) {
    map<string,ParamMap>::const_iterator pit;
    pit=parameters().find("modsel");
    if(pit==parameters().end()) return;
    ParamMap::const_iterator it;
    // nmssm or mssm
    it = pit->second.find(3);
    int inmssm = it!=pit->second.end() ? int(it->second) : 0;
    if(inmssm!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				    << " used but NMSSM read in " 
				    << Exception::runerror; 
    // RPV
    it = pit->second.find(4);
    int irpv = it!=pit->second.end() ? int(it->second) : 0;
    if(irpv!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				  << " used but RPV read in " 
				  << Exception::runerror; 
    // CPV
    it = pit->second.find(5);
    int icpv = it!=pit->second.end() ? int(it->second) : 0;
    if(icpv!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				  << " used but CPV read in " 
				  << Exception::runerror; 
    // flavour violation
    it = pit->second.find(6);
    int ifv = it!=pit->second.end() ? int(it->second) : 0;
    if(ifv!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				 << " used but flavour violation read in " 
				 << Exception::runerror;
  }
}
