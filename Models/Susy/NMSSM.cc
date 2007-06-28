// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSM class.
//

#include "NMSSM.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void NMSSM::persistentOutput(PersistentOStream & os) const {
  os << theHiggsAMix;
}

void NMSSM::persistentInput(PersistentIStream & is, int) {
  is >> theHiggsAMix;
}

ClassDescription<NMSSM> NMSSM::initNMSSM;
// Definition of the static class description member.

void NMSSM::Init() {

  static ClassDocumentation<NMSSM> documentation
    ("There is no documentation for the NMSSM class");

}

void NMSSM::extractParameters(bool checkmodel) {
  SusyBase::extractParameters(false);
  if(checkmodel) {
    map<string,ParamMap>::const_iterator pit;
    pit=parameters().find("modsel");
    if(pit==parameters().end()) return;
    ParamMap::const_iterator it;
    // nmssm or mssm
    it = pit->second.find(3);
    int inmssm = it!=pit->second.end() ? int(it->second) : 0;
    if(inmssm==0) throw Exception() << "R-parity, CP and flavour conserving NMSSM model"
				    << " used but MSSM read in " 
				    << Exception::runerror; 
    // RPV
    it = pit->second.find(4);
    int irpv = it!=pit->second.end() ? int(it->second) : 0;
    if(irpv!=0) throw Exception() << "NMSSM model does not support RPV"
				  << Exception::runerror; 
    // CPV
    it = pit->second.find(5);
    int icpv = it!=pit->second.end() ? int(it->second) : 0;
    if(icpv!=0) throw Exception() << "NMSSM model does not support CPV" 
				  << Exception::runerror; 
    // flavour violation
    it = pit->second.find(6);
    int ifv = it!=pit->second.end() ? int(it->second) : 0;
    if(ifv!=0) throw Exception() << "NMSSM model does not support "
				 << "flavour violation"
				 << Exception::runerror;
  }
}

void NMSSM::createMixingMatrices() {
  map<string,pair<pair<unsigned int,unsigned int>, 
    vector<MixingElement> > >::const_iterator it;
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    // pseudo-scalar higgs mixing
    if (name == "nmamix") {
      cerr << "testing amix\n";
      createMixingMatrix(theHiggsAMix,name,it->second.second,it->second.first);
    }
  }
  // base class for neutralinos and charginos
  MSSM::createMixingMatrices();
}
