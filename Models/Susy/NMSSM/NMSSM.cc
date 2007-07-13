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
  os << theHiggsAMix << _lambda << _kappa;
}

void NMSSM::persistentInput(PersistentIStream & is, int) {
  is >> theHiggsAMix >> _lambda >> _kappa;
}

ClassDescription<NMSSM> NMSSM::initNMSSM;
// Definition of the static class description member.

void NMSSM::Init() {

  static ClassDocumentation<NMSSM> documentation
    ("The NMSSM class is the base class for the NMSSM model");

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
  // get the NMSSM parameters
  map<string,ParamMap>::const_iterator pit;
  pit=parameters().find("extpar");
  _lambda=0.;
  _kappa =0.;
  if(pit!=parameters().end()) {
    ParamMap::const_iterator it = pit->second.find(61);
    if(it!=pit->second.end()) _lambda=it->second;
    it = pit->second.find(62);
    if(it!=pit->second.end()) _kappa=it->second;
  }
}

void NMSSM::createMixingMatrices() {
  map<string,pair<pair<unsigned int,unsigned int>, 
    vector<MixingElement> > >::const_iterator it;
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    // pseudo-scalar higgs mixing
    if (name == "nmamix") {
      createMixingMatrix(theHiggsAMix,name,it->second.second,it->second.first);
    }
  }
  // base class for neutralinos and charginos
  MSSM::createMixingMatrices();
}
