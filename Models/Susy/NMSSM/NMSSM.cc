// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSM class.
//

#include "NMSSM.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"

using namespace Herwig;

void NMSSM::persistentOutput(PersistentOStream & os) const {
  os << _lambda << _kappa << ounit(_theAlambda,GeV) 
     << ounit(_theAkappa, GeV) << ounit(_lambdaVEV, GeV)
     << ounit(_MQ3, GeV) << ounit(_MU2, GeV);
}

void NMSSM::persistentInput(PersistentIStream & is, int) {
  is >> _lambda >> _kappa >> iunit(_theAlambda,GeV) 
     >> iunit(_theAkappa, GeV) >> iunit(_lambdaVEV, GeV)
     >> iunit(_MQ3, GeV) >> iunit(_MU2, GeV);
}

ClassDescription<NMSSM> NMSSM::initNMSSM;
// Definition of the static class description member.

void NMSSM::Init() {

  static ClassDocumentation<NMSSM> documentation
    ("The NMSSM class is the base class for the NMSSM model");

}

void NMSSM::extractParameters(bool checkmodel) {
  MSSM::extractParameters(false);
  if(checkmodel) {
    map<string,ParamMap>::const_iterator pit;
    pit = parameters().find("modsel");
    if(pit == parameters().end()) return;
    ParamMap::const_iterator it;
    // nmssm or mssm
    it = pit->second.find(3);
    int inmssm = (it != pit->second.end()) ? int(it->second) : 0;
    if(inmssm == 0) 
      throw Exception() << "R-parity, CP and flavour conserving NMSSM model"
			<< " used but MSSM read in." << Exception::runerror; 
    // RPV
    it = pit->second.find(4);
    int irpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(irpv != 0) throw Exception() << "NMSSM model does not support RPV"
				  << Exception::runerror; 
    // CPV
    it = pit->second.find(5);
    int icpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(icpv != 0) throw Exception() << "NMSSM model does not support CPV" 
				  << Exception::runerror; 
    // flavour violation
    it = pit->second.find(6);
    int ifv = (it != pit->second.end()) ? int(it->second) : 0;
    if(ifv != 0) throw Exception() << "NMSSM model does not support "
				 << "flavour violation"
				 << Exception::runerror;
  }
  // get the NMSSM parameters
  map<string,ParamMap>::const_iterator pit;
  pit=parameters().find("msoft");
  if( pit != parameters().end() ) {
    ParamMap::const_iterator it;
    it = pit->second.find(43);
    if(it != pit->second.end()) _MQ3 = it->second*GeV;
    it = pit->second.find(46);
    if(it != pit->second.end()) _MU2 = it->second*GeV;
  }
  pit=parameters().find("nmssmrun");
  if( pit != parameters().end() ) {
    ParamMap::const_iterator it = pit->second.find(1);
    if(it != pit->second.end()) _lambda = it->second;
    it = pit->second.find(2);
    if(it != pit->second.end()) _kappa = it->second;
    it = pit->second.find(3);
    if(it != pit->second.end()) _theAlambda = it->second*GeV;
    it = pit->second.find(4);
    if(it != pit->second.end()) _theAkappa = it->second*GeV;
    it = pit->second.find(5);
    if(it != pit->second.end()) _lambdaVEV = it->second*GeV;
  }
  pit=parameters().find("extpar");
  if( pit != parameters().end() ) {
    ParamMap::const_iterator it = pit->second.find(61);
    if(_lambda==ZERO     && it != pit->second.end()) _lambda = it->second;
    it = pit->second.find(62);
    if(_kappa==ZERO      && it != pit->second.end()) _kappa = it->second;
    it = pit->second.find(63);
    if(_theAlambda==ZERO && it != pit->second.end()) _theAlambda = it->second*GeV;
    it = pit->second.find(64);
    if(_theAkappa==ZERO  && it != pit->second.end()) _theAkappa = it->second*GeV;
    it = pit->second.find(65);
    if(_lambdaVEV==ZERO  && it != pit->second.end()) _lambdaVEV = it->second*GeV;
    it = pit->second.find(43);
    if(_MQ3==ZERO        && it != pit->second.end()) _MQ3 = it->second*GeV;
    it = pit->second.find(46);
    if(_MU2==ZERO        && it != pit->second.end()) _MU2 = it->second*GeV;
  }
  else {
    throw Exception() << "NMSSM::extractParameters - There was no EXTPAR block "
		      << "in the extracted parameters list. The model cannot "
		      << "be used without these." << Exception::runerror;
  }
  pit=parameters().find("msoft");
  if( pit != parameters().end() ) {
    ParamMap::const_iterator it;
    if(_MQ3==ZERO) {
      it = pit->second.find(43);
      if(it != pit->second.end()) _MQ3 = it->second*GeV;
    }
    if(_MU2==ZERO) {
      it = pit->second.find(46);
      if(it != pit->second.end()) _MU2 = it->second*GeV;
    }
  }
}

void NMSSM::createMixingMatrices() {
  map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    // pseudo-scalar higgs mixing
    if (name == "nmamix") {
      MixingMatrixPtr temp;
      createMixingMatrix(temp,name,it->second.second,it->second.first);
      CPoddHiggsMix(temp);
    }
  }
  // base class for neutralinos and charginos
  MSSM::createMixingMatrices();
}
