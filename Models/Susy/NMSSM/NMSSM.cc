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
  os << theHiggsAMix << _lambda << _kappa << ounit(_theAlambda,GeV) 
     << ounit(_theAkappa, GeV) << ounit(_lambdaVEV, GeV) << _ffhvertex 
     << _wwhvertex << _whhvertex << _gogohvertex << _hhhvertex 
     << _hssvertex;
}

void NMSSM::persistentInput(PersistentIStream & is, int) {
  is >> theHiggsAMix >> _lambda >> _kappa >> iunit(_theAlambda,GeV) 
     >> iunit(_theAkappa, GeV) >> iunit(_lambdaVEV, GeV) >> _ffhvertex 
     >> _wwhvertex >> _whhvertex >> _gogohvertex >> _hhhvertex 
     >> _hssvertex;
}

ClassDescription<NMSSM> NMSSM::initNMSSM;
// Definition of the static class description member.

void NMSSM::Init() {

  static ClassDocumentation<NMSSM> documentation
    ("The NMSSM class is the base class for the NMSSM model");
  
  static Reference<NMSSM,FFSVertex> interfaceVertexNMSSMFFH
    ("Vertex/NMSSMFFH",
     "The higgs coupling to SM fermions in the NMSSM",
     &NMSSM::_ffhvertex, false, false, true, false, false);
  
  static Reference<NMSSM,VVSVertex> interfaceVertexNMSSMWWH
    ("Vertex/NMSSMWWH",
     "The coupling of 2 EW gauge bosons a higgs in the NMSSM",
     &NMSSM::_wwhvertex, false, false, true, false, false);

  static Reference<NMSSM,VSSVertex> interfaceVertexNMSSMWHH
    ("Vertex/NMSSMWHH",
     "The coupling of a pair of Higgs to EW gauge bosons",
     &NMSSM::_whhvertex, false, false, true, false, false);

  static Reference<NMSSM,FFSVertex> interfaceVertexNMSSMGOGOH
    ("Vertex/NMSSMGOGOH",
     "The coupling of a pair of gauginos to a higgs boson in the NMSSM",
     &NMSSM::_gogohvertex, false, false, true, false, false);

  static Reference<NMSSM,SSSVertex> interfaceVertexNMSSMHHH
    ("Vertex/NMSSMHHH",
     "The triple higgs coupling in the NMSSM",
     &NMSSM::_hhhvertex, false, false, true, false, false);

  static Reference<NMSSM,SSSVertex> interfaceVertexNMSSMHSS
    ("Vertex/NMSSMHSS",
     "The coupling of a pair of sfermions to a higgs in the NMSSM",
     &NMSSM::_hssvertex, false, false, true, false, false);

  static Reference<NMSSM,FFSVertex> interfaceVertexNMSSMNFS
    ("Vertex/NMSSMNFS",
     "The coupling of a neutralino to a fermion and an sfermion",
     &NMSSM::_nfsvertex, false, false, true, false, false);

  static Reference<NMSSM,FFVVertex> interfaceVertexNMSSMNNZ
    ("Vertex/NMSSMNNZ",
     "The coulping of a pair of neutralinos to a Z-boson in the NMSSM ",
     &NMSSM::_nnzvertex, false, false, true, false, false);
  
  static Reference<NMSSM,FFVVertex> interfaceVertexNMSSMCNW
    ("Vertex/NMSSMCNW",
     "The coupling of a chargino to a neutralino and W boson in the NMSSM",
     &NMSSM::_cnwvertex, false, false, true, false, false);
}

void NMSSM::extractParameters(bool checkmodel) {
  SusyBase::extractParameters(false);
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
  pit=parameters().find("extpar");
  if( pit != parameters().end() ) {
    ParamMap::const_iterator it = pit->second.find(61);
    if(it != pit->second.end()) _lambda = it->second;
    it = pit->second.find(62);
    if(it != pit->second.end()) _kappa = it->second;
    it = pit->second.find(63);
    if(it != pit->second.end()) _theAlambda = it->second*GeV;
    it = pit->second.find(64);
    if(it != pit->second.end()) _theAkappa = it->second*GeV;
    it = pit->second.find(65);
    if(it != pit->second.end()) _lambdaVEV = it->second*GeV;
  }
  else {
    throw Exception() << "NMSSM::extractParameters - There was no EXTPAR block "
		      << "in the extracted parameters list. The model cannot "
		      << "be used without these." << Exception::runerror;
  }
}

void NMSSM::createMixingMatrices() {
  map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
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
