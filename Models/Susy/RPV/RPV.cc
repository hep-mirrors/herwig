// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPV class.
//

#include "RPV.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RPV::persistentOutput(PersistentOStream & os ) const {
  os << lambdaLLE_ << lambdaLQD_ << lambdaUDD_ << ounit(vnu_,GeV)
     << LLEVertex_ << LQDVertex_ << UDDVertex_
     << HiggsAMix_ << HiggsPMix_;
}

void RPV::persistentInput(PersistentIStream & is, int) {
  is >> lambdaLLE_ >> lambdaLQD_ >> lambdaUDD_ >> iunit(vnu_,GeV)
     >> LLEVertex_ >> LQDVertex_ >> UDDVertex_
     >> HiggsAMix_ >> HiggsPMix_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<RPV,MSSM>
describeHerwigRPV("Herwig::RPV", "HwSusy.so HwRPV.so");

void RPV::Init() {

  static ClassDocumentation<RPV> documentation
    ("The RPV class is the base class for the implementation of the"
     " R-parity violating MSSM.");

  static Reference<RPV,AbstractFFSVertex> interfaceLLEVertex
    ("Vertex/LLE",
     "The vertex for the trillinear LLE interaction",
     &RPV::LLEVertex_, false, false, true, false, false);

  static Reference<RPV,AbstractFFSVertex> interfaceLQDVertex
    ("Vertex/LQD",
     "The vertex for the trillinear LQD interaction",
     &RPV::LQDVertex_, false, false, true, false, false);

  static Reference<RPV,AbstractFFSVertex> interfaceUDDVertex
    ("Vertex/UDD",
     "The vertex for the trillinear UDD interaction",
     &RPV::UDDVertex_, false, false, true, false, false);

}

void RPV::extractParameters(bool checkmodel) {
  MSSM::extractParameters(false);
  if(checkmodel) {
    map<string,ParamMap>::const_iterator pit;
    pit = parameters().find("modsel");
    if(pit == parameters().end()) return;
    ParamMap::const_iterator it;
    // nmssm or mssm
    it = pit->second.find(3);
    int inmssm = (it != pit->second.end()) ? int(it->second) : 0;
    if(inmssm != 0) 
      throw Exception() << "R-parity violating MSSM model"
			<< " used but NMSSM read in." << Exception::runerror; 
    // RPV
    it = pit->second.find(4);
    int irpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(irpv != 1) throw Exception() << "RPV model used but no RPV in input file"
				    << Exception::runerror; 
    // CPV
    it = pit->second.find(5);
    int icpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(icpv != 0) throw Exception() << "RPV model does not support CPV" 
				  << Exception::runerror; 
    // flavour violation
    it = pit->second.find(6);
    int ifv = (it != pit->second.end()) ? int(it->second) : 0;
    if(ifv != 0) throw Exception() << "RPV model does not support "
				 << "flavour violation"
				 << Exception::runerror;
  }
  // get the RPV parameters
  // lambda
  map<string,ParamMap>::const_iterator pit;
  pit=parameters().find("rvlamlle");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaLLE_[i][j][k] = it->second;
    }
  }
  // lambda'
  pit=parameters().find("rvlamlqd");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaLQD_[i][j][k] = it->second;
    }
  }
  // lambda''
  pit=parameters().find("rvlamudd");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaUDD_[i][j][k] = it->second;
    }
  }
  // sneutrino vevs
  pit=parameters().find("rvsnvev");
  vnu_.resize(3);
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first>0) {
	assert(it->first>=1&&it->first<=3);
	vnu_[it->first-1] = it->second*GeV;
      }
    }
  }
}

void RPV::createMixingMatrices() {
  map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    cerr << "testing in mixing loop " << name << "\n";
    // pseudo-scalar higgs mixing
    if (name == "rvamix") {
      createMixingMatrix(HiggsAMix_,name,it->second.second,it->second.first);
    }
    else if (name == "rvlmix") {
      createMixingMatrix(HiggsPMix_,name,it->second.second,it->second.first);
    }
  }
  // base class for neutralinos and charginos
  MSSM::createMixingMatrices();
}

void RPV::doinit() {
  MSSM::doinit();
  addVertex(LLEVertex_);
  addVertex(LQDVertex_);
  addVertex(UDDVertex_);
}
