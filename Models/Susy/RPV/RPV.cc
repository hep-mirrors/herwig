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
     << LLEVertex_ << LQDVertex_ << UDDVertex_;
}

void RPV::persistentInput(PersistentIStream & is, int) {
  is >> lambdaLLE_ >> lambdaLQD_ >> lambdaUDD_ >> iunit(vnu_,GeV)
     >> LLEVertex_ >> LQDVertex_ >> UDDVertex_;
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
    // pseudo-scalar higgs mixing
    if (name == "rvamix") {
      MixingMatrixPtr temp;
      createMixingMatrix(temp,name,it->second.second,it->second.first);
      CPoddHiggsMix(temp); 
      
    }
    else if (name == "rvlmix") {
      MixingMatrixPtr temp;
      createMixingMatrix(temp,name,it->second.second,it->second.first);
      ChargedHiggsMix(temp);
    }
  }
  // base class for neutralinos and charginos
  MSSM::createMixingMatrices();
  // now adjust the mixing matrices to have our structure
  // first bloodly SPHENO as it doesn't obey the SLHA  
  map<string,StringMap>::const_iterator sit = info().find("spinfo");
  string program;
  if(sit!=info().end()) {
    StringMap::const_iterator pit = sit->second.find(1);
    if(pit!=sit->second.end()) program = pit->second;
  }
  if(program=="SPheno") {
    map<string,ParamMap>::const_iterator fit=parameters().find("mass");
    if(fit==parameters().end()) 
      throw Exception() << "BLOCK MASS not found in input file"
			<< " can't set masses of SUSY particles"
			<< Exception::runerror;
    // adjust the charged scalars
    map<double,long> massMap;
    massMap[findValue(fit,     37,"mass",     "37")] =      37;
    massMap[findValue(fit,1000011,"mass","1000011")] = 1000011;
    massMap[findValue(fit,1000013,"mass","1000013")] = 1000013;
    massMap[findValue(fit,1000015,"mass","1000015")] = 1000015;
    massMap[findValue(fit,2000011,"mass","2000011")] = 2000011;
    massMap[findValue(fit,2000013,"mass","2000013")] = 2000013;
    massMap[findValue(fit,2000015,"mass","2000015")] = 2000015;
    vector<int> move(1,7);
    for(map<double,long>::iterator mit=massMap.begin();mit!=massMap.end();++mit) {
      if     (mit->second==     37) move.push_back(0);
      else if(mit->second==1000011) move.push_back(1);
      else if(mit->second==1000013) move.push_back(2);
      else if(mit->second==1000015) move.push_back(3);
      else if(mit->second==2000011) move.push_back(4);
      else if(mit->second==2000013) move.push_back(5);
      else if(mit->second==2000015) move.push_back(6);
    }
    CMatrix oldMat = ChargedHiggsMix()->getMatrix();
    CMatrix newMat(8,vector<Complex>(8,0.));
    for(unsigned int ix=0;ix<8;++ix) {
      for(unsigned int iy=0;iy<8;++iy)
	newMat[move[ix]][iy] = oldMat[ix][iy];
    }
    ChargedHiggsMix(new_ptr(MixingMatrix(newMat,ChargedHiggsMix()->getIds())));
    // adjust the pseudoscalars
    massMap.clear();
    massMap[findValue(fit,     36,"mass",     "36")] =      36;
    // extract the pseudoscalar masses and change the ids for the pseudoscalars
    // to those from the SLHA if needed
    if(fit->second.find(2000012)!=fit->second.end()) {
      massMap[findValue(fit,2000012,"mass","2000012")] = 1000017;
      idMap().insert(make_pair(2000012,1000017));
    }
    else {
      massMap[findValue(fit,1000017,"mass","1000017")] = 1000017;
    }
    if(fit->second.find(2000014)!=fit->second.end()) {
      massMap[findValue(fit,2000014,"mass","2000014")] = 1000018;
      idMap().insert(make_pair(2000014,1000018));
    }
    else {
      massMap[findValue(fit,1000018,"mass","1000018")] = 1000018;
    }
    if(fit->second.find(2000016)!=fit->second.end()) {
      massMap[findValue(fit,2000016,"mass","2000016")] = 1000019;
      idMap().insert(make_pair(2000016,1000019));
    }
    else {
      massMap[findValue(fit,1000019,"mass","1000019")] = 1000019;
    }
    move.clear(); move.push_back(4);
    for(map<double,long>::iterator mit=massMap.begin();mit!=massMap.end();++mit) {
      if     (mit->second==     36) move.push_back(0);
      else if(mit->second==1000017) move.push_back(1);
      else if(mit->second==1000018) move.push_back(2);
      else if(mit->second==1000019) move.push_back(3);
    }
    oldMat = CPoddHiggsMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[move[ix]][iy] = oldMat[ix][iy];
    }
    CPoddHiggsMix(new_ptr(MixingMatrix(newMat,CPoddHiggsMix()->getIds())));
    // adjust the neutral scalars
    massMap.clear();
    massMap[findValue(fit,     25,"mass",     "25")] =      25;
    massMap[findValue(fit,     35,"mass",     "35")] =      35;
    massMap[findValue(fit,1000012,"mass","1000012")] = 1000012;
    massMap[findValue(fit,1000014,"mass","1000014")] = 1000014;
    massMap[findValue(fit,1000016,"mass","1000016")] = 1000016;
    move.clear();
    for(map<double,long>::iterator mit=massMap.begin();mit!=massMap.end();++mit) {
      if     (mit->second==     25) move.push_back(0);
      else if(mit->second==     35) move.push_back(1);
      else if(mit->second==1000012) move.push_back(2);
      else if(mit->second==1000014) move.push_back(3);
      else if(mit->second==1000016) move.push_back(4);
    }
    oldMat = CPevenHiggsMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[move[ix]][iy] = oldMat[ix][iy];
    }
    CPevenHiggsMix(new_ptr(MixingMatrix(newMat,CPevenHiggsMix()->getIds())));
    // neutralino mixing
    move.resize(7);
    move[0] = 3; move[1] = 4; move[2] = 5; move[3] = 6;
    move[4] = 0; move[5] = 1; move[6] = 2;
    oldMat = neutralinoMix()->getMatrix();
    newMat = CMatrix(7,vector<Complex>(7,0.));
    for(unsigned int ix=0;ix<7;++ix) {
      for(unsigned int iy=0;iy<7;++iy)
	newMat[ix][move[iy]] = oldMat[ix][iy];
    }
    neutralinoMix(new_ptr(MixingMatrix(newMat,neutralinoMix()->getIds())));
    // chargino mixing
    move.resize(5);
    move[0] = 3; move[1] = 4;
    move[2] = 0; move[3] = 1; move[4] = 2;
    oldMat = charginoUMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[ix][move[iy]] = oldMat[ix][iy];
    }
    charginoUMix(new_ptr(MixingMatrix(newMat,charginoUMix()->getIds())));
    oldMat = charginoVMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[ix][move[iy]] = oldMat[ix][iy];
    }
    charginoVMix(new_ptr(MixingMatrix(newMat,charginoVMix()->getIds())));
  }
  // we don't want neutrinos first then neutralinos so swap them
  // neutralinos first then neutrinos
  vector<int> move(7);
  move[0] = 4; move[1] = 5; move[2] = 6;
  move[3] = 0; move[4] = 1; move[5] = 2; move[6] = 3;
  CMatrix oldMat = neutralinoMix()->getMatrix();
  CMatrix newMat(7,vector<Complex>(7,0.));
  for(unsigned int ix=0;ix<7;++ix) {
    for(unsigned int iy=0;iy<7;++iy)
      newMat[move[ix]][move[iy]] = oldMat[ix][iy];
  }
  neutralinoMix(new_ptr(MixingMatrix(newMat,neutralinoMix()->getIds())));
  // charginos the same
  move.resize(5);
  move[0] = 2; move[1] = 3; move[2] = 4;
  move[3] = 0; move[4] = 1;
  oldMat = charginoUMix()->getMatrix();
  newMat = CMatrix(5,vector<Complex>(5,0.));
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<5;++iy)
      newMat[move[ix]][move[iy]] = oldMat[ix][iy];
  }
  charginoUMix(new_ptr(MixingMatrix(newMat,charginoUMix()->getIds())));
  oldMat = charginoVMix()->getMatrix();
  newMat = CMatrix(5,vector<Complex>(5,0.));
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<5;++iy)
      newMat[move[ix]][move[iy]] = oldMat[ix][iy];
  }
  charginoVMix(new_ptr(MixingMatrix(newMat,charginoVMix()->getIds())));
}

void RPV::doinit() {
  MSSM::doinit();
  addVertex(LLEVertex_);
  addVertex(LQDVertex_);
  addVertex(UDDVertex_);
}
