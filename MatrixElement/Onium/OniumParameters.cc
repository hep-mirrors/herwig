// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OniumParameters class.
//

#include "OniumParameters.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/StringUtils.h"

using namespace Herwig;

IBPtr OniumParameters::clone() const {
  return new_ptr(*this);
}

IBPtr OniumParameters::fullclone() const {
  return new_ptr(*this);
}

void OniumParameters::persistentOutput(PersistentOStream & os) const {
  os << singletFromWaveFunction_ << ounit(R02_,GeV2*GeV)
     << ounit(Rp02_,GeV2*GeV2*GeV) << ounit(Rpp02_,GeV2*GeV2*GeV2*GeV)
     << ounit(O1_S_prod_,GeV2*GeV)
     << ounit(O1_P_prod_,GeV2*GeV2*GeV) << ounit(O1_D_prod_,GeV2*GeV2*GeV2*GeV)
     << ounit(O1_S_dec_,GeV2*GeV)
     << ounit(O1_P_dec_,GeV2*GeV2*GeV) << ounit(O1_D_dec_,GeV2*GeV2*GeV2*GeV)
     << ounit(O8_S_prod_,GeV2*GeV)
     << singletTripletMixing_;
}

void OniumParameters::persistentInput(PersistentIStream & is, int) {
  is >> singletFromWaveFunction_ >> iunit(R02_,GeV2*GeV)
     >> iunit(Rp02_,GeV2*GeV2*GeV) >> iunit(Rpp02_,GeV2*GeV2*GeV2*GeV)
     >> iunit(O1_S_prod_,GeV2*GeV)
     >> iunit(O1_P_prod_,GeV2*GeV2*GeV) >> iunit(O1_D_prod_,GeV2*GeV2*GeV2*GeV)
     >> iunit(O1_S_dec_,GeV2*GeV)
     >> iunit(O1_P_dec_,GeV2*GeV2*GeV) >> iunit(O1_D_dec_,GeV2*GeV2*GeV2*GeV)
     >> iunit(O8_S_prod_,GeV2*GeV)
     >> singletTripletMixing_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<OniumParameters,Interfaced>
describeHerwigOniumParameters("Herwig::OniumParameters", "HwOniumParameters.so");

void OniumParameters::Init() {

  static ClassDocumentation<OniumParameters> documentation
    ("The OniumParameters class stores the parameters for quarkonium production");

  static Command<OniumParameters> interfaceSetWaveFunction
    ("SetWaveFunction",
     "Set the value of the wavefunction, or its deriviatives at the origin",
     &OniumParameters::setWaveFunction, false);

  static Command<OniumParameters> interfaceSetSingletTripletMixing
    ("SetSingletTripletMixing",
     "Set the value of the singlet/triplet mixing for B_c states",
     &OniumParameters::setSingletTripletMixing, false);
  
  static Command<OniumParameters> interfaceSetOctetProductionMatrixElement
    ("SetOctetProductionMatrixElement",
     "Set the matrix element for the production of an octet state",
     &OniumParameters::setOctetProductionMatrixElement, false);

}

string OniumParameters::setWaveFunction(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg         = StringUtils::cdr(arg);
  // get the type of the wavefunction
  OniumState type;
  if(stype=="ccbar")
    type=ccbar;
  else if(stype=="bbbar")
    type=bbbar;
  else if(stype=="bcbar")
    type=bcbar;
  else if(stype=="cc")
    type=cc;
  else if(stype=="bb")
    type=bb;
  else if(stype=="bc")
    type=bc;
  else
    return "Error: Invalid string for wave function type " + stype;
  // get the principal quantum number and orbital AM quantum numbers
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  if(stype.size()!=2)
    return "Error: Invalid string for wave function level " + stype;
  unsigned int n = stoi(stype);
  // get the value
  double value = stof(arg);
  // now get the l value and set the wavefunction
  if(stype[1]=='S') {
    if(n>R02_[type].size()) R02_[type].resize(n,ZERO);
    R02_[type][n-1] = value*GeV2*GeV;
  }
  else if(stype[1]=='P') {
    if(n>Rp02_[type].size()) Rp02_[type].resize(n,ZERO);
    Rp02_[type][n-1] = value*GeV2*GeV2*GeV;
  }
  else if(stype[1]=='D') {
    if(n>Rpp02_[type].size()) Rpp02_[type].resize(n,ZERO);
    Rpp02_[type][n-1] = value*GeV2*GeV2*GeV2*GeV;
  }
  else
    return "Error: Invalid string for wave function l value " + stype;
  // success
  return "";
}

  
string OniumParameters::setSingletTripletMixing(string arg) {
  // get the principal quantum number and orbital AM quantum numbers
  string stype = StringUtils::car(arg);
  arg         = StringUtils::cdr(arg);
  if(stype.size()!=2)
    return "Error: Invalid string for mixing level " + stype;
  unsigned int n = stoi(stype);
  if(n>singletTripletMixing_.size()) {
    singletTripletMixing_.resize(n,vector<double>({0.,0.}));
  }
  // get the value
  double value = stof(arg);
  // now get the l value and set the wavefunction
  if(stype[1]=='S')
    return "No mixing for s-wave states";
  else if(stype[1]=='P')
    singletTripletMixing_[n-1][0] = value/180.*Constants::pi;
  else if(stype[1]=='D')
    singletTripletMixing_[n-1][1] = value/180.*Constants::pi;
  else
    return "Error: Invalid string for mixing l value " + stype;
  // success
  return "";
}

void OniumParameters::doinit() {
  Interfaced::doinit();
  if ( ! singletFromWaveFunction_ )
    throw Exception() << "Only calculate of singlet elements ffrom wavefunction currently supported\n";
  // s wave
  for(unsigned int type=0;type<R02_.size();++type) {
    if(O1_S_prod_[type].size()<R02_[type].size()) O1_S_prod_[type].resize(R02_[type].size(),{ZERO,ZERO});
    if(O1_S_dec_ [type].size()<R02_[type].size()) O1_S_dec_ [type].resize(R02_[type].size(),{ZERO,ZERO});
    for(unsigned int n=0;n<R02_[type].size();++n) {
      Energy3 Odec = 1.5/Constants::pi*R02_[type][n];
      O1_S_prod_[type][n][0] =    Odec;
      O1_S_prod_[type][n][1] = 3.*Odec;
      O1_S_dec_ [type][n][0] =    Odec;
      O1_S_dec_ [type][n][1] =    Odec;
    }
  }
  // p wave
  for(unsigned int type=0;type<Rp02_.size();++type) {
    if(O1_P_prod_[type].size()<Rp02_[type].size()) O1_P_prod_[type].resize(Rp02_[type].size(),{ZERO,ZERO,ZERO,ZERO});
    if(O1_P_dec_ [type].size()<Rp02_[type].size()) O1_P_dec_ [type].resize(Rp02_[type].size(),{ZERO,ZERO,ZERO,ZERO});
    for(unsigned int n=0;n<Rp02_[type].size();++n) {
      Energy5 Odec = 4.5/Constants::pi*Rp02_[type][n];
      O1_P_prod_[type][n][0] = 3.*Odec;
      O1_P_prod_[type][n][1] =    Odec;
      O1_P_prod_[type][n][2] = 3.*Odec;
      O1_P_prod_[type][n][3] = 5.*Odec;
      O1_P_dec_ [type][n][0] =    Odec;
      O1_P_dec_ [type][n][1] =    Odec;
      O1_P_dec_ [type][n][2] =    Odec;
      O1_P_dec_ [type][n][3] =    Odec;
    }
  }
  // d wave
  for(unsigned int type=0;type<Rpp02_.size();++type) {
    if(O1_D_prod_[type].size()<Rpp02_[type].size()) O1_D_prod_[type].resize(Rpp02_[type].size(),{ZERO,ZERO,ZERO,ZERO});
    if(O1_D_dec_ [type].size()<Rpp02_[type].size()) O1_D_dec_ [type].resize(Rpp02_[type].size(),{ZERO,ZERO,ZERO,ZERO});
    for(unsigned int n=0;n<Rpp02_[type].size();++n) {
      Energy7 Odec = 3.75*3./Constants::pi*Rpp02_[type][n];
      O1_D_prod_[type][n][0] = 5.*Odec;
      O1_D_prod_[type][n][1] = 3.*Odec;
      O1_D_prod_[type][n][2] = 5.*Odec;
      O1_D_prod_[type][n][3] = 7.*Odec;
      O1_D_dec_ [type][n][0] =    Odec;
      O1_D_dec_ [type][n][1] =    Odec;
      O1_D_dec_ [type][n][2] =    Odec;
      O1_D_dec_ [type][n][3] =    Odec;
    }
  }
}

string OniumParameters::setOctetProductionMatrixElement(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg         = StringUtils::cdr(arg);
  // get the type of the wavefunction
  OniumState type;
  if(stype=="ccbar")
    type=ccbar;
  else if(stype=="bbbar")
    type=bbbar;
  else
    return "Error: Invalid string for wave function type " + stype;
  // get the orbital AM and spin quantum numbers
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  if(stype.size()!=3)
    return "Error: Invalid string for spin and L quantum numbers " + stype;
  unsigned int s = (stoi(stype)-1)/2;
  unsigned int L=0;
  if     (stype[1]=='S') L=0;
  else if(stype[1]=='P') L=1;
  else if(stype[1]=='D') L=2;
  else 
    return "Error: Invalid string for spin and L quantum numbers " + stype;
  // get the pdg code and value
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  // get the pdg code
  long pid = stoi(stype);
  // get the value
  double value = stof(arg);
  // set the value
  if(L==0) {
    O8_S_prod_[type][s][pid] = value*GeV*GeV2;
  }
  else
    return "Error: Invalid string for wave function l value " + stype;
  // else if(stype[1]=='P') {
  //   if(n>Rp02_[type].size()) Rp02_[type].resize(n,ZERO);
  //   Rp02_[type][n-1] = value*GeV2*GeV2*GeV;
  // }
  // else if(stype[1]=='D') {
  //   if(n>Rpp02_[type].size()) Rpp02_[type].resize(n,ZERO);
  //   Rpp02_[type][n-1] = value*GeV2*GeV2*GeV2*GeV;
  // }
  // success
  return "";
}
