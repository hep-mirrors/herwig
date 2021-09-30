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
     << ounit(O1_S_,GeV2*GeV)
     << ounit(O1_P_,GeV2*GeV2*GeV) << ounit(O1_D_,GeV2*GeV2*GeV2*GeV);
}

void OniumParameters::persistentInput(PersistentIStream & is, int) {
  is >> singletFromWaveFunction_ >> iunit(R02_,GeV2*GeV)
     >> iunit(Rp02_,GeV2*GeV2*GeV) >> iunit(Rpp02_,GeV2*GeV2*GeV2*GeV)
     >> iunit(O1_S_,GeV2*GeV)
     >> iunit(O1_P_,GeV2*GeV2*GeV) >> iunit(O1_D_,GeV2*GeV2*GeV2*GeV);
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

void OniumParameters::doinit() {
  Interfaced::doinit();
  if ( ! singletFromWaveFunction_ )
    throw Exception() << "Only calculate of singlet elements ffrom wavefunction currently supported\n";
  for(unsigned int type=0;type<R02_.size();++type) {
    if(O1_S_[type].size()<R02_[type].size()) O1_S_[type].resize(R02_[type].size(),{ZERO,ZERO});
    for(unsigned int n=0;n<R02_[type].size();++n) {
      O1_S_[type][n][0] = 4.5/Constants::pi*R02_[type][n];
      O1_S_[type][n][1] = O1_S_[type][n][0];
    }
  }
}
