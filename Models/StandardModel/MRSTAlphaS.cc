// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MRSTAlphaS class.
//

#include "MRSTAlphaS.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

extern "C" {
  double mrstalphas_(double*,double*,double*,double*,double*,const int*);
}

MRSTAlphaS::MRSTAlphaS() 
  : _lambda(323.*MeV), _mcharm(1.43*GeV), _mbottom(4.3*GeV),
    _flavour(3), _iord(1)
{}

IBPtr MRSTAlphaS::clone() const {
  return new_ptr(*this);
}

IBPtr MRSTAlphaS::fullclone() const {
  return new_ptr(*this);
}

void MRSTAlphaS::persistentOutput(PersistentOStream & os) const {
  os << ounit(_lambda,GeV) << ounit(_mcharm,GeV) << ounit(_mbottom,GeV)
     << _flavour << _iord << ounit(_qsct,GeV2) << ounit(_qsbt,GeV2);
}

void MRSTAlphaS::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_lambda,GeV) >> iunit(_mcharm,GeV) >> iunit(_mbottom,GeV)
     >> _flavour >> _iord >> iunit(_qsct,GeV2) >> iunit(_qsbt,GeV2);
}

void MRSTAlphaS::doinit() {
  AlphaSBase::doinit();
  _qsct = 4.*sqr(_mcharm);
  _qsbt = 4.*sqr(_mbottom);
  for(Energy scale=1.*GeV;scale<201.*GeV;scale+=GeV) {
    cerr << scale/GeV << "\t" << value(sqr(scale),
				       *generator()->standardModel()) << "\n";
  }
}


ClassDescription<MRSTAlphaS> MRSTAlphaS::initMRSTAlphaS;
// Definition of the static class description member.

void MRSTAlphaS::Init() {

  static ClassDocumentation<MRSTAlphaS> documentation
    ("The MRSTAlphaS class implements the alphaS used by MRSt for testing");

  static Parameter<MRSTAlphaS,Energy> interfaceLambda
    ("Lambda",
     "Lambda QCd",
     &MRSTAlphaS::_lambda, GeV, 0.323*GeV, 0.25*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<MRSTAlphaS,Energy> interfaceCharmMass
    ("CharmMass",
     "The mass of the charm quark",
     &MRSTAlphaS::_mcharm, GeV, 1.43*GeV, 1.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<MRSTAlphaS,Energy> interfaceBottomMass
    ("BottomMass",
     "The mass of the bottom quark",
     &MRSTAlphaS::_mbottom, GeV, 4.3*GeV, 1.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<MRSTAlphaS,int> interfaceNumberOfFlavours
    ("NumberOfFlavours",
     "Number of flavours",
     &MRSTAlphaS::_flavour, 3, 1, 5,
     false, false, Interface::limited);

  static Switch<MRSTAlphaS,int> interfaceOrder
    ("Order",
     "Order",
     &MRSTAlphaS::_iord, 1, false, false);
  static SwitchOption interfaceOrderLeadingOrder
    (interfaceOrder,
     "LeadingOrder",
     "1st order",
     0);
  static SwitchOption interfaceOrderNLO
    (interfaceOrder,
     "NLO",
     "2nd order",
     1);
}

vector<Energy> MRSTAlphaS::LambdaQCDs() const {
  return vector<Energy>(6,_lambda);
}

vector<Energy2> MRSTAlphaS::flavourThresholds() const {
  vector<Energy2> thresholds(5,ZERO);
  thresholds[3] = _qsct;
  thresholds[4] = _qsbt;
  return thresholds;
}

double MRSTAlphaS::value(Energy2 scale, const StandardModelBase &) const {
  double t = log(scale/sqr(_lambda));
  double arg1=_lambda/GeV;
  double arg2=_flavour;
  double arg3=_qsct/GeV2;
  double arg4=_qsbt/GeV2;
  double output= mrstalphas_(&t,&arg1,&arg2,&arg4,&arg3,&_iord);
  return output;
}



