// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LittleHiggsModel class.
//

#include "LittleHiggsModel.h"
#include "ThePEG/Repository/BaseRepository.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LittleHiggsModel::persistentOutput(PersistentOStream & os) const {
  os << _g << _gp << _cott << _tantp << ounit(_v,GeV) << _lamratio
     << ounit(_mH,GeV) << _vacratio << ounit(_f,GeV) 
     << _lambda1 << _lambda2 << _s << _c << _sp << _cp;
}

void LittleHiggsModel::persistentInput(PersistentIStream & is, int) {
  is >> _g >> _gp >> _cott >> _tantp >> iunit(_v,GeV)  >> _lamratio
     >> iunit(_mH,GeV) >> _vacratio >> iunit(_f,GeV) 
     >> _lambda1 >> _lambda2 >> _s >> _c >> _sp >> _cp;
}

ClassDescription<LittleHiggsModel> LittleHiggsModel::initLittleHiggsModel;
// Definition of the static class description member.

void LittleHiggsModel::Init() {

  static ClassDocumentation<LittleHiggsModel> documentation
    ("The LittleHiggsModel class");

  static Parameter<LittleHiggsModel,double> interfacegCoupling
    ("gCoupling",
     "The g coupling",
     &LittleHiggsModel::_g,  sqrt(0.43), 0.0, 1.,
     false, false, Interface::limited);

  static Parameter<LittleHiggsModel,double> interfacegPrimeCoupling
    ("gPrimeCoupling",
     "The g' coupling",
     &LittleHiggsModel::_gp, sqrt(0.12), 0.0, 1.,
     false, false, Interface::limited);

  static Parameter<LittleHiggsModel,double> interfaceCotTheta
    ("CotTheta",
     "The cotangent of the theta mixing angle",
     &LittleHiggsModel::_cott, 1.0, 0.1, 10.,
     false, false, Interface::limited);

  static Parameter<LittleHiggsModel,double> interfaceTanThetaPrime
    ("TanThetaPrime",
     "The tangent of the theta' mixing angle",
     &LittleHiggsModel::_tantp, 1.0, 0.1, 10.0,
     false, false, Interface::limited);

  static Parameter<LittleHiggsModel,Energy> interfacef
    ("f",
     "The scale of the non-linear sigma-model",
     &LittleHiggsModel::_f, TeV, 3.*TeV, 0.0*TeV, 10.0*TeV,
     true, false, Interface::limited);

  static Parameter<LittleHiggsModel,double> interfaceLambdaRatio
    ("LambdaRatio",
     "The ratio lambda_2/lambda_1 of the top Yukawa couplings.",
     &LittleHiggsModel::_lamratio, 1.0, 0.01, 100.,
     false, false, Interface::limited);

  static Parameter<LittleHiggsModel,double> interfaceVEVRatio
    ("VEVRatio",
     "The ratio of the vacuum expection values v'/v",
     &LittleHiggsModel::_vacratio, 0.05, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<LittleHiggsModel,Energy> interfacemH
    ("mH",
     "The mass of the lightest Higgs",
     &LittleHiggsModel::_mH, GeV, 120.0*GeV, 100.0*GeV, 1000.0*GeV,
     false, false, Interface::limited);

}

void LittleHiggsModel::doinit() throw(InitException) {
  // compute the parameters of the model
  // W and Z masses
  Energy mw(getParticleData(ParticleID::Wplus)->mass()),
    mz(getParticleData(ParticleID::Z0)->mass());
  // vev
  _v = 2.*mw/_g;
  // cos and sin of mixing angles
  double theta (atan(1./_cott));
  _c = cos(theta );
  _s = sin(theta );
  double thetap(atan(_tantp  ));
  _cp = cos(thetap);
  _sp = sin(thetap);
  double sw2(sin2ThetaW()),cw2(1.-sin2ThetaW());
  // masses of the neutral gauge bosons
  double xH(2.5*_g*_gp*_s*_c*_sp*_cp*(sqr(_c*_sp)+sqr(_s*_cp))/
	    (5.*sqr(_g*_sp*_cp)-sqr(_gp*_s*_c)));
  double vf(sqr(_v/_f));
  Energy2 MZL2 = sqr(mz)*(1.-vf*(1./6.+0.25*sqr(sqr(_c)-sqr(_s))
				 +1.25*sqr(sqr(_cp)-sqr(_sp)))+8.*sqr(_vacratio));
  Energy2 MAH2 = sqr(mz)*sw2*(0.2/sqr(_sp*_cp)/vf-1.+0.25*xH*cw2/sqr(_s*_c)/sw2);
  Energy2 MZH2 = sqr(mw)*(1./sqr(_s*_c)/vf-1.-xH*sw2/sqr(_sp*_cp)/cw2);
  // mass of the charged gauge bosons
  Energy2 MWL2 = sqr(mw)*(1.-vf*(1./6.+0.25*sqr(sqr(_c)-sqr(_s))));
  Energy2 MWH2 = sqr(mw)*(1./sqr(_s*_c)/vf-1.);
  // top and heavy top yukawas
  Energy mt = getParticleData(ParticleID::t)->mass();
  Energy MT = _f/_v*(1.+sqr(_lamratio))/_lamratio*mt;
  _lambda1 = MT/sqrt(1.+sqr(_lamratio))/_f;
  _lambda2 = _lamratio*_lambda1;
  // masses of the Higgs bosons
  double r = 8.*_f/_v*_vacratio;
  double lamh = 2.*sqr(_mH/_v)/(1./r-0.25*r)/r;
  Energy2 MPhi2 = lamh*sqr(_f);
  // set the masses of the new particles
  resetMass( 32,sqrt(MAH2));
  resetMass( 33,sqrt(MZH2));
  resetMass( 34,sqrt(MWH2));
  resetMass(-34,sqrt(MWH2));
  resetMass(  8,MT);
  resetMass( -8,MT);
  resetMass( 35,sqrt(MPhi2));
  resetMass( 36,sqrt(MPhi2));
  resetMass( 37,sqrt(MPhi2));
  resetMass(-37,sqrt(MPhi2));
  resetMass( 38,sqrt(MPhi2));
  resetMass(-38,sqrt(MPhi2));
  // stuff from the base class
  StandardModel::doinit();
}

void LittleHiggsModel::resetMass(long id, Energy mass) {
  tPDPtr part = getParticleData(id);
  if(!part) return;
  part->init();
  const InterfaceBase * ifb = BaseRepository::FindInterface(part, "NominalMass");
  ostringstream os;
  os << abs(mass/GeV);
  ifb->exec(*part, "set", os.str());
}
