// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHModel class.
//

#include "LHModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHModel::LHModel() 
  : _cott(1.), _tantp(1.),
    _v(246.*GeV), _lamratio(1.), _mH(120.*GeV), _vacratio(0.05),
    _f(3.*TeV), _lambda1(0.), _lambda2(0.),
    _s(0.), _c(0.), _sp(0.), _cp(0.)
{}

void LHModel::persistentOutput(PersistentOStream & os) const {
  os << _cott << _tantp << ounit(_v,GeV) << _lamratio
     << ounit(_mH,GeV) << _vacratio << ounit(_f,GeV) 
     << _s0 << _c0 << _sP << _cP << _sPlus << _cPlus
     << _lambda1 << _lambda2 << _s << _c << _sp << _cp
     << WHHVertex_;
}

void LHModel::persistentInput(PersistentIStream & is, int) {
  is >> _cott >> _tantp >> iunit(_v,GeV)  >> _lamratio
     >> iunit(_mH,GeV) >> _vacratio >> iunit(_f,GeV)  
     >> _s0 >> _c0 >> _sP >> _cP >> _sPlus >> _cPlus
     >> _lambda1 >> _lambda2 >> _s >> _c >> _sp >> _cp
     >> WHHVertex_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHModel,StandardModel>
describeThePEGLHModel("Herwig::LHModel", "HwLHModel.so");

void LHModel::Init() {

  static ClassDocumentation<LHModel> documentation
    ("The LHModel class");

  static Parameter<LHModel,double> interfaceCotTheta
    ("CotTheta",
     "The cotangent of the theta mixing angle",
     &LHModel::_cott, 1.0, 0.1, 10.,
     false, false, Interface::limited);

  static Parameter<LHModel,double> interfaceTanThetaPrime
    ("TanThetaPrime",
     "The tangent of the theta' mixing angle",
     &LHModel::_tantp, 1.0, 0.1, 10.0,
     false, false, Interface::limited);

  static Parameter<LHModel,Energy> interfacef
    ("f",
     "The scale of the non-linear sigma-model",
     &LHModel::_f, TeV, 3.*TeV, 0.0*TeV, 100.0*TeV,
     true, false, Interface::limited);

  static Parameter<LHModel,double> interfaceLambdaRatio
    ("LambdaRatio",
     "The ratio lambda_1/lambda_2 of the top Yukawa couplings.",
     &LHModel::_lamratio, 1.0, 0.01, 100.,
     false, false, Interface::limited);

  static Parameter<LHModel,double> interfaceVEVRatio
    ("VEVRatio",
     "The ratio of the vacuum expection values v'/v",
     &LHModel::_vacratio, 0.05, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<LHModel,Energy> interfacemH
    ("mH",
     "The mass of the lightest Higgs",
     &LHModel::_mH, GeV, 120.0*GeV, 100.0*GeV, 1000.0*GeV,
     false, false, Interface::limited);

  static Reference<LHModel,AbstractVSSVertex> interfaceVertexWHH
    ("Vertex/WHH",
     "Pointer to the WHH vertex",
     &LHModel::WHHVertex_, false, false, true, false, false);

}

void LHModel::doinit() {
  // stuff from the base class
  StandardModel::doinit();
  // compute the parameters of the model
  // W and Z masses
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  Energy mz(getParticleData(ParticleID::Z0)   ->mass());
  // SM couplings
  double e  = sqrt(4.*Constants::pi*alphaEM(sqr(mz)));
  double sw2(sin2ThetaW()),cw2(1.-sin2ThetaW());
  double g  = e/sqrt(sw2);
  double gp = e/sqrt(cw2);
  // vev
  _v = 2.*mw/g;
  // cos and sin of mixing angles
  double theta (atan(1./_cott));
  _c = cos(theta );
  _s = sin(theta );
  double thetap(atan(_tantp  ));
  _cp = cos(thetap);
  _sp = sin(thetap);
  // xH (Eqn A35)
  double xH = 2.5*g*gp*_s*_c*_sp*_cp*(sqr(_c*_sp)+sqr(_s*_cp))/
    (5.*sqr(g*_sp*_cp)-sqr(gp*_s*_c));
  double vf(sqr(_v/_f));
  // masses of the neutral gauge bosons (Eqn 21,22,A37)
  Energy2 MAH2 = sqr(mz)*sw2*(0.2/sqr(_sp*_cp)/vf-1.+0.25*xH*cw2/sqr(_s*_c)/sw2);
  Energy2 MZH2 = sqr(mw)*(1./sqr(_s*_c)/vf-1.-xH*sw2/sqr(_sp*_cp)/cw2);
  // mass of the heavy charged gauge boson (Eqn. 19/A33) 
  Energy2 MWH2 = sqr(mw)*(1./sqr(_s*_c)/vf-1.);
  // top and heavy top yukawas (from Eqns 26,27)
  Energy mt = getParticleData(ParticleID::t)->mass();
  Energy MT = _f/_v*(1.+sqr(_lamratio))/_lamratio*mt;
  _lambda2 = MT/sqrt(1.+sqr(_lamratio))/_f;
  _lambda1 = _lamratio*_lambda2;
  // masses of the Higgs bosons (Eqns 12 and 13)
  double r = 8.*_f/_v*_vacratio;
  double lamh = 2.*sqr(_mH/_v)/(1.-0.25*sqr(r));
  if(lamh<0.) {
    throw Exception() << "Higgs trilinear coupling negative, reduce f or v'\n"
		      << Exception::runerror;
  }
  Energy2 MPhi2 = lamh*sqr(_f);
  // from Eqn A27
  _sP    = 2.*sqrt(2.)*_vacratio;
  _cP    = 1.-4.*sqr(_vacratio);
  _sPlus = 2.*_vacratio;
  _cPlus = 1.-2.*sqr(_vacratio);
  // from Eqn A28
  _s0 = 2.*sqrt(2.)*_vacratio;
  _c0 = 1.-4.*sqr(_vacratio);
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
}
