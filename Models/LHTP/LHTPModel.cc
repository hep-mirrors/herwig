// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPModel class.
//

#include "LHTPModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/BaseRepository.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "gsl/gsl_multiroots.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include <algorithm>

using namespace Herwig;
using namespace ThePEG;

// equations for top parameters for GSL
namespace {

// struct to provide the model parameters to the functions
struct tparams {
  Energy v;
  Energy f;
  Energy mt;
  double tan2a;
};

// equations defining tan 2alpha and mt expressed in form f(lambda1,lambda2)=0
// to be solved to give lambda1 and lambda2 by gsl
int top_equation(const gsl_vector * x, void *params, gsl_vector *f ) {
  // yukawa and check top mass
  const double lam1 = gsl_vector_get(x,0);
  const double lam2 = gsl_vector_get(x,1);
  Energy fs = ((struct tparams *) params)->f;
  Energy v = ((struct tparams *) params)->v;
  double sv = sin(sqrt(2.)*v/fs);
  double cv = cos(sqrt(2.)*v/fs);
  Energy mt    = ((struct tparams *) params)->mt;
  double tan2a = ((struct tparams *) params)->tan2a;
  double f1 = 4.*lam1*lam2*(1.+cv)/(4.*sqr(lam2)-sqr(lam1)*(2.*sqr(sv)+sqr(1.+cv)))
    -tan2a;
  double delta = 0.5*(sqr(lam2)+0.5*sqr(lam1)*(sqr(sv)+0.5*sqr(1.+cv)));
  double f2 = sqr(fs/mt)*delta*(1.-sqrt(1.-0.5*sqr(lam1*lam2*sv/delta)))-1.;
  if(lam1*lam2<0.) f1+=1e10;
  if(lam1*lam2<0.) f2+=1e10;
  gsl_vector_set(f,0,f1);
  gsl_vector_set(f,1,f2);
  return GSL_SUCCESS;
}

}

LHTPModel::LHTPModel()
  : _f(0.5*TeV), _salpha(sqrt(0.500)), _calpha(0.), _sbeta(0.), _cbeta(0.),
    _kappa(1.), _mh(120.*GeV), _v(246.*GeV),_g(sqrt(0.43)), _gp(sqrt(0.12)) {}

IBPtr LHTPModel::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPModel::fullclone() const {
  return new_ptr(*this);
}

void LHTPModel::persistentOutput(PersistentOStream & os) const {
  os << ounit(_f,TeV) << _salpha << _calpha << _sbeta << _cbeta
     << _kappa << ounit(_v,GeV) << _g << _gp << _sthetaH << _cthetaH;
}

void LHTPModel::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_f,TeV) >> _salpha >> _calpha >> _sbeta >> _cbeta
     >> _kappa >> iunit(_v,GeV) >> _g >> _gp >> _sthetaH >> _cthetaH;
}

ClassDescription<LHTPModel> 
LHTPModel::initLHTPModel;
// Definition of the static class description member.

void LHTPModel::Init() {

  static ClassDocumentation<LHTPModel> documentation
    ("The LHTPModel class implements the Little Higgs model"
     " with T-parity");

  static Parameter<LHTPModel,Energy> interfacef
    ("f",
     "The scale of the non-linear sigma-model",
     &LHTPModel::_f, TeV, 1.*TeV, 0.0*TeV, 10.0*TeV,
     true, false, Interface::limited);

  static Parameter<LHTPModel,double> interfaceSinAlpha
    ("SinAlpha",
     "The parameter controlling the mixing in the top quark sector of the model",
     &LHTPModel::_salpha, sqrt(0.5), 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<LHTPModel,double> interfaceKappa
    ("Kappa",
     "The parameter controlling the masses of the T-odd fermions",
     &LHTPModel::_kappa, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<LHTPModel,Energy> interfaceHiggsMass
    ("HiggsMass",
     "The mass of the lightest Higgs boson",
     &LHTPModel::_mh, GeV, 120.0*GeV, 100.0*GeV, 1000.0*GeV,
     false, false, Interface::limited);
}

void LHTPModel::resetMass(long id, Energy mass) {
  tPDPtr part = getParticleData(id);
  if(!part) return;
  //  part->init();
  const InterfaceBase * ifb = BaseRepository::FindInterface(part, "NominalMass");
  ostringstream os;
  os << abs(mass/GeV);
  ifb->exec(*part, "set", os.str());
}

void LHTPModel::doinit() {
  StandardModel::doinit();
  string name = CurrentGenerator::current().filename() +
    string("-BSMModelInfo.out");
  ofstream dummy(name.c_str());
  using Constants::pi;
  // compute the parameters of the model
  // W and Z masses
  Energy mw(getParticleData(ParticleID::Wplus)->mass()),
    mz(getParticleData(ParticleID::Z0)->mass());
  // couplings g and g'
  double ee = sqrt(4.*pi*alphaEM(sqr(mz)));
  double sw2(sin2ThetaW()),cw2(1.-sin2ThetaW());
  double sw(sqrt(sw2)),cw(sqrt(cw2));
  _g  = ee/sw;
  _gp = ee/cw;
  // vev
  _v = 2.*mw/_g;
  double vf(sqr(_v/_f));
  // calculate masses of the new particles from input
  // and SM parameters
  // masses of the new gauge bosons (MWH = MZH)
  Energy MAH = _gp*_f*sqrt(0.2)*(1.-0.625*vf);
  Energy MZH = _g *_f*          (1.-0.125*vf);
  // mixings
  _sthetaH = 1.25*_g*_gp/(5.*sqr(_g)-sqr(_gp))*vf;
  _cthetaH = sqrt(1.-sqr(_sthetaH));
  // masses of the new top quarks
  Energy MTp,MTm;
  topMixing(MTp,MTm);
  // masses of the T-odd fermions
  Energy Mdm = sqrt(2.)*_kappa*_f;
  Energy Mum = sqrt(2.)*_kappa*_f*(1.-0.125*vf);
  Energy Mlm = sqrt(2.)*_kappa*_f;
  Energy Mnm = sqrt(2.)*_kappa*_f*(1.-0.125*vf);
  // masses of the triplet higgs
  Energy MPhi = sqrt(2.)*_mh*_f/_v;
  // set the masses of the new particles
  // new gauge bosons
  resetMass( 32      , MAH  );
  resetMass( 33      , MZH  );
  resetMass( 34      , MZH  );
  resetMass(-34      , MZH  );
  // masses of the new top quarks
  resetMass(  8      , MTp  );
  resetMass( -8      , MTp  );
  resetMass(  4000008, MTm  );
  resetMass( -4000008, MTm  );
  //  masses of the Higgs bosons
  resetMass( 25      , _mh  );
  resetMass( 35      , MPhi );
  resetMass( 36      , MPhi );
  resetMass( 37      , MPhi );
  resetMass(-37      , MPhi );
  resetMass( 38      , MPhi );
  resetMass(-38      , MPhi );
  // masses of the T-odd quarks
  resetMass( 4000001, Mdm  );
  resetMass(-4000001, Mdm  );
  resetMass( 4000002, Mum  );
  resetMass(-4000002, Mum  );
  resetMass( 4000003, Mdm  );
  resetMass(-4000003, Mdm  );
  resetMass( 4000004, Mum  );
  resetMass(-4000004, Mum  );
  resetMass( 4000005, Mdm  );
  resetMass(-4000005, Mdm  );
  resetMass( 4000006, Mum  );
  resetMass(-4000006, Mum  );
  // masses of the T-odd leptons
  resetMass( 4000011, Mlm  );
  resetMass(-4000011, Mlm  );
  resetMass( 4000012, Mnm  );
  resetMass(-4000012, Mnm  );
  resetMass( 4000013, Mlm  );
  resetMass(-4000013, Mlm  );
  resetMass( 4000014, Mnm  );
  resetMass(-4000014, Mnm  );
  resetMass( 4000015, Mlm  );
  resetMass(-4000015, Mlm  );
  resetMass( 4000016, Mnm  );
  resetMass(-4000016, Mnm  );
}

void LHTPModel::topMixing(Energy & MTp, Energy & MTm) {
  Energy mt = getParticleData(ParticleID::t)->mass();
  _calpha = sqrt(1.-sqr(_salpha));
  double sv(sin(sqrt(2.)*_v/_f)),cv(cos(sqrt(2.)*_v/_f));
  // first guess for Yukawa's based on leading order in v/f expansion
  double lambda1(mt/_v/_calpha), lambda2(mt/_salpha/_v);
  MTp = lambda1/_salpha*_f;
  MTm = lambda1/_salpha*_calpha*_f;
  // special case where denominator of tan 2 alpha eqn is zero
  if(abs(_salpha-sqrt(0.5))<1e-4) {
    double a = 0.25*(2.*sqr(sv)+sqr(1.+cv));
    double b = 0.5*(a+0.5*(sqr(sv)+0.5*sqr(1.+cv)));
    lambda1 = mt/_f*sqrt(1./b/(1.-sqrt(1-0.5*a*sqr(sv/b))));
    lambda2 = sqrt(a)*lambda1;
  }
  // general case using GSL
  else {
    double ca = sqrt(1.-sqr(_salpha));
    double ta = _salpha/ca;
    double tan2a = 2.*ta/(1.-sqr(ta));
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    int status;
    size_t iter=0;
    const size_t n=2;
    struct tparams p = {_v,_f,mt,tan2a};
    gsl_multiroot_function f = {&top_equation, n, &p};
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector_set(x,0,lambda1);
    gsl_vector_set(x,1,lambda2);
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T,2);
    gsl_multiroot_fsolver_set(s, &f,x);
    do {
      iter++;
      status = gsl_multiroot_fsolver_iterate(s);
      if(status) break;
      status = gsl_multiroot_test_residual(s->f,1e-7);
    }
    while (status==GSL_CONTINUE && iter < 1000);
    gsl_multiroot_fsolver_free(s);
    lambda1 = gsl_vector_get(s->x,0);
    lambda2 = gsl_vector_get(s->x,1);
    gsl_vector_free(x);
  }
  // calculate the heavy top masses using full result
  double delta = 0.5*(sqr(lambda2)+0.5*sqr(lambda1)*(sqr(sv)+0.5*sqr(1.+cv)));
  double det = sqrt(1.-0.5*sqr(lambda1*lambda2*sv/delta));
  MTp = sqrt(sqr(_f)*delta*(1.+det));
  MTm = lambda2*_f;
  // beta mixing angle
  double beta = 0.5*atan(2.*sqrt(2.)*sqr(lambda1)*sv*(1.+cv)/
			 (4.*sqr(lambda2)+sqr(1.+cv)*sqr(lambda1)-2.*sqr(lambda1)*sv));
  _sbeta = sin(beta);
  _cbeta = cos(beta);
}
