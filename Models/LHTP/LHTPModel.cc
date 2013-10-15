// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPModel class.
//

#include "LHTPModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
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
  : f_(0.5*TeV), salpha_(sqrt(0.5)), calpha_(sqrt(0.5)), sbeta_(0.), cbeta_(0.), 
    sL_(0.), cL_(1.), sR_(0.), cR_(0.),
    kappaQuark_(1.), kappaLepton_(1.), mh_(125.*GeV), v_(246.*GeV),
    g_(sqrt(0.43)), gp_(sqrt(0.12)), approximate_(false)
{}

IBPtr LHTPModel::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPModel::fullclone() const {
  return new_ptr(*this);
}

void LHTPModel::persistentOutput(PersistentOStream & os) const {
  os << ounit(f_,TeV) << salpha_ << calpha_ << sbeta_ << cbeta_
     << kappaQuark_ << kappaLepton_ << ounit(v_,GeV) 
     << g_ << gp_ << sthetaH_ << cthetaH_ << approximate_
     << sL_ << cL_ << sR_ << cR_ << WHHVertex_;
}

void LHTPModel::persistentInput(PersistentIStream & is, int) {
  is >> iunit(f_,TeV) >> salpha_ >> calpha_ >> sbeta_ >> cbeta_
     >> kappaQuark_ >> kappaLepton_>> iunit(v_,GeV)
     >> g_ >> gp_ >> sthetaH_ >> cthetaH_ >> approximate_
     >> sL_ >> cL_ >> sR_ >> cR_ >> WHHVertex_;
}


// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPModel,StandardModel>
describeHerwigLHTPModel("Herwig::LHTPModel", "HwLHTPModel.so");

void LHTPModel::Init() {

  static ClassDocumentation<LHTPModel> documentation
    ("The LHTPModel class implements the Little Higgs model"
     " with T-parity");

  static Parameter<LHTPModel,Energy> interfacef
    ("f",
     "The scale of the non-linear sigma-model",
     &LHTPModel::f_, TeV, 1.*TeV, 0.0*TeV, 10.0*TeV,
     true, false, Interface::limited);

  static Parameter<LHTPModel,double> interfaceSinAlpha
    ("SinAlpha",
     "The parameter controlling the mixing in the top quark sector of the model",
     &LHTPModel::salpha_, sqrt(0.5), 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<LHTPModel,double> interfaceKappaQuark
    ("KappaQuark",
     "The parameter controlling the masses of the T-odd quarks",
     &LHTPModel::kappaQuark_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<LHTPModel,double> interfaceKappaLepton
    ("KappaLepton",
     "The parameter controlling the masses of the T-odd leptons",
     &LHTPModel::kappaLepton_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<LHTPModel,Energy> interfaceHiggsMass
    ("HiggsMass",
     "The mass of the lightest Higgs boson",
     &LHTPModel::mh_, GeV, 120.0*GeV, 100.0*GeV, 1000.0*GeV,
     false, false, Interface::limited);

  static Switch<LHTPModel,bool> interfaceApproximate
    ("Approximate",
     "Whether to use the full expression for the mases of the top quark"
     " and its partners or the second-order expansion in v/f.",
     &LHTPModel::approximate_, false, false, false);
  static SwitchOption interfaceApproximateYes
    (interfaceApproximate,
     "Yes",
     "Approximate",
     true);
  static SwitchOption interfaceApproximateNo
    (interfaceApproximate,
     "No",
     "Don't approximate",
     false);

  static Reference<LHTPModel,AbstractVSSVertex> interfaceVertexWHH
    ("Vertex/WHH",
     "Vertex for the interactions of the electroweak gauge"
     " bosons and two Higgs bosons.",
     &LHTPModel::WHHVertex_, false, false, true, false, false);

}

void LHTPModel::doinit() {
  addVertex(WHHVertex_);
  StandardModel::doinit();
  using Constants::pi;
  // compute the parameters of the model
  // W and Z masses
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  // Energy mz(getParticleData(ParticleID::Z0)->mass());
  // couplings g and g'
  // double ee = sqrt(4.*pi*alphaEM(sqr(mz)));
  double ee = sqrt(4.*pi*alphaEMMZ());
  double sw(sqrt(sin2ThetaW())),cw(sqrt(1.-sin2ThetaW()));
  g_  = ee/sw;
  gp_ = ee/cw;
  // vev
  v_ = 2.*mw/g_;
  double vf(sqr(v_/f_));
  // calculate masses of the new particles from input
  // and SM parameters
  // masses of the new gauge bosons (MWH = MZH)
  Energy MAH = gp_*f_*sqrt(0.2)*(1.-0.625*vf);
  Energy MZH = g_ *f_*          (1.-0.125*vf);
  // mixings
  sthetaH_ = 1.25*g_*gp_/(5.*sqr(g_)-sqr(gp_))*vf;
  cthetaH_ = sqrt(1.-sqr(sthetaH_));
  // masses of the new top quarks
  Energy MTp,MTm;
  topMixing(MTp,MTm);
  // mixings in the top sector
  sL_ = sqr(salpha_)*v_/f_;
  cL_ = sqrt(1.-sqr(sL_));
  sR_ = salpha_*(1.-0.5*sqr(calpha_)*(sqr(calpha_)-sqr(salpha_))*vf);
  cR_ = sqrt(1.-sqr(sR_));
  // masses of the T-odd fermions
  Energy Mdm = sqrt(2.)*kappaQuark_ *f_;
  Energy Mum = sqrt(2.)*kappaQuark_ *f_*(1.-0.125*vf);
  Energy Mlm = sqrt(2.)*kappaLepton_*f_;
  Energy Mnm = sqrt(2.)*kappaLepton_*f_*(1.-0.125*vf);
  // masses of the triplet higgs
  Energy MPhi = sqrt(2.)*mh_*f_/v_;
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
  resetMass( 25      , mh_  );
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
  double vf(sqr(v_/f_));
  Energy mt = getParticleData(ParticleID::t)->mass();
  calpha_ = sqrt(1.-sqr(salpha_));
  double sv(sin(sqrt(2.)*v_/f_)),cv(cos(sqrt(2.)*v_/f_));
  // first guess for Yukawa's based on second-order in v/f expansion
  double lambda1 = mt/v_/calpha_*(1.+(2.-3.*pow(salpha_,4))*vf/6.);
  double lambda2 = mt/v_/salpha_*(1.+(2.-3.*pow(calpha_,4))*vf/6.);
  // first guess for masses
  MTp = sqrt(sqr(lambda1)+sqr(lambda2))*f_*(1-0.5*vf*sqr(calpha_*salpha_));
  MTm = lambda2*f_;
  if(!approximate_) {
    // special case where denominator of tan 2 alpha eqn is zero
    if(abs(salpha_-sqrt(0.5))<1e-4) {
      double a = 0.25*(2.*sqr(sv)+sqr(1.+cv));
      double b = 0.5*(a+0.5*(sqr(sv)+0.5*sqr(1.+cv)));
      lambda1 = mt/f_*sqrt(1./b/(1.-sqrt(1.-0.5*a*sqr(sv/b))));
      lambda2 = sqrt(a)*lambda1;
    }
    // general case using GSL
    else {
      double ca = sqrt(1.-sqr(salpha_));
      double ta = salpha_/ca;
      double tan2a = 2.*ta/(1.-sqr(ta));
      const gsl_multiroot_fsolver_type *T;
      gsl_multiroot_fsolver *s;
      int status;
      size_t iter=0;
      const size_t n=2;
      struct tparams p = {v_,f_,mt,tan2a};
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
    MTp = sqrt(sqr(f_)*delta*(1.+det));
    MTm = lambda2*f_;
  }
  // beta mixing angle
  double beta = 0.5*atan(2.*sqrt(2.)*sqr(lambda1)*sv*(1.+cv)/
			 (4.*sqr(lambda2)+sqr(1.+cv)*sqr(lambda1)-2.*sqr(lambda1)*sv));
  sbeta_ = sin(beta);
  cbeta_ = cos(beta);
}
