// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KMatrix class.
//

#include "KMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "Herwig/Utilities/GSLIntegrator.h"

using namespace Herwig;

KMatrix::KMatrix(FlavourInfo flavour, vector<Channels> channels,
		 vector<Energy2> poles, vector<vector<Energy> > g)
  : flavour_(flavour), channels_(channels), poles_(poles), g_(g),
    initTable_(false), rho0_(0.), n_(3.),
    en_({0.311678*GeV2, 0.321678*GeV2, 0.331678*GeV2, 0.341678*GeV2, 0.351678*GeV2, 0.361678*GeV2, 0.371678*GeV2,
	 0.381678*GeV2, 0.391678*GeV2, 0.401678*GeV2, 0.411678*GeV2, 0.421678*GeV2, 0.431678*GeV2, 0.441678*GeV2,
	 0.451678*GeV2, 0.461678*GeV2, 0.471678*GeV2, 0.481678*GeV2, 0.491678*GeV2, 0.501678*GeV2, 0.511678*GeV2,
	 0.521678*GeV2, 0.531678*GeV2, 0.541678*GeV2, 0.551678*GeV2, 0.561678*GeV2, 0.571678*GeV2, 0.581678*GeV2,
	 0.591678*GeV2, 0.601678*GeV2, 0.611678*GeV2, 0.621678*GeV2, 0.631678*GeV2, 0.641678*GeV2, 0.651678*GeV2,
	 0.661678*GeV2, 0.671678*GeV2, 0.681678*GeV2, 0.691678*GeV2, 0.701678*GeV2, 0.711678*GeV2, 0.721678*GeV2,
	 0.731678*GeV2, 0.741678*GeV2, 0.751678*GeV2, 0.761678*GeV2, 0.771678*GeV2, 0.781678*GeV2, 0.791678*GeV2,
	 0.801678*GeV2, 0.811678*GeV2, 0.821678*GeV2, 0.831678*GeV2, 0.841678*GeV2, 0.851678*GeV2, 0.861678*GeV2,
	 0.871678*GeV2, 0.881678*GeV2, 0.891678*GeV2, 0.901678*GeV2, 0.911678*GeV2, 0.921678*GeV2, 0.931678*GeV2,
	 0.941678*GeV2, 0.951678*GeV2, 0.961678*GeV2, 0.971678*GeV2, 0.981678*GeV2, 0.991678*GeV2, 1*GeV2}),
    rho_({0, 1.64899e-13, 6.98821e-12, 6.09957e-11, 2.79118e-10, 8.97354e-10, 2.30904e-09, 5.09766e-09, 1.00643e-08,
	  1.82502e-08, 3.09555e-08, 4.97545e-08, 7.65082e-08, 1.13375e-07, 1.62821e-07, 2.27625e-07, 3.10888e-07,
	  4.16043e-07, 5.46856e-07, 7.07438e-07, 9.02248e-07, 1.13611e-06, 1.4142e-06, 1.74209e-06, 2.12573e-06,
	  2.57145e-06, 3.08602e-06, 3.67661e-06, 4.35082e-06, 5.11674e-06, 5.98288e-06, 6.95827e-06, 8.05243e-06,
	  9.27542e-06, 1.06378e-05, 1.21508e-05, 1.38263e-05, 1.56765e-05, 1.77145e-05, 1.99543e-05, 2.24101e-05,
	  2.50974e-05, 2.80322e-05, 3.12314e-05, 3.47128e-05, 3.84953e-05, 4.25987e-05, 4.70438e-05, 5.18526e-05,
	  5.70484e-05, 6.26555e-05, 6.86999e-05, 7.52087e-05, 8.22109e-05, 8.97367e-05, 9.78184e-05, 0.00010649,
	  0.000115787, 0.000125749, 0.000136414, 0.000147826, 0.000160031, 0.000173076, 0.000187012, 0.000201893,
	  0.000217776, 0.000234723, 0.000252798, 0.00027207, 0.000289072})
{}

void KMatrix::persistentOutput(PersistentOStream & os) const {
  os << ounit(poles_,GeV2) << ounit(g_,GeV)
     << ounit(mPiPlus_,GeV) << ounit(mPi0_,GeV)
     << ounit(mKPlus_,GeV) << ounit(mK0_,GeV) << ounit(mEta_,GeV)
     << ounit(mEtaPrime_,GeV)
     << initTable_ << rho0_ << n_ << ounit(en_,GeV2) << rho_;
}

void KMatrix::persistentInput(PersistentIStream & is, int) {
  is >> iunit(poles_,GeV2) >> iunit(g_,GeV)
     >> iunit(mPiPlus_,GeV) >> iunit(mPi0_,GeV)
     >> iunit(mKPlus_,GeV) >> iunit(mK0_,GeV) >> iunit(mEta_,GeV)
     >> iunit(mEtaPrime_,GeV)
     >> initTable_ >> rho0_ >> n_ >> iunit(en_,GeV2) >> rho_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<KMatrix,Interfaced>
describeHerwigKMatrix("Herwig::KMatrix", "Herwig.so");

void KMatrix::Init() {

  static ClassDocumentation<KMatrix> documentation
    ("The KMatrix class provides a base class for the implementation of "
     "K-matrix parameterizations in Herwig");

  static Switch<KMatrix,bool> interfaceInitializeTable
    ("InitializeTable",
     "Initialize the table for the four pion phase space",
     &KMatrix::initTable_, false, false, false);
  static SwitchOption interfaceInitializeTableYes
    (interfaceInitializeTable,
     "Yes",
     "Initialize the table",
     true);
  static SwitchOption interfaceInitializeTableNo
    (interfaceInitializeTable,
     "No",
     "Don't initialize the table",
     false);

    static Parameter<KMatrix,double> interfacePower
    ("Power",
     "Power for the 4 pion phase space",
     &KMatrix::n_, 3.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Command<KMatrix> interfaceSetPoles
    ("SetPoles",
     "Set the values of the poles",
     &KMatrix::setPoles, false);

  static Command<KMatrix> interfaceSetCouplings
    ("SetCouplings",
     "Set the values of the couplings",
     &KMatrix::setCouplings, false);

}

namespace {

class KMatrix4PiFunction {

public :

  KMatrix4PiFunction(Energy mrho, Energy Gamma, Energy mp) : M(mrho), mpi(mp) {
    gamma = Gamma/pow(1.-4.*sqr(mpi/M),1.5);
  }

  double operator() (Energy2 s2) const {
    Energy gam=gamma*pow(1.-4.*sqr(mpi)/s2,1.5);
    return M*gam/(sqr(s2-sqr(M))+sqr(M*gam))*sqrt(sqr(s+s1-s2)-4.*s*s1);
  }
  /** Argument type for the GSLIntegrator */
  typedef Energy2 ArgType;
  /** Return type for the GSLIntegrator */
  typedef double ValType;

  Energy M;
  Energy gamma;
  Energy mpi;
  Energy2 s1,s;
};
  
class KMatrix4PiIntegrand {
  
public :

  KMatrix4PiIntegrand(Energy mrho, Energy Gamma, Energy mp) :
    function(mrho,Gamma,mp)
  {}

  InvEnergy2 operator() (Energy2 s1) const {
    function.s  = s;
    function.s1 = s1;
    Energy2 upp = sqr(sqrt(s)-sqrt(s1));
    auto inner = gauss.value(function,4.*sqr(function.mpi),upp);
    Energy gam=function.gamma*pow(1.-4.*sqr(function.mpi)/s1,1.5);
    return function.M*gam/(sqr(s1-sqr(function.M))+sqr(function.M*gam))*inner/s/sqr(Constants::pi);
  }
  /** Argument type for the GSLIntegrator */
  typedef Energy2 ArgType;
  /** Return type for the GSLIntegrator */
  typedef InvEnergy2 ValType;

  GSLIntegrator gauss;
  mutable KMatrix4PiFunction function;
  Energy2 s;

};
}

void KMatrix::doinit() {
  Interfaced::doinit();
  // The charged pion mass
  mPiPlus_=getParticleData(ParticleID::piplus)->mass();
  // The neutral pion mass
  mPi0_=getParticleData(ParticleID::pi0)->mass();
  // The charged kaon mass
  mKPlus_=getParticleData(ParticleID::Kplus)->mass();
  // The neutral kaon mass
  mK0_=getParticleData(ParticleID::K0)->mass();
  // The eta mass
  mEta_=getParticleData(ParticleID::eta)->mass();
  // The eta' mass
  mEtaPrime_=getParticleData(ParticleID::etaprime)->mass();
  // set up of the four pion phase space
  if(initTable_) {
    tPDPtr rho = getParticleData(113);
    KMatrix4PiIntegrand integrand(rho->mass(),rho->width(),mPiPlus_);
    Energy2 s=16.*sqr(mPiPlus_);
    GSLIntegrator gauss;
    en_.clear();
    en_.reserve(200);
    rho_.clear();
    rho_.reserve(200);
    while (s<=GeV2) {
      integrand.s=s;
      double value = gauss.value(integrand,4.*sqr(mPiPlus_),sqr(sqrt(s)-2.*mPiPlus_));
      en_.push_back(s);
      rho_.push_back(value);
      s+=0.01*GeV2;
    }
    integrand.s=GeV2;
    double value = gauss.value(integrand,4.*sqr(mPiPlus_),sqr(GeV-2.*mPiPlus_));
    en_ .push_back(GeV2);
    rho_.push_back(value);
  }
  rho0_ = pow(1.-16.*sqr(mPiPlus_)/GeV2,n_)/rho_.back();
}

ublas::matrix<Complex> KMatrix::rho(Energy2 s) const {
  size_t msize = channels_.size();
  ublas::matrix<Complex> rho=ublas::zero_matrix<Complex>(msize,msize);
  for(unsigned int iChan=0;iChan<msize;++iChan) {
    double val(0);
    Energy m1,m2;
    switch (channels_[iChan]) {
    case PiPi:
      m1=m2=mPiPlus_;
      break;
    case KK:
      m1 = m2 = mKPlus_;
      break;
    case EtaEta:
      m1 = m2 = mEta_;
      break;
    case EtaEtaPrime:
      m1= mEta_;
      m2 = mEtaPrime_;
      break;
    case KPi:
      m1 = mKPlus_;
      m2 = mPiPlus_;
      break;
    case KEta:
      m1 = mKPlus_;
      m2 = mEta_;
      break;
    case KEtaPrime:
      m1 = mKPlus_;
      m2 = mEtaPrime_;
      break;
    case FourPi:
      m1 = m2 = 2.*mPiPlus_;
      break;
    default:
      assert(false);
    }
    val=(1.-sqr(m1+m2)/s);
    if (channels_[iChan]==FourPi) {
      if(!inter_) inter_ = make_InterpolatorPtr(rho_,en_,3);
      if(s<GeV2)
	rho(iChan,iChan) = rho0_*(*inter_)(s);
      else
	rho(iChan,iChan) = pow(val,n_);
    }
    else {
      if(s>sqr(m1+m2))
	rho(iChan,iChan) = Complex(sqrt(val),0.);
      else
	rho(iChan,iChan) = Complex(0.,sqrt(abs(val)));
    }
  }
  return rho;
}

ublas::vector<Complex> KMatrix::
amplitudes(Energy2 s, ublas::vector<Complex> pVector, bool multiplyByPoles) const {
  static const Complex ii(0.,1.);
  const ublas::identity_matrix<Complex> I(channels_.size());
  double fact(1.);
  if(multiplyByPoles) {
    for (Energy2 pole : poles_ ) fact *= 1.-s/pole;
  }
  // matrix for which we need the inverse
  ublas::matrix<Complex> m = fact*I-ii*prod(K(s,multiplyByPoles),rho(s));
  // inverse matrix
  ublas::matrix<Complex> inverse = ublas::identity_matrix<Complex>(m.size1());
  // 1x1 just a number
  if(m.size1()==1) {
    inverse(0,0) = 1./m(0,0);
  }
  // compute the inverse
  else {
    // create a permutation matrix for the LU-factorization
    ublas::permutation_matrix<std::size_t>  pm(m.size1());
    // perform LU-factorization
    int res = ublas::lu_factorize(m,pm);
    if( res != 0 ) {
      cerr << "problem with factorization\n";
      exit(1);
    }
    // back substitute to get the inverse
    ublas::lu_substitute(m, pm, inverse);
  }
  // compute the amplitudes
  return prod(inverse,pVector);
}

string KMatrix::setPoles(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  unsigned int npole = stoi(stype);
  poles_.resize(npole);
  for(unsigned int ix=0;ix<npole;++ix) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    poles_[ix] = stof(stype)*GeV2;
  }
  // success
  return "";
}

string KMatrix::setCouplings(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  unsigned int npole = stoi(stype);
  // parse second bit of the string
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  unsigned int nchannel = stoi(stype);
  g_.resize(npole);
  for(unsigned int ix=0;ix<npole;++ix) {
    g_[ix].resize(nchannel);
    for(unsigned int iy=0;iy<nchannel;++iy) {
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      g_[ix][iy] = stof(stype)*GeV;
    }
  }
  // success
  return "";
}
