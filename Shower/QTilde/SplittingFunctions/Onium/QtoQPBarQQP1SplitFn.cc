// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQPBarQQP1SplitFn class.
//

#include "QtoQPBarQQP1SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const double QtoQPBarQQP1SplitFn::pOver_ = 1.;

IBPtr QtoQPBarQQP1SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQPBarQQP1SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQPBarQQP1SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(R02_,GeV*GeV2) << fixedAlphaS_;
}

void QtoQPBarQQP1SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(R02_,GeV*GeV2) >> fixedAlphaS_;
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<QtoQPBarQQP1SplitFn,Sudakov1to2FormFactor>
  describeHerwigQtoQPBarQQP1SplitFn("Herwig::QtoQPBarQQP1SplitFn", "HwOniumShower.so");

void QtoQPBarQQP1SplitFn::Init() {

  static ClassDocumentation<QtoQPBarQQP1SplitFn> documentation
    ("The QtoQPBarQQP1SplitFn class implements the splitting function for q -> qbar' (qq')_1");

  static Parameter<QtoQPBarQQP1SplitFn,Energy3> interfaceR02
    ("R0Squared",
     "The radial wavefunction at the origin squared",
     &QtoQPBarQQP1SplitFn::R02_, GeV*GeV2, pow<3,1>(0.41*GeV), 0.0*GeV*GeV2, 10.0*GeV*GeV2,
     false, false, Interface::limited);
  
  static Parameter<QtoQPBarQQP1SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQPBarQQP1SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);

}

void QtoQPBarQQP1SplitFn::guesstz(Energy2 t1,unsigned int iopt,
				  const IdList &ids,
				  double enhance,bool ident,
				  double detune, 
				  Energy2 &t_main, double &z_main) {
  unsigned int pdfopt = iopt!=1 ? 0 : pdfFactor();
  double lower = integOverP(zLimits().first ,ids,pdfopt);
  double upper = integOverP(zLimits().second,ids,pdfopt);
  Energy M = ids[0]->mass()+ids[1]->mass();
  double a2 = ids[1]->mass()/M;
  double aS2 = fixedAlphaS_ < 0 ? sqr(alpha()->overestimateValue()) : sqr(fixedAlphaS_);
  Energy2 pre = 1./9.*aS2*R02_/sqr(a2)/M/Constants::pi;
  Energy2 c = (upper - lower) * colourFactor() * pre * enhance * detune;
  double r = UseRandom::rnd();
  assert(iopt<=2);
  if(iopt==1) {
    c *= pdfMax();
    //symmetry of FS gluon splitting
    if(ident) c*= 2;
  }
  else if(iopt==2) c*=-1.;
  // guess t
  t_main = t1/(1.-t1/c*log(r));
  // guess z
  z_main = invIntegOverP(lower + UseRandom::rnd()*(upper - lower),ids,pdfopt);
}

double QtoQPBarQQP1SplitFn::ratioP(const double z, const Energy2 t,
				   const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double W0 = z*(17.-22.*z+9.*sqr(z))/sqr(1.+z);
  double W1 = 8*(3.-8.*z+sqr(z))/(1+z);
  double W2 = -48.;
  double ratio =(W0+r*W1+sqr(r)*W2)/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQPBarQQP1SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQPBarQQP1SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQPBarQQP1SplitFn::matrixElement(const double z, const Energy2 t, 
					      const IdList & ids, const double phi, bool) {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double rz=sqrt(z);
  double r2=sqrt(2.);
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t-sqr(m*(1.+z)));
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  (*kernal)(0,1,0) = -2.*r2*phase*double(pT/m)*r/((1.-z)*rz);
  (*kernal)(0,1,1) = rz*((3.-z)/(1.+z)-8.*r/(1.-z));
  (*kernal)(0,1,2) = 2.*r2*double(pT/m)*r*rz/(phase*(1.-z));
  (*kernal)(0,0,0) = 2.*r2*(1.-z)/(1.+z)*(z + r*(1.+z))/rz;
  (*kernal)(0,0,1) = 0.;
  (*kernal)(0,0,2) = 0.;
  (*kernal)(1,1,0) = 0.;
  (*kernal)(1,1,1) = 0.;
  (*kernal)(1,1,2) = 2.*r2*(1.-z)/(1.+z)*(z + r*(1.+z))/rz;
  (*kernal)(1,0,0) = -2.*r2*r*phase*double(pT/m)*rz/(1.-z);  
  (*kernal)(1,0,1) = rz*((3.-z)/(1.+z)-8.*r/(1.-z));
  (*kernal)(1,0,2) = 4.*r*r2*0.5*double(pT/m)/(phase*(1.-z)*rz);
  return kernal;
}


