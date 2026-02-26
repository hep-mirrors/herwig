// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HalfHalfOneEWSplitFn class.
//

#include "HalfHalfOneEWSplitFn.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

using namespace Herwig;

IBPtr HalfHalfOneEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr HalfHalfOneEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void HalfHalfOneEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << gZ_ << gWL_ << _couplingValueLeftIm << _couplingValueLeftRe << _couplingValueRightIm
     << _couplingValueRightRe << longitudinalEWScheme_ << _cG
     << _yLeftIm << _yLeftRe << _yRightIm << _yRightRe;
}

void HalfHalfOneEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> gZ_ >> gWL_ >> _couplingValueLeftIm >> _couplingValueLeftRe >> _couplingValueRightIm
     >> _couplingValueRightRe >> longitudinalEWScheme_ >> _cG
     >> _yLeftIm >> _yLeftRe >> _yRightIm >> _yRightRe;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<HalfHalfOneEWSplitFn,Sudakov1to2FormFactor>
describeHerwigHalfHalfOneEWSplitFn("Herwig::HalfHalfOneEWSplitFn", "HwShower.so");

void HalfHalfOneEWSplitFn::Init() {

  static ClassDocumentation<HalfHalfOneEWSplitFn> documentation
    ("The HalfHalfOneEWSplitFn class implements the splitting q->qW and q->qZ");

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceCouplingValueLeftIm
    ("CouplingValue.Left.Im",
     "The numerical value (imaginary part) of the left-handed splitting coupling to be imported for BSM splittings",
     &HalfHalfOneEWSplitFn::_couplingValueLeftIm, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceCouplingValueLeftRe
    ("CouplingValue.Left.Re",
     "The numerical value (real part) of the left-handed splitting coupling to be imported for BSM splittings",
     &HalfHalfOneEWSplitFn::_couplingValueLeftRe, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceCouplingValueRightIm
    ("CouplingValue.Right.Im",
     "The numerical value (imaginary part) of the right-handed splitting coupling to be imported for BSM splittings",
     &HalfHalfOneEWSplitFn::_couplingValueRightIm, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceCouplingValueRightRe
    ("CouplingValue.Right.Re",
     "The numerical value (real part) of the right-handed splitting coupling to be imported for BSM splittings",
     &HalfHalfOneEWSplitFn::_couplingValueRightRe, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Switch<HalfHalfOneEWSplitFn,unsigned int> interfaceLongitudinalEWScheme
    ("LongitudinalEWScheme",
    "EW splitting scheme: 0 = Subtraction, 1 = GI",
    &HalfHalfOneEWSplitFn::longitudinalEWScheme_, 0, false, false);
  static SwitchOption interfaceSubtractionEWScheme
    (interfaceLongitudinalEWScheme,
    "Subtraction",
    "Using Subtraction picture in q->q'V EW splittings",
    0);
  static SwitchOption interfaceGIEWScheme
    (interfaceLongitudinalEWScheme,
    "GaugeInvariant",
    "Using gauge invariant picture in q->q'V EW splittings",
    1);

  static Parameter<HalfHalfOneEWSplitFn,double> interfaceCG
    ("GI.cG",
    "Relative weight cG multiplying the Goldstone piece in GI (V + cG * G)",
    &HalfHalfOneEWSplitFn::_cG, 1.0, -1.0e3, +1.0e3,
    false, false, Interface::limited);

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceGBYukawaLeftIm
    ("GBYukawa.Left.Im",
    "The numerical value (imaginary part) of the left-handed Goldstone Yukawa coupling",
    &HalfHalfOneEWSplitFn::_yLeftIm, 0.0, -1.0E6, +1.0E6,
    false, false, Interface::limited);

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceGBYukawaLeftRe
    ("GBYukawa.Left.Re",
    "The numerical value (real part) of the left-handed Goldstone Yukawa coupling",
    &HalfHalfOneEWSplitFn::_yLeftRe, 0.0, -1.0E6, +1.0E6,
    false, false, Interface::limited);

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceGBYukawaRightIm
    ("GBYukawa.Right.Im",
    "The numerical value (imaginary part) of the right-handed Goldstone Yukawa coupling",
    &HalfHalfOneEWSplitFn::_yRightIm, 0.0, -1.0E6, +1.0E6,
    false, false, Interface::limited);

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceGBYukawaRightRe
    ("GBYukawa.Right.Re",
    "The numerical value (real part) of the right-handed Goldstone Yukawa coupling",
    &HalfHalfOneEWSplitFn::_yRightRe, 0.0, -1.0E6, +1.0E6,
    false, false, Interface::limited);

}

void HalfHalfOneEWSplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sm->sin2ThetaW();
  // left-handled W coupling
  gWL_ = 1./sqrt(2.*sw2);
  // Z couplings
  double fact = 0.25/sqrt(sw2*(1.-sw2));
  for(int ix=1;ix<4;++ix) {
    gZ_[2*ix-1]  = make_pair(fact*(sm->vd()  + sm->ad()),
			     fact*(sm->vd()  - sm->ad() ));
    gZ_[2*ix  ]  = make_pair(fact*(sm->vu()  + sm->au() ),
			     fact*(sm->vu()  - sm->au() ));
    gZ_[2*ix+9 ] = make_pair(fact*(sm->ve()  + sm->ae() ),
			     fact*(sm->ve()  - sm->ae() ));
    gZ_[2*ix+10] = make_pair(fact*(sm->vnu() + sm->anu()),
			     fact*(sm->vnu() - sm->anu()));
  }
}

void HalfHalfOneEWSplitFn::getCouplings(Complex & gL, Complex & gR, const IdList & ids) const {
  if(_couplingValueLeftRe!=0||_couplingValueLeftIm!=0||_couplingValueRightRe!=0||_couplingValueRightIm!=0) {
    double e = sqrt(4.*Constants::pi*generator()->standardModel()
              ->alphaEM(sqr(getParticleData(ParticleID::Z0)->mass())));
    // a factor e is factored out since its already accounted for
    gL = Complex(_couplingValueLeftRe,_couplingValueLeftIm)/e;
    gR = Complex(_couplingValueRightRe,_couplingValueRightIm)/e;
  }
  else if(ids[2]->id()==ParticleID::Z0) {
    map<long,pair<double,double> >::const_iterator it = gZ_.find(abs(ids[0]->id()));
    assert(it!=gZ_.end());
    gL = Complex(0.,it->second.first);
    gR = Complex(0.,it->second.second);
  }
  else if(abs(ids[2]->id())==ParticleID::Wplus) {
    gL = Complex(0.,gWL_);
  }
  else
    assert(false);
  if(ids[0]->id()<0) swap(gL,gR);
}

void HalfHalfOneEWSplitFn::getGBYukawas(Complex & yL, Complex & yR, const IdList & ids) const {
  // Default: derive from SM masses/couplings unless override is set
  if (_yLeftRe!=0 || _yLeftIm!=0 || _yRightRe!=0 || _yRightIm!=0) {
    yL = Complex(_yLeftRe,  _yLeftIm);
    yR = Complex(_yRightRe, _yRightIm);
    return;
  }

  // SM inputs
  tcSMPtr sm = generator()->standardModel();
  const Energy MW = getParticleData(ParticleID::Wplus)->mass();
  const double sw2 = sm->sin2ThetaW();
  const double sw  = sqrt(sw2);
  // SU(2) coupling g from e = g * sin(thetaW)
  const double e   = sqrt(4.0*Constants::pi*sm->alphaEM(sqr(getParticleData(ParticleID::Z0)->mass())));
  const double gSU2 = e / sw;

  const long idV = ids[2]->id();
  const Energy m0 = ids[0]->mass();
  const Energy m1 = ids[1]->mass();

  if (abs(idV) == ParticleID::Wplus) {
    // q -> q' W:  yL =  g m_up / (sqrt(2) MW),   yR = - g m_down / (sqrt(2) MW)
    yL = Complex( gSU2 * (m0/MW) / sqrt(2.0), 0.0 );
    yR = Complex(-gSU2 * (m1/MW) / sqrt(2.0), 0.0 );
  }
  else if (idV == ParticleID::Z0) {
    // q -> q Z:   yL =  + g m_f / (2 MW),   yR = - g m_f / (2 MW)
    const double fac = gSU2 * (m0/MW) / 2.0;
    yL = Complex(+fac, 0.0);
    yR = Complex(-fac, 0.0);
  }
  else {
    // Photon or others
    yL = yR = Complex(0.0, 0.0);
  }
}

double HalfHalfOneEWSplitFn::ratioP(const double z, const Energy2 t,
				    const IdList & ids, const bool mass,
				    const RhoDMatrix & rho) const {
  Complex gL(0.,0.), gR(0.,0.), yL(0.,0.), yR(0.,0.);
  getCouplings(gL,gR,ids);
  getGBYukawas(yL,yR,ids);
  if (longitudinalEWScheme_ == 1) {
    yL *= _cG;
    yR *= _cG;
  }
  double gL2 = norm(gL), gR2 = norm(gR);
  double yL2 = norm(yL), yR2 = norm(yR);

  Energy m0 = 0.*GeV, m1 = 0.*GeV, m2 = 0.*GeV;
  double val = 0.0;

  if (longitudinalEWScheme_ == 0) {
    val = (gL2*abs(rho(0,0)) + gR2*abs(rho(1,1)))*(1.+sqr(z));

    if (mass) {
      m0 = ids[0]->mass();
      m1 = ids[1]->mass();
      m2 = ids[2]->mass();
      double m0t = m0/sqrt(t), m1t = m1/sqrt(t), m2t = m2/sqrt(t);

      val += (gL2*abs(rho(0,0)) + gR2*abs(rho(1,1)))*
               ((1.+sqr(z))*sqr(m0t) - (1.+z)*sqr(m1t) - (1.-z)*sqr(m2t))
           + (gR2*abs(rho(0,0)) + gL2*abs(rho(1,1)))*z*(1.-z)*sqr(m0t)
           - 2.*real(gR*conj(gL))*(abs(rho(1,1)) + abs(rho(0,0)))*
               (1.-z)*m0t*m1t;
    }

    val /= 2.0 * max(gL2,gR2);
  }
  // else {
  //   double overP = (2.*max(gL2,gR2) + max(yL2,yR2) )/(1.-z);
  //   val = ((1.+sqr(z))/(1.-z) * (gL2*abs(rho(0,0)) + gR2*abs(rho(1,1)))
  //       + sqr(1.-z)/(1.-z) * (yL2*abs(rho(0,0)) + yR2*abs(rho(1,1))));
  //
  //   if (mass) {
  //     m0 = ids[0]->mass();
  //     m1 = ids[1]->mass();
  //     m2 = ids[2]->mass();
  //     double m0t = m0/sqrt(t), m1t = m1/sqrt(t), m2t = m2/sqrt(t);
  //
  //     val += (-((-1 + z)*(-(yL*((pow(m0t,2) + pow(m1t,2) - pow(m2t,2))*yL
  //         + 2*m0t*m1t*yR)) - pow(m0t,2)*(gR2 - pow(yL,2)
  //         + pow(yR,2))*z)*abs(rho(0,0))) + 2*sqrt(gL2)*(-(sqrt(gR2)*m0t*m1t*(-1 + z))
  //         + sqrt(2)*m2t*(m1t*yL + m0t*yR*z))*abs(rho(0,0)) - 2*sqrt(gL2)*sqrt(gR2)*m0t*m1t*(-1 + z)*abs(rho(1,1))
  //         + (2*sqrt(2)*sqrt(gR2)*m2t*(m1t*yR + m0t*yL*z) + (-1 + z)*(yR*(2*m0t*m1t*yL
  //         + (pow(m0t,2) + pow(m1t,2) - pow(m2t,2))*yR)
  //         + pow(m0t,2)*(yL - yR)*(yL + yR)*z) + gR2*(-(pow(m2t,2)*(-1 + z))
  //         + pow(m1t,2)*(1 + z) - pow(m0t,2)*(1 + pow(z,2))))*abs(rho(1,1))
  //         + gL2*((-(pow(m2t,2)*(-1 + z)) + pow(m1t,2)*(1 + z)
  //         - pow(m0t,2)*(1 + pow(z,2)))*abs(rho(0,0)) + pow(m0t,2)*(-1 + z)*z*abs(rho(1,1))))
  //         /(pow(-1 + z,3)*pow(z,2));
  //   }
  //   val /= overP;
  // }
  else {
  double oneMz = 1.0 - z;
  double overP = (2.*max(gL2,gR2) + max(yL2,yR2))/oneMz;

  double r00 = abs(rho(0,0));
  double r11 = abs(rho(1,1));

  //cerr << "rho00=" << abs(rho(0,0)) << " rho11=" << abs(rho(1,1)) << "\n";

  // massless GI (cG already absorbed in yL,yR)
  val = ((1.+sqr(z))/oneMz) * (gL2*r00 + gR2*r11)
      + oneMz              * (yL2*r00 + yR2*r11);

  if (mass) {
    m0 = ids[0]->mass();
    m1 = ids[1]->mass();
    m2 = ids[2]->mass();
    double m0t = m0/sqrt(t), m1t = m1/sqrt(t), m2t = m2/sqrt(t);

    double Re_gLgR = real(gL*conj(gR));
    double Re_gLyL = real(gL*conj(yL));
    double Re_gLyR = real(gL*conj(yR));
    double Re_gRyL = real(gR*conj(yL));
    double Re_gRyR = real(gR*conj(yR));
    double Re_yLyR = real(yL*conj(yR));

    double At = sqr(m0t) + sqr(m1t) - sqr(m2t);

    // "minus" (rho00) and "plus" (rho11) pieces
    double Nm =
        oneMz * ( - (At*yL2 + 2.*m0t*m1t*Re_yLyR)
                  - sqr(m0t)*(gR2 - yL2 + yR2)*z )
      + gL2 * ( sqr(m2t)*oneMz - sqr(m1t)*(1.+z) + sqr(m0t)*(1.+sqr(z)) )
      - 2.*sqrt(2.)*m2t * ( m1t*Re_gLyL + m0t*z*Re_gLyR );

    double Np =
        gR2 * ( sqr(m2t)*oneMz - sqr(m1t)*(1.+z) + sqr(m0t)*(1.+sqr(z)) )
      - oneMz * ( 2.*m0t*m1t*Re_yLyR + At*yR2 + sqr(m0t)*(yL2 - yR2)*z )
      - 2.*sqrt(2.)*m2t * ( m1t*Re_gRyR + m0t*z*Re_gRyL )
      - gL2 * sqr(m0t) * oneMz * z;

    double Nint =
      2.*oneMz * m0t*m1t * Re_gLgR * (r00 + r11);

    // this matches your qt^2 = t*z*(1-z) rewriting
    val += (r00*Nm + r11*Np + Nint) / (pow(oneMz,3)*sqr(z));
  }

  val /= overP;
}

  return val;
}

double HalfHalfOneEWSplitFn::integOverP(const double z,
					const IdList & ids,
					unsigned int PDFfactor) const {
  Complex gL(0.,0.), gR(0.,0.), yL(0.,0.), yR(0.,0.);
  getCouplings(gL,gR,ids);

  double gL2 = norm(gL), gR2 = norm(gR);
  double coef;

  if (longitudinalEWScheme_ == 0) {
    coef = 2.*max(gL2,gR2);
  }
  else {
    getGBYukawas(yL,yR,ids);
    yL *= _cG;
    yR *= _cG;
    double yL2 = norm(yL), yR2 = norm(yR);
    coef = 2.*max(gL2,gR2) + max(yL2,yR2);
  }

  switch (PDFfactor) {
  case 0:
    return -coef*Math::log1m(z);
  case 1:
    return  coef*log(z/(1.-z));
  case 2:
    return  coef/(1.-z);
  case 3:
  default:
    throw Exception() << "HalfHalfOneEWSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double HalfHalfOneEWSplitFn::invIntegOverP(const double r, const IdList & ids,
					   unsigned int PDFfactor) const {
  Complex gL(0.,0.), gR(0.,0.), yL(0.,0.), yR(0.,0.);
  getCouplings(gL,gR,ids);

  double gL2 = norm(gL), gR2 = norm(gR);
  double coef;

  if (longitudinalEWScheme_ == 0) {
    coef = 2.*max(gL2,gR2);
  }
  else {
    getGBYukawas(yL,yR,ids);
    yL *= _cG;
    yR *= _cG;
    double yL2 = norm(yL), yR2 = norm(yR);
    coef = 2.*max(gL2,gR2) + max(yL2,yR2);
  }

  switch (PDFfactor) {
  case 0:
    return 1. - exp(-r/coef);
  case 1:
    return 1./(1. - exp(-r/coef));
  case 2:
    return 1. - coef/r;
  case 3:
  default:
    throw Exception() << "HalfHalfOneEWSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool HalfHalfOneEWSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[2]->iSpin()==PDT::Spin1 && !(_couplingValueLeftRe==0 && _couplingValueLeftIm==0 && _couplingValueRightRe==0 && _couplingValueRightIm==0)) {
    if(ids[0]->iCharge()!=ids[1]->iCharge()+ids[2]->iCharge()) return false;
    if((abs(ids[0]->id())>=1 && abs(ids[0]->id())<=6) && (abs(ids[1]->id())>=1 && abs(ids[1]->id())<=6)) return true;
  }
  else if(ids[2]->id()==ParticleID::Z0) {
    if(ids[0]->id()==ids[1]->id() &&
       ((ids[0]->id()>=1 && ids[0]->id()<=6) || (ids[0]->id()>=11&&ids[0]->id()<=16) )) return true;
  }
  else if(abs(ids[2]->id())==ParticleID::Wplus) {
    if(!((ids[0]->id()>=1 && ids[0]->id()<=6) || (ids[0]->id()>=11&&ids[0]->id()<=16) )) return false;
    if(!((ids[1]->id()>=1 && ids[1]->id()<=6) || (ids[1]->id()>=11&&ids[1]->id()<=16) )) return false;
    if(ids[0]->id()+1!=ids[1]->id() && ids[0]->id()-1!=ids[1]->id()) return false;
    int out = ids[1]->iCharge()+ids[2]->iCharge();
    if(ids[0]->iCharge()==out) return true;
  }
  return false;
}

DecayMEPtr HalfHalfOneEWSplitFn::matrixElement(const double z, const Energy2 t,
                                               const IdList & ids, const double phi,
                                               bool timelike) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));

  Energy m0 = timelike ? ids[0]->mass() : ZERO;
  Energy m1 = timelike ? ids[1]->mass() : ZERO;
  Energy m2 = ids[2]->mass();
  Complex gL(0.,0.), gR(0.,0.), yL(0.,0.), yR(0.,0.);
  getCouplings(gL,gR,ids);
  getGBYukawas(yL,yR,ids);

  Energy den = sqrt(2.)*sqrt(t);
  Energy pt  = sqrt(z*(1.-z)*t+z*(1.-z)*sqr(m0)-(1.-z)*sqr(m1)-z*sqr(m2));
  Complex phase  = exp(Complex(0.,1.)*phi);
  Complex cphase = conj(phase);

  (*kernal)(1,1,2) = sqrt(2.)*gR*pt/sqrt(z)/(1.-z)*cphase/den;
  (*kernal)(1,1,0) = -sqrt(2.*z)*gR*pt/(1.-z)*phase/den;
  (*kernal)(1,0,2) = -sqrt(2.)*(gL*z*m0-gR*m1)/sqrt(z)/den;
  (*kernal)(1,0,0) = 0;
  (*kernal)(0,1,2) = 0;
  (*kernal)(0,1,0) = -sqrt(2.)*(gR*z*m0-gL*m1)/sqrt(z)/den;
  (*kernal)(0,0,2) = sqrt(2.*z)*gL*pt/(1.-z)*cphase/den;
  (*kernal)(0,0,0) = -sqrt(2.)*gL*pt/sqrt(z)/(1.-z)*phase/den;

  (*kernal)(0,0,1) = -2.*gL*sqrt(z)*m2/(1.-z)/den;
  (*kernal)(0,1,1) = 0.;
  (*kernal)(1,0,1) = 0.;
  (*kernal)(1,1,1) = -2.*gR*sqrt(z)*m2/(1.-z)/den;

  if(longitudinalEWScheme_ == 1) {
    (*kernal)(0,0,1)+= _cG*sqrt(2.)*(m0*yR*z + m1*yL)/sqrt(z)/den;
    (*kernal)(0,1,1)+= -_cG*sqrt(2.)*cphase*yL*pt/sqrt(z)/den;
    (*kernal)(1,0,1)+= -_cG*sqrt(2.)*phase*yR*pt/sqrt(z)/den;
    (*kernal)(1,1,1)+= _cG*sqrt(2.)*(m0*yL*z + m1*yR)/sqrt(z)/den;
  }

  // return the answer
  return kernal;
}
