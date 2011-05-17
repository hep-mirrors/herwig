// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneLoopFFAWZVertex class.
//

#include "OneLoopFFAWZVertex.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Looptools/clooptools.h"

using namespace Herwig;
namespace LT = Looptools;

OneLoopFFAWZVertex::OneLoopFFAWZVertex() :
  order_(0), EWscheme_(2), mZ_(ZERO), mW_(ZERO),muZ2_(ZERO),
  muZ_(ZERO),muW2_(ZERO),muW_(ZERO),mH_(ZERO),mH2_(ZERO),
  cw2_(0.),cw_(0.),sw2_(0.),sw_(0.),alphaEW_(0.),e_(0.),
  dZ_AA_(0.),dZ_ZA_(0.),dZ_W_(0.),dZ_AZ_(0.),dZ_ZZ_(0.),
  dmuW2_(ZERO),dmuZ2_(ZERO),de_(0.),dsw_(0.) {
  kinematics(true);
}

IBPtr OneLoopFFAWZVertex::clone() const {
  return new_ptr(*this);
}

IBPtr OneLoopFFAWZVertex::fullclone() const {
  return new_ptr(*this);
}

void OneLoopFFAWZVertex::setCoupling(Energy2,tcPDPtr aa,tcPDPtr bb,tcPDPtr cc) {
  int iferm=abs(aa->id());
  assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16));
  // photon
  if(cc->id()==ParticleID::gamma) {
    norm(e_);
    left (-ef_[iferm]);
    right(-ef_[iferm]);
  }
  // Z boson
  else if(cc->id()==ParticleID::Z0) {
    norm(e_);
    left (gl_[iferm]);
    right(gr_[iferm]);
  }
  else
    assert(false);
  // higher order piece if needed
  if(order_!=0) {
    pair<Complex,Complex> vertexCorrection = 
      renormalisedVertex(invariant(2,2),aa,bb,cc);
    left ( left()*vertexCorrection.first );
    right(right()*vertexCorrection.second);
  }
}

void OneLoopFFAWZVertex::doinit() {
  Looptools::ltini();
  // PDG codes for the particles
  // the quarks
  for(int ix=1;ix<7;++ix) {
    addToList(-ix, ix, 22);
  }
  // the leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix, ix, 22);
  }
  // PDG codes for the particles
  // the quarks
  for(int ix=1;ix<7;++ix) {
    addToList(-ix, ix, 23);
  }
  // the leptons
  for(int ix=11;ix<17;++ix) {
    addToList(-ix, ix, 23);
  }
  FFVVertex::doinit();
  // parameters for the EW corrections
  // convert the mass and width into the scheme of DH (2.11)
  tcPDPtr Z0    = getParticleData(ThePEG::ParticleID::Z0);
  tcPDPtr Wplus = getParticleData(ParticleID::Wplus);
  mZ_ = Z0   ->mass()/sqrt(1.+sqr(Z0   ->width()/Z0   ->mass()));
  mW_ = Wplus->mass()/sqrt(1.+sqr(Wplus->width()/Wplus->mass()));
  tcPDPtr higgs = getParticleData(ParticleID::h0);
  mH_ = higgs->mass();
  mH2_ = sqr(mH_);
  // the coupling
  if(EWscheme_==0)
    alphaEW_ = generator()->standardModel()->alphaEM();
  else if(EWscheme_==1)
    alphaEW_ = generator()->standardModel()->alphaEMMZ();
  else if(EWscheme_==2)
    alphaEW_ = sqrt(2.)*generator()->standardModel()->fermiConstant()*
      sqr(mW_)*(1.-sqr(mW_/mZ_))/Constants::pi;
  // electric charge
  e_ = sqrt(4.*alphaEW_*Constants::pi);
  // \todo need consistent calculation of the width
  // for the moment take value from ParticeData object (experimental value)
  Energy gamZ = Z0   ->width()/sqrt(1.+sqr(Z0   ->width()/Z0   ->mass()));
  Energy gamW = Wplus->width()/sqrt(1.+sqr(Wplus->width()/Wplus->mass()));
  // complex paramters
  muZ2_ = complex<Energy2>(sqr(mZ_),-gamZ*mZ_);
  muZ_  = sqrt(Complex(muZ2_/GeV2))*GeV;
  muW2_ = complex<Energy2>(sqr(mW_),-gamW*mW_);
  muW_  = sqrt(Complex(muW2_/GeV2))*GeV;
  // weak mixings
  cw2_ = muW2_/muZ2_;
  cw_  = sqrt(cw2_);
  sw2_ = 1.-cw2_   ;
  sw_  = sqrt(sw2_);
  // charges
  ef_.resize(17);
  gl_.resize(17);
  gr_.resize(17);
  for(unsigned int ix=0;ix<3;++ix) {
    // electric
    ef_[2*ix+1 ] = -1./3.;
    ef_[2*ix+2 ] =  2./3.;
    ef_[2*ix+11] = -1.   ;
    ef_[2*ix+12] =  0.   ;
    // left
    gl_[2*ix+1 ] = (-0.5 - sw2_*ef_[2*ix+1 ])/sw_/cw_;
    gl_[2*ix+2 ] = ( 0.5 - sw2_*ef_[2*ix+2 ])/sw_/cw_;
    gl_[2*ix+11] = (-0.5 - sw2_*ef_[2*ix+11])/sw_/cw_;
    gl_[2*ix+12] = ( 0.5 - sw2_*ef_[2*ix+12])/sw_/cw_;
    // right
    gr_[2*ix+1 ] = -sw_/cw_*ef_[2*ix+1 ];
    gr_[2*ix+2 ] = -sw_/cw_*ef_[2*ix+2 ];
    gr_[2*ix+11] = -sw_/cw_*ef_[2*ix+11];
    gr_[2*ix+12] = -sw_/cw_*ef_[2*ix+12];
  }
  // masses
  if(fermionMasses_.empty()) {
    fermionMasses_.resize(17,ZERO);
    for(long ix=1;ix<7;++ix)
      fermionMasses_[ix] = getParticleData(ix)->mass();
    for(long ix=11;ix<17;++ix)
      fermionMasses_[ix] = getParticleData(ix)->mass();
    fermionMasses_[11] = 5.1099892E-004*GeV;
    fermionMasses_[ 1] = 6.6000000E-002*GeV;
    fermionMasses_[ 2] = 6.6000000E-002*GeV;
    fermionMasses_[13] = 0.105658369*GeV;
    fermionMasses_[ 3] = 6.6E-002 *GeV;
    fermionMasses_[ 4] = 1.2*GeV;
    fermionMasses_[15] =  1.77699*GeV;
    fermionMasses_[ 5] = 4.3*GeV;
    fermionMasses_[ 6] = 174.3*GeV;
  }
  // renormalisation constants
  // boson field renormalisation
  dZ_AA_ = -DSigma_T_AA(ZERO);
  dZ_ZA_ = 2.*Sigma_T_AZ(ZERO)/muZ2_;
  dZ_W_  = -DSigma_T_W(sqr(mW_));
  dZ_AZ_ = -2.*Sigma_T_AZ(sqr(mZ_))/sqr(mZ_)+(muZ2_/sqr(mZ_)-1.)*dZ_ZA_;
  dZ_ZZ_ = -DSigma_T_ZZ(sqr(mZ_));
  // counter terms for gauge boson masses
  dmuW2_ = Sigma_T_W (sqr(mW_))+(sqr(mW_)-muW2_)*dZ_W_ ;
  dmuZ2_ = Sigma_T_ZZ(sqr(mZ_))+(sqr(mZ_)-muZ2_)*dZ_ZZ_;
  dsw_ = -0.5*cw2_/sw2_*(dmuW2_/muW2_-dmuZ2_/muZ2_);
  // de/e
  de_  = -0.5*dZ_AA_-0.5*sw_/cw_*dZ_ZA_;
  // resummation of large logs in schemes other than alpha(0)
  if(EWscheme_==1) {
    de_ -= 0.5*deltaAlphaMZ();
  }
  else if(EWscheme_==2) {
    Complex dr = -dZ_AA_-2.*dsw_ +(Sigma_T_W(ZERO)-dmuW2_)/muW2_
      +cw_/sw_*dZ_ZA_+0.25*alphaEW_/sw2_/Constants::pi*
      (6.+0.5*(7.-4.*sw2_)/sw2_*log(cw2_));
    de_ -= 0.5*real(dr);
  }
}

void OneLoopFFAWZVertex::dofinish() {
  Looptools::ltexi();
  FFVVertex::dofinish();
}

void OneLoopFFAWZVertex::doinitrun() {
  Looptools::ltini();
  FFVVertex::doinitrun();
}

void OneLoopFFAWZVertex::persistentOutput(PersistentOStream & os) const {
  os << EWscheme_ << ounit(mZ_,GeV) << ounit(mW_,GeV) << ounit(muZ2_,GeV2)
     << ounit(muZ_,GeV) << ounit(muW2_,GeV2) << ounit(muW_,GeV)
     << ounit(mH_,GeV) << ounit(mH2_,GeV2) << cw2_ << cw_ << sw2_ 
     << sw_ << ef_ << gl_ << gr_ << ounit(fermionMasses_,GeV) 
     << alphaEW_ << e_ << dZ_AA_ << dZ_ZA_ << dZ_W_ 
     << dZ_AZ_ << dZ_ZZ_ << ounit(dmuW2_,GeV2) << ounit(dmuZ2_,GeV2) 
     << de_ << dsw_;
}

void OneLoopFFAWZVertex::persistentInput(PersistentIStream & is, int) {
  is >> EWscheme_ >> iunit(mZ_,GeV) >> iunit(mW_,GeV) >> iunit(muZ2_,GeV2)
     >> iunit(muZ_,GeV) >> iunit(muW2_,GeV2) >> iunit(muW_,GeV)
     >> iunit(mH_,GeV) >> iunit(mH2_,GeV2) >> cw2_ >> cw_ >> sw2_ 
     >> sw_ >> ef_ >> gl_ >> gr_ >> iunit(fermionMasses_,GeV) 
     >> alphaEW_ >> e_ >> dZ_AA_ >> dZ_ZA_ >> dZ_W_ 
     >> dZ_AZ_ >> dZ_ZZ_ >> iunit(dmuW2_,GeV2) >> iunit(dmuZ2_,GeV2) 
     >> de_ >> dsw_;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<OneLoopFFAWZVertex,FFVVertex>
  describeHerwigOneLoopFFAWZVertex("Herwig::OneLoopFFAWZVertex", "OneLoopFFAWZVertex.so");

void OneLoopFFAWZVertex::Init() {

  static ClassDocumentation<OneLoopFFAWZVertex> documentation
    ("There is no documentation for the OneLoopFFAWZVertex class");

  static Switch<OneLoopFFAWZVertex,unsigned int> interfaceEWScheme
    ("EWScheme",
     "The scheme for the electroweak corrections",
     &OneLoopFFAWZVertex::EWscheme_, 2, false, false);
  static SwitchOption interfaceEWSchemealpha0
    (interfaceEWScheme,
     "alpha0",
     "The alpha0 scheme",
     0);
  static SwitchOption interfaceEWSchemealphaMZ
    (interfaceEWScheme,
     "alphaMZ",
     "The alpha(MZ) scheme ",
     1);
  static SwitchOption interfaceEWSchemeGmuScheme
    (interfaceEWScheme,
     "GmuScheme",
     "The G_mu scheme",
     2);

}

// eqn B1 of Denner
complex<Energy2> OneLoopFFAWZVertex::Sigma_T_AA_Boson(Energy2 k2) {
  // bosonic piece
  Complex B0WW = LT::B0C(0.     ,muW2_/GeV2,muW2_/GeV2);
  Complex BkWW = LT::B0C(double(k2/GeV2),muW2_/GeV2,muW2_/GeV2);
  complex<Energy2> output = (3.*k2+4.*muW2_)*BkWW - 4.*muW2_*B0WW;
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B1 of Denner differentiated w.r.t. k2
Complex OneLoopFFAWZVertex::DSigma_T_AA_Boson(Energy2 k2) {
  // bosonic piece
  Complex  BkWW =  LT::B0C(double(k2/GeV2),muW2_/GeV2,muW2_/GeV2);
  complex<InvEnergy2> DBkWW = LT::DB0C(double(k2/GeV2),muW2_/GeV2,muW2_/GeV2)/GeV2;
  Complex output = 3.*BkWW+Complex((3.*k2+4.*muW2_)*DBkWW);
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B2 of Denner
complex<Energy2> OneLoopFFAWZVertex::Sigma_T_AZ_Boson(Energy2 k2) {
  // bosonic piece
  Complex B0WW = LT::B0C(0.     ,muW2_/GeV2,muW2_/GeV2);
  Complex BkWW = LT::B0C(double(k2/GeV2),muW2_/GeV2,muW2_/GeV2);
  complex<Energy2> output = -1./3./sw_/cw_*
    (((9.*cw2_+0.5)*k2+4.*(3.*cw2_+1.)*muW2_)*BkWW
     -2.*(6.*cw2_-1.)*muW2_*B0WW+k2/3.);
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B2 of Denner differentiated w.r.t. k2
Complex OneLoopFFAWZVertex::DSigma_T_AZ_Boson(Energy2 k2) {
  // bosonic piece
  Complex  BkWW =  LT::B0C(double(k2/GeV2),muW2_/GeV2,muW2_/GeV2);
  complex<InvEnergy2> DBkWW = LT::DB0C(double(k2/GeV2),muW2_/GeV2,muW2_/GeV2)/GeV2;
  Complex output = -1./3./sw_/cw_*
    (Complex(((9.*cw2_+0.5)*k2+4.*(3.*cw2_+1.)*muW2_)*DBkWW)
     +(9.*cw2_+0.5)*BkWW);
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B3 of Denner
complex<Energy2> OneLoopFFAWZVertex::Sigma_T_ZZ_Boson(Energy2 k2) {
  // bosonic piece
  complex<Energy2> output;
  Complex B0WW = LT::B0C(0.,muW2_/GeV2,muW2_/GeV2);
  Complex B0ZZ = LT::B0C(0.,muZ2_/GeV2,muZ2_/GeV2);
  Complex B0HH = LT::B0C(0.,double(mH2_ /GeV2),double(mH2_ /GeV2));
  Complex B0ZH = LT::B0C(0.             ,muZ2_/GeV2,double(mH2_/GeV2));
  if(k2!=ZERO) {
    Complex BkWW = LT::B0C(double(k2/GeV2),muW2_/GeV2,muW2_/GeV2);
    Complex BkZH = LT::B0C(double(k2/GeV2),muZ2_/GeV2,double(mH2_/GeV2));
    output = 
      1./6./sw2_/cw2_*
      (BkWW*((18.*sqr(cw2_)+2.*cw2_-0.5)*k2
	     +2.*muW2_*(12.*sqr(cw2_)+8.*cw2_-5.))
       -2.*(12.*sqr(cw2_)-4.*cw2_+1.)*muW2_*B0WW+(4.*cw2_-1.)*k2/3.)
      +1./12./sw2_/cw2_*
      ((2.*mH2_-10.*muZ2_-k2)*BkZH-2.*muZ2_*B0ZZ-2.*mH2_*B0HH
       -sqr(muZ2_-mH2_)/k2*(BkZH-B0ZH)-2.*k2/3.);
  }
  else {
    complex<InvEnergy2> DB0ZH = LT::DB0C(0.,muZ2_/GeV2,double(mH2_/GeV2))/GeV2;
    output = 
      2./sw2_/cw2_*B0WW*muW2_*(2.*cw2_-1.)
      +1./12./sw2_/cw2_*
      ((2.*mH2_-10.*muZ2_)*B0ZH-2.*muZ2_*B0ZZ-2.*mH2_*B0HH
       -sqr(muZ2_-mH2_)*DB0ZH);
  }
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B3 of Denner differentiated w.r.t. k2
Complex OneLoopFFAWZVertex::DSigma_T_ZZ_Boson(Energy2 k2) {
  assert(k2!=ZERO);
  Complex              BkWW =  LT::B0C(double(k2/GeV2),muW2_/GeV2,muW2_/GeV2);
  complex<InvEnergy2> DBkWW = LT::DB0C(double(k2/GeV2),muW2_/GeV2,muW2_/GeV2)/GeV2;
  Complex              BkZH =  LT::B0C(double(k2/GeV2),muZ2_/GeV2,double(mH2_/GeV2));
  complex<InvEnergy2> DBkZH = LT::DB0C(double(k2/GeV2),muZ2_/GeV2,double(mH2_/GeV2))/GeV2;
  Complex B0ZH = LT::B0C(0.             ,muZ2_/GeV2,double(mH2_/GeV2));
  Complex output = 
    1./6./sw2_/cw2_*
    (BkWW*(18.*sqr(cw2_)+2.*cw2_-0.5)+(4.*cw2_-1.)/3.
     +Complex(DBkWW*((18.*sqr(cw2_)+2.*cw2_-0.5)*k2+2.*muW2_*(12.*sqr(cw2_)+8.*cw2_-5.))))
    +1./12./sw2_/cw2_*
    (-2./3.-BkZH+Complex((2.*mH2_-10.*muZ2_-k2)*DBkZH)
     +Complex(sqr(muZ2_-mH2_)/k2*(BkZH/k2-DBkZH-B0ZH/k2)));
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B4 of Denner
complex<Energy2> OneLoopFFAWZVertex::Sigma_T_W_Boson(Energy2 k2) {
  // bosonic piece
  complex<Energy2> output;
  Complex B0WZ = LT::B0C(0.,muW2_/GeV2,muZ2_/GeV2);
  Complex B0WW = LT::B0C(0.,muW2_/GeV2,muW2_/GeV2);
  Complex B0ZZ = LT::B0C(0.,muZ2_/GeV2,muZ2_/GeV2);
  Complex B0Wl = LT::B0C(0.,muW2_/GeV2,0.);
  Complex B0HH = LT::B0C(0.,double(mH2_ /GeV2),double(mH2_ /GeV2));
  Complex B0WH = LT::B0C(0.             ,muW2_/GeV2,double(mH2_/GeV2));
  if(k2!=ZERO) {
    Complex BkWl = LT::B0C(double(k2/GeV2),muW2_/GeV2,0.);
    Complex BkWH = LT::B0C(double(k2/GeV2),muW2_/GeV2,double(mH2_/GeV2));
    Complex BkWZ = LT::B0C(double(k2/GeV2),muW2_/GeV2,muZ2_/GeV2);
    output = 
      2./3.*((2.*muW2_+5.*k2)*BkWl-2.*muW2_*B0WW
	     -sqr(muW2_)/k2*(BkWl-B0Wl)+k2/3.)
      +1./12./sw2_*
      (((40.*cw2_-1.)*k2+2.*(8.*cw2_+27.-5./cw2_)*muW2_)*BkWZ
       -2.*(8.*cw2_+1.)*(muW2_*B0WW+muZ2_*B0ZZ)+2.*(4.*cw2_-1.)*k2/3.
       -(8.*cw2_+1.)*sqr(muW2_-muZ2_)/k2*(BkWZ-B0WZ)) 
      +1./12./sw2_*
      ((2.*mH2_-10.*muW2_-k2)*BkWH-2.*muW2_*B0WW-2.*mH2_*B0HH-2.*k2/3.
       -sqr(muW2_-mH2_)/k2*(BkWH-B0WH));
  }
  else {
    Complex DB0Wl = LT::DB0C(0.,muW2_/GeV2,0.);
    Complex DB0WH = LT::DB0C(0.,muW2_/GeV2,double(mH2_/GeV2));
    Complex DB0WZ = LT::DB0C(0.,muW2_/GeV2,muZ2_/GeV2);
    output = 
      +4./3.*muW2_*(B0Wl-B0WW)
      -2./3.*sqr(muW2_)*DB0Wl/GeV2
      +1./12./sw2_*
      ((+2.*(8.*cw2_+27.-5./cw2_)*muW2_)*B0WZ
       -2.*(8.*cw2_+1.)*(muW2_*B0WW+muZ2_*B0ZZ)
       -(8.*cw2_+1.)*sqr(muW2_-muZ2_)*DB0WZ/GeV2) 
      +1./12./sw2_*
      ((2.*mH2_-10.*muW2_)*B0WH-2.*muW2_*B0WW-2.*mH2_*B0HH
       -sqr(muW2_-mH2_)*DB0WH/GeV2);
  }
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B4 of Denner
Complex OneLoopFFAWZVertex::DSigma_T_W_Boson(Energy2 k2) {
  assert(k2!=ZERO);
  Complex BkWl = LT::B0C(double(k2/GeV2),muW2_/GeV2,0.);
  Complex BkWH = LT::B0C(double(k2/GeV2),muW2_/GeV2,double(mH2_/GeV2));
  Complex BkWZ = LT::B0C(double(k2/GeV2),muW2_/GeV2,muZ2_/GeV2);
  complex<InvEnergy2> DBkWl = LT::DB0C(double(k2/GeV2),muW2_/GeV2,0.)/GeV2;
  complex<InvEnergy2> DBkWH = LT::DB0C(double(k2/GeV2),muW2_/GeV2,double(mH2_/GeV2))/GeV2;
  complex<InvEnergy2> DBkWZ = LT::DB0C(double(k2/GeV2),muW2_/GeV2,muZ2_/GeV2)/GeV2;
  Complex B0Wl = LT::B0C(0.,muW2_/GeV2,0.);
  Complex B0WZ = LT::B0C(0.,muW2_/GeV2,muZ2_/GeV2);
  Complex B0WH = LT::B0C(0.             ,muW2_/GeV2,double(mH2_/GeV2));
  // bosonic piece
  Complex output = 
    2./3.*(5.*BkWl+Complex((2.*muW2_+5.*k2)*DBkWl)
	   +Complex(sqr(muW2_)/k2*(BkWl/k2-DBkWl-B0Wl/k2))+1./3.)
    +1./12./sw2_*
    (Complex(((40.*cw2_-1.)*k2+2.*(8.*cw2_+27.-5./cw2_)*muW2_)*DBkWZ)
     +(40.*cw2_-1.)*BkWZ+2.*(4.*cw2_-1.)/3.
     +Complex((8.*cw2_+1.)*sqr(muW2_-muZ2_)/k2*(BkWZ/k2-DBkWZ-B0WZ/k2)))
    +1./12./sw2_*
    (-BkWH+Complex((2.*mH2_-10.*muW2_-k2)*DBkWH)-2./3.
     +Complex(sqr(muW2_-mH2_)/k2*(BkWH/k2-DBkWH-B0WH/k2)));
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B1 of Denner (fermionic piece)
complex<Energy2> OneLoopFFAWZVertex::Sigma_T_AA_Fermion(Energy2 k2) {
  complex<Energy2> output(ZERO);
  for(unsigned int ix=1;ix<17;++ix) {
    if(ix==7) ix=11;
    if(ix>11&&ix%2==0) continue;
    Energy2 mf2 = sqr(fermionMasses_[ix]);
    double Nc = ix<7 ? 3. : 1.;
    Complex B0 = LT::B0(0.             ,mf2/GeV2,mf2/GeV2);
    Complex Bk = LT::B0(double(k2/GeV2),mf2/GeV2,mf2/GeV2);
    output += Nc*sqr(ef_[ix])*(-(k2+2.*mf2)*Bk+2.*mf2*B0+k2/3.);
  }
  // multiply by prefactor and return
  return -alphaEW_/3./Constants::pi*output;
}

// eqn B1 of Denner (fermionic piece) differentiated w.r.t. k2
Complex OneLoopFFAWZVertex::DSigma_T_AA_Fermion(Energy2 k2) {
  Complex output(0.);
  for(unsigned int ix=1;ix<17;++ix) {
    if(ix==7) ix=11;
    if(ix>11&&ix%2==0) continue;
    Energy2 mf2 = sqr(fermionMasses_[ix]);
    double Nc = ix<7 ? 3. : 1.;
    Complex              Bk =  LT::B0(double(k2/GeV2),mf2/GeV2,mf2/GeV2);
    complex<InvEnergy2> DBk = LT::DB0(double(k2/GeV2),mf2/GeV2,mf2/GeV2)/GeV2;
    output += Nc*sqr(ef_[ix])*(-Bk-Complex((k2+2.*mf2)*DBk)+1./3.);
  }
  // multiply by prefactor and return
  return -alphaEW_/3./Constants::pi*output;
}

// eqn B2 (fermionic piece) of Denner
complex<Energy2> OneLoopFFAWZVertex::Sigma_T_AZ_Fermion(Energy2 k2) {
  complex<Energy2> output(ZERO);
  for(unsigned int ix=1;ix<17;++ix) {
    if(ix==7) ix=11;
    if(ix>11&&ix%2==0) continue;
    Energy2 mf2 = sqr(fermionMasses_[ix]);
    double Nc = ix<7 ? 3. : 1.;
    Complex B0 = LT::B0(0.             ,mf2/GeV2,mf2/GeV2);
    Complex Bk = LT::B0(double(k2/GeV2),mf2/GeV2,mf2/GeV2);
    output += -2./3.*Nc*ef_[ix]*(gl_[ix]+gr_[ix])*
      (-(k2+2.*mf2)*Bk+2.*mf2*B0+k2/3.);
  }
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B2 (fermionic piece) of Denner differentiated w.r.t. k2
Complex OneLoopFFAWZVertex::DSigma_T_AZ_Fermion(Energy2 k2) {
  Complex output(0.);
  for(unsigned int ix=1;ix<17;++ix) {
    if(ix==7) ix=11;
    if(ix>11&&ix%2==0) continue;
    Energy2 mf2 = sqr(fermionMasses_[ix]);
    double Nc = ix<7 ? 3. : 1.;
    Complex              Bk =  LT::B0(double(k2/GeV2),mf2/GeV2,mf2/GeV2);
    complex<InvEnergy2> DBk = LT::DB0(double(k2/GeV2),mf2/GeV2,mf2/GeV2)/GeV2;
    output += -2./3.*Nc*ef_[ix]*(gl_[ix]+gr_[ix])*
      (-Bk-Complex((k2+2.*mf2)*DBk)+1./3.);
  }
  // multiply by prefactor and return
 return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B3 (fermionic piece) of Denner
complex<Energy2> OneLoopFFAWZVertex::Sigma_T_ZZ_Fermion(Energy2 k2) {
  complex<Energy2> output(ZERO);
  for(unsigned int ix=1;ix<17;++ix) {
    if(ix==7) ix=11;
    Energy2 mf2 = sqr(fermionMasses_[ix]);
    double Nc = ix<7 ? 3. : 1.;
    Complex B0 = LT::B0(0.             ,mf2/GeV2,mf2/GeV2);
    Complex Bk = LT::B0(double(k2/GeV2),mf2/GeV2,mf2/GeV2);
    output += 2./3.*Nc*
      ((sqr(gl_[ix])+sqr(gr_[ix]))*(-(k2+2.*mf2)*Bk+2.*mf2*B0+k2/3.)
       +0.75/sw2_/cw2_*mf2*Bk);
  }
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B3 (fermionic piece) of Denner differentiated w.r.t. k2
Complex OneLoopFFAWZVertex::DSigma_T_ZZ_Fermion(Energy2 k2) {
  Complex output(0.);
  for(unsigned int ix=1;ix<17;++ix) {
    if(ix==7) ix=11;
    Energy2 mf2 = sqr(fermionMasses_[ix]);
    double Nc = ix<7 ? 3. : 1.;
    Complex              Bk =  LT::B0(double(k2/GeV2),mf2/GeV2,mf2/GeV2);
    complex<InvEnergy2> DBk = LT::DB0(double(k2/GeV2),mf2/GeV2,mf2/GeV2)/GeV2;
    output += 2./3.*Nc*
      ((sqr(gl_[ix])+sqr(gr_[ix]))*(-Bk-Complex((k2+2.*mf2)*DBk)+1./3.)
       +0.75/sw2_/cw2_*Complex(mf2*DBk));
  }
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B4 (fermionic piece) of Denner
complex<Energy2> OneLoopFFAWZVertex::Sigma_T_W_Fermion(Energy2 k2) {
  complex<Energy2> output(ZERO);
  // leptons
  for(unsigned int ix=11;ix<17;ix+=2) {
    Energy2 ml2 = sqr(fermionMasses_[ix]);
    Complex B0ll = LT::B0(0.             ,ml2/GeV2,ml2/GeV2);
    Complex B00l = LT::B0(0.             ,0.      ,ml2/GeV2);
    if(k2!=ZERO) {
      Complex Bk0l = LT::B0(double(k2/GeV2),0.      ,ml2/GeV2);
      output += 1./3./sw2_*(-(k2-0.5*ml2)*Bk0l+k2/3.
			    +ml2*B0ll+0.5*sqr(ml2)/k2*(Bk0l-B00l));
    }
    else {
      complex<InvEnergy2> DBk0l = LT::DB0(0.,0.      ,ml2/GeV2)/GeV2;
      output += 1./3./sw2_*(0.5*ml2*B00l
			    +ml2*B0ll+0.5*sqr(ml2)*DBk0l);
    }
  }
  // quarks (NO CKM mixing included)
  for(unsigned int ix=1;ix<7;ix+=2) {
    Energy2 md2 = sqr(fermionMasses_[ix  ]);
    Energy2 mu2 = sqr(fermionMasses_[ix+1]);
    Complex B0uu = LT::B0(0.             ,mu2/GeV2,mu2/GeV2);
    Complex B0dd = LT::B0(0.             ,md2/GeV2,md2/GeV2);
    Complex B0ud = LT::B0(0.             ,mu2/GeV2,md2/GeV2);
    if(k2!=ZERO) {
      Complex Bkud = LT::B0(double(k2/GeV2),mu2/GeV2,md2/GeV2);
      output += 1./sw2_*(-(k2-0.5*(mu2+md2))*Bkud+k2/3.
			 +md2*B0dd+mu2*B0uu+0.5*sqr(mu2-md2)/k2*(Bkud-B0ud));
    }
    else {
      complex<InvEnergy2> DBkud = LT::DB0(0.,mu2/GeV2,md2/GeV2)/GeV2;
      output += 1./sw2_*(+0.5*(mu2+md2)*B0ud
			 +md2*B0dd+mu2*B0uu+0.5*sqr(mu2-md2)*DBkud);
    }
  }
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

// eqn B4 of Denner
Complex OneLoopFFAWZVertex::DSigma_T_W_Fermion(Energy2 k2) {
  assert(k2!=ZERO);
  Complex output(ZERO);
  // leptons
  for(unsigned int ix=11;ix<17;ix+=2) {
    Energy2 ml2 = sqr(fermionMasses_[ix]);
    Complex              Bk0l =  LT::B0(double(k2/GeV2),0.,ml2/GeV2);
    complex<InvEnergy2> DBk0l = LT::DB0(double(k2/GeV2),0.,ml2/GeV2)/GeV2;
    Complex B00l = LT::B0(0.             ,0.      ,ml2/GeV2);
    output += 1./3./sw2_*(-Bk0l-(k2-0.5*ml2)*DBk0l+1./3.
			  -0.5*sqr(ml2)/k2*(Bk0l/k2-DBk0l-B00l/k2));
  }
  // quarks (NO CKM mixing included)
  for(unsigned int ix=1;ix<7;ix+=2) {
    Energy2 md2 = sqr(fermionMasses_[ix  ]);
    Energy2 mu2 = sqr(fermionMasses_[ix+1]);
    Complex              Bkud =  LT::B0(double(k2/GeV2),mu2/GeV2,md2/GeV2);
    complex<InvEnergy2> DBkud = LT::DB0(double(k2/GeV2),mu2/GeV2,md2/GeV2)/GeV2;
    Complex B0ud = LT::B0(0.             ,mu2/GeV2,md2/GeV2);
    output += 1./sw2_*(-Bkud-(k2-0.5*(mu2+md2))*DBkud+1./3.
		       -0.5*sqr(mu2-md2)/k2*(Bkud/k2-DBkud-B0ud/k2));
  }
  // multiply by prefactor and return
  return -0.25*alphaEW_/Constants::pi*output;
}

pair<Complex,Complex> OneLoopFFAWZVertex::
renormalisedVertex(Energy2 sHat,tcPDPtr aa,tcPDPtr,tcPDPtr cc) {
  // get the gauge boson and the fermion
  int iferm = abs(aa->id());
  assert((iferm>=1&&iferm<=5) || (iferm>=11&&iferm<=16));
  Complex Qf  = ef_[iferm];
  double I3f = iferm % 2 == 0 ? 0.5 : -0.5;
  int ibos  = cc->id();
  Complex B00Z = LT::B0C(0.       ,0.,muZ2_/GeV2);
  Complex B00W = LT::B0C(0.       ,0.,muW2_/GeV2);
  Complex BS00 = LT::B0C(double(sHat/GeV2),0.,0.        );
  Complex BSWW = LT::B0C(double(sHat/GeV2),muW2_/GeV2,muW2_/GeV2);
  complex<InvEnergy2> C0Z0 = LT::C0C(0.,0.,double(sHat/GeV2),0.,muZ2_/GeV2,0.)/GeV2;
  complex<InvEnergy2> C0W0 = LT::C0C(0.,0.,double(sHat/GeV2),0.,muW2_/GeV2,0.)/GeV2;
  complex<InvEnergy2> CW0W = LT::C0C(0.,0.,double(sHat/GeV2),
				     muW2_/GeV2,0.,muW2_/GeV2)/GeV2;
  // fermion field renormalisation is common
  Complex B1Z = LT::B1C(0.,0.,muZ2_/GeV2);
  Complex dZ_F_L;
  if(iferm!=5) {
    Complex B1W = LT::B1C(0.,0.,muW2_/GeV2);
    dZ_F_L = 0.25*alphaEW_/Constants::pi*(sqr(gl_[iferm])*(2.*B1Z+1.)
					  +0.5/sw2_*( 2.           *B1W+1.));
  }
  else {
    Energy2 mt2=sqr(fermionMasses_[6]);
    Complex B1W = LT::B1C(0.,double(mt2/GeV2),muW2_/GeV2);
    dZ_F_L = 0.25*alphaEW_/Constants::pi*(sqr(gl_[iferm])*(2.*B1Z+1.)
					  +0.5/sw2_*((2.+mt2/muW2_)*B1W+1.));
  }
  Complex dZ_F_R = 0.25*alphaEW_/Constants::pi* sqr(gr_[iferm])*(2.*B1Z+1.);
  // vertex correction
  pair<Complex,Complex> output = make_pair(0.,0.);
  // photon
  if(ibos==ParticleID::gamma) {
    // first the unrenormalised vertex
    output.second = -0.25*alphaEW_*sqr(ef_[iferm])*sw2_/cw2_/Constants::pi*
      (2.-2./sHat*(muZ2_+2.*sHat)*B00Z+(3.*sHat+2.*muZ2_)/sHat*BS00
       +2.*Complex(sqr(muZ2_+sHat)/sHat*C0Z0));
    output.first = 0.125*alphaEW_/Constants::pi/sHat/sw2_*
      (1./Qf*(2.*Qf*(2.*sHat+muW2_)*B00W
	      +(2.*I3f-Qf)*(2.*sHat+(3.*sHat+2.*muW2_)*BS00
			    +2.*sqr(sHat+muW2_)*C0W0)
	      -2.*I3f*((sHat+2.*muW2_)*BSWW-2.*muW2_*(2.*sHat+muW2_)*CW0W))
       +sqr(I3f-Qf*sw2_)/cw2_*
       (-4.*sHat+4.*(2.*sHat+muZ2_)*B00Z-2.*(3.*sHat+2.*muZ2_)*BS00
	-4.*sqr(sHat+muZ2_)*C0Z0));
    if(iferm==5) {
      Energy2 mt2 = sqr(fermionMasses_[6]);
      Complex B0tW = LT::B0C(0.       ,double(mt2/GeV2),muW2_/GeV2);
      Complex BStt = LT::B0C(double(sHat/GeV2),double(mt2/GeV2),double(mt2/GeV2));
      complex<InvEnergy2> CtWt = LT::C0C(0.,0.,double(sHat/GeV2),
					 double(mt2/GeV2),muW2_/GeV2,double(mt2/GeV2))/GeV2;
      complex<InvEnergy2> CWtW = LT::C0C(0.,0.,double(sHat/GeV2),
					 muW2_/GeV2,double(mt2/GeV2),muW2_/GeV2)/GeV2;
      output.first += 0.125*alphaEW_/Constants::pi/sHat/sw2_*
	(-2.*((2.*sHat+muW2_)*(B00W-B0tW)+(3.*sHat+2.*muW2_)*(BS00-BStt)
	      +2.*sqr(sHat+muW2_)*(C0W0-CtWt)
	      +3.*muW2_*(2.*sHat+muW2_)*(CW0W-CWtW))
	 -mt2/muW2_*((mt2+muW2_)*(B0tW+2.*BStt-3.*BSWW)
		     +sHat*(BStt-1.5*BSWW-2.5)
		     +(muW2_*(2.*sHat+3.*muW2_)-mt2*(mt2+sHat))*(2.*CtWt+3.*CWtW)));
    }
    // counter terms
    output.first  += de_ + 0.5*dZ_AA_ + dZ_F_L - 0.5*gl_[iferm]/Qf*dZ_ZA_;
    output.second += de_ + 0.5*dZ_AA_ + dZ_F_R - 0.5*gr_[iferm]/Qf*dZ_ZA_;
  }
  // Z0
  else if(ibos==ParticleID::Z0) {    
    output.second = -0.25*alphaEW_*sqr(ef_[iferm])*sw2_/cw2_/Constants::pi*
      (2.-2./sHat*(muZ2_+2.*sHat)*B00Z+(3.*sHat+2.*muZ2_)/sHat*BS00
       +2.*Complex(sqr(muZ2_+sHat)/sHat*C0Z0));
    output.first = 0.125*alphaEW_/Constants::pi/sHat/sw2_*
      (2.*(2.*sHat+muW2_)*B00W
       +1./(I3f-Qf*sw2_)*((I3f*cw2_-I3f*sw2_+Qf*sw2_)*(2.*sHat+(3.*sHat+2.*muW2_)*BS00
						       +2.*sqr(sHat+muW2_)*C0W0)
			  -2.*cw2_*I3f*((sHat+2.*muW2_)*BSWW-2.*muW2_*(2.*sHat+muW2_)*CW0W))   
       +sqr(I3f-Qf*sw2_)/cw2_*(-4.*sHat+4.*(2.*sHat+muZ2_)*B00Z-2.*(3.*sHat+2.*muZ2_)*BS00
			       -4.*sqr(sHat+muZ2_)*C0Z0));
    if(iferm==5) {
      Energy2 mt2 = sqr(fermionMasses_[6]);
      Complex B0tW = LT::B0C(0.       ,double(mt2/GeV2),muW2_/GeV2);
      Complex BStt = LT::B0C(double(sHat/GeV2),double(mt2/GeV2),double(mt2/GeV2));
      complex<InvEnergy2> CtWt = LT::C0C(0.,0.,double(sHat/GeV2),
					 double(mt2/GeV2),muW2_/GeV2,double(mt2/GeV2))/GeV2;
      complex<InvEnergy2> CWtW = LT::C0C(0.,0.,double(sHat/GeV2),
					 muW2_/GeV2,double(mt2/GeV2),muW2_/GeV2)/GeV2;
      output.first -= 0.125*alphaEW_/Constants::pi/sHat/sw2_/(2.*sw2_-3.)*
	(2.*(2.*sHat+muW2_)*(2.*sw2_-3.)*(B00W-B0tW)
	 +(3.*sHat+2.*muW2_)*(4.*sw2_-3.)*(BS00-BStt)
	 -12.*cw2_*muW2_*(2.*sHat+muW2_)*(CW0W-CWtW)
	 +2.*(4.*sw2_-3.)*sqr(sHat+muW2_)*(C0W0-CtWt)
	 +mt2/muW2_*((mt2+muW2_)*(4.*BStt*sw2_-3.*BSWW*(sw2_-cw2_)+B0tW*(2.*sw2_-3.))
		     +2.*(sHat*sw2_-3.*muW2_)*BStt+sHat*(1.5-5.*sw2_)
		     +(6.*muW2_-1.5*sHat*(sw2_-cw2_))*BSWW-12.*sHat*muW2_*CtWt
		     -3.*(4.*sqr(muW2_)-2.*mt2*muW2_-mt2*sHat)*(CtWt+CWtW)
		     -3.*(3.*sqr(muW2_)-sqr(mt2))*CWtW
		     +2.*sw2_*(muW2_*(2.*sHat+3.*muW2_)-mt2*(mt2+sHat))*
		     (2.*CtWt+3.*CWtW)));
    }
    // counter terms
    output.first  += de_ + dsw_/cw2_ + dZ_F_L + 0.5*dZ_ZZ_ - 0.5*Qf*dZ_AZ_/gl_[iferm] 
      - 2.*I3f*dsw_/sw_/cw_/gl_[iferm] ;
    output.second += de_ + dsw_/cw2_ + dZ_F_R + 0.5*dZ_ZZ_ + 0.5*dZ_AZ_/sw_*cw_;
  }
  else
    assert(false);
  return output;
}

Complex OneLoopFFAWZVertex::deltaAlphaMZ() {
  Complex Pi0(0.),PiMZ(0.);
  for(unsigned int ix=1;ix<17;++ix) {
    if(ix==6) ix=11;
    if(ix>11&&ix%2==0) continue;
    Energy2 mf2 = sqr(fermionMasses_[ix]);
    double Nc = ix<7 ? 3. : 1.;
    // Pi(0)
    Complex              B0 =  LT::B0(0.,mf2/GeV2,mf2/GeV2);
    complex<InvEnergy2> DB0 = LT::DB0(0.,mf2/GeV2,mf2/GeV2)/GeV2;
    Pi0 += Nc*sqr(ef_[ix])*(-B0-Complex(2.*mf2*DB0)+1./3.);
    // Pi(MZ)
    Complex Bk = LT::B0(sqr(mZ_/GeV),mf2/GeV2,mf2/GeV2);
    Complex mu2 = double(mf2/sqr(mZ_));
    PiMZ += Nc*sqr(ef_[ix])*(-(1.+2.*mu2)*Bk+2.*mu2*B0+1./3.);
  }
  // multiply by prefactor and return
  return -alphaEW_/3./Constants::pi*(Pi0-real(PiMZ));
}

VectorWaveFunction 
OneLoopFFAWZVertex::selfEnergyCorrection(tcPDPtr particle,
					 const VectorWaveFunction & old) {
  assert(old.direction()==intermediate);
  long oldId = old.particle()->id();
  long newId = particle->id();
  Energy2 scale = old.m2();
  Complex fact(0.);
  if(abs(oldId)==ParticleID::Wplus) {
    assert(oldId==newId);
    fact = -SigmaHat_T_W(old.m2())/(old.m2()-muW2_);
    particle = particle->CC();
  }
  else if(oldId==ParticleID::gamma&&
	  newId==ParticleID::gamma) {
    fact = -SigmaHat_T_AA(old.m2())/old.m2();
  }
  else if(oldId==ParticleID::gamma&&
	  newId==ParticleID::Z0) {
    fact = -SigmaHat_T_AZ(old.m2())/(old.m2()-muZ2_);
  }
  else if(oldId==ParticleID::Z0&&
	  newId==ParticleID::gamma) {
    fact = -SigmaHat_T_AZ(old.m2())/old.m2();
  }
  else if(oldId==ParticleID::Z0&&
	  newId==ParticleID::Z0) {
    fact = -SigmaHat_T_ZZ(old.m2())/(old.m2()-muZ2_);
  }
  else 
    assert(false);
  return VectorWaveFunction(old.momentum(),particle,
			    fact*old.wave(),old.direction());
}

complex<InvEnergy2> 
OneLoopFFAWZVertex::neutralCurrentBT(int hel1, int hel2,
				     complex<Energy2> mV2,
				     complex<Energy2> mVp2, 
				     Energy2 mQ2, Energy2 mq2, Energy2 ml2,
				     Energy2 sHat, Energy2 tHat, Energy2 uHat) {
  assert(abs(hel1)==1&&abs(hel2)==1);
  if(hel1==-hel2) {
    complex<InvEnergy2> Cl = LT::C0C(double(ml2/GeV2),double(ml2/GeV2),double(sHat/GeV2),
				     mV2/GeV2,double(ml2/GeV2),mVp2/GeV2)/GeV2;
    complex<InvEnergy2> Cq = LT::C0C(double(mq2/GeV2),double(mq2/GeV2),double(sHat/GeV2),
				     mV2/GeV2,double(mQ2/GeV2),mVp2/GeV2)/GeV2;
    complex<InvEnergy4> D = LT::D0C(double(mq2/GeV2),double(mq2/GeV2),
				    double(ml2/GeV2),double(ml2/GeV2),
				    double(sHat/GeV2),double(tHat/GeV2),
				    mV2/GeV2,double(mQ2/GeV2),mVp2/GeV2,
				    double(ml2/GeV2))/GeV2/GeV2;
    return -2.*(Cl+Cq-(tHat-mQ2)*D);
  }
  else {
    return 1./sqr(uHat)*
      (2.*uHat*(LT::B0C(double(sHat/GeV2),       mV2/GeV2 ,mVp2/GeV2)-
		LT::B0C(double(tHat/GeV2),double(mQ2/GeV2),0.))
       -tHat*(mQ2-mV2-mVp2-tHat+uHat)*
       (LT::C0C(double(ml2/GeV2),double(mq2/GeV2),double(tHat/GeV2),
		0.,mV2 /GeV2,double(mQ2/GeV2))+
	LT::C0C(double(ml2/GeV2),double(mq2/GeV2),double(tHat/GeV2),
		0.,mVp2/GeV2,double(mQ2/GeV2)))/GeV2
       -(sqr(tHat)+sqr(uHat)+sHat*(mQ2-mV2-mVp2))*
       (LT::C0C(double(ml2/GeV2),double(ml2/GeV2),double(sHat/GeV2),
		mV2/GeV2,double(ml2/GeV2),mVp2/GeV2)+
	LT::C0C(double(mq2/GeV2),double(mq2/GeV2),double(sHat/GeV2),
		mV2/GeV2,double(mQ2/GeV2),mVp2/GeV2))/GeV2
       +(tHat*sqr(mV2+mVp2-mQ2-2.*sHat)
	 +(uHat*(2.*uHat-mQ2)-2.*sqr(sHat))*(mV2+mVp2-mQ2-2.*sHat)
	 -2.*uHat*(complex<Energy4>(sqr(uHat))-mV2*mVp2)
	 + uHat*sHat*(complex<Energy2>(sHat)-mQ2)
	 - sHat*sqr(sHat) )*
       LT::D0C(double(mq2/GeV2),double(mq2/GeV2),
	       double(ml2/GeV2),double(ml2/GeV2),
	       double(sHat/GeV2),double(tHat/GeV2),
	       mV2/GeV2,double(mQ2/GeV2),mVp2/GeV2,double(ml2/GeV2))/GeV2/GeV2
       );
  }
}

vector<vector<complex<InvEnergy2> > > 
OneLoopFFAWZVertex::neutralCurrentFBox(tcPDPtr q1, tcPDPtr q2,
				       tcPDPtr l1, tcPDPtr l2,
				       Energy2 sHat, Energy2 tHat, Energy2 uHat) {
  assert(q1->id()==-q2->id());
  int iq = abs(q1->id());
  assert(iq<=5);
  assert(l1->id()==-l2->id());
  int il = abs(l1->id());
  assert(il>=11 && il<=16);
vector<vector<complex<InvEnergy2> > > output(2,vector<complex<InvEnergy2> >());
  // B functions for the ZZ boxes 
  complex<InvEnergy2> btZs = neutralCurrentBT(1, 1,muZ2_,muZ2_,ZERO,ZERO,ZERO,
					      sHat,tHat,uHat);
  complex<InvEnergy2> btZo = neutralCurrentBT(1,-1,muZ2_,muZ2_,ZERO,ZERO,ZERO,
					      sHat,tHat,uHat);
  complex<InvEnergy2> buZs = neutralCurrentBU(1, 1,muZ2_,muZ2_,ZERO,ZERO,ZERO,
					      sHat,tHat,uHat);
  complex<InvEnergy2> buZo = neutralCurrentBU(1,-1,muZ2_,muZ2_,ZERO,ZERO,ZERO,
					      sHat,tHat,uHat);
  // and for WW
  complex<InvEnergy2> bW;
  if(iq%2==0)
    bW = neutralCurrentBU(1,1,muW2_,muW2_,ZERO,ZERO,ZERO,sHat,tHat,uHat);
  else if(iq!=5)
    bW = neutralCurrentBT(1,1,muW2_,muW2_,ZERO,ZERO,ZERO,sHat,tHat,uHat);
  else
    bW = neutralCurrentBT(1,1,muW2_,muW2_,sqr(fermionMasses_[6]),ZERO,ZERO,
			  sHat,tHat,uHat);
  // loop to calculate the coefficients 
  for(int sigma=0;sigma<2;++sigma) {
    Complex gq = sigma==0 ? gl_[iq] : gr_[iq];
    for(int tau=0;tau<2;++tau) {
      Complex gl = tau==0 ? gl_[il] : gr_[il];
      // Z component
      complex<InvEnergy2> value = sigma==tau ? btZs + buZs : btZo + buZo;   
      value *= sqr(alphaEW_*gq*gl);
      // W component
      if(sigma==0&&tau==0)
	value += 0.25*sqr(alphaEW_/sw2_)*bW;
      // store the total
      output[sigma].push_back(value);
    }
  }
  return output;
}

void OneLoopFFAWZVertex::neutralCurrentME(tcPDPtr q1, tcPDPtr q2,
					  tcPDPtr l1, tcPDPtr l2,
					  Energy2 sHat, Energy2 tHat, Energy2 uHat) {
  assert(q1->id()==-q2->id());
  int iq = abs(q1->id());
  assert(iq<=5);
  assert(l1->id()==-l2->id());
  int il = abs(l1->id());
  assert(il>=11 && il<=16);
  // vertex correction factors
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0    = getParticleData(ParticleID::Z0);
  pair<Complex,Complex> dlGamma = renormalisedVertex(sHat,l1,l2,gamma);
  pair<Complex,Complex> dqGamma = renormalisedVertex(sHat,q1,q2,gamma);
  pair<Complex,Complex> dlZ0    = renormalisedVertex(sHat,l1,l2,Z0   );
  pair<Complex,Complex> dqZ0    = renormalisedVertex(sHat,q1,q2,Z0   );
  Complex dlg[2] = {dlGamma.first,dlGamma.second};
  Complex dlZ[2] = {dlZ0   .first,dlZ0   .second};
  Complex dqg[2] = {dqGamma.first,dqGamma.second};
  Complex dqZ[2] = {dqZ0   .first,dqZ0   .second};
  // box pieces
  vector<vector<complex<InvEnergy2> > > boxCoeff = 
    neutralCurrentFBox(q1,q2,l1,l2,sHat,tHat,uHat);
  // full result
  Complex loSum(0.),EWSum(0.);
  for(int sigma=0;sigma<2;++sigma) {
    Complex gq = sigma==0 ? gl_[iq] : gr_[iq];
    for(int tau=0;tau<2;++tau) {
      Complex gl = tau==0 ? gl_[il] : gr_[il];
      // helicity amplitude
      Energy2 amp = sigma==tau ? 2.*uHat : 2.*tHat; 
      Complex lo = ef_[iq]*ef_[il]+Complex(sHat/(sHat-muZ2_))*gq*gl;
      // LO piece
      Complex loAmp = -sqr(e_)*lo*amp/sHat;
      // self energy piece
      Complex self = sqr(e_)*amp*(ef_[iq]*ef_[il]        /sqr(  sHat      )*SigmaHat_T_AA(sHat)+
				  gq     *gl             /sqr(  sHat-muZ2_)*SigmaHat_T_ZZ(sHat)-
				  (ef_[il]*gq+ef_[iq]*gl)/sHat/(sHat-muZ2_)*SigmaHat_T_AZ(sHat));
      // vertex piece
      Complex vertex = -sqr(e_)*amp*(ef_[iq]*ef_[il]*(dqg[sigma]+dlg[tau])/sHat +
				     gq     *gl     *(dqZ[sigma]+dlZ[tau])/(sHat-muZ2_));
      // box
      Complex box = boxCoeff[sigma][tau]*amp;
      Complex loop = self+vertex+box;
      EWSum += loAmp*conj(loop)+loop*conj(loAmp);
      loSum += std::norm(loAmp);
    }
  }
  cerr << "TEST OF EW CORRECTION LO: " << loSum/12. << " NLO :" << EWSum/12. << "\n";
  if(abs(EWSum/loSum)>0.01) cerr << "testing ratio " << EWSum/loSum << "\n";
}

void OneLoopFFAWZVertex::clearCache() {
  LT::clearcache();
}
