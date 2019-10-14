// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHZPVertex class.
//

#include "SMHZPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SMHZPVertex::SMHZPVertex()
  :_couplast(0.),_q2last(),_mw(),_mz(),_massopt(1),
   _minloop(6),_maxloop(6) {
  orderInGs(0);
  orderInGem(3);
  kinematics(true);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr SMHZPVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SMHZPVertex::fullclone() const {
  return new_ptr(*this);
}

void SMHZPVertex::doinit() {
  GeneralVVSVertex::doinit();
  //PDG codes for particles at vertices
  addToList(23,22,25);
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) 
    throw InitException() 
      << "SMHGGVertex::doinit() - The pointer to the SM object is null."
      << Exception::abortnow;
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  _mz = getParticleData(ThePEG::ParticleID::Z0)->mass();
  GeneralVVSVertex::doinit();
}

void SMHZPVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << ounit(_mz,GeV) << _massopt
     << _minloop << _maxloop;
}

void SMHZPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw, GeV) >> iunit(_mz,GeV) >> _massopt
     >> _minloop >> _maxloop;
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<SMHZPVertex,GeneralVVSVertex>
describeHerwigSMHZPVertex("Herwig::SMHZPVertex", "libHerwig.so");

void SMHZPVertex::Init() {

  static ClassDocumentation<SMHZPVertex> documentation
    ("The SMHZPVertex class provides a simple implementation of the "
     "Higgs-Z-Photon loop looping to allow the calculation of the "
     "associated Higgs decay mode H -> Z gamma.");


  static Parameter<SMHZPVertex,int> interfaceMinQuarkInLoop
    ("MinQuarkInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &SMHZPVertex::_minloop, 6, 1, 6,
     false, false, Interface::limited);

  static Parameter<SMHZPVertex,int> interfaceMaxQuarkInLoop
    ("MaxQuarkInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &SMHZPVertex::_maxloop, 6, 1, 6,
     false, false, Interface::limited);

  static Switch<SMHZPVertex,unsigned int> interfaceMassOption
    ("LoopMassScheme",
     "Switch for the treatment of the masses in the loops ",
     &SMHZPVertex::_massopt, 2, false, false);
  static SwitchOption interfaceHeavyMass
    (interfaceMassOption,
     "PoleMasses",
     "The loop is calculcated with the pole quark masses",
     1);
  static SwitchOption interfaceNormalMass
    (interfaceMassOption,
     "RunningMasses",
     "running quark masses are taken in the loop",
     2);
}

void SMHZPVertex::setCoupling(Energy2 q2, tcPDPtr part2,
			      tcPDPtr part3,
#ifndef NDEBUG
			      tcPDPtr part1) {
#else
			      tcPDPtr) {
#endif
  if(part3->id()==ParticleID::Z0) swap(part2,part3);
  assert( part1->id() == ParticleID::h0 &&
 	  part2->id() == ParticleID::Z0 && part3->id() == ParticleID::gamma );
  int Qminloop = _minloop;
  int Qmaxloop = _maxloop;
  if (_maxloop < _minloop) swap(Qmaxloop,Qminloop);
  double cw = sqrt(1.-sin2ThetaW()),sw=sqrt(sin2ThetaW()),tw = sw/cw;
  if(q2 != _q2last||_couplast==0.) {
    double g   = weakCoupling(q2);
    double e2 = sqr(electroMagneticCoupling(q2));
    _couplast = UnitRemoval::E * e2 * g / 16. / _mw/ sqr(Constants::pi);
    _q2last = q2;
  }
  norm(_couplast);
  Complex loop(0.);
  // quark loops
  for ( int i = Qminloop; i <= Qmaxloop; ++i ) {
    tcPDPtr qrk = getParticleData(i);
    Energy mass = (2 == _massopt) ? _theSM->mass(q2,qrk) : qrk->mass();
    double charge = i%2==0 ?
      generator()->standardModel()->eu() : generator()->standardModel()->ed();
    double gv = i%2==0 ?
      generator()->standardModel()->vu() : generator()->standardModel()->vd();
    double tau = 0.25*invariant(0,0)/sqr(mass), lambda(0.25*sqr(_mz/mass));
    loop += 3.*charge*gv *(I1(tau,lambda)-I2(tau,lambda))/(sw*cw);
  }
  // lepton loops
  int Lminloop = 3; // still fixed value
  int Lmaxloop = 3; // still fixed value
  for (int i = Lminloop; i <= Lmaxloop; ++i) {
    tcPDPtr lpt = getParticleData(9 + 2*i);
    Energy mass = (2 == _massopt) ? _theSM->mass(q2,lpt) : lpt->mass();
    double charge = generator()->standardModel()->ee();
    double gv     = generator()->standardModel()->ve();
    double tau = 0.25*invariant(0,0)/sqr(mass), lambda(0.25*sqr(_mz/mass));
    loop += charge*gv*(I1(tau,lambda)-I2(tau,lambda))/(sw*cw);
  }
  // W loop
  double tau = 0.25*invariant(0,0)/sqr(_mw), lambda(0.25*sqr(_mz/_mw));
  loop += ( 4.*(3.-sqr(tw))*I2(tau,lambda) +
	    ((1.+2.*tau)*sqr(tw)-(5.+2.*tau))*I1(tau,lambda))/tw;
  a00(loop);
  a11(0.0);
  a12(0.0);
  a21(-loop);
  a22(0.0);
  aEp(0.0);
  // test of the width calculation
  // Energy mh = getParticleData(25)->mass();
  // Energy pre = sqr(weakCoupling(q2))*pow(electroMagneticCoupling(q2),4)*mh*sqr(mh/_mw)
  //   /128./16./pow(Constants::pi,5)*pow(double(1.-sqr(_mz/mh)),3)*
  //   std::real(std::norm(loop));
}

Complex SMHZPVertex::I1(double tau,double lambda) const {
  return (-0.5+0.5/(tau-lambda)*(f(tau)-f(lambda))+
	  lambda/(tau-lambda)*(g(tau)-g(lambda)))/(tau-lambda);
}

Complex SMHZPVertex::I2(double tau,double lambda) const {
  return 0.5/(tau-lambda)*(f(tau)-f(lambda));
}

Complex SMHZPVertex::f(double tau) const {
  if(tau>0 && tau<= 1.) {
    return sqr(asin(sqrt(tau)));
  }
  else if(tau>1.) {
    double lx = log(sqrt(tau)+sqrt(tau-1));
    return -sqr(lx)+0.25*sqr(Constants::pi)+Complex(0.,1.)*Constants::pi*lx;
  }
  else {
    assert(false);
    return 0.;
  }
}

Complex SMHZPVertex::g(double tau) const {
  if(tau>0 && tau<= 1.) {
    return sqrt((1.-tau)/tau)*asin(sqrt(tau));
  }
  else if(tau>1.) {
    double lx = log(sqrt(tau)+sqrt(tau-1));
    double root = sqrt((tau-1.)/tau);
    return root*(lx-0.5*Complex(0,1)*Constants::pi);
  }
  else {
    assert(false);
    return 0.;
  }
}
