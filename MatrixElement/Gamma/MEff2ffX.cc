// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2ffX class.
//

#include "MEff2ffX.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

void MEff2ffX::getDiagrams() const {
  vector<DiagPtr> diags = amp_->getDiagrams(1);
  for(DiagPtr diag : diags) add(diag);
}

Energy2 MEff2ffX::scale() const {
  return sHat();
}

int MEff2ffX::nDim() const {
  return 4+amp_->nDim(1);
}

void MEff2ffX::setKinematics() {
  HwMEBase::setKinematics(); // Always call the base class method first.
}

bool MEff2ffX::generateKinematics(const double * r) {
  // initialise the jacobian
  jacobian(1.);
  // roots
  Energy rS(sqrt(sHat()));
  Energy W(ZERO);
  Energy m = mePartonData()[0]->mass();
  if(amp_->nDim(1)==0) {
    W = mePartonData()[4]->mass();
  }
  else {
    Energy2 jacW(ZERO);
    Energy Wmax = rS-2.*m;
    W = amp_->generateW(r[4],tcPDVector(mePartonData().begin()+4,mePartonData().end()),mHatMin_,Wmax,jacW,sHat());
    jacobian(jacW*jacobian()/sHat());
  }
  Energy2 m2=sqr(m), W2=sqr(W);
  double beta = sqrt(1.-4.*m2/sHat());
  // initialise the jacobian
  jacobian(jacobian()*0.25*sqr(Constants::pi)/beta);
  // outer integral is t2
  Energy2 deltat2 = beta*sqrt((sHat()-W2)*(sHat()-sqr(W+2.*m)));
  Energy2 t2min = -0.5*(sHat()-W2-2.*m*W-4.*m2+deltat2);
  Energy2 t2max = sqr(m*W*(W+2.*m))/t2min/sHat();
  Energy2 thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[2]);
  if(thmin>abs(t2max)) t2max=-thmin;
  if(abs(t2max)<Q2_2min_) t2max= -Q2_2min_;
  if(abs(t2min)>Q2_2max_) t2min= -Q2_2max_;
  if(t2max<t2min) return false;
  t2_ = t2min*pow(t2max/t2min,r[0]);
  double y2 = sqrt(1.-4.*m2/t2_);
  jacobian(jacobian()*log(t2max/t2min)*t2_/sHat());
  // then t1
  Energy2 Q  = (t2_*sHat() - m2*t2_ - sqr(m*W) -sHat()*t2_*beta*y2)/m2;
  Energy2 a1 = 2.*(Q+t2_+2.*m2+W2);
  Energy4 b1 = sqr(Q)-8.*m2*t2_-sqr(t2_)-8.*sqr(m*W)+2.*t2_*W2-pow<4,1>(W);
  Energy6 c1 = 4.*sqr(m*(W2-t2_));
  Energy8 delta1 = (Q+t2_-4.*m*W-W2)*(Q+t2_+4.*m*W-W2)*(sqr(Q)-16.*m2*t2_-2.*Q*t2_+sqr(t2_)+2.*Q*W2-2.*t2_*W2+pow<4,1>(W));
  Energy2 t1min = -0.5*(b1+sqrt(delta1))/a1;
  Energy2 t1max = c1/a1/t1min;
  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[3]);
  if(thmin>abs(t1max)) t1max=-thmin;
  if(abs(t1max)<Q2_1min_) t1max= -Q2_1min_;
  if(abs(t1min)>Q2_1max_) t1min= -Q2_1max_;
  if(t1max<t1min) return false;
  t1_ = t1min*pow(t1max/t1min,r[1]);
  double y1 = sqrt(1.-4.*m2/t1_);
  jacobian(jacobian()*log(t1max/t1min)*t1_/sHat());
  // then s1
  Energy2 nu = 0.5*(W2-t1_-t2_);
  Energy  K  = sqrt(sqr(nu)-t1_*t2_)/W;
  double d1 = log(sHat()/(nu+K*W)*sqr(1.+beta)/((1.+y1)*(1.+y2)));
  jacobian(jacobian()*d1);
  Energy2 X1 = (nu+K*W)*(1.+y1)*exp(d1*r[2]);
  Energy2 s1 = 0.5*X1+m2+t2_+2.*m2*t2_/X1;
  // finally s2
  Energy6 G3 = pow<6,1>(m) - 2.*sqr(m2)*s1 + sHat()*t2_*(sHat()-s1+t2_) + m2*(sqr(s1)-3.*sHat()*t2_);
  Energy6 G4 = sqr(m2)*t1_ + t1_*(sqr(s1)+t2_*W2+s1*(t1_-t2_-W2))
    - m2*(2.*s1*t1_ - sqr(t2_-W2) + t1_*(t2_ + W2));
  auto Delta = 16.*G3*G4;
  Energy6 rDelta = sqrt(Delta);
  Energy4 a = sqr(s1-t2_-m2)-4.*m2*t2_;
  Energy6 b = -2.*pow<6,1>(m) + 4.*sqr(m2)*s1 - 2.*m2*sqr(s1) + 2*sqr(m2)*t1_ - 2.*m2*sHat()*t1_
    + 2.*sHat()*s1*t1_ - 2*sqr(s1)*t1_ + 8.*sqr(m2)*t2_ - 2*m2*sHat()*t2_ + 2*sHat()*s1*t2_ - 2*m2*t1_*t2_ + 2*sHat()*t1_*t2_
    + 2*s1*t1_*t2_ - 2*m2*sqr(t2_) - 2*sHat()*sqr(t2_) - 4*sqr(m2)*W2 + 2*m2*sHat()*W2 + 4*m2*s1*W2
    - 2*sHat()*s1*W2 + 2*sHat()*t2_*W2;
  Energy2 s2 = -0.5/a*(b+rDelta*cos(Constants::pi*r[3]));
  // now generate the momenta
  // Energies and maginitudes of the momenta
  Energy E1 = 0.5*(sHat()+m2-s2)/rS; Energy P1 = sqrt(sqr(E1)-m2);
  Energy E2 = 0.5*(sHat()+m2-s1)/rS; Energy P2 = sqrt(sqr(E2)-m2);
  Energy EX = 0.5*(s1+s2-2*m2)/rS  ; Energy PX = sqrt(sqr(EX)-W2);
  // angle wrt beams
  // e-
  Energy6 D1 = 0.25*(-pow<6,1>(m)+ 2*sqr(m2)*s2-sHat()*t1_*(sHat()-s2+t1_)-sqr(m)*(sqr(s2) - 3.*sHat()*t1_));
  double sTheta1 = 2.*sqrt(D1)/sHat()/beta/P1;
  double cTheta1 = 0.5*(sHat()-s2+2.*t1_-3.*m2)/beta/rS/P1;
  double theta1 = asin(sTheta1);
  cHalf1_ = cos(0.5*theta1);
  sHalf1_ = sin(0.5*theta1);
  if(cTheta1<0.) {
    swap(cHalf1_,sHalf1_);
    cTheta1=-sqrt(1.-sqr(sTheta1));
  }
  else
    cTheta1=sqrt(1.-sqr(sTheta1));
  // e+
  Energy6 D3 = -0.25*G3;
  double sTheta2 = 2.*sqrt(D3)/sHat()/beta/P2;
  double cTheta2 =-0.5*(sHat()-s1+2.*t2_-3.*m2)/beta/rS/P2;
  double theta2 = asin(sTheta2);
  cHalf2_ = cos(0.5*theta2);
  sHalf2_ = sin(0.5*theta2);
  if(cTheta2<0.) {
    swap(cHalf2_,sHalf2_);
    cTheta2=-sqrt(1.-sqr(sTheta2));
  }
  else
    cTheta2=sqrt(1.-sqr(sTheta2));
  // Energy6 D2 = 0.25*(-sqr(m2)*t2_ + m2*(2*s2*t2_ - sqr(t1_-W2) + t2_*(t1_ + W2)) - t2_*(sqr(s2) + t1_*W2 - s2*(t1_ - t2_ + W2)));
  // hadronic system
  Energy6 D6 = -0.25*m2*(s1-m2)*(s2-m2)
    + 0.125*sHat()*((-m2 + s2 - t1_)*(-m2 + s1 - t2_) + t1_*t2_ - (sHat()-4*m2)*(W2-t1_-t2_));
  Energy6 D5 = D1+D3+2.*D6;
  if(D5<ZERO) return false;
  double sThetaX = 2.*sqrt(max(ZERO,D5))/sHat()/beta/PX;
  double cThetaX = 0.5*(s2-s1+2.*t2_-2.*t1_)/beta/rS/PX;
  // azimuths
  auto Delta4 = -Delta/64./a*sqr(sin(r[3]*Constants::pi));
  if(Delta4>ZERO) return false;
  double sPhi1 = 2.*sqrt(max(ZERO,-Delta4))/sHat()/beta/PX/sThetaX/P1/sTheta1;
  double cPhi1 = (D1+D6)/sqrt(D1*D5);
  double sPhi2 =-2.*sqrt(-Delta4)/sHat()/beta/PX/sThetaX/P2/sTheta2;
  double cPhi2 = (D3+D6)/sqrt(D3*D5);
  // rotation of the system
  double phi = Constants::twopi*UseRandom::rnd();
  double cPhi = cos(phi), sPhi = sin(phi);
  phi1_ = atan2(sPhi1,cPhi1)+phi;
  phi2_ = atan2(sPhi2,cPhi2)+phi;
  // set the momenta
  meMomenta()[2].setE(E1);
  meMomenta()[3].setE(E2);
  meMomenta()[2].setX( P1*sTheta1*(cPhi1*cPhi-sPhi1*sPhi));
  meMomenta()[3].setX( P2*sTheta2*(cPhi2*cPhi-sPhi2*sPhi));
  meMomenta()[2].setY( P1*sTheta1*(cPhi1*sPhi+sPhi1*cPhi));
  meMomenta()[3].setY( P2*sTheta2*(cPhi2*sPhi+sPhi2*cPhi));
  meMomenta()[2].setZ(P1*cTheta1);
  meMomenta()[3].setZ(P2*cTheta2);
  Lorentz5Momentum pX(-PX*sThetaX*cPhi,-PX*sThetaX*sPhi, PX*cThetaX,EX,W);
  // momentum conservation test
  Lorentz5Momentum ptotal;
  for(unsigned int ix=0;ix<4;++ix) {
    if(ix<2) ptotal -=meMomenta()[ix];
    else     ptotal +=meMomenta()[ix];
  }
  ptotal+=pX;
  //double test = (sqr(ptotal.x())+sqr(ptotal.y())+sqr(ptotal.z())+sqr(ptotal.t()))/sHat();
  if(meMomenta().size()==5) {
    meMomenta()[4] = pX;
  }
  else {
    vector<Lorentz5Momentum> pout(meMomenta().size()-4);
    tcPDVector tout(mePartonData().begin()+4,mePartonData().end());
    double jac = amp_->generateKinematics(r+5,W2,pout,tout);
    Boost bv = pX.boostVector();
    for(unsigned int ix=0;ix<pout.size();++ix)
      meMomenta()[ix+4] = pout[ix].boost(bv);
    jacobian(pow(Constants::twopi,3)*jac*jacobian());
  }
  // check the cuts
  vector<LorentzMomentum> out(meMomenta().size()-2);
  tcPDVector tout(meMomenta().size()-2);
  for(unsigned int ix=0;ix<meMomenta().size()-2;++ix) {
    out [ix] = meMomenta()   [2+ix];
    tout[ix] = mePartonData()[2+ix];
  }
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}

double MEff2ffX::me2() const {
  // calculate the leptonic currents
  vector<VectorWaveFunction> current1 =  firstCurrent(mePartonData()[0],meMomenta()[0],meMomenta()[2]);
  vector<VectorWaveFunction> current2 = secondCurrent(mePartonData()[1],meMomenta()[1],meMomenta()[3]);
  DVector save;
  double output = amp_->me2(current1,current2,t1_,t2_,sHat(),
			    vector<Lorentz5Momentum>(meMomenta().begin()+4,meMomenta().end()),
			    cPDVector(mePartonData().begin()+4,mePartonData().end()),save);
  meInfo(save);
  return output;
}

CrossSection MEff2ffX::dSigHatDR() const {
  return 0.5/pow(Constants::twopi,5)*sqr(hbarc)*me2()*jacobian()/sHat()*sqr(sHat()*UnitRemoval::InvE2);
}

Selector<MEBase::DiagramIndex>
MEff2ffX::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    unsigned int id = abs(diags[i]->id())-1;
    if(id<meInfo().size())
      sel.insert(meInfo()[0], i);
    else if ( meInfo().empty() ) sel.insert(1., i);
    else
      assert(false);
  }
  return sel;
}

Selector<const ColourLines *>
MEff2ffX::colourGeometries(tcDiagPtr diag) const {
  return amp_->colourGeometries(1,mePartonData(),diag);
}

IBPtr MEff2ffX::clone() const {
  return new_ptr(*this);
}

IBPtr MEff2ffX::fullclone() const {
  return new_ptr(*this);
}

void MEff2ffX::persistentOutput(PersistentOStream & os) const {
  os << FFPVertex_ << gamma_ << amp_ << currentMode_
     << ounit(Q2_1min_,GeV2) << ounit(Q2_1max_,GeV2)
     << ounit(Q2_2min_,GeV2) << ounit(Q2_2max_,GeV2)
     << ounit(mHatMin_,GeV) << formFactor_;
}

void MEff2ffX::persistentInput(PersistentIStream & is, int) {
  is >> FFPVertex_ >> gamma_ >> amp_ >> currentMode_
     >> iunit(Q2_1min_,GeV2) >> iunit(Q2_1max_,GeV2)
     >> iunit(Q2_2min_,GeV2) >> iunit(Q2_2max_,GeV2) >> iunit(mHatMin_,GeV)
     >> formFactor_;
}

void MEff2ffX::doinit() {
  HwMEBase::doinit();
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  FFPVertex_ = hwsm->vertexFFP();
  gamma_ = getParticleData(ParticleID::gamma);

  vector<Energy2> q2 = {0.*GeV2,1.*GeV2,2.*GeV2,3.*GeV2,4.*GeV2,100.*GeV2};
  vector<double>  ff = {1.,1.,1.,1.,1.,1.};
  formFactor_ = make_InterpolatorPtr(ff,q2,3);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEff2ffX,HwMEBase>
describeHerwigMEff2ffX("Herwig::MEff2ffX", "HwMEGammaGamma.so");

void MEff2ffX::Init() {

  static ClassDocumentation<MEff2ffX> documentation
    ("The MEff2ffX class implements e+e- -> e+e- gamma gamma processes with"
     " gamma gamma-X via the GammaGammaAmplitude");

  static Parameter<MEff2ffX,Energy2> interfaceQ2_1Min
    ("Q2_1Min",
     "The minimum value of Q2 for the off-shell photon from the electron",
     &MEff2ffX::Q2_1min_, GeV2, 0.*GeV2, 0.0*GeV2, Constants::MaxEnergy2,
     false, false, Interface::limited);

  static Parameter<MEff2ffX,Energy2> interfaceQ2_1Max
    ("Q2_1Max",
     "The maximum value of Q2 for the off-shell photon from the electron",
     &MEff2ffX::Q2_1max_, GeV2, 0.*GeV2, 0.0*GeV2, Constants::MaxEnergy2,
     false, false, Interface::limited);

  static Parameter<MEff2ffX,Energy2> interfaceQ2_2Min
    ("Q2_2Min",
     "The minimum value of Q2 for the off-shell photon from the electron",
     &MEff2ffX::Q2_2min_, GeV2, 0.*GeV2, 0.0*GeV2, Constants::MaxEnergy2,
     false, false, Interface::limited);

  static Parameter<MEff2ffX,Energy2> interfaceQ2_2Max
    ("Q2_2Max",
     "The maximum value of Q2 for the off-shell photon from the electron",
     &MEff2ffX::Q2_2max_, GeV2, 0.*GeV2, 0.0*GeV2, Constants::MaxEnergy2,
     false, false, Interface::limited);

  static Reference<MEff2ffX,GammaGammaAmplitude> interfaceAmplitude
    ("Amplitude",
     "The gamma gamma -> X amplitude",
     &MEff2ffX::amp_, false, false, true, false, false);

  static Switch<MEff2ffX,unsigned int> interfaceCurrentMode
    ("CurrentMode",
     "Approximation in which to calculate the electron and position currents",
     &MEff2ffX::currentMode_, 0, false, false);
  static SwitchOption interfaceCurrentModeFull
    (interfaceCurrentMode,
     "Full",
     "Use the full result, calculated to be numerically stable in the collineart limit",
     0);
  static SwitchOption interfaceCurrentModeEquivalentPhoton
    (interfaceCurrentMode,
     "EquivalentPhoton",
     "Use the equivalent photon approximation",
     1);
  static SwitchOption interfaceCurrentModeEquivalentPhotonNoSpin
    (interfaceCurrentMode,
     "EquivalentPhotonNoSpin",
     "Use the equivalent photon approximation, neglecting spin correlations",
     2);

  static Parameter<MEff2ffX,Energy> interfaceMHatMin
    ("MHatMin",
     "The minimum Mhat for the core process",
     &MEff2ffX::mHatMin_, GeV, ZERO, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

}

vector<VectorWaveFunction> MEff2ffX::firstCurrent(tcPDPtr inPart,
						  const Lorentz5Momentum & pin,
						  const Lorentz5Momentum & pout) const {
  double ee = FFPVertex_->electroMagneticCoupling(ZERO);
  Lorentz5Momentum pGamma = pin-pout;
  // form factors
  Energy m = inPart->mass();
  double F1(0.),F2(0.);
  if(abs(inPart->id())==ParticleID::eminus) {
    F1 = 1.;
  }
  else if(abs(inPart->id())==ParticleID::pplus) {
    // defining the form factors F1 and F2 from Ginzburg
    Energy2 q02 = 0.71*GeV2;
    double Ge = 1/sqr(1-t1_/q02);
    double Gm = Ge*sqrt(7.78);
    double tau = -t1_/4./sqr(m);
    F1 = (Ge + tau*Gm)/(1.+ tau);
    F2 = (Gm - Ge)/(1 + tau);
  }
  else {
    LorentzPolarizationVector current = (*formFactor_)(t1_)*UnitRemoval::E/t1_*ee*(inPart->iCharge()/3)*(pin+pout);
    vector<VectorWaveFunction> output;
    output.push_back(VectorWaveFunction(pGamma,gamma_,current));
    return output;
  }
  // calculation of the current
  Complex II(0.,1.);
  Complex phase = exp(II*phi1_);
  Energy Ea = pin.t();
  // full current, safe in small x limit
  vector<LorentzPolarizationVector> current(4);
  if(currentMode_==0) {
    double mr2 = sqr(m/Ea);
    double x = pout.t()/Ea;
    double b = sqrt(1.-mr2), b1 = sqrt(1-mr2/sqr(x)), bb1=(1.+b)*(1.+b1);
    double fact1 = (F1+F2)*2.*Ea*sqrt(x)*sHalf1_/sqrt(bb1)/t1_*UnitRemoval::E*ee;
    double fact2 = (F1+F2)*m/sqrt(x)            /sqrt(bb1)/t1_*UnitRemoval::E*ee;
    // F1 piece
    current[0] = fact1*LorentzPolarizationVector(0.5*phase*(bb1-mr2/x)+b1*sqr(cHalf1_)*(mr2+bb1*x)*cos(phi1_)/(1.-x),
						 -0.5*II*phase*(bb1-mr2/x)+b1*sqr(cHalf1_)*(mr2+bb1*x)*sin(phi1_)/(1.-x),
						 cHalf1_*mr2*(1.+b1-(1.+b)/x)/(sHalf1_*(1.-x))-b1*cHalf1_*sHalf1_*(mr2+bb1*x)/(1.-x),0);
    current[1] = fact2*LorentzPolarizationVector(-cHalf1_*(1 +b-(1 +b1)*x)+2.*b1*cHalf1_*sqr(sHalf1_)*x*(1.+b+(1.+b1)*x)*cos(phi1_)/(phase*(1.-x)),
						 II*cHalf1_*(1.+b-(1.+b1)*x)+2.*b1*cHalf1_*sqr(sHalf1_)*x*(1.+b+(1.+b1)*x)*sin(phi1_)/(phase*(1.-x)),
						 -2.*sHalf1_*x*((1.+b)*(1.+b1*sqr(sHalf1_))-(1.+b1)*(1.-b1*sqr(sHalf1_))*x)/(phase*(1.-x)),0);
    current[2] = fact2*LorentzPolarizationVector(cHalf1_   *(1.+b-(1.+b1)*x)-2.*b1*cHalf1_*phase*sqr(sHalf1_)*x*(1.+b+(1.+b1)*x)*cos(phi1_)/(1.-x),
						 II*cHalf1_*(1.+b-(1.+b1)*x)-2.*b1*cHalf1_*phase*sqr(sHalf1_)*x*(1.+b+(1.+b1)*x)*sin(phi1_)/(1.-x),
						 2.*phase*sHalf1_*x*((1.+b)*(1.+b1*sqr(sHalf1_))-(1.+b1)*(1.-b1*sqr(sHalf1_))*x)/(1.-x),0);
    current[3] = fact1*LorentzPolarizationVector(0.5*(bb1-mr2/x)/phase+b1*sqr(cHalf1_)*(mr2+bb1*x)*cos(phi1_)/(1.-x),
						 0.5*II*(bb1-mr2/x)/phase+b1*sqr(cHalf1_)*(mr2+bb1*x)*sin(phi1_)/(1.-x),
						 cHalf1_*mr2*(1.+b1-(1.+b)/x)/(sHalf1_*(1.-x))-b1*cHalf1_*sHalf1_*(mr2+bb1*x)/(1.-x),0);
    // F2 piece if non-zero
    if(F2!=0) {
      double fact3 = F2*2.*    Ea *sqrt(x)*cHalf1_/sqrt(bb1)*(1.+b+x*(1.+b1))/(1.-x)/t1_*UnitRemoval::E*ee;
      double fact4 = F2*   sqr(Ea)*(mr2+bb1*x)/m/sqrt(x)*sHalf1_/sqrt(bb1)/(1.-x)/t1_*UnitRemoval::E*ee;
      current[0] -= fact3*LorentzPolarizationVector(b1*cHalf1_*sHalf1_*cos(phi1_),b1*cHalf1_*sHalf1_*sin(phi1_),
						    -b1*sqr(sHalf1_)-0.5*mr2*(1.-x)*(1.+x)/((b+b1)*sqr(x)),0.);
      current[1] -= fact4/phase*LorentzPolarizationVector(2.*b1*cHalf1_*sHalf1_*x*cos(phi1_),
							  2.*b1*cHalf1_*sHalf1_*x*sin(phi1_),
							  -x*(2.*b1*sqr(sHalf1_)+mr2*(1.-x)*(1.+x)/((b+b1)*sqr(x))),0);
      current[2] -= fact4*phase*LorentzPolarizationVector(-2.*b1*cHalf1_*sHalf1_*x*cos(phi1_),
							  -2.*b1*cHalf1_*sHalf1_*x*sin(phi1_),
							  x*(2.*b1*sqr(sHalf1_)+mr2*(1.-x)*(1.+x)/((b+b1)*sqr(x))),0.);
      current[3] -= fact3*LorentzPolarizationVector(b1*cHalf1_*sHalf1_*cos(phi1_),b1*cHalf1_*sHalf1_*sin(phi1_),
						    -b1*sqr(sHalf1_)-0.5*mr2*(1.-x)*(1.+x)/((b+b1)*sqr(x)),0.);
    }
  }
  // approximate modes
  else if(currentMode_<=2) {
    Lorentz5Momentum p = pin;
    Lorentz5Momentum n(ZERO,ZERO,-Ea,Ea);
    double z = (n*pout)/(n*p);
    Energy2 pT2 = -z*t1_-sqr(1-z)*sqr(m);
    if(pT2<ZERO) pT2=ZERO;
    Energy pT = sqrt(pT2);
    // Equivalent Photon (with spin correlations)
    if(currentMode_==1) {
      double fact = ee*Ea*UnitRemoval::E/t1_;
      // normal piece
      double fact1 = fact*(F1+F2)*pT/Ea/sqrt(z)/(1.-z);
      double fact2 = fact*(F1+F2)* m/Ea/sqrt(z)*(1.-z);
      current[0] = fact1*LorentzPolarizationVector(phase+z/phase,-II*(phase-z/phase),0., 0.);
      current[1] = fact2*LorentzPolarizationVector(-1.,II,0.,0.);
      current[2] = fact2*LorentzPolarizationVector( 1.,II,0.,0.);
      current[3] = fact1*LorentzPolarizationVector(1./phase+z*phase,II*(1./phase-phase*z),0.,0.);
      // F2 piece if non-zero
      if(F2!=0) {
	double fact3 = fact*F2*    pT /Ea  *(1.+z)/(1.-z)/sqrt(z);
	double fact4 = fact*F2*sqr(pT)/Ea/m       /(1.-z)/sqrt(z);
	current[0] -= fact3      *LorentzPolarizationVector( cos(phi1_), sin(phi1_),0,0);
	current[1] -= fact4/phase*LorentzPolarizationVector( cos(phi1_), sin(phi1_),0,0);
	current[2] -= fact4*phase*LorentzPolarizationVector(-cos(phi1_),-sin(phi1_),0,0);
	current[3] -= fact3      *LorentzPolarizationVector( cos(phi1_), sin(phi1_),0,0);
      }
    }
    // no spin correlations
    else if(currentMode_==2) {
      Energy2 m2=sqr(m);
      double ort=sqrt(0.5);
      double fact = ort*ee*UnitRemoval::E/t1_*
	sqrt(2.*pT2*((4.*sqr(F1))/sqr(1.-z) + (2.*sqr(F1)+4.*F1*F2+3.*sqr(F2))/z)
	     +2.*sqr(F2)*sqr(pT2)/(m2*sqr(1.-z)*z)
	     +4.*sqr(F1+F2)*m2*sqr(1.-z)/z);
      current.resize(2);
      current[0] = fact*ort*LorentzPolarizationVector( 1,-II,0.,0.);
      current[1] = fact*ort*LorentzPolarizationVector(-1,-II,0.,0.);
    }
  }
  // no other option
  else {
    assert(false);
  }
  // test code
  // SpinorWaveFunction    ein (pin,inPart,incoming);
  // SpinorBarWaveFunction eout(pout,inPart,outgoing);
  // vector<SpinorWaveFunction>    fin;
  // vector<SpinorBarWaveFunction> fout;
  // for(unsigned int ix=0;ix<2;++ix) {
  //   ein .reset(ix); fin .push_back(ein );
  //   eout.reset(ix); fout.push_back(eout);
  // }
  // cerr << "testing +z dirn " << inPart->PDGName() << " " << currentMode_ << " " << F1 << " " << F2 << "\n";
  // for(unsigned int ih1=0;ih1<2;++ih1) {
  //   for(unsigned int ih2=0;ih2<2;++ih2) {
  //     LorentzPolarizationVector test  = (F1+F2)*ee/t1_*UnitRemoval::E2*fin[ih1].wave().vectorCurrent(fout[ih2].wave());
  //     test -= F2*(pin+pout)/2./m*ee/t1_*UnitRemoval::E2*fin[ih1].wave().scalar(fout[ih2].wave());
  //     test -= test.t()*pGamma/pGamma.t();
  //     cerr << "testing in current\n"
  //   	   << ih1 << " " << ih2 << " " << test.x() << " " << test.y() << " " << test.z() << " " << test.t() << "\n"
  //   	   << ih1 << " " << ih2 << " " << current[2*ih1+ih2].x() << " " << current[2*ih1+ih2].y() << " "
  //   	   << current[2*ih1+ih2].z() << " " << current[2*ih1+ih2].t() << "\n"
  //    	   << (test.x()-current[2*ih1+ih2].x())/(test.x()+current[2*ih1+ih2].x()) << " "
  //    	   << (test.y()-current[2*ih1+ih2].y())/(test.y()+current[2*ih1+ih2].y()) << " "
  //    	   << (test.z()-current[2*ih1+ih2].z())/(test.z()+current[2*ih1+ih2].z()) << " ";
  //     if(current[2*ih1+ih2].t()!=0.)
  //   	cerr << (test.t()-current[2*ih1+ih2].t())/(test.t()+current[2*ih1+ih2].t());
  //     cerr << "\n";
  //   }
  // }
  // return the currents
  vector<VectorWaveFunction> output; output.reserve(4);
  for(unsigned int ix=0;ix<4;++ix)
    output.push_back(VectorWaveFunction(pGamma,gamma_,current[ix]));
  return output;
}

vector<VectorWaveFunction> MEff2ffX::secondCurrent(tcPDPtr inPart,
						   const Lorentz5Momentum & pin,
						   const Lorentz5Momentum & pout) const {
  Lorentz5Momentum pGamma = pin-pout;
  double ee = FFPVertex_->electroMagneticCoupling(ZERO);
  // form factors
  Energy m = inPart->mass();
  double F1(0.),F2(0.);
  if(abs(inPart->id())==ParticleID::eminus) {
    F1 = 1.;
  }
  else if(abs(inPart->id())==ParticleID::pplus) {
    //defining the form factors F1 and F2 from Ginzburg
    Energy2 q02 = 0.71*GeV2;
    double Ge = 1/sqr(1-t2_/q02);
    double Gm = Ge*sqrt(7.78);
    double tau = -t2_/4./sqr(m);
    F1 = (Ge + tau*Gm)/(1.+ tau);
    F2 = (Gm - Ge)/(1 + tau);
  }
  else {
    LorentzPolarizationVector current = (*formFactor_)(t2_)*UnitRemoval::E/t2_*ee*(inPart->iCharge()/3)*(pin+pout);
    vector<VectorWaveFunction> output;
    output.push_back(VectorWaveFunction(pGamma,gamma_,current));
    return output;
  }
  // calculation of the current
  Complex II(0.,1.);
  Complex phase = exp(II*phi2_);
  Energy Ea = pin.t();
  vector<LorentzPolarizationVector> current(4);
  // full result stable for small pT and x->1
  if(currentMode_==0) {
    double mr2 = sqr(m/Ea);
    double x = pout.t()/Ea;
    double b = sqrt(1.-mr2), b1 = sqrt(1-mr2/sqr(x)), bb1=(1.+b)*(1.+b1);
    double fact1 = (F1+F2)*2.*Ea*sqrt(x)*cHalf2_/sqrt(bb1)/t2_*UnitRemoval::E*ee;
    double fact2 = (F1+F2)*m/sqrt(x)            /sqrt(bb1)/t2_*UnitRemoval::E*ee;
    current[0] = fact1*LorentzPolarizationVector(0.5*(bb1-mr2/x)   +b1*phase*sqr(sHalf2_)*(mr2+bb1*x)*cos(phi2_)/(1.-x),
						 0.5*II*(bb1-mr2/x)+b1*phase*sqr(sHalf2_)*(mr2+bb1*x)*sin(phi2_)/(1.-x),
						 phase*(mr2*sHalf2_*(1.+b-(1.+b1)*x)/cHalf2_/x+b1*cHalf2_*sHalf2_*(mr2+bb1*x))/(1.-x),0);
    current[1] = fact2*LorentzPolarizationVector(sHalf2_*(   (1.+b-(1.+b1)*x)/phase-2.*b1*sqr(cHalf2_)*x*(1.+b+(1.+b1)*x)*cos(phi2_)/(1.-x)),
						 sHalf2_*(II*(1.+b-(1.+b1)*x)/phase-2.*b1*sqr(cHalf2_)*x*(1.+b+(1.+b1)*x)*sin(phi2_)/(1.-x)),
						 -2.*cHalf2_*x*((1.+b)*(1.+b1*sqr(cHalf2_))-(1.+b1)*(1.-b1*sqr(cHalf2_))*x)/(1.-x),0.);
    current[2] = fact2*LorentzPolarizationVector(sHalf2_*(  -phase*(1.+b-(1.+b1)*x)+2.*b1*sqr(cHalf2_)*x*(1.+b+(1.+b1)*x)*cos(phi2_)/(1.-x)),
						 sHalf2_*(II*phase*(1.+b-(1.+b1)*x)+2.*b1*sqr(cHalf2_)*x*(1.+b+(1.+b1)*x)*sin(phi2_)/(1.-x)),
						 2.*cHalf2_*x*((1.+b)*(1.+b1*sqr(cHalf2_))-(1.+b1)*(1.-b1*sqr(cHalf2_))*x)/(1.-x),0);
    current[3] = fact1*LorentzPolarizationVector(0.5*(bb1-mr2/x)+b1*sqr(sHalf2_)*(mr2+bb1*x)*cos(phi2_)/phase/(1.-x),
						 -0.5*II*(bb1-mr2/x)+b1*sqr(sHalf2_)*(mr2+bb1*x)*sin(phi2_)/phase/(1.-x),
						 (mr2*sHalf2_*(1.+b-(1.+b1)*x)/x/cHalf2_+b1*cHalf2_*sHalf2_*(mr2+bb1*x))/phase/(1.-x),0.);
    // F2 piece if non-zero
    if(F2!=0) {
      double fact3 = F2*2.*    Ea   *sqrt(x)*sHalf2_/sqrt(bb1)*(1.+b+x*(1.+b1))/(1.-x)/t2_*UnitRemoval::E*ee;
      double fact4 = F2*   sqr(Ea)/m/sqrt(x)*cHalf2_/sqrt(bb1)*(mr2+bb1*x)     /(1.-x)/t2_*UnitRemoval::E*ee;
      current[0] -= fact3*phase*LorentzPolarizationVector(b1*cHalf2_*sHalf2_*cos(phi2_),b1*cHalf2_*sHalf2_*sin(phi2_),
							  b1*sqr(cHalf2_)+0.5*mr2*(1.-x)*(1.+x)/((b+b1)*sqr(x)),0);
      current[1] -= fact4*      LorentzPolarizationVector(-2.*b1*cHalf2_*sHalf2_*x*cos(phi2_),
							  -2.*b1*cHalf2_*sHalf2_*x*sin(phi2_),
							  -x*(2.*b1*sqr(cHalf2_)+mr2*(1.-x)*(1.+x)/((b+b1)*sqr(x))),0.);
      current[2] -= fact4*      LorentzPolarizationVector(2.*b1*cHalf2_*sHalf2_*x*cos(phi2_),
							  2.*b1*cHalf2_*sHalf2_*x*sin(phi2_),
							  x*(2.*b1*sqr(cHalf2_)+mr2*(1.-x)*(1.+x)/((b+b1)*sqr(x))),0.);
      current[3] -= fact3/phase*LorentzPolarizationVector(b1*cHalf2_*sHalf2_*cos(phi2_),b1*cHalf2_*sHalf2_*sin(phi2_),
							  b1*sqr(cHalf2_)+0.5*mr2*(1.-x)*(1.+x)/((b+b1)*sqr(x)),0);
    }
  }
  // approximate modes
  else if(currentMode_<=2) {
    Lorentz5Momentum p = pin;
    Lorentz5Momentum n(ZERO,ZERO,Ea,Ea);
    double z = (n*pout)/(n*p);
    Energy2 pT2 = -z*t2_-sqr(1-z)*sqr(m);
    if(pT2<ZERO) pT2=ZERO;
    Energy pT = sqrt(pT2);
    // Equivalent Photon with spin correlations
    if(currentMode_<=1) {
      double fact = ee*Ea*UnitRemoval::E/t2_;
      double fact1 = fact*(F1+F2)*pT/Ea/sqrt(z)/(1.-z);
      double fact2 = fact*(F1+F2)*m /Ea/sqrt(z)*(1.-z);
      current[0] = fact1      *LorentzPolarizationVector(1.+z*sqr(phase), II*(1.-z*sqr(phase)),0., 0.);
      current[1] = fact2/phase*LorentzPolarizationVector( 1.,II,0.,0.);
      current[2] = fact2*phase*LorentzPolarizationVector(-1.,II,0.,0.);
      current[3] = fact1      *LorentzPolarizationVector(1.+z/sqr(phase),-II*(1.-z/sqr(phase)),0.,0.);
      // F2 piece if non-zero
      if(F2!=0) {
	double fact1 = fact*F2*    pT *(1.+z)/(Ea  *(1.-z)*sqrt(z));
	double fact2 = fact*F2*sqr(pT)       /(Ea*m*(1.-z)*sqrt(z));
	current[0] -= fact1*phase*LorentzPolarizationVector( cos(phi2_), sin(phi2_),0,0);
	current[1] -= fact2      *LorentzPolarizationVector(-cos(phi2_),-sin(phi2_),0,0);
	current[2] -= fact2      *LorentzPolarizationVector( cos(phi2_), sin(phi2_),0,0);
	current[3] -= fact1/phase*LorentzPolarizationVector( cos(phi2_), sin(phi2_),0,0);
      }
    }
    // no spin correlations
    else if(currentMode_==2) {
      double ort=sqrt(0.5);
      Energy2 m2=sqr(m);
      double fact = ort*ee*UnitRemoval::E/t2_*
	sqrt(2.*pT2*((4.*sqr(F1))/sqr(1.-z) + (2.*sqr(F1)+4.*F1*F2+3.*sqr(F2))/z)
	     +2.*sqr(F2)*sqr(pT2)/(m2*sqr(1.-z)*z)
	     +4.*sqr(F1+F2)*m2*sqr(1.-z)/z);
      current.resize(2);
      current[0] = fact*ort*LorentzPolarizationVector(-1,-II,0.,0.);
      current[1] = fact*ort*LorentzPolarizationVector( 1,-II,0.,0.);
    }
  }
  // no other option
  else {
    assert(false);
  }
  // test code
  // SpinorWaveFunction    ein (pin,inPart,incoming);
  // SpinorBarWaveFunction eout(pout,inPart,outgoing);
  // vector<SpinorWaveFunction>    fin;
  // vector<SpinorBarWaveFunction> fout;
  // for(unsigned int ix=0;ix<2;++ix) {
  //   ein .reset(ix); fin .push_back(ein );
  //   eout.reset(ix); fout.push_back(eout);
  // }
  // cerr << "testing -z dirn " << inPart->PDGName() << " " << currentMode_ << " " << F1 << " " << F2 << "\n";
  // for(unsigned int ih1=0;ih1<2;++ih1) {
  //   for(unsigned int ih2=0;ih2<2;++ih2) {
  //     LorentzPolarizationVector test  = (F1+F2)*ee/t2_*UnitRemoval::E2*fin[ih1].wave().vectorCurrent(fout[ih2].wave());
  //     test -= F2*(pin+pout)/2./m*ee/t2_*UnitRemoval::E2*fin[ih1].wave().scalar(fout[ih2].wave());
  //     test -= test.t()*pGamma/pGamma.t();
  //     cerr << "testing in current\n"
  //   	   << ih1 << " " << ih2 << " " << test.x() << " " << test.y() << " " << test.z() << " " << test.t() << "\n"
  //   	   << ih1 << " " << ih2 << " " << current[2*ih1+ih2].x() << " " << current[2*ih1+ih2].y() << " "
  //   	   << current[2*ih1+ih2].z() << " " << current[2*ih1+ih2].t() << "\n"
  //    	   << (test.x()-current[2*ih1+ih2].x())/(test.x()+current[2*ih1+ih2].x()) << " "
  //    	   << (test.y()-current[2*ih1+ih2].y())/(test.y()+current[2*ih1+ih2].y()) << " "
  //    	   << (test.z()-current[2*ih1+ih2].z())/(test.z()+current[2*ih1+ih2].z()) << " ";
  //     if(current[2*ih1+ih2].t()!=0.)
  //   	cerr << (test.t()-current[2*ih1+ih2].t())/(test.t()+current[2*ih1+ih2].t());
  //     cerr << "\n";
  //   }
  // }
  // return the answer
  vector<VectorWaveFunction> output; output.reserve(4);
  for(unsigned int ix=0;ix<4;++ix)
    output.push_back(VectorWaveFunction(pGamma,gamma_,current[ix]));
  return output;
}

void MEff2ffX::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  tParticleVector hard;
  hard.reserve(sub->outgoing().size()+2);
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  for(unsigned int ix=0;ix<sub->outgoing().size();++ix)
    hard.push_back(sub->outgoing()[ix]);
  // calculate the fermionic currents
  vector<VectorWaveFunction> current1 =  firstCurrent(hard[0]->dataPtr(),hard[0]->momentum(),hard[2]->momentum());
  vector<VectorWaveFunction> current2 = secondCurrent(hard[1]->dataPtr(),hard[1]->momentum(),hard[3]->momentum());
  // wavefunctions for the fermions
  vector<SpinorWaveFunction>    f1,f2;
  vector<SpinorBarWaveFunction> a1,a2;
  if(hard[0]->dataPtr()->iSpin()==PDT::Spin1Half) {
    if(hard[0]->id()>0) {
      SpinorWaveFunction   (f1,hard[0],incoming,false,true);
      SpinorBarWaveFunction(a1,hard[2],outgoing,true,true);
    }
    else {
      SpinorWaveFunction   (f1,hard[2],outgoing,true,true);
      SpinorBarWaveFunction(a1,hard[0],incoming,false,true);
    }
  }
  else {
    ScalarWaveFunction(hard[0],incoming,false);
    ScalarWaveFunction(hard[2],outgoing,true );
  }
  if(hard[1]->dataPtr()->iSpin()==PDT::Spin1Half) {
    if(hard[1]->id()>0) {
      SpinorWaveFunction   (f2,hard[1],incoming,false,true);
      SpinorBarWaveFunction(a2,hard[3],outgoing,true,true);
    }
    else {
      SpinorWaveFunction   (f2,hard[3],outgoing,true,true);
      SpinorBarWaveFunction(a2,hard[1],incoming,false,true);
    }
  }
  else {
    ScalarWaveFunction(hard[1],incoming,false);
    ScalarWaveFunction(hard[3],outgoing,true );
  }
  tParticleVector pTemp(hard.begin()+4,hard.end());
  ProductionMatrixElement me = amp_->me(current1,current2,pTemp);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<hard.size();++ix) {
    hard[ix]->spinInfo()->productionVertex(hardvertex);
  }
}
