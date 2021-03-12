// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2eeX class.
//

#include "MEee2eeX.h"
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

using namespace Herwig;

void MEee2eeX::getDiagrams() const {
  vector<DiagPtr> diags = amp_->getDiagrams(1);
  for(DiagPtr diag : diags) add(diag);
}

Energy2 MEee2eeX::scale() const {
  return sHat();
}

int MEee2eeX::nDim() const {
  return 4+amp_->nDim(1);
}

void MEee2eeX::setKinematics() {
  HwMEBase::setKinematics(); // Always call the base class method first.
}

bool MEee2eeX::generateKinematics(const double * r) {
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
    W = amp_->generateW(r[4],tcPDVector(mePartonData().begin()+4,mePartonData().end()),rS-2.*m,jacW);
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
  double test = (sqr(ptotal.x())+sqr(ptotal.y())+sqr(ptotal.z())+sqr(ptotal.t()))/sHat();
  if(test>1e-20) cerr << "testing the sum " << ptotal/GeV << " " << test << "\n";
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

double MEee2eeX::me2() const {
  // calculate the leptonic currents
  vector<VectorWaveFunction> eCurrent = electronCurrent();
  vector<VectorWaveFunction> pCurrent = positronCurrent();
  DVector save;
  double output = amp_->me2(eCurrent,pCurrent,t1_,t2_,sHat(),
			    vector<Lorentz5Momentum>(meMomenta().begin()+4,meMomenta().end()),
			    cPDVector(mePartonData().begin()+4,mePartonData().end()),save);
  meInfo(save);
  return output;
}

CrossSection MEee2eeX::dSigHatDR() const {
  return 0.5/pow(Constants::twopi,5)*sqr(hbarc)*me2()*jacobian()/sHat()*sqr(sHat()*UnitRemoval::InvE2);
}

Selector<MEBase::DiagramIndex>
MEee2eeX::diagrams(const DiagramVector & diags) const {
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
MEee2eeX::colourGeometries(tcDiagPtr diag) const {
  return amp_->colourGeometries(1,mePartonData(),diag);
}

IBPtr MEee2eeX::clone() const {
  return new_ptr(*this);
}

IBPtr MEee2eeX::fullclone() const {
  return new_ptr(*this);
}

void MEee2eeX::persistentOutput(PersistentOStream & os) const {
  os << FFPVertex_ << gamma_ << amp_ << currentMode_
     << ounit(Q2_1min_,GeV2) << ounit(Q2_1max_,GeV2)
     << ounit(Q2_2min_,GeV2) << ounit(Q2_2max_,GeV2);
}

void MEee2eeX::persistentInput(PersistentIStream & is, int) {
  is >> FFPVertex_ >> gamma_ >> amp_ >> currentMode_
     >> iunit(Q2_1min_,GeV2) >> iunit(Q2_1max_,GeV2)
     >> iunit(Q2_2min_,GeV2) >> iunit(Q2_2max_,GeV2);
}

void MEee2eeX::doinit() {
  HwMEBase::doinit();
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  FFPVertex_ = hwsm->vertexFFP();
  gamma_ = getParticleData(ParticleID::gamma);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEee2eeX,HwMEBase>
describeHerwigMEee2eeX("Herwig::MEee2eeX", "HwMEGammaGamma.so");

void MEee2eeX::Init() {

  static ClassDocumentation<MEee2eeX> documentation
    ("The MEee2eeX class implements e+e- -> e+e- gamma gamma processes with"
     " gamma gamma-X via the GammaGammaAmplitude");

  static Parameter<MEee2eeX,Energy2> interfaceQ2_1Min
    ("Q2_1Min",
     "The minimum value of Q2 for the off-shell photon from the electron",
     &MEee2eeX::Q2_1min_, GeV2, 0.*GeV2, 0.0*GeV2, Constants::MaxEnergy2,
     false, false, Interface::limited);

  static Parameter<MEee2eeX,Energy2> interfaceQ2_1Max
    ("Q2_1Max",
     "The maximum value of Q2 for the off-shell photon from the electron",
     &MEee2eeX::Q2_1max_, GeV2, 0.*GeV2, 0.0*GeV2, Constants::MaxEnergy2,
     false, false, Interface::limited);

  static Parameter<MEee2eeX,Energy2> interfaceQ2_2Min
    ("Q2_2Min",
     "The minimum value of Q2 for the off-shell photon from the electron",
     &MEee2eeX::Q2_2min_, GeV2, 0.*GeV2, 0.0*GeV2, Constants::MaxEnergy2,
     false, false, Interface::limited);

  static Parameter<MEee2eeX,Energy2> interfaceQ2_2Max
    ("Q2_2Max",
     "The maximum value of Q2 for the off-shell photon from the electron",
     &MEee2eeX::Q2_2max_, GeV2, 0.*GeV2, 0.0*GeV2, Constants::MaxEnergy2,
     false, false, Interface::limited);
  
  static Reference<MEee2eeX,GammaGammaAmplitude> interfaceAmplitude
    ("Amplitude",
     "The gamma gamma -> X amplitude",
     &MEee2eeX::amp_, false, false, true, false, false);

  static Switch<MEee2eeX,unsigned int> interfaceCurrentMode
    ("CurrentMode",
     "Approximation in which to calculate the electron and position currents",
     &MEee2eeX::currentMode_, 0, false, false);
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

}

vector<VectorWaveFunction> MEee2eeX::electronCurrent() const {
  Lorentz5Momentum pGamma = meMomenta()[0]-meMomenta()[2];
  double ee = FFPVertex_->electroMagneticCoupling(ZERO);
  Complex II(0.,1.);
  Complex phase = exp(II*phi1_);
  Energy m = mePartonData()[0]->mass();
  vector<VectorWaveFunction> output; output.reserve(4);
  // fuill current, safe in small x limit
  if(currentMode_==0) {
    double mr2 = sqr(m/meMomenta()[0].t());
    double x = meMomenta()[2].t()/meMomenta()[0].t();
    double b = sqrt(1.-mr2), b1 = sqrt(1-mr2/sqr(x)), bb1=(1.+b)*(1.+b1);
    double fact1 = 2.*meMomenta()[0].t()*sqrt(x)/t1_*UnitRemoval::E*ee*sHalf1_/sqrt(bb1);
    double fact2 = m/sqrt(x)/t1_*UnitRemoval::E*ee/sqrt(bb1);
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact1*LorentzPolarizationVector(0.5*phase*(bb1-mr2/x)+
										       b1*(mr2+bb1*x)*cos(phi1_)*sqr(cHalf1_)/(1.-x),
										       -0.5*II*phase*(bb1-mr2/x)
										       +b1*(mr2+bb1*x)*sqr(cHalf1_)*sin(phi1_)/(1.-x),
										       (1.+b1-(1.+b)/x)*mr2*cHalf1_/sHalf1_/(1.-x)
										       -b1*(mr2+bb1*x)*cHalf1_*sHalf1_/(1.-x),0)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact2*LorentzPolarizationVector(-(1.+b-(1.+b1)*x)*cHalf1_,
										       II*(1.+b-(1.+b1)*x)*cHalf1_,
										       (1.+b-(1.+b1)*x)*sHalf1_/phase,
										       (1.+b+(1.+b1)*x)*sHalf1_/phase)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact2*LorentzPolarizationVector((1.+b-(1.+b1)*x)*cHalf1_,
										       II*(1.+b-(1.+b1)*x)*cHalf1_,
										       -phase*(1.+b-(1.+b1)*x)*sHalf1_,
										       -phase*(1.+b+(1.+b1)*x)*sHalf1_)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact1*LorentzPolarizationVector(0.5/phase*(bb1-mr2/x)+ 
										       b1*(mr2+bb1*x)*cos(phi1_)*sqr(cHalf1_)/(1.-x),
										       0.5/phase*II*(bb1-mr2/x)+ 
										       b1*(mr2+bb1*x)*sqr(cHalf1_)*sin(phi1_)/(1.-x),
										       mr2*(1.+b1-(1.+b)/x)*cHalf1_/sHalf1_/(1.-x) + 
										       -b1*(mr2+bb1*x)*cHalf1_*sHalf1_/(1.-x),0)));
    // test code
    // SpinorWaveFunction    ein (meMomenta()[0],mePartonData()[0],incoming);
    // SpinorBarWaveFunction eout(meMomenta()[2],mePartonData()[2],outgoing);
    // vector<SpinorWaveFunction>    fin;
    // vector<SpinorBarWaveFunction> fout;
    // for(unsigned int ix=0;ix<2;++ix) {
    //   ein .reset(ix); fin .push_back(ein );
    //   eout.reset(ix); fout.push_back(eout);
    // }
    // for(unsigned int ih1=0;ih1<2;++ih1) {
    //   for(unsigned int ih2=0;ih2<2;++ih2) {
    //     LorentzPolarizationVector test = FFPVertex_->evaluate(ZERO,1,gamma_,fin[ih1],fout[ih2]).wave();
    //     if(ih1==ih2) test -= test.t()*pGamma/pGamma.t();
    //     cerr << "testing in current\n"
    // 	   << ih1 << " " << ih2 << " " << test.x() << " " << test.y() << " " << test.z() << " " << test.t() << "\n"
    // 	   << ih1 << " " << ih2 << " " << output[2*ih1+ih2].x() << " " << output[2*ih1+ih2].y() << " "
    // 	   << output[2*ih1+ih2].z() << " " << output[2*ih1+ih2].t() << "\n"
    // 	   << (test.x()-output[2*ih1+ih2].x())/(test.x()+output[2*ih1+ih2].x()) << " "
    // 	   << (test.y()-output[2*ih1+ih2].y())/(test.y()+output[2*ih1+ih2].y()) << " "
    // 	   << (test.z()-output[2*ih1+ih2].z())/(test.z()+output[2*ih1+ih2].z()) << " ";
    //     if(output[2*ih1+ih2].t()!=0.)
    // 	cerr << (test.t()-output[2*ih1+ih2].t())/(test.t()+output[2*ih1+ih2].t());
    //     cerr << "\n";
    //   }
    // }
  }
  // Equivalent Photon
  else if(currentMode_==1) {
    Lorentz5Momentum p = meMomenta()[0];
    Lorentz5Momentum n(ZERO,ZERO,-meMomenta()[0].t(),meMomenta()[0].t());
    double z = (n*meMomenta()[2])/(n*p);
    Energy2 pT2 = -z*t1_-sqr(1-z)*sqr(m);
    if(pT2<ZERO) pT2=ZERO;
    Energy pT = sqrt(pT2);
    Energy Ea = meMomenta()[0].t();
    double fact = ee*Ea*UnitRemoval::E/t1_;
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact*LorentzPolarizationVector(  pT/Ea/phase*(sqr(phase)+z)/sqrt(z)/(1.-z),
											-II*pT/Ea/phase*(sqr(phase)-z)/sqrt(z)/(1.-z),0., 0.)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact*LorentzPolarizationVector(-Complex(m/Ea/sqrt(z)*(1.-z)),II*m/Ea*(1.-z)/sqrt(z),0.,0.)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact*LorentzPolarizationVector( Complex(m/Ea/sqrt(z)*(1.-z)),II*m/Ea*(1.-z)/sqrt(z),0.,0.)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact*LorentzPolarizationVector( pT/Ea/sqrt(z)/(1.-z)/phase*(1.+z*sqr(phase)),
										       -II*pT/Ea/sqrt(z)/(1.-z)*(1.-sqr(phase)*z)/phase,0.,0.)));
  }
  // no other option
  else {
    assert(false);
  }
  return output;
}

vector<VectorWaveFunction> MEee2eeX::positronCurrent() const {
  Lorentz5Momentum pGamma = meMomenta()[1]-meMomenta()[3];
  double ee = FFPVertex_->electroMagneticCoupling(ZERO);
  Complex II(0.,1.);
  Complex phase = exp(II*phi2_);
  Energy m = mePartonData()[1]->mass();
  vector<VectorWaveFunction> output; output.reserve(4);
  if(currentMode_==0) {
    double mr2 = sqr(m/meMomenta()[1].t());
    double x = meMomenta()[3].t()/meMomenta()[1].t();
    double b = sqrt(1.-mr2), b2 = sqrt(1-mr2/sqr(x)), bb2=(1.+b)*(1.+b2);
    double fact1 = 2.*meMomenta()[1].t()*sqrt(x)/t2_*UnitRemoval::E*ee*cHalf2_/sqrt(bb2);
    double fact2 = m/sqrt(x)/t2_*UnitRemoval::E*ee/sqrt(bb2);
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact1*LorentzPolarizationVector(0.5*(bb2-mr2/x) +
										       b2*phase*(mr2+bb2*x)*cos(phi2_)*
										       sqr(sHalf2_)/(1.-x),
										       0.5*II*(bb2-mr2/x) + 
										       b2*phase*(mr2+bb2*x)*sin(phi2_)*
										       sqr(sHalf2_)/(1.-x),
										       phase*(b2*(mr2+bb2*x)*cHalf2_*sHalf2_ + 
											      mr2*(1.+b-(1.+b2)*x)*sHalf2_/cHalf2_/x)/(1.-x),0)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact2*LorentzPolarizationVector((1.+b-(1.+b2)*x)*sHalf2_/phase,
										       II*(1.+b-(1.+b2)*x)*sHalf2_/phase,
										       (1.+b-(1.+b2)*x)*cHalf2_,
										       -(1.+b+(1.+b2)*x)*cHalf2_)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact2*LorentzPolarizationVector(-phase*(1.+b-(1.+b2)*x)*sHalf2_,
										       II*phase*(1.+b-(1.+b2)*x)*sHalf2_,
										       -(1.+b-(1.+b2)*x)*cHalf2_,
										       (1.+b+(1.+b2)*x)*cHalf2_)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact1*LorentzPolarizationVector(0.5*(bb2-mr2/x) + 
										       b2/phase*(mr2+bb2*x)*cos(phi2_)*
										       sqr(sHalf2_)/(1.-x),
										       -0.5*II*(bb2-mr2/x) + 
										       b2/phase*(mr2+bb2*x)*sin(phi2_)*
										       sqr(sHalf2_)/(1.-x),
										       (b2*(mr2+bb2*x)*cHalf2_*sHalf2_ + 
											(mr2*(1.+b-(1.+b2)*x)*sHalf2_/cHalf2_)/x)/phase/(1.-x),0.)));
    // test code
    // SpinorBarWaveFunction pin (meMomenta()[1],mePartonData()[1],incoming);
    // SpinorWaveFunction    pout(meMomenta()[3],mePartonData()[3],outgoing);
    // vector<SpinorBarWaveFunction> fin;
    // vector<SpinorWaveFunction>    fout;
    // for(unsigned int ix=0;ix<2;++ix) {
    //   pin .reset(ix); fin .push_back(pin );
    //   pout.reset(ix); fout.push_back(pout);
    // }
    // for(unsigned int ih1=0;ih1<2;++ih1) {
    //   for(unsigned int ih2=0;ih2<2;++ih2) {
    //     LorentzPolarizationVector test = FFPVertex_->evaluate(ZERO,1,gamma_,fout[ih2],fin[ih1]).wave();
    //     if(ih1==ih2) test -= test.t()*pGamma/pGamma.t();
    //     cerr << "testing in current\n"
    // 	   << ih1 << " " << ih2 << " " << test.x() << " " << test.y() << " " << test.z() << " " << test.t() << "\n"
    // 	   << ih1 << " " << ih2 << " " << output[2*ih1+ih2].x() << " " << output[2*ih1+ih2].y() << " "
    // 	   << output[2*ih1+ih2].z() << " " << output[2*ih1+ih2].t() << "\n"
    // 	   << (test.x()-output[2*ih1+ih2].x())/(test.x()+output[2*ih1+ih2].x()) << " "
    // 	   << (test.y()-output[2*ih1+ih2].y())/(test.y()+output[2*ih1+ih2].y()) << " "
    // 	   << (test.z()-output[2*ih1+ih2].z())/(test.z()+output[2*ih1+ih2].z()) << " ";
    //     if(output[2*ih1+ih2].t()!=0.) cerr << (test.t()-output[2*ih1+ih2].t())/(test.t()+output[2*ih1+ih2].t());
    //     cerr << "\n";
    //   }
    // }
  }
  // Equivalent Photon
  else if(currentMode_==1) {
    Lorentz5Momentum p = meMomenta()[1];
    Lorentz5Momentum n(ZERO,ZERO,meMomenta()[1].t(),meMomenta()[1].t());
    double z = (n*meMomenta()[3])/(n*p);
    Energy2 pT2 = -z*t2_-sqr(1-z)*sqr(m);
    if(pT2<ZERO) pT2=ZERO;
    Energy pT = sqrt(pT2);
    Energy Eb = meMomenta()[1].t();
    double fact = ee*Eb*UnitRemoval::E/t2_;
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact*LorentzPolarizationVector(   pT/Eb/(1.-z)/sqrt(z)*(1.+z*sqr(phase)),
										      II*pT/Eb/(1.-z)/sqrt(z)*(1.-z*sqr(phase)),0., 0.)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact*LorentzPolarizationVector(   m/Eb*(1.-z)/sqrt(z)/phase,
										      II*m/Eb*(1.-z)/sqrt(z)/phase,0.,0.)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact*LorentzPolarizationVector(-  m/Eb*phase*(1.-z)/sqrt(z),
										      II*m/Eb*phase*(1.-z)/sqrt(z),0.,0.)));
    output.push_back(VectorWaveFunction(pGamma,gamma_, fact*LorentzPolarizationVector(    pT/Eb/(1.-z)/sqrt(z)*(1.+z/sqr(phase)),
										      -II*pT/Eb/(1.-z)/sqrt(z)*(1.-z/sqr(phase)),0.,0.)));
  }
  // no other option
  else {
    assert(false);
  }
  return output;
}





// Energy2 MEee2eePseudoScalar::scale() const {
//   return sHat();
// }



// Selector<MEBase::DiagramIndex>
// MEee2eePseudoScalar::diagrams(const DiagramVector & diags) const {
//   Selector<DiagramIndex> sel;
//   for ( DiagramIndex i = 0; i < diags.size(); ++i )
//     sel.insert(1.0, i);
//   return sel;
// }




// void MEee2eePseudoScalar::constructVertex(tSubProPtr ) {
// //   // extract the particles in the hard process
// //   ParticleVector hard;
// //   hard.push_back(sub->incoming().first);
// //   hard.push_back(sub->incoming().second);
// //   hard.push_back(sub->outgoing()[0]);
// //   hard.push_back(sub->outgoing()[1]);
// //   hard.push_back(sub->outgoing()[2]);
// //   // ensure right order
// //   if(hard[0]->id()<0) swap(hard[0],hard[1]);
// //   if(hard[3]->id()==ParticleID::h0) swap(hard[2],hard[3]);
// //   if(hard[4]->id()==ParticleID::h0) swap(hard[2],hard[4]);
// //   if(hard[3]->id()<0) swap(hard[3],hard[4]);
// //   vector<SpinorWaveFunction>    fin,aout;
// //   vector<SpinorBarWaveFunction> ain,fout;
// //   SpinorWaveFunction(   fin ,hard[0],incoming,false,true);
// //   SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
// //   ScalarWaveFunction(        hard[2],outgoing,true,true);
// //   SpinorBarWaveFunction(fout,hard[3],outgoing,true ,true);
// //   SpinorWaveFunction(   aout,hard[4],outgoing,true ,true);
// //   helicityME(fin,ain,fout,aout,true);
// //   // construct the vertex
// //   HardVertexPtr hardvertex=new_ptr(HardVertex());
// //   // set the matrix element for the vertex
// //   hardvertex->ME(_me);
// //   // set the pointers and to and from the vertex
// //   for(unsigned int ix=0;ix<5;++ix) {
// //     hard[ix]->spinInfo()->productionVertex(hardvertex);
// //   }
// }
