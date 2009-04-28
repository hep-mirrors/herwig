// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2HiggsVBFPowheg class.
//

#include "MEPP2HiggsVBFPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/PDT/StandardMatchers.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"

using namespace Herwig;

MEPP2HiggsVBFPowheg::MEPP2HiggsVBFPowheg() 
  : scaleOpt_(1),  muF_(100.*GeV), scaleFact_(1.), contrib_(1), power_(0.1)
{}

int MEPP2HiggsVBFPowheg::nDim() const {
  return MEPP2HiggsVBF::nDim()+3;
}

bool MEPP2HiggsVBFPowheg::generateKinematics(const double * r) {
  // Born kinematics
  if(!MEPP2HiggsVBFPowheg::generateKinematics(r)) return false;

  int in;
  // hadron and momentum fraction
  if(UseRandom::rnd()< 0.5) {
    in = 0;
    _hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr());
    _xB = lastX1();
  }
  else {
    in = 1;
    _hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr());
    _xB = lastX2();
  }
  // set Q2 process momenta
  if(in == 0){
    _q2 = -(meMomenta()[0]-meMomenta()[2]).m2();
    _pa = meMomenta()[5]+meMomenta()[2]-meMomenta()[0];
    _pb = meMomenta()[0];
    _pc = meMomenta()[2];
  }
  else {
    _q2 = -(meMomenta()[1]-meMomenta()[3]).m2();
    _pa = meMomenta()[5]+meMomenta()[3]-meMomenta()[1];
    _pb = meMomenta()[1];
    _pc = meMomenta()[3];
  }

  // xp
  int ndim=nDim();
  double rhomin = pow(1.-_xB,1.-power_); 
  double rho = r[ndim-1]*rhomin;
  _xp = 1.-pow(rho,1./(1.-power_));
  // zp 
  //double lambda = r[ndim-2]*rhomin;
  //_zp = 1.-pow(lambda,1./(1.-power_));
  _zp = r[ndim-2];
  // phi
  _phi = r[ndim-3]*Constants::twopi;

  jac_  = rhomin/(1.-power_)*pow(1.-_xp,power_);
  //jac_ *= pow(1.-_zp,power_)/(1.-power_);
  jac_ *= 1./Constants::twopi;
  jacobian(jacobian()*jac_);
  return true;
}

IBPtr MEPP2HiggsVBFPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2HiggsVBFPowheg::fullclone() const {
  return new_ptr(*this);
}


void MEPP2HiggsVBFPowheg::persistentOutput(PersistentOStream & os) const {
  os << ounit(muF_,GeV) << scaleFact_ << scaleOpt_ << contrib_
     << _sinW << _cosW << ounit(_mz2,GeV2) << power_;
}

void MEPP2HiggsVBFPowheg::persistentInput(PersistentIStream & is, int) {
  is >> iunit(muF_,GeV) >> scaleFact_ >> scaleOpt_ >> contrib_
     >> _sinW >> _cosW >> iunit(_mz2,GeV2) >> power_;
}

ClassDescription<MEPP2HiggsVBFPowheg>
 MEPP2HiggsVBFPowheg::initMEPP2HiggsVBFPowheg;
// Definition of the static class description member.

void MEPP2HiggsVBFPowheg::Init() {

  static ClassDocumentation<MEPP2HiggsVBFPowheg> documentation
    ("The MENeutralCurrentDISPowheg class implements the NLO matrix element"
     " for neutral current DIS in the Powheg scheme.");

  static Switch<MEPP2HiggsVBFPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEPP2HiggsVBFPowheg::contrib_, 1, false, false);
  static SwitchOption interfaceContributionLeadingOrder
    (interfaceContribution,
     "LeadingOrder",
     "Just generate the leading order cross section",
     0);
  static SwitchOption interfaceContributionPositiveNLO
    (interfaceContribution,
     "PositiveNLO",
     "Generate the positive contribution to the full NLO cross section",
     1);
  static SwitchOption interfaceContributionNegativeNLO
    (interfaceContribution,
     "NegativeNLO",
     "Generate the negative contribution to the full NLO cross section",
     2);

  static Switch<MEPP2HiggsVBFPowheg,unsigned int> interfaceScaleOption
    ("ScaleOption",
     "Option for the choice of factorization (and renormalization) scale",
     &MEPP2HiggsVBFPowheg::scaleOpt_, 1, false, false);
  static SwitchOption interfaceDynamic
    (interfaceScaleOption,
     "Dynamic",
     "Dynamic factorization scale equal to the current sqrt(sHat())",
     1);
  static SwitchOption interfaceFixed
    (interfaceScaleOption,
     "Fixed",
     "Use a fixed factorization scale set with FactorizationScaleValue",
     2);

  static Parameter<MEPP2HiggsVBFPowheg,Energy> interfaceFactorizationScale
    ("FactorizationScale",
     "Value to use in the event of a fixed factorization scale",
     &MEPP2HiggsVBFPowheg::muF_, GeV, 100.0*GeV, 1.0*GeV, 500.0*GeV,
     true, false, Interface::limited);

  static Parameter<MEPP2HiggsVBFPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before Q2 if using a running scale",
     &MEPP2HiggsVBFPowheg::scaleFact_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<MEPP2HiggsVBFPowheg,double> interfaceSamplingPower
    ("SamplingPower",
     "Power for the sampling of xp",
     &MEPP2HiggsVBFPowheg::power_, 0.6, 0.0, 1.,
     false, false, Interface::limited);

}

Energy2 MEPP2HiggsVBFPowheg::scale() const {
  return scaleOpt_ == 1 ? 
    sqr(scaleFact_)*MEPP2HiggsVBF::scale() : sqr(scaleFact_*muF_);
}

CrossSection MEPP2HiggsVBFPowheg::dSigHatDR() const {
  return NLOWeight()*MEPP2HiggsVBF::dSigHatDR();
}

double MEPP2HiggsVBFPowheg::NLOWeight() const {
  // If only leading order is required return 1:
  if(contrib_==0) return 1.;
  // scale and prefactors
  Energy2 mu2(scale());
  double aS = 2.*SM().alphaS(mu2);
  double CFfact = 4./3.*aS/Constants::twopi;
  double TRfact = 1./2.*aS/Constants::twopi;

  // Boost
  Axis axis(_pa.vect().unit());
  LorentzRotation rot;
  double sinth(sqr(axis.x())+sqr(axis.y()));
  if(axis.perp2()>1e-20) {
    rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    rot.rotateX(Constants::pi);
  }
  if(abs(1.-_pa.e()/_pa.vect().mag())>1e-6) rot.boostZ(_pa.e()/_pa.vect().mag());
  _pb *= rot;
  if(_pb.perp2()/GeV2>1e-20) {
    Boost trans = -1./_pb.e()*_pb.vect();
    trans.setZ(0.);
    rot.boost(trans);
  }
  
  // LO + dipole subtracted virtual + collinear quark bit with LO pdf
  double virt = 1.+CFfact*(-4.5-1./3.*sqr(Constants::pi)+1.5*log(_q2/mu2/(1.-_xB))
			   +2.*log(1.-_xB)*log(_q2/mu2)+sqr(log(1.-_xB)));
  virt /= jac_;
  // PDF from leading-order
  double loPDF = _hadron->pdf()->xfx(_hadron,mePartonData()[1],mu2,_xB)/_xB;
  // NLO gluon PDF
  tcPDPtr gluon = getParticleData(ParticleID::g);
  double gPDF   = _hadron->pdf()->xfx(_hadron,gluon,mu2,_xB/_xp)*_xp/_xB;
  // NLO quark PDF
  double qPDF   = _hadron->pdf()->xfx(_hadron,mePartonData()[1],mu2,_xB/_xp)*_xp/_xB;
  // collinear counterterms
  // gluon
  double collg = 
    TRfact/_xp*gPDF*(2.*_xp*(1.-_xp)+(sqr(_xp)+sqr(1.-_xp))*log((1.-_xp)*_q2/_xp/mu2));
  // quark
  double collq = 
    CFfact/_xp*qPDF*(1-_xp-2./(1.-_xp)*log(_xp)-(1.+_xp)*log((1.-_xp)/_xp*_q2/mu2))+
    CFfact/_xp*(qPDF-_xp*loPDF)*(2./(1.-_xp)*log(_q2*(1.-_xp)/mu2)-1.5/(1.-_xp));
  // cacluate lepton kinematic variables
  Lorentz5Momentum  q = _pa;
  Lorentz5Momentum p1 = _pb;
  Lorentz5Momentum p2 = _pc;
  Lorentz5Momentum p3 = meMomenta()[5];
  double x1 = 2.*p1*q/(q*q);
  double x2 = 2.*p2*q/(q*q);
  double x2perp = 4.*(1-_xp)*(1.-_zp)*_zp/_xp;
  double sintheta1 = sqrt(x2perp/(sqr(x1)+x2perp));
  double sintheta2 = sqrt(x2perp/(sqr(x2)+x2perp));
  // dipoles kinematc variables
  double x = 1.-p3*p2/(p1*p3+p1*p2);
  double z = p2*p1/(p1*p3+p1*p2);
  // q -> qg term
  double real1 = CFfact*qPDF/loPDF*1./((1.-_xp)*(1.-_zp))*
    (1.+sqr(_xp)*(sqr(x2)+1.5*x2perp)*sqr(1.-sintheta2*cos(_phi))/(1.+0.5*sqr(sintheta2)));
  double dipole1 = CFfact*16.*sqr(Constants::pi)*qPDF/loPDF*1./_q2*(sqr(x)+sqr(z))/((1.-x)*(1.-z)); 
  double realq = real1-dipole1;
  // g -> q qbar term
  double real2 = TRfact*gPDF/loPDF*1./(1.-_zp)*
    (sqr(_xp)*(sqr(x1)+1.5*x2perp)*sqr(1.+sintheta1*cos(_phi))/(1.+0.5*sqr(sintheta1))+sqr(_xp)*(sqr(x2)+1.5*x2perp)*sqr(1.-sintheta2*cos(_phi))/(1.+0.5*sqr(sintheta2)));
  double dipole2 = TRfact*16.*sqr(Constants::pi)*gPDF/loPDF*1./_q2*(sqr(x)+sqr(1-x))/(1-z);
  double realg = real2-dipole2;
  // return the full result
  double wgt = virt+((collq+collg)/loPDF+realq+realg);
 
  // boost back 
  rot.invert();

  return contrib_ == 1 ? max(0.,wgt) : max(0.,-wgt);
}

void MEPP2HiggsVBFPowheg::doinit() {
  MEPP2HiggsVBF::doinit();
  // electroweak parameters
  _sinW = generator()->standardModel()->sin2ThetaW();
  _cosW = sqrt(1.-_sinW);
  _sinW = sqrt(_sinW);
  _mz2 = sqr(getParticleData(ParticleID::Z0)->mass());



  tPDPtr proton=getParticleData(2212);
  _hadron = dynamic_ptr_cast<tcBeamPtr>(proton);
  tcPDPtr gluon = getParticleData(ParticleID::g);
  double x = 0.1365;
  Energy2 Q2=320.*GeV2;
  // alphaS bits
  double aS = SM().alphaS(Q2);
  double CFfact = 4./3.*aS/Constants::twopi;
  double TRfact = 1./2.*aS/Constants::twopi;
  // test of F2
  double F2[6]={0.,0.,0.,0.,0.},F2e[6]={0.,0.,0.,0.,0.};
  for(int ix=-5;ix<=5;++ix) {
    if(ix==0) continue;
    tPDPtr parton = getParticleData(ix);
    double eq2 = sqr(double(getParticleData(ix)->iCharge())/3.);
    double loPDF = _hadron->pdf()->xfx(_hadron,parton,Q2,x)/x;
    double wsum(0.),wsqsum(0.),wgt(0.);
    unsigned int npoint(1000000);
    for(unsigned int iy=0;iy<npoint;++iy) {
      // generate z
      double rhomin = pow(1.-x,1.-power_); 
      double rho = UseRandom::rnd()*rhomin;
      double z = 1.-pow(rho,1./(1.-power_));
      double jac = rhomin/(1.-power_)*pow(1.-z,power_);
//       double z = x+UseRandom::rnd()*(1.-x);
//       double jac =(1.-x);
      double gPDF = _hadron->pdf()->xfx(_hadron,gluon ,Q2,x/z)*z/x;
      double qPDF = _hadron->pdf()->xfx(_hadron,parton,Q2,x/z)*z/x;
      double f2g = gPDF/z*TRfact*((sqr(1-z)+sqr(z))*log((1-z)/z)+8*z*(1.-z)-1.);
      double f2q = 
	loPDF/jac*(1.+CFfact*(-1.5*log(1.-x)+sqr(log(1.-x))-sqr(Constants::pi)/3.-4.5))
 	+qPDF          *CFfact/z*(3.+2.*z-(1.+z)*log(1.-z)-(1.+sqr(z))/(1.-z)*log(z))
 	+(qPDF-z*loPDF)*CFfact/z*(2.*log(1.-z)-1.5)/(1.-z);
      // NLO weight
      wgt = (f2g+f2q);
      // LO weight
      //wgt = loPDF/jac;
      wgt *= x*eq2*jac;
      wsum += wgt;
      wsqsum += sqr(wgt);
    }
    wsum /= double(npoint);
     cerr << "testing average " << wsum << " " << 0.5*(1.-sqr(x)) << "\n";
//     cerr << "testing average " << wsum << " " << (1.-x) << "\n";
    wsqsum /= double(npoint);
    F2 [abs(ix)] += wsum;
    F2e[abs(ix)] += max((wsqsum-sqr(wsum))/double(npoint),0.);
  }
  double total(0.),error(0.);
  for(unsigned int ix=1;ix<=5;++ix) {
    total += F2[ix];
    error += F2e[ix];
    F2e[ix] = sqrt(F2e[ix]);
    cerr << "testing contribution " << F2[ix] << " +/- " << F2e[ix] << "\n";
  }
  error = sqrt(error);
  cerr << "total " << total << " +/- " << error << "\n";
}



