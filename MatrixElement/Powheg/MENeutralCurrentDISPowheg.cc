// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MENeutralCurrentDISPowheg class.
//

#include "MENeutralCurrentDISPowheg.h"
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

MENeutralCurrentDISPowheg::MENeutralCurrentDISPowheg() 
  : scaleOpt_(1),  muF_(100.*GeV), scaleFact_(1.), contrib_(1), power_(0.1)
{}

int MENeutralCurrentDISPowheg::nDim() const {
  return MENeutralCurrentDIS::nDim()+1;
}

bool MENeutralCurrentDISPowheg::generateKinematics(const double * r) {
  // Born kinematics
  if(!MENeutralCurrentDIS::generateKinematics(r)) return false;
  // hadron and momentum fraction
  if(HadronMatcher::Check(*lastParticles().first->dataPtr())) {
    _hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr());
    _xB = lastX1();
  }
  else {
    _hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr());
    _xB = lastX2();
  }
  // Q2
  _q2 = -(meMomenta()[0]-meMomenta()[2]).m2();
  // xp
  int ndim=nDim();
  double rhomin = pow(1.-_xB,1.-power_); 
  double rho = r[ndim-1]*rhomin;
  _xp = 1.-pow(rho,1./(1.-power_));
  jac_ = rhomin/(1.-power_)*pow(1.-_xp,power_);
  jacobian(jacobian()*jac_);
  return true; 
}

IBPtr MENeutralCurrentDISPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MENeutralCurrentDISPowheg::fullclone() const {
  return new_ptr(*this);
}

void MENeutralCurrentDISPowheg::persistentOutput(PersistentOStream & os) const {
  os << ounit(muF_,GeV) << scaleFact_ << scaleOpt_ << contrib_
     << _sinW << _cosW << ounit(_mz2,GeV2) << power_;
}

void MENeutralCurrentDISPowheg::persistentInput(PersistentIStream & is, int) {
  is >> iunit(muF_,GeV) >> scaleFact_ >> scaleOpt_ >> contrib_
     >> _sinW >> _cosW >> iunit(_mz2,GeV2) >> power_;
}

ClassDescription<MENeutralCurrentDISPowheg> 
MENeutralCurrentDISPowheg::initMENeutralCurrentDISPowheg;
// Definition of the static class description member.

void MENeutralCurrentDISPowheg::Init() {

  static ClassDocumentation<MENeutralCurrentDISPowheg> documentation
    ("The MENeutralCurrentDISPowheg class implements the NLO matrix element"
     " for neutral current DIS in the Powheg scheme.");

  static Switch<MENeutralCurrentDISPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MENeutralCurrentDISPowheg::contrib_, 1, false, false);
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

  static Switch<MENeutralCurrentDISPowheg,unsigned int> interfaceScaleOption
    ("ScaleOption",
     "Option for the choice of factorization (and renormalization) scale",
     &MENeutralCurrentDISPowheg::scaleOpt_, 1, false, false);
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

  static Parameter<MENeutralCurrentDISPowheg,Energy> interfaceFactorizationScale
    ("FactorizationScale",
     "Value to use in the event of a fixed factorization scale",
     &MENeutralCurrentDISPowheg::muF_, GeV, 100.0*GeV, 1.0*GeV, 500.0*GeV,
     true, false, Interface::limited);

  static Parameter<MENeutralCurrentDISPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before Q2 if using a running scale",
     &MENeutralCurrentDISPowheg::scaleFact_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<MENeutralCurrentDISPowheg,double> interfaceSamplingPower
    ("SamplingPower",
     "Power for the sampling of xp",
     &MENeutralCurrentDISPowheg::power_, 0.6, 0.0, 1.,
     false, false, Interface::limited);

}

Energy2 MENeutralCurrentDISPowheg::scale() const {
  return scaleOpt_ == 1 ? 
    sqr(scaleFact_)*MENeutralCurrentDIS::scale() : sqr(scaleFact_*muF_);
}

CrossSection MENeutralCurrentDISPowheg::dSigHatDR() const {
  return NLOWeight()*MENeutralCurrentDIS::dSigHatDR();
}

double MENeutralCurrentDISPowheg::NLOWeight() const {
  // If only leading order is required return 1:
  if(contrib_==0) return 1./jac_;
  // scale and prefactors
  Energy2 mu2(scale());
  double aS = SM().alphaS(mu2);
  double CFfact = 4./3.*aS/Constants::twopi;
  double TRfact = 1./2.*aS/Constants::twopi;
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
  // calculate the A coefficient for the real pieces
  double a(A());
  // cacluate lepton kinematic variables
  Lorentz5Momentum q = meMomenta()[0]-meMomenta()[2];
  double  yB = (q*meMomenta()[1])/(meMomenta()[0]*meMomenta()[1]);
  double l = 2./yB-1.;
  // q -> qg term
  double realq = CFfact/_xp/(1.+a*l+sqr(l))*qPDF/loPDF*
    (2.+2.*sqr(l)-_xp+3.*_xp*sqr(l)+a*l*(2.*_xp+1.));
  // g -> q qbar term
  double realg =-TRfact/_xp/(1.+a*l+sqr(l))*gPDF/loPDF*
    ((1.+sqr(l)+2.*(1.-3.*sqr(l))*_xp*(1.-_xp))
     +2.*a*l*(1.-2.*_xp*(1.-_xp)));
  // return the full result
  double wgt = virt+((collq+collg)/loPDF+realq+realg); 
  //   double f2g = gPDF/_xp*TRfact*((sqr(1-_xp)+sqr(_xp))*log((1-_xp)/_xp)+
  // 				8*_xp*(1.-_xp)-1.);
  //   double f2q = 
  //     loPDF/jac_*(1.+CFfact*(-1.5*log(1.-_xB)+sqr(log(1.-_xB))
  // 			   -sqr(Constants::pi)/3.-4.5))
  //     +qPDF            *CFfact/_xp*(3.+2.*_xp-(1.+_xp)*log(1.-_xp)
  // 				  -(1.+sqr(_xp))/(1.-_xp)*log(_xp))
  //     +(qPDF-_xp*loPDF)*CFfact/_xp*(2.*log(1.-_xp)/(1.-_xp)-1.5/(1.-_xp));
  //   double wgt = (f2g+f2q)/loPDF;
  return contrib_ == 1 ? max(0.,wgt) : max(0.,-wgt);
}

void MENeutralCurrentDISPowheg::doinit() {
  MENeutralCurrentDIS::doinit();
  // electroweak parameters
  _sinW = generator()->standardModel()->sin2ThetaW();
  _cosW = sqrt(1.-_sinW);
  _sinW = sqrt(_sinW);
  _mz2 = sqr(getParticleData(ParticleID::Z0)->mass());
}

double MENeutralCurrentDISPowheg::A() const {
  double output;
  double factZ = (gammaZOption()==0||gammaZOption()==2) ? 
    0.25*double(_q2/(_q2+_mz2))/_sinW/_cosW : 0;
  double factG = (gammaZOption()==0||gammaZOption()==1) ? 
    1 : 0;
  double cvl,cal,cvq,caq;
  if(abs(mePartonData()[0]->id())%2==0) {
    cvl = generator()->standardModel()->vnu()*factZ+
          generator()->standardModel()->enu()*factG;
    cal = generator()->standardModel()->anu()*factZ;
  }
  else {
    cvl = generator()->standardModel()->ve()*factZ+
          generator()->standardModel()->ee()*factG;
    cal = generator()->standardModel()->ae()*factZ;
  }
  if(abs(mePartonData()[1]->id())%2==0) {
    cvq = generator()->standardModel()->vu()*factZ+
          generator()->standardModel()->eu()*factG;
    caq = generator()->standardModel()->au()*factZ;
  }
  else {
    cvq = generator()->standardModel()->vd()*factZ+
          generator()->standardModel()->ed()*factG;
    caq = generator()->standardModel()->ad()*factZ;
  }
  output = 8.*cvl*cal*cvq*caq/(sqr(cvl)+sqr(cal))/(sqr(cvq)+sqr(caq));
  if(mePartonData()[1]->id()<0) output *= -1.;
  if(mePartonData()[0]->id()<0) output *= -1;
  return output;
}
