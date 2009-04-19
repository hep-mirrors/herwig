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
  return MENeutralCurrentDIS::nDim()+2;
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
  double rho = r[ndim-2]*rhomin;
  _xp = 1.-pow(rho,1./(1.-power_));
  jac_ = rhomin/(1.-power_)*pow(1.-_xp,power_);
  jacobian(jacobian()*jac_);
  // zp
  _zp = r[ndim-1];
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
  if(contrib_==0) return 1.;
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
//   double realq = 2.*CFfact/_xp/(1.+a*l+sqr(l))*qPDF/loPDF*
//     (1.+sqr(l)-_xp*_zp+3.*_xp*_zp*sqr(l)+a*l*(_xp+_zp));
  //   double realq = 2.*CFfact/_xp*qPDF/loPDF*(1.+3.*_xp*_zp);
  // g -> q qbar term
  double realg =-TRfact/_xp/(1.+a*l+sqr(l))*gPDF/loPDF*
    ((1.+sqr(l)+2.*(1.-3.*sqr(l))*_xp*(1.-_xp))
     +2.*a*l*(1.-2.*_xp*(1.-_xp)));
//   double realg =-2.*TRfact/_xp/(1.+a*l+sqr(l))*gPDF/loPDF*
//     (_zp*(1.+sqr(l)+2.*(1.-3.*sqr(l))*_xp*(1.-_xp))
//      +a*l*(1.-2.*_xp+2.*sqr(_xp)));
  //   double realg =-2.*TRfact/_xp*gPDF/loPDF*(_zp*(1.-6.*_xp*(1.-_xp)));
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

  // kinematically variables for the real emission terms
//   double x1 = -1./_xp;
//   double x2 = 1.-(1.-_zp)/_xp;
//   double x3 = 2.+x1-x2;
//   double xT = sqrt(4.*(1.-_xp)*(1.-_zp)*_zp/_xp);
//   double cos2 = x2/sqrt(sqr(x2)+sqr(xT));
//   double cos3 = x3/sqrt(sqr(x3)+sqr(xT));
//   double R2 = 0.5*(3.*sqr(cos2)-1.+2.*a*cos2*l+3.*sqr(l)-sqr(l*cos2))/(1.+a*l+sqr(l));
//   double qme = 1.+sqr(_xp)*(sqr(x2)+sqr(xT))*R2;
//   qme *= CFfact/_xp/(1.-_xp)/(1.-_zp);
//   double dipoleQG = CFfact*(-2.+(sqr(x1)+sqr(x2))*sqr(_xp)/(1.-_xp)/(1.-_zp))/_xp;
//   double testq = (qme-dipoleQG)*qPDF/loPDF;
//   double R3 = 0.5*(3.*sqr(cos3)-1.-2.*a*cos3*l+3.*sqr(l)-sqr(l*cos3))/(1.+a*l+sqr(l));
//   double gme = TRfact*_xp/(1.-_zp)*((sqr(x2)+sqr(xT))*R2+(sqr(x3)+sqr(xT))*R3);
//   double dipoleQQ = TRfact/_xp/(1.-_zp)*(1.-2.*_xp+2.*sqr(_xp));
//   double realg=(gme-dipoleQQ)*gPDF/loPDF;

//   double testg2 = 2.*TRfact/_xp/(1.+a*l+sqr(l))*gPDF*
//     (-_zp*(1+sqr(l))+6*sqr(l)*_xp*_zp*(1-_xp)-2*_xp*_zp*(1-_xp)+a*l*(-1+2*_xp*(1-_xp)));


//   cerr << "QQBAR TEST RATIO WITH SAME PDF " 
//        << (realg-testg/loPDF)/(realg+testg/loPDF) << "\n";


//   cerr << "testing cos3 " << -(_xp-_zp)/(-_xp-_zp+2*_xp*_zp) << " " << cos3 << "\n";




//   cerr << "testing A "
//        <<  TRfact*_xp/(1.-_zp)*gPDF/loPDF*
//     ((sqr(x2)+sqr(xT))*(0.5*(3.*sqr(cos2)-1.+2.*a*cos2*l+3.*sqr(l)-sqr(l*cos2))/(1.+a*l+sqr(l)))+
//      (sqr(x3)+sqr(xT))*(0.5*(3.*sqr(cos3)-1.-2.*a*cos3*l+3.*sqr(l)-sqr(l*cos3))/(1.+a*l+sqr(l))))
//        << " " << gme << "\n";
//   cerr << "testing B "
//        << TRfact/_xp/(1.-_zp)*(1.-2.*_xp+2.*sqr(_xp))*gPDF/loPDF
//        << " " << dipoleQQ*gPDF/loPDF << "\n";

//   cerr << "testing C "
//        <<  (TRfact*_xp/(1.-_zp)*gPDF/loPDF*
//     ((sqr(x2)+sqr(xT))*(0.5*(3.*sqr(cos2)-1.+2.*a*cos2*l+3.*sqr(l)-sqr(l*cos2))/(1.+a*l+sqr(l)))+
//      (sqr(x3)+sqr(xT))*(0.5*(3.*sqr(cos3)-1.-2.*a*cos3*l+3.*sqr(l)-sqr(l*cos3))/(1.+a*l+sqr(l)))))
//     -(TRfact/_xp/(1.-_zp)*(1.-2.*_xp+2.*sqr(_xp))*gPDF/loPDF) << "\n";
//   cerr << "testing D "
//        <<  TRfact/(1.-_zp)*gPDF/loPDF/_xp*
//     ((sqr(_xp)*
//       ((sqr(x2)+sqr(xT))*(0.5*(3.*sqr(cos2)-1.+2.*a*cos2*l+3.*sqr(l)-sqr(l*cos2))/(1.+a*l+sqr(l)))+
//        (sqr(x3)+sqr(xT))*(0.5*(3.*sqr(cos3)-1.-2.*a*cos3*l+3.*sqr(l)-sqr(l*cos3))/(1.+a*l+sqr(l)))))
//      -(1.-2.*_xp+2.*sqr(_xp))) << "\n";


//  cerr << "testing M "
//       <<  0.5*TRfact/(1.-_zp)*gPDF/loPDF/_xp/(1.+a*l+sqr(l))*
//    (-4.*(2.*_xp+1.-2.*sqr(_xp))*_zp*(1.-_zp)

//     +2.*a*l*(-2*_xp+1+2*sqr(_xp))*(2*_zp-1)

//     +sqr(l)*(+24*_xp*_zp+4*sqr(_zp)-4*_zp+24*sqr(_xp)*sqr(_zp)-24*sqr(_xp)*_zp-24*_xp*sqr(_zp))
   

//    -2.*(1.-2.*_xp+2.*sqr(_xp))*(a*l)) << "\n";

//  cerr << "testing N "
//       <<  0.5*TRfact*gPDF/loPDF/_xp/(1.+a*l+sqr(l))*
//    (-4.*(2.*_xp+1.-2.*sqr(_xp))*_zp

//     -4.*a*l*(-2*_xp+1+2*sqr(_xp))

//     +sqr(l)*(+24*_xp*(1.-_xp)-4)*_zp
   

//    ) << "\n";

//  cerr << "testing 0 "
//       <<  2.*TRfact*gPDF/loPDF/_xp/(1.+a*l+sqr(l))*
//    (-(2.*_xp+1.-2.*sqr(_xp))*_zp

//     -a*l*(-2*_xp+1+2*sqr(_xp))

//     +sqr(l)*(+6*_xp*(1.-_xp)-1)*_zp
   

//    ) << "\n";
//   exit(0);

//   //


//   // tests of the matrix elements
//   // the momenta
//   Energy Q(sqrt(_q2));
//   Lorentz5Momentum p1( 0.5*Q*xT,  ZERO, -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
//   Lorentz5Momentum p2(-0.5*Q*xT,  ZERO , -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
//   Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
//   Lorentz5Momentum qq(ZERO,ZERO,-Q,ZERO);
//   // tests of the real q -> q g channel
//   // subtracted with same pdf in dipole and me
// //   double testq = 2.*CFfact/_xp/(1.+a*l+sqr(l))*qPDF/loPDF*
// //     (1.+sqr(l)-_xp*_zp+3.*_xp*_zp*sqr(l)+a*l*(_xp+_zp));
//   cerr << "QG TEST RATIO WITH SAME PDF " 
//        << (qme-dipoleQG*qPDF/loPDF)/testq << "\n";
//   // momenta for the FS dipole
//   Lorentz5Momentum pj(p2),pa(p0),pi(p1);
//   // computations of the FS dipole
//   double xFS = (pi*pa+pj*pa-pi*pj)/(pi*pa+pj*pa);
//   double zFS = (pi*pa)/(pi*pa+pj*pa);
//   // in terms of momenta
//   InvEnergy2 DFSA =  0.5/(pi*pj)/xFS*8.*Constants::pi*aS*4./3.*(2./(1.-zFS+1.-xFS)-1.-zFS);
//   // in terms of x etc
//   InvEnergy2 DFSB =  8.*Constants::pi*aS*4./3./_q2/(1.+x1)*(2.*sqr(x1)/(x1+x2)-(1.+2.*x1-x2));
//   cerr << "FS dipoles " << (DFSA-DFSB)/(DFSA+DFSB) << "\n";
//   // computations of the IS dipole
//   Lorentz5Momentum pk(p1);
//   pi=p2;
//   double xIS = (pk*pa+pi*pa-pi*pk)/(pk*pa+pi*pa);
//   double uIS = (pi*pa)/(pi*pa+pk*pa);
//   // in terms of momenta
//   InvEnergy2 DISA =  0.5/(pa*pi)/xIS*8.*Constants::pi*aS*4./3.*(2./(1.-xIS+uIS)-(1.+xIS));
//   // in terms of x etc
//   InvEnergy2 DISB =-8.*Constants::pi*aS*4./3./_q2/(1.-x2)*(2.*sqr(x1)/(x1+x2)+(1.-x1));
//   cerr << "IS dipoles " << (DISA-DISB)/(DISA+DISB) << "\n";
//   double qgtest = 0.25*(DISA+DFSA)/sqr(Constants::twopi)*_q2/_xp;
//   cerr << "full dipole " << (qgtest-dipoleQG)/(qgtest+dipoleQG) << "\n";



//   // tests of the real g -> q qbar channel
//   // dipole
//   // in terms of momenta
//   InvEnergy2 DQQA =  0.5/(pa*pi)/xIS*8.*Constants::pi*aS*0.5*(1.-2.*xIS*(1.-xIS));
//   // in terms of x etc
//   InvEnergy2 DQQB =-8.*Constants::pi*aS*0.5/_q2/(1.-x2)/x1*(1.+sqr(1.+x1));
//   cerr << "testing QQ dipole " << (DQQA-DQQB)/(DQQA+DQQB) << "\n";
//   double qqtest = 0.25*DQQA/sqr(Constants::twopi)*_q2/_xp;
//   cerr << "full dipole " << (qqtest-dipoleQQ)/(qqtest+dipoleQQ) << "\n";




//   exit(0);
