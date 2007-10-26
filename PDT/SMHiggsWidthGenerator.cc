// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsWidthGenerator class.
//

#include "SMHiggsWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace Herwig;

void SMHiggsWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << _widthopt << _branchingopt << ounit(_mw,GeV) << ounit(_mz,GeV) 
     << ounit(_gamw,GeV) << ounit(_gamz,GeV) << ounit(_qmass,GeV) 
     << ounit(_lmass,GeV);
}

void SMHiggsWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _widthopt >> _branchingopt >> iunit(_mw,GeV) >> iunit(_mz,GeV) 
     >> iunit(_gamw,GeV) >> iunit(_gamz,GeV) >> iunit(_qmass,GeV) 
     >> iunit(_lmass,GeV);
}

ClassDescription<SMHiggsWidthGenerator> SMHiggsWidthGenerator::initSMHiggsWidthGenerator;
// Definition of the static class description member.

void SMHiggsWidthGenerator::Init() {

  static ClassDocumentation<SMHiggsWidthGenerator> documentation
    ("The SMHiggsWidthGenerator class calculates the running Higgs width as in"
     "hep-ph/9505211.",
     "The Higgs width was calculated as in \\cite{Seymour:1995qg}.",
     "\\bibitem{Seymour:1995qg} M.~H.~Seymour,\n"
     "Phys.\\ Lett.\\  B {\\bf 354} (1995) 409 [arXiv:hep-ph/9505211].\n"
     "%%CITATION = PHLTA,B354,409;%%\n");

  static Switch<SMHiggsWidthGenerator,unsigned int> interfaceWidthOption
    ("WidthScheme",
     "Option for the treatment of the Higss Width calculation",
     &SMHiggsWidthGenerator::_widthopt, 1, false, false);
  static SwitchOption interfaceFixedWidth
    (interfaceWidthOption,
     "Fixed",
     "Fixed Higgs width, taken from ThePEGParticles.in",
     1);
  static SwitchOption interfaceNLLWidth
    (interfaceWidthOption,
     "NLLcorrected",
     "NLL corrected Higgs width (a-la FORTRAN HERWIG)",
     2);
  static SwitchOption interfaceLOWidthOption
    (interfaceWidthOption,
     "LO",
     "LO Higgs width (formula taken from The \"Higgs Hunter's Guide\")",
     3);

  static Switch<SMHiggsWidthGenerator,unsigned int> interfaceBranchingOption
    ("Branching",
     "Option to switch on/off branchings in the total Higgs width calculation",
     &SMHiggsWidthGenerator::_branchingopt, 3, false, false);
  static SwitchOption interfaceFermionOnlyBranchings
    (interfaceBranchingOption,
     "Fermion",
     "Fermion branchings in the full Higgs width",
     1);
  static SwitchOption interfacePlusWWZZBranchings
    (interfaceBranchingOption,
     "FermionWZ",
     "Fermion and WW/ZZ branchings in the full Higgs width",
     2);
  static SwitchOption interfacePlusGammaBranching
    (interfaceBranchingOption,
     "FermionWZGamma",
     "Fermion 2gamma, and WW/ZZ branchings in the full Higgs width",
     3);
  static SwitchOption interfacePlausGluonBranching
    (interfaceBranchingOption,
     "FermionWZGammaGluon",
     "Fermion 2gamma, 2gluons, and WW/ZZ branchings in the full Higgs width",
     4);
}

bool SMHiggsWidthGenerator::accept(const ParticleData & in) const {
  return in.id()==ParticleID::h0;
}

Energy SMHiggsWidthGenerator::width(const ParticleData & in, Energy m) const {
  switch (_widthopt) {
  case 1:
    return in.width();
  case 2: 
    return calcNLLRunningWidth(m);
  case 3: 
    return  calcLORunningWidth(m);
  default: 
    throw Exception() << "Unknown width option in SMHiggsWidthGenerator::width()" 
		      << Exception::runerror;
  }
}

DecayMap SMHiggsWidthGenerator::rate(const ParticleData & p) const {
  return p.decaySelector();
}

// Taken from HERWIG 6510.
Complex SMHiggsWidthGenerator::HwW2(double tau) const {
  using Constants::pi;
  if (tau > 1.0) {
    return sqr(asin(1.0/sqrt(tau)));
  } 
  else if (tau < 1.0) {
    double FNsqr = sqrt(1-tau);
    double FNlog = log((1+FNsqr)/(1-FNsqr));
    return Complex(-0.25 * (sqr(FNlog)-pi*pi),0.5*pi*FNlog);
  } 
  else {
    return sqr(0.5*pi);
  }
}

// Taken from HERWIG 6510 with some simplifications.
Energy SMHiggsWidthGenerator::calcNLLRunningWidth(Energy Mh) const {
  using Constants::pi;
  Energy2 q2 = sqr(Mh);
  Energy QCDLambda = generator()->standardModel()->LambdaQCD(q2);
//  double alphaEM   = generator()->standardModel()->alphaEM(q2);
  double alphaEM   = generator()->standardModel()->alphaEM();
  double alphaS    = generator()->standardModel()->alphaS(q2);
  double sw2       = generator()->standardModel()->sin2ThetaW();
  _partialwidths=vector<Energy>(14,0.*GeV);

  double ncolour = generator()->standardModel()->Nc();
  double Ca = ncolour;
  double Cf = (sqr(ncolour)-1.0)/(2.0*Ca);
  const unsigned int minflavour = 1;
  const unsigned int maxflavour = 6;

// All calculation are being done for Monte-Carlo QCD Lambda, except Higgs width...
  double bcoeff4=(11.*Ca-10.)/(12.*pi);
  double kfac=Ca*(67./18.-sqr(pi)/6.)-25./9.;
  QCDLambda /= exp(kfac/(4.*pi*bcoeff4))/sqrt(2.);

// H->fermion pair
  double nflavour=minflavour-1;
  for (unsigned int i = minflavour; i <= maxflavour; ++i) {
    if (2.0*_qmass[i] < Mh) nflavour+=1.;
  }
  double k1 = 5./sqr(pi);
  double k0 = 3./(4.*sqr(pi));
  double beta0 = (11.*Ca-2.0*nflavour)/3.;
  double beta1 = (34.*sqr(Ca)-(10.*Ca+6.*Cf)*nflavour)/3.;
  double gam0 = -8.;
  double gam1 = -404./3.+40.*nflavour/9.;
  double SClog = log(sqr(Mh/QCDLambda));
  double Cd = 1.+(k1/k0-2.*gam0+gam0*beta1/sqr(beta0)*log(SClog)+(gam0*beta1-gam1*beta0)/sqr(beta0))/(beta0*SClog);
  Energy2 GFermiINV = 8.*sw2*sqr(_mw)/alphaEM;

// quarks: partialw[1-6]
  for (unsigned int i = minflavour; i <= maxflavour; ++i) {
    Energy mf = _qmass[i];
    double xf = sqr(mf/Mh);
    if (mf > QCDLambda) mf *= pow(log(Mh/QCDLambda)/log(mf/QCDLambda),gam0/(2.0*beta0));
    if (xf < 0.25) _partialwidths[i] = Ca*Mh*sqr(mf)*pow(1.0-4.0*xf,1.5)*Cd/GFermiINV;
  }

// leptons: partialw[7-9]
  for (unsigned int i = 0; i < 3; ++i) {
    Energy mf = _lmass[i];
    double xf = sqr(mf/Mh);
    if (xf < 0.25) _partialwidths[7+i] = Mh*sqr(mf)*pow(1.0-4.0*xf,1.5)/GFermiINV;
  }

// H->W*W*/Z*Z*: partialw[10,11]
  if (_branchingopt > 1) {
    double xw = sqr(_mw/Mh);
    double xz = sqr(_mz/Mh);
    double yw = _mw*_gamw/sqr(Mh);
    double yz = _mz*_gamz/sqr(Mh);
    _partialwidths[10] = Mh*Mh*Mh*HwDoubleBW(xw,yw)/2./GFermiINV;
    _partialwidths[11] = Mh*Mh*Mh*HwDoubleBW(xz,yz)/4./GFermiINV;
  }

// H->gamma,gamma: partialw[12]
  if (_branchingopt > 2) {
    double taut = sqr(2.0*_qmass[ParticleID::t]/Mh);
    double tauw = sqr(2.0*_mw/Mh);
    Complex ftaut = HwW2 (taut);
    Complex ftauw = HwW2 (tauw);
    double re = 4.0/3.0*(-2.0*taut*(1.0+(1.0-taut)*ftaut.real()))+(2.0+3.0*tauw*(1+(2.0-tauw)*ftauw.real()));
    double im = 4.0/3.0*(-2.0*taut*(    (1.0-taut)*ftaut.imag()))+(   3.0*tauw*(  (2.0-tauw)*ftauw.imag()));
    _partialwidths[12] = sqr(alphaEM/pi)*Mh*Mh*Mh*(sqr(re)+ sqr(im))/32./GFermiINV;
  }

// H->gluon,gluon: partialw[13]
  if (_branchingopt > 3) {
    double tau = sqr(2.0*_qmass[ParticleID::t]/Mh);
    Complex ftau = HwW2(tau);
    double re = 1+(1.0-tau)*ftau.real();
    double im =   (1.0-tau)*ftau.imag();
    _partialwidths[13] = sqr(alphaS/pi)*Mh*Mh*Mh*sqr(tau)*(sqr(re)+ sqr(im))/4./GFermiINV;
  }

  Energy higgswidth = Energy();
  for (unsigned int i = 1; i < 14; ++i) {
    higgswidth += _partialwidths[i];
  }
  return higgswidth;
}

// Taken from HERWIG 6510 with simplifications.
Energy SMHiggsWidthGenerator::calcLORunningWidth(Energy Mh) const {
  Energy2 q2 = sqr(Mh);
  Energy QCDLambda = generator()->standardModel()->LambdaQCD(q2);
//  double alphaEM   = generator()->standardModel()->alphaEM(q2);
  double alphaEM   = generator()->standardModel()->alphaEM();
  double alphaS    = generator()->standardModel()->alphaS(q2);
  double sw2       = generator()->standardModel()->sin2ThetaW();
  _partialwidths = vector<Energy>(14,0.*GeV);

  double ncolour = generator()->standardModel()->Nc();
  double Ca = ncolour;
  const unsigned int minflavour = 1;
  const unsigned int maxflavour = 6;
  double pi = Constants::pi;

// All calculation are being done for Monte-Carlo QCD Lambda, except Higgs width...
  double bcoeff4=(11.*Ca-10.)/(12.*pi);
  double kfac(Ca*(67./18.-sqr(pi)/6.)-25./9.);
  QCDLambda /= exp(kfac/(4.*pi*bcoeff4))/sqrt(2.);

// H->fermion pair
  double nflavour = minflavour-1;
  for (unsigned int i = minflavour; i <= maxflavour; ++i) {
    if (2.0*_qmass[i] < Mh) nflavour+=1.;
  }
  Energy2 GFermiINV = 8.*sw2*sqr(_mw)/alphaEM;

// quarks: partialw[1-6]
  for (unsigned int i = minflavour; i <= maxflavour; ++i) {
    Energy mf = _qmass[i];
    double xf = sqr(mf/Mh);
    if (xf < 0.25) _partialwidths[i] = Ca*Mh*sqr(mf)*pow(1.0-4.0*xf,1.5)/GFermiINV;
  }

// leptons: partialw[7-9]
  for (unsigned int i = 0; i < 3; ++i) {
    Energy mf = _lmass[i];
    double xf = sqr(mf/Mh);
    if (xf < 0.25) _partialwidths[7+i] = Mh*sqr(mf)*pow(1.0-4.0*xf,1.5)/GFermiINV;
  }

// H->W*W*/Z*Z*: partialw[10,11]
  if (_branchingopt > 1) {
    double xfw = sqr(_mw/Mh);
    double xfz = sqr(_mz/Mh);
    if (2.0*_mw < Mh) _partialwidths[10] = Mh*Mh*Mh*sqrt(1-4.0*xfw)*(1-xfw+0.75*sqr(xfw))/2./GFermiINV;
    if (2.0*_mz < Mh) _partialwidths[11] = Mh*Mh*Mh*sqrt(1-4.0*xfz)*(1-xfz+0.75*sqr(xfz))/4./GFermiINV;
  }

// H->gamma,gamma: partialw[12]
  if (_branchingopt > 2) {
    double taut = sqr(2.0*_qmass[ParticleID::t]/Mh);
    double tauw = sqr(2.0*_mw/Mh);
    Complex ftaut = HwW2(taut);
    Complex ftauw = HwW2(tauw);
    double re = 4.0/3.0*(-2.0*taut*(1.0+(1.0-taut)*ftaut.real()))+(2.0+3.0*tauw*(1+(2.0-tauw)*ftauw.real()));
    double im = 4.0/3.0*(-2.0*taut*(    (1.0-taut)*ftaut.imag()))+(   3.0*tauw*(  (2.0-tauw)*ftauw.imag()));
    _partialwidths[12] = sqr(alphaEM/pi)*Mh*Mh*Mh*(sqr(re)+ sqr(im))/32./GFermiINV;
  }

// H->gluon,gluon: partialw[13]
  if (_branchingopt > 3) {
    double tau = sqr(2.0*_qmass[ParticleID::t]/Mh);
    Complex ftau = HwW2 (tau);
    double re = 1+(1.0-tau)*ftau.real();
    double im =   (1.0-tau)*ftau.imag();
    _partialwidths[13] = sqr(alphaS/pi)*Mh*Mh*Mh*sqr(tau)*(sqr(re)+ sqr(im))/4./GFermiINV;
  }

  Energy higgswidth = Energy();
  for (unsigned int i = 1; i < 13; ++i) {
    higgswidth += _partialwidths[i];
  }
  return higgswidth;
}

// Taken from HERWIG 6510.
double SMHiggsWidthGenerator::HwDoubleBW(double x, double y) const {
// Calculate the Double Breit-Wigner Integral
//  x=(mw/mh)**2, y=mw*gw/mh**2
  double limit = 0.425;
  double nbin = 25;
  double itgerl = 0.0;
  if (y < 0.0) return itgerl;
  if (x > limit) {
// Direct Integration
    double fac1 = 0.25/nbin;
    for (unsigned int i = 0; i < nbin; ++i) {
      double x1 = (i+0.5)*fac1;
      double fac2 = (sqr(1-sqrt(x1))-x1)/nbin;
      double sq = 1.0;
      int j = 0;
      while (j < nbin && 0.0 < sq) {
        double x2 = (j+0.5)*fac2+x1;
        sq = 1.0+x1*x1+x2*x2-2*(x1+x2+x1*x2);
        itgerl += 2.0*(sqr(1-x1-x2)+8.0*x1*x2)*sqrt(sq)/(sqr(x1-x)+y*y)*y/(sqr(x2-x)+y*y)*y*fac1*fac2;
        ++j;
      }
    }
  } else {
// Integration using tan theta substitution
    double th1low = atan ((0.0-x)/y);
    double th1high = atan ((1.0-x)/y);
    double fac1 = (th1high-th1low)/nbin;
    for (unsigned int i = 0; i < nbin; ++i) {
      double th1 = (i+0.5)*fac1+th1low;
      double x1 = y*tan(th1)+x;
      double x2max = min(x1,sqr(1-sqrt(x1)));
      double th2low = atan ((0-x)/y);
      double th2high = atan ((x2max-x)/y);
      double fac2 = (th2high-th2low)/nbin;
      double sq = 1.0;
      int j = 0;
      while (j < nbin && 0.0 < sq) {
        double th2 = (j+0.5)*fac2+th2low;
        double x2 = y*tan(th2)+x;
        double sq = 1.0+x1*x1+x2*x2-2*(x1+x2+x1*x2);
        itgerl += 2.0*(sqr(1-x1-x2)+8*x1*x2)*sqrt(sq)*fac1*fac2;
        ++j;
      }
    }
  }
  itgerl *= 1/sqr(Constants::pi);
  return itgerl;
}

void SMHiggsWidthGenerator::doinit() throw(InitException) {
  WidthGenerator::doinit();
  // extract W and Z mass and width
  tPDPtr w = getParticleData(ParticleID::Wplus);
  tPDPtr z = getParticleData(ParticleID::Z0);
  _mw = w->mass();
  _mz = z->mass();
  _gamw = w->width();
  _gamz = z->width();
  for(unsigned int ix=1;ix<7;++ix) {
    tcPDPtr q = getParticleData(ix);
    _qmass[ix] = q->mass();
  }
  for(unsigned int ix=0;ix<3;++ix) {
    tcPDPtr lepton = getParticleData(11+2*ix);
    _lmass[ix] = lepton->mass();
  }
}
