// -*- C++ -*-
//
// SMHiggsWidthGenerator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace Herwig;

void SMHiggsWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << _widthopt << _offshell << ounit(_mw,GeV) << ounit(_mz,GeV) 
     << ounit(_gamw,GeV) << ounit(_gamz,GeV) << ounit(_qmass,GeV) 
     << ounit(_lmass,GeV) << _sw2 << _ca << _cf;
}

void SMHiggsWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _widthopt >> _offshell >> iunit(_mw,GeV) >> iunit(_mz,GeV) 
     >> iunit(_gamw,GeV) >> iunit(_gamz,GeV) >> iunit(_qmass,GeV) 
     >> iunit(_lmass,GeV) >> _sw2 >> _ca >> _cf;
}

ClassDescription<SMHiggsWidthGenerator> SMHiggsWidthGenerator::initSMHiggsWidthGenerator;
// Definition of the static class description member.

void SMHiggsWidthGenerator::Init() {

  static ClassDocumentation<SMHiggsWidthGenerator> documentation
    ("The SMHiggsWidthGenerator class calculates the running Higgs width as in"
     "hep-ph/9505211.",
     "The Higgs width was calculated as in \\cite{Seymour:1995qg}.",
     "%\\cite{Seymour:1995qg}\n"
     "\\bibitem{Seymour:1995qg}\n"
     "  M.~H.~Seymour,\n"
     "  %``The Higgs boson line shape and perturbative unitarity,''\n"
     "  Phys.\\ Lett.\\  B {\\bf 354}, 409 (1995)\n"
     "  [arXiv:hep-ph/9505211].\n"
     "  %%CITATION = PHLTA,B354,409;%%\n"
     );

  static Switch<SMHiggsWidthGenerator,unsigned int> interfaceWidthOption
    ("WidthScheme",
     "Option for the treatment of the Higss Width calculation",
     &SMHiggsWidthGenerator::_widthopt, 2, false, false);
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

  static Parameter<SMHiggsWidthGenerator,double> interfaceOffShell
    ("OffShell",
     "Number of times the width the Higgs is allowed to be off-shell",
     &SMHiggsWidthGenerator::_offshell, 10., 0.01, 100.0,
     false, false, Interface::limited);

}

bool SMHiggsWidthGenerator::accept(const ParticleData & in) const {
  return in.id()==ParticleID::h0;
}

Energy SMHiggsWidthGenerator::width(const ParticleData & in, Energy m) const {
  if(_widthopt==1) {
    return in.width();
  }
  else if(_widthopt <=3 ) {
    Energy higgswidth = ZERO;
    for (unsigned int i = 1; i < 14; ++i) higgswidth += partialWidth(m,i);
    return higgswidth;
  }
  else  
    throw Exception() << "Unknown width option in SMHiggsWidthGenerator::width()" 
		      << Exception::runerror;
}

DecayMap SMHiggsWidthGenerator::rate(const ParticleData & p) const {
  if(_mw==ZERO) return p.decaySelector();
  else            return branching(p.mass(),p);
}

DecayMap SMHiggsWidthGenerator::rate(const Particle & p) {
  return branching(p.mass(),p.data());
}

DecayMap SMHiggsWidthGenerator::branching(Energy scale, const ParticleData & p) const {
  // if not using running width return original
  if(_widthopt==1) return p.decaySelector();
  // calculate the partial widths
  vector<Energy> partial;
  Energy total(ZERO);
  for(unsigned int ix=0;ix<14;++ix) {
    partial.push_back(partialWidth(scale,ix));
    total+=partial.back();
  }
  // produce the new decay selector
  DecayMap newdm;
  for(DecaySet::const_iterator it=p.decayModes().begin();
      it!=p.decayModes().end();++it) {
    tDMPtr mode=*it;
    if(mode->orderedProducts().size()!=2||!mode->on()) continue;
    // particle antiparticle
    long id=abs(mode->orderedProducts()[0]->id());
    double br=0;
    if(mode->orderedProducts()[0]->id()==-mode->orderedProducts()[1]->id()) {
      // quarks
      if(id<=6)                            br = partial[id      ]/total;
      // leptons
      else if(id>=11&&id<=15&&(id-9)%2==0) br = partial[(id+3)/2]/total;
      // WW
      else if(id==ParticleID::Wplus)       br = partial[10      ]/total;
      // unknown mode
      else continue;
    }
    else if(mode->orderedProducts()[0]->id()==mode->orderedProducts()[1]->id()) {
      // gamma gamma
      if(id==ParticleID::Z0)         br = partial[11]/total;
      else if(id==ParticleID::gamma) br = partial[12]/total;
      else if(id==ParticleID::g)     br = partial[13]/total;
      // unknown mode
      else continue;      
    }
    // unknown mode
    else continue;
    // insert the mode into the new selector
    newdm.insert(br,mode);
  }
  return newdm;
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

Energy SMHiggsWidthGenerator::partialWidth(Energy Mh,unsigned int imode) const {
  useMe();
  using Constants::pi;
  if(Mh!=_qlast) {
    _qlast = Mh;
    Energy2 q2 = sqr(_qlast);
    // couplings
    _lambdaQCD = generator()->standardModel()->LambdaQCD(q2);
    _alphaEM   = generator()->standardModel()->alphaEM();
    _alphaS    = generator()->standardModel()->alphaS(q2);
    // QCD correction factors for H -> f fbar
    double nflavour=0.;
    for (unsigned int i =1; i <= 6; ++i) if (2.0*_qmass[i] < Mh) nflavour+=1.;
    // All calculation are being done for Monte-Carlo QCD Lambda, except Higgs width...
    // not needed in C++ as should be normal lambda
    // MC only in shower
    //double bcoeff4=(11.*_ca-10.)/(12.*pi);
    //double kfac=_ca*(67./18.-sqr(pi)/6.)-25./9.;
    //_lambdaQCD /= exp(kfac/(4.*pi*bcoeff4))/sqrt(2.);
    double k1 = 5./sqr(pi);
    double k0 = 3./(4.*sqr(pi));
    _beta0 = (11.*_ca-2.0*nflavour)/3.;
    double beta1 = (34.*sqr(_ca)-(10.*_ca+6.*_cf)*nflavour)/3.;
    _gam0 = -8.;
    double gam1 = -404./3.+40.*nflavour/9.;
    double SClog = log(sqr(Mh/_lambdaQCD));
    if(SClog<=0.)
      _cd = 1.;
    else
      _cd = 1.+(k1/k0-2.*_gam0+_gam0*beta1/sqr(_beta0)*log(SClog)+
		(_gam0*beta1-gam1*_beta0)/sqr(_beta0))/(_beta0*SClog);
    _gfermiinv = 8.*_sw2*sqr(_mw)/_alphaEM;
  }
  // output value
  Energy3 output(ZERO);
  // quark modes
  if(imode<=6) {
    Energy mf = _qmass[imode];
    double xf = sqr(mf/Mh);
    if( xf >= 0.25 ) return ZERO;
    if(_widthopt==2) {
      if (mf > _lambdaQCD) mf *= pow(log(Mh/_lambdaQCD)/log(mf/_lambdaQCD),
				     _gam0/(2.0*_beta0));
      output = _ca*Mh*sqr(mf)*pow(1.0-4.0*xf,1.5)*_cd;
    }
    else {
      output = _ca*Mh*sqr(mf)*pow(1.0-4.0*xf,1.5);
    }
  }
  // lepton modes
  else if(imode<=9) {
    Energy mf = _lmass[imode-7];
    double xf = sqr(mf/Mh);
    if (xf < 0.25) output = Mh*sqr(mf)*pow(1.0-4.0*xf,1.5);
  }
  // H->W*W*
  else if(imode==10) {
    if(_widthopt==2) {
      double xw = sqr(_mw/Mh);
      double yw = _mw*_gamw/sqr(Mh);
      output = pow<3,1>(Mh)*HwDoubleBW(xw,yw)/2.;
    }
    else {
      double xfw = sqr(_mw/Mh);
      if(2.0*_mw < Mh) output = pow<3,1>(Mh)*sqrt(1-4.0*xfw)*(1-xfw+0.75*sqr(xfw))/2.;
    }
  }
  // H->Z*Z*
  else if(imode==11) {
    if(_widthopt==2) {
      double xz = sqr(_mz/Mh);
      double yz = _mz*_gamz/sqr(Mh);
      output = pow<3,1>(Mh)*HwDoubleBW(xz,yz)/4.;
    }
    else {
      double xfz = sqr(_mz/Mh);
      if (2.0*_mz < Mh) output = pow<3,1>(Mh)*sqrt(1-4.0*xfz)*(1-xfz+0.75*sqr(xfz))/4.;
    }
  }
  // H->gamma gamma
  else if(imode==12) {
    double taut = sqr(2.0*_qmass[ParticleID::t]/Mh);
    double tauw = sqr(2.0*_mw/Mh);
    Complex ftaut = HwW2 (taut);
    Complex ftauw = HwW2 (tauw);
    double re = 4.0/3.0*(-2.0*taut*(1.0+(1.0-taut)*ftaut.real()))+(2.0+3.0*tauw*(1+(2.0-tauw)*ftauw.real()));
    double im = 4.0/3.0*(-2.0*taut*(    (1.0-taut)*ftaut.imag()))+(   3.0*tauw*(  (2.0-tauw)*ftauw.imag()));
    output = sqr(_alphaEM/pi)*pow<3,1>(Mh)*(sqr(re)+ sqr(im))/32.;
  }
  // H -> gg
  else if(imode==13) {
    double tau = sqr(2.0*_qmass[ParticleID::t]/Mh);
    Complex ftau = HwW2(tau);
    double re = 1+(1.0-tau)*ftau.real();
    double im =   (1.0-tau)*ftau.imag();
    output = sqr(_alphaS/pi)*pow<3,1>(Mh)*sqr(tau)*(sqr(re)+ sqr(im))/4.;
  }
  return output/_gfermiinv;
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

void SMHiggsWidthGenerator::doinit() {
  WidthGenerator::doinit();
  // extract W and Z mass and width
  tPDPtr w = getParticleData(ParticleID::Wplus);
  tPDPtr z = getParticleData(ParticleID::Z0);
  _mw = w->mass();
  _mz = z->mass();
  _gamw = w->width();
  _gamz = z->width();
  // quark masses
  for(int ix=1;ix<7;++ix) {
    tcPDPtr q = getParticleData(ix);
    _qmass[ix] = q->mass();
  }
  // lepton masses
  for ( int ix=0; ix<3; ++ix ) {
    tcPDPtr lepton = getParticleData(11+2*ix);
    _lmass[ix] = lepton->mass();
  }
  // sin2 theta_w
  _sw2       = generator()->standardModel()->sin2ThetaW();
  // casmirs
  double ncolour = generator()->standardModel()->Nc();
  _ca = ncolour;
  _cf = (sqr(ncolour)-1.0)/(2.0*_ca);
  // reset the width and set the limits
  tPDPtr h0=getParticleData(ParticleID::h0);
  Energy wid = width(*h0,h0->mass());
  h0->width(wid);
  h0->widthCut(_offshell*wid);
}

pair<Energy,Energy> SMHiggsWidthGenerator::width(Energy scale, const ParticleData & p) const {
  if(_widthopt==1) return make_pair(p.width(),p.width());
  // calculate the partial widths
  vector<Energy> partial;
  Energy total(ZERO);
  for(unsigned int ix=0;ix<14;++ix) {
    partial.push_back(partialWidth(scale,ix));
    total+=partial.back();
  }
  // sum for the on modes
  Energy partialon(ZERO);
  for(DecaySet::const_iterator it=p.decayModes().begin();it!=p.decayModes().end();++it) {
    tDMPtr mode=*it;
    if(!mode->on()||mode->orderedProducts().size()!=2) continue;
    // particle antiparticle
    long id=abs(mode->orderedProducts()[0]->id());
    if(mode->orderedProducts()[0]->id()==-mode->orderedProducts()[1]->id()) {
      // quarks
      if(id<=6)                            partialon += partial[id      ];
      // leptons
      else if(id>=11&&id<=15&&(id-9)%2==0) partialon += partial[(id+3)/2];
      // WW
      else if(id==ParticleID::Wplus)       partialon += partial[10      ];
      // unknown mode
      else continue;
    }
    else if(mode->orderedProducts()[0]->id()==mode->orderedProducts()[1]->id()) {
      // gamma gamma
      if(id==ParticleID::Z0)         partialon += partial[11];
      else if(id==ParticleID::gamma) partialon += partial[12];
      else if(id==ParticleID::g)     partialon += partial[13];
      // unknown mode
      else continue;      
    }
    // unknown mode
    else continue;
  }
  return make_pair(partialon,total);
}
