// -*- C++ -*-
//
// SMHiggsWidthGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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

using namespace Herwig;

IBPtr SMHiggsWidthGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr SMHiggsWidthGenerator::fullclone() const {
  return new_ptr(*this);
}

void SMHiggsWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << widthopt_ << offshell_ << ounit(mw_,GeV) << ounit(mz_,GeV) 
     << ounit(gamw_,GeV) << ounit(gamz_,GeV) << ounit(qmass_,GeV) 
     << ounit(lmass_,GeV) << sw2_ << ca_ << cf_ << locMap_;
}

void SMHiggsWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> widthopt_ >> offshell_ >> iunit(mw_,GeV) >> iunit(mz_,GeV) 
     >> iunit(gamw_,GeV) >> iunit(gamz_,GeV) >> iunit(qmass_,GeV) 
     >> iunit(lmass_,GeV) >> sw2_ >> ca_ >> cf_ >> locMap_;
}

ClassDescription<SMHiggsWidthGenerator> SMHiggsWidthGenerator::initSMHiggsWidthGenerator;
// Definition of the static class description member.

void SMHiggsWidthGenerator::Init() {

  static ClassDocumentation<SMHiggsWidthGenerator> documentation
    ("The SMHiggsWidthGenerator class calculates the running Higgs width as in "
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
     &SMHiggsWidthGenerator::widthopt_, 2, false, false);
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
     &SMHiggsWidthGenerator::offshell_, 10., 0.01, 100.0,
     false, false, Interface::limited);

}

Energy SMHiggsWidthGenerator::width(const ParticleData & in, Energy m) const {
  if(widthopt_==1) {
    return in.width();
  }
  else if(widthopt_ <=3 ) {
    Energy higgswidth = ZERO;
    for (unsigned int i = 0; i < decayModes().size(); ++i)
      higgswidth += partialWidth(i,m);
    return higgswidth;
  }
  else  
    throw Exception() << "Unknown width option in SMHiggsWidthGenerator::width()" 
		      << Exception::runerror;
}

DecayMap SMHiggsWidthGenerator::rate(const ParticleData & p) const {
  if(mw_==ZERO) return p.decaySelector();
  else          return branching(p.mass(),p);
}

DecayMap SMHiggsWidthGenerator::rate(const Particle & p) {
  return branching(p.mass(),p.data());
}

DecayMap SMHiggsWidthGenerator::branching(Energy scale, 
					  const ParticleData & p) const {
  // if not using running width return original
  if(widthopt_==1) return p.decaySelector();
  // calculate the partial widths
  vector<Energy> partial;
  Energy total(ZERO);
  for(unsigned int ix=0;ix<decayModes().size();++ix) {
    partial.push_back(partialWidth(ix,scale));
    total+=partial.back();
  }
  // produce the new decay selector
  DecayMap newdm;
  for(unsigned int ix=0;ix<decayModes().size();++ix) {
    tDMPtr mode = decayModes()[ix];
    if(mode->orderedProducts().size()!=2||!mode->on()) continue;
    double br = partial[ix]/total;
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

Energy SMHiggsWidthGenerator::partialWidth(int iloc,Energy Mh) const {
  useMe();
  using Constants::pi;
  if(Mh!=qlast_) {
    qlast_ = Mh;
    Energy2 q2 = sqr(qlast_);
    // couplings
    lambdaQCD_ = generator()->standardModel()->LambdaQCD(q2);
    alphaEM_   = generator()->standardModel()->alphaEM();
    alphaS_    = generator()->standardModel()->alphaS(q2);
    // QCD correction factors for H -> f fbar
    double nflavour=0.;
    for (unsigned int i =1; i <= 6; ++i) 
      if (2.0*qmass_[i] < Mh) nflavour+=1.;
    // All calculation are being done for Monte-Carlo QCD Lambda, except Higgs width...
    // not needed in C++ as should be normal lambda
    // MC only in shower
    //double bcoeff4=(11.*ca_-10.)/(12.*pi);
    //double kfac=ca_*(67./18.-sqr(pi)/6.)-25./9.;
    //lambdaQCD_ /= exp(kfac/(4.*pi*bcoeff4))/sqrt(2.);
    double k1 = 5./sqr(pi);
    double k0 = 3./(4.*sqr(pi));
    beta0_ = (11.*ca_-2.0*nflavour)/3.;
    double beta1 = (34.*sqr(ca_)-(10.*ca_+6.*cf_)*nflavour)/3.;
    gam0_ = -8.;
    double gam1 = -404./3.+40.*nflavour/9.;
    double SClog = log(sqr(Mh/lambdaQCD_));
    if(SClog<=0.)
      cd_ = 1.;
    else
      cd_ = 1.+(k1/k0-2.*gam0_+gam0_*beta1/sqr(beta0_)*log(SClog)+
 		(gam0_*beta1-gam1*beta0_)/sqr(beta0_))/(beta0_*SClog);
    gfermiinv_ = 8.*sw2_*sqr(mw_)/alphaEM_;
  }
  // output value
  Energy3 output(ZERO);
  // work out which mode
  map<int,int>::const_iterator it = locMap_.find(iloc);
  if(it==locMap_.end()) return ZERO;
  int imode = it->second;
  // quark modes
  if(imode<=6) {
    Energy mf = qmass_[imode];
    double xf = sqr(mf/Mh);
    if( xf >= 0.25 ) return ZERO;
    if(widthopt_==2) {
      if (mf > lambdaQCD_) mf *= pow(log(Mh/lambdaQCD_)/log(mf/lambdaQCD_),
 				     gam0_/(2.0*beta0_));
      output = ca_*Mh*sqr(mf)*pow(1.0-4.0*xf,1.5)*cd_;
    }
    else {
      output = ca_*Mh*sqr(mf)*pow(1.0-4.0*xf,1.5);
    }
  }
  // lepton modes
  else if(imode<=9) {
    Energy mf = lmass_[imode-7];
    double xf = sqr(mf/Mh);
    if (xf < 0.25) output = Mh*sqr(mf)*pow(1.0-4.0*xf,1.5);
  }
  // H->W*W*
  else if(imode==10) {
    if(widthopt_==2) {
      double xw = sqr(mw_/Mh);
      double yw = mw_*gamw_/sqr(Mh);
      output = pow<3,1>(Mh)*HwDoubleBW(xw,yw)/2.;
    }
    else {
      double xfw = sqr(mw_/Mh);
      if(2.0*mw_ < Mh) output = pow<3,1>(Mh)*sqrt(1-4.0*xfw)*(1-xfw+0.75*sqr(xfw))/2.;
    }
  }
  // H->Z*Z*
  else if(imode==11) {
    if(widthopt_==2) {
      double xz = sqr(mz_/Mh);
      double yz = mz_*gamz_/sqr(Mh);
      output = pow<3,1>(Mh)*HwDoubleBW(xz,yz)/4.;
    }
    else {
      double xfz = sqr(mz_/Mh);
      if (2.0*mz_ < Mh) output = pow<3,1>(Mh)*sqrt(1-4.0*xfz)*(1-xfz+0.75*sqr(xfz))/4.;
    }
  }
  // H->gamma gamma
  else if(imode==12) {
    double taut = sqr(2.0*qmass_[ParticleID::t]/Mh);
    double tauw = sqr(2.0*mw_/Mh);
    Complex ftaut = HwW2 (taut);
    Complex ftauw = HwW2 (tauw);
    double re = 4.0/3.0*(-2.0*taut*(1.0+(1.0-taut)*ftaut.real())) + 
      (2.0+3.0*tauw*(1+(2.0-tauw)*ftauw.real()));
    double im = 4.0/3.0*(-2.0*taut*(    (1.0-taut)*ftaut.imag())) +
      (   3.0*tauw*(  (2.0-tauw)*ftauw.imag()));
    output = sqr(alphaEM_/pi)*pow<3,1>(Mh)*(sqr(re)+ sqr(im))/32.;
  }
  // H -> gg
  else if(imode==13) {
    double tau = sqr(2.0*qmass_[ParticleID::t]/Mh);
    Complex ftau = HwW2(tau);
    double re = 1+(1.0-tau)*ftau.real();
    double im =   (1.0-tau)*ftau.imag();
    output = sqr(alphaS_/pi)*pow<3,1>(Mh)*sqr(tau)*(sqr(re)+ sqr(im))/4.;
  }
  return output/gfermiinv_;
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
        itgerl += 2.0*(sqr(1-x1-x2)+8.0*x1*x2)*sqrt(sq)
	  /(sqr(x1-x)+y*y)*y/(sqr(x2-x)+y*y)*y*fac1*fac2;
        ++j;
      }
    }
  } 
  else {
    // Integration using tan theta substitution
    double th1low = atan2 ((0.0-x),y);
    double th1high = atan2 ((1.0-x),y);
    double fac1 = (th1high-th1low)/nbin;
    for (unsigned int i = 0; i < nbin; ++i) {
      double th1 = (i+0.5)*fac1+th1low;
      double x1 = y*tan(th1)+x;
      double x2max = min(x1,sqr(1-sqrt(x1)));
      double th2low = atan2 ((0-x),y);
      double th2high = atan2 ((x2max-x),y);
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
  // extract W and Z mass and width
  tPDPtr w = getParticleData(ParticleID::Wplus);
  tPDPtr z = getParticleData(ParticleID::Z0);
  mw_   = w->mass();
  mz_   = z->mass();
  gamw_ = w->width();
  gamz_ = z->width();
  // quark masses
  for(int ix=1;ix<7;++ix) {
    tcPDPtr q = getParticleData(ix);
    qmass_[ix] = q->mass();
  }
  // lepton masses
  for ( int ix=0; ix<3; ++ix ) {
    tcPDPtr lepton = getParticleData(11+2*ix);
    lmass_[ix] = lepton->mass();
  }
  // sin2 theta_w
  sw2_ = generator()->standardModel()->sin2ThetaW();
  // casmirs
  double ncolour = generator()->standardModel()->Nc();
  ca_ = ncolour;
  cf_ = (sqr(ncolour)-1.0)/(2.0*ca_);
  GenericWidthGenerator::doinit();
  if(particle()->widthGenerator()!=this) return;
  // construct the map
  for(unsigned int ix=0;ix<decayModes().size();++ix) {
    tDMPtr mode = decayModes()[ix];
    if(mode->orderedProducts().size()!=2) continue;
    // particle antiparticle
    long id=abs(mode->orderedProducts()[0]->id());
    if(mode->orderedProducts()[0]->id()==-mode->orderedProducts()[1]->id()) {
      // leptons
      if(id>=11&&id<=15&&(id-9)%2==0) 
	id = (id+3)/2;
      // WW
      else if(id==ParticleID::Wplus)       
	id = 10;
      // unknown mode
      else if(id>6) 
	continue;
    }
    else if(mode->orderedProducts()[0]->id()==mode->orderedProducts()[1]->id()) {
      // gamma gamma
      if(id==ParticleID::Z0)         id = 11;
      else if(id==ParticleID::gamma) id = 12;
      else if(id==ParticleID::g)     id = 13;
      // unknown mode
      else continue;      
    }
    // unknown mode
    else continue;
    // set pu the map
    locMap_[ix] = id;
  }
  // reset the width and set the limits
  if(particle()->id() != ParticleID::h0)
    throw Exception() << "Must be the Standard Model Higgs boson "
		      << "in SMHiggsWidthGenerator::doinit()"
		      << Exception::runerror;
  Energy wid = width(*particle(),particle()->mass());
  particle()->width   (wid);
  particle()->widthCut(offshell_*wid);
}

pair<Energy,Energy> SMHiggsWidthGenerator::width(Energy scale,
						 const ParticleData & p) const {
  if(widthopt_==1) return make_pair(p.width(),p.width());
  // calculate the partial widths
  vector<Energy> partial;
  Energy total(ZERO);
  for(unsigned int ix=0;ix<decayModes().size();++ix) {
    partial.push_back(partialWidth(ix,scale));
    total += partial.back();
  }
  // sum for the on modes
  Energy partialon(ZERO);
  for(unsigned int ix=0;ix<decayModes().size();++ix) {
    tDMPtr mode= decayModes()[ix];
    if(!mode->on()||mode->orderedProducts().size()!=2) continue;
    partialon += partial[ix];
  }
  return make_pair(partialon,total);
}
