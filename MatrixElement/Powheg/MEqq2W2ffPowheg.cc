// -*- C++ -*-
//
// MEqq2W2ffPowheg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2W2ffPowheg class.
//

#include "MEqq2W2ffPowheg.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Utilities/Maths.h"

using namespace Herwig;
using Herwig::Math::ReLi2;

MEqq2W2ffPowheg::MEqq2W2ffPowheg() : 
  _gluon(), TR_(0.5), CF_(4./3.),
  _contrib(1)    ,_nlo_alphaS_opt(0), _fixed_alphaS(0.115895),
  _a(0.5)        ,_p(0.7)           , _eps(1.0e-8), _scaleopt(1),
  _fixedScale(100.*GeV), _scaleFact(1.) {
  massOption(vector<unsigned int>(2,1));
}

void MEqq2W2ffPowheg::doinit() {
  // gluon ParticleData object
  _gluon = getParticleData(ParticleID::g);
  MEqq2W2ff::doinit();
}

Energy2 MEqq2W2ffPowheg::scale() const {
  return _scaleopt == 0 ? sqr(_fixedScale) : _scaleFact*sHat();
}

void MEqq2W2ffPowheg::persistentOutput(PersistentOStream & os) const {
  os << _contrib << _nlo_alphaS_opt << _fixed_alphaS << _a << _p << _gluon
     << _scaleopt << ounit(_fixedScale,GeV) << _scaleFact;
}

void MEqq2W2ffPowheg::persistentInput(PersistentIStream & is, int) { 
  is >> _contrib >> _nlo_alphaS_opt >> _fixed_alphaS >> _a >> _p >> _gluon
     >> _scaleopt >> iunit(_fixedScale,GeV) >> _scaleFact;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEqq2W2ffPowheg,MEqq2W2ff>
describeHerwigMEqq2W2ffPowheg("Herwig::MEqq2W2ffPowheg", "HwMEHadron.so HwPowhegMEHadron.so");

void MEqq2W2ffPowheg::Init() {

    static ClassDocumentation<MEqq2W2ffPowheg> documentation
    ("The MEqq2W2ffPowheg class implements the matrix element for"
     "q qbar to Standard Model fermions via W exchange using helicity amplitude"
     "techniques including the NLO correction in the POWHEG formalism",
     "The qq$\\to$W$\\to$ff POWHEG matrix element is described in \\cite{Hamilton:2008pd}.",
     "%\\cite{Hamilton:2008pd}\n"
     "\\bibitem{Hamilton:2008pd}\n"
     "  K.~Hamilton, P.~Richardson and J.~Tully,\n"
     "  ``A Positive-Weight Next-to-Leading Order Monte Carlo Simulation of Drell-Yan\n"
     "  Vector Boson Production,''\n"
     "  JHEP {\\bf 0810} (2008) 015\n"
     "  [arXiv:0806.0290 [hep-ph]].\n"
     "  %%CITATION = JHEPA,0810,015;%%\n");
   static Switch<MEqq2W2ffPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEqq2W2ffPowheg::_contrib, 1, false, false);
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

  static Switch<MEqq2W2ffPowheg,unsigned int> interfaceNLOalphaSopt
    ("NLOalphaSopt",
     "Whether to use a fixed or a running QCD coupling for the NLO weight",
     &MEqq2W2ffPowheg::_nlo_alphaS_opt, 0, false, false);
  static SwitchOption interfaceNLOalphaSoptRunningAlphaS
    (interfaceNLOalphaSopt,
     "RunningAlphaS",
     "Use the usual running QCD coupling evaluated at scale scale()",
     0);
  static SwitchOption interfaceNLOalphaSoptFixedAlphaS
    (interfaceNLOalphaSopt,
     "FixedAlphaS",
     "Use a constant QCD coupling for comparison/debugging purposes",
     1);

  static Parameter<MEqq2W2ffPowheg,double> interfaceFixedNLOalphaS
    ("FixedNLOalphaS",
     "The value of alphaS to use for the nlo weight if _nlo_alphaS_opt=1",
     &MEqq2W2ffPowheg::_fixed_alphaS, 0.115895, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<MEqq2W2ffPowheg,double> interfaceCorrectionCoefficient
    ("CorrectionCoefficient",
     "The magnitude of the correction term to reduce the negative contribution",
     &MEqq2W2ffPowheg::_a, 0.5, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<MEqq2W2ffPowheg,double> interfaceCorrectionPower
    ("CorrectionPower",
     "The power of the correction term to reduce the negative contribution",
     &MEqq2W2ffPowheg::_p, 0.7, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<MEqq2W2ffPowheg,unsigned int> interfaceScaleOption
    ("ScaleOption",
     "Option for the scale to be used",
     &MEqq2W2ffPowheg::_scaleopt, 1, false, false);
  static SwitchOption interfaceScaleOptionFixed
    (interfaceScaleOption,
     "Fixed",
     "Use a fixed scale",
     0);
  static SwitchOption interfaceScaleOptionsHat
    (interfaceScaleOption,
     "Dynamic",
     "Use the off-shell vector boson mass as the scale",
     1);

  static Parameter<MEqq2W2ffPowheg,Energy> interfaceFixedScale
    ("FixedScale",
     "The fixed scale to use if required",
     &MEqq2W2ffPowheg::_fixedScale, GeV, 100.0*GeV, 10.0*GeV, 1000.0*GeV,
     false, false, Interface::limited);

  static Parameter<MEqq2W2ffPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before sHat if using a running scale",
     &MEqq2W2ffPowheg::_scaleFact, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
}

int MEqq2W2ffPowheg::nDim() const {
  return 3;
}

bool MEqq2W2ffPowheg::generateKinematics(const double * r) {
  _xt=*(r+1);
  _v =*(r+2);
  return MEqq2W2ff::generateKinematics(r);
}

CrossSection MEqq2W2ffPowheg::dSigHatDR() const {
  // Get Born momentum fractions xbar_a and xbar_b:
  _xb_a = lastX1();
  _xb_b = lastX2();
  return MEqq2W2ff::dSigHatDR()*NLOweight();
}

double MEqq2W2ffPowheg::NLOweight() const {
  // If only leading order is required return 1:
  if(_contrib==0) return 1.;
  useMe();
  // Get particle data for QCD particles:
  _parton_a=mePartonData()[0];
  _parton_b=mePartonData()[1];
  // get BeamParticleData objects for PDF's
  _hadron_A=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().first->dataPtr());
  _hadron_B=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().second->dataPtr());
  // If necessary swap the particle data vectors so that _xb_a, 
  // mePartonData[0], beam[0] relate to the inbound quark: 
  if(!(lastPartons().first ->dataPtr()==_parton_a&&
       lastPartons().second->dataPtr()==_parton_b)) {
    swap(_xb_a    ,_xb_b);
    swap(_hadron_A,_hadron_B);
  }
  // calculate the PDF's for the Born process
  _oldq    = _hadron_A->pdf()->xfx(_hadron_A,_parton_a,scale(),_xb_a)/_xb_a;
  _oldqbar = _hadron_B->pdf()->xfx(_hadron_B,_parton_b,scale(),_xb_b)/_xb_b;
  // Calculate alpha_S
  _alphaS2Pi = _nlo_alphaS_opt==1 ? _fixed_alphaS : SM().alphaS(scale());
  _alphaS2Pi /= 2.*Constants::pi;
  // Calculate the invariant mass of the dilepton pair
  _mll2 = sHat();
  _mu2  = scale();
  // Calculate the integrand
  // q qbar contribution
  double wqqvirt      = Vtilde_qq();
  double wqqcollin    = Ctilde_qq(x(_xt,1.),1.) + Ctilde_qq(x(_xt,0.),0.);
  double wqqreal      = Ftilde_qq(_xt,_v);
  double wqq          = wqqvirt+wqqcollin+wqqreal;
  // q g contribution
  double wqgcollin    = Ctilde_qg(x(_xt,0.),0.);
  double wqgreal      = Ftilde_qg(_xt,_v);
  double wqg          = wqgreal+wqgcollin;
  // g qbar contribution
  double wgqbarcollin = Ctilde_gq(x(_xt,1.),1.);
  double wgqbarreal   = Ftilde_gq(_xt,_v);
  double wgqbar       = wgqbarreal+wgqbarcollin;
  // total
  double wgt          = 1.+(wqq+wqg+wgqbar);
  //trick to try and reduce neg wgt contribution
  if(_xt<1.-_eps)
    wgt += _a*(1./pow(1.-_xt,_p)-(1.-pow(_eps,1.-_p))/(1.-_p)/(1.-_eps));
  // return the answer
  assert(isfinite(wgt));
  return _contrib==1 ? max(0.,wgt) : max(0.,-wgt);
}

double MEqq2W2ffPowheg::x(double xt, double v) const {
  double x0(xbar(v));
  return x0+(1.-x0)*xt;
}

double MEqq2W2ffPowheg::x_a(double x, double v) const {
  if(x==1.) return _xb_a;
  if(v==0.) return _xb_a;
  if(v==1.) return _xb_a/x;
  return (_xb_a/sqrt(x))*sqrt((1.-(1.-x)*(1.-v))/(1.-(1.-x)*v));
}

double MEqq2W2ffPowheg::x_b(double x, double v) const {
  if(x==1.) return _xb_b;
  if(v==0.) return _xb_b/x;
  if(v==1.) return _xb_b;
  return (_xb_b/sqrt(x))*sqrt((1.-(1.-x)*v)/(1.-(1.-x)*(1.-v)));
}

double MEqq2W2ffPowheg::xbar(double v) const {
  double xba2(sqr(_xb_a)), xbb2(sqr(_xb_b)), omv(-999.);
  double xbar1(-999.), xbar2(-999.);
  if(v==1.) return _xb_a;
  if(v==0.) return _xb_b;
  omv = 1.-v;
  xbar1=4.*  v*xba2/
    (sqrt(sqr(1.+xba2)*4.*sqr(omv)+16.*(1.-2.*omv)*xba2)+2.*omv*(1.-_xb_a)*(1.+_xb_a));
  xbar2=4.*omv*xbb2/
    (sqrt(sqr(1.+xbb2)*4.*sqr(  v)+16.*(1.-2.*  v)*xbb2)+2.*  v*(1.-_xb_b)*(1.+_xb_b));
  return max(xbar1,xbar2);
}

double MEqq2W2ffPowheg::Ltilde_qq(double x, double v) const {
  if(x==1.) return 1.;
  double xa(x_a(x,v)),xb(x_b(x,v));
  double newq    = (_hadron_A->pdf()->xfx(_hadron_A,_parton_a,scale(),   xa)/   xa);
  double newqbar = (_hadron_B->pdf()->xfx(_hadron_B,_parton_b,scale(),   xb)/   xb);
  return( newq * newqbar / _oldq / _oldqbar );
}

double MEqq2W2ffPowheg::Ltilde_qg(double x, double v) const {
  double xa(x_a(x,v)),xb(x_b(x,v));
  double newq    = (_hadron_A->pdf()->xfx(_hadron_A,_parton_a,scale(),   xa)/   xa);
  double newg2   = (_hadron_B->pdf()->xfx(_hadron_B,_gluon   ,scale(),   xb)/   xb);
  return( newq * newg2 / _oldq / _oldqbar );
}

double MEqq2W2ffPowheg::Ltilde_gq(double x, double v) const {
  double xa(x_a(x,v)),xb(x_b(x,v));
  double newg1   = (_hadron_A->pdf()->xfx(_hadron_A,_gluon   ,scale(),   xa)/   xa);
  double newqbar = (_hadron_B->pdf()->xfx(_hadron_B,_parton_b,scale(),   xb)/   xb);
  return( newg1 * newqbar / _oldq / _oldqbar );
}

double MEqq2W2ffPowheg::Vtilde_qq() const {
  return _alphaS2Pi*CF_*(-3.*log(_mu2/_mll2)+(2.*sqr(Constants::pi)/3.)-8.);
}

double MEqq2W2ffPowheg::Ccalbar_qg(double x) const {
  return (sqr(x)+sqr(1.-x))*(log(_mll2/(_mu2*x))+2.*log(1.-x))+2.*x*(1.-x);
}

double MEqq2W2ffPowheg::Ctilde_qg(double x, double v) const {
  return  _alphaS2Pi*TR_ * ((1.-xbar(v))/x) * Ccalbar_qg(x)*Ltilde_qg(x,v);
}

double MEqq2W2ffPowheg::Ctilde_gq(double x, double v) const {
  return  _alphaS2Pi*TR_ * ((1.-xbar(v))/x) * Ccalbar_qg(x)*Ltilde_gq(x,v);
}

double MEqq2W2ffPowheg::Ctilde_qq(double x, double v) const {
  double wgt 
    = ((1.-x)/x+(1.+x*x)/(1.-x)/x*(2.*log(1.-x)-log(x)))*Ltilde_qq(x,v)
    -  4.*log(1.-x)/(1.-x)
    +  2./(1.-xbar(v))*log(1.-xbar(v))*log(1.-xbar(v))
    + (2./(1.-xbar(v))*log(1.-xbar(v))-2./(1.-x)+(1.+x*x)/x/(1.-x)*Ltilde_qq(x,v))
    *log(_mll2/_mu2);
  return _alphaS2Pi*CF_*(1.-xbar(v))*wgt;    
}

double MEqq2W2ffPowheg::Fcal_qq(double x, double v) const {
  return (sqr(1.-x)*(1.-2.*v*(1.-v))+2.*x)/x*Ltilde_qq(x,v);
}

double MEqq2W2ffPowheg::Fcal_qg(double x, double v) const {
  return ((1.-xbar(v))/x)*
    (2.*x*(1.-x)*v+sqr((1.-x)*v)+sqr(x)+sqr(1.-x))*Ltilde_qg(x,v);
}

double MEqq2W2ffPowheg::Fcal_gq(double x, double v) const {
  return ((1.-xbar(v))/x)*
    (2.*x*(1.-x)*(1.-v)+sqr((1.-x)*(1.-v))+sqr(x)+sqr(1.-x))*Ltilde_gq(x,v);
}

double MEqq2W2ffPowheg::Ftilde_qg(double xt, double v) const {
  return _alphaS2Pi*TR_*
    ( Fcal_qg(x(xt,v),v) - Fcal_qg(x(xt,0.),0.) )/v;
}

double MEqq2W2ffPowheg::Ftilde_gq(double xt, double v) const {
  return _alphaS2Pi*TR_*
    ( Fcal_gq(x(xt,v),v) - Fcal_gq(x(xt,1.),1.) )/(1.-v);
}

double MEqq2W2ffPowheg::Ftilde_qq(double xt, double v) const {
  double eps(1e-10);
  // is emission into regular or singular region?
  if(xt>=0. && xt<1.-eps && v>eps && v<1.-eps) { 
    // x<1, v>0, v<1 (regular emission, neither soft or collinear):
    return _alphaS2Pi*CF_*
      (( ( Fcal_qq(x(xt, v), v) - Fcal_qq(x(xt,1.),1.) ) / (1.-v)+
	 ( Fcal_qq(x(xt, v), v) - Fcal_qq(x(xt,0.),0.) ) / v )/(1.-xt)
       + ( log(1.-xbar(v)) - log(1.-_xb_a))*2./(1.-v)
       + ( log(1.-xbar(v)) - log(1.-_xb_b))*2./v);
  } 
  else {
    // make sure emission is actually in the allowed phase space:
    if(!(v>=0. && v<=1. && xt>=0. && xt<=1.)) {
      ostringstream s;
      s << "MEqq2W2ffPowheg::Ftilde_qq : \n" << "xt(" << xt << ") and / or v(" 
	<< v << ") not in the phase space.";
      generator()->logWarning(Exception(s.str(),Exception::warning)); 
      return 0.;
    }
    // is emission soft singular?
    if(xt>=1.-eps) {
      // x=1:
      if(v<=eps) {
	// x==1, v=0 (soft and collinear with particle b):
	return _alphaS2Pi*CF_*
	  (   ( log(1.-xbar(v)) - log(1.-_xb_a))*2./(1.-v)
	      );
      } else if(v>=1.-eps) {
	// x==1, v=1 (soft and collinear with particle a):
	return _alphaS2Pi*CF_*
	  (   ( log(1.-xbar(v)) - log(1.-_xb_b))*2./v
	      );
      } else {
	// x==1, 0<v<1 (soft wide angle emission):
	return _alphaS2Pi*CF_*
	  (   ( log(1.-xbar(v)) - log(1.-_xb_a))*2./(1.-v)
	      + ( log(1.-xbar(v)) - log(1.-_xb_b))*2./v
	      );
      }
    } else {
      // x<1:
      if(v<=eps) {
	// x<1 but v=0 (collinear with particle b, but not soft):
	return _alphaS2Pi*CF_*
	  ( ( ( Fcal_qq(x(xt, v), v) - Fcal_qq(x(xt,1.),1.) ) / (1.-v)
	      )/(1.-xt)
	    + ( log(1.-xbar(v)) - log(1.-_xb_a))*2./(1.-v) 
	    );
      } else if(v>=1.-eps) {
	// x<1 but v=1 (collinear with particle a, but not soft):
	return _alphaS2Pi*CF_*
	  ( ( ( Fcal_qq(x(xt, v), v) - Fcal_qq(x(xt,0.),0.) ) / v 
	      )/(1.-xt)
	    + ( log(1.-xbar(v)) - log(1.-_xb_b))*2./v 
	    );
      }
    }
  }
  return 0.;
}
