// -*- C++ -*-
//
// MEPP2HiggsPowheg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2HiggsPowheg class.
//

#include "MEPP2HiggsPowheg.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
 
using namespace Herwig;

MEPP2HiggsPowheg::MEPP2HiggsPowheg() : 
  CF_(4./3.)  ,  CA_(3.)            , TR_(0.5)        , nlf_(5)          ,
  beta0_((11.*CA_/3. - 4.*TR_*nlf_/3.)/(4.*Constants::pi))                 ,
  contrib_(1) ,  nlo_alphaS_opt_(0) , fixed_alphaS_(0.118109485),
  scaleopt_(1),  mu_F_(100.*GeV)    ,  mu_UV_(100.*GeV) , scaleFact_(1.) 
{}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2HiggsPowheg,MEPP2Higgs>
describeHerwigMEPP2HiggsPowheg("Herwig::MEPP2HiggsPowheg", "HwMEHadron.so HwPowhegMEHadron.so");

void MEPP2HiggsPowheg::persistentOutput(PersistentOStream & os) const {
  os << contrib_       << nlo_alphaS_opt_  << fixed_alphaS_ 
     << scaleopt_      << ounit(mu_F_,GeV) << ounit(mu_UV_,GeV)   
     << scaleFact_     ;
}

void MEPP2HiggsPowheg::persistentInput(PersistentIStream & is, int) {
  is >> contrib_       >> nlo_alphaS_opt_  >> fixed_alphaS_ 
     >> scaleopt_      >> iunit(mu_F_,GeV) >> iunit(mu_UV_,GeV)  
     >> scaleFact_     ;
}

void MEPP2HiggsPowheg::Init() {

  static ClassDocumentation<MEPP2HiggsPowheg> documentation
    ("The MEPP2HiggsPowheg class implements the matrix elements for"
     " Higgs production (with decay H->W-W+) in hadron-hadron collisions.",
     "The PP$\\to$Higgs POWHEG matrix element is described in \\cite{Hamilton:2009za}.",
     "%\\cite{Hamilton:2009za}\n"
     "\\bibitem{Hamilton:2009za}\n"
     "  K.~Hamilton, P.~Richardson and J.~Tully,\n"
     "  %``A Positive-Weight Next-to-Leading Order Monte Carlo Simulation for Higgs\n"
     "  %Boson Production,''\n"
     "  JHEP {\\bf 0904} (2009) 116\n"
     "  [arXiv:0903.4345 [hep-ph]].\n"
     "  %%CITATION = JHEPA,0904,116;%%\n"
     );

  static Switch<MEPP2HiggsPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEPP2HiggsPowheg::contrib_, 1, false, false);
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

  static Switch<MEPP2HiggsPowheg,unsigned int> interfaceNLOalphaSopt
    ("NLOalphaSopt",
     "Whether to use a fixed or a running QCD coupling for the NLO weight",
     &MEPP2HiggsPowheg::nlo_alphaS_opt_, 0, false, false);
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

  static Parameter<MEPP2HiggsPowheg,double> interfaceFixedNLOalphaS
    ("FixedNLOalphaS",
     "The value of alphaS to use for the nlo weight if nlo_alphaS_opt_=1",
     &MEPP2HiggsPowheg::fixed_alphaS_, 0.11803463, 0., 1.0,
     false, false, Interface::limited);

  static Switch<MEPP2HiggsPowheg,unsigned int> interfaceFactorizationScaleOption
    ("FactorizationScaleOption",
     "Option for the choice of factorization (and renormalization) scale",
     &MEPP2HiggsPowheg::scaleopt_, 1, false, false);
  static SwitchOption interfaceDynamic
    (interfaceFactorizationScaleOption,
     "Dynamic",
     "Dynamic factorization scale equal to the current sqrt(sHat())",
     1);
  static SwitchOption interfaceFixed
    (interfaceFactorizationScaleOption,
     "Fixed",
     "Use a fixed factorization scale set with FactorizationScaleValue",
     2);

  static Parameter<MEPP2HiggsPowheg,Energy> interfaceFactorizationScaleValue
    ("FactorizationScaleValue",
     "Value to use in the event of a fixed factorization scale",
     &MEPP2HiggsPowheg::mu_F_, GeV, 100.0*GeV, 50.0*GeV, 500.0*GeV,
     true, false, Interface::limited);

  static Parameter<MEPP2HiggsPowheg,Energy> interfaceRenormalizationScaleValue
    ("RenormalizationScaleValue",
     "Value to use for the (UV) renormalization scale",
     &MEPP2HiggsPowheg::mu_UV_, GeV, 100.0*GeV, 50.0*GeV, 500.0*GeV,
     true, false, Interface::limited);

  static Parameter<MEPP2HiggsPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before sHat if using a running scale",
     &MEPP2HiggsPowheg::scaleFact_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

}

Energy2 MEPP2HiggsPowheg::scale() const {
  return scaleopt_ == 1 ?  sqr(scaleFact_)*sHat() : sqr(scaleFact_*mu_F_);
}

int MEPP2HiggsPowheg::nDim() const {
  return 2;
}

bool MEPP2HiggsPowheg::generateKinematics(const double * r) {
  // Generate the radiative integration variables:
  xt_= *r;
  y_ = *(r+1) * 2. - 1.;
  // Continue with lo matrix element code:
  return MEPP2Higgs::generateKinematics(r);
}

double MEPP2HiggsPowheg::me2() const {
  useMe();
  double output = MEPP2Higgs::me2();
  if (mePartonData()[0]->id() == ParticleID::g && 
      mePartonData()[1]->id() == ParticleID::g) {
    get_born_variables();
    // NB - lo_ggME_ equals sqr(alphaS/(pi*vev))*
    // sqr(p2_)/576. _ALL_IN_MeVs_!
    lo_ggME_ = output;
    if(output==0.) return 0.;
    output *= NLOweight();
  }
  return output;
}

double MEPP2HiggsPowheg::NLOweight() const {
  // If only leading order is required return 1:
  if(contrib_==0) return 1.;

  // ATTENTION!!! - for consistency with LO matrix element ALL
  // energy dimensions should be understood as MeVs!

  // If necessary swap the particle data vectors so that xbp_, 
  // mePartonData[0], beam[0] relate to the inbound particle a. 
  // [ Irrelevant for gg collisions!! Kept as an aide memoire. ]
  if(!(lastPartons().first ->dataPtr()==a_lo_&&
       lastPartons().second->dataPtr()==b_lo_)) {
    swap(xbp_     ,xbm_     );
    swap(etabarp_ ,etabarm_ );
    swap(hadron_A_,hadron_B_);
  }

  // calculate the PDF's for the Born process
  lo_lumi_ = hadron_A_->pdf()->xfx(hadron_A_,a_lo_,scale(),xbp_)/xbp_
           * hadron_B_->pdf()->xfx(hadron_B_,b_lo_,scale(),xbm_)/xbm_;

  // Calculate alpha_S and alpha_S/(2*pi)
  alphaS_ = nlo_alphaS_opt_==1 ? fixed_alphaS_ : SM().alphaS(scale());
  double alsOn2pi(alphaS_/2./Constants::pi);

  // Particle data objects for the new plus and minus colliding partons.
  tcPDPtr a_nlo, b_nlo;

  // Calculate the integrand
  double wgt(0.);

  // g g contribution
  a_nlo=getParticleData(ParticleID::g);
  b_nlo=getParticleData(ParticleID::g);
  double wggvirt      = Vtilde_universal() + M_V_regular()/lo_ggME_;
  double wggcollin    = alsOn2pi
                      * ( Ctilde_Ltilde_gg_on_x(a_nlo,b_nlo,xt_, 1.) 
                        + Ctilde_Ltilde_gg_on_x(a_nlo,b_nlo,xt_,-1.));
  double wggreal      = alsOn2pi
                      * Rtilde_Ltilde_gg_on_x(a_nlo,b_nlo,xt_,y_);
  double wgg          = wggvirt + wggcollin + wggreal;
  // g q + g qbar contributions
  a_nlo=getParticleData(ParticleID::g);
  double wgqcollin(0.)   , wgqreal(0.)   , wgq(0.)   ;
  for(int ix=1; ix<=nlf_; ++ix) {
    b_nlo=getParticleData( ix);
    wgqcollin         = alsOn2pi*Ctilde_Ltilde_qg_on_x(a_nlo,b_nlo,xt_,-1.);
    wgqreal           = alsOn2pi*Rtilde_Ltilde_gq_on_x(a_nlo,b_nlo,xt_,y_);
    wgq              += wgqreal+wgqcollin;

    b_nlo=getParticleData(-ix);
    wgqcollin         = alsOn2pi*Ctilde_Ltilde_qg_on_x(a_nlo,b_nlo,xt_,-1.);
    wgqreal           = alsOn2pi*Rtilde_Ltilde_gq_on_x(a_nlo,b_nlo,xt_,y_);
    wgq              += wgqreal+wgqcollin;
  }
  // q g + qbar g contributions
  b_nlo=getParticleData(ParticleID::g);
  double wqgcollin(0.)   , wqgreal(0.)   , wqg(0.)   ;
  for(int ix=1; ix<=nlf_; ++ix) {
    a_nlo=getParticleData( ix);
    wqgcollin         = alsOn2pi*Ctilde_Ltilde_qg_on_x(a_nlo,b_nlo,xt_, 1.);
    wqgreal           = alsOn2pi*Rtilde_Ltilde_qg_on_x(a_nlo,b_nlo,xt_,y_);
    wqg              += wqgreal+wqgcollin;

    a_nlo=getParticleData(-ix);
    wqgcollin         = alsOn2pi*Ctilde_Ltilde_qg_on_x(a_nlo,b_nlo,xt_, 1.);
    wqgreal           = alsOn2pi*Rtilde_Ltilde_qg_on_x(a_nlo,b_nlo,xt_,y_);
    wqg              += wqgreal+wqgcollin;
  }
  // q qbar + qbar q contributions
  double wqqbarreal(0.), wqqbar(0.);
  for(int ix=1; ix<=nlf_; ++ix) {
    a_nlo=getParticleData( ix);
    b_nlo=getParticleData(-ix);
    wqqbarreal    = alsOn2pi*Rtilde_Ltilde_qqbar_on_x(a_nlo,b_nlo,xt_,y_);
    wqqbar       += wqqbarreal;
    a_nlo=getParticleData(-ix);
    b_nlo=getParticleData( ix);
    wqqbarreal    = alsOn2pi*Rtilde_Ltilde_qbarq_on_x(a_nlo,b_nlo,xt_,y_);
    wqqbar       += wqqbarreal;
  }
  // total
  wgt                 = 1.+(wgg+wgq+wqg+wqqbar);
  return contrib_==1 ? max(0.,wgt) : max(0.,-wgt);
}

void MEPP2HiggsPowheg::get_born_variables() const {
  assert(meMomenta().size()==3);
  // Particle data for QCD particles:
  a_lo_=mePartonData()[0];
  b_lo_=mePartonData()[1];
  // BeamParticleData objects for PDF's
  hadron_A_=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().first->dataPtr());
  hadron_B_=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().second->dataPtr());
  // Leading order momentum fractions and associated etabar's
  xbp_ = lastX1(); etabarp_ = sqrt(1.-xbp_);
  xbm_ = lastX2(); etabarm_ = sqrt(1.-xbm_);
  // Assign Born variables
  p2_    = sHat();
  s2_    = p2_;
  return;
}

double MEPP2HiggsPowheg::xbar(double y) const {
  double xbp2(sqr(xbp_)), xbm2(sqr(xbm_));
  double omy(-999.)     , opy( 999.)     ;
  double xbar1(-999.)   , xbar2(-999.)   ;
  if(y== 1.) return xbp_;
  if(y==-1.) return xbm_;
  omy = 1.-y; 
  opy = 1.+y;
  xbar1=2.*opy*xbp2/
    (sqrt(sqr(1.+xbp2)*sqr(omy)+16.*y*xbp2)+omy*(1.-xbp_)*(1.+xbp_));
  xbar2=2.*omy*xbm2/
    (sqrt(sqr(1.+xbm2)*sqr(opy)-16.*y*xbm2)+opy*(1.-xbm_)*(1.+xbm_));
  return max(xbar1,xbar2);
}

double MEPP2HiggsPowheg::etabar(double y) const {
  return sqrt(1.-xbar(y));
}

double MEPP2HiggsPowheg::xp(double x, double y) const {
  if(x== 1.) return xbp_  ;
  if(y==-1.) return xbp_  ;
  if(y== 1.) return xbp_/x;
  return (xbp_/sqrt(x))*sqrt((2.-(1.-x)*(1.-y))/(2.-(1.-x)*(1.+y)));
}

double MEPP2HiggsPowheg::xm(double x, double y) const {
  if(x== 1.) return xbm_  ;
  if(y==-1.) return xbm_/x;
  if(y== 1.) return xbm_  ;
  return (xbm_/sqrt(x))*sqrt((2.-(1.-x)*(1.+y))/(2.-(1.-x)*(1.-y)));
}

double MEPP2HiggsPowheg::Lhat_ab(tcPDPtr a, tcPDPtr b, double x, double y) const {
  double nlo_lumi(-999.);
  double xp_x_y(xp(x,y)),xm_x_y(xm(x,y));
  nlo_lumi = (hadron_A_->pdf()->xfx(hadron_A_,a,scale(),xp_x_y)/xp_x_y)
           * (hadron_B_->pdf()->xfx(hadron_B_,b,scale(),xm_x_y)/xm_x_y);
  return nlo_lumi / lo_lumi_;
}

double MEPP2HiggsPowheg::Vtilde_universal() const {
  return  alphaS_/2./Constants::pi*CA_ 
        * ( log(p2_/scale())*( 2.*(2.*Constants::pi*beta0_/CA_)
	     	          + 4.*log(etabarp_)+4.*log(etabarm_))
	                  + 8.*sqr(log(etabarp_)) + 8.*sqr(log(etabarm_))
	                  - 2.*sqr(Constants::pi)/3.
	  )
        + alphaS_/2./Constants::pi*CA_ 
        * ( 8./(1.+y_)*log(etabar(y_)/etabarm_)
     	  + 8./(1.-y_)*log(etabar(y_)/etabarp_)
          );
}

double MEPP2HiggsPowheg::Ctilde_Ltilde_gg_on_x(tcPDPtr a, tcPDPtr b, 
					       double xt, double y ) const {
  if(y!= 1.&&y!=-1.) { cout << "\nCtilde_gg::y value not allowed."; }
  if(y== 1.&&a->id()!=21) 
    cout << "\nCtilde_gg::for Cgg^plus  a must be a gluon! id = " 
	 << a->id() << "\n";
  if(y==-1.&&b->id()!=21) 
    cout << "\nCtilde_gg::for Cgg^minus b must be a gluon! id = " 
	 << b->id() << "\n";
  double x_pm      = x(xt,y);
  double etabar_pm = y == 1. ? etabarp_ : etabarm_ ;
  return ( ( (1./(1.-xt))*log(p2_/scale()/x_pm)+4.*log(etabar_pm)/(1.-xt)
       	   + 2.*log(1.-xt)/(1.-xt)
           )*2.*CA_*(x_pm+sqr(1.-x_pm)/x_pm+x_pm*sqr(1.-x_pm))

	 )*Lhat_ab(a,b,x_pm,y) / x_pm
       - ( ( (1./(1.-xt))*log(p2_/scale()     )+4.*log(etabar_pm)/(1.-xt)
	   + 2.*log(1.-xt)/(1.-xt)
	   )*2.*CA_
	 );
}

double MEPP2HiggsPowheg::Ctilde_Ltilde_qg_on_x(tcPDPtr a, tcPDPtr b, 
					       double xt, double y ) const {
  if(y!= 1.&&y!=-1.) { cout << "\nCtilde_qg::y value not allowed."; }
  if(y== 1.&&!(abs(a->id())>0&&abs(a->id())<7)) 
    cout << "\nCtilde_qg::for Cqg^plus  a must be a quark! id = " 
	 << a->id() << "\n";
  if(y==-1.&&!(abs(b->id())>0&&abs(b->id())<7)) 
    cout << "\nCtilde_qg::for Cqg^minus b must be a quark! id = "
	 << b->id() << "\n";
  double x_pm      = x(xt,y);
  double etabar_pm = y == 1. ? etabarp_ : etabarm_ ;
  return ( ( (1./(1.-xt))*log(p2_/scale()/x_pm)+4.*log(etabar_pm)/(1.-xt)
       	   + 2.*log(1.-xt)/(1.-xt)
           )*(1.-x_pm)*CF_*(1.+sqr(1.-x_pm))/x_pm
	 + sqr(etabar_pm)*CF_*x_pm
	 )*Lhat_ab(a,b,x_pm,y) / x_pm;
}

double MEPP2HiggsPowheg::Ctilde_Ltilde_gq_on_x(tcPDPtr a, tcPDPtr b, 
					       double xt, double y ) const {
  if(y!= 1.&&y!=-1.) { cout << "\nCtilde_gq::y value not allowed."; }
  if(y== 1.&&a->id()!=21)
    cout << "\nCtilde_gq::for Cgq^plus  a must be a gluon! id = " 
	 << a->id() << "\n";
  if(y==-1.&&b->id()!=21) 
    cout << "\nCtilde_gq::for Cgq^minus b must be a gluon! id = " 
	 << b->id() << "\n";
  double x_pm      = x(xt,y);
  double etabar_pm = y == 1. ? etabarp_ : etabarm_ ;
  return ( ( (1./(1.-xt))*log(p2_/scale()/x_pm)+4.*log(etabar_pm)/(1.-xt)
       	   + 2.*log(1.-xt)/(1.-xt)
           )*(1.-x_pm)*TR_*(sqr(x_pm)+sqr(1.-x_pm))
	 + sqr(etabar_pm)*TR_*2.*x_pm*(1.-x_pm)
	 )*Lhat_ab(a,b,x_pm,y) / x_pm;
}

double MEPP2HiggsPowheg::M_V_regular() const {
  // If the we are using a dynamic scale we use the same scale
  // for the renormalization scale as for the factorization scale
  Energy2 mu_UV2 = scaleopt_ == 1 ? scale() : sqr(mu_UV_);
  return alphaS_/2./Constants::pi*CA_*
                        (  11./3.
			+  4.*sqr(Constants::pi)/3.
			- (4.*Constants::pi*beta0_/CA_)*log(p2_/mu_UV2)
			)*lo_ggME_;
}

Energy2 MEPP2HiggsPowheg::t_u_M_R_qqbar(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*32./9./sqr(p2_)/s(xt,y)*tk(xt,y)*uk(xt,y)
           *( sqr(tk(xt,y)) + sqr(uk(xt,y))
            )*lo_ggME_;
}

Energy2 MEPP2HiggsPowheg::t_u_M_R_qbarq(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*32./9./sqr(p2_)/s(xt,y)*uk(xt,y)*tk(xt,y)
           *( sqr(uk(xt,y)) + sqr(tk(xt,y))
            )*lo_ggME_;
}

Energy2 MEPP2HiggsPowheg::t_u_M_R_gg(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*3./sqr(p2_)/s(xt,y)
           *( sqr(sqr(p2_    )) + sqr(sqr(s( xt,y)))
	    + sqr(sqr(tk(xt,y))) + sqr(sqr(uk(xt,y)))
            )*lo_ggME_;
}

Energy2 MEPP2HiggsPowheg::t_u_M_R_qg(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*-4./3./sqr(p2_)*uk(xt,y)
           *( sqr(s(xt,y)) + sqr(uk(xt,y))
            )*lo_ggME_;
}

Energy2 MEPP2HiggsPowheg::t_u_M_R_gq(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*-4./3./sqr(p2_)*tk(xt,y)
           *( sqr(s(xt,y)) + sqr(tk(xt,y))
            )*lo_ggME_;
}

double MEPP2HiggsPowheg::Rtilde_Ltilde_qqbar_on_x(tcPDPtr a , tcPDPtr b,
					     double  xt, double y ) const {
  return ( ( 
	     1./s(xt ,y  )
	   * t_u_M_R_qqbar(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt , 1.)
	   * t_u_M_R_qqbar(xt , 1.)*Lhat_ab(a,b,x(xt , 1.), 1.)

	   - 1./s( 1.,y  ) * t_u_M_R_qqbar( 1.,y  )

	   + 1./s( 1., 1.) * t_u_M_R_qqbar( 1., 1.)

           )*2./(1.-y)/(1.-xt)
	 + ( 
	     1./s(xt ,y  )
	   * t_u_M_R_qqbar(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt ,-1.)
	   * t_u_M_R_qqbar(xt ,-1.)*Lhat_ab(a,b,x(xt ,-1.),-1.)

	   - 1./s( 1.,y  ) * t_u_M_R_qqbar( 1.,y  )

	   + 1./s( 1.,-1.) * t_u_M_R_qqbar( 1.,-1.)

           )*2./(1.+y)/(1.-xt)
	 ) / lo_ggME_ / 8. / Constants::pi / alphaS_;
}

double MEPP2HiggsPowheg::Rtilde_Ltilde_qbarq_on_x(tcPDPtr a , tcPDPtr b,
					     double  xt, double y ) const {
  return ( ( 
	     1./s(xt ,y  )
	   * t_u_M_R_qbarq(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt , 1.)
	   * t_u_M_R_qbarq(xt , 1.)*Lhat_ab(a,b,x(xt , 1.), 1.)

	   - 1./s( 1.,y  ) * t_u_M_R_qbarq( 1.,y  )

	   + 1./s( 1., 1.) * t_u_M_R_qbarq( 1., 1.)

           )*2./(1.-y)/(1.-xt)
	 + ( 
	     1./s(xt ,y  )
	   * t_u_M_R_qbarq(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt ,-1.)
	   * t_u_M_R_qbarq(xt ,-1.)*Lhat_ab(a,b,x(xt ,-1.),-1.)

	   - 1./s( 1.,y  ) * t_u_M_R_qbarq( 1.,y  )

	   + 1./s( 1.,-1.) * t_u_M_R_qbarq( 1.,-1.)

           )*2./(1.+y)/(1.-xt)
	 ) / lo_ggME_ / 8. / Constants::pi / alphaS_;
}

double MEPP2HiggsPowheg::Rtilde_Ltilde_gg_on_x(tcPDPtr a , tcPDPtr b, 
					     double  xt, double y ) const {
  return ( ( 
	     1./s(xt ,y  )
	   * t_u_M_R_gg(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt , 1.)
	   * t_u_M_R_gg(xt , 1.)*Lhat_ab(a,b,x(xt , 1.), 1.)

	   - 1./s( 1.,y  ) * t_u_M_R_gg( 1.,y  )

	   + 1./s( 1., 1.) * t_u_M_R_gg( 1., 1.)

           )*2./(1.-y)/(1.-xt)
	 + ( 
	     1./s(xt ,y  )
	   * t_u_M_R_gg(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt ,-1.)
	   * t_u_M_R_gg(xt ,-1.)*Lhat_ab(a,b,x(xt ,-1.),-1.)

	   - 1./s( 1.,y  ) * t_u_M_R_gg( 1.,y  )

	   + 1./s( 1.,-1.) * t_u_M_R_gg( 1.,-1.)

           )*2./(1.+y)/(1.-xt)
	 ) / lo_ggME_ / 8. / Constants::pi / alphaS_;
}

double MEPP2HiggsPowheg::Rtilde_Ltilde_gq_on_x(tcPDPtr a , tcPDPtr b, 
					     double  xt, double y ) const {
  return ( ( 
	     1./s(xt ,y  )
	   * t_u_M_R_gq(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt , 1.)
	   * t_u_M_R_gq(xt , 1.)*Lhat_ab(a,b,x(xt , 1.), 1.)

	   - 1./s( 1.,y  ) * t_u_M_R_gq( 1.,y  )

	   + 1./s( 1., 1.) * t_u_M_R_gq( 1., 1.)

           )*2./(1.-y)/(1.-xt)
	 + ( 
	     1./s(xt ,y  )
	   * t_u_M_R_gq(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt ,-1.)
	   * t_u_M_R_gq(xt ,-1.)*Lhat_ab(a,b,x(xt ,-1.),-1.)

	   - 1./s( 1.,y  ) * t_u_M_R_gq( 1.,y  )

	   + 1./s( 1.,-1.) * t_u_M_R_gq( 1.,-1.)

           )*2./(1.+y)/(1.-xt)
	 ) / lo_ggME_ / 8. / Constants::pi / alphaS_;
}

double MEPP2HiggsPowheg::Rtilde_Ltilde_qg_on_x(tcPDPtr a , tcPDPtr b, 
					     double  xt, double y ) const {
  return ( ( 
	     1./s(xt ,y  )
	   * t_u_M_R_qg(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt , 1.)
	   * t_u_M_R_qg(xt , 1.)*Lhat_ab(a,b,x(xt , 1.), 1.)

	   - 1./s( 1.,y  ) * t_u_M_R_qg( 1.,y  )

	   + 1./s( 1., 1.) * t_u_M_R_qg( 1., 1.)

           )*2./(1.-y)/(1.-xt)
	 + ( 
	     1./s(xt ,y  )
	   * t_u_M_R_qg(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt ,-1.)
	   * t_u_M_R_qg(xt ,-1.)*Lhat_ab(a,b,x(xt ,-1.),-1.)

	   - 1./s( 1.,y  ) * t_u_M_R_qg( 1.,y  )

	   + 1./s( 1.,-1.) * t_u_M_R_qg( 1.,-1.)

           )*2./(1.+y)/(1.-xt)
         ) / lo_ggME_ / 8. / Constants::pi / alphaS_;
} 
