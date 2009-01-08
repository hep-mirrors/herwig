// -*- C++ -*-
//
// MEPP2VVPowheg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2VVPowheg class.
//

#include "MEPP2VVPowheg.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/MatrixElement/HardVertex.h"

using namespace Herwig;

MEPP2VVPowheg::MEPP2VVPowheg() :  
    CF_(4./3.),    TR_(0.5),
    contrib_(1),   nlo_alphaS_opt_(0) , fixed_alphaS_(0.118109485),
    removebr_(1)
{}

void MEPP2VVPowheg::persistentOutput(PersistentOStream & os) const {
  os << contrib_    << nlo_alphaS_opt_  << fixed_alphaS_
     << removebr_ ;
}

void MEPP2VVPowheg::persistentInput(PersistentIStream & is, int) {
  is >> contrib_    >> nlo_alphaS_opt_  >> fixed_alphaS_
     >> removebr_ ;
}

ClassDescription<MEPP2VVPowheg> MEPP2VVPowheg::initMEPP2VVPowheg;
// Definition of the static class description member.

void MEPP2VVPowheg::Init() {

  static Switch<MEPP2VVPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEPP2VVPowheg::contrib_, 1, false, false);
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

  static Switch<MEPP2VVPowheg,unsigned int> interfaceNLOalphaSopt
    ("NLOalphaSopt",
     "Whether to use a fixed or a running QCD coupling for the NLO weight",
     &MEPP2VVPowheg::nlo_alphaS_opt_, 0, false, false);
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

  static Parameter<MEPP2VVPowheg,double> interfaceFixedNLOalphaS
    ("FixedNLOalphaS",
     "The value of alphaS to use for the nlo weight if nlo_alphaS_opt_=1",
     &MEPP2VVPowheg::fixed_alphaS_, 0.11803463, 0., 1.0,
     false, false, Interface::limited);

  static Switch<MEPP2VVPowheg,unsigned int> interfaceremovebr
    ("removebr",
     "Whether to multiply the event weights by the MCFM branching ratios",
     &MEPP2VVPowheg::removebr_, 1, false, false);
  static SwitchOption interfaceProductionCrossSection
    (interfaceremovebr,
     "true",
     "Do not multiply in the branching ratios (default running)",
     1);
  static SwitchOption interfaceIncludeBRs
    (interfaceremovebr,
     "false",
     "Multiply by MCFM branching ratios for comparison/debugging purposes",
     0);

}

int MEPP2VVPowheg::nDim() const {
  int output = MEPP2VV::nDim(); 
  if(contrib_>0) output += 2;
  return output;
}

bool MEPP2VVPowheg::generateKinematics(const double * r) {
  if(contrib_>0) {
    // Generate the radiative integration variables:
    xt_= *r;
    y_ = *(r+1) * 2. - 1.;
  }
  // Continue with lo matrix element code:
  return MEPP2VV::generateKinematics(r);
}

double MEPP2VVPowheg::me2() const {
  double output(0.0);
  useMe();
  output = MEPP2VV::me2();
  double mcfm_brs(1.);
  if(!removebr_) {
    switch(MEPP2VV::process()) {
    case 1: // W+(->e+,nu_e) W-(->e-,nu_ebar) (MCFM: 61 [nproc])
      mcfm_brs *= 0.109338816;
      mcfm_brs *= 0.109338816;
      break;
    case 2: // W+/-(mu+,nu_mu / mu-,nu_mubar) Z(nu_e,nu_ebar) 
            // (MCFM: 72+77 [nproc])
      mcfm_brs *= 0.109338816;
      mcfm_brs *= 0.06839002;
      break;
    case 3: // Z(mu-,mu+) Z(e-,e+) (MCFM: 86 [nproc])
      mcfm_brs *= 0.034616433;
      mcfm_brs *= 0.034616433;
      mcfm_brs *= 2.;  // as identical particle factor 1/2 is now obsolete.
      break;
    case 4: // W+(mu+,nu_mu) Z(nu_e,nu_ebar) (MCFM: 72 [nproc])
      mcfm_brs *= 0.109338816;
      mcfm_brs *= 0.06839002;
      break;
    case 5: // W-(mu-,nu_mubar) Z(nu_e,nu_ebar) (MCFM: 77 [nproc])
      mcfm_brs *= 0.109338816;
      mcfm_brs *= 0.06839002;
      break;
    }
  }
  get_LO_data(output);
  output *= NLOweight();
  output *= mcfm_brs;
  return output;
}

void MEPP2VVPowheg::get_LO_data(double lo_me2) const {

  // In this member we want to get the lo_lumi_ and the lo_me2_ 
  // as these are denominators in the NLO weight. We want sb_,
  // tb_, ub_, k12b_, k22b_ to evaluate the expression for the
  // virtual correction. We also want Yb_, theta1b_, theta2b_, 
  // as these, together with the diboson invariant mass (sb_) 
  // are preserved in the projection from the 2->3 real emission 
  // process to the 2->2 LO/virtual ones i.e. they are the Born 
  // variables which will be needed to construct the 2->3 along
  // with the radiative variables.

  // First a few sanity checks (these can be removed when the code is done):
  if(mePartonData()[0]->id()<0)
    cout << "Error in get_born_variables:\n" 
	 << "mePartonData()[0] is an antiquark, id=" 
	 << ab_->PDGName() << endl;
  if(mePartonData()[1]->id()>0)
    cout << "Error in get_born_variables:\n" 
	 << "mePartonData()[1] is an quark, id=" 
	 << bb_->PDGName() << endl;
  bool alarm(false);
  bool wminus_first(false);
  switch(MEPP2VV::process()) {
  case 1: // W+(->e+,nu_e) W-(->e-,nu_ebar) (MCFM: 61 [nproc])
    if(abs(mePartonData()[2]->id())!=24||abs(mePartonData()[3]->id())!=24) 
      alarm=true;
    if(mePartonData()[2]->id()<0) wminus_first=true;
    break;
  case 2: // W+/-(mu+,nu_mu / mu-,nu_mubar) Z(nu_e,nu_ebar) 
    // (MCFM: 72+77 [nproc])
    if(abs(mePartonData()[2]->id())!=24||mePartonData()[3]->id()!=23) 
      alarm=true;
    break;
  case 3: // Z(mu-,mu+) Z(e-,e+) (MCFM: 86 [nproc])
    if(mePartonData()[2]->id()!= 23||mePartonData()[3]->id()!=23) 
      alarm=true;
    break;
  case 4: // W+(mu+,nu_mu) Z(nu_e,nu_ebar) (MCFM: 72 [nproc])
    if(mePartonData()[2]->id()!= 24||mePartonData()[3]->id()!=23) 
      alarm=true;
    break;
  case 5: // W-(mu-,nu_mubar) Z(nu_e,nu_ebar) (MCFM: 77 [nproc])
    if(mePartonData()[2]->id()!=-24||mePartonData()[3]->id()!=23) 
      alarm=true;
    break;
  }
  if(alarm) { 
    cout << "Error in get_born_variables: unexpected final state labelling.\n";
    cout << "mePartonData()[2] = " << mePartonData()[2]->PDGName() << endl; 
    cout << "mePartonData()[3] = " << mePartonData()[3]->PDGName() << endl; 
  }

  // Now get all data on the LO process needed for the NLO computation:

  // Store the value of the leading order squared matrix element:
  lo_me2_ = lo_me2;

  // BeamParticleData objects for PDF's
  hadron_A_=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().first->dataPtr());
  hadron_B_=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().second->dataPtr());

  // Leading order momentum fractions and associated etabar's:
  xpb_   = lastX1(); 
  etapb_ = sqrt(1.-xpb_);
  xmb_   = lastX2(); 
  etamb_ = sqrt(1.-xmb_);

  // Particle data for incoming QCD particles:
  ab_=mePartonData()[0];  // This is the quark in MEPP2VV.cc
  bb_=mePartonData()[1];  // This is the antiquark in MEPP2VV.cc

  // If necessary swap the particle data vectors so that xpb_, 
  // mePartonData[0], hadron_A_ relate to the inbound particle a
  // (the quark) and xmb_, mePartonData[1], hadron_B_ relate to 
  // particle b (the antiquark). See MEPP2VV.cc for more info. 
  flipped_ = false;
  if(!(lastPartons().first ->dataPtr()==ab_&&
       lastPartons().second->dataPtr()==bb_)) {
    swap(xpb_     ,xmb_     );
    swap(etapb_   ,etamb_   );
    swap(hadron_A_,hadron_B_);
    flipped_ = true;
  }

  // Now get the partonic flux for the Born process:
  lo_lumi_ = hadron_A_->pdf()->xfx(hadron_A_,ab_,scale(),xpb_)/xpb_
           * hadron_B_->pdf()->xfx(hadron_B_,bb_,scale(),xmb_)/xmb_;

  // The Born momenta according to the notation of the FMNR papers,
  // in the diboson centre of mass frame:
  p1b_ = meMomenta()[0];
  p2b_ = meMomenta()[1];
  k1b_ = meMomenta()[2];
  k2b_ = meMomenta()[3];
  // For W+W- events make sure k1 corresponds to the W+ momentum:
  if(wminus_first) swap(k1b_,k2b_);
  // Boost the momenta so they are in the diboson centre of mass frame:
  Lorentz5Momentum ktot(k1b_+k2b_);
  Boost CMSBoostb(-ktot.boostVector());
  p1b_.boost(CMSBoostb);
  p2b_.boost(CMSBoostb);
  k1b_.boost(CMSBoostb);
  k2b_.boost(CMSBoostb);

  // The diboson invariant mass / shat, that and uhat:
  sb_ = sHat();
  tb_ = (p1b_-k1b_).m2();
  ub_ = (p1b_-k2b_).m2();

  // The diboson rapidity:
  Yb_ = 0.5*log(xpb_/xmb_);
  // Note Yb_ = + lastY() if flipped = false 
  // but  Yb_ = - lastY() if flipped = true.
  // Yb_ is always defined with the quark travelling in the +z direction!

  // Masses of the final state bosons:
  k12b_    = k1b_.m2();
  k22b_    = k2b_.m2();

  // Polar and azimuthal angles of the dibosons in their rest frame:
  theta1b_ = acos(k1b_.z()/k1b_.vect().mag());
  theta2b_ = atan2(k1b_.x(),k1b_.y());

  // The masses and angles can be used to reconstruct the momenta in 
  // the diboson centre of mass frame as follows (as a sanity check):
  //  Energy4 lambda = sb_*sb_+k12b_*k12b_+k22b_*k22b_
  //      -2.*sb_*k12b_- 2.*sb_*k22b_ - 2.*k12b_*k22b_;
  //  Energy p  = sqrt(lambda)/2./sqrt(sb_);
  //  Energy pt = p*sin(theta1b_); 
  //  Energy pl = p*cos(theta1b_); 
  //  Lorentz5Momentum p1b(0.*GeV,0.*GeV, p1b_.t(),p1b_.t(),0.*GeV);
  //  Lorentz5Momentum p2b(0.*GeV,0.*GeV,-p2b_.t(),p2b_.t(),0.*GeV);
  //  Lorentz5Momentum k1b( pt*sin(theta2b_), pt*cos(theta2b_), pl,
  //		       sqrt(p*p+k12b_),sqrt(k12b_));
  //  Lorentz5Momentum k2b(-pt*sin(theta2b_),-pt*cos(theta2b_),-pl,
  //		       sqrt(p*p+k22b_),sqrt(k22b_));
  // Checks showed these momenta agreed with p1b_,p2b_,k1b_,k2b_ to within 
  // 0.01 eV. These therefore correspond to theta's of Frixione et al.

  return;
}

double MEPP2VVPowheg::NLOweight() const {
  // If only leading order is required return 1:
  if(contrib_==0) return 1.;

  // Calculate alpha_S and alpha_S/(2*pi)
  alphaS_ = nlo_alphaS_opt_==1 ? fixed_alphaS_ : SM().alphaS(scale());
  double alsOn2pi(alphaS_/2./Constants::pi);

  // Particle data objects for the new plus and minus colliding partons.
  tcPDPtr a_nlo, b_nlo;

  // Calculate the integrand
  double wgt(0.);

  // q qb  contribution
  a_nlo=ab_;
  b_nlo=bb_;
  double wqqbvirt   = Vtilde_universal() + M_V_regular()/lo_me2_;
  double wqqbcollin = alsOn2pi
                    * ( Ctilde_Ltilde_qq_on_x(a_nlo,b_nlo,xt_, 1.) 
                      + Ctilde_Ltilde_qq_on_x(a_nlo,b_nlo,xt_,-1.));
  double wqqbreal   = alsOn2pi
                    * Rtilde_Ltilde_qqb_on_x(a_nlo,b_nlo,xt_,y_);
  double wqqb       = wqqbvirt + wqqbcollin + wqqbreal;
  // q g   contribution
  a_nlo=ab_;
  b_nlo=getParticleData(ParticleID::g);
  double wqgcollin  = alsOn2pi*Ctilde_Ltilde_gq_on_x(a_nlo,b_nlo,xt_,-1.);
  double wqgreal    = alsOn2pi*Rtilde_Ltilde_qg_on_x(a_nlo,b_nlo,xt_,y_);
  double wqg        = wqgreal+wqgcollin;
  // g qb  contribution
  a_nlo=getParticleData(ParticleID::g);
  b_nlo=bb_;
  double wgqbcollin = alsOn2pi*Ctilde_Ltilde_gq_on_x(a_nlo,b_nlo,xt_, 1.);
  double wgqbreal   = alsOn2pi*Rtilde_Ltilde_gqb_on_x(a_nlo,b_nlo,xt_,y_);
  double wgqb       = wgqbreal+wgqbcollin;
  // total contribution
  wgt               = 1.+(wqqb+wgqb+wqg);
  return contrib_==1 ? max(0.,wgt) : max(0.,-wgt);
}

double MEPP2VVPowheg::xbar(double y) const {
  double xbp2(sqr(xpb_)), xbm2(sqr(xmb_));
  double omy(-999.)     , opy( 999.)     ;
  double xbar1(-999.)   , xbar2(-999.)   ;
  if(y== 1.) return xpb_;
  if(y==-1.) return xmb_;
  omy = 1.-y; 
  opy = 1.+y;
  xbar1=2.*opy*xbp2/
    (sqrt(sqr(1.+xbp2)*sqr(omy)+16.*y*xbp2)+omy*(1.-xpb_)*(1.+xpb_));
  xbar2=2.*omy*xbm2/
    (sqrt(sqr(1.+xbm2)*sqr(opy)-16.*y*xbm2)+opy*(1.-xmb_)*(1.+xmb_));
  return max(xbar1,xbar2);
}

double MEPP2VVPowheg::etab(double y) const {
  return sqrt(1.-xbar(y));
}

double MEPP2VVPowheg::xp(double x, double y) const {
  if(x== 1.) return xpb_  ;
  if(y==-1.) return xpb_  ;
  if(y== 1.) return xpb_/x;
  return (xpb_/sqrt(x))*sqrt((2.-(1.-x)*(1.-y))/(2.-(1.-x)*(1.+y)));
}

double MEPP2VVPowheg::xm(double x, double y) const {
  if(x== 1.) return xmb_  ;
  if(y==-1.) return xmb_/x;
  if(y== 1.) return xmb_  ;
  return (xmb_/sqrt(x))*sqrt((2.-(1.-x)*(1.+y))/(2.-(1.-x)*(1.-y)));
}

double MEPP2VVPowheg::Lhat_ab(tcPDPtr a, tcPDPtr b, double x, double y) const {
  double nlo_lumi(-999.);
  double xp_x_y(xp(x,y)),xm_x_y(xm(x,y));
  nlo_lumi = (hadron_A_->pdf()->xfx(hadron_A_,a,scale(),xp_x_y)/xp_x_y)
           * (hadron_B_->pdf()->xfx(hadron_B_,b,scale(),xm_x_y)/xm_x_y);
  return nlo_lumi / lo_lumi_;
}

double MEPP2VVPowheg::Vtilde_universal() const {
  return  alphaS_/2./Constants::pi*CF_ 
        * (   log(sb_/sqr(mu_F()))
	    * (3. + 4.*log(etapb_)+4.*log(etamb_))
	    + 8.*sqr(log(etapb_)) + 8.*sqr(log(etamb_))
	    - 2.*sqr(Constants::pi)/3.
	  )
        + alphaS_/2./Constants::pi*CF_ 
        * ( 8./(1.+y_)*log(etab(y_)/etamb_)
     	  + 8./(1.-y_)*log(etab(y_)/etapb_)
          );
}

double MEPP2VVPowheg::Ctilde_Ltilde_qq_on_x(tcPDPtr a, tcPDPtr b, 
					       double xt, double y ) const {
  if(y!= 1.&&y!=-1.) { cout << "\nCtilde_qq::y value not allowed."; }
  if(y== 1.&&!(abs(a->id())>0&&abs(a->id()<7))) 
    cout << "\nCtilde_qq::for Cqq^plus  a must be a quark! id = " 
	 << a->id() << "\n";
  if(y==-1.&&!(abs(b->id())>0&&abs(b->id()<7))) 
    cout << "\nCtilde_qq::for Cqq^minus b must be a quark! id = " 
	 << b->id() << "\n";
  double x_pm    = x(xt,y);
  double etab_pm = y == 1. ? etapb_ : etamb_ ;
  return ( ( (1./(1.-xt))*log(sb_/sqr(mu_F())/x_pm)+4.*log(etab_pm)/(1.-xt)
       	   + 2.*log(1.-xt)/(1.-xt)
           )*CF_*(1.+sqr(x_pm)) 
	 + sqr(etab_pm)*CF_*(1.-x_pm)
	 )*Lhat_ab(a,b,x_pm,y) / x_pm
       - ( ( (1./(1.-xt))*log(sb_/sqr(mu_F())     )+4.*log(etab_pm)/(1.-xt)
	   + 2.*log(1.-xt)/(1.-xt)
	   )*CF_*2. 
	 );
}

double MEPP2VVPowheg::Ctilde_Ltilde_gq_on_x(tcPDPtr a, tcPDPtr b, 
					       double xt, double y ) const {
  if(y!= 1.&&y!=-1.) { cout << "\nCtilde_gq::y value not allowed."; }
  if(y== 1.&&a->id()!=21)
    cout << "\nCtilde_gq::for Cgq^plus  a must be a gluon! id = " 
	 << a->id() << "\n";
  if(y==-1.&&b->id()!=21) 
    cout << "\nCtilde_gq::for Cgq^minus b must be a gluon! id = " 
	 << b->id() << "\n";
  double x_pm    = x(xt,y);
  double etab_pm = y == 1. ? etapb_ : etamb_ ;
  return ( ( (1./(1.-xt))*log(sb_/sqr(mu_F())/x_pm)+4.*log(etab_pm)/(1.-xt)
	           + 2.*log(1.-xt)/(1.-xt)
	    )*(1.-x_pm)*TR_*(sqr(x_pm)+sqr(1.-x_pm))
	 + sqr(etab_pm)*TR_*2.*x_pm*(1.-x_pm)
	 )*Lhat_ab(a,b,x_pm,y) / x_pm;
}

double MEPP2VVPowheg::Rtilde_Ltilde_qqb_on_x(tcPDPtr a , tcPDPtr b,
				           double  xt, double y ) const {
  return ( ( 
	     1./s(xt ,y  )
	   * t_u_M_R_qqb(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt , 1.)
	   * t_u_M_R_qqb(xt , 1.)*Lhat_ab(a,b,x(xt , 1.), 1.)

	   - 1./s( 1.,y  ) * t_u_M_R_qqb( 1.,y  )

	   + 1./s( 1., 1.) * t_u_M_R_qqb( 1., 1.)

           )*2./(1.-y)/(1.-xt)
	 + ( 
	     1./s(xt ,y  )
	   * t_u_M_R_qqb(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt ,-1.)
	   * t_u_M_R_qqb(xt ,-1.)*Lhat_ab(a,b,x(xt ,-1.),-1.)

	   - 1./s( 1.,y  ) * t_u_M_R_qqb( 1.,y  )

	   + 1./s( 1.,-1.) * t_u_M_R_qqb( 1.,-1.)

           )*2./(1.+y)/(1.-xt)
	 ) / lo_me2_ / 8. / Constants::pi / alphaS_;
}

double MEPP2VVPowheg::Rtilde_Ltilde_gqb_on_x(tcPDPtr a , tcPDPtr b, 
					     double  xt, double y ) const {
  return ( ( 
	     1./s(xt ,y  )
	   * t_u_M_R_gqb(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt , 1.)
	   * t_u_M_R_gqb(xt , 1.)*Lhat_ab(a,b,x(xt , 1.), 1.)

	   - 1./s( 1.,y  ) * t_u_M_R_gqb( 1.,y  )

	   + 1./s( 1., 1.) * t_u_M_R_gqb( 1., 1.)

           )*2./(1.-y)/(1.-xt)
	 + ( 
	     1./s(xt ,y  )
	   * t_u_M_R_gqb(xt ,y  )*Lhat_ab(a,b,x(xt ,y  ),y  )

  	   - 1./s(xt ,-1.)
	   * t_u_M_R_gqb(xt ,-1.)*Lhat_ab(a,b,x(xt ,-1.),-1.)

	   - 1./s( 1.,y  ) * t_u_M_R_gqb( 1.,y  )

	   + 1./s( 1.,-1.) * t_u_M_R_gqb( 1.,-1.)

           )*2./(1.+y)/(1.-xt)
	 ) / lo_me2_ / 8. / Constants::pi / alphaS_;
}

double MEPP2VVPowheg::Rtilde_Ltilde_qg_on_x(tcPDPtr a , tcPDPtr b, 
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
         ) / lo_me2_ / 8. / Constants::pi / alphaS_;
} 

double MEPP2VVPowheg::M_V_regular() const {
  return alphaS_/2./Constants::pi*CF_*
                        (  11./3.
			+  4.*sqr(Constants::pi)/3.
			- (4.*Constants::pi*CF_/CF_)*log(sb_/sqr(mu_UV()))
			)*lo_me2_;
}

Energy2 MEPP2VVPowheg::t_u_M_R_qqb(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*32./9./sqr(sb_)/s(xt,y)*tk(xt,y)*uk(xt,y)
           *( sqr(tk(xt,y)) + sqr(uk(xt,y))
            )*lo_me2_;
}

Energy2 MEPP2VVPowheg::t_u_M_R_qg(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*-4./3./sqr(sb_)*uk(xt,y)
           *( sqr(s(xt,y)) + sqr(uk(xt,y))
            )*lo_me2_;
}

Energy2 MEPP2VVPowheg::t_u_M_R_gqb(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*-4./3./sqr(sb_)*tk(xt,y)
           *( sqr(s(xt,y)) + sqr(tk(xt,y))
            )*lo_me2_;
}




