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
  double xt(-999.);
  double y( -999.);
  if(contrib_>0) {
    // Generate the radiative integration variables:
    xt = *r;
    y  = *(r+1) * 2. - 1.;
  }

  // Continue with lo matrix element code:
  bool output(MEPP2VV::generateKinematics(r));

  // Work out the kinematics for the leading order / virtual process 
  // and also get the leading order luminosity function:
  getKinematics(xt,y);

  return output;
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

  // Store the value of the leading order squared matrix element:
  lo_me2_ = output;
  output *= NLOweight();
  output *= mcfm_brs;
  return output;
}

void MEPP2VVPowheg::getKinematics(double xt, double y) {

  // In this member we want to get the lo_lumi_ as this is a 
  // common denominator in the NLO weight. We want also the 
  // born2to2Kinematics object and all of the real2to3Kinematics
  // objects needed for the NLO weight.

  // First a few sanity checks (these can be removed when the code is done):
  if(mePartonData()[0]->id()<0)
    cout << "Error in get_born_variables:\n" 
	 << "mePartonData()[0] is an antiquark, id=" 
	 << mePartonData()[0]->PDGName() << endl;
  if(mePartonData()[1]->id()>0)
    cout << "Error in get_born_variables:\n" 
	 << "mePartonData()[1] is an quark, id=" 
	 << mePartonData()[1]->PDGName() << endl;
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

  // Should be the hadron containing particle a (the quark):
  hadron_A_=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().first->dataPtr());
  // Should be the hadron containing particle b (the anti-quark):
  hadron_B_=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().second->dataPtr());

  // Leading order momentum fractions:
  double xa(lastX1()); // Should be the quark momentum fraction. 
  double xb(lastX2()); // Should be the anti-quark momentum fraction. 

  // Particle data for incoming QCD particles:
  ab_ = mePartonData()[0];  // This is the quark in MEPP2VV.cc
  bb_ = mePartonData()[1];  // This is the antiquark in MEPP2VV.cc

  // If the lastPartons.first() and lastPartons.second() are 
  // not a quark and antiquark respectively, swap xa<->xb and
  // hadron_A_<->hadron_B_, as xa and xb are defined to be the
  // quark and anti-quark momentum fractions respectively, also,
  // hadron_A_ and hadron_B_, are defined to be the hadrons 
  // containing the colliding partons a and b respectively. 
  // See MEPP2VV.cc for more info. 
  flipped_ = false;
  if(!(lastPartons().first ->dataPtr()==ab_&&
       lastPartons().second->dataPtr()==bb_)) {
    swap(xa       ,xb       );
    swap(hadron_A_,hadron_B_);
    flipped_ = true;
  }

  // Now get the partonic flux for the Born process:
  lo_lumi_ = hadron_A_->pdf()->xfx(hadron_A_,ab_,scale(),xa)/xa
           * hadron_B_->pdf()->xfx(hadron_B_,bb_,scale(),xb)/xb;

  // For W+W- events make sure k1 corresponds to the W+ momentum:
  if(MEPP2VV::process()==1&&wminus_first) swap(meMomenta()[2],meMomenta()[3]);

  // Create the object containing all 2->2 __kinematic__ information:
  B_ = born2to2Kinematics(meMomenta(),xa,xb);

  // Revert momentum swap in case meMomenta and mePartonData correlation
  // needs preserving for other things.
  if(MEPP2VV::process()==1&&wminus_first) swap(meMomenta()[2],meMomenta()[3]);

  // Check the Born kinematics objects is internally consistent:
  B_.sanityCheck();

  // If we are going beyond leading order then lets calculate all of 
  // the necessary real emission kinematics.
  if(contrib_>0) {
    // Soft limit of the 2->3 real emission kinematics:
    S_   = real2to3Kinematics(B_, 1.,  y);
    // Soft-collinear limit of the 2->3 kinematics (emission in +z direction):
    SCp_ = real2to3Kinematics(B_, xt, 1.);
    // Soft-collinear limit of the 2->3 kinematics (emission in -z direction):
    SCm_ = real2to3Kinematics(B_, xt,-1.);
    // Collinear limit of the 2->3 kinematics (emission in +z direction):
    Cp_  = real2to3Kinematics(B_, xt, 1.);
    // Collinear limit of the 2->3 kinematics (emission in -z direction):
    Cm_  = real2to3Kinematics(B_, xt,-1.);
    // The resolved 2->3 real emission kinematics:
    H_   = real2to3Kinematics(B_, xt,  y);

    // Check all the real kinematics objects are internally consistent:
    S_.sanityCheck();
    SCp_.sanityCheck();
    SCm_.sanityCheck();
    Cp_.sanityCheck();
    Cm_.sanityCheck();
    H_.sanityCheck();
  }

  return;
}


double MEPP2VVPowheg::NLOweight() const {
  // If only leading order is required return 1:
  if(contrib_==0) return 1.;

  // Calculate alpha_S and alpha_S/(2*pi)
  alphaS_ = nlo_alphaS_opt_==1 ? fixed_alphaS_ : SM().alphaS(scale());
  double alsOn2pi(alphaS_/2./Constants::pi);

  // Particle data objects for the new plus and minus colliding partons.
  tcPDPtr a_nlo, b_nlo, gluon;
  gluon = getParticleData(ParticleID::g);

  // Calculate the integrand
  double wgt(0.);

  // q qb  contribution
  a_nlo=ab_;
  b_nlo=bb_;
  double wqqbvirt   = Vtilde_universal(S_) + M_V_regular(B_)/lo_me2_;
  double wqqbcollin = Ctilde_Ltilde_qq_on_x(a_nlo,b_nlo,Cp_) 
                    + Ctilde_Ltilde_qq_on_x(a_nlo,b_nlo,Cm_);
  double wqqbreal   = Rtilde_Ltilde_qqb_on_x(a_nlo,b_nlo);
  double wqqb       = wqqbvirt + wqqbcollin + wqqbreal;
  // q g   contribution
  a_nlo=ab_;
  b_nlo=gluon;
  double wqgcollin  = Ctilde_Ltilde_gq_on_x(a_nlo,b_nlo,Cm_);
  double wqgreal    = Rtilde_Ltilde_qg_on_x(a_nlo,b_nlo);
  double wqg        = wqgreal+wqgcollin;
  // g qb  contribution
  a_nlo=gluon;
  b_nlo=bb_;
  double wgqbcollin = Ctilde_Ltilde_gq_on_x(a_nlo,b_nlo,Cp_);
  double wgqbreal   = Rtilde_Ltilde_gqb_on_x(a_nlo,b_nlo);
  double wgqb       = wgqbreal+wgqbcollin;
  // total contribution
  wgt               = 1.+alsOn2pi*(wqqb+wgqb+wqg);
  return contrib_==1 ? max(0.,wgt) : max(0.,-wgt);
}

double MEPP2VVPowheg::Lhat_ab(tcPDPtr a, tcPDPtr b, 
			      real2to3Kinematics Kinematics) const {
    double nlo_lumi(-999.);
  double xp(Kinematics.xpr()),xm(Kinematics.xmr());
  nlo_lumi = (hadron_A_->pdf()->xfx(hadron_A_,a,scale(),xp)/xp)
           * (hadron_B_->pdf()->xfx(hadron_B_,b,scale(),xm)/xm);
  return nlo_lumi / lo_lumi_;
}

double MEPP2VVPowheg::Vtilde_universal(real2to3Kinematics S) const {
  double xbar_y = S.xbar();
  double y = S.y();
  double etapb(S.bornVariables().etapb());
  double etamb(S.bornVariables().etamb());
  Energy2 sb(S.s2r());
  return  alphaS_/2./Constants::pi*CF_ 
        * (   log(sb/sqr(mu_F()))
	    * (3. + 4.*log(etapb)+4.*log(etamb))
	    + 8.*sqr(log(etapb)) +8.*sqr(log(etamb))
	    - 2.*sqr(Constants::pi)/3.
	  )
        + alphaS_/2./Constants::pi*CF_ 
        * ( 8./(1.+y)*log(sqrt(1.-xbar_y)/etamb)
     	  + 8./(1.-y)*log(sqrt(1.-xbar_y)/etapb)
          );
}

double MEPP2VVPowheg::Ctilde_Ltilde_qq_on_x(tcPDPtr a, tcPDPtr b, 
					    real2to3Kinematics C) const {
  if(C.y()!= 1.&&C.y()!=-1.) 
    cout << "\nCtilde_qq::y value not allowed.";
  if(C.y()== 1.&&!(abs(a->id())>0&&abs(a->id()<7))) 
    cout << "\nCtilde_qq::for Cqq^plus  a must be a quark! id = " 
	 << a->id() << "\n";
  if(C.y()==-1.&&!(abs(b->id())>0&&abs(b->id()<7))) 
    cout << "\nCtilde_qq::for Cqq^minus b must be a quark! id = " 
	 << b->id() << "\n";
  double xt = C.xt();
  double x  = C.xr();
  double etab = C.y() == 1. ? C.bornVariables().etapb() 
                            : C.bornVariables().etamb() ;
  Energy2 sb(C.s2r());
  return ( ( (1./(1.-xt))*log(sb/sqr(mu_F())/x)+4.*log(etab)/(1.-xt)
       	   + 2.*log(1.-xt)/(1.-xt)
           )*CF_*(1.+sqr(x)) 
	 + sqr(etab)*CF_*(1.-x)
	 )*Lhat_ab(a,b,C) / x
       - ( ( (1./(1.-xt))*log(sb/sqr(mu_F())  )+4.*log(etab)/(1.-xt)
	   + 2.*log(1.-xt)/(1.-xt)
	   )*CF_*2. 
	 );
}

double MEPP2VVPowheg::Ctilde_Ltilde_gq_on_x(tcPDPtr a, tcPDPtr b, 
					    real2to3Kinematics C) const {
  if(C.y()!= 1.&&C.y()!=-1.) 
    cout << "\nCtilde_gq::y value not allowed.";
  if(C.y()== 1.&&a->id()!=21)
    cout << "\nCtilde_gq::for Cgq^plus  a must be a gluon! id = " 
	 << a->id() << "\n";
  if(C.y()==-1.&&b->id()!=21) 
    cout << "\nCtilde_gq::for Cgq^minus b must be a gluon! id = " 
	 << b->id() << "\n";
  double xt = C.xt();
  double x  = C.xr();
  double etab = C.y() == 1. ? C.bornVariables().etapb() 
                            : C.bornVariables().etamb() ;
  Energy2 sb(C.s2r());
  return ( ( (1./(1.-xt))*log(sb/sqr(mu_F())/x)+4.*log(etab)/(1.-xt)
	   + 2.*log(1.-xt)/(1.-xt)
	   )*(1.-x)*TR_*(sqr(x)+sqr(1.-x))
	 + sqr(etab)*TR_*2.*x*(1.-x)
	 )*Lhat_ab(a,b,C) / x;
}

double MEPP2VVPowheg::Rtilde_Ltilde_qqb_on_x(tcPDPtr a , tcPDPtr b) const {
  double xt(H_.xt()); 
  double y(H_.y()); 
  Energy2 s(H_.sr());
  Energy2 s2(H_.s2r());
  return 
     ( ( (t_u_M_R_qqb(H_)*Lhat_ab(a,b,H_)-t_u_M_R_qqb(Cp_)*Lhat_ab(a,b,Cp_))/s
       - (t_u_M_R_qqb(S_) - t_u_M_R_qqb(SCp_))/s2

       )*2./(1.-y)/(1.-xt)

     + ( (t_u_M_R_qqb(H_)*Lhat_ab(a,b,H_)-t_u_M_R_qqb(Cm_)*Lhat_ab(a,b,Cm_))/s
       - (t_u_M_R_qqb(S_) - t_u_M_R_qqb(SCm_))/s2

       )*2./(1.+y)/(1.-xt)

     ) / lo_me2_ / 8. / Constants::pi / alphaS_;
}

double MEPP2VVPowheg::Rtilde_Ltilde_gqb_on_x(tcPDPtr a , tcPDPtr b) const {
  double xt(H_.xt()); 
  double y(H_.y()); 
  Energy2 s(H_.sr());
  Energy2 s2(H_.s2r());
  return 
     ( ( (t_u_M_R_gqb(H_)*Lhat_ab(a,b,H_)-t_u_M_R_gqb(Cp_)*Lhat_ab(a,b,Cp_))/s
       - (t_u_M_R_gqb(S_) - t_u_M_R_gqb(SCp_)                              )/s2

       )*2./(1.-y)/(1.-xt)

     + ( (t_u_M_R_gqb(H_)*Lhat_ab(a,b,H_)-t_u_M_R_gqb(Cm_)*Lhat_ab(a,b,Cm_))/s
       - (t_u_M_R_gqb(S_) - t_u_M_R_gqb(SCm_)                              )/s2

       )*2./(1.+y)/(1.-xt)

     ) / lo_me2_ / 8. / Constants::pi / alphaS_;
}

double MEPP2VVPowheg::Rtilde_Ltilde_qg_on_x(tcPDPtr a , tcPDPtr b) const {
  double xt(H_.xt()); 
  double y(H_.y()); 
  Energy2 s(H_.sr());
  Energy2 s2(H_.s2r());
  return 
     ( ( (t_u_M_R_qg(H_)*Lhat_ab(a,b,H_) - t_u_M_R_qg(Cp_)*Lhat_ab(a,b,Cp_))/s
       - (t_u_M_R_qg(S_)                 - t_u_M_R_qg(SCp_)                )/s2

       )*2./(1.-y)/(1.-xt)

     + ( (t_u_M_R_qg(H_)*Lhat_ab(a,b,H_) - t_u_M_R_qg(Cm_)*Lhat_ab(a,b,Cm_))/s
       - (t_u_M_R_qg(S_)                 - t_u_M_R_qg(SCm_)                )/s2

       )*2./(1.+y)/(1.-xt)

     ) / lo_me2_ / 8. / Constants::pi / alphaS_;
}

double MEPP2VVPowheg::M_V_regular(born2to2Kinematics S) const {
  Energy2 s(S.sb());
  Energy2 t(S.tb());
  Energy2 u(S.ub());
  Energy2 mW2(S.k12b());
  Energy2 mZ2(S.k22b());
  return alphaS_/2./Constants::pi*CF_*
                        (  11./3.
			+  4.*sqr(Constants::pi)/3.
			- (4.*Constants::pi*CF_/CF_)*log(s/sqr(mu_UV()))
			)*lo_me2_;
}

Energy2 MEPP2VVPowheg::t_u_M_R_qqb(real2to3Kinematics theKinematics) const {
  // First the Born variables:
  Energy2 s2(theKinematics.s2r());
  Energy2 k12(theKinematics.k12r());
  Energy2 k22(theKinematics.k22r());
  double  theta1(theKinematics.theta1r());
  double  theta2(theKinematics.theta2r());
  // Then the rest:
  Energy2 s(theKinematics.sr());
  Energy2 tk(theKinematics.tkr());
  Energy2 uk(theKinematics.ukr());
  double  cpsi(theKinematics.cpsir());
  double  cpsipr(theKinematics.cpsipr());
  double  betax(theKinematics.betaxr());
  double  v1(theKinematics.v1r());
  double  v2(theKinematics.v2r());
  Energy2 q1(theKinematics.q1r());
  Energy2 q2(theKinematics.q2r());
  Energy2 q1hat(theKinematics.q1hatr());
  Energy2 q2hat(theKinematics.q2hatr());
  Energy2 w1(theKinematics.w1r());
  Energy2 w2(theKinematics.w2r());

  return 0.*GeV2;
}

Energy2 MEPP2VVPowheg::t_u_M_R_qg(real2to3Kinematics theKinematics) const {
  // First the Born variables:
  Energy2 s2(theKinematics.s2r());
  Energy2 k12(theKinematics.k12r());
  Energy2 k22(theKinematics.k22r());
  double  theta1(theKinematics.theta1r());
  double  theta2(theKinematics.theta2r());
  // Then the rest:
  Energy2 s(theKinematics.sr());
  Energy2 tk(theKinematics.tkr());
  Energy2 uk(theKinematics.ukr());
  double  cpsi(theKinematics.cpsir());
  double  cpsipr(theKinematics.cpsipr());
  double  betax(theKinematics.betaxr());
  double  v1(theKinematics.v1r());
  double  v2(theKinematics.v2r());
  Energy2 q1(theKinematics.q1r());
  Energy2 q2(theKinematics.q2r());
  Energy2 q1hat(theKinematics.q1hatr());
  Energy2 q2hat(theKinematics.q2hatr());
  Energy2 w1(theKinematics.w1r());
  Energy2 w2(theKinematics.w2r());

  return 0.*GeV2;
}

Energy2 MEPP2VVPowheg::t_u_M_R_gqb(real2to3Kinematics theKinematics) const {
  // First the Born variables:
  Energy2 s2(theKinematics.s2r());
  Energy2 k12(theKinematics.k12r());
  Energy2 k22(theKinematics.k22r());
  double  theta1(theKinematics.theta1r());
  double  theta2(theKinematics.theta2r());
  // Then the rest:
  Energy2 s(theKinematics.sr());
  Energy2 tk(theKinematics.tkr());
  Energy2 uk(theKinematics.ukr());
  double  cpsi(theKinematics.cpsir());
  double  cpsipr(theKinematics.cpsipr());
  double  betax(theKinematics.betaxr());
  double  v1(theKinematics.v1r());
  double  v2(theKinematics.v2r());
  Energy2 q1(theKinematics.q1r());
  Energy2 q2(theKinematics.q2r());
  Energy2 q1hat(theKinematics.q1hatr());
  Energy2 q2hat(theKinematics.q2hatr());
  Energy2 w1(theKinematics.w1r());
  Energy2 w2(theKinematics.w2r());

  return 0.*GeV2;
}




