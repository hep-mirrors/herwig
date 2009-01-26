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
    CF_(4./3.),    TR_(0.5),  NC_(3.),
    contrib_(1),   nlo_alphaS_opt_(0) , fixed_alphaS_(0.118109485),
    removebr_(1) {  
    massOption(true ,1);
    massOption(false,1);
}

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
    xt = (*(r+1));
    y  = (*(r+2)) * 2. - 1.;
//    xt = UseRandom::rnd();
//    y  = UseRandom::rnd() * 2. -1.;
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
    SCp_ = real2to3Kinematics(B_, 1., 1.);
    // Soft-collinear limit of the 2->3 kinematics (emission in -z direction):
    SCm_ = real2to3Kinematics(B_, 1.,-1.);
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

  // Calculate alpha_S and alpha_S/(2*pi).
  alphaS_ = nlo_alphaS_opt_==1 ? fixed_alphaS_ : SM().alphaS(scale());
  double alsOn2pi(alphaS_/2./pi);

  // Particle data objects for the new plus and minus colliding partons.
  tcPDPtr a_nlo, b_nlo, gluon;
  gluon = getParticleData(ParticleID::g);

  // Get the all couplings.
  double gW(sqrt(4.0*pi*SM().alphaEM(scale())/SM().sin2ThetaW()));
  double sin2ThetaW(SM().sin2ThetaW());
  double cosThetaW(sqrt(1.-sin2ThetaW));
  guL_ = gW/2./cosThetaW*( 1.-4./3.*sin2ThetaW);
  gdL_ = gW/2./cosThetaW*(-1.+2./3.*sin2ThetaW);
  eZ_  = gW*cosThetaW;

  // Get the CKM entry. Note that this code was debugged 
  // considerably; the call to CKM(particle,particle) 
  // did not appear to work, so we extract the elements
  // as follows below. The right numbers now appear to 
  // to be associated with the right quarks.
  double Kij(-999.);
  int up_id(-999),dn_id(-999);
  if(abs(ab_->id())%2==0&&abs(bb_->id())%2==1) {
    up_id = abs(ab_->id());  
    dn_id = abs(bb_->id());  
  }
  else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==0) {
    up_id = abs(bb_->id());  
    dn_id = abs(ab_->id());  
  }
  else {
    cout << "MEPP2VVPowheg::NLOweight" << endl;
    cout << "No quarks in the call to the CKM matrix!" << endl;
  }
  up_id /= 2;
  up_id -= 1;
  dn_id -= 1;
  dn_id /= 2;
  Kij = sqrt(SM().CKM(up_id,dn_id));

  Fij2_ = sqr(gW/2./sqrt(2.)*Kij);

  // Calculate the integrand
  double wgt(0.);

  // q qb  contribution
  a_nlo=ab_;
  b_nlo=bb_;

  ///////////////////////////////////////
  ///////////////////////////////////////
  // DEBUGGING - switching off TGCs    //
  // eZ_=0.;
  // guL_=0.;
  ///////////////////////////////////////
  ///////////////////////////////////////

  double wqqbvirt   = Vtilde_universal(S_) + M_V_regular(S_)
                                           / lo_me2_;
  double wqqbcollin = alsOn2pi*( Ctilde_Ltilde_qq_on_x(a_nlo,b_nlo,Cp_) 
                               + Ctilde_Ltilde_qq_on_x(a_nlo,b_nlo,Cm_) );
  double wqqbreal   = alsOn2pi*Rtilde_Ltilde_qqb_on_x(a_nlo,b_nlo);
  double wqqb       = wqqbvirt + wqqbcollin + wqqbreal;
  // q g   contribution
  a_nlo=ab_;
  b_nlo=gluon;
  double wqgcollin  = alsOn2pi*Ctilde_Ltilde_gq_on_x(a_nlo,b_nlo,Cm_);
  double wqgreal    = alsOn2pi*Rtilde_Ltilde_qg_on_x(a_nlo,b_nlo);
  double wqg        = wqgreal + wqgcollin;
  // g qb  contribution
  a_nlo=gluon;
  b_nlo=bb_;
  double wgqbcollin = alsOn2pi*Ctilde_Ltilde_gq_on_x(a_nlo,b_nlo,Cp_);
  double wgqbreal   = alsOn2pi*Rtilde_Ltilde_gqb_on_x(a_nlo,b_nlo);
  double wgqb       = wgqbreal+wgqbcollin;
  // total contribution
  wgt               = 1.+(wqqb+wgqb+wqg);
  // Temporary warning in case of nans & infs (not getting any so far!):
  if(isnan(wgt)||isinf(wgt)) { 
    cout << "NAN / INF detected: wgt = " << wgt << endl;
    cout << "Resetting wgt to zero." << endl;
    wgt = 0.;
  }
//  Debugging output:
//     cout << "\n\n\n";
//     cout << ab_->PDGName() << ", " << bb_->PDGName() << endl;
//     cout << "lo_me2_ " << lo_me2_ << ",  " << "M_Born " << M_Born(B_) << "  ratio " << M_Born(B_)/lo_me2_ << endl;
//     cout << "xt = " << H_.xt() << endl;
//     cout << "xr = " << H_.xr() << ", y = " << H_.y() << endl;
//     cout << "sb + tkb + ukb = "
// 	 << B_.sb()/GeV2  << " + "
// 	 << B_.tb()/GeV2 << " + "
// 	 << B_.ub()/GeV2 << " = "
// 	 << (B_.sb()+B_.tb()+B_.ub())/GeV2 << endl;
//     cout << "sr + tkr + ukr = "
// 	 << H_.sr()/GeV2  << " + "
// 	 << H_.tkr()/GeV2 << " + "
// 	 << H_.ukr()/GeV2 << " = "
// 	 << (H_.sr()+H_.tkr()+H_.ukr())/GeV2 << endl;
//     cout << "s2r        " << H_.s2r()/GeV2 << endl;
//     cout << "sqrt(k12)  " << sqrt(H_.k12r())/GeV << endl;
//     cout << "sqrt(k22)  " << sqrt(H_.k22r())/GeV << endl;
//     cout << "gW^2       " << gW*gW      << endl;
//     cout << "sin2ThetaW " << sin2ThetaW << endl;
//     cout << "sqr(Kij)   " << Kij*Kij    << endl;
//     cout << "wqqbvirt   " << wqqbvirt   << endl;
//     cout << "wqqbcollin " << wqqbcollin << endl;
//     cout << "wqqbreal   " << wqqbreal   << endl;
//     cout << "wqqb       " << wqqb       << endl;
//     cout << "t_u_M_R_qqb(H_)    " << t_u_M_R_qqb(H_)  /GeV2 << endl;
//     cout << "t_u_M_R_qqb(Cp_)   " << t_u_M_R_qqb(Cp_) /GeV2 << endl;
//     cout << "t_u_M_R_qqb(Cm_)   " << t_u_M_R_qqb(Cm_) /GeV2 << endl;
//     cout << "t_u_M_R_qqb(S_)    " << t_u_M_R_qqb(S_)  /GeV2 << endl;
//     cout << "t_u_M_R_qqb(SCp_)  " << t_u_M_R_qqb(SCp_)/GeV2 << endl;
//     cout << "t_u_M_R_qqb(SCm_)  " << t_u_M_R_qqb(SCm_)/GeV2 << endl;
//     cout << "Splitting approx plus:  " 
// 	 << Cp_.sr()*8.*pi*alphaS_/Cp_.xr()
// 	      *CF_*(1.+sqr(Cp_.xr()))*M_Born(B_)  /GeV2
// 	 << endl;
//     cout << "Splitting approx minus:  " 
// 	 << Cm_.sr()*8.*pi*alphaS_/Cm_.xr()
// 	      *CF_*(1.+sqr(Cm_.xr()))*M_Born(B_)  /GeV2
// 	 << endl;
//     cout << "( ( (t_u_M_R_qqb(H_)*Lhat_ab(a,b,H_) - t_u_M_R_qqb(Cp_)*Lhat_ab(ab_,bb_,Cp_))/s)*2./(1.-y)/(1.-xt) ) / lo_me2_ / 8. / pi / alphaS_ \n"
// 	 << ( ( (t_u_M_R_qqb(H_)*Lhat_ab(ab_,bb_,H_) - t_u_M_R_qqb(Cp_)*Lhat_ab(ab_,bb_,Cp_))/H_.sr())*2./(1.-H_.y())/(1.-H_.xt()) ) / lo_me2_ / 8. / pi / alphaS_ 
// 	 << endl;
//     cout << "( ( - (t_u_M_R_qqb(S_)                 - t_u_M_R_qqb(SCp_)                )/s2)*2./(1.-y)/(1.-xt) ) / lo_me2_ / 8. / pi / alphaS_ \n"
// 	 << ( ( - (t_u_M_R_qqb(S_)                 - t_u_M_R_qqb(SCp_)                )/H_.s2r())*2./(1.-H_.y())/(1.-H_.xt()) ) / lo_me2_ / 8. / pi / alphaS_
// 	 << endl;
//     cout << "( ( (t_u_M_R_qqb(H_)*Lhat_ab(ab_,bb_,H_) - t_u_M_R_qqb(Cm_)*Lhat_ab(ab_,bb_,Cm_))/s)*2./(1.+y)/(1.-xt) ) / lo_me2_ / 8. / pi / alphaS_ \n"
// 	 << ( ( (t_u_M_R_qqb(H_)*Lhat_ab(ab_,bb_,H_) - t_u_M_R_qqb(Cm_)*Lhat_ab(ab_,bb_,Cm_))/H_.sr())*2./(1.+H_.y())/(1.-H_.xt()) ) / lo_me2_ / 8. / pi / alphaS_ 
// 	 << endl;
//     cout << "( ( - (t_u_M_R_qqb(S_)                 - t_u_M_R_qqb(SCm_)                )/s2)*2./(1.+y)/(1.-xt) ) / lo_me2_ / 8. / pi / alphaS_ \n"
// 	 << ( ( - (t_u_M_R_qqb(S_)                 - t_u_M_R_qqb(SCm_)                )/H_.s2r())*2./(1.+H_.y())/(1.-H_.xt()) ) / lo_me2_ / 8. / pi / alphaS_
// 	 << endl;
//     cout << "wqgcollin  " << wqgcollin  << endl;
//     cout << "wqgreal    " << wqgreal    << endl;
//     cout << "wqg        " << wqg        << endl;
//     cout << "wgqbcollin " << wgqbcollin << endl;
//     cout << "wgqbreal   " << wgqbreal   << endl;
//     cout << "wgqb       " << wgqb       << endl;
//     cout << "wgt        " << wgt        << endl;
  return contrib_==1 ? max(0.,wgt) : max(0.,-wgt);
}

double MEPP2VVPowheg::Lhat_ab(tcPDPtr a, tcPDPtr b, 
			      real2to3Kinematics Kinematics) const {
  if(!(abs(a->id())<=6||a->id()==21)||!(abs(b->id())<=6||b->id()==21))
    cout << "MEPP2VVPowheg::Lhat_ab: Error," 
         << "particle a = " << a->PDGName() << ", "
         << "particle b = " << b->PDGName() << endl;
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
  return  alphaS_/2./pi*CF_ 
        * (   log(sb/sqr(mu_F()))
	    * (3. + 4.*log(etapb)+4.*log(etamb))
	    + 8.*sqr(log(etapb)) +8.*sqr(log(etamb))
	    - 2.*sqr(pi)/3.
	  )
        + alphaS_/2./pi*CF_ 
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
  if(!(abs(a->id())<=6||a->id()==21)||!(abs(b->id())<=6||b->id()==21))
    cout << "MEPP2VVPowheg::Rtilde_Ltilde_qqb_on_x: Error," 
         << "particle a = " << a->PDGName() << ", "
         << "particle b = " << b->PDGName() << endl;
  double xt(H_.xt()); 
  double y(H_.y()); 
  Energy2 s(H_.sr());
  Energy2 s2(H_.s2r());
  return 
   ( ( (t_u_M_R_qqb(H_)*Lhat_ab(a,b,H_) - t_u_M_R_qqb(Cp_)*Lhat_ab(a,b,Cp_))/s
     - (t_u_M_R_qqb(S_)                 - t_u_M_R_qqb(SCp_)                )/s2

     )*2./(1.-y)/(1.-xt)

   + ( (t_u_M_R_qqb(H_)*Lhat_ab(a,b,H_) - t_u_M_R_qqb(Cm_)*Lhat_ab(a,b,Cm_))/s
     - (t_u_M_R_qqb(S_)                 - t_u_M_R_qqb(SCm_)                )/s2

     )*2./(1.+y)/(1.-xt)

   ) / lo_me2_ / 8. / pi / alphaS_;
}

double MEPP2VVPowheg::Rtilde_Ltilde_gqb_on_x(tcPDPtr a , tcPDPtr b) const {
  if(!(abs(a->id())<=6||a->id()==21)||!(abs(b->id())<=6||b->id()==21))
    cout << "MEPP2VVPowheg::Rtilde_Ltilde_gqb_on_x: Error," 
         << "particle a = " << a->PDGName() << ", "
         << "particle b = " << b->PDGName() << endl;
  double xt(H_.xt()); 
  double y(H_.y()); 
  Energy2 s(H_.sr());
  Energy2 s2(H_.s2r());
  return 
   ( ( (t_u_M_R_gqb(H_)*Lhat_ab(a,b,H_) - t_u_M_R_gqb(Cp_)*Lhat_ab(a,b,Cp_))/s
     - (t_u_M_R_gqb(S_)                 - t_u_M_R_gqb(SCp_)                )/s2

     )*2./(1.-y)/(1.-xt)

   + ( (t_u_M_R_gqb(H_)*Lhat_ab(a,b,H_) - t_u_M_R_gqb(Cm_)*Lhat_ab(a,b,Cm_))/s
     - (t_u_M_R_gqb(S_)                 - t_u_M_R_gqb(SCm_)                )/s2

     )*2./(1.+y)/(1.-xt)

   ) / lo_me2_ / 8. / pi / alphaS_;
}

double MEPP2VVPowheg::Rtilde_Ltilde_qg_on_x(tcPDPtr a , tcPDPtr b) const {
  if(!(abs(a->id())<=6||a->id()==21)||!(abs(b->id())<=6||b->id()==21))
    cout << "MEPP2VVPowheg::Rtilde_Ltilde_qg_on_x: Error," 
         << "particle a = " << a->PDGName() << ", "
         << "particle b = " << b->PDGName() << endl;
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

   ) / lo_me2_ / 8. / pi / alphaS_;
}

/***************************************************************************/
// The following three functions are identically  \tilde{I}_{4,t}, 
// \tilde{I}_{3,WZ} and \tilde{I}_{3,W} given in Eqs. B.8,B.9,B.10 
// of NPB 383(1992)3-44, respectively. They are related to / derived 
// from the loop integrals in Eqs. A.3, A.5 and A.8 of the same paper.
InvEnergy4 TildeI4t(Energy2 s, Energy2 t, Energy2 mW2, Energy2 mZ2);
InvEnergy2 TildeI3WZ(Energy2 s, Energy2 mW2, Energy2 mZ2, double beta);
InvEnergy2 TildeI3W(Energy2 s, Energy2 t, Energy2 mW2);

/***************************************************************************/
// The following six functions are identically  I_{dd}^{(1)}, I_{ud}^{(1)}, 
// I_{uu}^{(1)}, F_{u}^{(1)}, F_{d}^{(1)}, H^{(1)} from Eqs. B.4, B.5, B.3,
// B.3, B.6, B.7 of NPB 383(1992)3-44, respectively. They make up the 
// one-loop matrix element. Ixx functions correspond to the graphs
// with no TGC, Fx functions are due to non-TGC graphs interfering 
// with TGC graphs, while the H function is due purely to TGC graphs. 
double Idd1(Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2,double beta);
double Iud1(Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2,double beta);
double Iuu1(Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2,double beta);
Energy2 Fu1(Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2,double beta);
Energy2 Fd1(Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2,double beta);
Energy4 H1 (Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2);

/***************************************************************************/
// M_V_Regular is the regular part of the one-loop matrix element 
// exactly as defined in Eqs. B.1 and B.2 of of NPB 383(1992)3-44.
double MEPP2VVPowheg::M_V_regular(real2to3Kinematics S) const {
  Energy2 s(S.bornVariables().sb());
  Energy2 t(S.bornVariables().tb());
  Energy2 u(S.bornVariables().ub());
  Energy2 mW2(S.k12r()); // N.B. the diboson masses are preserved in getting
  Energy2 mZ2(S.k22r()); // the 2->2 from the 2->3 kinematics.
  double  beta(S.betaxr()); // N.B. for x=1 \beta_x=\beta in NPB 383(1992)3-44.
  
  return 4.*pi*alphaS_*Fij2_*CF_*(1./sqr(4.*pi))/NC_
       * ( gdL_*gdL_*Idd1(s,t,u,mW2,mZ2,beta)
	 + gdL_*guL_*Iud1(s,t,u,mW2,mZ2,beta)
	 + guL_*guL_*Iuu1(s,t,u,mW2,mZ2,beta) 
	 - eZ_/(s-mW2) * ( gdL_*Fd1(s,t,u,mW2,mZ2,beta)
	                 - guL_*Fu1(s,t,u,mW2,mZ2,beta)
	                 )
         + sqr(eZ_/(s-mW2)) * H1(s,t,u,mW2,mZ2)
	 );
}

/***************************************************************************/
InvEnergy4 TildeI4t(Energy2 s, Energy2 t, Energy2 mW2, Energy2 mZ2) {
  double sqrBrackets;
  sqrBrackets = (  sqr(log(-t/mW2))/2.+log(-t/mW2)*log(-t/mZ2)/2.
	        - 2.*log(-t/mW2)*log((mW2-t)/mW2)-2.*ReLi2(t/mW2)
	        );

  swap(mW2,mZ2);
  sqrBrackets+= (  sqr(log(-t/mW2))/2.+log(-t/mW2)*log(-t/mZ2)/2.
	        - 2.*log(-t/mW2)*log((mW2-t)/mW2)-2.*ReLi2(t/mW2)
	        );
  swap(mW2,mZ2);

  return sqrBrackets/s/t;
}

InvEnergy2 TildeI3WZ(Energy2 s, Energy2 mW2, Energy2 mZ2, double beta) {
  Energy2 sig(mZ2+mW2);
  Energy2 del(mZ2-mW2);
  double  sqrBrackets ;
  sqrBrackets  = ( ReLi2(2.*mW2/(sig-del*(del/s+beta)))
	 	 + ReLi2((1.-del/s+beta)/2.)
	         + sqr(log((1.-del/s+beta)/2.))/2.
		 + log((1.-del/s-beta)/2.)*log((1.+del/s-beta)/2.)
                 );

  beta *= -1;
  sqrBrackets -= ( ReLi2(2.*mW2/(sig-del*(del/s+beta)))
	 	 + ReLi2((1.-del/s+beta)/2.)
	         + sqr(log((1.-del/s+beta)/2.))/2.
		 + log((1.-del/s-beta)/2.)*log((1.+del/s-beta)/2.)
                 );
  beta *= -1;

  swap(mW2,mZ2);
  del  *= -1.;
  sqrBrackets += ( ReLi2(2.*mW2/(sig-del*(del/s+beta)))
	 	 + ReLi2((1.-del/s+beta)/2.)
	         + sqr(log((1.-del/s+beta)/2.))/2.
		 + log((1.-del/s-beta)/2.)*log((1.+del/s-beta)/2.)
                 );
  swap(mW2,mZ2);
  del  *= -1.;

  beta *= -1;
  swap(mW2,mZ2);
  del  *= -1.;
  sqrBrackets -= ( ReLi2(2.*mW2/(sig-del*(del/s+beta)))
	 	 + ReLi2((1.-del/s+beta)/2.)
	         + sqr(log((1.-del/s+beta)/2.))/2.
		 + log((1.-del/s-beta)/2.)*log((1.+del/s-beta)/2.)
                 );
  beta *= -1;
  swap(mW2,mZ2);
  del  *= -1.;

  return sqrBrackets/s/beta;
}

InvEnergy2 TildeI3W(Energy2 s, Energy2 t, Energy2 mW2) {
  return 
    1./(mW2-t)*(sqr(log(mW2/s))/2.-sqr(log(-t/s))/2.-sqr(pi)/2.);
}

/***************************************************************************/
double  Idd1(Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2, double beta) { 
    Energy2 sig(mZ2+mW2);
    Energy2 del(mZ2-mW2);
    double Val(0.);

    Val +=    2.*(22.*t*t+t*(19.*s+18.*sig)+18.*mW2*mZ2)/t/t     
            - 8.*(u*t+2*s*sig)/mW2/mZ2
            - 2.*sqr(t-u)/t/s/sqr(beta);

    Val += +( 2.*(8.*t*t+4.*t*(s-3.*sig)+4*sqr(sig)-5.*s*sig+s*s)/t/s/sqr(beta)
	    + 4.*(t*(3.*u+s)-3.*mW2*mZ2)/t/t
	    + 6.*(t+u)*sqr(t-u)/t/s/s/sqr(sqr(beta))
	    )*log(-t/s);

    Val += +( ( 8.*t*t*(-2.*s+del)+8.*t*(-s*s+3.*s*sig-2.*del*sig)
	       - 2.*(s-sig)*(s*s-4.*s*sig+3.*del*sig)
	      )/t/s/s/beta/beta
	      + 16.*s*(t-mZ2)/(t*(u+s)-mW2*mZ2)
	      + 2.*(4.*t*t+t*(10.*s-3.*mZ2-9.*mW2)+12.*mW2*mZ2)/t/t
	      -6.*(s-del)*(t+u)*sqr(t-u)/t/s/s/s/sqr(sqr(beta))
            )*log(-t/mW2);

    Val +=  ( - ( 4.*t*t*(2.*sig-3.*s)
	        - 4.*t*(s-sig)*(2.*s-3.*sig)
	        - 2.*(s-2.*sig)*sqr(s-sig)
	        )/t/s/beta/beta
	      + ( 4.*sig*t-3.*s*s+4.*s*sig
	        - 4.*(mW2*mW2+mZ2*mZ2)
	        )/t
	      - 3.*sqr(t*t-u*u)/t/s/s/sqr(sqr(beta))
	    )*TildeI3WZ(s,mW2,mZ2,beta);

    Val += +( 4.*(t*u+2.*s*sig)/3./mW2/mZ2 - 4.*(t-2.*u)/3./t
	    )*pi*pi;

    Val += -( 4.*s*(t*u-2.*mW2*mZ2)/t
	    )*TildeI4t(s,t,mW2,mZ2);

    Val +=  ( 8.*(t-mW2)*(u*t-2.*mW2*mZ2)/t/t
	    )*TildeI3W(s,t,mW2);

    swap(mW2,mZ2);
    del *= -1;
    Val +=    2.*(22.*t*t+t*(19.*s+18.*sig)+18.*mW2*mZ2)/t/t     
            - 8.*(u*t+2*s*sig)/mW2/mZ2
            - 2.*sqr(t-u)/t/s/sqr(beta);

    Val += +( 2.*(8.*t*t+4.*t*(s-3.*sig)+4*sqr(sig)-5.*s*sig+s*s)/t/s/sqr(beta)
	    + 4.*(t*(3.*u+s)-3.*mW2*mZ2)/t/t
	    + 6.*(t+u)*sqr(t-u)/t/s/s/sqr(sqr(beta))
	    )*log(-t/s);

    Val += +( ( 8.*t*t*(-2.*s+del)+8.*t*(-s*s+3.*s*sig-2.*del*sig)
	       - 2.*(s-sig)*(s*s-4.*s*sig+3.*del*sig)
	      )/t/s/s/beta/beta
	      + 16.*s*(t-mZ2)/(t*(u+s)-mW2*mZ2)
	      + 2.*(4.*t*t+t*(10.*s-3.*mZ2-9.*mW2)+12.*mW2*mZ2)/t/t
	      -6.*(s-del)*(t+u)*sqr(t-u)/t/s/s/s/sqr(sqr(beta))
            )*log(-t/mW2);

    Val +=  ( - ( 4.*t*t*(2.*sig-3.*s)
	        - 4.*t*(s-sig)*(2.*s-3.*sig)
	        - 2.*(s-2.*sig)*sqr(s-sig)
	        )/t/s/beta/beta
	      + ( 4.*sig*t-3.*s*s+4.*s*sig
	        - 4.*(mW2*mW2+mZ2*mZ2)
	        )/t
	      - 3.*sqr(t*t-u*u)/t/s/s/sqr(sqr(beta))
	    )*TildeI3WZ(s,mW2,mZ2,beta);

    Val += +( 4.*(t*u+2.*s*sig)/3./mW2/mZ2 - 4.*(t-2.*u)/3./t
	    )*pi*pi;

    Val += -( 4.*s*(t*u-2.*mW2*mZ2)/t
	    )*TildeI4t(s,t,mW2,mZ2);

    Val +=  ( 8.*(t-mW2)*(u*t-2.*mW2*mZ2)/t/t
	    )*TildeI3W(s,t,mW2);
    swap(mW2,mZ2);
    del *= -1;

    return Val;
}

/***************************************************************************/
double  Iud1(Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2, double beta) { 
    Energy2 sig(mZ2+mW2);
    Energy2 del(mZ2-mW2);
    double Val(0.);

    Val +=    2.*(4.*t*t+t*(9.*s-4.*sig)-18.*s*sig)/t/u     
            + 8.*(t*u+2.*s*sig)/mW2/mZ2
            + 4.*s*s*(2.*t-sig)/u/(mW2*mZ2-t*(u+s))
	    - 2.*sqr(t-u)/u/s/sqr(beta);

    Val +=  ( 2.*(8.*t*t-4.*t*(s+3.*sig)-(s-sig)*(3.*s+4.*sig))/u/s/sqr(beta)
	    + 6.*(t+u)*sqr(t-u)/u/s/s/sqr(sqr(beta))
	    - 12.*s*(t-sig)/t/u
	    )*log(-t/s);
    
    Val +=  ( (2./u/s/s/sqr(beta))*( 4.*t*t*(-2.*s+del)
				   + 4.*t*(s*s+s*(mZ2+5.*mW2)-2.*sig*del)
				   + (s-sig)*(3.*s*s+8.*mW2*s-3.*sig*del)
                                   )
	    + (2.*t*(18.*s+3.*mW2+mZ2)-24.*s*sig)/t/u
	    - 8.*s*(2.*t*t-t*(3.*s+4.*mZ2+2.*mW2)+2.*mZ2*(s+sig))
                  /u/(mW2*mZ2-t*(u+s))
	    - 8.*s*s*t*(2.*t-sig)*(t-mZ2)/u/sqr(mW2*mZ2-t*(u+s))
	    + 6.*(s-del)*(s-sig)*sqr(t-u)/u/s/s/s/sqr(sqr(beta))
	    )*log(-t/mW2);

    Val +=  ( -2.*(2.*t*t*(2.*sig-3.*s)+6.*sig*t*(s-sig)+sqr(s-sig)*(s+2.*sig))
	         /u/s/sqr(beta)
	      +3.*s*(4.*t-4.*sig-s)/u
	      -3.*sqr(s-sig)*sqr(t-u)/u/s/s/sqr(sqr(beta))
	    )*TildeI3WZ(s,mW2,mZ2,beta);
    
    Val +=  ( 4.*(u+4.*s)/3./u - 4.*(u*t+2.*s*sig)/3./mW2/mZ2 
	    )*pi*pi;

    Val += -( 16.*s*(t-sig)*(t-mW2)/t/u
	    )*TildeI3W(s,t,mW2);

    Val +=  ( 8.*s*s*(t-sig)/u
	    )*TildeI4t(s,t,mW2,mZ2);

    swap(t,u);
    Val +=    2.*(4.*t*t+t*(9.*s-4.*sig)-18.*s*sig)/t/u     
            + 8.*(t*u+2.*s*sig)/mW2/mZ2
            + 4.*s*s*(2.*t-sig)/u/(mW2*mZ2-t*(u+s))
	    - 2.*sqr(t-u)/u/s/sqr(beta);

    Val +=  ( 2.*(8.*t*t-4.*t*(s+3.*sig)-(s-sig)*(3.*s+4.*sig))/u/s/sqr(beta)
	    + 6.*(t+u)*sqr(t-u)/u/s/s/sqr(sqr(beta))
	    - 12.*s*(t-sig)/t/u
	    )*log(-t/s);
    
    Val +=  ( (2./u/s/s/sqr(beta))*( 4.*t*t*(-2.*s+del)
				   + 4.*t*(s*s+s*(mZ2+5.*mW2)-2.*sig*del)
				   + (s-sig)*(3.*s*s+8.*mW2*s-3.*sig*del)
                                   )
	    + (2.*t*(18.*s+3.*mW2+mZ2)-24.*s*sig)/t/u
	    - 8.*s*(2.*t*t-t*(3.*s+4.*mZ2+2.*mW2)+2.*mZ2*(s+sig))
                  /u/(mW2*mZ2-t*(u+s))
	    - 8.*s*s*t*(2.*t-sig)*(t-mZ2)/u/sqr(mW2*mZ2-t*(u+s))
	    + 6.*(s-del)*(s-sig)*sqr(t-u)/u/s/s/s/sqr(sqr(beta))
	    )*log(-t/mW2);

    Val +=  ( -2.*(2.*t*t*(2.*sig-3.*s)+6.*sig*t*(s-sig)+sqr(s-sig)*(s+2.*sig))
	         /u/s/sqr(beta)
	      +3.*s*(4.*t-4.*sig-s)/u
	      -3.*sqr(s-sig)*sqr(t-u)/u/s/s/sqr(sqr(beta))
	    )*TildeI3WZ(s,mW2,mZ2,beta);
    
    Val +=  ( 4.*(u+4.*s)/3./u - 4.*(u*t+2.*s*sig)/3./mW2/mZ2 
	    )*pi*pi;

    Val += -( 16.*s*(t-sig)*(t-mW2)/t/u
	    )*TildeI3W(s,t,mW2);

    Val +=  ( 8.*s*s*(t-sig)/u
	    )*TildeI4t(s,t,mW2,mZ2);
    swap(t,u);

    swap(mW2,mZ2);
    del *= -1;
    Val +=    2.*(4.*t*t+t*(9.*s-4.*sig)-18.*s*sig)/t/u     
            + 8.*(t*u+2.*s*sig)/mW2/mZ2
            + 4.*s*s*(2.*t-sig)/u/(mW2*mZ2-t*(u+s))
	    - 2.*sqr(t-u)/u/s/sqr(beta);

    Val +=  ( 2.*(8.*t*t-4.*t*(s+3.*sig)-(s-sig)*(3.*s+4.*sig))/u/s/sqr(beta)
	    + 6.*(t+u)*sqr(t-u)/u/s/s/sqr(sqr(beta))
	    - 12.*s*(t-sig)/t/u
	    )*log(-t/s);
    
    Val +=  ( (2./u/s/s/sqr(beta))*( 4.*t*t*(-2.*s+del)
				   + 4.*t*(s*s+s*(mZ2+5.*mW2)-2.*sig*del)
				   + (s-sig)*(3.*s*s+8.*mW2*s-3.*sig*del)
                                   )
	    + (2.*t*(18.*s+3.*mW2+mZ2)-24.*s*sig)/t/u
	    - 8.*s*(2.*t*t-t*(3.*s+4.*mZ2+2.*mW2)+2.*mZ2*(s+sig))
                  /u/(mW2*mZ2-t*(u+s))
	    - 8.*s*s*t*(2.*t-sig)*(t-mZ2)/u/sqr(mW2*mZ2-t*(u+s))
	    + 6.*(s-del)*(s-sig)*sqr(t-u)/u/s/s/s/sqr(sqr(beta))
	    )*log(-t/mW2);

    Val +=  ( -2.*(2.*t*t*(2.*sig-3.*s)+6.*sig*t*(s-sig)+sqr(s-sig)*(s+2.*sig))
	         /u/s/sqr(beta)
	      +3.*s*(4.*t-4.*sig-s)/u
	      -3.*sqr(s-sig)*sqr(t-u)/u/s/s/sqr(sqr(beta))
	    )*TildeI3WZ(s,mW2,mZ2,beta);
    
    Val +=  ( 4.*(u+4.*s)/3./u - 4.*(u*t+2.*s*sig)/3./mW2/mZ2 
	    )*pi*pi;

    Val += -( 16.*s*(t-sig)*(t-mW2)/t/u
	    )*TildeI3W(s,t,mW2);

    Val +=  ( 8.*s*s*(t-sig)/u
	    )*TildeI4t(s,t,mW2,mZ2);
    swap(mW2,mZ2);
    del *= -1;

    swap(t,u);
    swap(mW2,mZ2);
    del *= -1;
    Val +=    2.*(4.*t*t+t*(9.*s-4.*sig)-18.*s*sig)/t/u     
            + 8.*(t*u+2.*s*sig)/mW2/mZ2
            + 4.*s*s*(2.*t-sig)/u/(mW2*mZ2-t*(u+s))
	    - 2.*sqr(t-u)/u/s/sqr(beta);

    Val +=  ( 2.*(8.*t*t-4.*t*(s+3.*sig)-(s-sig)*(3.*s+4.*sig))/u/s/sqr(beta)
	    + 6.*(t+u)*sqr(t-u)/u/s/s/sqr(sqr(beta))
	    - 12.*s*(t-sig)/t/u
	    )*log(-t/s);
    
    Val +=  ( (2./u/s/s/sqr(beta))*( 4.*t*t*(-2.*s+del)
				   + 4.*t*(s*s+s*(mZ2+5.*mW2)-2.*sig*del)
				   + (s-sig)*(3.*s*s+8.*mW2*s-3.*sig*del)
                                   )
	    + (2.*t*(18.*s+3.*mW2+mZ2)-24.*s*sig)/t/u
	    - 8.*s*(2.*t*t-t*(3.*s+4.*mZ2+2.*mW2)+2.*mZ2*(s+sig))
                  /u/(mW2*mZ2-t*(u+s))
	    - 8.*s*s*t*(2.*t-sig)*(t-mZ2)/u/sqr(mW2*mZ2-t*(u+s))
	    + 6.*(s-del)*(s-sig)*sqr(t-u)/u/s/s/s/sqr(sqr(beta))
	    )*log(-t/mW2);

    Val +=  ( -2.*(2.*t*t*(2.*sig-3.*s)+6.*sig*t*(s-sig)+sqr(s-sig)*(s+2.*sig))
	         /u/s/sqr(beta)
	      +3.*s*(4.*t-4.*sig-s)/u
	      -3.*sqr(s-sig)*sqr(t-u)/u/s/s/sqr(sqr(beta))
	    )*TildeI3WZ(s,mW2,mZ2,beta);
    
    Val +=  ( 4.*(u+4.*s)/3./u - 4.*(u*t+2.*s*sig)/3./mW2/mZ2 
	    )*pi*pi;

    Val += -( 16.*s*(t-sig)*(t-mW2)/t/u
	    )*TildeI3W(s,t,mW2);

    Val +=  ( 8.*s*s*(t-sig)/u
	    )*TildeI4t(s,t,mW2,mZ2);
    swap(t,u);
    swap(mW2,mZ2);
    del *= -1;

    return Val;
}

/***************************************************************************/
double  Iuu1(Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2, double beta) {
    double Val(Idd1(s,u,t,mW2,mZ2,beta));
    return Val;
}

/***************************************************************************/
Energy2 Fd1 (Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2, double beta) {
    Energy2 sig(mZ2+mW2);
    Energy2 del(mZ2-mW2);
    Energy2 Val(0.*GeV2);

    Val +=    4.*(17.*t*t+t*(11.*s-13.*sig)+17.*(s*sig+mW2*mZ2))/t     
	    + 16.*(s-sig)*(t*u+2.*s*sig)/mW2/mZ2
	    + 4*s*s*(2.*t-sig)/(t*(u+s)-mW2*mZ2);

    Val +=  ( 8.*(t-u)/sqr(beta)
	    - 4.*(3.*t*t-t*(s+3.*sig)+3.*(s*sig+mW2*mZ2))/t
	    )*log(-t/s);

    Val +=  ( 8.*(t*t-t*(2.*s+3.*mW2+mZ2)+3.*(s*sig+mW2*mZ2))/t
	    + 8.*s*(t*(3.*s+2.*sig)-2.*mZ2*(s+sig))/(t*(u+s)-mW2*mZ2)
	    + 8.*s*s*t*(2.*t-sig)*(t-mZ2)/sqr(t*(u+s)-mW2*mZ2)
	    - 8.*(s-del)*(t-u)/s/sqr(beta)  
	    )*log(-t/mW2);

    Val +=  ( 4.*(s-sig)*(t-u)/sqr(beta)
	    + 4.*(sig-3.*s)*t
	    + 4.*(4.*s*sig-mZ2*mZ2-mW2*mW2)
   	    )*TildeI3WZ(s,mW2,mZ2,beta);

    Val += -( 8.*(3.*t*t+2.*t*(2.*s-sig)+2.*(s*sig+mW2*mZ2))/3./t
	    + 8.*(s-sig)*(t*u+2.*s*sig)/3./mW2/mZ2
	    )*pi*pi;

    Val +=  ( 4.*(s*t*t-s*(s+sig)*t+2.*s*(s*sig+mW2*mZ2))
   	    )*TildeI4t(s,t,mW2,mZ2);

    Val += -( 8.*(t-mW2)*(t*t-t*(s+sig)+2.*(s*sig+mW2*mZ2))/t
   	    )*TildeI3W(s,t,mW2);


    swap(mW2,mZ2);
    del *= -1;
    Val +=    4.*(17.*t*t+t*(11.*s-13.*sig)+17.*(s*sig+mW2*mZ2))/t     
	    + 16.*(s-sig)*(t*u+2.*s*sig)/mW2/mZ2
	    + 4*s*s*(2.*t-sig)/(t*(u+s)-mW2*mZ2);

    Val +=  ( 8.*(t-u)/sqr(beta)
	    - 4.*(3.*t*t-t*(s+3.*sig)+3.*(s*sig+mW2*mZ2))/t
	    )*log(-t/s);

    Val +=  ( 8.*(t*t-t*(2.*s+3.*mW2+mZ2)+3.*(s*sig+mW2*mZ2))/t
	    + 8.*s*(t*(3.*s+2.*sig)-2.*mZ2*(s+sig))/(t*(u+s)-mW2*mZ2)
	    + 8.*s*s*t*(2.*t-sig)*(t-mZ2)/sqr(t*(u+s)-mW2*mZ2)
	    - 8.*(s-del)*(t-u)/s/sqr(beta)  
	    )*log(-t/mW2);

    Val +=  ( 4.*(s-sig)*(t-u)/sqr(beta)
	    + 4.*(sig-3.*s)*t
	    + 4.*(4.*s*sig-mZ2*mZ2-mW2*mW2)
   	    )*TildeI3WZ(s,mW2,mZ2,beta);

    Val += -( 8.*(3.*t*t+2.*t*(2.*s-sig)+2.*(s*sig+mW2*mZ2))/3./t
	    + 8.*(s-sig)*(t*u+2.*s*sig)/3./mW2/mZ2
	    )*pi*pi;

    Val +=  ( 4.*(s*t*t-s*(s+sig)*t+2.*s*(s*sig+mW2*mZ2))
   	    )*TildeI4t(s,t,mW2,mZ2);

    Val += -( 8.*(t-mW2)*(t*t-t*(s+sig)+2.*(s*sig+mW2*mZ2))/t
   	    )*TildeI3W(s,t,mW2);
    swap(mW2,mZ2);
    del *= -1;

    return Val;
}

/***************************************************************************/
Energy2 Fu1 (Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2, double beta) { 
    Energy2 Val(Fd1(s,u,t,mW2,mZ2,beta));
    return Val;
}

/***************************************************************************/
Energy4 H1  (Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2) { 
    Energy2 sig(mZ2+mW2);
    Energy2 del(mZ2-mW2);
    Energy4 Val(0.*GeV2*GeV2);
    Val  =   8.*t*t+8.*t*(s-sig)+s*s+6.*s*sig+mZ2*mZ2+10.*mW2*mZ2+mW2*mW2
	   - sqr(s-sig)*(t*u+2.*s*sig)/mW2/mZ2;
    Val *= ( 16.-8.*pi*pi/3.);
    return Val;
}

Energy2 t_u_Rdd(Energy2 s  , Energy2 tk , Energy2 uk , Energy2 q1 , Energy2 q2,
		Energy2 mW2, Energy2 mZ2);
Energy2 t_u_Rud(Energy2 s  , Energy2 tk , Energy2 uk , Energy2 q1 , Energy2 q2,
		Energy2 q1h, Energy2 q2h, Energy2 mW2, Energy2 mZ2);
Energy2 t_u_Ruu(Energy2 s  , Energy2 tk , Energy2 uk, Energy2 q1h, Energy2 q2h,
		Energy2 mW2, Energy2 mZ2);
Energy4 t_u_RZds(Energy2 s ,Energy2 tk , Energy2 uk , Energy2 q1, Energy2 q2,
		 Energy2 s2,Energy2 mW2, Energy2 mZ2);
Energy4 t_u_RZda(Energy2 s , Energy2 tk , Energy2 uk , Energy2 q1, Energy2 q2,
		 Energy2 s2, Energy2 mW2, Energy2 mZ2);
Energy4 t_u_RZd(Energy2 s  , Energy2 tk , Energy2 uk , Energy2 q1 , Energy2 q2 ,
		Energy2 s2 , Energy2 mW2, Energy2 mZ2);
Energy4 t_u_RZu(Energy2 s  , Energy2 tk , Energy2 uk , Energy2 q1h, Energy2 q2h,
		Energy2 s2 , Energy2 mW2, Energy2 mZ2);
Energy6 t_u_RZs(Energy2 s , Energy2 tk , Energy2 uk , Energy2 q1, Energy2 q2,
		Energy2 s2, Energy2 mW2, Energy2 mZ2);
Energy6 t_u_RZa(Energy2 s , Energy2 tk , Energy2 uk , Energy2 q1, Energy2 q2,
		Energy2 s2, Energy2 mW2, Energy2 mZ2);
Energy6 t_u_RZ(Energy2 s  , Energy2 tk , Energy2 uk , Energy2 q1, Energy2 q2,
	       Energy2 s2 , Energy2 mW2, Energy2 mZ2);

/***************************************************************************/
// t_u_M_R_qqb is the real emission q + qb -> n + g matrix element 
// exactly as defined in Eqs. C.1 of NPB 383(1992)3-44, multiplied by
// tk * uk!
Energy2 MEPP2VVPowheg::t_u_M_R_qqb(real2to3Kinematics R) const {
  // First the Born variables:
  Energy2 s2(R.s2r());
  Energy2 mW2(R.k12r());
  Energy2 mZ2(R.k22r());
  // Then the rest:
  Energy2 s(R.sr());
  Energy2 tk(R.tkr());
  Energy2 uk(R.ukr());
  Energy2 q1(R.q1r());
  Energy2 q2(R.q2r());
  Energy2 q1h(R.q1hatr());
  Energy2 q2h(R.q2hatr());
  Energy2 w1(R.w1r());
  Energy2 w2(R.w2r());

  return -2.*pi*alphaS_*Fij2_*CF_/NC_
       * (    gdL_*gdL_*t_u_Rdd(s,tk,uk,q1,q2,mW2,mZ2)
	 + 2.*gdL_*guL_*t_u_Rud(s,tk,uk,q1,q2,q1h,q2h,mW2,mZ2)
	 +    guL_*guL_*t_u_Ruu(s,tk,uk,q1h,q2h,mW2,mZ2)
	 - 2.*eZ_/(s2-mW2) * ( gdL_
			      *t_u_RZd(s,tk,uk,q1 ,q2 ,s2,mW2,mZ2)
	                     - guL_
			      *t_u_RZu(s,tk,uk,q1h,q2h,s2,mW2,mZ2)
	                     )
         + sqr(eZ_/(s2-mW2)) *t_u_RZ(s,tk,uk,q1,q2,s2,mW2,mZ2)
	 );
}

Energy2 t_u_Rdd(Energy2 s ,Energy2 tk ,Energy2 uk ,Energy2 q1,Energy2 q2,
            Energy2 mW2, Energy2 mZ2) {
  Energy2 Val(0.*GeV2);

  Val +=   4.*(q2*(uk+2.*s+q2)+q1*(s+q1))/mW2/mZ2*uk
         + 16.*(uk+s)/q2*uk
         - 4.*(2.*uk+4.*s+q2)/mW2*uk
         - 4.*(2.*uk+5.*s+q2+2.*q1-mW2)/mZ2*uk
         + 4.*q1*s*(s+q1)/mW2/mZ2
         + 16.*s*(s+q2-mZ2-mW2)/q1
         - 4.*s*(4.*s+q2+q1)/mW2
         + 16.*mW2*mZ2*s/q1/q2
         + 4.*s
         + 16.*mZ2*(tk-2.*mW2)/q1/q2/q2*tk*uk
         + 16.*(2.*mZ2+mW2-tk)/q1/q2*tk*uk
         + 16.*mW2*(s-mZ2-mW2)/q1/q2*uk
         + 16.*mZ2*(q1-2.*mW2)/q2/q2*uk
         + 32.*mW2*mW2*mZ2/q1/q2/q2*uk
         + 16.*mW2/q1*uk
         + 4.*uk
         + 8./q2*tk*uk
         + 4.*q1/mW2/mZ2*tk*uk
         - 24./q1*tk*uk
         - 4./mW2*tk*uk;

  swap(mW2,mZ2);
  swap(q1,q2);
  swap(tk,uk);
  Val +=   4.*(q2*(uk+2.*s+q2)+q1*(s+q1))/mW2/mZ2*uk
         + 16.*(uk+s)/q2*uk
         - 4.*(2.*uk+4.*s+q2)/mW2*uk
         - 4.*(2.*uk+5.*s+q2+2.*q1-mW2)/mZ2*uk
         + 4.*q1*s*(s+q1)/mW2/mZ2
         + 16.*s*(s+q2-mZ2-mW2)/q1
         - 4.*s*(4.*s+q2+q1)/mW2
         + 16.*mW2*mZ2*s/q1/q2
         + 4.*s
         + 16.*mZ2*(tk-2.*mW2)/q1/q2/q2*tk*uk
         + 16.*(2.*mZ2+mW2-tk)/q1/q2*tk*uk
         + 16.*mW2*(s-mZ2-mW2)/q1/q2*uk
         + 16.*mZ2*(q1-2.*mW2)/q2/q2*uk
         + 32.*mW2*mW2*mZ2/q1/q2/q2*uk
         + 16.*mW2/q1*uk
         + 4.*uk
         + 8./q2*tk*uk
         + 4.*q1/mW2/mZ2*tk*uk
         - 24./q1*tk*uk
         - 4./mW2*tk*uk;
  swap(mW2,mZ2);
  swap(q1,q2);
  swap(tk,uk);

  return Val;
}
Energy2 t_u_Rud(Energy2 s ,Energy2 tk ,Energy2 uk ,Energy2 q1,Energy2 q2,
	    Energy2 q1h,Energy2 q2h,Energy2 mW2, Energy2 mZ2) {
  Energy2 Val(0.*GeV2); 
  
  Val +=   (uk*s*(uk+3.*s+q1h)+s*s*(s+mZ2)-(s+uk)*(2.*mZ2*s+3.*mW2*s+mW2*q1h)
           ) * 8./q1/q2h/q2*uk  
         - (uk*(uk+3.*s+q1h-mW2)-(q2+s)*(q2-s)+s*(q2-mW2)+q1h*(q2-mW2)+mW2*q2
	   ) * 4.*s/mZ2/q1/q2h*uk
         - 4.*((s+uk+q2h-2.*mZ2)*(s+q1h-mZ2)-mZ2*q1)/mW2/q2*uk
         + 4.*(2.*s*uk+2.*mW2*uk+5.*s*s+2.*q1h*s-2.*mZ2*s)/q1/q2h*uk
         + 4.*(2.*s*uk-s*s-2.*q1h*s+2.*mW2*s+2.*mW2*q1h)/q1/q2h/q2*tk*uk
         + ((2.*uk+s)*(s+q1h)+s*(q2+q2h)+2.*q2*(s+q2h)-q1*s+q1*q2+q1h*q2h
	   ) /mW2/mZ2*uk
         + 8.*s*(uk-q1h+mZ2)/q1/q2*uk
         + 4.*s*(-uk+s-q2+q1+q1h)/mZ2/q2h*uk
         + 4.*s*(-uk-q2+q1h)/mZ2/q1*uk
         + 8.*(mZ2*uk-s*s+mW2*s-2.*mZ2*q1-2.*mZ2*q1h)/q2h/q2*uk
         + 2.*(-uk-9.*s-4.*q2-5.*q2h-3.*q1-4.*q1h+8.*mZ2)/mW2*uk
         + 2.*(-4.*uk+3.*s+5.*q1+4.*q1h)/q2h*uk
         + 2.*(s*tk+q2*tk+s*s-q2*q2+q1h*q2)/mW2/mZ2*tk
         - 8.*s*(tk+s+q1h)/mW2/q2*tk
         + 2.*(-tk+3.*s+q2-q1h)/mW2*tk
         - 8.*s*s*s/q1h/q2
         - 2.*s*q2*(s+q2)/mW2/mZ2
         + 2.*s*(2.*s+q2)/mZ2
         + 2.*s*(2.*s+q2)/mW2
         - 16.*s*s/q1h
         - 2.*s
         - 16.*s*s/q1h/q2*tk
         - 8.*s/q2*tk
         - 16.*s/q1h*tk
         + 6.*s/mZ2*tk
         + 4.*s/q1*uk
         + 4.*s/mZ2*uk
         + 12.*uk
         + 4.*s*(tk+q1h-mW2)/mZ2/q1/q2h*tk*uk
         + 2.*(s+4.*q1+5.*q1h-4.*mZ2)/q2*uk
         - 4.*s*s*s/q1h/q1/q2h/q2*tk*uk
         - 4.*s*s/q1h/q2h/q2*tk*uk
         - 4.*s*s/q1h/q1/q2*tk*uk
         + 8.*s*s/mW2/q1h/q2*tk*uk
         - 4.*s*s/q1h/q1/q2h*tk*uk
         + 4.*(s+mZ2)/mW2/q2*tk*uk
         - 4.*s/q1h/q2h*tk*uk
         - 4.*s/q1h/q1*tk*uk
         + 12.*s/mW2/q1h*tk*uk
         - (s+4.*q2)/mW2/mZ2*tk*uk
         - 4.*(s+2.*mZ2)/q2h/q2*tk*uk
         - 4.*(3.*s+2.*q1h)/q1/q2*tk*uk
         - 8.*mW2/q1/q2h*tk*uk
         + 8./q2h*tk*uk
         + 8./q1*tk*uk;

  swap(mW2,mZ2);
  swap(q1,q2);
  swap(tk,uk);
  swap(q1h,q2h); // Note this swap is done in accordance with MC@NLO.
                 // It is not in NPB 383(1992)3-44 Eq.C.4!
  Val +=   (uk*s*(uk+3.*s+q1h)+s*s*(s+mZ2)-(s+uk)*(2.*mZ2*s+3.*mW2*s+mW2*q1h)
           ) * 8./q1/q2h/q2*uk  
         - (uk*(uk+3.*s+q1h-mW2)-(q2+s)*(q2-s)+s*(q2-mW2)+q1h*(q2-mW2)+mW2*q2
	   ) * 4.*s/mZ2/q1/q2h*uk
         - 4.*((s+uk+q2h-2.*mZ2)*(s+q1h-mZ2)-mZ2*q1)/mW2/q2*uk
         + 4.*(2.*s*uk+2.*mW2*uk+5.*s*s+2.*q1h*s-2.*mZ2*s)/q1/q2h*uk
         + 4.*(2.*s*uk-s*s-2.*q1h*s+2.*mW2*s+2.*mW2*q1h)/q1/q2h/q2*tk*uk
         + ((2.*uk+s)*(s+q1h)+s*(q2+q2h)+2.*q2*(s+q2h)-q1*s+q1*q2+q1h*q2h
	   ) /mW2/mZ2*uk
         + 8.*s*(uk-q1h+mZ2)/q1/q2*uk
         + 4.*s*(-uk+s-q2+q1+q1h)/mZ2/q2h*uk
         + 4.*s*(-uk-q2+q1h)/mZ2/q1*uk
         + 8.*(mZ2*uk-s*s+mW2*s-2.*mZ2*q1-2.*mZ2*q1h)/q2h/q2*uk
         + 2.*(-uk-9.*s-4.*q2-5.*q2h-3.*q1-4.*q1h+8.*mZ2)/mW2*uk
         + 2.*(-4.*uk+3.*s+5.*q1+4.*q1h)/q2h*uk
         + 2.*(s*tk+q2*tk+s*s-q2*q2+q1h*q2)/mW2/mZ2*tk
         - 8.*s*(tk+s+q1h)/mW2/q2*tk
         + 2.*(-tk+3.*s+q2-q1h)/mW2*tk
         - 8.*s*s*s/q1h/q2
         - 2.*s*q2*(s+q2)/mW2/mZ2
         + 2.*s*(2.*s+q2)/mZ2
         + 2.*s*(2.*s+q2)/mW2
         - 16.*s*s/q1h
         - 2.*s
         - 16.*s*s/q1h/q2*tk
         - 8.*s/q2*tk
         - 16.*s/q1h*tk
         + 6.*s/mZ2*tk
         + 4.*s/q1*uk
         + 4.*s/mZ2*uk
         + 12.*uk
         + 4.*s*(tk+q1h-mW2)/mZ2/q1/q2h*tk*uk
         + 2.*(s+4.*q1+5.*q1h-4.*mZ2)/q2*uk
         - 4.*s*s*s/q1h/q1/q2h/q2*tk*uk
         - 4.*s*s/q1h/q2h/q2*tk*uk
         - 4.*s*s/q1h/q1/q2*tk*uk
         + 8.*s*s/mW2/q1h/q2*tk*uk
         - 4.*s*s/q1h/q1/q2h*tk*uk
         + 4.*(s+mZ2)/mW2/q2*tk*uk
         - 4.*s/q1h/q2h*tk*uk
         - 4.*s/q1h/q1*tk*uk
         + 12.*s/mW2/q1h*tk*uk
         - (s+4.*q2)/mW2/mZ2*tk*uk
         - 4.*(s+2.*mZ2)/q2h/q2*tk*uk
         - 4.*(3.*s+2.*q1h)/q1/q2*tk*uk
         - 8.*mW2/q1/q2h*tk*uk
         + 8./q2h*tk*uk
         + 8./q1*tk*uk;
  swap(mW2,mZ2);
  swap(q1,q2);
  swap(tk,uk);
  swap(q1h,q2h); // Note this swap is done in accordance with MC@NLO.
                 // It is not in NPB 383(1992)3-44 Eq.C.4!

  swap(tk,uk);
  swap(q1,q2h);
  swap(q2,q1h);
  Val +=   (uk*s*(uk+3.*s+q1h)+s*s*(s+mZ2)-(s+uk)*(2.*mZ2*s+3.*mW2*s+mW2*q1h)
           ) * 8./q1/q2h/q2*uk  
         - (uk*(uk+3.*s+q1h-mW2)-(q2+s)*(q2-s)+s*(q2-mW2)+q1h*(q2-mW2)+mW2*q2
	   ) * 4.*s/mZ2/q1/q2h*uk
         - 4.*((s+uk+q2h-2.*mZ2)*(s+q1h-mZ2)-mZ2*q1)/mW2/q2*uk
         + 4.*(2.*s*uk+2.*mW2*uk+5.*s*s+2.*q1h*s-2.*mZ2*s)/q1/q2h*uk
         + 4.*(2.*s*uk-s*s-2.*q1h*s+2.*mW2*s+2.*mW2*q1h)/q1/q2h/q2*tk*uk
         + ((2.*uk+s)*(s+q1h)+s*(q2+q2h)+2.*q2*(s+q2h)-q1*s+q1*q2+q1h*q2h
	   ) /mW2/mZ2*uk
         + 8.*s*(uk-q1h+mZ2)/q1/q2*uk
         + 4.*s*(-uk+s-q2+q1+q1h)/mZ2/q2h*uk
         + 4.*s*(-uk-q2+q1h)/mZ2/q1*uk
         + 8.*(mZ2*uk-s*s+mW2*s-2.*mZ2*q1-2.*mZ2*q1h)/q2h/q2*uk
         + 2.*(-uk-9.*s-4.*q2-5.*q2h-3.*q1-4.*q1h+8.*mZ2)/mW2*uk
         + 2.*(-4.*uk+3.*s+5.*q1+4.*q1h)/q2h*uk
         + 2.*(s*tk+q2*tk+s*s-q2*q2+q1h*q2)/mW2/mZ2*tk
         - 8.*s*(tk+s+q1h)/mW2/q2*tk
         + 2.*(-tk+3.*s+q2-q1h)/mW2*tk
         - 8.*s*s*s/q1h/q2
         - 2.*s*q2*(s+q2)/mW2/mZ2
         + 2.*s*(2.*s+q2)/mZ2
         + 2.*s*(2.*s+q2)/mW2
         - 16.*s*s/q1h
         - 2.*s
         - 16.*s*s/q1h/q2*tk
         - 8.*s/q2*tk
         - 16.*s/q1h*tk
         + 6.*s/mZ2*tk
         + 4.*s/q1*uk
         + 4.*s/mZ2*uk
         + 12.*uk
         + 4.*s*(tk+q1h-mW2)/mZ2/q1/q2h*tk*uk
         + 2.*(s+4.*q1+5.*q1h-4.*mZ2)/q2*uk
         - 4.*s*s*s/q1h/q1/q2h/q2*tk*uk
         - 4.*s*s/q1h/q2h/q2*tk*uk
         - 4.*s*s/q1h/q1/q2*tk*uk
         + 8.*s*s/mW2/q1h/q2*tk*uk
         - 4.*s*s/q1h/q1/q2h*tk*uk
         + 4.*(s+mZ2)/mW2/q2*tk*uk
         - 4.*s/q1h/q2h*tk*uk
         - 4.*s/q1h/q1*tk*uk
         + 12.*s/mW2/q1h*tk*uk
         - (s+4.*q2)/mW2/mZ2*tk*uk
         - 4.*(s+2.*mZ2)/q2h/q2*tk*uk
         - 4.*(3.*s+2.*q1h)/q1/q2*tk*uk
         - 8.*mW2/q1/q2h*tk*uk
         + 8./q2h*tk*uk
         + 8./q1*tk*uk;
  swap(tk,uk);
  swap(q1,q2h);
  swap(q2,q1h);

  swap(mW2,mZ2);
  swap(q1,q1h);
  swap(q2,q2h);
  Val +=   (uk*s*(uk+3.*s+q1h)+s*s*(s+mZ2)-(s+uk)*(2.*mZ2*s+3.*mW2*s+mW2*q1h)
           ) * 8./q1/q2h/q2*uk  
         - (uk*(uk+3.*s+q1h-mW2)-(q2+s)*(q2-s)+s*(q2-mW2)+q1h*(q2-mW2)+mW2*q2
	   ) * 4.*s/mZ2/q1/q2h*uk
         - 4.*((s+uk+q2h-2.*mZ2)*(s+q1h-mZ2)-mZ2*q1)/mW2/q2*uk
         + 4.*(2.*s*uk+2.*mW2*uk+5.*s*s+2.*q1h*s-2.*mZ2*s)/q1/q2h*uk
         + 4.*(2.*s*uk-s*s-2.*q1h*s+2.*mW2*s+2.*mW2*q1h)/q1/q2h/q2*tk*uk
         + ((2.*uk+s)*(s+q1h)+s*(q2+q2h)+2.*q2*(s+q2h)-q1*s+q1*q2+q1h*q2h
	   ) /mW2/mZ2*uk
         + 8.*s*(uk-q1h+mZ2)/q1/q2*uk
         + 4.*s*(-uk+s-q2+q1+q1h)/mZ2/q2h*uk
         + 4.*s*(-uk-q2+q1h)/mZ2/q1*uk
         + 8.*(mZ2*uk-s*s+mW2*s-2.*mZ2*q1-2.*mZ2*q1h)/q2h/q2*uk
         + 2.*(-uk-9.*s-4.*q2-5.*q2h-3.*q1-4.*q1h+8.*mZ2)/mW2*uk
         + 2.*(-4.*uk+3.*s+5.*q1+4.*q1h)/q2h*uk
         + 2.*(s*tk+q2*tk+s*s-q2*q2+q1h*q2)/mW2/mZ2*tk
         - 8.*s*(tk+s+q1h)/mW2/q2*tk
         + 2.*(-tk+3.*s+q2-q1h)/mW2*tk
         - 8.*s*s*s/q1h/q2
         - 2.*s*q2*(s+q2)/mW2/mZ2
         + 2.*s*(2.*s+q2)/mZ2
         + 2.*s*(2.*s+q2)/mW2
         - 16.*s*s/q1h
         - 2.*s
         - 16.*s*s/q1h/q2*tk
         - 8.*s/q2*tk
         - 16.*s/q1h*tk
         + 6.*s/mZ2*tk
         + 4.*s/q1*uk
         + 4.*s/mZ2*uk
         + 12.*uk
         + 4.*s*(tk+q1h-mW2)/mZ2/q1/q2h*tk*uk
         + 2.*(s+4.*q1+5.*q1h-4.*mZ2)/q2*uk
         - 4.*s*s*s/q1h/q1/q2h/q2*tk*uk
         - 4.*s*s/q1h/q2h/q2*tk*uk
         - 4.*s*s/q1h/q1/q2*tk*uk
         + 8.*s*s/mW2/q1h/q2*tk*uk
         - 4.*s*s/q1h/q1/q2h*tk*uk
         + 4.*(s+mZ2)/mW2/q2*tk*uk
         - 4.*s/q1h/q2h*tk*uk
         - 4.*s/q1h/q1*tk*uk
         + 12.*s/mW2/q1h*tk*uk
         - (s+4.*q2)/mW2/mZ2*tk*uk
         - 4.*(s+2.*mZ2)/q2h/q2*tk*uk
         - 4.*(3.*s+2.*q1h)/q1/q2*tk*uk
         - 8.*mW2/q1/q2h*tk*uk
         + 8./q2h*tk*uk
         + 8./q1*tk*uk;
  swap(mW2,mZ2);
  swap(q1,q1h);
  swap(q2,q2h);

  return Val;
 }

Energy2 t_u_Ruu(Energy2 s ,Energy2 tk ,Energy2 uk ,Energy2 q1h,Energy2 q2h,
		Energy2 mW2, Energy2 mZ2) {
  return t_u_Rdd(s,tk,uk,q1h,q2h,mZ2,mW2);
}

Energy4 t_u_RZds(Energy2 s ,Energy2 tk ,Energy2 uk ,Energy2 q1,Energy2 q2,
	    Energy2 s2, Energy2 mW2, Energy2 mZ2) {
  Energy4 Val(0.*GeV2*GeV2); 
  Energy2 sig(mZ2+mW2);

  Val +=   ( q1*q2*(5./2.*s*s+5.*s*tk+3.*tk*tk)+(tk*uk*uk+q1*q1*q2)*(tk+s)
	   + q1*(tk*tk*uk+s*uk*uk-s*s*tk+s*s*uk)+q1*q1*q1*(uk+s)-q1*q1*s*s2 
           ) * 8./q1/q2 
         - ( tk*tk*(4.*uk+s+q1-2.*q2)+tk*(sqr(q1+q2)-q1*s-3.*q2*s-2.*q1*q1)
	   - q1*s*(4.*s-2.*q1-q2)+tk*uk*(q1+3.*s)
	   ) * 4.*sig/q1/q2 
         - 4.*sig*sig*(s*(2.*s+q1)+tk*(uk+5./2.*tk+5.*s+q1+q2)
                      )/mW2/mZ2 
         + 2.*sig*s2*(4.*sqr(s+tk)+tk*(uk+s+4.*q1+2.*q2)+2.*q1*(2.*s+q1)
                     )/mW2/mZ2 
         + 4.*sig*sig*(s2+s-q1+q2)/q1/q2*tk 
         - 16.*mW2*mZ2*(tk*uk/2.+q2*tk-q1*s)/q1/q2 
         - 4.*s2*s2*q1*(tk+s+q1)/mW2/mZ2 
         + sig*sig*sig*(uk+tk)/mW2/mZ2 
         + 4.*mW2*mZ2*sig*(uk+tk)/q1/q2;

  swap(mW2,mZ2);
  swap(q1,q2);
  swap(tk,uk);
  Val +=   ( q1*q2*(5./2.*s*s+5.*s*tk+3.*tk*tk)+(tk*uk*uk+q1*q1*q2)*(tk+s)
	   + q1*(tk*tk*uk+s*uk*uk-s*s*tk+s*s*uk)+q1*q1*q1*(uk+s)-q1*q1*s*s2 
           ) * 8./q1/q2 
         - ( tk*tk*(4.*uk+s+q1-2.*q2)+tk*(sqr(q1+q2)-q1*s-3.*q2*s-2.*q1*q1)
	   - q1*s*(4.*s-2.*q1-q2)+tk*uk*(q1+3.*s)
	   ) * 4.*sig/q1/q2 
         - 4.*sig*sig*(s*(2.*s+q1)+tk*(uk+5./2.*tk+5.*s+q1+q2)
                      )/mW2/mZ2 
         + 2.*sig*s2*(4.*sqr(s+tk)+tk*(uk+s+4.*q1+2.*q2)+2.*q1*(2.*s+q1)
                     )/mW2/mZ2 
         + 4.*sig*sig*(s2+s-q1+q2)/q1/q2*tk 
         - 16.*mW2*mZ2*(tk*uk/2.+q2*tk-q1*s)/q1/q2 
         - 4.*s2*s2*q1*(tk+s+q1)/mW2/mZ2 
         + sig*sig*sig*(uk+tk)/mW2/mZ2 
         + 4.*mW2*mZ2*sig*(uk+tk)/q1/q2;
  swap(mW2,mZ2);
  swap(q1,q2);
  swap(tk,uk);

  return Val;
}
Energy4 t_u_RZda(Energy2 s ,Energy2 tk ,Energy2 uk ,Energy2 q1,Energy2 q2,
		 Energy2 s2, Energy2 mW2, Energy2 mZ2) {
  Energy4 Val(0.*GeV2*GeV2);

  Val +=   4.*mZ2*(2.*uk*uk-s*tk+q1*(uk-tk-s+q1+0.5*q2)+q2*(s-3.*q2)
                  ) /q1/q2*tk
         - 4.*mZ2*mZ2*(q1-tk-2.*s-q2)/q1/q2*tk
         - 2.*mZ2*(tk+2.*s+2.*q2)/mW2*tk
         - 2.*s2*(s+2.*q2)/mZ2*tk
         + 8.*mW2*mZ2*mZ2/q1/q2*tk
         + 2.*mZ2*mZ2/mW2*tk;

  swap(mW2,mZ2); // N.B. Here we subtract!
  Val -=   4.*mZ2*(2.*uk*uk-s*tk+q1*(uk-tk-s+q1+0.5*q2)+q2*(s-3.*q2)
                  ) /q1/q2*tk
         - 4.*mZ2*mZ2*(q1-tk-2.*s-q2)/q1/q2*tk
         - 2.*mZ2*(tk+2.*s+2.*q2)/mW2*tk
         - 2.*s2*(s+2.*q2)/mZ2*tk
         + 8.*mW2*mZ2*mZ2/q1/q2*tk
         + 2.*mZ2*mZ2/mW2*tk;
  swap(mW2,mZ2);

  swap(q1,q2); // N.B. Here we subtract!
  swap(tk,uk);
  Val -=   4.*mZ2*(2.*uk*uk-s*tk+q1*(uk-tk-s+q1+0.5*q2)+q2*(s-3.*q2)
                  ) /q1/q2*tk
         - 4.*mZ2*mZ2*(q1-tk-2.*s-q2)/q1/q2*tk
         - 2.*mZ2*(tk+2.*s+2.*q2)/mW2*tk
         - 2.*s2*(s+2.*q2)/mZ2*tk
         + 8.*mW2*mZ2*mZ2/q1/q2*tk
         + 2.*mZ2*mZ2/mW2*tk;
  swap(q1,q2);
  swap(tk,uk);

  swap(mW2,mZ2); // N.B. Here we add!
  swap(q1,q2); 
  swap(tk,uk);
  Val +=   4.*mZ2*(2.*uk*uk-s*tk+q1*(uk-tk-s+q1+0.5*q2)+q2*(s-3.*q2)
                  ) /q1/q2*tk
         - 4.*mZ2*mZ2*(q1-tk-2.*s-q2)/q1/q2*tk
         - 2.*mZ2*(tk+2.*s+2.*q2)/mW2*tk
         - 2.*s2*(s+2.*q2)/mZ2*tk
         + 8.*mW2*mZ2*mZ2/q1/q2*tk
         + 2.*mZ2*mZ2/mW2*tk;
  swap(mW2,mZ2);
  swap(q1,q2);
  swap(tk,uk);

  return Val;
}
Energy4 t_u_RZd(Energy2 s , Energy2 tk , Energy2 uk , Energy2 q1 , Energy2 q2 ,
		Energy2 s2, Energy2 mW2, Energy2 mZ2) {
  Energy4 Val(0.*GeV2*GeV2); 

  Val = t_u_RZds(s,tk,uk,q1,q2,s2,mW2,mZ2)
      + t_u_RZda(s,tk,uk,q1,q2,s2,mW2,mZ2);

  return Val;
}
Energy4 t_u_RZu(Energy2 s , Energy2 tk , Energy2 uk , Energy2 q1h, Energy2 q2h,
		Energy2 s2, Energy2 mW2, Energy2 mZ2) {
  Energy4 Val(0.*GeV2*GeV2);

  Val = t_u_RZd(s,tk,uk,q1h,q2h,s2,mZ2,mW2);

  return Val;
}
Energy6 t_u_RZs(Energy2 s,Energy2 tk,Energy2 uk,Energy2 q1,Energy2 q2,
		Energy2 s2,Energy2 mW2,Energy2 mZ2) {
  Energy6 Val(0.*GeV2*GeV2*GeV2); 
  Energy2 sig(mZ2+mW2);

  Val +=   2.*sig*sig*s2*( tk*(3.*uk+9.*tk+19.*s+6.*q1+4.*q2)+8.*s*s+6.*q1*s
	                 + 2.*q1*q1
                         )/mW2/mZ2
         - 2.*sig*sig*sig*(tk*(3.*uk+6.*tk+11.*s+2.*q1+2.*q2)+2.*s*(2.*s+q1))
           / mW2/mZ2
         - 2.*sig*s2*s2*(tk*(uk+4.*tk+9.*s+6.*q1+2.*q2)+4.*sqr(s+q1)-2.*q1*s)
           /mW2/mZ2
         - 16.*sig*(2.*tk*(uk/2.-tk-s+q1+q2)-s*(3.*s/2.-2.*q1))
         + 8.*s2*(s*(s/2.+tk)+4.*q1*(tk+s+q1))
         + 4.*s2*s2*s2*q1*(tk+s+q1)/mW2/mZ2
         + 8.*sig*sig*(2.*tk+s/2.)
         + 2.*sig*sig*sig*sig*tk/mW2/mZ2
         + 32.*mW2*mZ2*s;

  swap(mW2,mZ2);
  swap(q1,q2);
  swap(tk,uk);
  Val +=   2.*sig*sig*s2*( tk*(3.*uk+9.*tk+19.*s+6.*q1+4.*q2)+8.*s*s+6.*q1*s
	                 + 2.*q1*q1
                         )/mW2/mZ2
         - 2.*sig*sig*sig*(tk*(3.*uk+6.*tk+11.*s+2.*q1+2.*q2)+2.*s*(2.*s+q1))
           / mW2/mZ2
         - 2.*sig*s2*s2*(tk*(uk+4.*tk+9.*s+6.*q1+2.*q2)+4.*sqr(s+q1)-2.*q1*s)
           /mW2/mZ2
         - 16.*sig*(2.*tk*(uk/2.-tk-s+q1+q2)-s*(3.*s/2.-2.*q1))
         + 8.*s2*(s*(s/2.+tk)+4.*q1*(tk+s+q1))
         + 4.*s2*s2*s2*q1*(tk+s+q1)/mW2/mZ2
         + 8.*sig*sig*(2.*tk+s/2.)
         + 2.*sig*sig*sig*sig*tk/mW2/mZ2
         + 32.*mW2*mZ2*s;
  swap(mW2,mZ2);
  swap(q1,q2);
  swap(tk,uk);
  
  return Val;
}
Energy6 t_u_RZa(Energy2 s,Energy2 tk,Energy2 uk,Energy2 q1,Energy2 q2,
		Energy2 s2,Energy2 mW2,Energy2 mZ2) {
  Energy6 Val(0.*GeV2*GeV2*GeV2);
  Energy2 sig(mZ2+mW2);

  Val += - 2.*mZ2*(2.*tk+11.*s+18.*q2)*tk
         - 2.*mZ2*mZ2*(2.*tk+3.*s+2.*q2)/mW2*tk
         + 2.*mZ2*s2*(tk+3.*s+4.*q2)/mW2*tk
         - 2.*s2*s2*(s+2.*q2)/mW2*tk
         + 2.*mZ2*mZ2*mZ2/mW2*tk
         + 20.*mZ2*mZ2*tk;

  swap(mW2,mZ2);
  Val -= - 2.*mZ2*(2.*tk+11.*s+18.*q2)*tk
         - 2.*mZ2*mZ2*(2.*tk+3.*s+2.*q2)/mW2*tk
         + 2.*mZ2*s2*(tk+3.*s+4.*q2)/mW2*tk
         - 2.*s2*s2*(s+2.*q2)/mW2*tk
         + 2.*mZ2*mZ2*mZ2/mW2*tk
         + 20.*mZ2*mZ2*tk;
  swap(mW2,mZ2);

  swap(q1,q2);
  swap(tk,uk);
  Val -= - 2.*mZ2*(2.*tk+11.*s+18.*q2)*tk
         - 2.*mZ2*mZ2*(2.*tk+3.*s+2.*q2)/mW2*tk
         + 2.*mZ2*s2*(tk+3.*s+4.*q2)/mW2*tk
         - 2.*s2*s2*(s+2.*q2)/mW2*tk
         + 2.*mZ2*mZ2*mZ2/mW2*tk
         + 20.*mZ2*mZ2*tk;
  swap(q1,q2);
  swap(tk,uk);

  swap(mW2,mZ2);
  swap(q1,q2);
  swap(tk,uk);
  Val += - 2.*mZ2*(2.*tk+11.*s+18.*q2)*tk
         - 2.*mZ2*mZ2*(2.*tk+3.*s+2.*q2)/mW2*tk
         + 2.*mZ2*s2*(tk+3.*s+4.*q2)/mW2*tk
         - 2.*s2*s2*(s+2.*q2)/mW2*tk
         + 2.*mZ2*mZ2*mZ2/mW2*tk
         + 20.*mZ2*mZ2*tk;
  swap(mW2,mZ2);
  swap(q1,q2);
  swap(tk,uk);

  return Val;
}
Energy6 t_u_RZ(Energy2 s , Energy2 tk , Energy2 uk , Energy2 q1, Energy2 q2,
	       Energy2 s2, Energy2 mW2, Energy2 mZ2) {
  Energy6 Val(0.*GeV2*GeV2*GeV2); 

  Val = t_u_RZs(s,tk,uk,q1,q2,s2,mW2,mZ2)
      + t_u_RZa(s,tk,uk,q1,q2,s2,mW2,mZ2);

  return Val;
}
/***************************************************************************/
// t_u_M_R_qg is the real emission q + qb -> n + g matrix element 
// exactly as defined in Eqs. C.9 of NPB 383(1992)3-44, multiplied by
// tk * uk!
Energy2 MEPP2VVPowheg::t_u_M_R_qg(real2to3Kinematics R) const {
  // First the Born variables:
  Energy2 s2(R.s2r());
  Energy2 mW2(R.k12r());
  Energy2 mZ2(R.k22r());
  // Then the rest:
  Energy2 s(R.sr());
  Energy2 tk(R.tkr());
  Energy2 uk(R.ukr());
  Energy2 q1(R.q1r());
  Energy2 q2(R.q2r());
  Energy2 q1h(R.q1hatr());
  Energy2 q2h(R.q2hatr());
  Energy2 w1(R.w1r());
  Energy2 w2(R.w2r());

  Energy2 Val(0.*GeV2);

  swap(s,tk);
  swap(q2,w2);
  swap(q2h,w1);
  Val  =  -2.*pi*alphaS_*Fij2_*CF_/NC_
        * (    gdL_*gdL_*t_u_Rdd(s,tk,uk,q1,q2,mW2,mZ2)
	  + 2.*gdL_*guL_*t_u_Rud(s,tk,uk,q1,q2,q1h,q2h,mW2,mZ2)
	  +    guL_*guL_*t_u_Ruu(s,tk,uk,q1h,q2h,mW2,mZ2)
	  - 2.*eZ_/(s2-mW2) * ( gdL_
		 	       *t_u_RZd(s,tk,uk,q1 ,q2 ,s2,mW2,mZ2)
	                      - guL_
			       *t_u_RZu(s,tk,uk,q1h,q2h,s2,mW2,mZ2)
	                      )
          + sqr(eZ_/(s2-mW2)) *t_u_RZ(s,tk,uk,q1,q2,s2,mW2,mZ2)
	  );
  swap(s,tk);
  swap(q2,w2);
  swap(q2h,w1);

  Val *= -tk/s * TR_/CF_; 

  return Val;
}
/***************************************************************************/
// t_u_M_R_gqb is the real emission g + qb -> n + q matrix element 
// exactly as defined in Eqs. C.9 of NPB 383(1992)3-44, multiplied by
// tk * uk!
Energy2 MEPP2VVPowheg::t_u_M_R_gqb(real2to3Kinematics R) const {
  // First the Born variables:
  Energy2 s2(R.s2r());
  Energy2 mW2(R.k12r());
  Energy2 mZ2(R.k22r());
  // Then the rest:
  Energy2 s(R.sr());
  Energy2 tk(R.tkr());
  Energy2 uk(R.ukr());
  Energy2 q1(R.q1r());
  Energy2 q2(R.q2r());
  Energy2 q1h(R.q1hatr());
  Energy2 q2h(R.q2hatr());
  Energy2 w1(R.w1r());
  Energy2 w2(R.w2r());

  Energy2 Val(0.*GeV2);

  swap(s,uk);
  swap(q1,w1);
  swap(q1h,w2);
  Val  =  -2.*pi*alphaS_*Fij2_*CF_/NC_
        * (    gdL_*gdL_*t_u_Rdd(s,tk,uk,q1,q2,mW2,mZ2)
	  + 2.*gdL_*guL_*t_u_Rud(s,tk,uk,q1,q2,q1h,q2h,mW2,mZ2)
	  +    guL_*guL_*t_u_Ruu(s,tk,uk,q1h,q2h,mW2,mZ2)
	  - 2.*eZ_/(s2-mW2) * ( gdL_
		 	       *t_u_RZd(s,tk,uk,q1 ,q2 ,s2,mW2,mZ2)
	                      - guL_
			       *t_u_RZu(s,tk,uk,q1h,q2h,s2,mW2,mZ2)
	                      )
          + sqr(eZ_/(s2-mW2)) *t_u_RZ(s,tk,uk,q1,q2,s2,mW2,mZ2)
	  );
  swap(s,uk);
  swap(q1,w1);
  swap(q1h,w2);

  Val *= -uk/s * TR_/CF_; 

  return Val;
}

/***************************************************************************/
// The following six functions are I_{dd}^{(0)}, I_{ud}^{(0)}, 
// I_{uu}^{(0)}, F_{u}^{(0)}, F_{d}^{(0)}, H^{(0)} from Eqs. 3.9 - 3.14
// They make up the Born matrix element. Ixx functions correspond to the 
// graphs with no TGC, Fx functions are due to non-TGC graphs interfering 
// with TGC graphs, while the H function is due purely to TGC graphs. 
double Idd0(Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2);
double Iud0(Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2);
double Iuu0(Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2);
Energy2 Fu0(Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2);
Energy2 Fd0(Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2);
Energy4 H0 (Energy2 s,Energy2 t,Energy2 u,Energy2 mW2,Energy2 mZ2);

/***************************************************************************/
// M_V_Regular is the regular part of the one-loop matrix element 
// exactly as defined in Eqs. B.1 and B.2 of of NPB 383(1992)3-44.
double MEPP2VVPowheg::M_Born(born2to2Kinematics B) const {
  Energy2 s(B.sb());
  Energy2 t(B.tb());
  Energy2 u(B.ub());
  Energy2 mW2(B.k12b()); // N.B. the diboson masses are preserved in getting
  Energy2 mZ2(B.k22b()); // the 2->2 from the 2->3 kinematics.
  
  return Fij2_/2./NC_
       * (    gdL_*gdL_*Idd0(s,t,u,mW2,mZ2)
	 + 2.*gdL_*guL_*Iud0(s,t,u,mW2,mZ2)
	 +    guL_*guL_*Iuu0(s,t,u,mW2,mZ2) 
	 - 2.*eZ_/(s-mW2) * ( gdL_*Fd0(s,t,u,mW2,mZ2)
	                    - guL_*Fu0(s,t,u,mW2,mZ2)
	                    )
         + sqr(eZ_/(s-mW2)) * H0(s,t,u,mW2,mZ2)
	 );
}

/***************************************************************************/
double  Idd0(Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2) { 
    return    8.*((u*t/mW2/mZ2-1.)/4.+s/2.*(mW2+mZ2)/mW2/mZ2)
      	    + 8.*(u/t-mW2*mZ2/t/t);
}

/***************************************************************************/
double  Iud0(Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2) { 
    return  - 8.*((u*t/mW2/mZ2-1.)/4.+s/2.*(mW2+mZ2)/mW2/mZ2)
	    + 8.*s/t/u*(mW2+mZ2);
}

/***************************************************************************/
double  Iuu0(Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2) {
    return Idd0(s,u,t,mW2,mZ2);
}

/***************************************************************************/
Energy2 Fd0 (Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2) {
    return  - 8.*s*( (u*t/mW2/mZ2-1.)*(1.-(mW2+mZ2)/s-4.*mW2*mZ2/s/t)/4.
		   + (mW2+mZ2)/2./mW2/mZ2*(s-mW2-mZ2+2.*mW2*mZ2/t)
	           );
}

/***************************************************************************/
Energy2 Fu0 (Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2) { 
    return Fd0(s,u,t,mW2,mZ2);
}

/***************************************************************************/
Energy4 H0  (Energy2 s, Energy2 t, Energy2 u, Energy2 mW2, Energy2 mZ2) { 
    return    8.*s*s*(u*t/mW2/mZ2-1.)*( 1./4.-(mW2+mZ2)/2./s
				      + (sqr(mW2+mZ2)+8.*mW2*mZ2)/4./s/s
	                              )
	    + 8.*s*s*(mW2+mZ2)/mW2/mZ2*(s/2.-mW2-mZ2+sqr(mW2-mZ2)/2./s);
}
