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

void MEPP2VVPowheg::doinit() throw(InitException) {
  MEPP2VV::doinit();
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
  gW_ = sqrt(4.0*pi*SM().alphaEM(scale())/SM().sin2ThetaW());
  sin2ThetaW_ = SM().sin2ThetaW();
  double cosThetaW(sqrt(1.-sin2ThetaW_));
  guL_ = gW_/2./cosThetaW*( 1.-4./3.*sin2ThetaW_);
  gdL_ = gW_/2./cosThetaW*(-1.+2./3.*sin2ThetaW_);
  guR_ = gW_/2./cosThetaW*(   -4./3.*sin2ThetaW_);
  gdR_ = gW_/2./cosThetaW*(   +2./3.*sin2ThetaW_);
  eZ_  = gW_*cosThetaW;
  eZ2_ = sqr(eZ_);

  // If the process is W-Z instead of W+Z we must transform these
  // couplings as follows, according to NPB 383(1992)3-44 Eq.3.23
  if(mePartonData()[2]->id()==-24&&mePartonData()[3]->id()==23) { 
    swap(guL_,gdL_);
    eZ_ *= -1.;
  }

  // Get the CKM entry. Note that this code was debugged 
  // considerably; the call to CKM(particle,particle) 
  // did not appear to work, so we extract the elements
  // as follows below. The right numbers now appear to 
  // to be associated with the right quarks.
  double Kij(-999.);  
  // W+Z / W-Z
  if(abs(mePartonData()[2]->id())==24&&mePartonData()[3]->id()==23) { 
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
      cout << "MEPP2VVPowheg:" << endl;
      cout << "WZ needs an up and a down type quark as incoming!" << endl;
    }
    up_id /= 2;
    up_id -= 1;
    dn_id -= 1;
    dn_id /= 2;
    Kij = sqrt(SM().CKM(up_id,dn_id));
  } 
  // W+W-
  else if(abs(mePartonData()[2]->id())==24&&abs(mePartonData()[3]->id())==24) {
    if(!MEPP2VV::mixingInWW()) {
      Kij = 1.0;
    }
    else {
      if(abs(ab_->id())%2==0&&abs(bb_->id())%2==0) {
	int up_ida(abs(ab_->id())/2-1);
	int up_idb(abs(bb_->id())/2-1);
	Kij  = sqrt(std::norm( MEPP2VV::CKM(up_ida,0)*MEPP2VV::CKM(up_idb,0)
			     + MEPP2VV::CKM(up_ida,1)*MEPP2VV::CKM(up_idb,1)
			     + MEPP2VV::CKM(up_ida,2)*MEPP2VV::CKM(up_idb,2)));
      }
      else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) {
	int dn_ida((abs(ab_->id())-1)/2);
	int dn_idb((abs(bb_->id())-1)/2);
	Kij  = sqrt(std::norm( MEPP2VV::CKM(0,dn_ida)*MEPP2VV::CKM(0,dn_idb)
			     + MEPP2VV::CKM(1,dn_ida)*MEPP2VV::CKM(1,dn_idb)
			     + MEPP2VV::CKM(2,dn_ida)*MEPP2VV::CKM(2,dn_idb)));
      }
      else {
	cout << "MEPP2VVPowheg:" << endl;
	cout << "WW needs 2 down-type / 2 up-type!" << endl;
      }
    }
  }
  // ZZ 
  else if(mePartonData()[2]->id()==23&&mePartonData()[3]->id()==23) {
    Kij  = 2.*sqrt(2.)/gW_;
  }
  else {
    cout << "MEPP2VVPowheg: incompatible final state particles!" << endl;
  }

  Fij2_ = sqr(gW_/2./sqrt(2.)*Kij);

  // Get the leading order matrix element
  M_Born_      = M_Born_WZ(B_);
  // Get the regular part of the virtual correction
  M_V_regular_ = M_V_regular(S_);
  // Get the q + qbar real emission matrix element
  t_u_M_R_qqb_ = t_u_M_R_qqb(H_);

  // Calculate the integrand
  double wgt(0.);

  // q qb  contribution
  a_nlo=ab_;
  b_nlo=bb_;

  double wqqbvirt   = Vtilde_universal(S_) + M_V_regular(S_)/lo_me2_;
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

  // Debugging output:
  if(sanityCheck()) {
    //  sanityCheck();
    cout << "MEPP2VV::mixingInWW() " << MEPP2VV::mixingInWW() << endl;
    cout << ab_->PDGName() << ", " 
	 << bb_->PDGName() << ", " 
         << mePartonData()[2]->PDGName() << ", " 
	 << mePartonData()[3]->PDGName() << endl; 
    cout << "lo_me2_ - M_Born_ (rel) = " 
	 <<  lo_me2_-M_Born_                << "  (" 
	 << (lo_me2_-M_Born_)/M_Born_       << ")\n";
    cout << "lo_me2_, M_Born_    " << lo_me2_ << ", " << M_Born_    << endl;
    cout << "xr = " << H_.xr() << ", y = " << H_.y() << endl;
    cout << "root(sb) = " << sqrt(B_.sb())/GeV << endl;
    cout << "sb+tb+ub = " 
	 << B_.sb()/GeV2 << " + " 
	 << B_.tb()/GeV2 << " + " << B_.ub()/GeV2 << endl;
    cout << "sqrt(k12)  " << sqrt(H_.k12r())/GeV << endl;
    cout << "sqrt(k22)  " << sqrt(H_.k22r())/GeV << endl;
    cout << "sqr(Kij)   " << Kij*Kij    << endl;
    cout << "wqqbvirt   " << wqqbvirt   << endl;
    cout << "wqqbcollin " << wqqbcollin << endl;
    cout << "wqqbreal   " << wqqbreal   << endl;
    cout << "wqqb       " << wqqb       << endl;
    cout << "wqgcollin  " << wqgcollin  << endl;
    cout << "wqgreal    " << wqgreal    << endl;
    cout << "wqg        " << wqg        << endl;
    cout << "wgqbcollin " << wgqbcollin << endl;
    cout << "wgqbreal   " << wgqbreal   << endl;
    cout << "wgqb       " << wgqb       << endl;
    cout << "wgt        " << wgt        << endl;
  }

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

  Energy2 t_u_M_R_qqb_H (t_u_M_R_qqb(H_));
  Energy2 t_u_M_R_qqb_Cp(8.*pi*alphaS_*Cp_.sr()/Cp_.xr()
			*CF_*(1.+sqr(Cp_.xr()))*M_Born_);
  Energy2 t_u_M_R_qqb_Cm(8.*pi*alphaS_*Cm_.sr()/Cm_.xr()
			*CF_*(1.+sqr(Cm_.xr()))*M_Born_);

  return 
   ( ( (t_u_M_R_qqb_H*Lhat_ab(a,b,H_) - t_u_M_R_qqb_Cp*Lhat_ab(a,b,Cp_))/s
     - (t_u_M_R_qqb(S_)               - t_u_M_R_qqb(SCp_)              )/s2

     )*2./(1.-y)/(1.-xt)

   + ( (t_u_M_R_qqb_H*Lhat_ab(a,b,H_) - t_u_M_R_qqb_Cm*Lhat_ab(a,b,Cm_))/s
     - (t_u_M_R_qqb(S_)               - t_u_M_R_qqb(SCm_)              )/s2

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

  Energy2 t_u_M_R_gqb_H (t_u_M_R_gqb(H_));
  Energy2 t_u_M_R_gqb_Cp(8.*pi*alphaS_*Cp_.sr()/Cp_.xr()*(1.-Cp_.xr())
			*TR_*(sqr(Cp_.xr())+sqr(1.-Cp_.xr()))*M_Born_);
  Energy2 t_u_M_R_gqb_Cm(t_u_M_R_gqb(Cm_));

  return 
   ( ( (t_u_M_R_gqb_H*Lhat_ab(a,b,H_) - t_u_M_R_gqb_Cp*Lhat_ab(a,b,Cp_))/s
     - (t_u_M_R_gqb(S_)               - t_u_M_R_gqb(SCp_)              )/s2

     )*2./(1.-y)/(1.-xt)

   + ( (t_u_M_R_gqb_H*Lhat_ab(a,b,H_) - t_u_M_R_gqb_Cm*Lhat_ab(a,b,Cm_))/s
     - (t_u_M_R_gqb(S_)               - t_u_M_R_gqb(SCm_)              )/s2

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

  Energy2 t_u_M_R_qg_H (t_u_M_R_qg(H_));
  Energy2 t_u_M_R_qg_Cp(t_u_M_R_qg(Cp_));
  Energy2 t_u_M_R_qg_Cm(8.*pi*alphaS_*Cm_.sr()/Cm_.xr()*(1.-Cm_.xr())
		       *TR_*(sqr(Cm_.xr())+sqr(1.-Cm_.xr()))*M_Born_);

  return 
   ( ( (t_u_M_R_qg_H*Lhat_ab(a,b,H_) - t_u_M_R_qg_Cp*Lhat_ab(a,b,Cp_))/s
     - (t_u_M_R_qg(S_)               - t_u_M_R_qg(SCp_)                )/s2

     )*2./(1.-y)/(1.-xt)

   + ( (t_u_M_R_qg_H*Lhat_ab(a,b,H_) - t_u_M_R_qg_Cm*Lhat_ab(a,b,Cm_))/s
     - (t_u_M_R_qg(S_)               - t_u_M_R_qg(SCm_)                )/s2

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
  
  double cosThetaW(sqrt(1.-sin2ThetaW_));
  double eZ2(eZ2_);
  double eZ(eZ_);
  double gdL(gdL_);
  double guL(guL_);
  double gdR(gdR_);
  double guR(guR_);

  // W+W-
  if(abs(mePartonData()[2]->id())==24&&abs(mePartonData()[3]->id())==24) {
    double e2(sqr(gW_)*sin2ThetaW_);
    if(abs(ab_->id())%2==0&&abs(bb_->id())%2==0) {
      // N.B. OLD eZ used to calculate new eZ2 *then* new eZ is set!
      if(ab_->id()==-bb_->id()) {
        eZ2 = 1./2.*sqr(s-mW2)/Fij2_
	    * (e2*e2/s/s*(sqr( 2./3.+eZ*(guL+guR)/2./e2*s/(s-mW2/sqr(cosThetaW)))
		         +sqr(       eZ*(guL-guR)/2./e2*s/(s-mW2/sqr(cosThetaW))))
              );
        eZ  = -1./2./Fij2_/(gW_*gW_/4./sqrt(Fij2_))*(s-mW2)
            * (gW_*gW_*e2/4./s *( 2./3.+2.*eZ*guL/2./e2*s/(s-mW2/sqr(cosThetaW))));
      } else {
	eZ2 =0.;
	eZ  =0.;
      }
      gdL = gW_/sqrt(2.);
      guL = 0.;
    }
    else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) {
      // N.B. OLD eZ used to calculate new eZ2 *then* new eZ is set!
      if(ab_->id()==-bb_->id()) {
	eZ2 = 1./2.*sqr(s-mW2)/Fij2_
	    * (e2*e2/s/s*(sqr(-1./3.+eZ*(gdL+gdR)/2./e2*s/(s-mW2/sqr(cosThetaW)))
		         +sqr(       eZ*(gdL-gdR)/2./e2*s/(s-mW2/sqr(cosThetaW))))
	      );
	eZ  = -1./2./Fij2_/(gW_*gW_/4./sqrt(Fij2_))*(s-mW2)
            * (gW_*gW_*e2/4./s *(-1./3.+2.*eZ*gdL/2./e2*s/(s-mW2/sqr(cosThetaW))));
      } else {
	eZ2 =0.;
	eZ  =0.;
      }
      guL = gW_/sqrt(2.);
      gdL = 0.;
    }
  }
  // ZZ 
  else if(mePartonData()[2]->id()==23&&mePartonData()[3]->id()==23) {
    eZ  = 0.;
    eZ2 = 0.;
    double gV2,gA2;
    gV2 = sqr(guL/2.-gW_/2./cosThetaW*2./3.*sin2ThetaW_);
    gA2 = sqr(guL/2.+gW_/2./cosThetaW*2./3.*sin2ThetaW_);
    guL = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
    gV2 = sqr(gdL/2.+gW_/2./cosThetaW*1./3.*sin2ThetaW_);
    gA2 = sqr(gdL/2.-gW_/2./cosThetaW*1./3.*sin2ThetaW_);
    gdL = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
    if(abs(ab_->id())%2==0&&abs(bb_->id())%2==0)      gdL = guL;
    else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) guL = gdL;
    else {
      cout << "MEPP2VVPowheg:" << endl;
      cout << "ZZ needs 2 down-type / 2 up-type!" << endl;
    }
  }

  return 4.*pi*alphaS_*Fij2_*CF_*(1./sqr(4.*pi))/NC_
       * ( gdL*gdL*Idd1(s,t,u,mW2,mZ2,beta)
	 + gdL*guL*Iud1(s,t,u,mW2,mZ2,beta)
	 + guL*guL*Iuu1(s,t,u,mW2,mZ2,beta) 
	 - eZ/(s-mW2) * ( gdL*Fd1(s,t,u,mW2,mZ2,beta)
	                - guL*Fu1(s,t,u,mW2,mZ2,beta)
	                 )
         + eZ2/sqr(s-mW2) * H1(s,t,u,mW2,mZ2)
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

    Val +=    2.*(22.*t*t+t*(19.*s-18.*sig)+18.*mW2*mZ2)/t/t     
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
    Val +=    2.*(22.*t*t+t*(19.*s-18.*sig)+18.*mW2*mZ2)/t/t     
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

  double cosThetaW(sqrt(1.-sin2ThetaW_));
  double eZ2(eZ2_);
  double eZ(eZ_);
  double gdL(gdL_);
  double guL(guL_);
  double gdR(gdR_);
  double guR(guR_);

  // W+W-
  if(abs(mePartonData()[2]->id())==24&&abs(mePartonData()[3]->id())==24) {
    double e2(sqr(gW_)*sin2ThetaW_);
    if(abs(ab_->id())%2==0&&abs(bb_->id())%2==0) {
      // N.B. OLD eZ used to calculate new eZ2 *then* new eZ is set!
      if(ab_->id()==-bb_->id()) {
        eZ2 = 1./2.*sqr(s2-mW2)/Fij2_
	    * (e2*e2/s2/s2*(sqr( 2./3.+eZ*(guL+guR)/2./e2*s2/(s2-mW2/sqr(cosThetaW)))
		           +sqr(       eZ*(guL-guR)/2./e2*s2/(s2-mW2/sqr(cosThetaW))))
              );
        eZ  = -1./2./Fij2_/(gW_*gW_/4./sqrt(Fij2_))*(s2-mW2)
            * (gW_*gW_*e2/4./s2 *( 2./3.+2.*eZ*guL/2./e2*s2/(s2-mW2/sqr(cosThetaW))));
      } else {
	eZ2 =0.;
	eZ  =0.;
      }
      gdL = gW_/sqrt(2.);
      guL = 0.;
    }
    else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) {
      // N.B. OLD eZ used to calculate new eZ2 *then* new eZ is set!
      if(ab_->id()==-bb_->id()) {
	eZ2 = 1./2.*sqr(s2-mW2)/Fij2_
	    * (e2*e2/s2/s2*(sqr(-1./3.+eZ*(gdL+gdR)/2./e2*s2/(s2-mW2/sqr(cosThetaW)))
		           +sqr(       eZ*(gdL-gdR)/2./e2*s2/(s2-mW2/sqr(cosThetaW))))
	      );
	eZ  = -1./2./Fij2_/(gW_*gW_/4./sqrt(Fij2_))*(s2-mW2)
            * (gW_*gW_*e2/4./s2 *(-1./3.+2.*eZ*gdL/2./e2*s2/(s2-mW2/sqr(cosThetaW))));
      } else {
	eZ2 =0.;
	eZ  =0.;
      }
      guL = gW_/sqrt(2.);
      gdL = 0.;
    }
  }
  // ZZ 
  else if(mePartonData()[2]->id()==23&&mePartonData()[3]->id()==23) {
    eZ  = 0.;
    eZ2 = 0.;
    double gV2,gA2;
    gV2 = sqr(guL/2.-gW_/2./cosThetaW*2./3.*sin2ThetaW_);
    gA2 = sqr(guL/2.+gW_/2./cosThetaW*2./3.*sin2ThetaW_);
    guL = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
    gV2 = sqr(gdL/2.+gW_/2./cosThetaW*1./3.*sin2ThetaW_);
    gA2 = sqr(gdL/2.-gW_/2./cosThetaW*1./3.*sin2ThetaW_);
    gdL = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
    if(abs(ab_->id())%2==0&&abs(bb_->id())%2==0)      gdL = guL;
    else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) guL = gdL;
    else {
      cout << "MEPP2VVPowheg:" << endl;
      cout << "ZZ needs 2 down-type / 2 up-type!" << endl;
    }
  }

  return -2.*pi*alphaS_*Fij2_*CF_/NC_
       * (    gdL*gdL*t_u_Rdd(s,tk,uk,q1,q2,mW2,mZ2)
	 + 2.*gdL*guL*t_u_Rud(s,tk,uk,q1,q2,q1h,q2h,mW2,mZ2)
	 +    guL*guL*t_u_Ruu(s,tk,uk,q1h,q2h,mW2,mZ2)
	 - 2.*eZ/(s2-mW2) * ( gdL
			    * t_u_RZd(s,tk,uk,q1 ,q2 ,s2,mW2,mZ2)
	                    - guL
			    * t_u_RZu(s,tk,uk,q1h,q2h,s2,mW2,mZ2)
	                    )
         + eZ2/sqr(s2-mW2) *t_u_RZ(s,tk,uk,q1,q2,s2,mW2,mZ2)
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

  double cosThetaW(sqrt(1.-sin2ThetaW_));
  double eZ2(eZ2_);
  double eZ(eZ_);
  double gdL(gdL_);
  double guL(guL_);
  double gdR(gdR_);
  double guR(guR_);

  // W+W-
  if(abs(mePartonData()[2]->id())==24&&abs(mePartonData()[3]->id())==24) {
    double e2(sqr(gW_)*sin2ThetaW_);
    if(abs(ab_->id())%2==0&&abs(bb_->id())%2==0) {
      // N.B. OLD eZ used to calculate new eZ2 *then* new eZ is set!
      if(ab_->id()==-bb_->id()) {
        eZ2 = 1./2.*sqr(s2-mW2)/Fij2_
	    * (e2*e2/s2/s2*(sqr( 2./3.+eZ*(guL+guR)/2./e2*s2/(s2-mW2/sqr(cosThetaW)))
		           +sqr(       eZ*(guL-guR)/2./e2*s2/(s2-mW2/sqr(cosThetaW))))
              );
        eZ  = -1./2./Fij2_/(gW_*gW_/4./sqrt(Fij2_))*(s2-mW2)
            * (gW_*gW_*e2/4./s2 *( 2./3.+2.*eZ*guL/2./e2*s2/(s2-mW2/sqr(cosThetaW))));
      } else {
	eZ2 =0.;
	eZ  =0.;
      }
      gdL = gW_/sqrt(2.);
      guL = 0.;
    }
    else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) {
      // N.B. OLD eZ used to calculate new eZ2 *then* new eZ is set!
      if(ab_->id()==-bb_->id()) {
	eZ2 = 1./2.*sqr(s2-mW2)/Fij2_
	    * (e2*e2/s2/s2*(sqr(-1./3.+eZ*(gdL+gdR)/2./e2*s2/(s2-mW2/sqr(cosThetaW)))
		           +sqr(       eZ*(gdL-gdR)/2./e2*s2/(s2-mW2/sqr(cosThetaW))))
	      );
	eZ  = -1./2./Fij2_/(gW_*gW_/4./sqrt(Fij2_))*(s2-mW2)
            * (gW_*gW_*e2/4./s2 *(-1./3.+2.*eZ*gdL/2./e2*s2/(s2-mW2/sqr(cosThetaW))));
      } else {
	eZ2 =0.;
	eZ  =0.;
      }
      guL = gW_/sqrt(2.);
      gdL = 0.;
    }
  }
  // ZZ 
  else if(mePartonData()[2]->id()==23&&mePartonData()[3]->id()==23) {
    eZ  = 0.;
    eZ2 = 0.;
    double gV2,gA2;
    gV2 = sqr(guL/2.-gW_/2./cosThetaW*2./3.*sin2ThetaW_);
    gA2 = sqr(guL/2.+gW_/2./cosThetaW*2./3.*sin2ThetaW_);
    guL = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
    gV2 = sqr(gdL/2.+gW_/2./cosThetaW*1./3.*sin2ThetaW_);
    gA2 = sqr(gdL/2.-gW_/2./cosThetaW*1./3.*sin2ThetaW_);
    gdL = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
    if(abs(ab_->id())%2==0&&abs(bb_->id())%2==0)      gdL = guL;
    else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) guL = gdL;
    else {
      cout << "MEPP2VVPowheg:" << endl;
      cout << "ZZ needs 2 down-type / 2 up-type!" << endl;
    }
  }

  Energy2 Val(0.*GeV2);

  swap(s,tk);
  swap(q2,w2);
  swap(q2h,w1);
  Val  =  -2.*pi*alphaS_*Fij2_*CF_/NC_
        * (    gdL*gdL*t_u_Rdd(s,tk,uk,q1,q2,mW2,mZ2)
	  + 2.*gdL*guL*t_u_Rud(s,tk,uk,q1,q2,q1h,q2h,mW2,mZ2)
	  +    guL*guL*t_u_Ruu(s,tk,uk,q1h,q2h,mW2,mZ2)
	  - 2.*eZ/(s2-mW2) * ( gdL
		 	     * t_u_RZd(s,tk,uk,q1 ,q2 ,s2,mW2,mZ2)
	                     - guL
			     * t_u_RZu(s,tk,uk,q1h,q2h,s2,mW2,mZ2)
	                     )
          + eZ2/sqr(s2-mW2) *t_u_RZ(s,tk,uk,q1,q2,s2,mW2,mZ2)
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

  double cosThetaW(sqrt(1.-sin2ThetaW_));
  double eZ2(eZ2_);
  double eZ(eZ_);
  double gdL(gdL_);
  double guL(guL_);
  double gdR(gdR_);
  double guR(guR_);

  // W+W-
  if(abs(mePartonData()[2]->id())==24&&abs(mePartonData()[3]->id())==24) {
    double e2(sqr(gW_)*sin2ThetaW_);
    if(abs(ab_->id())%2==0&&abs(bb_->id())%2==0) {
      // N.B. OLD eZ used to calculate new eZ2 *then* new eZ is set!
      if(ab_->id()==-bb_->id()) {
        eZ2 = 1./2.*sqr(s2-mW2)/Fij2_
	    * (e2*e2/s2/s2*(sqr( 2./3.+eZ*(guL+guR)/2./e2*s2/(s2-mW2/sqr(cosThetaW)))
		           +sqr(       eZ*(guL-guR)/2./e2*s2/(s2-mW2/sqr(cosThetaW))))
              );
        eZ  = -1./2./Fij2_/(gW_*gW_/4./sqrt(Fij2_))*(s2-mW2)
            * (gW_*gW_*e2/4./s2 *( 2./3.+2.*eZ*guL/2./e2*s2/(s2-mW2/sqr(cosThetaW))));
      } else {
	eZ2 =0.;
	eZ  =0.;
      }
      gdL = gW_/sqrt(2.);
      guL = 0.;
    }
    else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) {
      // N.B. OLD eZ used to calculate new eZ2 *then* new eZ is set!
      if(ab_->id()==-bb_->id()) {
	eZ2 = 1./2.*sqr(s2-mW2)/Fij2_
	    * (e2*e2/s2/s2*(sqr(-1./3.+eZ*(gdL+gdR)/2./e2*s2/(s2-mW2/sqr(cosThetaW)))
		           +sqr(       eZ*(gdL-gdR)/2./e2*s2/(s2-mW2/sqr(cosThetaW))))
	      );
	eZ  = -1./2./Fij2_/(gW_*gW_/4./sqrt(Fij2_))*(s2-mW2)
            * (gW_*gW_*e2/4./s2 *(-1./3.+2.*eZ*gdL/2./e2*s2/(s2-mW2/sqr(cosThetaW))));
      } else {
	eZ2 =0.;
	eZ  =0.;
      }
      guL = gW_/sqrt(2.);
      gdL = 0.;
    }
  }
  // ZZ 
  else if(mePartonData()[2]->id()==23&&mePartonData()[3]->id()==23) {
    eZ  = 0.;
    eZ2 = 0.;
    double gV2,gA2;
    gV2 = sqr(guL/2.-gW_/2./cosThetaW*2./3.*sin2ThetaW_);
    gA2 = sqr(guL/2.+gW_/2./cosThetaW*2./3.*sin2ThetaW_);
    guL = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
    gV2 = sqr(gdL/2.+gW_/2./cosThetaW*1./3.*sin2ThetaW_);
    gA2 = sqr(gdL/2.-gW_/2./cosThetaW*1./3.*sin2ThetaW_);
    gdL = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
    if(abs(ab_->id())%2==0&&abs(bb_->id())%2==0)      gdL = guL;
    else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) guL = gdL;
    else {
      cout << "MEPP2VVPowheg:" << endl;
      cout << "ZZ needs 2 down-type / 2 up-type!" << endl;
    }
  }

  Energy2 Val(0.*GeV2);

  swap(s,uk);
  swap(q1,w1);
  swap(q1h,w2);
  Val  =  -2.*pi*alphaS_*Fij2_*CF_/NC_
        * (    gdL*gdL*t_u_Rdd(s,tk,uk,q1,q2,mW2,mZ2)
	  + 2.*gdL*guL*t_u_Rud(s,tk,uk,q1,q2,q1h,q2h,mW2,mZ2)
	  +    guL*guL*t_u_Ruu(s,tk,uk,q1h,q2h,mW2,mZ2)
	  - 2.*eZ/(s2-mW2) * ( gdL
		 	     * t_u_RZd(s,tk,uk,q1 ,q2 ,s2,mW2,mZ2)
	                     - guL
			     * t_u_RZu(s,tk,uk,q1h,q2h,s2,mW2,mZ2)
	                     )
          + eZ2/sqr(s2-mW2) *t_u_RZ(s,tk,uk,q1,q2,s2,mW2,mZ2)
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
// M_Born_WZ is the Born matrix element exactly as defined in Eqs. 3.3-3.14
// of of NPB 383(1992)3-44.
double MEPP2VVPowheg::M_Born_WZ(born2to2Kinematics B) const {
  Energy2 s(B.sb());
  Energy2 t(B.tb());
  Energy2 u(B.ub());
  Energy2 mW2(B.k12b()); // N.B. the diboson masses are preserved in getting
  Energy2 mZ2(B.k22b()); // the 2->2 from the 2->3 kinematics.
  
  double cosThetaW(sqrt(1.-sin2ThetaW_));
  double eZ2(eZ2_);
  double eZ(eZ_);
  double gdL(gdL_);
  double guL(guL_);
  double gdR(gdR_);
  double guR(guR_);

  // W+W-
  if(abs(mePartonData()[2]->id())==24&&abs(mePartonData()[3]->id())==24) {
    double e2(sqr(gW_)*sin2ThetaW_);
    if(abs(ab_->id())%2==0&&abs(bb_->id())%2==0) {
      // N.B. OLD eZ used to calculate new eZ2 *then* new eZ is set!
      if(ab_->id()==-bb_->id()) {
        eZ2 = 1./2.*sqr(s-mW2)/Fij2_
	    * (e2*e2/s/s*(sqr( 2./3.+eZ*(guL+guR)/2./e2*s/(s-mW2/sqr(cosThetaW)))
		         +sqr(       eZ*(guL-guR)/2./e2*s/(s-mW2/sqr(cosThetaW))))
              );
        eZ  = -1./2./Fij2_/(gW_*gW_/4./sqrt(Fij2_))*(s-mW2)
            * (gW_*gW_*e2/4./s *( 2./3.+2.*eZ*guL/2./e2*s/(s-mW2/sqr(cosThetaW))));
      } else {
	eZ2 =0.;
	eZ  =0.;
      }
      gdL = gW_/sqrt(2.);
      guL = 0.;
    }
    else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) {
      // N.B. OLD eZ used to calculate new eZ2 *then* new eZ is set!
      if(ab_->id()==-bb_->id()) {
	eZ2 = 1./2.*sqr(s-mW2)/Fij2_
	    * (e2*e2/s/s*(sqr(-1./3.+eZ*(gdL+gdR)/2./e2*s/(s-mW2/sqr(cosThetaW)))
		         +sqr(       eZ*(gdL-gdR)/2./e2*s/(s-mW2/sqr(cosThetaW))))
	      );
	eZ  = -1./2./Fij2_/(gW_*gW_/4./sqrt(Fij2_))*(s-mW2)
            * (gW_*gW_*e2/4./s *(-1./3.+2.*eZ*gdL/2./e2*s/(s-mW2/sqr(cosThetaW))));
      } else {
	eZ2 =0.;
	eZ  =0.;
      }
      guL = gW_/sqrt(2.);
      gdL = 0.;
    }
  }
  // ZZ 
  else if(mePartonData()[2]->id()==23&&mePartonData()[3]->id()==23) {
    eZ  = 0.;
    eZ2 = 0.;
    double gV2,gA2;
    gV2 = sqr(guL/2.-gW_/2./cosThetaW*2./3.*sin2ThetaW_);
    gA2 = sqr(guL/2.+gW_/2./cosThetaW*2./3.*sin2ThetaW_);
    guL = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
    gV2 = sqr(gdL/2.+gW_/2./cosThetaW*1./3.*sin2ThetaW_);
    gA2 = sqr(gdL/2.-gW_/2./cosThetaW*1./3.*sin2ThetaW_);
    gdL = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
    if(abs(ab_->id())%2==0&&abs(bb_->id())%2==0)      gdL = guL;
    else if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) guL = gdL;
    else {
      cout << "MEPP2VVPowheg:" << endl;
      cout << "ZZ needs 2 down-type / 2 up-type!" << endl;
    }
  }

  return Fij2_/2./NC_
       * (    
	      gdL*gdL*Idd0(s,t,u,mW2,mZ2)
	 + 2.*gdL*guL*Iud0(s,t,u,mW2,mZ2)
	 +    guL*guL*Iuu0(s,t,u,mW2,mZ2) 
	 - 2.*eZ/(s-mW2) * ( gdL*Fd0(s,t,u,mW2,mZ2)
	                   - guL*Fu0(s,t,u,mW2,mZ2)
	                   )
         + eZ2/sqr(s-mW2) * H0(s,t,u,mW2,mZ2)
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

/***************************************************************************/
bool MEPP2VVPowheg::sanityCheck() const {

  bool alarm(false);

  Energy2 prefacs(8.*pi*alphaS_*S_.sr() /S_.xr() );
  Energy2 prefacsp(8.*pi*alphaS_*SCp_.sr() /SCp_.xr() );
  Energy2 prefacsm(8.*pi*alphaS_*SCm_.sr() /SCm_.xr() );
  Energy2 prefacp(8.*pi*alphaS_*Cp_.sr()/Cp_.xr());
  Energy2 prefacm(8.*pi*alphaS_*Cm_.sr()/Cm_.xr());

  double xp(Cp_.xr());
  double xm(Cm_.xr());

  double M_B_WW(M_Born_WW(B_));
  double M_B_ZZ(M_Born_ZZ(B_));
  double M_V_reg_WW(M_V_regular_WW(S_));
  double M_V_reg_ZZ(M_V_regular_ZZ(S_));
  Energy2 t_u_qqb_WW(t_u_M_R_qqb_WW(H_));
  Energy2 t_u_qqb_ZZ(t_u_M_R_qqb_ZZ(H_));

  // Check that the native leading order Herwig++ matrix 
  // element is equivalent to the WZ leading order matrix 
  // element in NPB 383 (1992) 3-44, with the relevant WZ->WW
  // WZ->ZZ transformation applied (M_Born_).
//   if(fabs((lo_me2_ - M_Born_)/M_Born_)>1.e-2) {
//     alarm=true;
//     cout << "lo_me2_ - M_Born_ (rel) = " 
// 	 <<  lo_me2_ - M_Born_          << "  (" 
// 	 << (lo_me2_ - M_Born_)/M_Born_ << ")\n";
//   }

  // Check that the transformation from NPB 383 (1992) 3-44 WZ 
  // matrix elements to WW matrix elements actually works, by
  // comparing them to the explicit WW matrix elements in 
  // NPB 410 (1993) 280-324.
  if(abs(mePartonData()[2]->id())==24&&abs(mePartonData()[3]->id())==24) {
    if(fabs((M_Born_     -M_B_WW    )/M_B_WW    )>1.e-6) {
	alarm=true;
	cout << "WZ->WW transformation error!\n"; 
	cout << "M_Born_ - M_B_WW (rel) = " 
	     <<  M_Born_ - M_B_WW         << "  (" 
	     << (M_Born_ - M_B_WW)/M_B_WW << ")\n";
	cout << "M_Born_ = " << M_Born_   << endl;
	cout << "M_B_WW  = " << M_B_WW    << endl;
    }
    if(fabs((M_V_regular_-M_V_reg_WW)/M_V_reg_WW)>1.e-6) {
	alarm=true;
	cout << "WZ->WW transformation error!\n"; 
	cout << "M_V_regular_ - M_V_reg_WW (rel) = " 
	     <<  M_V_regular_ - M_V_reg_WW             << "  (" 
	     << (M_V_regular_ - M_V_reg_WW)/M_V_reg_WW << ")\n";
	cout << "M_V_regular_   = " << M_V_regular_    << endl;
	cout << "M_V_reg_WW     = " << M_V_reg_WW      << endl;
    }
    if(fabs((t_u_M_R_qqb_-t_u_qqb_WW)/t_u_qqb_WW)>1.e-6) {
	alarm=true;
	cout << "WZ->WW transformation error!\n"; 
	cout << "t_u_M_R_qqb_ - t_u_qqb_WW (rel) = " 
	     << (t_u_M_R_qqb_ - t_u_qqb_WW)/GeV2       << "  (" 
	     << (t_u_M_R_qqb_ - t_u_qqb_WW)/t_u_qqb_WW << ")\n";
	cout << "t_u_M_R_qqb_ = " << t_u_M_R_qqb_/GeV2 << endl;
	cout << "t_u_qqb_WW   = " << t_u_qqb_WW  /GeV2 << endl;
    }
  }

  // Check that the transformation from NPB 383 (1992) 3-44 WZ 
  // matrix elements to ZZ matrix elements actually works, by
  // comparing them to the explicit ZZ matrix elements in 
  // NPB 357 (1991) 409-438.
  if(abs(mePartonData()[2]->id())==23&&abs(mePartonData()[3]->id())==23) {
    if(fabs((M_Born_     -M_B_ZZ    )/M_B_ZZ    )>1.e-6) {
	alarm=true;
	cout << "WZ->ZZ transformation error!\n"; 
	cout << "M_Born_ - M_B_ZZ (rel) = " 
	     <<  M_Born_ - M_B_ZZ         << "  (" 
	     << (M_Born_ - M_B_ZZ)/M_B_ZZ << ")\n";
	cout << "M_Born_ = " << M_Born_   << endl;
	cout << "M_B_ZZ  = " << M_B_ZZ    << endl;
    }
    if(fabs((M_V_regular_-M_V_reg_ZZ)/M_V_reg_ZZ)>1.e-6) {
	alarm=true;
	cout << "WZ->ZZ transformation error!\n"; 
	cout << "M_V_regular_ - M_V_reg_ZZ (rel) = " 
	     <<  M_V_regular_ - M_V_reg_ZZ             << "  (" 
	     << (M_V_regular_ - M_V_reg_ZZ)/M_V_reg_ZZ << ")\n";
	cout << "M_V_regular_ = " << M_V_regular_      << endl;
	cout << "M_V_reg_ZZ   = " << M_V_reg_ZZ        << endl;
    }
    if(fabs((t_u_M_R_qqb_-t_u_qqb_ZZ)/t_u_qqb_ZZ)>1.e-6) {
	alarm=true;
	cout << "WZ->ZZ transformation error!\n"; 
	cout << "t_u_M_R_qqb_ - t_u_qqb_ZZ (rel) = " 
	     << (t_u_M_R_qqb_ - t_u_qqb_ZZ)/GeV2       << "  (" 
	     << (t_u_M_R_qqb_ - t_u_qqb_ZZ)/t_u_qqb_ZZ << ")\n";
	cout << "t_u_M_R_qqb_ = " << t_u_M_R_qqb_/GeV2 << endl;
	cout << "t_u_qqb_ZZ   = " << t_u_qqb_ZZ  /GeV2 << endl;
    }
  }

  // Check the soft limit of the q + qbar matrix element.
  Energy2 absDiff_qqbs 
      = t_u_M_R_qqb(S_) - prefacs*2.*CF_*M_Born_;
  double  relDiff_qqbs = absDiff_qqbs / t_u_M_R_qqb(S_);
  if(fabs(relDiff_qqbs)>1.e-7) {
    alarm=true;
    cout << "\n";
    cout << "t_u_M_R_qqb(S_)    " << t_u_M_R_qqb(S_)  /GeV2 << endl;
    cout << "t_u_M_R_qqb(S_)-8*pi*alphaS*sHat/x*2*Cab*M_Born_ (rel):\n"
	 << absDiff_qqbs / GeV2 << "   (" << relDiff_qqbs << ")\n";
  }

  // Check the positive soft-collinearlimit of the q + qbar matrix element.
  Energy2 absDiff_qqbsp 
      = t_u_M_R_qqb(SCp_) - prefacsp*2.*CF_*M_Born_;
  double  relDiff_qqbsp = absDiff_qqbsp / t_u_M_R_qqb(SCp_);
  if(fabs(relDiff_qqbsp)>1.e-7) {
    alarm=true;
    cout << "\n";
    cout << "t_u_M_R_qqb(SCp_)  " << t_u_M_R_qqb(SCp_)/GeV2 << endl;
    cout << "t_u_M_R_qqb(SCp_)-8*pi*alphaS*sHat/x*2*Cab*M_Born_ (rel):\n"
	 << absDiff_qqbsp / GeV2 << "   (" << relDiff_qqbsp << ")\n";
  }

  // Check the negative soft-collinearlimit of the q + qbar matrix element.
  Energy2 absDiff_qqbsm 
      = t_u_M_R_qqb(SCm_) - prefacsm*2.*CF_*M_Born_;
  double  relDiff_qqbsm = absDiff_qqbsm / t_u_M_R_qqb(SCm_);
  if(fabs(relDiff_qqbsm)>1.e-7) {
    alarm=true;
    cout << "\n";
    cout << "t_u_M_R_qqb(SCm_)  " << t_u_M_R_qqb(SCm_)/GeV2 << endl;
    cout << "t_u_M_R_qqb(SCm_)-8*pi*alphaS*sHat/x*2*Cab*M_Born_ (rel):\n"
	 << absDiff_qqbsm / GeV2 << "   (" << relDiff_qqbsm << ")\n";
  }

  // Check the positive collinearlimit of the q + qbar matrix element.
  Energy2 absDiff_qqbp 
      = t_u_M_R_qqb(Cp_) - prefacp*CF_*(1.+xp*xp)*M_Born_;
  double  relDiff_qqbp = absDiff_qqbp / t_u_M_R_qqb(Cp_);
  if(fabs(relDiff_qqbp)>1.e-7) {
    alarm=true;
    cout << "\n";
    cout << "t_u_M_R_qqb(Cp_)   " << t_u_M_R_qqb(Cp_) /GeV2 << endl;
    cout << "t_u_M_R_qqb(Cp_)-8*pi*alphaS*sHat/x*(1-x)*Pqq*M_Born_ (rel):\n"
	 << absDiff_qqbp / GeV2 << "   (" << relDiff_qqbp << ")\n";
  }

  // Check the negative collinearlimit of the q + qbar matrix element.
  Energy2 absDiff_qqbm 
      = t_u_M_R_qqb(Cm_) - prefacm*CF_*(1.+xm*xm)*M_Born_;
  double  relDiff_qqbm = absDiff_qqbm / t_u_M_R_qqb(Cm_);
  if(fabs(relDiff_qqbm)>1.e-7) {
    alarm=true;
    cout << "\n";
    cout << "t_u_M_R_qqb(Cm_)   " << t_u_M_R_qqb(Cm_) /GeV2 << endl;
    cout << "t_u_M_R_qqb(Cm_)-8*pi*alphaS*sHat/x*(1-x)*Pqq*M_Born_ (rel):\n"
	 << absDiff_qqbm / GeV2 << "   (" << relDiff_qqbm << ")\n";
  }

  // Check the positive collinear limit of the g + qbar matrix element.
  Energy2 absDiff_gqbp
      = t_u_M_R_gqb(Cp_) - prefacp*(1.-xp)*TR_*(xp*xp+sqr(1.-xp))*M_Born_;
  double  relDiff_gqbp =  absDiff_gqbp/ t_u_M_R_gqb(Cp_);
  if(fabs(relDiff_gqbp)>1.e-7) {
    alarm=true;
    cout << "\n";
    cout << "t_u_M_R_gqb(Cp_)   " << t_u_M_R_gqb(Cp_) /GeV2 << endl;
    cout << "t_u_M_R_gqb(Cp_)-8*pi*alphaS*sHat/x*(1-x)*Pgq*M_Born_ (rel):\n"
	 << absDiff_gqbp / GeV2 << "   (" << relDiff_gqbp << ")\n";
  }

  // Check the negative collinear limit of the q + g matrix element.
  Energy2 absDiff_qgm
      = t_u_M_R_qg(Cm_)  - prefacm*(1.-xm)*TR_*(xm*xm+sqr(1.-xm))*M_Born_;
  double  relDiff_qgm  =  absDiff_qgm / t_u_M_R_qg(Cm_);
  if(fabs(relDiff_qgm)>1.e-7) {
    alarm=true;
    cout << "\n";
    cout << "t_u_M_R_qg(Cm_)   " << t_u_M_R_qg(Cm_) /GeV2 << endl;
    cout << "t_u_M_R_qg(Cm_)-8*pi*alphaS*sHat/x*(1-x)*Pgq*M_Born_ (rel):\n"
	 << absDiff_qgm  / GeV2 << "   (" << relDiff_qgm  << ")\n";
  }

  return alarm;
}

/***************************************************************************/
// M_Born_ZZ is the Born matrix element exactly as defined in Eqs. 2.18-2.19
// of of NPB 357(1991)409-438.
double MEPP2VVPowheg::M_Born_ZZ(born2to2Kinematics B) const {
  Energy2 s(B.sb());
  Energy2 t(B.tb());
  Energy2 u(B.ub());
  Energy2 mW2(B.k12b()); // N.B. the diboson masses are preserved in getting
  Energy2 mZ2(B.k22b()); // the 2->2 from the 2->3 kinematics.
  double cosThetaW(sqrt(1.-sin2ThetaW_));
  
  double gV2,gA2,gX,gY,gZ;
  gV2  = sqr(guL_/2.-gW_/2./cosThetaW*2./3.*sin2ThetaW_);
  gA2  = sqr(guL_/2.+gW_/2./cosThetaW*2./3.*sin2ThetaW_);
  gX   = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
  gV2  = sqr(gdL_/2.+gW_/2./cosThetaW*1./3.*sin2ThetaW_);
  gA2  = sqr(gdL_/2.-gW_/2./cosThetaW*1./3.*sin2ThetaW_);
  gY   = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
  gZ   = gX;
  if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) gZ = gY;
  
  return 1./NC_*sqr(gZ*2.)*(t/u+u/t+4.*mZ2*s/t/u-mZ2*mZ2*(1./t/t+1./u/u));

}

/***************************************************************************/
// M_V_regular_ZZ is the one-loop ZZ matrix element exactly as defined in 
// Eqs. B.1 & B.2 of NPB 357(1991)409-438.
double MEPP2VVPowheg::M_V_regular_ZZ(real2to3Kinematics S) const {
  Energy2 s(S.bornVariables().sb());
  Energy2 t(S.bornVariables().tb());
  Energy2 u(S.bornVariables().ub());
  Energy2 mW2(S.k12r()); // N.B. the diboson masses are preserved in getting
  Energy2 mZ2(S.k22r()); // the 2->2 from the 2->3 kinematics.
  double  beta(S.betaxr()); // N.B. for x=1 \beta_x=\beta in NPB 383(1992)3-44.
  double cosThetaW(sqrt(1.-sin2ThetaW_));
  
  double gV2,gA2,gX,gY,gZ;
  gV2  = sqr(guL_/2.-gW_/2./cosThetaW*2./3.*sin2ThetaW_);
  gA2  = sqr(guL_/2.+gW_/2./cosThetaW*2./3.*sin2ThetaW_);
  gX   = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
  gV2  = sqr(gdL_/2.+gW_/2./cosThetaW*1./3.*sin2ThetaW_);
  gA2  = sqr(gdL_/2.-gW_/2./cosThetaW*1./3.*sin2ThetaW_);
  gY   = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
  gZ   = gX;
  if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) gZ = gY;
  
  double M_V_reg(0.);
  M_V_reg  = 2.*s*sqr(gZ*2.)*4.*pi*alphaS_*CF_/NC_/sqr(4.*pi)/2.
    *(   2.*sqr(t+mZ2)/sqr(beta)/s/t/u + 4.*s/(t-mZ2)/u
       - ( 16.*t*t*t+(28.*s-68.*mZ2)*t*t+(18.*s*s-36.*mZ2*s+88.*mZ2*mZ2)*t
	 + 18.*mZ2*mZ2*s-36.*mZ2*mZ2*mZ2
	 )/t/t/s/u
       + ( 12.*s/(t-mZ2)/u-4.*mZ2*s/sqr(t-mZ2)/u+2.*(t+4.*s)/s/u
	 - 6.*(s*s+mZ2*mZ2)/s/t/u+6.*mZ2*mZ2*(2.*mZ2-s)/t/t/s/u 
	 )*log(-t/mZ2)
       + ( - ( 5.*t*t*t+(8.*s-18.*mZ2)*t*t+(6.*s*s+25.*mZ2*mZ2)*t
	     + 6.*mZ2*mZ2*s-12.*mZ2*mZ2*mZ2
	     )/t/t/s/u
	   - 12.*mZ2*sqr(t+mZ2)/sqr(sqr(beta))/s/s/t/u
	   + ( 3.*t*t-26.*mZ2*t-25.*mZ2*mZ2)/sqr(beta)/s/t/u
	 )*log(s/mZ2)
       + ( (-2.*t*t+8.*mZ2*t-2.*s*s-12.*mZ2*mZ2)/u + 4.*mZ2*mZ2*(2.*mZ2-s)/t/u)
       / (s*t)
       * ( 2.*sqr(log(-t/mZ2))-4.*log(-t/mZ2)*log((mZ2-t)/mZ2)-4.*ReLi2(t/mZ2))
       + ( 4.*(t*t-5.*mZ2*t+s*s+10.*mZ2*mZ2)/s/u
         + 4.*mZ2*(-s*s+2.*mZ2*s-10.*mZ2*mZ2)/s/t/u
         + 8.*mZ2*mZ2*mZ2*(2.*mZ2-s)/t/t/s/u
	 )
       / (t-mZ2)
       * (pi*pi/2.+log(-t/mZ2)*log(-t/s)-1./2.*sqr(log(-t/mZ2)))
       + ( ( (2.*s-3.*mZ2)*t*t+(6.*mZ2*mZ2-8.*mZ2*s)*t+2.*s*s*s-4.*mZ2*s*s
	   + 12.*mZ2*mZ2*s-3.*mZ2*mZ2*mZ2
	   ) /s/t/u
	 + 12.*mZ2*mZ2*sqr(t+mZ2)/sqr(sqr(beta))/s/s/t/u
	 - (mZ2*t*t-30.*mZ2*mZ2*t-27.*mZ2*mZ2*mZ2)/beta/beta/s/t/u
	 )
       / (beta*s)
       * (pi*pi/3.+sqr(log((1.-beta)/(1.+beta)))+4.*ReLi2(-(1.-beta)/(1.+beta)))
       + (4.*(t+4.*s-4.*mZ2)/3./s/u+4.*sqr(s-2.*mZ2)/3./s/t/u)*pi*pi
     );

  swap(t,u);
  M_V_reg += 2.*s*sqr(gZ*2.)*4.*pi*alphaS_*CF_/NC_/sqr(4.*pi)/2.
    *(   2.*sqr(t+mZ2)/sqr(beta)/s/t/u + 4.*s/(t-mZ2)/u
       - ( 16.*t*t*t+(28.*s-68.*mZ2)*t*t+(18.*s*s-36.*mZ2*s+88.*mZ2*mZ2)*t
	 + 18.*mZ2*mZ2*s-36.*mZ2*mZ2*mZ2
	 )/t/t/s/u
       + ( 12.*s/(t-mZ2)/u-4.*mZ2*s/sqr(t-mZ2)/u+2.*(t+4.*s)/s/u
	 - 6.*(s*s+mZ2*mZ2)/s/t/u+6.*mZ2*mZ2*(2.*mZ2-s)/t/t/s/u 
	 )*log(-t/mZ2)
       + ( - ( 5.*t*t*t+(8.*s-18.*mZ2)*t*t+(6.*s*s+25.*mZ2*mZ2)*t
	     + 6.*mZ2*mZ2*s-12.*mZ2*mZ2*mZ2
	     )/t/t/s/u
	   - 12.*mZ2*sqr(t+mZ2)/sqr(sqr(beta))/s/s/t/u
	   + ( 3.*t*t-26.*mZ2*t-25.*mZ2*mZ2)/sqr(beta)/s/t/u
	 )*log(s/mZ2)
       + ( (-2.*t*t+8.*mZ2*t-2.*s*s-12.*mZ2*mZ2)/u + 4.*mZ2*mZ2*(2.*mZ2-s)/t/u)
       / (s*t)
       * ( 2.*sqr(log(-t/mZ2))-4.*log(-t/mZ2)*log((mZ2-t)/mZ2)-4.*ReLi2(t/mZ2))
       + ( 4.*(t*t-5.*mZ2*t+s*s+10.*mZ2*mZ2)/s/u
         + 4.*mZ2*(-s*s+2.*mZ2*s-10.*mZ2*mZ2)/s/t/u
         + 8.*mZ2*mZ2*mZ2*(2.*mZ2-s)/t/t/s/u
	 )
       / (t-mZ2)
       * (pi*pi/2.+log(-t/mZ2)*log(-t/s)-1./2.*sqr(log(-t/mZ2)))
       + ( ( (2.*s-3.*mZ2)*t*t+(6.*mZ2*mZ2-8.*mZ2*s)*t+2.*s*s*s-4.*mZ2*s*s
	   + 12.*mZ2*mZ2*s-3.*mZ2*mZ2*mZ2
	   ) /s/t/u
	 + 12.*mZ2*mZ2*sqr(t+mZ2)/sqr(sqr(beta))/s/s/t/u
	 - (mZ2*t*t-30.*mZ2*mZ2*t-27.*mZ2*mZ2*mZ2)/beta/beta/s/t/u
	 )
       / (beta*s)
       * (pi*pi/3.+sqr(log((1.-beta)/(1.+beta)))+4.*ReLi2(-(1.-beta)/(1.+beta)))
       + (4.*(t+4.*s-4.*mZ2)/3./s/u+4.*sqr(s-2.*mZ2)/3./s/t/u)*pi*pi
     );

  return M_V_reg;
}

/***************************************************************************/
// t_u_M_R_qqb_ZZ is the real emission q + qb -> n + g matrix element 
// exactly as defined in Eqs. C.1 of NPB 357(1991)409-438, multiplied by
// tk * uk!
Energy2 MEPP2VVPowheg::t_u_M_R_qqb_ZZ(real2to3Kinematics R) const {
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

  double cosThetaW(sqrt(1.-sin2ThetaW_));
  
  double gV2,gA2,gX,gY,gZ;
  gV2  = sqr(guL_/2.-gW_/2./cosThetaW*2./3.*sin2ThetaW_);
  gA2  = sqr(guL_/2.+gW_/2./cosThetaW*2./3.*sin2ThetaW_);
  gX   = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
  gV2  = sqr(gdL_/2.+gW_/2./cosThetaW*1./3.*sin2ThetaW_);
  gA2  = sqr(gdL_/2.-gW_/2./cosThetaW*1./3.*sin2ThetaW_);
  gY   = sqrt(gV2*gV2+gA2*gA2+6.*gA2*gV2)/2.;
  gZ   = gX;
  if(abs(ab_->id())%2==1&&abs(bb_->id())%2==1) gZ = gY;

  Energy2 t_u_qqb(0.*GeV2);
  t_u_qqb  = (2.*s)*sqr(gZ*2.)*4.*pi*alphaS_*CF_/NC_/2.
    * ( - ( tk*uk*uk+2.*s*uk*uk-tk*tk*uk
	  - 2.*s*tk*uk+mZ2*(tk*tk-uk*uk+2.*s*uk-2.*s*tk-2.*s*s)
	  )/q1h/q1/q2h/s*tk
	+ 2.*(tk*uk*uk-mZ2*uk*(s+3.*tk)+mZ2*mZ2*(2.*uk-s))/q1/q2/s
	+ ( tk*uk*(uk+s)-mZ2*(uk*uk+3.*tk*uk+3.*s*uk+s*tk)
	  + 2.*mZ2*mZ2*(uk+tk+2.*s)
	  )/q1h/q1/q2/s*tk
	+ ( tk*(uk*uk+tk*uk-s*s)+mZ2*(4.*s*uk-3.*tk*uk-tk*tk+4.*s*s) 
	  )/q1h/q2/s
	- ( tk*uk+s*uk-s*tk-s*s+2.*mZ2*(s-tk) ) /q1h/q1/s*tk
	+ q2*(tk*uk-s*uk-2.*s*tk-2.*s*s)/q1/q2h/s
	+ 2.*(tk*uk-tk*tk-s*tk-s*s+mZ2*(2.*s-uk))/q1/s
	- 2.*mZ2*(uk*uk-2.*mZ2*uk+2.*mZ2*mZ2)/q1/q1/q2/s*tk
	+ (2.*s*uk+tk*tk+3.*s*tk+2*s*s)/q1h/s
	+ q1*(uk+s)*(uk+tk)/q1h/q2h/s
	+ (tk*uk+s*uk+3.*s*tk+2.*s*s-mZ2*(uk+tk+2.*s))/q1h/q2h/s*uk
	+ (uk-tk)/2./q1h/q2h/s*(q1*(uk+s)/q2/tk-q2*(tk+s)/q1/uk)*tk*uk
	+ (tk-2.*mZ2)*(uk-2.*mZ2)/q1h/q1/q2h/q2*tk*uk
	- (q1*q1+q2*q2)/q1/q2
	- 2.*mZ2*(q2-2.*mZ2)/q1/q1/s*tk
      );

  swap(tk ,uk );
  swap(q1 ,q2 );
  swap(q1h,q2h);
  t_u_qqb += (2.*s)*sqr(gZ*2.)*4.*pi*alphaS_*CF_/NC_/2.
    * ( - ( tk*uk*uk+2.*s*uk*uk-tk*tk*uk
	  - 2.*s*tk*uk+mZ2*(tk*tk-uk*uk+2.*s*uk-2.*s*tk-2.*s*s)
	  )/q1h/q1/q2h/s*tk
	+ 2.*(tk*uk*uk-mZ2*uk*(s+3.*tk)+mZ2*mZ2*(2.*uk-s))/q1/q2/s
	+ ( tk*uk*(uk+s)-mZ2*(uk*uk+3.*tk*uk+3.*s*uk+s*tk)
	  + 2.*mZ2*mZ2*(uk+tk+2.*s)
	  )/q1h/q1/q2/s*tk
	+ ( tk*(uk*uk+tk*uk-s*s)+mZ2*(4.*s*uk-3.*tk*uk-tk*tk+4.*s*s) 
	  )/q1h/q2/s
	- ( tk*uk+s*uk-s*tk-s*s+2.*mZ2*(s-tk) ) /q1h/q1/s*tk
	+ q2*(tk*uk-s*uk-2.*s*tk-2.*s*s)/q1/q2h/s
	+ 2.*(tk*uk-tk*tk-s*tk-s*s+mZ2*(2.*s-uk))/q1/s
	- 2.*mZ2*(uk*uk-2.*mZ2*uk+2.*mZ2*mZ2)/q1/q1/q2/s*tk
	+ (2.*s*uk+tk*tk+3.*s*tk+2*s*s)/q1h/s
	+ q1*(uk+s)*(uk+tk)/q1h/q2h/s
	+ (tk*uk+s*uk+3.*s*tk+2.*s*s-mZ2*(uk+tk+2.*s))/q1h/q2h/s*uk
	+ (uk-tk)/2./q1h/q2h/s*(q1*(uk+s)/q2/tk-q2*(tk+s)/q1/uk)*tk*uk
	+ (tk-2.*mZ2)*(uk-2.*mZ2)/q1h/q1/q2h/q2*tk*uk
	- (q1*q1+q2*q2)/q1/q2
	- 2.*mZ2*(q2-2.*mZ2)/q1/q1/s*tk
      );
  swap(tk ,uk );
  swap(q1 ,q2 );
  swap(q1h,q2h);

  swap(q1 ,q1h);
  swap(q2 ,q2h);
  t_u_qqb += (2.*s)*sqr(gZ*2.)*4.*pi*alphaS_*CF_/NC_/2.
    * ( - ( tk*uk*uk+2.*s*uk*uk-tk*tk*uk
	  - 2.*s*tk*uk+mZ2*(tk*tk-uk*uk+2.*s*uk-2.*s*tk-2.*s*s)
	  )/q1h/q1/q2h/s*tk
	+ 2.*(tk*uk*uk-mZ2*uk*(s+3.*tk)+mZ2*mZ2*(2.*uk-s))/q1/q2/s
	+ ( tk*uk*(uk+s)-mZ2*(uk*uk+3.*tk*uk+3.*s*uk+s*tk)
	  + 2.*mZ2*mZ2*(uk+tk+2.*s)
	  )/q1h/q1/q2/s*tk
	+ ( tk*(uk*uk+tk*uk-s*s)+mZ2*(4.*s*uk-3.*tk*uk-tk*tk+4.*s*s) 
	  )/q1h/q2/s
	- ( tk*uk+s*uk-s*tk-s*s+2.*mZ2*(s-tk) ) /q1h/q1/s*tk
	+ q2*(tk*uk-s*uk-2.*s*tk-2.*s*s)/q1/q2h/s
	+ 2.*(tk*uk-tk*tk-s*tk-s*s+mZ2*(2.*s-uk))/q1/s
	- 2.*mZ2*(uk*uk-2.*mZ2*uk+2.*mZ2*mZ2)/q1/q1/q2/s*tk
	+ (2.*s*uk+tk*tk+3.*s*tk+2*s*s)/q1h/s
	+ q1*(uk+s)*(uk+tk)/q1h/q2h/s
	+ (tk*uk+s*uk+3.*s*tk+2.*s*s-mZ2*(uk+tk+2.*s))/q1h/q2h/s*uk
	+ (uk-tk)/2./q1h/q2h/s*(q1*(uk+s)/q2/tk-q2*(tk+s)/q1/uk)*tk*uk
	+ (tk-2.*mZ2)*(uk-2.*mZ2)/q1h/q1/q2h/q2*tk*uk
	- (q1*q1+q2*q2)/q1/q2
	- 2.*mZ2*(q2-2.*mZ2)/q1/q1/s*tk
      );
  swap(q1 ,q1h);
  swap(q2 ,q2h);

  swap(tk ,uk );
  swap(q1 ,q2h);
  swap(q2 ,q1h);
  t_u_qqb += (2.*s)*sqr(gZ*2.)*4.*pi*alphaS_*CF_/NC_/2.
    * ( - ( tk*uk*uk+2.*s*uk*uk-tk*tk*uk
	  - 2.*s*tk*uk+mZ2*(tk*tk-uk*uk+2.*s*uk-2.*s*tk-2.*s*s)
	  )/q1h/q1/q2h/s*tk
	+ 2.*(tk*uk*uk-mZ2*uk*(s+3.*tk)+mZ2*mZ2*(2.*uk-s))/q1/q2/s
	+ ( tk*uk*(uk+s)-mZ2*(uk*uk+3.*tk*uk+3.*s*uk+s*tk)
	  + 2.*mZ2*mZ2*(uk+tk+2.*s)
	  )/q1h/q1/q2/s*tk
	+ ( tk*(uk*uk+tk*uk-s*s)+mZ2*(4.*s*uk-3.*tk*uk-tk*tk+4.*s*s) 
	  )/q1h/q2/s
	- ( tk*uk+s*uk-s*tk-s*s+2.*mZ2*(s-tk) ) /q1h/q1/s*tk
	+ q2*(tk*uk-s*uk-2.*s*tk-2.*s*s)/q1/q2h/s
	+ 2.*(tk*uk-tk*tk-s*tk-s*s+mZ2*(2.*s-uk))/q1/s
	- 2.*mZ2*(uk*uk-2.*mZ2*uk+2.*mZ2*mZ2)/q1/q1/q2/s*tk
	+ (2.*s*uk+tk*tk+3.*s*tk+2*s*s)/q1h/s
	+ q1*(uk+s)*(uk+tk)/q1h/q2h/s
	+ (tk*uk+s*uk+3.*s*tk+2.*s*s-mZ2*(uk+tk+2.*s))/q1h/q2h/s*uk
	+ (uk-tk)/2./q1h/q2h/s*(q1*(uk+s)/q2/tk-q2*(tk+s)/q1/uk)*tk*uk
	+ (tk-2.*mZ2)*(uk-2.*mZ2)/q1h/q1/q2h/q2*tk*uk
	- (q1*q1+q2*q2)/q1/q2
	- 2.*mZ2*(q2-2.*mZ2)/q1/q1/s*tk
      );
  swap(tk ,uk );
  swap(q1 ,q2h);
  swap(q2 ,q1h);
  
  return t_u_qqb;
}

/***************************************************************************/
// M_B_WW is the Born matrix element exactly as defined in Eqs. 3.2-3.8
// of of NPB 410(1993)280-384.
double MEPP2VVPowheg::M_Born_WW(born2to2Kinematics B) const {
  Energy2 s(B.sb());
  Energy2 t(B.tb());
  Energy2 u(B.ub());
  Energy2 mW2(B.k12b()); // N.B. the diboson masses are preserved in getting
  Energy2 mZ2(B.k22b()); // the 2->2 from the 2->3 kinematics.
    
  bool up_type = abs(ab_->id())%2==0 ? true : false;
  double Qi    = up_type ? 2./3.    : -1./3. ; 
  double giL   = up_type ? guL_/2.  : gdL_/2.; 
  double giR   = up_type ? guR_/2.  : gdR_/2.; 
  double e2    = sqr(gW_)*sin2ThetaW_;
  
  double cos2ThetaW(1.-sin2ThetaW_);

  double ctt_i(gW_*gW_*gW_*gW_/16.);
  InvEnergy2 cts_i(gW_*gW_*e2/4./s *(Qi+2.*eZ_*giL/e2*s/(s-mW2/cos2ThetaW)));
  InvEnergy4 css_i(e2*e2/s/s*(sqr(Qi+eZ_*(giL+giR)/e2*s/(s-mW2/cos2ThetaW))
		      	     +sqr(   eZ_*(giL-giR)/e2*s/(s-mW2/cos2ThetaW)))
                  );

  if(!MEPP2VV::mixingInWW()&&ab_->id()!=-bb_->id()) {
    return 0.;
  }

  if(MEPP2VV::mixingInWW()) {
    ctt_i *= 8.*Fij2_/gW_/gW_;
    cts_i *= sqrt(8.*Fij2_/gW_/gW_);
    if(ab_->id()!=-bb_->id()) {
      cts_i = 0./GeV2;
      css_i = 0./GeV2/GeV2;
    }
  }

  if(!up_type) swap(t,u);
  double signf = up_type ? 1. : -1.;

  return 1./4./NC_ 
      * ( 
	  ctt_i*( 16.*(u*t/mW2/mW2-1.)*(1./4.+mW2*mW2/t/t)+16.*s/mW2)
        - cts_i*( 16.*(u*t/mW2/mW2-1.)*(s/4.-mW2/2.-mW2*mW2/t)
		+ 16.*s*(s/mW2-2.+2.*mW2/t)
	        )
	       *signf
        + 
          css_i*(  8.*(u*t/mW2/mW2-1.)*(s*s/4.-s*mW2+3.*mW2*mW2)
 		+  8.*s*s*(s/mW2-4.)
 	        )
	);
}

/***************************************************************************/
// M_V_regular_WW is the regular part of the one-loop WW matrix element 
// exactly as defined in Eqs. C.1 - C.7 of of NPB 410(1993)280-324 ***
// modulo a factor 1/(2s) ***, which is a flux factor that those authors 
// absorb in the matrix element. 
double MEPP2VVPowheg::M_V_regular_WW(real2to3Kinematics S) const {
  Energy2 s(S.bornVariables().sb());
  Energy2 t(S.bornVariables().tb());
  Energy2 u(S.bornVariables().ub());
  Energy2 mW2(S.k12r()); // N.B. the diboson masses are preserved in getting
  Energy2 mZ2(S.k22r()); // the 2->2 from the 2->3 kinematics.
  double  beta(S.betaxr()); // N.B. for x=1 \beta_x=\beta in NPB 383(1992)3-44.
    
  bool up_type = abs(ab_->id())%2==0 ? true : false;
  double Qi    = up_type ? 2./3. : -1./3.; 
  double giL   = up_type ? guL_/2.  : gdL_/2.; 
  double giR   = up_type ? guR_/2.  : gdR_/2.; 
  double e2    = sqr(gW_)*sin2ThetaW_;
  
  double cos2ThetaW(1.-sin2ThetaW_);

  double ctt_i(gW_*gW_*gW_*gW_/16.);
  InvEnergy2 cts_i(gW_*gW_*e2/4./s *(Qi+2.*eZ_*giL/e2*s/(s-mW2/cos2ThetaW)));
  InvEnergy4 css_i(e2*e2/s/s*(sqr(Qi+eZ_*(giL+giR)/e2*s/(s-mW2/cos2ThetaW))
		      	     +sqr(   eZ_*(giL-giR)/e2*s/(s-mW2/cos2ThetaW)))
                  );

  if(!MEPP2VV::mixingInWW()&&ab_->id()!=-bb_->id()) {
    return 0.;
  }

  if(MEPP2VV::mixingInWW()) {
    ctt_i *= 8.*Fij2_/gW_/gW_;
    cts_i *= sqrt(8.*Fij2_/gW_/gW_);
    if(ab_->id()!=-bb_->id()) {
      cts_i = 0./GeV2;
      css_i = 0./GeV2/GeV2;
    }
  }

  if(!up_type) swap(t,u);
  double signf = up_type ? 1. : -1.;

  InvEnergy4 TildeI4  = ( 2.*sqr(log(-t/mW2))-4.*log((mW2-t)/mW2)*log(-t/mW2)
			- 4.*ReLi2(t/mW2) )/s/t;
  InvEnergy2 TildeI3t = 1./(mW2-t)
                          *(sqr(log(mW2/s))/2.-sqr(log(-t/s))/2.-pi*pi/2.);
  InvEnergy2 TildeI3l = 1./s/beta*( 4.*ReLi2((beta-1.)/(beta+1.))
				  + sqr(log((1.-beta)/(1.+beta)))
				  + pi*pi/3.);
  double Fup1_st(0.);
  Fup1_st = 4.*(80.*t*t+73.*s*t-140.*mW2*t+72.*mW2*mW2)/t/t
          - 4.*sqr(4.*t+s)/s/beta/beta/t
          - 128.*(t+2.*s)/mW2
          + 64.*t*(t+s)/mW2/mW2
          - (32.*(t*t-3.*s*t-3.*mW2*mW2)/t/t+128.*s/(t-mW2))*log(-t/mW2)
          + ( 8.*(6.*t*t+8.*s*t-19.*mW2*t+12.*mW2*mW2)/t/t
	    - (32.*t*t-128.*s*t-26.*s*s)/s/beta/beta/t
	    + 6.*sqr(4.*t+s)/s/sqr(sqr(beta))/t
	    )*log(s/mW2)
          + 32.*s*(2.*mW2*mW2/t-u)*TildeI4
          - 64.*(t-mW2)*(2.*mW2*mW2/t/t-u/t)*TildeI3t
          + ( (16.*t*(4.*mW2-u)-49.*s*s+72.*mW2*s-48.*mW2*mW2)/2./t
	    + 2.*(8.*t*t-14.*s*t-3.*s*s)/beta/beta/t
	    - 3.*sqr(4.*t+s)/2./sqr(sqr(beta))/t
	    )*TildeI3l
          + 32./3.*( 2.*(t+2.*s)/mW2 
		   - (3.*t+2.*s-4.*mW2)/t
		   - t*(t+s)/mW2/mW2
	           )*pi*pi;
  Energy2 Jup1_st(0.*GeV2);
  Jup1_st = -128.*(t*t+2.*s*t+2.*s*s)/mW2 
          - 16.*(t*t-21.*s*t-26.*mW2*t+34.*mW2*s+17.*mW2*mW2)/t
          + 64.*s*t*(t+s)/mW2/mW2 +32.*s*s/(t-mW2)
          + ( 16.*(t-5.*s+2.*mW2)-48.*mW2*(2.*s+mW2)/t
	    + 64.*s*(2.*t+s)/(t-mW2) - 32.*s*s*t/sqr(t-mW2)
	    )*log(-t/mW2)
          + ( 16.*(4.*t+s)/beta/beta 
	    - 16.*(3.*t-2.*s)
	    + 48.*mW2*(2.*t-2.*s-mW2)/t
	    )*log(s/mW2)
          + 16.*s*(t*(2.*s+u)-2.*mW2*(2.*s+mW2))*TildeI4
          + 32.*(t-mW2)*(2.*mW2*(2.*s+mW2)/t-2.*s-u)*TildeI3t
          + ( 32.*s*t-12.*s*s+32.*mW2*mW2
	    - 16.*mW2*(2.*t+7.*s)-4.*s*(4.*t+s)/beta/beta
	    )*TildeI3l
          + 32./3.*( 2.*(t*t+2.*s*t+2.*s*s)/mW2
                   - s*t*(t+s)/mW2/mW2-2.*mW2*(2.*t-2.*s-mW2)/t-t-4.*s
	           )*pi*pi;
  Energy4 Kup1_st(0.*GeV2*GeV2);
  Kup1_st = 16.*( 12.*t*t+20.*s*t-24.*mW2*t+17.*s*s-4.*mW2*s+12.*mW2*mW2
		+ s*s*t*(t+s)/mW2/mW2-2.*s*(2.*t*t+3.*s*t+2.*s*s)/mW2)
                       *(2.-pi*pi/3.);

  return pi*alphaS_*CF_/NC_/(sqr(4.*pi)) 
      * ( ctt_i*Fup1_st - cts_i*Jup1_st*signf + css_i*Kup1_st );
}

/***************************************************************************/
// t_u_M_R_qqb is the real emission q + qb -> n + g matrix element 
// exactly as defined in Eqs. C.1 of NPB 383(1992)3-44, multiplied by
// tk * uk!
Energy2 MEPP2VVPowheg::t_u_M_R_qqb_WW(real2to3Kinematics R) const {
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

  bool up_type = abs(ab_->id())%2==0 ? true : false;
  double Qi    = up_type ? 2./3. : -1./3.; 
  double giL   = up_type ? guL_/2.  : gdL_/2.; 
  double giR   = up_type ? guR_/2.  : gdR_/2.; 
  double e2    = sqr(gW_)*sin2ThetaW_;
  
  double cos2ThetaW(1.-sin2ThetaW_);

  double ctt_i(gW_*gW_*gW_*gW_/16.);
  InvEnergy2 cts_i(gW_*gW_*e2/4./s2*(Qi+2.*eZ_*giL/e2*s2/(s2-mW2/cos2ThetaW)));
  InvEnergy4 css_i(e2*e2/s2/s2*(sqr(Qi+eZ_*(giL+giR)/e2*s2/(s2-mW2/cos2ThetaW))
		      	       +sqr(   eZ_*(giL-giR)/e2*s2/(s2-mW2/cos2ThetaW)))
                  );

  if(!MEPP2VV::mixingInWW()&&ab_->id()!=-bb_->id()) {
    return 0.*GeV2;
  }

  if(MEPP2VV::mixingInWW()) {
    ctt_i *= 8.*Fij2_/gW_/gW_;
    cts_i *= sqrt(8.*Fij2_/gW_/gW_);
    if(ab_->id()!=-bb_->id()) {
      cts_i = 0./GeV2;
      css_i = 0./GeV2/GeV2;
    }
  }

  if(!up_type) { 
    swap(q1,q1h);
    swap(q2,q2h);
  }
  double signf = up_type ? 1. : -1.;

  Energy2 t_u_Xup(0.*GeV2);
  Energy4 t_u_Yup(0.*GeV2*GeV2);
  Energy6 t_u_Zup(0.*GeV2*GeV2*GeV2);

  t_u_Xup  = 32.*mW2*(tk*uk+3.*q2*uk+q2*s+q1*q2)/q1/q2/q2*tk
           + 32.*mW2*q1/q2/q2*uk
           - 64.*mW2*s/q2
           - 32.*tk*(uk-q2)/q1/q2*tk
           + 64.*mW2*mW2*mW2/q1/q1/q2*tk
           - 16.*(2.*tk-2.*s-q2)/q2*uk
           + 16.*s*(2.*s+2.*q1+q2/2.)/q2
           - 8.*(4.*tk+uk+9.*s+2.*q2+2.*q1)/mW2*tk
           - 16.*s*(2.*s+q1)/mW2
           - 64.*mW2*mW2*(tk*uk+q2*tk+q1*uk-q2*s/2.)/q1/q2/q2
           + 8.*s2*q1*(tk+s+q1)/mW2/mW2;
  swap(tk,uk);
  swap(q1,q2);
  t_u_Xup += 32.*mW2*(tk*uk+3.*q2*uk+q2*s+q1*q2)/q1/q2/q2*tk
           + 32.*mW2*q1/q2/q2*uk
           - 64.*mW2*s/q2
           - 32.*tk*(uk-q2)/q1/q2*tk
           + 64.*mW2*mW2*mW2/q1/q1/q2*tk
           - 16.*(2.*tk-2.*s-q2)/q2*uk
           + 16.*s*(2.*s+2.*q1+q2/2.)/q2
           - 8.*(4.*tk+uk+9.*s+2.*q2+2.*q1)/mW2*tk
           - 16.*s*(2.*s+q1)/mW2
           - 64.*mW2*mW2*(tk*uk+q2*tk+q1*uk-q2*s/2.)/q1/q2/q2
           + 8.*s2*q1*(tk+s+q1)/mW2/mW2;
  swap(tk,uk);
  swap(q1,q2);

  t_u_Yup  = - 16.*tk*(uk*(uk+s+q1)+q2*(s-2.*q1))/q1/q2*tk
             - 32.*mW2*mW2*s/q2
             - 32.*mW2*mW2*mW2/q1/q2*tk
             + 16.*(2.*q2*uk+s*s+q1*s+5.*q2*s+q1*q2+2.*q2*q2)/q2*tk
             - 16.*(q2*q2+s*s-q2*s)/q1*tk
             + 16.*s*(q1*s+3./2.*q2*s+q1*q2-q1*q1)/q2
             + 16.*mW2*tk*(4.*uk+s+q1-2.*q2)/q1/q2*tk
             + 16.*mW2*(3.*s*uk+q1*uk-q1*s-3.*q2*s-q1*q1+q2*q2)/q1/q2*tk
             + 16.*mW2*s*(q2-4.*s+2.*q1)/q2
             - 8.*s2*(4.*tk+uk+9.*s+4.*q1+2.*q2)/mW2*tk
             - 16.*s2*(2.*s*s+2.*q1*s+q1*q1)/mW2
             - 32.*mW2*mW2*(tk+uk/2.+2.*s-q1)/q1/q2*tk
             + 8.*s2*s2*q1*(tk+s+q1)/mW2/mW2;
  swap(tk,uk);
  swap(q1,q2);
  t_u_Yup += - 16.*tk*(uk*(uk+s+q1)+q2*(s-2.*q1))/q1/q2*tk
             - 32.*mW2*mW2*s/q2
             - 32.*mW2*mW2*mW2/q1/q2*tk
             + 16.*(2.*q2*uk+s*s+q1*s+5.*q2*s+q1*q2+2.*q2*q2)/q2*tk
             - 16.*(q2*q2+s*s-q2*s)/q1*tk
             + 16.*s*(q1*s+3./2.*q2*s+q1*q2-q1*q1)/q2
             + 16.*mW2*tk*(4.*uk+s+q1-2.*q2)/q1/q2*tk
             + 16.*mW2*(3.*s*uk+q1*uk-q1*s-3.*q2*s-q1*q1+q2*q2)/q1/q2*tk
             + 16.*mW2*s*(q2-4.*s+2.*q1)/q2
             - 8.*s2*(4.*tk+uk+9.*s+4.*q1+2.*q2)/mW2*tk
             - 16.*s2*(2.*s*s+2.*q1*s+q1*q1)/mW2
             - 32.*mW2*mW2*(tk+uk/2.+2.*s-q1)/q1/q2*tk
             + 8.*s2*s2*q1*(tk+s+q1)/mW2/mW2;
  swap(tk,uk);
  swap(q1,q2);

  t_u_Zup  =   8.*s2*(9.*tk+3.*uk+20.*s+10.*q1+4.*q2)*tk
             + 8.*s2*(17./2.*s*s+10.*q1*s+6.*q1*q1)
             - 4.*s2*s2*(4.*tk+uk+9.*s+6.*q1+2.*q2)/mW2*tk
             - 8.*s2*s2*(2.*s*s+3.*q1*s+2.*q1*q1)/mW2
             - 16.*mW2*(2.*tk+5.*uk+7.*s+6.*q1+6.*q2)*tk
             - 16.*mW2*s*(s+6.*q1)
             + 4.*s2*s2*s2*q1*(tk+s+q1)/mW2/mW2
             + 48.*mW2*mW2*s2;
  swap(tk,uk);
  swap(q1,q2);
  t_u_Zup +=   8.*s2*(9.*tk+3.*uk+20.*s+10.*q1+4.*q2)*tk
             + 8.*s2*(17./2.*s*s+10.*q1*s+6.*q1*q1)
             - 4.*s2*s2*(4.*tk+uk+9.*s+6.*q1+2.*q2)/mW2*tk
             - 8.*s2*s2*(2.*s*s+3.*q1*s+2.*q1*q1)/mW2
             - 16.*mW2*(2.*tk+5.*uk+7.*s+6.*q1+6.*q2)*tk
             - 16.*mW2*s*(s+6.*q1)
             + 4.*s2*s2*s2*q1*(tk+s+q1)/mW2/mW2
             + 48.*mW2*mW2*s2;
  swap(tk,uk);
  swap(q1,q2);

  return -pi*alphaS_*CF_/NC_ 
      * ( ctt_i*t_u_Xup - cts_i*t_u_Yup*signf + css_i*t_u_Zup );

}

