// -*- C++ -*-
//
// MEPP2HiggsPowheg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2HiggsPowheg class.
//

#include "MEPP2HiggsPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig++/MatrixElement/General/HardVertex.h"

using namespace Herwig;

MEPP2HiggsPowheg::MEPP2HiggsPowheg() : 
  contrib_(1) ,  nlo_alphaS_opt_(0) , fixed_alphaS_(0.11803463),
  widthopt_(1),  usersWidth_(0.00468456293*GeV)         ,
  scaleopt_(1),  mu_F_(100.*GeV)    , scaleFact_(1.)    ,
  shapeopt_(2),  processopt_(1)     ,  minflavouropt_(4), maxflavouropt_(5), 
  mh_(0.*GeV) ,  wh_(0.*GeV)
{}

ClassDescription<MEPP2HiggsPowheg> MEPP2HiggsPowheg::initMEPP2HiggsPowheg;
// Definition of the static class description member.

void MEPP2HiggsPowheg::persistentOutput(PersistentOStream & os) const {
  os << hggvertex      << ffhvertex        << theSM 
     << contrib_       << nlo_alphaS_opt_  << fixed_alphaS_ 
     << widthopt_      << ounit(usersWidth_,GeV) 
     << scaleopt_      << ounit(mu_F_,GeV) << scaleFact_     << shapeopt_      
     << processopt_    << minflavouropt_   << maxflavouropt_ << hmass_        
     << ounit(mh_,GeV) << ounit(wh_,GeV);
}

void MEPP2HiggsPowheg::persistentInput(PersistentIStream & is, int) {
  is >> hggvertex      >> ffhvertex        >> theSM 
     >> contrib_       >> nlo_alphaS_opt_  >> fixed_alphaS_ 
     >> widthopt_      >> iunit(usersWidth_,GeV) 
     >> scaleopt_      >> iunit(mu_F_,GeV) >> scaleFact_     >> shapeopt_      
     >> processopt_    >> minflavouropt_   >> maxflavouropt_ >> hmass_     
     >> iunit(mh_,GeV) >> iunit(wh_,GeV);
}

void MEPP2HiggsPowheg::Init() {

  static ClassDocumentation<MEPP2HiggsPowheg> documentation
    ("The MEPP2HiggsPowheg class implements the matrix elements for"
     " Higgs production (with decay H->W-W+) in hadron-hadron collisions.");

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

  static Switch<MEPP2HiggsPowheg,unsigned int> interfaceWidthOption
    ("WidthOption",
     "Option to allow user to specify the width",
     &MEPP2HiggsPowheg::widthopt_, 1, false, false);
  static SwitchOption interfaceWidthGenerator
    (interfaceWidthOption,
     "WidthGenerator",
     "Herwig++/ThePEG calculates the width (default)",
     1);
  static SwitchOption interfaceSpecifyWidth
    (interfaceWidthOption,
     "SpecifyWidth",
     "The user can set the width manually setting interface MyHiggsWidth",
     2);

  static Parameter<MEPP2HiggsPowheg,Energy> interfaceMyHiggsWidth
    ("MyHiggsWidth",
     "Value to use when overriding widthGenerator with interface WidthOption",
     &MEPP2HiggsPowheg::usersWidth_, GeV, 0.00468456293*GeV, 0.0*GeV, 10.0*GeV,
     true, false, Interface::limited);

  static Switch<MEPP2HiggsPowheg,unsigned int> interfaceFactorizationScaleOption
    ("FactorizationScaleOption",
     "Option for the choice of factorization scale",
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

  static Parameter<MEPP2HiggsPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before sHat if using a running scale",
     &MEPP2HiggsPowheg::scaleFact_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Switch<MEPP2HiggsPowheg,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &MEPP2HiggsPowheg::shapeopt_, 1, false, false);
  static SwitchOption interfaceStandardShapeFixed
    (interfaceShapeOption,
     "FixedBreitWigner",
     "Breit-Wigner s-channel resonanse",
     1);
  static SwitchOption interfaceStandardShapeRunning
    (interfaceShapeOption,
     "MassGenerator",
     "Use the mass generator to give the shape",
     2);

  static Switch<MEPP2HiggsPowheg,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2HiggsPowheg::processopt_, 1, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     1);
  static SwitchOption interfaceProcess1
    (interfaceProcess,
     "qqbar",
     "Only include the incoming q qbar subprocess",
     2);
  static SwitchOption interfaceProcessgg
    (interfaceProcess,
     "gg",
     "Only include the incoming gg subprocess",
     3);

  static Parameter<MEPP2HiggsPowheg,unsigned int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The minimum flavour of the incoming quarks in the hard process",
     &MEPP2HiggsPowheg::minflavouropt_, 4, 3, 5,
     false, false, Interface::limited);

  static Parameter<MEPP2HiggsPowheg,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the incoming quarks in the hard process",
     &MEPP2HiggsPowheg::maxflavouropt_, 5, 3, 5,
     false, false, Interface::limited);
}

void MEPP2HiggsPowheg::doinit() throw(InitException) {
  // Set colour factors:
  CA_    = 3. ;  CF_  = 4./3. ;  TR_  = 0.5 ;  nlf_ = 5. ;
  beta0_ = (11.*CA_/3. - 4.*TR_*nlf_/3.)/(4.*Constants::pi);
  // MEPP2HiggsPowheg.cc::doinit() resumes here:
  MEBase::doinit();
  // get the vertex pointers from the SM object
  theSM = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!theSM) {
    throw InitException() << "Wrong type of StandardModel object in MEPP2HiggsPowheg::doinit(),"
                          << " the Herwig++ version must be used" 
                          << Exception::runerror;
  }
  hggvertex = theSM->vertexHGG();
  ffhvertex = theSM->vertexFFH();
  // get the mass generator for the higgs
  PDPtr h0 = getParticleData(ParticleID::h0);
  mh_ = h0->mass();
  wh_ = widthopt_==2 ? usersWidth_ : h0->generateWidth(mh_);
  if(h0->massGenerator()) {
    hmass_=dynamic_ptr_cast<SMHiggsMassGeneratorPtr>(h0->massGenerator());
  }
  if(shapeopt_==2&&!hmass_) throw InitException()
    << "If using the mass generator for the line shape in MEPP2HiggsPowheg::doinit()"
    << "the mass generator must be an instance of the SMHiggsMassGenerator class"
    << Exception::runerror;
}

unsigned int MEPP2HiggsPowheg::orderInAlphaS() const {
  return 2;
}

unsigned int MEPP2HiggsPowheg::orderInAlphaEW() const {
  return 1;
}

Energy2 MEPP2HiggsPowheg::scale() const {
  return scaleopt_ == 1 ?  scaleFact_*sHat() : sqr(mu_F_);
}

int MEPP2HiggsPowheg::nDim() const {
  return 2;
}

bool MEPP2HiggsPowheg::generateKinematics(const double * r) {
  xt_=*(r  );
  y_ =*(r+1) * 2. - 1.;
  Lorentz5Momentum pout = meMomenta()[0] + meMomenta()[1];
  pout.rescaleMass();
  meMomenta()[2].setMass(pout.mass());
  meMomenta()[2] = LorentzMomentum(pout.x(),pout.y(),pout.z(),pout.t());
  jacobian(1.0);
  // check whether it passes all the cuts: returns true if it does
  vector<LorentzMomentum> out(1,meMomenta()[2]);
  tcPDVector tout(1,mePartonData()[2]);
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}

void MEPP2HiggsPowheg::getDiagrams() const {
  tcPDPtr h0=getParticleData(ParticleID::h0);
  // gg -> H process
  if(processopt_==1||processopt_==3) {
    tcPDPtr g=getParticleData(ParticleID::g);
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, h0, -1)));
  }
  // q qbar -> H processes
  if(processopt_==1||processopt_==2) {
    for (unsigned int i = minflavouropt_; i <= maxflavouropt_; ++i) {
      tcPDPtr q = getParticleData(i);
      tcPDPtr qb = q->CC();
    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, h0, -2)));
    }
  }
}

CrossSection MEPP2HiggsPowheg::dSigHatDR() const {
  // Get Born momentum fractions xbp_ and xbm_:
  get_born_variables();
  using Constants::pi;
  InvEnergy2 bwfact;
  if(widthopt_==2) {
    bwfact = wh_*sqrt(sHat())/pi/(sqr(sHat()-sqr(mh_))+sqr(mh_*wh_));
    double cs = me2()*jacobian()*pi*double(UnitRemoval::E4 * bwfact/sHat());
    return UnitRemoval::InvE2 * sqr(hbarc) * cs;
  }  
  if(shapeopt_==1) {
    bwfact = mePartonData()[2]->generateWidth(sqrt(sHat()))*sqrt(sHat())/pi/
      (sqr(sHat()-sqr(mh_))+sqr(mh_*wh_));
  } else {
    bwfact = hmass_->BreitWignerWeight(sqrt(sHat()),0);
  }
  double cs = me2() * jacobian() * pi * double(UnitRemoval::E4 * bwfact/sHat());
  return UnitRemoval::InvE2 * sqr(hbarc) * cs;
}

double MEPP2HiggsPowheg::me2() const {
  double output(0.0);
  useMe();
  ScalarWaveFunction hout(meMomenta()[2],mePartonData()[2],outgoing);

// Safety code to garantee the reliable behaviour of Higgs shape limits 
// (important for heavy and broad Higgs resonance).
  Energy hmass = meMomenta()[2].m();
  tcPDPtr h0 = mePartonData()[2];
  Energy mass = h0->mass();
  Energy halfmass = .5*mass;
  if (.0*GeV > hmass) return 0.0;
  // stricly speaking the condition is applicable if h0->widthUpCut() == h0->widthLoCut()...
  if (h0->widthLoCut() > halfmass) {
    if ((mass + h0->widthUpCut() < hmass || mass - h0->widthLoCut() > hmass)) return 0.0;
  } else {
    if (mass + halfmass < hmass || halfmass > hmass) return 0.0;
  }

  if (mePartonData()[0]->id() == ParticleID::g && 
      mePartonData()[1]->id() == ParticleID::g) {
    VectorWaveFunction gin1(meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction gin2(meMomenta()[1],mePartonData()[1],incoming);

    vector<VectorWaveFunction> g1,g2;
    for(unsigned int i = 0; i < 2; ++i) {
      gin1.reset(2*i);
      g1.push_back(gin1);
      gin2.reset(2*i);
      g2.push_back(gin2);
    }
    output = ggME(g1,g2,hout,false);
  } else {
    if (mePartonData()[0]->id() == -mePartonData()[1]->id()) {
      SpinorWaveFunction    qin (meMomenta()[0],mePartonData()[0],incoming);
      SpinorBarWaveFunction qbin(meMomenta()[1],mePartonData()[1],incoming);

      vector<SpinorWaveFunction> fin;
      vector<SpinorBarWaveFunction> ain;
      for (unsigned int i = 0; i < 2; ++i) {
        qin.reset(i);
        fin.push_back(qin);
        qbin.reset(i);
        ain.push_back(qbin);
      }
      output = qqME(fin,ain,hout,false);
    }
    else {
    throw Exception() << "Unknown subprocess in MEPP2HiggsPowheg::me2()" 
		      << Exception::runerror;
    }
  }
  return output;
}

Selector<MEBase::DiagramIndex> MEPP2HiggsPowheg::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for (DiagramIndex i = 0; i < diags.size(); ++i)
    sel.insert(1.0, i);
  return sel;
}

Selector<const ColourLines *> MEPP2HiggsPowheg::colourGeometries(tcDiagPtr diag) const {
  // colour lines
  static const ColourLines line1("1 -2,2 -1");
  static const ColourLines line2("1 -2");
  // select the colour flow
  Selector<const ColourLines *> sel;
  if (diag->id() == -1) {
    sel.insert(1.0, &line1);
  } else {
    sel.insert(1.0, &line2);
  }
  // return the answer
  return sel;
}

void MEPP2HiggsPowheg::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  if(hard[0]->id() < hard[1]->id()) {
    swap(hard[0],hard[1]);
  }
  // identify the process and calculate the matrix element
  if(hard[0]->id() == ParticleID::g && hard[1]->id() == ParticleID::g) {
    vector<VectorWaveFunction> g1,g2;
    vector<SpinorBarWaveFunction> q;
    vector<SpinorWaveFunction> qbar;
    VectorWaveFunction (g1,hard[0],incoming,false,true,true);
    VectorWaveFunction (g2,hard[1],incoming,false,true,true);
    ScalarWaveFunction hout(hard[2],outgoing,true,true);
    g1[1] = g1[2];
    g2[1] = g2[2];
    ggME(g1,g2,hout,true);
  } else {
    vector<SpinorWaveFunction>    q1;
    vector<SpinorBarWaveFunction> q2;
    SpinorWaveFunction    (q1,hard[0],incoming,false,true);
    SpinorBarWaveFunction (q2,hard[1],incoming,false,true);
    ScalarWaveFunction     hout(hard[2],outgoing,true,true);
    qqME(q1,q2,hout,true);
  }
  // construct the vertex
  HardVertexPtr hardvertex = new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me_);
  // set the pointers and to and from the vertex
  for(unsigned int i = 0; i < 3; ++i) {
    dynamic_ptr_cast<SpinfoPtr>(hard[i]->spinInfo())->setProductionVertex(hardvertex);
  }
}

double MEPP2HiggsPowheg::ggME(vector<VectorWaveFunction> g1, 
                          vector<VectorWaveFunction> g2, 
                          ScalarWaveFunction & in, 
                          bool calc) const {
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1,PDT::Spin0);
  Energy2 s(sHat());
  double me2(0.0);
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < 2; ++j) {
      Complex diag = hggvertex->evaluate(s,g1[i],g2[j],in);
      me2 += norm(diag);
      if(calc) newme(2*i, 2*j, 0) = diag;
    }
  }
  if(calc) me_.reset(newme);
  // initial colour and spin factors: colour -> (8/64) and spin -> (1/4)
  lo_ggME_ = me2/32.;
  return lo_ggME_;
}


double MEPP2HiggsPowheg::qqME(vector<SpinorWaveFunction> & fin, 
                          vector<SpinorBarWaveFunction> & ain, 
                          ScalarWaveFunction & in, 
                          bool calc) const {
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0);
  Energy2 s(scale());
  double me2(0.0);
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < 2; ++j) {
      Complex diag = ffhvertex->evaluate(s,fin[i],ain[j],in);
      me2+=norm(diag);
      if(calc) newme(i, j, 0) = diag;
    }
  }
  if(calc) me_.reset(newme);
  // final colour/spin factors
  return me2/12.;
}

double MEPP2HiggsPowheg::NLOweight() const {
  // If only leading order is required return 1:
  if(contrib_==0) return 1.;
  // Get particle data for QCD particles:
  a_lo_=mePartonData()[0];
  b_lo_=mePartonData()[1];
  // get BeamParticleData objects for PDF's
  hadron_A_=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().first->dataPtr());
  hadron_B_=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().second->dataPtr());
  // If necessary swap the particle data vectors so that xbp_, 
  // mePartonData[0], beam[0] relate to the inbound gluon: 
  if(!(lastPartons().first ->dataPtr()==a_lo_&&
       lastPartons().second->dataPtr()==b_lo_)) {
    swap(xbp_     ,xbm_     );
    swap(hadron_A_,hadron_B_);
  }
  // calculate the PDF's for the Born process
  oldlumi_ = hadron_A_->pdf()->xfx(hadron_A_,a_lo_,scale(),xbp_)/xbp_
           * hadron_B_->pdf()->xfx(hadron_B_,b_lo_,scale(),xbm_)/xbm_;
  // Calculate alpha_S
  alphaS_ = nlo_alphaS_opt_==1 ? fixed_alphaS_ : SM().alphaS(scale());
  // Calculate the integrand
  double wgt(0.);
//   // gg contribution
//   double wggvirt      = Vtilde_gg();
//   double wggcollin    = Ctilde_gg(x(_xt, 1.), 1.) 
//                       + Ctilde_gg(x(_xt,-1.),-1.);
//   double wggreal      = Ftilde_gg(_xt,_y);
//   double wgg          = wggvirt+wggcollin+wggreal;
//   // g q contribution
//   double wgqcollin    = Ctilde_gq(x(_xt,0.),0.);
//   double wgqreal      = Ftilde_gq(_xt,_y);
//   double wgq          = wgqreal+wgqcollin;
//   // q qbar contribution
//   double wqqreal      = Ftilde_qq(_xt,_v);
//   double wqq          = wqqreal+wqqcollin;
//   // total
//   wgt                 = 1.+(wgg+wgq+wqq);
  // return the answer
  return contrib_==1 ? max(0.,wgt) : max(0.,-wgt);
}

void MEPP2HiggsPowheg::get_born_variables() const {
  xbp_ = lastX1(); etabarp_ = sqrt(1.-xbp_);
  xbm_ = lastX2(); etabarm_ = sqrt(1.-xbm_);
  p2_    = sHat();
  s2_    = p2_;
  if(meMomenta().size()==4) {
    p12_    = meMomenta()[1].m2();
    p22_    = meMomenta()[2].m2();
    theta1_= acos(meMomenta()[2].z()/meMomenta()[2].vect().mag());
    theta2_= atan(meMomenta()[2].x()/meMomenta()[2].y());
  }
  return;
}
Energy2 MEPP2HiggsPowheg::s(double xt, double y) const {
  return  p2_/x(xt,y);
}

Energy2 MEPP2HiggsPowheg::tk(double xt, double y) const {
  double  x_xt_y(x(xt,y));
  return -0.5*p2_/x_xt_y*(1.- x_xt_y)*(1.-y);
}

Energy2 MEPP2HiggsPowheg::uk(double xt, double y) const {
  double  x_xt_y(x(xt,y));
  return -0.5*p2_/x_xt_y*(1.- x_xt_y)*(1.+y);
}

double MEPP2HiggsPowheg::betax(double xt, double y) const {
  double  x_xt_y(x(xt,y));
  Energy2 s_xt_y(s(xt,y));
  return   sqrt(1.-sqr(sqrt(p12_)+sqrt(p22_))/x_xt_y/s_xt_y)
	 * sqrt(1.-sqr(sqrt(p12_)-sqrt(p22_))/x_xt_y/s_xt_y);
}

double MEPP2HiggsPowheg::v1(double xt, double y) const {
  return   betax(xt,y)/(1.-(p22_-p12_)/x(xt,y)/s(xt,y));
}

double MEPP2HiggsPowheg::v2(double xt, double y) const {
  return   betax(xt,y)/(1.+(p22_-p12_)/x(xt,y)/s(xt,y));
}

double MEPP2HiggsPowheg::cpsi(double xt, double y) const {
  Energy2 s_xt_y(s(xt,y));
  return 1.- s_xt_y 
           / 2.
           / ((s_xt_y+tk(xt,y))/2./sqrt(s2_)*(s_xt_y +uk(xt,y))/2./sqrt(s2_));
}

double MEPP2HiggsPowheg::cpsipr(double xt, double y) const {
  Energy2 tk_xt_y(tk(xt,y));
  return 1.+ tk_xt_y
           / 2.
           / ((s(xt,y)+tk_xt_y)/2./sqrt(s2_)*-(tk_xt_y+uk(xt,y))/2./sqrt(s2_));
}

Energy2 MEPP2HiggsPowheg::q1(double xt, double y) const {
  return p12_ - 0.5*(s(xt,y)+tk(xt,y))
                   *betax(xt,y)/v1(xt,y)*(1.-v1(xt,y)*cos(theta1_));
}

Energy2 MEPP2HiggsPowheg::q2(double xt, double y) const {
  return p22_ - 0.5*(s(xt,y)+uk(xt,y))
                   *betax(xt,y)/v2(xt,y)
                   *(1.+ v2(xt,y)*cos(theta2_)*sin(theta1_)
                                *sqrt(1.-sqr(cpsi(xt,y)))
		       + v2(xt,y)*cos(theta1_)*cpsi(xt,y)
 	            );
}

Energy2 MEPP2HiggsPowheg::q1hat(double xt, double y) const {
  return p12_ + p22_ - s(xt,y)  - tk(xt,y) - q1(xt,y);
}

Energy2 MEPP2HiggsPowheg::q2hat(double xt, double y) const {
  return p12_ + p22_ - s(xt,y)  - uk(xt,y) - q2(xt,y);
}

Energy2 MEPP2HiggsPowheg::w1(double xt, double y) const {
  return p12_ - q1(xt,y)  + q2(xt,y) - tk(xt,y);
}

Energy2 MEPP2HiggsPowheg::w2(double xt, double y) const {
  return p22_ + q1(xt,y)  - q2(xt,y) - uk(xt,y);
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

double MEPP2HiggsPowheg::x(double xt, double y) const {
  double x0(xbar(y));
  return x0+(1.-x0)*xt;
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
  double newlumi(-999.);
  double xp_x_y(xp(x,y)),xm_x_y(xm(x,y));
  newlumi = (hadron_A_->pdf()->xfx(hadron_A_,a,scale(),xp_x_y)/xp_x_y)
          * (hadron_B_->pdf()->xfx(hadron_B_,b,scale(),xm_x_y)/xm_x_y);
  return newlumi / oldlumi_;
}

double MEPP2HiggsPowheg::Vtilde_universal() const {
  return  alphaS_/2./Constants::pi*CA_ 
        * ( log(p2_/sqr(mu_F_))*( 2.*(2.*Constants::pi*beta0_/CA_)
	     	          + 4.*(log(etabarp_)+log(etabarm_)))
	                  + 8.*sqr(log(etabarp_)) + 8.*sqr(log(etabarm_))
	                  - 2.*sqr(Constants::pi)/3.
	  );
}

double MEPP2HiggsPowheg::Ctilde_Ltilde_qq_on_x(tcPDPtr a, tcPDPtr b, 
					       double xt, double y ) const {
  if(y!= 1.&&y!=-1.) { cout << "\nCtilde_qq::y value not allowed."; }
  if(y== 1.&&!(abs(a->id())>0&&abs(a->id()<7))) 
    cout << "\nCtilde_qq::for Cqq^plus  a must be a quark! id = ",a->id();
  if(y==-1.&&!(abs(b->id())>0&&abs(b->id()<7))) 
    cout << "\nCtilde_qq::for Cqq^minus b must be a quark! id = ",b->id();
  double x_pm      = y == 1  ? xbp_     : xbm_     ;
  double etabar_pm = y == 1. ? etabarp_ : etabarm_ ;
  return ( ( ( (1./(1.-xt))*log(p2_/sqr(mu_F_)/x_pm)+4.*log(etabar_pm)/(1.-xt)
       	     + 2.*log(1.-xt)/(1.-xt)
             )*CF_*(1.+sqr(x_pm)) 
	   + sqr(etabar_pm)*CF_*(1.-x_pm)
	   )*Lhat_ab(a,b,x_pm,y)
	 - ( ( (1./(1.-xt))*log(p2_/sqr(mu_F_)     )+4.*log(etabar_pm)/(1.-xt)
       	     + 2.*log(1.-xt)/(1.-xt)
             )*CF_*2. 
	   )
         )/x_pm;
}

double MEPP2HiggsPowheg::Ctilde_Ltilde_gg_on_x(tcPDPtr a, tcPDPtr b, 
					       double xt, double y ) const {
  if(y!= 1.&&y!=-1.) { cout << "\nCtilde_gg::y value not allowed."; }
  if(y== 1.&&a->id()!=21) 
    cout << "\nCtilde_gg::for Cgg^plus  a must be a gluon! id = ",a->id();
  if(y==-1.&&b->id()!=21) 
    cout << "\nCtilde_gg::for Cgg^minus b must be a gluon! id = ",b->id();
  double x_pm      = y == 1  ? xbp_     : xbm_     ;
  double etabar_pm = y == 1. ? etabarp_ : etabarm_ ;
  return ( ( ( (1./(1.-xt))*log(p2_/sqr(mu_F_)/x_pm)+4.*log(etabar_pm)/(1.-xt)
       	     + 2.*log(1.-xt)/(1.-xt)
             )*2.*CA_*(x_pm+sqr(1.-x_pm)/x_pm+x_pm*sqr(1.-x_pm))

	   )*Lhat_ab(a,b,x_pm,y)
         - ( ( (1./(1.-xt))*log(p2_/sqr(mu_F_)     )+4.*log(etabar_pm)/(1.-xt)
	     + 2.*log(1.-xt)/(1.-xt)
	     )*2.*CA_
	   )
         ) / x_pm;
}

double MEPP2HiggsPowheg::Ctilde_Ltilde_qg_on_x(tcPDPtr a, tcPDPtr b, 
					       double xt, double y ) const {
  if(y!= 1.&&y!=-1.) { cout << "\nCtilde_qg::y value not allowed."; }
  if(y== 1.&&!(abs(a->id())>0&&abs(a->id()<7))) 
    cout << "\nCtilde_qg::for Cqg^plus  a must be a quark! id = ",a->id();
  if(y==-1.&&!(abs(b->id())>0&&abs(b->id()<7))) 
    cout << "\nCtilde_qg::for Cqg^minus b must be a quark! id = ",b->id();
  double x_pm      = y == 1  ? xbp_     : xbm_     ;
  double etabar_pm = y == 1. ? etabarp_ : etabarm_ ;
  return ( ( ( (1./(1.-xt))*log(p2_/sqr(mu_F_)/x_pm)+4.*log(etabar_pm)/(1.-xt)
       	     + 2.*log(1.-xt)/(1.-xt)
             )*(1.-x_pm)*CF_*(1.+sqr(1.-x_pm))/x_pm
	   + sqr(etabar_pm)*CF_*x_pm
	   )*Lhat_ab(a,b,x_pm,y)
         ) / x_pm;
}

double MEPP2HiggsPowheg::Ctilde_Ltilde_gq_on_x(tcPDPtr a, tcPDPtr b, 
					       double xt, double y ) const {
  if(y!= 1.&&y!=-1.) { cout << "\nCtilde_gq::y value not allowed."; }
  if(y== 1.&&a->id()!=21)
    cout << "\nCtilde_gq::for Cgq^plus  a must be a gluon! id = ",a->id();
  if(y==-1.&&b->id()!=21)
    cout << "\nCtilde_gq::for Cgq^minus b must be a gluon! id = ",b->id();
  double x_pm      = y == 1  ? xbp_     : xbm_     ;
  double etabar_pm = y == 1. ? etabarp_ : etabarm_ ;
  return ( ( ( (1./(1.-xt))*log(p2_/sqr(mu_F_)/x_pm)+4.*log(etabar_pm)/(1.-xt)
       	     + 2.*log(1.-xt)/(1.-xt)
             )*(1.-x_pm)*TR_*(sqr(x_pm)+sqr(1.-x_pm))
	   + sqr(etabar_pm)*TR_*2.*x_pm*(1.-x_pm)
	   )*Lhat_ab(a,b,x_pm,y)
         ) / x_pm;
}

double MEPP2HiggsPowheg::M_V_regular() const {
  return alphaS_/2./Constants::pi*CA_*
                        (  11./3.
			+  4.*sqr(Constants::pi)/3.
			- (4.*Constants::pi*beta0_/CA_)*log(p2_/sqr(mu_UV_))
			)*lo_ggME_;
}

double MEPP2HiggsPowheg::M_R_qq(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*
                        ( 0.
                        )*lo_ggME_;
}

double MEPP2HiggsPowheg::M_R_gg(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*
                        ( 0.
			)*lo_ggME_;
}

double MEPP2HiggsPowheg::M_R_qg(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*
                        ( 0.
			)*lo_ggME_;
}

double MEPP2HiggsPowheg::M_R_gq(double xt, double y) const {
  return 8.*Constants::pi*alphaS_*
                        ( 0.
			)*lo_ggME_;
}

double MEPP2HiggsPowheg::Rcal_Ltilde_qq_on_x(tcPDPtr a , tcPDPtr b,
					     double  xt, double y ) const {
  return   M_R_qq(xt ,y)*Lhat_ab(a,b,x(xt,y),y)
	 - M_R_qq( 1.,y)*Lhat_ab(a,b,     1.,y);
}

double MEPP2HiggsPowheg::Rcal_Ltilde_gg_on_x(tcPDPtr a , tcPDPtr b, 
					     double  xt, double y ) const {
  return   M_R_gg(xt ,y)*Lhat_ab(a,b,x(xt,y),y)
	 - M_R_gg( 1.,y)*Lhat_ab(a,b,     1.,y);
}

double MEPP2HiggsPowheg::Rcal_Ltilde_gq_on_x(tcPDPtr a , tcPDPtr b, 
					     double  xt, double y ) const {
  return   M_R_gq(xt ,y)*Lhat_ab(a,b,x(xt,y),y)
	 - M_R_gq( 1.,y)*Lhat_ab(a,b,     1.,y);
}

double MEPP2HiggsPowheg::Rcal_Ltilde_qg_on_x(tcPDPtr a , tcPDPtr b, 
					     double  xt, double y ) const {
  return   M_R_qg(xt ,y)*Lhat_ab(a,b,x(xt,y),y)
	 - M_R_qg( 1.,y)*Lhat_ab(a,b,     1.,y);
}




