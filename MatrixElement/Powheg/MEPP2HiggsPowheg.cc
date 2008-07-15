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
  _contrib(1) , _nlo_alphaS_opt(0) , _fixed_alphaS(0.11803463),
  _a(0.5)     , _p(0.7)            , _eps(1.0e-8)      ,
  widthopt(1) ,  usersWidth(0.00468456293*GeV)         ,
  scaleopt(1) ,  mu_F(100.*GeV)    , _scaleFact(1.)    ,
  shapeopt(2) ,  processopt(1)     ,  minflavouropt(4) , maxflavouropt(5), 
  _mh(0.*GeV) , _wh(0.*GeV)
{}

ClassDescription<MEPP2HiggsPowheg> MEPP2HiggsPowheg::initMEPP2HiggsPowheg;
// Definition of the static class description member.

void MEPP2HiggsPowheg::persistentOutput(PersistentOStream & os) const {
  os << hggvertex      << ffhvertex       << theSM 
     << _contrib       << _nlo_alphaS_opt << _fixed_alphaS 
     << _a             << _p              << _eps
     << widthopt       << ounit(usersWidth,GeV) 
     << scaleopt       << ounit(mu_F,GeV) << _scaleFact    << shapeopt      
     << processopt     << minflavouropt   << maxflavouropt << _hmass        
     << ounit(_mh,GeV) << ounit(_wh,GeV);
}

void MEPP2HiggsPowheg::persistentInput(PersistentIStream & is, int) {
  is >> hggvertex      >> ffhvertex       >> theSM 
     >> _contrib       >> _nlo_alphaS_opt >> _fixed_alphaS 
     >> _a             >> _p              >> _eps
     >> widthopt       >> iunit(usersWidth,GeV) 
     >> scaleopt       >> iunit(mu_F,GeV) >> _scaleFact    >> shapeopt      
     >> processopt     >> minflavouropt   >> maxflavouropt >> _hmass        
     >> iunit(_mh,GeV) >> iunit(_wh,GeV);
}

void MEPP2HiggsPowheg::Init() {

  static ClassDocumentation<MEPP2HiggsPowheg> documentation
    ("The MEPP2HiggsPowheg class implements the matrix elements for"
     " Higgs production (with decay H->W-W+) in hadron-hadron collisions.");

  static Switch<MEPP2HiggsPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEPP2HiggsPowheg::_contrib, 1, false, false);
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
     &MEPP2HiggsPowheg::_nlo_alphaS_opt, 0, false, false);
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
     "The value of alphaS to use for the nlo weight if _nlo_alphaS_opt=1",
     &MEPP2HiggsPowheg::_fixed_alphaS, 0.11803463, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<MEPP2HiggsPowheg,double> interfaceCorrectionCoefficient
    ("CorrectionCoefficient",
     "The magnitude of the correction term to reduce the negative contribution",
     &MEPP2HiggsPowheg::_a, 0.5, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<MEPP2HiggsPowheg,double> interfaceCorrectionPower
    ("CorrectionPower",
     "The power of the correction term to reduce the negative contribution",
     &MEPP2HiggsPowheg::_p, 0.7, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<MEPP2HiggsPowheg,unsigned int> interfaceWidthOption
    ("WidthOption",
     "Option to allow user to specify the width",
     &MEPP2HiggsPowheg::widthopt, 1, false, false);
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
     &MEPP2HiggsPowheg::usersWidth, GeV, 0.00468456293*GeV, 0.0*GeV, 10.0*GeV,
     true, false, Interface::limited);

  static Switch<MEPP2HiggsPowheg,unsigned int> interfaceFactorizationScaleOption
    ("FactorizationScaleOption",
     "Option for the choice of factorization scale",
     &MEPP2HiggsPowheg::scaleopt, 1, false, false);
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
     &MEPP2HiggsPowheg::mu_F, GeV, 100.0*GeV, 50.0*GeV, 500.0*GeV,
     true, false, Interface::limited);

  static Parameter<MEPP2HiggsPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before sHat if using a running scale",
     &MEPP2HiggsPowheg::_scaleFact, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Switch<MEPP2HiggsPowheg,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &MEPP2HiggsPowheg::shapeopt, 1, false, false);
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
     &MEPP2HiggsPowheg::processopt, 1, false, false);
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
     &MEPP2HiggsPowheg::minflavouropt, 4, 3, 5,
     false, false, Interface::limited);

  static Parameter<MEPP2HiggsPowheg,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the incoming quarks in the hard process",
     &MEPP2HiggsPowheg::maxflavouropt, 5, 3, 5,
     false, false, Interface::limited);
}

void MEPP2HiggsPowheg::doinit() throw(InitException) {
  // 2 thing for the NLO correction
  _CF = 4./3.; 
  _TR = 0.5;
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
  _mh = h0->mass();
  _wh = widthopt==2 ? usersWidth : h0->generateWidth(_mh);
  if(h0->massGenerator()) {
    _hmass=dynamic_ptr_cast<SMHiggsMassGeneratorPtr>(h0->massGenerator());
  }
  if(shapeopt==2&&!_hmass) throw InitException()
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
  return scaleopt == 1 ?  _scaleFact*sHat() : sqr(mu_F);
}

int MEPP2HiggsPowheg::nDim() const {
  return 0;
}

bool MEPP2HiggsPowheg::generateKinematics(const double *) {
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
  if(processopt==1||processopt==3) {
    tcPDPtr g=getParticleData(ParticleID::g);
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, h0, -1)));
  }
  // q qbar -> H processes
  if(processopt==1||processopt==2) {
    for (unsigned int i = minflavouropt; i <= maxflavouropt; ++i) {
      tcPDPtr q = getParticleData(i);
      tcPDPtr qb = q->CC();
    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, h0, -2)));
    }
  }
}

CrossSection MEPP2HiggsPowheg::dSigHatDR() const {
  using Constants::pi;
  InvEnergy2 bwfact;
  if(widthopt==2) {
    bwfact = _wh*sqrt(sHat())/pi/(sqr(sHat()-sqr(_mh))+sqr(_mh*_wh));
    double cs = me2()*jacobian()*pi*double(UnitRemoval::E4 * bwfact/sHat());
    return UnitRemoval::InvE2 * sqr(hbarc) * cs;
  }  
  if(shapeopt==1) {
    bwfact = mePartonData()[2]->generateWidth(sqrt(sHat()))*sqrt(sHat())/pi/
      (sqr(sHat()-sqr(_mh))+sqr(_mh*_wh));
  } else {
    bwfact = _hmass->BreitWignerWeight(sqrt(sHat()),0);
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
  hardvertex->ME(_me);
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
  if(calc) _me.reset(newme);
  // initial colour and spin factors: colour -> (8/64) and spin -> (1/4)
  return me2/32.;
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
  if(calc) _me.reset(newme);
  // final colour/spin factors
  return me2/12.;
}

double MEPP2HiggsPowheg::NLOweight() const {
  // If only leading order is required return 1:
  if(_contrib==0) return 1.;
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
  return _contrib==1 ? max(0.,wgt) : max(0.,-wgt);
}

double MEPP2HiggsPowheg::x(double xt, double v) const {
  double x0(xbar(v));
  return x0+(1.-x0)*xt;
}

double MEPP2HiggsPowheg::x_a(double x, double v) const {
  if(x==1.) return _xb_a;
  if(v==0.) return _xb_a;
  if(v==1.) return _xb_a/x;
  return (_xb_a/sqrt(x))*sqrt((1.-(1.-x)*(1.-v))/(1.-(1.-x)*v));
}

double MEPP2HiggsPowheg::x_b(double x, double v) const {
  if(x==1.) return _xb_b;
  if(v==0.) return _xb_b/x;
  if(v==1.) return _xb_b;
  return (_xb_b/sqrt(x))*sqrt((1.-(1.-x)*v)/(1.-(1.-x)*(1.-v)));
}

double MEPP2HiggsPowheg::xbar(double v) const {
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

double MEPP2HiggsPowheg::Ltilde_qq(double x, double v) const {
  if(x==1.) return 1.;
  double xa(x_a(x,v)),xb(x_b(x,v));
  double newq    = (_hadron_A->pdf()->xfx(_hadron_A,_parton_a,scale(),   xa)/   xa);
  double newqbar = (_hadron_B->pdf()->xfx(_hadron_B,_parton_b,scale(),   xb)/   xb);
  return( newq * newqbar / _oldq / _oldqbar );
}

double MEPP2HiggsPowheg::Ltilde_qg(double x, double v) const {
  double xa(x_a(x,v)),xb(x_b(x,v));
  double newq    = (_hadron_A->pdf()->xfx(_hadron_A,_parton_a,scale(),   xa)/   xa);
  double newg2   = (_hadron_B->pdf()->xfx(_hadron_B,_gluon   ,scale(),   xb)/   xb);
  return( newq * newg2 / _oldq / _oldqbar );
}

double MEPP2HiggsPowheg::Ltilde_gq(double x, double v) const {
  double xa(x_a(x,v)),xb(x_b(x,v));
  double newg1   = (_hadron_A->pdf()->xfx(_hadron_A,_gluon   ,scale(),   xa)/   xa);
  double newqbar = (_hadron_B->pdf()->xfx(_hadron_B,_parton_b,scale(),   xb)/   xb);
  return( newg1 * newqbar / _oldq / _oldqbar );
}

double MEPP2HiggsPowheg::Vtilde_qq() const {
  return _alphaS2Pi*_CF*(-3.*log(_mu2/_mll2)+(2.*sqr(Constants::pi)/3.)-8.);
}

double MEPP2HiggsPowheg::Ccalbar_qg(double x) const {
  return (sqr(x)+sqr(1.-x))*(log(_mll2/(_mu2*x))+2.*log(1.-x))+2.*x*(1.-x);
}

double MEPP2HiggsPowheg::Ctilde_qg(double x, double v) const {
  return  _alphaS2Pi*_TR * ((1.-xbar(v))/x) * Ccalbar_qg(x)*Ltilde_qg(x,v);
}

double MEPP2HiggsPowheg::Ctilde_gq(double x, double v) const {
  return  _alphaS2Pi*_TR * ((1.-xbar(v))/x) * Ccalbar_qg(x)*Ltilde_gq(x,v);
}

double MEPP2HiggsPowheg::Ctilde_qq(double x, double v) const {
  double wgt 
    = ((1.-x)/x+(1.+x*x)/(1.-x)/x*(2.*log(1.-x)-log(x)))*Ltilde_qq(x,v)
    -  4.*log(1.-x)/(1.-x)
    +  2./(1.-xbar(v))*log(1.-xbar(v))*log(1.-xbar(v))
    + (2./(1.-xbar(v))*log(1.-xbar(v))-2./(1.-x)+(1.+x*x)/x/(1.-x)*Ltilde_qq(x,v))
    *log(_mll2/_mu2);
  return _alphaS2Pi*_CF*(1.-xbar(v))*wgt;    
}

double MEPP2HiggsPowheg::Fcal_qq(double x, double v) const {
  return (sqr(1.-x)*(1.-2.*v*(1.-v))+2.*x)/x*Ltilde_qq(x,v);
}

double MEPP2HiggsPowheg::Fcal_qg(double x, double v) const {
  return ((1.-xbar(v))/x)*
    (2.*x*(1.-x)*v+sqr((1.-x)*v)+sqr(x)+sqr(1.-x))*Ltilde_qg(x,v);
}

double MEPP2HiggsPowheg::Fcal_gq(double x, double v) const {
  return ((1.-xbar(v))/x)*
    (2.*x*(1.-x)*(1.-v)+sqr((1.-x)*(1.-v))+sqr(x)+sqr(1.-x))*Ltilde_gq(x,v);
}

double MEPP2HiggsPowheg::Ftilde_qg(double xt, double v) const {
  return _alphaS2Pi*_TR*
    ( Fcal_qg(x(xt,v),v) - Fcal_qg(x(xt,0.),0.) )/v;
}

double MEPP2HiggsPowheg::Ftilde_gq(double xt, double v) const {
  return _alphaS2Pi*_TR*
    ( Fcal_gq(x(xt,v),v) - Fcal_gq(x(xt,1.),1.) )/(1.-v);
}

double MEPP2HiggsPowheg::Ftilde_qq(double xt, double v) const {
  return _alphaS2Pi*_CF*
    (( ( Fcal_qq(x(xt, v), v) - Fcal_qq(x(xt,1.),1.) ) / (1.-v)+
       ( Fcal_qq(x(xt, v), v) - Fcal_qq(x(xt,0.),0.) ) / v )/(1.-xt)
     + ( log(1.-xbar(v)) - log(1.-_xb_a))*2./(1.-v)
     + ( log(1.-xbar(v)) - log(1.-_xb_b))*2./v);
}
