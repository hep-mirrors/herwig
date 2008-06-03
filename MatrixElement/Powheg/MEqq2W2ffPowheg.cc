// -*- C++ -*-
//
// MEqq2W2ffPowheg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2W2ffPowheg class.
//

#include "MEqq2W2ffPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/MatrixElement/General/HardVertex.h"
#include <cassert>
#include "ThePEG/PDF/BeamParticleData.h"
#include "Herwig++/Utilities/Maths.h"
#include "Herwig++/Utilities/Histogram.h"
#include "Herwig++/MatrixElement/General/HardVertex.h"

using namespace Herwig;

void MEqq2W2ffPowheg::doinit() throw(InitException) {
  ME2to2Base::doinit();
  _wp=getParticleData(ThePEG::ParticleID::Wplus);
  _wm=getParticleData(ThePEG::ParticleID::Wminus);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) 
    _theFFWVertex = hwsm->vertexFFW();
  else
    throw InitException() << "Must be the Herwig++ StandardModel class in "
			  << "MEqq2W2ffPowheg::doinit" << Exception::abortnow;
}

void MEqq2W2ffPowheg::doinitrun() {  
  maxy = 0.;
  miny = 0.;
  no_wgts = 0;
  no_negwgts = 0;

}



void MEqq2W2ffPowheg::getDiagrams() const {
  // which intgermediates to include
  bool wplus =_plusminus==0||_plusminus==1;
  bool wminus=_plusminus==0||_plusminus==2;
  // possible incoming and outgoing particles
  typedef std::vector<pair<long,long> > Pairvector;
  // possible parents
  Pairvector parentpair;
  parentpair.reserve(6);
  // don't even think of putting 'break' in here!
  switch(_maxflavour) {
  case 5:
    parentpair.push_back(make_pair(ParticleID::b, ParticleID::cbar));
    parentpair.push_back(make_pair(ParticleID::b, ParticleID::ubar));
  case 4:
    parentpair.push_back(make_pair(ParticleID::s, ParticleID::cbar));
    parentpair.push_back(make_pair(ParticleID::d, ParticleID::cbar));
  case 3:
    parentpair.push_back(make_pair(ParticleID::s, ParticleID::ubar));
  case 2:
    parentpair.push_back(make_pair(ParticleID::d, ParticleID::ubar));
  default:
    ;
  }
  // possible children
  Pairvector childpair;
  childpair.reserve(9);
  childpair.push_back(make_pair(ParticleID::eminus,   ParticleID::nu_ebar));
  childpair.push_back(make_pair(ParticleID::muminus,  ParticleID::nu_mubar));
  childpair.push_back(make_pair(ParticleID::tauminus, ParticleID::nu_taubar));
  childpair.push_back(make_pair(ParticleID::d, ParticleID::ubar));
  childpair.push_back(make_pair(ParticleID::s, ParticleID::ubar));
  childpair.push_back(make_pair(ParticleID::b, ParticleID::ubar));
  childpair.push_back(make_pair(ParticleID::d, ParticleID::cbar));
  childpair.push_back(make_pair(ParticleID::s, ParticleID::cbar));
  childpair.push_back(make_pair(ParticleID::b, ParticleID::cbar));
  // loop over the children
  bool lepton,quark;
  Pairvector::const_iterator child = childpair.begin();
  for (; child != childpair.end(); ++child) {
    assert(child->first > 0 && child->second < 0);
    // allowed leptonic decay
    lepton = child->first > 10 && (_process==0 
				   || _process==2 
				   || (child->first - 5)/2 == int(_process));
    // allowed quark decay
    quark  = child->first < 10 && (_process==0 
				   || _process==1 
				   || (child->second == -2 && (child->first + 11)/2 == int(_process)) 
				   || (child->second == -4 && (child->first + 17)/2 == int(_process))
				   );
    // if decay not allowed skip
    if(!(quark || lepton)) continue;
    // decay products
    tcPDPtr lNeg1 = getParticleData(child->first);
    tcPDPtr lNeg2 = getParticleData(child->second);
    tcPDPtr lPos1 = lNeg2->CC();
    tcPDPtr lPos2 = lNeg1->CC();
    Pairvector::const_iterator parent = parentpair.begin();
    for (; parent != parentpair.end(); ++parent) {
      // parents
      tcPDPtr qNeg1 = getParticleData(parent->first);
      tcPDPtr qNeg2 = getParticleData(parent->second);
      tcPDPtr qPos1 = qNeg2->CC();
      tcPDPtr qPos2 = qNeg1->CC();
      // diagrams
      if(wminus) add(new_ptr((Tree2toNDiagram(2), 
			      qNeg1, qNeg2, 
			      1, _wm, 
			      3, lNeg1, 3, lNeg2, -1)));

      if(wplus)  add(new_ptr((Tree2toNDiagram(2), 
			      qPos1, qPos2, 
			      1, _wp, 
			      3, lPos1, 3, lPos2, -2)));
    }
  }
}

Energy2 MEqq2W2ffPowheg::scale() const {
  // Checked sqrt(sHat)-sqrt(lastX1()*lastX2()*lastS())<1 MeV for 10K events. 
  return 10000.0*GeV2;
//  return sHat();
}

double MEqq2W2ffPowheg::me2() const {
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction       q(meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction qbar(meMomenta()[1],mePartonData()[1],incoming);
  SpinorBarWaveFunction    f(meMomenta()[2],mePartonData()[2],outgoing);
  SpinorWaveFunction    fbar(meMomenta()[3],mePartonData()[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q.reset(ix)   ; fin.push_back(q);
    qbar.reset(ix); ain.push_back(qbar);
    f.reset(ix)   ;fout.push_back(f);
    fbar.reset(ix);aout.push_back(fbar);
  }
  return qqbarME(fin,ain,fout,aout,false);
}

unsigned int MEqq2W2ffPowheg::orderInAlphaS() const {
  return 0;
}

unsigned int MEqq2W2ffPowheg::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEqq2W2ffPowheg::diagrams(const DiagramVector &) const {
  Selector<DiagramIndex> sel;
  sel.insert(1.0, 0);
  return sel;
}

Selector<const ColourLines *>
MEqq2W2ffPowheg::colourGeometries(tcDiagPtr) const {
  static const ColourLines c1("1 -2");
  static const ColourLines c2("1 -2,4 -5");
  Selector<const ColourLines *> sel;
  if(abs(mePartonData()[2]->id())<=6) sel.insert(1.0, &c2);
  else                                sel.insert(1.0, &c1);
  return sel;
}

void MEqq2W2ffPowheg::persistentOutput(PersistentOStream & os) const {
  os << _maxflavour << _plusminus << _process << _theFFWVertex 
     << _wp << _wm << _contrib << _nlo_alphaS_opt << _fixed_alphaS
     << _a << _p;
}

void MEqq2W2ffPowheg::persistentInput(PersistentIStream & is, int) { 
  is >> _maxflavour >> _plusminus >> _process >> _theFFWVertex  
     >> _wp >> _wm >> _contrib >> _nlo_alphaS_opt >> _fixed_alphaS
     >> _a >> _p; 
}

ClassDescription<MEqq2W2ffPowheg> MEqq2W2ffPowheg::initMEqq2W2ffPowheg;
// Definition of the static class description member.

void MEqq2W2ffPowheg::Init() {

    static ClassDocumentation<MEqq2W2ffPowheg> documentation
    ("The MEqq2W2ffPowheg class implements the matrix element for"
     "q qbar to Standard Model fermions via W exchange using helicity amplitude"
     "techniques");

  static Parameter<MEqq2W2ffPowheg,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEqq2W2ffPowheg::_maxflavour, 5, 0, 5, false, false, true);

  static Switch<MEqq2W2ffPowheg,unsigned int> interfacePlusMinus
    ("Wcharge",
     "Which intermediate W bosons to include",
     &MEqq2W2ffPowheg::_plusminus, 0, false, false);
  static SwitchOption interfacePlusMinusAll
    (interfacePlusMinus,
     "Both",
     "Include W+ and W-",
     0);
  static SwitchOption interfacePlusMinusPlus
    (interfacePlusMinus,
     "Plus",
     "Only include W+",
     1);
  static SwitchOption interfacePlusMinusMinus
    (interfacePlusMinus,
     "Minus",
     "Only include W-",
     2);

  static Switch<MEqq2W2ffPowheg,unsigned int> interfaceProcess
    ("Process",
     "Which processes to include",
     &MEqq2W2ffPowheg::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all SM fermions as outgoing particles",
     0);
  static SwitchOption interfaceProcessQuarks
    (interfaceProcess,
     "Quarks",
     "Only include outgoing quarks",
     1);
  static SwitchOption interfaceProcessLeptons
    (interfaceProcess,
     "Leptons",
     "All include outgoing leptons",
     2);
  static SwitchOption interfaceProcessElectron
    (interfaceProcess,
     "Electron",
     "Only include outgoing e nu_e",
     3);
  static SwitchOption interfaceProcessMuon
    (interfaceProcess,
     "Muon",
     "Only include outgoing mu nu_mu",
     4);
  static SwitchOption interfaceProcessTau
    (interfaceProcess,
     "Tau",
     "Only include outgoing tauu nu_tau",
     5);
  static SwitchOption interfaceProcessUpDown
    (interfaceProcess,
     "UpDown",
     "Only include outgoing u dbar/ d ubar",
     6);
  static SwitchOption interfaceProcessUpStrange
    (interfaceProcess,
     "UpStrange",
     "Only include outgoing u sbar/ s ubar",
     7);
  static SwitchOption interfaceProcessUpBottom
    (interfaceProcess,
     "UpBottom",
     "Only include outgoing u bbar/ b ubar",
     8);
  static SwitchOption interfaceProcessCharmDown
    (interfaceProcess,
     "CharmDown",
     "Only include outgoing c dbar/ d cbar",
     9);
  static SwitchOption interfaceProcessCharmStrange
    (interfaceProcess,
     "CharmStrange",
     "Only include outgoing c sbar/ s cbar",
     10);
  static SwitchOption interfaceProcessCharmBottom
    (interfaceProcess,
     "CharmBottom",
     "Only include outgoing c bbar/ b cbar",
     11);
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

}

double MEqq2W2ffPowheg::qqbarME(vector<SpinorWaveFunction>    & fin ,
			  vector<SpinorBarWaveFunction> & ain ,
			  vector<SpinorBarWaveFunction> & fout,
			  vector<SpinorWaveFunction>    & aout,
			  bool calc) const {
  // matrix element to be stored
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  // positive or negative W boson
  bool positive = mePartonData()[0]->iCharge() 
    + mePartonData()[1]->iCharge() > 0;
  unsigned int ihel1,ihel2,ohel1,ohel2;
  // sum over helicities to get the matrix element
  double me = 0.;
  VectorWaveFunction inter;
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      for(ohel1=0;ohel1<2;++ohel1) {
	for(ohel2=0;ohel2<2;++ohel2) {
	  Complex diag = 0.0;
	  if (positive) {
	    // the Wp exchange
	    inter = _theFFWVertex->evaluate(scale(),1,_wp,fin[ihel1],ain[ihel2]);
	    diag = _theFFWVertex->evaluate(scale(),aout[ohel2],fout[ohel1],inter);
	  } 
	  else {
	    // the Wm exchange
	    inter = _theFFWVertex->evaluate(scale(),1,_wm,fin[ihel1],ain[ihel2]);
	    diag = _theFFWVertex->evaluate(scale(),aout[ohel2],fout[ohel1],inter);
	  }
	  // sum over helicities
	  me += real(diag*conj(diag));
	  if(calc) newme(ihel1,ihel2,ohel1,ohel2) = diag;
	}
      }
    }
  }
  // results
  // spin and colour factor
  double colspin=1./12.;
  if(abs(fout[0].id())<=6) colspin*=3.;
  if(calc) _me.reset(newme);
  
  // Get Born momentum fractions xbar_a and xbar_b:
  _xb_a = lastX1();
  _xb_b = lastX2();
  return me*NLOweight()*colspin;
}

void MEqq2W2ffPowheg::constructVertex(tSubProPtr sub) {
  SpinfoPtr spin[4];
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);hard.push_back(sub->outgoing()[1]);
  // order of particles
  unsigned int order[4]={0,1,2,3};
  if(hard[0]->id()<0){order[0]=1;order[1]=0;}
  if(hard[2]->id()<0){order[2]=3;order[3]=2;}
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction(   fin ,hard[order[0]],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[order[1]],incoming,false,true);
  SpinorBarWaveFunction(fout,hard[order[2]],outgoing,true ,true);
  SpinorWaveFunction(   aout,hard[order[3]],outgoing,true ,true);
  qqbarME(fin,ain,fout,aout,true);
  // get the spin info objects
  for(unsigned int ix=0;ix<4;++ix)
    {spin[ix]=dynamic_ptr_cast<SpinfoPtr>(hard[order[ix]]->spinInfo());}
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix){spin[ix]->setProductionVertex(hardvertex);}
}

int MEqq2W2ffPowheg::nDim() const {
  return 3;
}

bool MEqq2W2ffPowheg::generateKinematics(const double * r) {
  _xt=*(r+1);
  _v =*(r+2);
  //_xt=UseRandom::rnd();
  //_v=UseRandom::rnd();
  return ME2to2Base::generateKinematics(r);
}

double MEqq2W2ffPowheg::NLOweight() const {

  // If only leading order is required return 1:
  if(_contrib==0) return 1.;
  using Constants::pi;

  // Get particle data for QCD particles:
  _parton_a=mePartonData()[0];
  _parton_b=mePartonData()[1];
  _gluon=getParticleData(ParticleID::g);
  _hadron_A=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>(lastParticles().first->dataPtr());
  _hadron_B=dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>(lastParticles().second->dataPtr());
  // If necessary swap the particle data vectors so that _xb_a, 
  // mePartonData[0], beam[0] relate to the inbound quark: 

  if(!(lastPartons().first ->dataPtr()==_parton_a&&
       lastPartons().second->dataPtr()==_parton_b)) {
    swap(_xb_a    ,_xb_b);
    swap(_hadron_A,_hadron_B);
  }

  // Calculate alpha_S, CF, TF etc:
  if(_nlo_alphaS_opt==1) _alphaS = _fixed_alphaS;
  else _alphaS = SM().alphaS(scale());
  _CF = 4./3.; _TF = 0.5;
  // Calculate the invariant mass of the dilepton pair
  _mll2 = sHat();
  _mu2  = scale();

  // Calculate the integrand:
  double wqqvirt      = Vtilde_qq();
  double wqqcollin    = Ctilde_qq(x(_xt,1.),1.) + Ctilde_qq(x(_xt,0.),0.);
  double wqqreal      = Ftilde_qq(_xt,_v);
  double wqq          = wqqvirt+wqqcollin+wqqreal;

  double wqgcollin    = Ctilde_qg(x(_xt,0.),0.);
//  double wqgcollin    = Ctilde_qg_trick(x(_xt,0.),0.);
  double wqgreal      = Ftilde_qg(_xt,_v);
  double wqg          = wqgreal+wqgcollin;

  double wgqbarcollin = Ctilde_gq(x(_xt,1.),1.);
//  double wgqbarcollin = Ctilde_gq_trick(x(_xt,1.),1.);
  double wgqbarreal   = Ftilde_gq(_xt,_v);
  double wgqbar       = wgqbarreal+wgqbarcollin;

  double wgt          = 1.+(wqq+wqg+wgqbar);

  //trick to try and reduce neg wgt contribution
  if(_xt<1.-_eps) {
    wgt += _a*(1./pow(1.-_xt,_p)-(1.-pow(_eps,1.-_p))/(1.-_p)/(1.-_eps));
  }
  
  no_wgts ++;
  if( wgt < 0. ) no_negwgts ++;

  // Fill histogram bins
  int x_bin        = int(100.*_xt);
  int v_bin        = int(100.*_v);
  int xba_bin      = int(100.*_xb_a);
  int xbb_bin      = int(100.*_xb_b);
  int shatovrs_bin = int(100.*_xb_a*_xb_b);
  x_h[x_bin]     += wgt;
  v_h[v_bin]     += wgt;
  xba_h[xba_bin] += wgt;
  xbb_h[xbb_bin] += wgt;
  shatovrs_h[shatovrs_bin] += wgt;
  if ( wgt > 0. ) {
    x_pos_h[x_bin] += wgt;
    v_pos_h[v_bin] += wgt;
    xba_pos_h[xba_bin] += wgt;
    xbb_pos_h[xbb_bin] += wgt;
    shatovrs_pos_h[shatovrs_bin] += wgt;
  } 
  else if( wgt < 0. ) {
    x_neg_h[x_bin] += wgt;
    v_neg_h[v_bin] += wgt;
    xba_neg_h[xba_bin] += wgt;
    xbb_neg_h[xbb_bin] += wgt;
    shatovrs_neg_h[shatovrs_bin] += wgt;
  }

  if(lastY()>maxy) maxy=lastY();
  if(lastY()<miny) miny=lastY();

  if(wgt > _max_wgt){
    // cerr<<"maxwgt = "<<wgt<<" at xt = "<<_xt<<", vt = "<< _v<<"\n";
    _max_wgt = wgt;
  }
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

  double xa(x_a(x,v));
  double xb(x_b(x,v));

  double newq    = (_hadron_A->pdf()->xfx(_hadron_A,_parton_a,scale(),   xa)/   xa);
  double newqbar = (_hadron_B->pdf()->xfx(_hadron_B,_parton_b,scale(),   xb)/   xb);

  double oldq    = (_hadron_A->pdf()->xfx(_hadron_A,_parton_a,scale(),_xb_a)/_xb_a);
  double oldqbar = (_hadron_B->pdf()->xfx(_hadron_B,_parton_b,scale(),_xb_b)/_xb_b);

  return( newq * newqbar / oldq / oldqbar );

}
double MEqq2W2ffPowheg::Ltilde_qg(double x, double v) const {
  double xa(x_a(x,v));
  double xb(x_b(x,v));
  
  double newq    = (_hadron_A->pdf()->xfx(_hadron_A,_parton_a,scale(),   xa)/   xa);
  double newg2   = (_hadron_B->pdf()->xfx(_hadron_B,_gluon   ,scale(),   xb)/   xb);
    
  double oldq    = (_hadron_A->pdf()->xfx(_hadron_A,_parton_a,scale(),_xb_a)/_xb_a);
  double oldqbar = (_hadron_B->pdf()->xfx(_hadron_B,_parton_b,scale(),_xb_b)/_xb_b);

  return( newq * newg2 / oldq / oldqbar );

}
double MEqq2W2ffPowheg::Ltilde_gq(double x, double v) const {
  double xa(x_a(x,v));
  double xb(x_b(x,v));

  double newg1   = (_hadron_A->pdf()->xfx(_hadron_A,_gluon   ,scale(),   xa)/   xa);
  double newqbar = (_hadron_B->pdf()->xfx(_hadron_B,_parton_b,scale(),   xb)/   xb);

  double oldq    = (_hadron_A->pdf()->xfx(_hadron_A,_parton_a,scale(),_xb_a)/_xb_a);
  double oldqbar = (_hadron_B->pdf()->xfx(_hadron_B,_parton_b,scale(),_xb_b)/_xb_b);

  return( newg1 * newqbar / oldq / oldqbar );
}
double MEqq2W2ffPowheg::Vtilde_qq() const {
    return (_alphaS*_CF/(2.*pi))
          *(-3.*log(_mu2/_mll2)+(2.*pi*pi/3.)-8.);
}
double MEqq2W2ffPowheg::Ccalbar_qg(double x) const {
    return (sqr(x)+sqr(1.-x))*(log(_mll2/(_mu2*x))+2.*log(1.-x))+2.*x*(1.-x);
}
double MEqq2W2ffPowheg::Ctilde_qg(double x, double v) const {
    return  (_alphaS*_TF/(2.*pi))
	  * ((1.-xbar(v))/x)
	  * Ccalbar_qg(x)*Ltilde_qg(x,v);
}
double MEqq2W2ffPowheg::Ctilde_gq(double x, double v) const {
    return  (_alphaS*_TF/(2.*pi))
	  * ((1.-xbar(v))/x)
          * Ccalbar_qg(x)*Ltilde_gq(x,v);
}
double MEqq2W2ffPowheg::Ctilde_qq(double x, double v) const {
    double wgt 
       = ((1.-x)/x+(1.+x*x)/(1.-x)/x*(2.*log(1.-x)-log(x)))*Ltilde_qq(x,v)
       -  4.*log(1.-x)/(1.-x)
       +  2./(1.-xbar(v))*log(1.-xbar(v))*log(1.-xbar(v))
       + (2./(1.-xbar(v))*log(1.-xbar(v))-2./(1.-x)+(1.+x*x)/x/(1.-x)*Ltilde_qq(x,v))
	 *log(_mll2/_mu2);
    return (_alphaS*_CF/(2.*pi))*(1.-xbar(v))*wgt;    
}
double MEqq2W2ffPowheg::Fcal_qq(double x, double v) const {
    double tmp = (sqr(1.-x)*(1.-2.*v*(1.-v))+2.*x)/x;
    return tmp*Ltilde_qq(x,v);
}
double MEqq2W2ffPowheg::Fcal_qg(double x, double v) const {
    double tmp = 2.*x*(1.-x)*v+sqr((1.-x)*v)+sqr(x)+sqr(1.-x);
    return ((1.-xbar(v))/x)*tmp*Ltilde_qg(x,v);
}
double MEqq2W2ffPowheg::Fcal_gq(double x, double v) const {
    double tmp = 2.*x*(1.-x)*(1.-v)+sqr((1.-x)*(1.-v))+sqr(x)+sqr(1.-x);
    return ((1.-xbar(v))/x)*tmp*Ltilde_gq(x,v);
}
double MEqq2W2ffPowheg::Ftilde_qg(double xt, double v) const {
    double tmp = ( Fcal_qg(x(xt,v),v) - Fcal_qg(x(xt,0.),0.)
	         )/v;
    return (_alphaS*_TF/(2.*pi))*tmp;
}
double MEqq2W2ffPowheg::Ftilde_gq(double xt, double v) const {
    double tmp = ( Fcal_gq(x(xt,v),v) - Fcal_gq(x(xt,1.),1.)
	         )/(1.-v);
    return (_alphaS*_TF/(2.*pi))*tmp;
}
double MEqq2W2ffPowheg::Ftilde_qq(double xt, double v) const {
    double tmp = 
        ( ( Fcal_qq(x(xt, v), v) - Fcal_qq(x(xt,1.),1.) ) / (1.-v)
      +   ( Fcal_qq(x(xt, v), v) - Fcal_qq(x(xt,0.),0.) ) / v
	)/(1.-xt)
      + ( log(1.-xbar(v)) - log(1.-_xb_a))*2./(1.-v)
      + ( log(1.-xbar(v)) - log(1.-_xb_b))*2./v;
    return (_alphaS*_CF/(2.*pi))*tmp;
}


double MEqq2W2ffPowheg::Ctilde_qg_trick(double x, double v) const {
    double ctilde_qg_minus_ct = (_alphaS*_TF/(2.*pi))*((1.-xbar(v))/x) 
	                      * Ccalbar_qg(x)*(Ltilde_qg(x,v)-Ltilde_qg(1.,0.));
    double integrated_ct = (_alphaS*_TF/(2.*pi))
	                 * ( ((_xb_b-2.)*_xb_b+0.5*log(_xb_b))*log(_xb_b)
			   + 2.*ReLi2(_xb_b)
			   + 0.5-pi*pi/3.
			   - 2.*sqr(1.-_xb_b)*log(1.-_xb_b)
			   - 0.5*(4.-3.*_xb_b)*_xb_b
			   - (sqr(1.-_xb_b)+log(_xb_b))*log(_mll2/_mu2)
                           )*Ltilde_qg(1.,0.);
    return ctilde_qg_minus_ct + integrated_ct;
}
double MEqq2W2ffPowheg::Ctilde_gq_trick(double x, double v) const {
    double ctilde_gq_minus_ct = (_alphaS*_TF/(2.*pi))*((1.-xbar(v))/x)
	                      * Ccalbar_qg(x)*(Ltilde_gq(x,v)-Ltilde_gq(1.,1.));
    double integrated_ct = (_alphaS*_TF/(2.*pi))
	                 * ( ((_xb_a-2.)*_xb_a+0.5*log(_xb_a))*log(_xb_a)
			   + 2.*ReLi2(_xb_a)
			   + 0.5-pi*pi/3.
			   - 2.*sqr(1.-_xb_a)*log(1.-_xb_a)
			   - 0.5*(4.-3.*_xb_a)*_xb_a
			   - (sqr(1.-_xb_a)+log(_xb_a))*log(_mll2/_mu2)
                           )*Ltilde_gq(1.,1.);
    return ctilde_gq_minus_ct + integrated_ct;
}

void MEqq2W2ffPowheg::dofinish() {
  cerr << "\n";
  cerr << "alphaS = " << _alphaS << "\n";
  cerr << "a = " << _a << "\n";
  cerr << "Y ranged from "    << miny << " to " << maxy << "\n";
  cerr << "total wgts = " << no_wgts    << "\n"
       << "neg   wgts = " << no_negwgts << "\n";
  cerr << "percentage neg wgts = " << double(no_negwgts)/double(no_wgts)*100. << "% \n"; 

  ME2to2Base::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  
  //output xt histogram
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"xt distribution: all wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"xt\"\n";
  double step = 1. / 100.;
  double x=0.;
  for( unsigned int ix = 0; ix < x_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << x_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"xt distribution: pos wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"xt\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < x_pos_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << x_pos_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"xt distribution: neg wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"xt\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < x_neg_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << x_neg_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  //output v histogram
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"v distribution: all wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < v_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << v_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"v distribution: pos wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < v_pos_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << v_pos_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"v distribution: neg wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < v_neg_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << v_neg_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  //output _xb_a histogram
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"x_a (born) distribution: all wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x_a (born)\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < xba_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << xba_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"x_a (born) distribution: pos wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x_a (born)\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < xba_pos_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << xba_pos_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"x_a (born) distribution: neg wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x_a (born)\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < xba_neg_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << xba_neg_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  //output _xb_a histogram
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"x_b (born) distribution: all wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x_b (born)\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < xbb_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << xbb_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"x_b (born) distribution: pos wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x_b (born)\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < xbb_pos_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << xbb_pos_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"x_b (born) distribution: neg wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x_b (born)\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < xbb_neg_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << xbb_neg_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  //output shatovrs histogram
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"x_a*x_b (born) distribution: all wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x_a*x_b (born)\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < shatovrs_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << shatovrs_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"x_a*x_b (born) distribution: pos wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x_a*x_b (born)\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < shatovrs_pos_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << shatovrs_pos_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"x_a*x_b (born) distribution: neg wgts\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x_a*x_b (born)\"\n";
  step = 1. / 100.;
  x=0.;
  for( unsigned int ix = 0; ix < shatovrs_neg_h.size(); ++ix ) {
    x += step;
    outfile << x << "\t" << shatovrs_neg_h[ix].mean() << "\n";
  }
  outfile << "HIST\n";

  outfile.close();
}

