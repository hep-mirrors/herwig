// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISBase class.
//

#include "DISBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Utilities/Maths.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "Herwig/PDT/StandardMatchers.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include <numeric>
#include "Herwig/Shower/Base/ShowerProgenitor.h"
#include "Herwig/Shower/Base/ShowerTree.h"
#include "Herwig/Shower/Base/Branching.h"
#include "Herwig/Shower/Base/HardTree.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

// namespace {
// using namespace Herwig;
// using namespace ThePEG::Helicity;
// 
// void debuggingMatrixElement(bool BGF,const Lorentz5Momentum & pin,
// 			    const Lorentz5Momentum & p1,
// 			    const Lorentz5Momentum & p2,
// 			    tcPDPtr gluon,
// 			    const Lorentz5Momentum & pl1,
// 			    const Lorentz5Momentum & pl2,
// 			    const Lorentz5Momentum & pq1,
// 			    const Lorentz5Momentum & pq2,
// 			    tcPDPtr lepton1,tcPDPtr lepton2,
// 			    tcPDPtr quark1 ,tcPDPtr quark2,
// 			    Energy2 Q2,double phi, double x2, double x3, 
// 			    double xperp, double zp, double xp,
// 			    const vector<double> & azicoeff, 
// 			    bool normalize) {
//   tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>
//     (CurrentGenerator::current().standardModel());
//   assert(hwsm);
//   vector<AbstractFFVVertexPtr> weakVertex;
//   vector<PDPtr> bosons;
//   AbstractFFVVertexPtr strongVertex = hwsm->vertexFFG();
//   if(lepton1->id()==lepton2->id()) {
//     weakVertex.push_back(hwsm->vertexFFZ());
//     bosons.push_back(hwsm->getParticleData(ParticleID::Z0));
//     weakVertex.push_back(hwsm->vertexFFP());
//     bosons.push_back(hwsm->getParticleData(ParticleID::gamma));
//   }
//   else {
//     weakVertex.push_back(hwsm->vertexFFW());
//     bosons.push_back(hwsm->getParticleData(ParticleID::Wplus));
//   }
//   if(!BGF) {
//     SpinorWaveFunction    l1,q1,qp1;
//     SpinorBarWaveFunction l2,q2,qp2;
//     VectorWaveFunction    gl(p2,gluon,outgoing);
//     if(lepton1->id()>0) {
//       l1  = SpinorWaveFunction   (pl1,lepton1,incoming);
//       l2  = SpinorBarWaveFunction(pl2,lepton2,outgoing);
//     }
//     else {
//       l1  = SpinorWaveFunction   (pl2,lepton2,outgoing);
//       l2  = SpinorBarWaveFunction(pl1,lepton1,incoming);
//     }
//     if(quark1->id()>0) {
//       q1  = SpinorWaveFunction   (pq1,quark1,incoming);
//       q2  = SpinorBarWaveFunction(pq2,quark2,outgoing);
//       qp1 = SpinorWaveFunction   (pin,quark1,incoming);
//       qp2 = SpinorBarWaveFunction(p1 ,quark2,outgoing);
//     }
//     else {
//       q1  = SpinorWaveFunction   (pq2,quark2,outgoing);
//       q2  = SpinorBarWaveFunction(pq1,quark1,incoming);
//       qp1 = SpinorWaveFunction   (p1 ,quark2,outgoing);
//       qp2 = SpinorBarWaveFunction(pin,quark1,incoming);
//     }
//     double lome(0.),realme(0.);
//     for(unsigned int lhel1=0;lhel1<2;++lhel1) {
//       l1.reset(lhel1);
//       for(unsigned int lhel2=0;lhel2<2;++lhel2) { 
// 	l2.reset(lhel2);
// 	for(unsigned int qhel1=0;qhel1<2;++qhel1) {
// 	  q1.reset(qhel1);
// 	  qp1.reset(qhel1);
// 	  for(unsigned int qhel2=0;qhel2<2;++qhel2) {
// 	    q2.reset(qhel2);
// 	    qp2.reset(qhel2);
// 	    // leading order matrix element 
// 	    Complex diagLO(0.);
// 	    for(unsigned int ix=0;ix<weakVertex.size();++ix) {
// 	      VectorWaveFunction inter = 
// 		weakVertex[ix]->evaluate(Q2,3,bosons[ix],l1,l2);
// 	      diagLO += weakVertex[ix]->evaluate(Q2,q1,q2,inter);
// 	    }
// 	    lome   += norm(diagLO);
// 	    // real emission matrix element
// 	    for(unsigned int ghel=0;ghel<2;++ghel) {
//  	      gl.reset(2*ghel);
// 	      Complex diagReal(0.);
// 	      for(unsigned int ix=0;ix<weakVertex.size();++ix) {
// 		VectorWaveFunction inter = 
// 		  weakVertex[ix]->evaluate(Q2,3,bosons[ix],l1,l2);
// 		SpinorWaveFunction off1 = 
// 		  strongVertex->evaluate(Q2,5,qp1.particle(),qp1,gl);
// 		Complex diag1 = weakVertex[ix]->evaluate(Q2,off1,qp2,inter);
// 		SpinorBarWaveFunction off2 = 
// 		  strongVertex->evaluate(Q2,5,qp2.particle(),qp2,gl);
// 		Complex diag2 = weakVertex[ix]->evaluate(Q2,qp1,off2,inter);
// 		diagReal += diag1+diag2;
// 	      }
// 	      realme += norm(diagReal);
// 	    }
// 	  }
// 	}
//       }
//     }
//     double test1 = realme/lome/hwsm->alphaS(Q2)*Q2*UnitRemoval::InvE2;
//     double cphi(cos(phi));
//     double test2;
//     if(normalize) {
//       test2 = 8.*Constants::pi/(1.-xp)/(1.-zp)*
// 	(azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi))*
// 	(1.+sqr(xp)*(sqr(x2)+1.5*sqr(xperp)));
//     }
//     else {
//       test2 = 8.*Constants::pi/(1.-xp)/(1.-zp)*
// 	(azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi));
//     }
//     cerr << "testing RATIO A  " << test1/test2 << "\n";
//   }
//   else {
//     SpinorWaveFunction    l1,q1,qp1;
//     SpinorBarWaveFunction l2,q2,qp2;
//     VectorWaveFunction    gl(pin,gluon,incoming);
//     if(lepton1->id()>0) {
//       l1  = SpinorWaveFunction   (pl1,lepton1,incoming);
//       l2  = SpinorBarWaveFunction(pl2,lepton2,outgoing);
//     }
//     else {
//       l1  = SpinorWaveFunction   (pl2,lepton2,outgoing);
//       l2  = SpinorBarWaveFunction(pl1,lepton1,incoming);
//     }
//     if(quark1->id()>0) {
//       q1  = SpinorWaveFunction   (pq1,quark1      ,incoming);
//       q2  = SpinorBarWaveFunction(pq2,quark2      ,outgoing);
//       qp2 = SpinorBarWaveFunction(p1    ,quark2      ,outgoing);
//       qp1 = SpinorWaveFunction   (p2    ,quark1->CC(),outgoing);
//     }
//     else {
//       q1  = SpinorWaveFunction   (pq2,quark2      ,outgoing);
//       q2  = SpinorBarWaveFunction(pq1,quark1      ,incoming);
//       qp2 = SpinorBarWaveFunction(p2    ,quark1->CC(),outgoing);
//       qp1 = SpinorWaveFunction   (p1    ,quark2      ,outgoing);
//     }
//     double lome(0.),realme(0.);
//     for(unsigned int lhel1=0;lhel1<2;++lhel1) {
//       l1.reset(lhel1);
//       for(unsigned int lhel2=0;lhel2<2;++lhel2) { 
// 	l2.reset(lhel2);
// 	for(unsigned int qhel1=0;qhel1<2;++qhel1) {
// 	  q1.reset(qhel1);
// 	  qp1.reset(qhel1);
// 	  for(unsigned int qhel2=0;qhel2<2;++qhel2) {
// 	    q2.reset(qhel2);
// 	    qp2.reset(qhel2);
// 	    // leading order matrix element 
// 	    Complex diagLO(0.);
// 	    for(unsigned int ix=0;ix<weakVertex.size();++ix) {
// 	      VectorWaveFunction inter = 
// 		weakVertex[ix]->evaluate(Q2,3,bosons[ix],l1,l2);
// 	      diagLO += weakVertex[ix]->evaluate(Q2,q1,q2,inter);
// 	    }
// 	    lome   += norm(diagLO);
// 	    // real emission matrix element
// 	    for(unsigned int ghel=0;ghel<2;++ghel) {
//   	      gl.reset(2*ghel);
// 	      Complex diagReal(0.);
// 	      for(unsigned int ix=0;ix<weakVertex.size();++ix) {
// 		VectorWaveFunction inter = 
// 		  weakVertex[ix]->evaluate(Q2,3,bosons[ix],l1,l2);
// 		SpinorWaveFunction off1 = 
// 		  strongVertex->evaluate(Q2,5,qp1.particle(),qp1,gl);
// 		Complex diag1 = weakVertex[ix]->evaluate(Q2,off1,qp2,inter);
// 		SpinorBarWaveFunction off2 = 
// 		  strongVertex->evaluate(Q2,5,qp2.particle(),qp2,gl);
// 		Complex diag2 = weakVertex[ix]->evaluate(Q2,qp1,off2,inter);
// 		diagReal += diag1+diag2;
// 	      }
// 	      realme += norm(diagReal);
// 	    }
// 	  }
// 	}
//       }
//     }
//     double test1 = realme/lome/hwsm->alphaS(Q2)*Q2*UnitRemoval::InvE2;
//     double cphi(cos(phi));
//     double test2;
//     if(normalize) {
//       test2 = 8.*Constants::pi/zp/(1.-zp)*
// 	(azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi))*
// 	sqr(xp)*(sqr(x3)+sqr(x2)+3.*sqr(xperp));
//     }
//     else {
//       test2 = 8.*Constants::pi/zp/(1.-zp)*
// 	(azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi));
//     }
//     cerr << "testing RATIO B " << test1/test2 << "\n";
//   }
// }
// 
// }
 
DISBase::DISBase()  : initial_(6.), final_(3.),
		      procProb_(0.35),
		      comptonInt_(0.), bgfInt_(0.),
		      comptonWeight_(50.), BGFWeight_(150.), 
		      pTmin_(0.1*GeV), 
		      scaleOpt_(1),  muF_(100.*GeV), scaleFact_(1.),
		      contrib_(0), power_(0.1)
{}

DISBase::~DISBase() {}

void DISBase::persistentOutput(PersistentOStream & os) const {
  os << comptonInt_ << bgfInt_ << procProb_ << initial_ << final_ << alpha_
     << ounit(pTmin_,GeV) << comptonWeight_ << BGFWeight_ << gluon_
     << ounit(muF_,GeV) << scaleFact_ << scaleOpt_ << contrib_<< power_;
}

void DISBase::persistentInput(PersistentIStream & is, int) {
  is >> comptonInt_ >> bgfInt_ >> procProb_  >> initial_ >> final_ >> alpha_
     >> iunit(pTmin_,GeV) >> comptonWeight_ >> BGFWeight_ >> gluon_
     >> iunit(muF_,GeV) >> scaleFact_ >> scaleOpt_ >> contrib_ >> power_;
}

AbstractClassDescription<DISBase> DISBase::initDISBase;
// Definition of the static class description member.

void DISBase::Init() {
  
  static ClassDocumentation<DISBase> documentation
    ("The DISBase class provides the base class for the "
     "implementation of DIS type processes including the "
     "hard corrections in either the old-fashioned matrix "
     "element correction of POWHEG approaches");

  static Parameter<DISBase,double> interfaceProcessProbability
    ("ProcessProbability",
     "The probabilty of the QCD compton process for the process selection",
     &DISBase::procProb_, 0.3, 0.0, 1.,
     false, false, Interface::limited);

  static Reference<DISBase,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &DISBase::alpha_, false, false, true, false, false);
  
  static Parameter<DISBase,Energy> interfacepTMin
    ("pTMin",
     "The minimum pT",
     &DISBase::pTmin_, GeV, 1.*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DISBase,double> interfaceComptonWeight
    ("ComptonWeight",
     "Weight for the overestimate ofthe compton channel",
     &DISBase::comptonWeight_, 50.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<DISBase,double> interfaceBGFWeight
    ("BGFWeight",
     "Weight for the overestimate of the BGF channel",
     &DISBase::BGFWeight_, 100.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Switch<DISBase,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &DISBase::contrib_, 0, false, false);
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

  static Switch<DISBase,unsigned int> interfaceScaleOption
    ("ScaleOption",
     "Option for the choice of factorization (and renormalization) scale",
     &DISBase::scaleOpt_, 1, false, false);
  static SwitchOption interfaceDynamic
    (interfaceScaleOption,
     "Dynamic",
     "Dynamic factorization scale equal to the current sqrt(sHat())",
     1);
  static SwitchOption interfaceFixed
    (interfaceScaleOption,
     "Fixed",
     "Use a fixed factorization scale set with FactorizationScaleValue",
     2);

  static Parameter<DISBase,Energy> interfaceFactorizationScale
    ("FactorizationScale",
     "Value to use in the event of a fixed factorization scale",
     &DISBase::muF_, GeV, 100.0*GeV, 1.0*GeV, 500.0*GeV,
     true, false, Interface::limited);

  static Parameter<DISBase,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before Q2 if using a running scale",
     &DISBase::scaleFact_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DISBase,double> interfaceSamplingPower
    ("SamplingPower",
     "Power for the sampling of xp",
     &DISBase::power_, 0.6, 0.0, 1.,
     false, false, Interface::limited);
}

void DISBase::doinit() {
  HwMEBase::doinit();
  // integrals of me over phase space
  double r5=sqrt(5.),darg((r5-1.)/(r5+1.)),ath(0.5*log((1.+1./r5)/(1.-1./r5)));
  comptonInt_ = 2.*(-21./20.-6./(5.*r5)*ath+sqr(Constants::pi)/3.
		    -2.*Math::ReLi2(1.-darg)-2.*Math::ReLi2(1.-1./darg));
  bgfInt_ = 121./9.-56./r5*ath;
  // extract the gluon ParticleData objects
  gluon_ = getParticleData(ParticleID::g);
}

void DISBase::initializeMECorrection(ShowerTreePtr tree, double & initial,
				     double & final) {
  initial = initial_;
  final   = final_;
  // incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data())) {
      partons_[0] = cit->first->progenitor()->dataPtr();
      pq_[0] = cit->first->progenitor()->momentum();
    }
    else if(LeptonMatcher::Check(cit->first->progenitor()->data())) {
      leptons_[0] = cit->first->progenitor()->dataPtr();
      pl_[0] = cit->first->progenitor()->momentum();
    }
  }
  // outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data())) {
      partons_[1] = cit->first->progenitor()->dataPtr();
      pq_[1] = cit->first->progenitor()->momentum();
    }
    else if(LeptonMatcher::Check(cit->first->progenitor()->data())) {
      leptons_[1] = cit->first->progenitor()->dataPtr();
      pl_[1] = cit->first->progenitor()->momentum();
    }
  }
  // extract the born variables
  q_ =pl_[0]-pl_[1];
  q2_ = -q_.m2();
  double  yB = (q_*pq_[0])/(pl_[0]*pq_[0]); 
  l_ = 2./yB-1.;
  // calculate the A coefficient for the correlations
  acoeff_ = A(leptons_[0],leptons_[1],
	      partons_[0],partons_[1],q2_);
}

void DISBase::applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  static const double eps=1e-6;
  // find the incoming and outgoing quarks and leptons
  ShowerParticlePtr quark[2],lepton[2];
  PPtr hadron;
  // incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data())) {
      hadron = cit->first->original()->parents()[0];
      quark [0] = cit->first->progenitor();
      beam_ = cit->first->beam();
    }
    else if(LeptonMatcher::Check(cit->first->progenitor()->data())) {
      lepton[0] = cit->first->progenitor();
    }
  }
  pdf_ = beam_->pdf();
  assert(beam_&&pdf_&&quark[0]&&lepton[0]);
  // outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data()))
      quark [1] = cit->first->progenitor();
    else if(LeptonMatcher::Check(cit->first->progenitor()->data())) {
      lepton[1] = cit->first->progenitor();
    }
  }
  // momentum fraction
  assert(quark[1]&&lepton[1]);
  xB_ = quark[0]->x();
  // calculate the matrix element
  vector<double> azicoeff;
  // select the type of process
  bool BGF = UseRandom::rnd()>procProb_;
  double xp,zp,wgt,x1,x2,x3,xperp;
  // generate a QCD compton process
  if(!BGF) {
    wgt = generateComptonPoint(xp,zp);
    if(xp<eps) return;
    // common pieces
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1.);
    wgt *= 2./3./Constants::pi*alpha_->value(scale)/procProb_;
    // PDF piece
    wgt *= pdf_->xfx(beam_,quark[0]->dataPtr(),scale,xB_/xp)/
           pdf_->xfx(beam_,quark[0]->dataPtr(),q2_  ,xB_);
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = ComptonME(xp,x2,xperp,true);
  }
  // generate a BGF process
  else {
    wgt = generateBGFPoint(xp,zp);
    if(xp<eps) return;
    // common pieces 
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1);
    wgt *= 0.25/Constants::pi*alpha_->value(scale)/(1.-procProb_);
    // PDF piece
    wgt *= pdf_->xfx(beam_,gluon_              ,scale,xB_/xp)/
           pdf_->xfx(beam_,quark[0]->dataPtr(),q2_  ,xB_);
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = BGFME(xp,x2,x3,xperp,true);
  }
  // compute the azimuthal average of the weight
  wgt *= (azicoeff[0]+0.5*azicoeff[2]);
  // decide whether or not to accept the weight
  if(UseRandom::rnd()>wgt) return;
  // if generate generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi));
    phiwgt = azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi);
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in DISMECorrection"
				  << "::applyHardMatrixElementCorrection() to"
				  << " generate phi" << Exception::eventerror;
  // construct lorentz transform from lab to breit frame
  Lorentz5Momentum phadron =  hadron->momentum();
  phadron.setMass(0.*GeV);
  phadron.rescaleEnergy();
  Lorentz5Momentum pcmf = phadron+0.5/xB_*q_;
  pcmf.rescaleMass();
  LorentzRotation rot(-pcmf.boostVector());
  Lorentz5Momentum pbeam = rot*phadron;
  Axis axis(pbeam.vect().unit());
  double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
  rot.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  Lorentz5Momentum pl    = rot*pl_[0];
  rot.rotateZ(-atan2(pl.y(),pl.x()));
  pl_[0] *= rot;
  pl_[1] *= rot;
  pq_[0] *= rot;
  pq_[1] *= rot;
  // compute the new incoming and outgoing momenta
  Energy Q(sqrt(q2_));
  Lorentz5Momentum p1 = Lorentz5Momentum( 0.5*Q*xperp*cos(phi), 0.5*Q*xperp*sin(phi),
					  -0.5*Q*x2,0.*GeV,0.*GeV);
  p1.rescaleEnergy();
  Lorentz5Momentum p2 = Lorentz5Momentum(-0.5*Q*xperp*cos(phi),-0.5*Q*xperp*sin(phi),
					 -0.5*Q*x3,0.*GeV,0.*GeV);
  p2.rescaleEnergy();
  Lorentz5Momentum pin(0.*GeV,0.*GeV,-0.5*x1*Q,-0.5*x1*Q,0.*GeV);
//   debuggingMatrixElement(BGF,pin,p1,p2,gluon_,pl_[0],pl_[1],pq_[0],pq_[1],
// 			 lepton[0]->dataPtr(),lepton[1]->dataPtr(),
// 			 quark [0]->dataPtr(),quark [1]->dataPtr(),
// 			 q2_,phi,x2,x3,xperp,zp,xp,azicoeff,true);
  // we need the Lorentz transform back to the lab
  rot.invert();
  // transform the momenta to lab frame
  pin *= rot;
  p1  *= rot;
  p2  *= rot;
  // test to ensure outgoing particles can be put on-shell
  if(!BGF) {
    if(p1.e()<quark[1]->dataPtr()->constituentMass()) return;
    if(p2.e()<gluon_              ->constituentMass()) return;
  }
  else {
    if(p1.e()<quark[1]->dataPtr()      ->constituentMass()) return;
    if(p2.e()<quark[0]->dataPtr()->CC()->constituentMass()) return;
  }
  // create the new particles and add to ShowerTree
  bool isquark = quark[0]->colourLine();
  if(!BGF) {
    PPtr newin  = new_ptr(Particle(*quark[0]));
    newin->set5Momentum(pin);
    PPtr newg   = gluon_              ->produceParticle(p2 );
    PPtr newout = quark[1]->dataPtr()->produceParticle(p1 ); 
    ColinePtr col=isquark ? 
      quark[0]->colourLine() : quark[0]->antiColourLine();
    ColinePtr newline=new_ptr(ColourLine());
    // final-state emission
    if(xp>zp) {
      col->removeColoured(newout,!isquark);
      col->addColoured(newin,!isquark);
      col->addColoured(newg,!isquark);
      newline->addColoured(newg,isquark);
      newline->addColoured(newout,!isquark);
    }
    // initial-state emission
    else {
      col->removeColoured(newin ,!isquark);
      col->addColoured(newout,!isquark);
      col->addColoured(newg,isquark);
      newline->addColoured(newg,!isquark);
      newline->addColoured(newin,!isquark);
    }
    PPtr orig;
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	  cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[0]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newin);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
      cit->first->progenitor(sp);
      tree->incomingLines()[cit->first]=sp;
      sp->x(xB_/xp);
      cit->first->perturbative(xp>zp);
      if(xp<=zp) orig=cit->first->original();
    }
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[1]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newout);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newout,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(xp<=zp);
      if(xp>zp) orig=cit->first->original();
    }
    assert(orig);
    // add the gluon
    ShowerParticlePtr sg=new_ptr(ShowerParticle(*newg,1,true));
    ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(orig,newg,sg));
    gluon->perturbative(false);
    tree->outgoingLines().insert(make_pair(gluon,sg));
    tree->hardMatrixElementCorrection(true);
  }
  else {
    PPtr newin   = gluon_                   ->produceParticle(pin);
    PPtr newqbar = quark[0]->dataPtr()->CC()->produceParticle(p2 );
    PPtr newout  = quark[1]->dataPtr()      ->produceParticle(p1 );
    ColinePtr col=isquark ? quark[0]->colourLine() : quark[0]->antiColourLine();
    ColinePtr newline=new_ptr(ColourLine()); 
    col    ->addColoured(newin  ,!isquark);
    newline->addColoured(newin  , isquark);
    col    ->addColoured(newout ,!isquark);
    newline->addColoured(newqbar, isquark);
    PPtr orig;
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	  cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[0]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newin);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
      cit->first->progenitor(sp);
      tree->incomingLines()[cit->first]=sp;
      sp->x(xB_/xp);
      cit->first->perturbative(false);
      orig=cit->first->original();
    }
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[1]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newout);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newout,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(true);
    }
    assert(orig);
    // add the (anti)quark
    ShowerParticlePtr sqbar=new_ptr(ShowerParticle(*newqbar,1,true));
    ShowerProgenitorPtr qbar=new_ptr(ShowerProgenitor(orig,newqbar,sqbar));
    qbar->perturbative(false);
    tree->outgoingLines().insert(make_pair(qbar,sqbar));
    tree->hardMatrixElementCorrection(true);
  }
}

bool DISBase::softMatrixElementVeto(ShowerProgenitorPtr initial,
				    ShowerParticlePtr parent, Branching br) {
  bool veto = !UseRandom::rndbool(parent->isFinalState() ? 1./final_ : 1./initial_);
  // check if me correction should be applied
  long id[2]={initial->id(),parent->id()};
  if(id[0]!=id[1]||id[1]==ParticleID::g) return veto;
  // get the pT
  Energy pT=br.kinematics->pT();
  // check if hardest so far
  if(pT<initial->highestpT()) return veto;
  double kappa(sqr(br.kinematics->scale())/q2_),z(br.kinematics->z());
  double zk((1.-z)*kappa);
  // final-state
  double wgt(0.);
  if(parent->isFinalState()) {
    double zp=z,xp=1./(1.+z*zk);
    double xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    double x2 = 1.-(1.-zp)/xp;
    vector<double> azicoeff = ComptonME(xp,x2,xperp,false);
    wgt = (azicoeff[0]+0.5*azicoeff[2])*xp/(1.+sqr(z))/final_;
    if(wgt<.0||wgt>1.) {
      ostringstream wstring;
      wstring << "Soft ME correction weight too large or "
	      << "negative for FSR in DISBase::"
	      << "softMatrixElementVeto() soft weight " 
	      << " xp = " << xp << " zp = " << zp
	      << " weight = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(), 
					 Exception::warning) );
    }
  }
  else {
    double xp = 2.*z/(1.+zk+sqrt(sqr(1.+zk)-4.*z*zk));
    double zp = 0.5* (1.-zk+sqrt(sqr(1.+zk)-4.*z*zk));
    double xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    double x1 = -1./xp, x2 = 1.-(1.-zp)/xp, x3 = 2.+x1-x2;
    // compton
    if(br.ids[0]!=ParticleID::g) {
      vector<double> azicoeff = ComptonME(xp,x2,xperp,false);
      wgt = (azicoeff[0]+0.5*azicoeff[2])*xp*(1.-z)/(1.-xp)/(1.+sqr(z))/
	(1.-zp+xp-2.*xp*(1.-zp));
    }
    // BGF
    else {
      vector<double> azicoeff = BGFME(xp,x2,x3,xperp,true);
      wgt = (azicoeff[0]+0.5*azicoeff[2])*xp/(1.-zp+xp-2.*xp*(1.-zp))/(sqr(z)+sqr(1.-z));
    }
    wgt /=initial_;
    if(wgt<.0||wgt>1.) {
      ostringstream wstring;
      wstring << "Soft ME correction weight too large or "
	      << "negative for ISR in DISBase::"
	      << "softMatrixElementVeto() soft weight " 
	      << " xp = " << xp << " zp = " << zp
	      << " weight = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(), 
					 Exception::warning) );
    }
  }
  // if not vetoed
  if(UseRandom::rndbool(wgt)) return false;
  // otherwise
  parent->vetoEmission(br.type,br.kinematics->scale());
  return true;
}

double DISBase::generateComptonPoint(double &xp, double & zp) {
  static const double maxwgt = 1.;
  double wgt;
  do {
    xp  = UseRandom::rnd();
    double zpmin = xp, zpmax = 1./(1.+xp*(1.-xp));
    zp = 1.-pow((1.-zpmin)/(1.-zpmax),UseRandom::rnd())*(1.-zpmax);
    wgt = log((1.-zpmin)/(1.-zpmax))*(1.-zp);
    if(UseRandom::rndbool()) swap(xp,zp);
    double xperp2 = 4.*(1.-xp)*(1.-zp)*zp/xp,x2=1.-(1.-zp)/xp;
    wgt *= 2.*(1.+sqr(xp)*(sqr(x2)+1.5*xperp2))/(1.-xp)/(1.-zp);
    if(wgt>maxwgt) {
      ostringstream wstring;
      wstring << "DISBase::generateComptonPoint "
	      << "Weight greater than maximum "
	      << "wgt = " << wgt << " maxwgt = 1\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
  }
  while(wgt<UseRandom::rnd()*maxwgt);
  return comptonInt_;
}

double DISBase::generateBGFPoint(double &xp, double & zp) {
  static const double maxwgt = 25.;
  double wgt;
  do {
    xp = UseRandom::rnd();
    double zpmax = 1./(1.+xp*(1.-xp)), zpmin = 1.-zpmax;
    zp = 1.-pow((1.-zpmin)/(1.-zpmax),UseRandom::rnd())*(1.-zpmax);
    wgt = log((1.-zpmin)/(1.-zpmax))*(1.-zp);
    double x1 = -1./xp;
    double x2 = 1.-(1.-zp)/xp;
    double x3 = 2.+x1-x2;
    double xperp2 = 4.*(1.-xp)*(1.-zp)*zp/xp;
    wgt *= sqr(xp)/(1.-zp)*(sqr(x3)+sqr(x2)+3.*xperp2);
    if(wgt>maxwgt) {
      ostringstream wstring;
      wstring << "DISBase::generateBGFPoint "
	      << "Weight greater than maximum "
	      << "wgt = " << wgt << " maxwgt = 1\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
  }
  while(wgt<UseRandom::rnd()*maxwgt);
  return bgfInt_;
//   static const double maxwgt = 2.,npow=0.34,ac=1.0;
//   double wgt;
//   do {
//     double rho = UseRandom::rnd();
//     xp = 1.-pow(rho,1./(1.-npow));
//     wgt = (sqr(xp)+ac+sqr(1.-xp));
//     if(wgt>1.+ac) cerr << "testing violates BGF maxA " << wgt << "\n";
//   }
//   while(wgt<UseRandom::rnd()*(1.+ac));
//   double xpwgt = -((6.-5.*npow+sqr(npow))*ac-3.*npow+sqr(npow)+4) 
//     /(sqr(npow)*(npow-6.)+11.*npow-6.);
//   xpwgt *= pow(1.-xp,npow)/wgt;
//   double xp2(sqr(xp)),lxp(log(xp)),xp4(sqr(xp2)),lxp1(log(1.-xp));
//   double zpwgt = (2.*xp4*(lxp+lxp1-3.)+4.*xp*xp2*(3.-lxp-lxp1)
// 		  +xp2*(-13.+lxp+lxp1)+xp*(+7.+lxp+lxp1)-lxp-lxp1-1.)/(1.+xp-xp2);
//   do {
//     double zpmax = 1./(1.+xp*(1.-xp)), zpmin = 1.-zpmax;
//     zp = 1.-pow((1.-zpmin)/(1.-zpmax),UseRandom::rnd())*(1.-zpmax);
//     wgt = log((1.-zpmin)/(1.-zpmax))*(1.-zp);
//     double x1 = -1./xp;
//     double x2 = 1.-(1.-zp)/xp;
//     double x3 = 2.+x1-x2;
//     double xperp2 = 4.*(1.-xp)*(1.-zp)*zp/xp;
//     wgt *= sqr(xp)/(1.-zp)*(sqr(x3)+sqr(x2)+3.*xperp2);
//     if(wgt>maxwgt*zpwgt) cerr << "testing violates BGF maxB " << wgt/xpwgt << "\n";
//   }
//   while(wgt<UseRandom::rnd()*maxwgt);
//   return zpwgt*xpwgt;
}

vector<double> DISBase::ComptonME(double xp, double x2, double xperp,
				  bool norm) {
  vector<double> output(3,0.);
  double cos2 =   x2 /sqrt(sqr(x2)+sqr(xperp));
  double sin2 = xperp/sqrt(sqr(x2)+sqr(xperp));
  double root = sqrt(sqr(l_)-1.);
  output[0] = sqr(cos2)+acoeff_*cos2*l_+sqr(l_);
  output[1] = -acoeff_*cos2*root*sin2-2.*l_*root*sin2;
  output[2] = sqr(root)*sqr(sin2);
  double lo(1+acoeff_*l_+sqr(l_));
  double denom = norm ? 1.+sqr(xp)*(sqr(x2)+1.5*sqr(xperp)) : 1.;
  double fact  = sqr(xp)*(sqr(x2)+sqr(xperp))/lo;
  for(unsigned int ix=0;ix<output.size();++ix) 
    output[ix] = ((ix==0 ? 1. : 0.) +fact*output[ix])/denom;
  return output;
}

vector<double> DISBase::BGFME(double xp, double x2, double x3, 
			      double xperp, bool norm) {
  vector<double> output(3,0.);
  double cos2  =   x2 /sqrt(sqr(x2)+sqr(xperp));
  double sin2  = xperp/sqrt(sqr(x2)+sqr(xperp));
  double fact2 = sqr(xp)*(sqr(x2)+sqr(xperp));
  double cos3  =   x3 /sqrt(sqr(x3)+sqr(xperp));
  double sin3  = xperp/sqrt(sqr(x3)+sqr(xperp));
  double fact3 = sqr(xp)*(sqr(x3)+sqr(xperp));
  double root = sqrt(sqr(l_)-1.);
  output[0] = fact2*(sqr(cos2)+acoeff_*cos2*l_+sqr(l_)) +
              fact3*(sqr(cos3)-acoeff_*cos3*l_+sqr(l_));
  output[1] = - fact2*(acoeff_*cos2*root*sin2+2.*l_*root*sin2)
              - fact3*(acoeff_*cos3*root*sin3-2.*l_*root*sin3);
  output[2] = fact2*(sqr(root)*sqr(sin2)) +
              fact3*(sqr(root)*sqr(sin3));
  double lo(1+acoeff_*l_+sqr(l_));
  double denom = norm ? sqr(xp)*(sqr(x3)+sqr(x2)+3.*sqr(xperp))*lo : lo;
  for(unsigned int ix=0;ix<output.size();++ix) output[ix] /= denom;
  return output;
}

HardTreePtr DISBase::generateHardest(ShowerTreePtr tree,
				     vector<ShowerInteraction::Type> inter) {
  bool found = false;
  // check if generating QCD radiation
  for(unsigned int ix=0;ix<inter.size();++ix) {
    found |= inter[ix]==ShowerInteraction::QCD;
  }
  if(!found) return HardTreePtr();
  ShowerParticlePtr quark[2],lepton[2];
  PPtr hadron;
  // incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data())) {
      hadron = cit->first->original()->parents()[0];
      quark [0] = cit->first->progenitor();
      beam_ = cit->first->beam();
    }
    else if(LeptonMatcher::Check(cit->first->progenitor()->data())) {
      lepton[0] = cit->first->progenitor();
      leptons_[0] = lepton[0]->dataPtr();
    }
  }
  pdf_=beam_->pdf();
  assert(beam_&&pdf_&&quark[0]&&lepton[0]);
  // outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data()))
      quark [1] = cit->first->progenitor();
    else if(LeptonMatcher::Check(cit->first->progenitor()->data())) {
      lepton[1] = cit->first->progenitor();
      leptons_[1] = lepton[1]->dataPtr();
    }
  }
  assert(quark[1]&&lepton[1]);
  // Particle data objects
  for(unsigned int ix=0;ix<2;++ix) partons_[ix] = quark[ix]->dataPtr();
  // extract the born variables
  q_ =lepton[0]->momentum()-lepton[1]->momentum();
  q2_ = -q_.m2();
  xB_ = quark[0]->x();
  double  yB = 
    (                   q_*quark[0]->momentum())/
    (lepton[0]->momentum()*quark[0]->momentum()); 
  l_ = 2./yB-1.;
  // construct lorentz transform from lab to breit frame
  Lorentz5Momentum phadron =  hadron->momentum();
  phadron.setMass(0.*GeV);
  phadron.rescaleRho();
  Lorentz5Momentum pb     = quark[0]->momentum();
  Axis axis(q_.vect().unit());
  double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
  LorentzRotation rot_ = LorentzRotation();
  if(axis.perp2()>1e-20) {
    rot_.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    rot_.rotateX(Constants::pi);
  }
  if(abs(1.-q_.e()/q_.vect().mag())>1e-6) rot_.boostZ( q_.e()/q_.vect().mag());
  pb *= rot_;
  if(pb.perp2()/GeV2>1e-20) {
    Boost trans = -1./pb.e()*pb.vect();
    trans.setZ(0.);
    rot_.boost(trans);
  }
  Lorentz5Momentum pl    = rot_*lepton[0]->momentum();
  rot_.rotateZ(-atan2(pl.y(),pl.x()));
  // momenta of the particles
  pl_[0]=rot_*lepton[0]->momentum();
  pl_[1]=rot_*lepton[1]->momentum();
  pq_[0]=rot_* quark[0]->momentum();
  pq_[1]=rot_* quark[1]->momentum();
  q_ *= rot_;
  // coefficient for the matrix elements
  acoeff_ = A(lepton[0]->dataPtr(),lepton[1]->dataPtr(),
	      quark [0]->dataPtr(),quark [1]->dataPtr(),q2_);
  // generate a compton point
  generateCompton();
  generateBGF();
  // no valid emission, return
  if(pTCompton_<ZERO&&pTBGF_<ZERO) return HardTreePtr();
  // type of emission, pick highest pT
  bool isCompton=pTCompton_>pTBGF_;
//   // find the sudakov for the branching
//   SudakovPtr sudakov;
//   // ISR
//   if(ComptonISFS_||!isCompton) {
//     BranchingList branchings=evolver()->splittingGenerator()->initialStateBranchings();
//     long index = abs(partons_[0]->id());
//     IdList br(3);
//     if(isCompton) {
//       br[0] = index;
//       br[1] = index;
//       br[2] = ParticleID::g;
//     }
//     else {
//       br[0] = ParticleID::g;
//       br[1] =  abs(partons_[0]->id());
//       br[2] = -abs(partons_[0]->id());
//     }
//     for(BranchingList::const_iterator cit = branchings.lower_bound(index); 
// 	cit != branchings.upper_bound(index); ++cit ) {
//       IdList ids = cit->second.second;
//       if(ids[0]==br[0]&&ids[1]==br[1]&&ids[2]==br[2]) {
// 	sudakov=cit->second.first;
// 	break;
//       }
//     }
//   }
//   // FSR
//   else {
//     BranchingList branchings = 
//       evolver()->splittingGenerator()->finalStateBranchings();
//     long index=abs(partons_[1]->id());
//     for(BranchingList::const_iterator cit = branchings.lower_bound(index);
// 	cit != branchings.upper_bound(index); ++cit ) {
//       IdList ids = cit->second.second;
//       if(ids[0]==index&&ids[1]==index&&ids[2]==ParticleID::g) {
// 	sudakov = cit->second.first;
// 	break; 	    
//       }
//     }
//   }
//   if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
// 				 << "DISBase::generateHardest()" 
// 				 << Exception::runerror;
  // add the leptons
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  spaceBranchings.push_back(new_ptr(HardBranching(lepton[0],SudakovPtr(),
						  HardBranchingPtr(),
						  HardBranching::Incoming)));
  allBranchings.push_back(spaceBranchings.back());
  allBranchings.push_back(new_ptr(HardBranching(lepton[1],SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  // compton hardest
  if(isCompton) {
    rot_.invert();
    for(unsigned int ix=0;ix<ComptonMomenta_.size();++ix) {
      ComptonMomenta_[ix].transform(rot_);
    }
    ShowerParticlePtr newqout (new_ptr(ShowerParticle(partons_[1],true)));
    newqout->set5Momentum(ComptonMomenta_[1]);
    ShowerParticlePtr newg(new_ptr(ShowerParticle(gluon_,true)));
    newg->set5Momentum(ComptonMomenta_[2]);
    ShowerParticlePtr newqin   (new_ptr(ShowerParticle(partons_[0],false )));
    newqin->set5Momentum(ComptonMomenta_[0]);
    if(ComptonISFS_) {
      ShowerParticlePtr newspace(new_ptr(ShowerParticle(partons_[0],false)));
      newspace->set5Momentum(ComptonMomenta_[0]-ComptonMomenta_[2]);
      HardBranchingPtr spaceBranch(new_ptr(HardBranching(newqin,SudakovPtr(),
							 HardBranchingPtr(),
							 HardBranching::Incoming)));
      HardBranchingPtr offBranch(new_ptr(HardBranching(newspace,SudakovPtr(),
						       spaceBranch,
						       HardBranching::Incoming)));
      spaceBranch->addChild(offBranch);
      HardBranchingPtr g(new_ptr(HardBranching(newg,SudakovPtr(),spaceBranch,
					       HardBranching::Outgoing)));
      spaceBranch->addChild(g);
      spaceBranch->type(offBranch->branchingParticle()->id()>0 ? 
			ShowerPartnerType::QCDColourLine : ShowerPartnerType::QCDAntiColourLine);
      HardBranchingPtr outBranch(new_ptr(HardBranching(newqout,SudakovPtr(),
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
      spaceBranchings.push_back(spaceBranch);
      allBranchings.push_back(offBranch);
      allBranchings.push_back(outBranch);
      ColinePtr newin(new_ptr(ColourLine())),newout(new_ptr(ColourLine()));
      newin ->addColoured(newqin  ,partons_[0]->id()<0);
      newin ->addColoured(newg    ,partons_[0]->id()<0);
      newout->addColoured(newspace,partons_[0]->id()<0);
      newout->addColoured(newqout ,partons_[1]->id()<0);
      newout->addColoured(newg    ,partons_[1]->id()>0);
    }
    else {
      ShowerParticlePtr newtime(new_ptr(ShowerParticle(partons_[1],true)));
      newtime->set5Momentum(ComptonMomenta_[1]+ComptonMomenta_[2]);
      HardBranchingPtr spaceBranch(new_ptr(HardBranching(newqin,SudakovPtr(),
							 HardBranchingPtr(),
							 HardBranching::Incoming)));
      HardBranchingPtr offBranch(new_ptr(HardBranching(newtime,SudakovPtr(),
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
      HardBranchingPtr g(new_ptr(HardBranching(newg,SudakovPtr(),offBranch,
					       HardBranching::Outgoing)));
      HardBranchingPtr outBranch(new_ptr(HardBranching(newqout,SudakovPtr(),offBranch,
						       HardBranching::Outgoing)));
      offBranch->addChild(outBranch);
      offBranch->addChild(g);
      offBranch->type(offBranch->branchingParticle()->id()>0 ? 
		       ShowerPartnerType::QCDColourLine : ShowerPartnerType::QCDAntiColourLine);
      spaceBranchings.push_back(spaceBranch);
      allBranchings.push_back(spaceBranch);
      allBranchings.push_back(offBranch);	 
      ColinePtr newin(new_ptr(ColourLine())),newout(new_ptr(ColourLine()));
      newin ->addColoured(newqin ,newqin->dataPtr()->iColour()!=PDT::Colour3);
      newin ->addColoured(newtime,newqin->dataPtr()->iColour()!=PDT::Colour3);
      newin ->addColoured(newg   ,newqin->dataPtr()->iColour()!=PDT::Colour3);
      newout->addColoured(newg   ,newqin->dataPtr()->iColour()==PDT::Colour3);
      newout->addColoured(newqout,newqin->dataPtr()->iColour()!=PDT::Colour3);
    }
  }
  // BGF hardest
  else {
    rot_.invert();
    for(unsigned int ix=0;ix<BGFMomenta_.size();++ix) {
      BGFMomenta_[ix].transform(rot_);
    }
    ShowerParticlePtr newq   (new_ptr(ShowerParticle(partons_[1],true)));
    newq->set5Momentum(BGFMomenta_[1]);
    ShowerParticlePtr newqbar(new_ptr(ShowerParticle(partons_[0]->CC(),true)));
    newqbar->set5Momentum(BGFMomenta_[2]);
    ShowerParticlePtr newg   (new_ptr(ShowerParticle(gluon_,false)));
    newg->set5Momentum(BGFMomenta_[0]);
    ShowerParticlePtr newspace(new_ptr(ShowerParticle(partons_[0],false)));
    newspace->set5Momentum(BGFMomenta_[0]-BGFMomenta_[2]);
    HardBranchingPtr spaceBranch(new_ptr(HardBranching(newg,SudakovPtr(),HardBranchingPtr(),
						       HardBranching::Incoming)));
    HardBranchingPtr offBranch(new_ptr(HardBranching(newspace,SudakovPtr(),spaceBranch,
						     HardBranching::Incoming)));
    HardBranchingPtr qbar(new_ptr(HardBranching(newqbar,SudakovPtr(),spaceBranch,
						HardBranching::Outgoing)));
    spaceBranch->addChild(offBranch);
    spaceBranch->addChild(qbar);
    spaceBranch->type(offBranch->branchingParticle()->id()>0 ? 
		     ShowerPartnerType::QCDColourLine : ShowerPartnerType::QCDAntiColourLine);
    HardBranchingPtr outBranch(new_ptr(HardBranching(newq,SudakovPtr(),
						     HardBranchingPtr(),
						     HardBranching::Outgoing)));
    spaceBranchings.push_back(spaceBranch);
    allBranchings.push_back(offBranch);
    allBranchings.push_back(outBranch); 	
    ColinePtr newin(new_ptr(ColourLine())),newout(new_ptr(ColourLine()));
    newout->addColoured(newspace,newspace->dataPtr()->iColour()!=PDT::Colour3);
    newout->addColoured(newq    ,newspace->dataPtr()->iColour()!=PDT::Colour3);
    newout->addColoured(newg    ,newspace->dataPtr()->iColour()!=PDT::Colour3);
    newin ->addColoured(newg    ,newspace->dataPtr()->iColour()==PDT::Colour3);
    newin ->addColoured(newqbar ,newspace->dataPtr()->iColour()==PDT::Colour3);
  }
  allBranchings[2]->colourPartner(allBranchings[3]);
  allBranchings[3]->colourPartner(allBranchings[2]);
  HardTreePtr newTree(new_ptr(HardTree(allBranchings,spaceBranchings,
				       ShowerInteraction::QCD)));
  // Set the maximum pt for all other emissions and connect hard and shower tree
  Energy pT = isCompton ? pTCompton_ : pTBGF_;
  // incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    // set maximum pT
    if(QuarkMatcher::Check(cit->first->progenitor()->data()))
      cit->first->maximumpT(pT,ShowerInteraction::QCD);
    for(set<HardBranchingPtr>::iterator cjt=newTree->branchings().begin();
	cjt!=newTree->branchings().end();++cjt) {
      if(!(*cjt)->branchingParticle()->isFinalState()&&
	 (*cjt)->branchingParticle()->id()==cit->first->progenitor()->id()) {
	newTree->connect(cit->first->progenitor(),*cjt);
	tPPtr beam =cit->first->original();
	if(!beam->parents().empty()) beam=beam->parents()[0];
	(*cjt)->beam(beam);
	HardBranchingPtr parent=(*cjt)->parent();
	while(parent) {
	  parent->beam(beam);
	  parent=parent->parent();
	};
      }
    }
  }
  // outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    // set maximum pT
    if(QuarkMatcher::Check(cit->first->progenitor()->data()))
      cit->first->maximumpT(pT,ShowerInteraction::QCD);
    for(set<HardBranchingPtr>::iterator cjt=newTree->branchings().begin();
	cjt!=newTree->branchings().end();++cjt) {
      if((*cjt)->branchingParticle()->isFinalState()&&
	 (*cjt)->branchingParticle()->id()==cit->first->progenitor()->id()) {
	newTree->connect(cit->first->progenitor(),*cjt);
      }
    }
  }
  return newTree;
}

void DISBase::generateCompton() {
  // maximum value of the xT
  double xT = sqrt((1.-xB_)/xB_);
  double xTMin = 2.*pTmin_/sqrt(q2_);
  double zp;
  // prefactor
  double a = alpha_->overestimateValue()*comptonWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.);
  vector<double> azicoeff;
  do {
    wgt = 0.;
    // intergration variables dxT/xT^3
    xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
    // zp
    zp = UseRandom::rnd();
    xp = 1./(1.+0.25*sqr(xT)/zp/(1.-zp));
    // check allowed
    if(xp<xB_||xp>1.) continue;
    // phase-space piece of the weight
    wgt = 8.*(1.-xp)*zp/comptonWeight_;
    // PDF piece of the weight
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1.);
    wgt *= pdf_->xfx(beam_,partons_[0],scale,xB_/xp)/
           pdf_->xfx(beam_,partons_[0],q2_  ,xB_);
    // me piece of the weight
    double x2 = 1.-(1.-zp)/xp;
    azicoeff = ComptonME(xp,x2,xT,false);
    wgt *= 4./3.*alpha_->ratio(0.25*q2_*sqr(xT))*(azicoeff[0]+0.5*azicoeff[2]);
    if(wgt>1.||wgt<0.) {
      ostringstream wstring;
      wstring << "DISBase::generateCompton() "
	      << "Weight greater than one or less than zero"
	      << "wgt = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
  }
  while(xT>xTMin&&UseRandom::rnd()>wgt);
  if(xT<=xTMin) {
    pTCompton_=-GeV;
    return;
  }
  // generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi));
    phiwgt = azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi);
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in DISMECorrection"
				  << "::generateCompton() to"
				  << " generate phi" << Exception::eventerror;
  // momenta for the configuration
  Energy Q(sqrt(q2_));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		       -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  pTCompton_ = 0.5*Q*xT;
  ComptonMomenta_.resize(3);
  ComptonMomenta_[0] = p0;
  ComptonMomenta_[1] = p1;
  ComptonMomenta_[2] = p2;
  ComptonISFS_ = zp>xp;
//   debuggingMatrixElement(false,p0,p1,p2,gluon_,pl_[0],pl_[1],pq_[0],pq_[1],
// 			 leptons_[0],leptons_[1],
// 			 partons_[0],partons_[1],
// 			 q2_,phi,x2,x3,xT,zp,xp,azicoeff,false);
}

void DISBase::generateBGF() {
  // maximum value of the xT
  double xT = (1.-xB_)/xB_;
  double xTMin = 2.*max(pTmin_,pTCompton_)/sqrt(q2_);
  double zp;
  // prefactor
  double a = alpha_->overestimateValue()*BGFWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.);
  vector<double> azicoeff;
  do {
    wgt = 0.;
    // intergration variables dxT/xT^3
    xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
    // zp
    zp = UseRandom::rnd();
    xp = 1./(1.+0.25*sqr(xT)/zp/(1.-zp));
    // check allowed
    if(xp<xB_||xp>1.) continue;
    // phase-space piece of the weight
    wgt = 8.*sqr(1.-xp)*zp/BGFWeight_;
    // PDF piece of the weight
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1.);
    wgt *= pdf_->xfx(beam_,gluon_     ,scale,xB_/xp)/
           pdf_->xfx(beam_,partons_[0],q2_  ,xB_);
    // me piece of the weight
    double x1 = -1./xp;
    double x2 = 1.-(1.-zp)/xp;
    double x3 = 2.+x1-x2;
    azicoeff = BGFME(xp,x2,x3,xT,false);
    wgt *= 0.5*alpha_->ratio(0.25*q2_*sqr(xT))*
      (azicoeff[0]+0.5*azicoeff[2]);
    if(wgt>1.||wgt<0.) {
      ostringstream wstring;
      wstring << "DISBase::generateBGF() "
	      << "Weight greater than one or less than zero"
	      << "wgt = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
  }
  while(xT>xTMin&&UseRandom::rnd()>wgt);
  if(xT<=xTMin) {
    pTBGF_=-GeV;
    return;
  }
  // generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi));
    phiwgt = azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi);
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in DISMECorrection"
				  << "::generateBGF() to"
				  << " generate phi" << Exception::eventerror;
  // momenta for the configuration
  Energy Q(sqrt(q2_));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		       -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  pTBGF_=0.5*Q*xT;
  BGFMomenta_.resize(3);
  BGFMomenta_[0]=p0;
  BGFMomenta_[1]=p1;
  BGFMomenta_[2]=p2;
//   debuggingMatrixElement(true,p0,p1,p2,gluon_,pl_[0],pl_[1],pq_[0],pq_[1],
// 			 leptons_[0],leptons_[1],
// 			 partons_[0],partons_[1],
// 			 q2_,phi,x2,x3,xT,zp,xp,azicoeff,false);
}

int DISBase::nDim() const {
  return HwMEBase::nDim() + (contrib_>0 ? 1 : 0 );
}

bool DISBase::generateKinematics(const double * r) {
  // Born kinematics
  if(!HwMEBase::generateKinematics(r)) return false;
  if(contrib_!=0) {
    // hadron and momentum fraction
    if(HadronMatcher::Check(*lastParticles().first->dataPtr())) {
      hadron_ = dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr());
      xB_ = lastX1();
    }
    else {
      hadron_ = dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr());
      xB_ = lastX2();
    }
    // Q2
    q2_ = -(meMomenta()[0]-meMomenta()[2]).m2();
    // xp
    int ndim=nDim();
    double rhomin = pow(1.-xB_,1.-power_); 
    double rho = r[ndim-1]*rhomin;
    xp_ = 1.-pow(rho,1./(1.-power_));
    jac_ = rhomin/(1.-power_)*pow(1.-xp_,power_);
    jacobian(jacobian()*jac_);
  }
  return true; 
}

Energy2 DISBase::scale() const {
  return scaleOpt_ == 1 ? 
    -sqr(scaleFact_)*tHat() : sqr(scaleFact_*muF_);
}

CrossSection DISBase::dSigHatDR() const {
  return NLOWeight()*HwMEBase::dSigHatDR();
}

double DISBase::NLOWeight() const {
  // If only leading order is required return 1:
  if(contrib_==0) return 1.;
  // scale and prefactors
  Energy2 mu2(scale());
  double aS = SM().alphaS(mu2);
  double CFfact = 4./3.*aS/Constants::twopi;
  double TRfact = 1./2.*aS/Constants::twopi;
  // LO + dipole subtracted virtual + collinear quark bit with LO pdf
  double virt = 1.+CFfact*(-4.5-1./3.*sqr(Constants::pi)+1.5*log(q2_/mu2/(1.-xB_))
			   +2.*log(1.-xB_)*log(q2_/mu2)+sqr(log(1.-xB_)));
  virt /= jac_;
  // PDF from leading-order
  double loPDF = hadron_->pdf()->xfx(hadron_,mePartonData()[1],mu2,xB_)/xB_;
  // NLO gluon PDF
  tcPDPtr gluon = getParticleData(ParticleID::g);
  double gPDF   = hadron_->pdf()->xfx(hadron_,gluon,mu2,xB_/xp_)*xp_/xB_;
  // NLO quark PDF
  double qPDF   = hadron_->pdf()->xfx(hadron_,mePartonData()[1],mu2,xB_/xp_)*xp_/xB_;
  // collinear counterterms
  // gluon
  double collg = 
    TRfact/xp_*gPDF*(2.*xp_*(1.-xp_)+(sqr(xp_)+sqr(1.-xp_))*log((1.-xp_)*q2_/xp_/mu2));
  // quark
  double collq = 
    CFfact/xp_*qPDF*(1-xp_-2./(1.-xp_)*log(xp_)-(1.+xp_)*log((1.-xp_)/xp_*q2_/mu2))+
    CFfact/xp_*(qPDF-xp_*loPDF)*(2./(1.-xp_)*log(q2_*(1.-xp_)/mu2)-1.5/(1.-xp_));
  // calculate the A coefficient for the real pieces
  double a(A(mePartonData()[0],mePartonData()[2],
	     mePartonData()[1],mePartonData()[3],q2_));
  // cacluate lepton kinematic variables
  Lorentz5Momentum q = meMomenta()[0]-meMomenta()[2];
  double  yB = (q*meMomenta()[1])/(meMomenta()[0]*meMomenta()[1]);
  double l = 2./yB-1.;
  // q -> qg term
  double realq = CFfact/xp_/(1.+a*l+sqr(l))*qPDF/loPDF*
    (2.+2.*sqr(l)-xp_+3.*xp_*sqr(l)+a*l*(2.*xp_+1.));
  // g -> q qbar term
  double realg =-TRfact/xp_/(1.+a*l+sqr(l))*gPDF/loPDF*
    ((1.+sqr(l)+2.*(1.-3.*sqr(l))*xp_*(1.-xp_))
     +2.*a*l*(1.-2.*xp_*(1.-xp_)));
  // return the full result
  double wgt = virt+((collq+collg)/loPDF+realq+realg); 
  //   double f2g = gPDF/xp_*TRfact*((sqr(1-xp_)+sqr(xp_))*log((1-xp_)/xp_)+
  // 				8*xp_*(1.-xp_)-1.);
  //   double f2q = 
  //     loPDF/jac_*(1.+CFfact*(-1.5*log(1.-xB_)+sqr(log(1.-xB_))
  // 			   -sqr(Constants::pi)/3.-4.5))
  //     +qPDF            *CFfact/xp_*(3.+2.*xp_-(1.+xp_)*log(1.-xp_)
  // 				  -(1.+sqr(xp_))/(1.-xp_)*log(xp_))
  //     +(qPDF-xp_*loPDF)*CFfact/xp_*(2.*log(1.-xp_)/(1.-xp_)-1.5/(1.-xp_));
  //   double wgt = (f2g+f2q)/loPDF;
  return contrib_ == 1 ? max(0.,wgt) : max(0.,-wgt);
}
