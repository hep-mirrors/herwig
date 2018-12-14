// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2QQH class.
//

#include "GeneralQQHiggs.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

GeneralQQHiggs::GeneralQQHiggs() : quarkFlavour_(6), process_(0), shapeOpt_(2),
			       mh_(), wh_(), alpha_(1.1)
{}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GeneralQQHiggs,HwMEBase>
describeHerwigGeneralQQHiggs("Herwig::GeneralQQHiggs", "Herwig.so");

void GeneralQQHiggs::Init() {

  static ClassDocumentation<GeneralQQHiggs> documentation
    ("The GeneralQQHiggs class implements the matrix elements for the "
     "production of the Higgs boson in association with a heavy quark-antiquark pair");

  static Switch<GeneralQQHiggs,unsigned int> interfaceQuarkType
    ("QuarkType",
     "The type of quark",
     &GeneralQQHiggs::quarkFlavour_, 6, false, false);
  static SwitchOption interfaceQuarkTypeBottom
    (interfaceQuarkType,
     "Bottom",
     "Produce bottom-antibottom",
     5);
  static SwitchOption interfaceQuarkTypeTop
    (interfaceQuarkType,
     "Top",
     "Produce top-antitop",
     6);

  static Switch<GeneralQQHiggs,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &GeneralQQHiggs::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcess1
    (interfaceProcess,
     "gg",
     "Include only gg -> QQbarH processes",
     1);
  static SwitchOption interfaceProcessqbarqbarqbarqbar
    (interfaceProcess,
     "qqbar",
     "Include only qbar qbar -> QQbarH processes",
     2);

  static Switch<GeneralQQHiggs,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &GeneralQQHiggs::shapeOpt_, 2, false, false);
  static SwitchOption interfaceStandardShapeFixed
    (interfaceShapeOption,
     "FixedBreitWigner",
     "Breit-Wigner s-channel resonance",
     1);
  static SwitchOption interfaceStandardShapeRunning
    (interfaceShapeOption,
     "MassGenerator",
     "Use the mass generator to give the shape",
     2);
  static SwitchOption interfaceStandardShapeOnShell
    (interfaceShapeOption,
     "OnShell",
     "Produce an on-shell Higgs boson",
     0);

  static Parameter<GeneralQQHiggs,double> interfaceAlpha
    ("Alpha",
     "Power for the generation of the tranverse mass in the pT mapping",
     &GeneralQQHiggs::alpha_, 1.1, 0.0, 10.0,
     false, false, Interface::limited);

}

Energy2 GeneralQQHiggs::scale() const {
  return sHat();
}

int GeneralQQHiggs::nDim() const {
  return 4 + int(shapeOpt_>0);
}

unsigned int GeneralQQHiggs::orderInAlphaS() const {
  return 2;
}

unsigned int GeneralQQHiggs::orderInAlphaEW() const {
  return 1;
}

IBPtr GeneralQQHiggs::clone() const {
  return new_ptr(*this);
}

IBPtr GeneralQQHiggs::fullclone() const {
  return new_ptr(*this);
}

void GeneralQQHiggs::setKinematics() {
  HwMEBase::setKinematics();
}

void GeneralQQHiggs::persistentOutput(PersistentOStream & os) const {
  os << quarkFlavour_ << process_ << shapeOpt_
     << ounit(mh_,GeV) << ounit(wh_,GeV) << hmass_
     << GGGVertex_ << QQGVertex_ << QQHVertex_
     << gluon_ << higgs_ << quark_ << antiquark_
     << alpha_;
}

void GeneralQQHiggs::persistentInput(PersistentIStream & is, int) {
  is >> quarkFlavour_ >> process_ >> shapeOpt_
     >> iunit(mh_,GeV) >> iunit(wh_,GeV) >> hmass_
     >> GGGVertex_ >> QQGVertex_ >> QQHVertex_
     >> gluon_ >> higgs_ >> quark_ >> antiquark_
     >> alpha_;
}

void GeneralQQHiggs::doinit() {
  HwMEBase::doinit();
  // stuff for the higgs mass
  mh_ = higgs_->mass();
  wh_ = higgs_->width();
  if(higgs_->massGenerator()) {
    hmass_=dynamic_ptr_cast<GenericMassGeneratorPtr>(higgs_->massGenerator());
  }
  if(shapeOpt_==2&&!hmass_) 
    throw InitException()
      << "If using the mass generator for the line shape in GeneralQQHiggs::doinit()"
      << "the mass generator must be an instance of the GenericMassGenerator class"
      << Exception::runerror;
  // get the vertex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!hwsm) throw InitException() << "Wrong type of StandardModel object in "
				  << "GeneralQQHiggs::doinit() the Herwig"
				  << " version must be used" 
				  << Exception::runerror;
  GGGVertex_ = hwsm->vertexGGG();
  QQGVertex_ = hwsm->vertexFFG();
  // get the particle data objects
  gluon_=getParticleData(ParticleID::g);
  for(int ix=1;ix<=6;++ix) {
    quark_.push_back(    getParticleData( ix));
    antiquark_.push_back(getParticleData(-ix));
  }
}

bool GeneralQQHiggs::generateKinematics(const double * r) {
  jacobian(1.);
  // CMS energy
  Energy rs = sqrt(sHat());
  // quark mass
  Energy mq(quark_[quarkFlavour_-1]->mass());
  // generate the higgs mass
  Energy mh(mh_);
  if(shapeOpt_!=0) {
    Energy mhmax = min(rs-2.*mq,higgs_->massMax());
    Energy mhmin = max(ZERO    ,higgs_->massMin());
    if(mhmax<=mhmin) return false;
    double rhomin = atan2((sqr(mhmin)-sqr(mh_)), mh_*wh_);
    double rhomax = atan2((sqr(mhmax)-sqr(mh_)), mh_*wh_);
    mh = sqrt(mh_*wh_*tan(rhomin+r[4]*(rhomax-rhomin))+sqr(mh_));
    jacobian(jacobian()*(rhomax-rhomin));
  }
  if(rs<mh+2.*mq) return false;
  // limits for virtual quark mass
  Energy2 mmin(sqr(mq+mh)),mmax(sqr(rs-mq));
  double rhomin,rhomax;
  if(alpha_==0.) {
    rhomin = mmin/sqr(mq);
    rhomax = mmax/sqr(mq);
  }
  else if(alpha_==1.) {
    rhomax = log((mmax-sqr(mq))/sqr(mq));
    rhomin = log((mmin-sqr(mq))/sqr(mq));
  }
  else {
    rhomin = pow((mmax-sqr(mq))/sqr(mq),1.-alpha_);
    rhomax = pow((mmin-sqr(mq))/sqr(mq),1.-alpha_);
    jacobian(jacobian()/(alpha_-1.));
  }
  // branch for mass smoothing
  Energy2 m132,m232;
  Energy p1,p2;
  // first branch
  if(r[1]<=0.5) {
    double rtemp = 2.*r[1];
    double rho = rhomin+rtemp*(rhomax-rhomin);
    if(alpha_==0) 
      m132 = sqr(mq)*rho;
    else if(alpha_==1)
      m132 = sqr(mq)*(exp(rho)+1.);
    else
      m132 = sqr(mq)*(pow(rho,1./(1.-alpha_))+1.);
    Energy m13 = sqrt(m132);
    try {
      p1 = SimplePhaseSpace::getMagnitude(sHat(), m13, mq);
      p2 = SimplePhaseSpace::getMagnitude(m132,mq,mh);
    } catch ( ImpossibleKinematics & e ) {
      return false;
    }
    Energy ptmin = lastCuts().minKT(mePartonData()[3]);
    double ctmin = -1.0, ctmax = 1.0;
    if ( ptmin > ZERO ) {
      double ctm = 1.0 - sqr(ptmin/p1);
      if ( ctm <= 0.0 ) return false;
      ctmin = max(ctmin, -sqrt(ctm));
      ctmax = min(ctmax, sqrt(ctm));
    }
    double cos1 = getCosTheta(ctmin,ctmax,r[0]);
    double sin1(sqrt(1.-sqr(cos1)));
    double phi1 = Constants::twopi*UseRandom::rnd();
    Lorentz5Momentum p13(sin1*p1*cos(phi1),sin1*p1*sin(phi1),cos1*p1,
			 sqrt(sqr(p1)+m132),m13);
    meMomenta()[3].setVect(Momentum3(-sin1*p1*cos(phi1),-sin1*p1*sin(phi1),-cos1*p1));
    meMomenta()[3].setMass(mq);
    meMomenta()[3].rescaleEnergy();
    bool test=Kinematics::twoBodyDecay(p13,mq,mh,-1.+2*r[2],r[3]*Constants::twopi,
				       meMomenta()[2],meMomenta()[4]);
    if(!test) return false;
    m232 = (meMomenta()[3]+meMomenta()[4]).m2();
    double D = 2./(pow(sqr(mq)/(m132-sqr(mq)),alpha_)+
		   pow(sqr(mq)/(m232-sqr(mq)),alpha_));
    jacobian(0.5*jacobian()*rs/m13*sqr(mq)*D*(rhomax-rhomin)/sHat());
  }
  // second branch
  else {
    double rtemp = 2.*(r[1]-0.5);
    double rho = rhomin+rtemp*(rhomax-rhomin);
    if(alpha_==0) 
      m232 = sqr(mq)*rho;
    else if(alpha_==1)
      m232 = sqr(mq)*(exp(rho)+1.);
    else
      m232 = sqr(mq)*(pow(rho,1./(1.-alpha_))+1.);
    Energy m23 = sqrt(m232);
    try {
      p1 = SimplePhaseSpace::getMagnitude(sHat(), m23, mq);
      p2 = SimplePhaseSpace::getMagnitude(m232,mq,mh);
    } catch ( ImpossibleKinematics & e ) {
      return false;
    }
    Energy ptmin = lastCuts().minKT(mePartonData()[2]);
    double ctmin = -1.0, ctmax = 1.0;
    if ( ptmin > ZERO ) {
      double ctm = 1.0 - sqr(ptmin/p1);
      if ( ctm <= 0.0 ) return false;
      ctmin = max(ctmin, -sqrt(ctm));
      ctmax = min(ctmax, sqrt(ctm));
    }
    double cos1 = getCosTheta(ctmin,ctmax,r[0]);
    double sin1(sqrt(1.-sqr(cos1)));
    double phi1 = Constants::twopi*UseRandom::rnd();
    Lorentz5Momentum p23(-sin1*p1*cos(phi1),-sin1*p1*sin(phi1),-cos1*p1,
			 sqrt(sqr(p1)+m232),m23);
    meMomenta()[2].setVect(Momentum3(sin1*p1*cos(phi1),sin1*p1*sin(phi1),cos1*p1));
    meMomenta()[2].setMass(mq);
    meMomenta()[2].rescaleEnergy();
    bool test=Kinematics::twoBodyDecay(p23,mq,mh,-1.+2*r[2],r[3]*Constants::twopi,
				       meMomenta()[3],meMomenta()[4]);
    if(!test) return false;
    m132 = (meMomenta()[2]+meMomenta()[4]).m2();
    double D = 2./(pow(sqr(mq)/(m132-sqr(mq)),alpha_)+
		   pow(sqr(mq)/(m232-sqr(mq)),alpha_));
    jacobian(0.5*jacobian()*rs/m23*sqr(mq)*D*(rhomax-rhomin)/sHat());
  }
  // calculate jacobian
  jacobian(0.125*jacobian()*p1*p2/sHat());
  // check cuts
  vector<LorentzMomentum> out;
  tcPDVector tout;
  for(unsigned int ix=2;ix<5;++ix) {
    out .push_back(meMomenta()[ix]);
    tout.push_back(mePartonData()[ix]);
  }
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}

CrossSection GeneralQQHiggs::dSigHatDR() const {
  using Constants::pi;
  // jacobian factor for the higgs
  InvEnergy2 bwfact(ZERO);
  Energy moff = meMomenta()[4].mass();
  if(shapeOpt_==1) {
    bwfact = mePartonData()[4]->generateWidth(moff)*moff/pi/
      (sqr(sqr(moff)-sqr(mh_))+sqr(mh_*wh_));
  }
  else if(shapeOpt_==2) {
    bwfact = hmass_->BreitWignerWeight(moff);
  }
  double jac1 = shapeOpt_==0 ? 
    1. : double(bwfact*(sqr(sqr(moff)-sqr(mh_))+sqr(mh_*wh_))/(mh_*wh_));
  return sqr(hbarc)*me2()*jacobian()*jac1/sHat()/pow(Constants::twopi,3);
}

void GeneralQQHiggs::getDiagrams() const {
  tPDPtr Q  = quark_[quarkFlavour_-1];
  tPDPtr QB = antiquark_[quarkFlavour_-1];
  // gg -> q qbar h0 subprocesses
  if(process_==0||process_==1) {
    // first t-channel
    add(new_ptr((Tree2toNDiagram(3), gluon_, QB, gluon_,
		 1, Q, 4, Q , 2, QB, 4, higgs_, -1)));
    add(new_ptr((Tree2toNDiagram(4), gluon_, QB, QB, gluon_,
		 1, Q, 3, QB, 2, higgs_, -2)));
    add(new_ptr((Tree2toNDiagram(3),gluon_,QB,gluon_,
		 1, Q, 2, QB, 5, QB, 5, higgs_, -3)));
    // interchange
    add(new_ptr((Tree2toNDiagram(3),gluon_,Q,gluon_,
		 2, Q, 4, Q , 1, QB, 4, higgs_, -4)));
    add(new_ptr((Tree2toNDiagram(4),gluon_,Q,Q,gluon_,
		 3, Q, 1, QB, 2, higgs_, -5)));
    add(new_ptr((Tree2toNDiagram(3),gluon_,Q,gluon_,
		 2, Q, 1, QB, 5, QB, 5, higgs_, -6)));
    // s-channel
    add(new_ptr((Tree2toNDiagram(2),gluon_,gluon_, 1, gluon_,
		 3, Q, 4, Q,  3, QB, 4, higgs_, -7)));
    add(new_ptr((Tree2toNDiagram(2),gluon_,gluon_, 1, gluon_,
		 3,Q, 3, QB, 5, QB, 5, higgs_, -8)));
  }
  // q qbar -> q qbar
  if(process_==0||process_==2) {
    for(unsigned int ix=1;ix<5;++ix) {
      // gluon s-channel
      add(new_ptr((Tree2toNDiagram(2),quark_[ix-1], antiquark_[ix-1],
		   1, gluon_, 3, Q, 4, Q , 3, QB, 4, higgs_,  -9)));
      add(new_ptr((Tree2toNDiagram(2),quark_[ix-1], antiquark_[ix-1],
		   1, gluon_, 3, Q, 3, QB, 5, QB, 5, higgs_, -10)));
    }
  }
}

double GeneralQQHiggs::me2() const {
  // total matrix element
  double me(0.);
  // gg initiated processes
  if(mePartonData()[0]->id()==ParticleID::g) {
    VectorWaveFunction      g1w(meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction      g2w(meMomenta()[1],mePartonData()[1],incoming);
    SpinorBarWaveFunction    qw(meMomenta()[2],mePartonData()[2],outgoing);
    SpinorWaveFunction    qbarw(meMomenta()[3],mePartonData()[3],outgoing);
    ScalarWaveFunction    higgs(meMomenta()[4],mePartonData()[4],1.,outgoing);
    vector<VectorWaveFunction> g1,g2;
    vector<SpinorBarWaveFunction> q;
    vector<SpinorWaveFunction> qbar;
    for(unsigned int ix=0;ix<2;++ix) {
      g1w.reset(2*ix);g1.push_back(g1w);
      g2w.reset(2*ix);g2.push_back(g2w);
      qw.reset(ix);q.push_back(qw);
      qbarw.reset(ix);qbar.push_back(qbarw);
    }
    // calculate the matrix element
    me=ggME(g1,g2,q,qbar,higgs,0);
  }
  // q qbar initiated
  else {
    SpinorWaveFunction    q1w(meMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction q2w(meMomenta()[1],mePartonData()[1],incoming);
    SpinorBarWaveFunction q3w(meMomenta()[2],mePartonData()[2],outgoing);
    SpinorWaveFunction    q4w(meMomenta()[3],mePartonData()[3],outgoing);
    ScalarWaveFunction    higgs(meMomenta()[4],mePartonData()[4],1.,outgoing);
    vector<SpinorWaveFunction>    q1,q4;
    vector<SpinorBarWaveFunction> q2,q3;
    for(unsigned int ix=0;ix<2;++ix) {
      q1w.reset(ix);q1.push_back(q1w);
      q2w.reset(ix);q2.push_back(q2w);
      q3w.reset(ix);q3.push_back(q3w);
      q4w.reset(ix);q4.push_back(q4w);
    }
    // calculate the matrix element
    me = qqME(q1,q2,q3,q4,higgs,0);
  }
  return me*sHat()*UnitRemoval::InvE2;
}

double GeneralQQHiggs::ggME(vector<VectorWaveFunction> &g1,
			  vector<VectorWaveFunction> &g2,
			  vector<SpinorBarWaveFunction> & q,
			  vector<SpinorWaveFunction> & qbar,
			  ScalarWaveFunction & hwave,
			  unsigned int iflow) const {
  // scale
  Energy2 mt(scale());
  Energy mass = q[0].mass();
  // matrix element to be stored
  if(iflow!=0) me_.reset(ProductionMatrixElement(PDT::Spin1,PDT::Spin1,
						 PDT::Spin1Half,PDT::Spin1Half,
						 PDT::Spin0));
  // calculate the matrix element
  double output(0.),sumflow[2]={0.,0.};
  double sumdiag[8]={0.,0.,0.,0.,0.,0.,0.,0.};
  Complex diag[8],flow[2];
  VectorWaveFunction interv;
  SpinorWaveFunction inters,QBoff;
  SpinorBarWaveFunction intersb,Qoff;

  for(unsigned int ihel1=0;ihel1<2;++ihel1) { 
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      interv = GGGVertex_->evaluate(mt,5,gluon_,g1[ihel1],g2[ihel2]);
      for(unsigned int ohel1=0;ohel1<2;++ohel1) { 
	Qoff = QQHVertex_->evaluate(mt,3,q[ohel1].particle(),
				    q[ohel1],hwave,mass);
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  QBoff = QQHVertex_->evaluate(mt,3,qbar[ohel2].particle(),
				       qbar[ohel2],hwave,mass);
	  // 1st diagram
	  inters  = QQGVertex_->evaluate(mt,1,qbar[ohel2].particle(),
					 qbar[ohel2],g2[ihel2],mass);
	  diag[0] = QQGVertex_->evaluate(mt,inters,Qoff,g1[ihel1]);
	  // 2nd diagram
	  intersb = QQGVertex_->evaluate(mt,1,q[ohel1].particle(),
					 q[ohel1],g1[ihel1],mass);
	  diag[1] = QQHVertex_->evaluate(mt,inters,intersb,hwave);
	  // 3rd diagram
	  diag[2] = QQGVertex_->evaluate(mt,QBoff,intersb,g2[ihel2]);
	  // 4th diagram
	  inters  = QQGVertex_->evaluate(mt,1,qbar[ohel2].particle(),
					 qbar[ohel2],g1[ihel1],mass);
	  diag[3] = QQGVertex_->evaluate(mt,inters,Qoff,g2[ihel2]);
	  // 5th diagram
	  intersb = QQGVertex_->evaluate(mt,1,q[ohel1].particle(),
					 q[ohel1],g2[ihel2],mass);
	  diag[4] = QQHVertex_->evaluate(mt,inters,intersb,hwave);
	  // 6th diagram
	  diag[5] = QQGVertex_->evaluate(mt,QBoff,intersb,g1[ihel1]);
	  // 7th diagram
	  diag[6] = QQGVertex_->evaluate(mt,qbar[ohel2],Qoff    ,interv);
	  // 8th diagram
	  diag[7] = QQGVertex_->evaluate(mt,QBoff      ,q[ohel1],interv);
	  // colour flows
	  flow[0]=diag[0]+diag[1]+diag[2]+(diag[6]+diag[7]);
	  flow[1]=diag[3]+diag[4]+diag[5]-(diag[6]+diag[7]);
	  // sums
	  for(unsigned int ix=0;ix<8;++ix) sumdiag[ix] += norm(diag[ix]);
	  for(unsigned int ix=0;ix<2;++ix) sumflow[ix] += norm(flow[ix]);
	  // total
	  output +=real(flow[0]*conj(flow[0])+flow[1]*conj(flow[1])
			-0.25*flow[0]*conj(flow[1]));
	  // store the me if needed
	  if(iflow!=0) me_(2*ihel1,2*ihel2,ohel1,ohel2,0)=flow[iflow-1];
	}
      }
    }
  }
  // select a colour flow
  flow_ = 1 + UseRandom::rnd2(sumflow[0],sumflow[1]);
  if(flow_==1) sumdiag[0]=sumdiag[1]=sumdiag[2]=0.;
  else         sumdiag[3]=sumdiag[4]=sumdiag[5]=0.;
  // select a diagram from that flow
  double prob = UseRandom::rnd();
  for(unsigned int ix=0;ix<8;++ix) {
    if(prob<=sumdiag[ix]) {
      diagram_=1+ix;
      break;
    }
    prob -= sumdiag[ix];
  }				   
  // final part of colour and spin factors
  return output/48.;
}

double GeneralQQHiggs::qqME(vector<SpinorWaveFunction> & q1,
			  vector<SpinorBarWaveFunction> & q2,
			  vector<SpinorBarWaveFunction>    & q3,
			  vector<SpinorWaveFunction>    & q4,
			  ScalarWaveFunction & hwave,
			  unsigned int iflow) const {
  // scale
  Energy2 mt(scale());
  Energy mass = q3[0].mass();
  // matrix element to be stored
  if(iflow!=0) me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
						 PDT::Spin1Half,PDT::Spin1Half,
						 PDT::Spin0));
  // calculate the matrix element
  double output(0.),sumdiag[2]={0.,0.};
  Complex diag[2];
  VectorWaveFunction interv;
  SpinorWaveFunction QBoff;
  SpinorBarWaveFunction Qoff;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) { 
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      interv = QQGVertex_->evaluate(mt,5,gluon_,q1[ihel1],q2[ihel2]);
      for(unsigned int ohel1=0;ohel1<2;++ohel1) { 
	Qoff = QQHVertex_->evaluate(mt,3,q3[ohel1].particle(),
				    q3[ohel1],hwave,mass);
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  QBoff = QQHVertex_->evaluate(mt,3,q4[ohel2].particle(),
				       q4[ohel2],hwave,mass);
	  // 1st diagram
	  diag[0] = QQGVertex_->evaluate(mt,q4[ohel2],Qoff,interv);
	  // 2nd diagram
	  diag[1] = QQGVertex_->evaluate(mt,QBoff,q3[ohel1],interv);
	  // sum of diagrams
	  for(unsigned int ix=0;ix<2;++ix) sumdiag[ix] += norm(diag[ix]);
	  diag[0] += diag[1];
	  output += norm(diag[0]);
	  if(iflow!=0) me_(ihel1,ihel2,ohel1,ohel2,0) = diag[0];
	}
      }
    }
  }
  // only 1 colour flow
  flow_=1;
  // select a diagram
  diagram_ = 9+UseRandom::rnd2(sumdiag[0],sumdiag[1]);
  // final part of colour and spin factors
  return output/18.;
}

Selector<const ColourLines *>
GeneralQQHiggs::colourGeometries(tcDiagPtr diag) const {
  // colour lines for gg -> Q Qbar H
  static const ColourLines cgg[10]=
    {ColourLines("1 4 5, -1 -2  3   , -3 -6   "),
     ColourLines("1 5  , -1 -2 -3  4, -4 -6   "),
     ColourLines("1 4  , -1 -2  3   , -3 -5 -6"),
     ColourLines("3 4 5,  1  2 -3   , -1 -6   "),
     ColourLines("4 5  ,  1  2  3 -4, -1 -6"),
     ColourLines("3 4  ,  1  2 -3   , -1 -5 -6"),
     ColourLines("1 3 4 5, -1  2, -2 -3 -6"),
     ColourLines("2 3 4 5,  1 -2, -1 -3 -6"),
     ColourLines("1 3 4, -1  2, -2 -3 -5 -6"),
     ColourLines("2 3 4,  1 -2, -1 -3 -5 -6")};
  // colour lines for q qbar -> Q Qbar H
  static const ColourLines cqq[2]=
    {ColourLines("1 3 4 5, -2 -3 -6"),
     ColourLines("1 3 4  , -2 -3 -5 -6")};
  // select the colour flow (as all ready picked just insert answer)
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
    // gg -> q qbar subprocess
  case 1: case 2: case 3: case 4: case 5: case 6:
    sel.insert(1.0, &cgg[abs(diag->id())-1]);
    break;
  case 7:
    sel.insert(1.0, &cgg[5 + flow_]);
    break;
  case 8:
    sel.insert(1.0, &cgg[7 + flow_]);
    break;
    // q qbar -> q qbar subprocess
  case 9: case 10:
    sel.insert(1.0, &cqq[abs(diag->id())-9]);
    break;
  }
  return sel;
}

Selector<MEBase::DiagramIndex>
GeneralQQHiggs::diagrams(const DiagramVector & diags) const {
  // select the diagram, this is easy for us as we have already done it
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if(diags[i]->id()==-int(diagram_)) sel.insert(1.0, i);
    else sel.insert(0., i);
  }
  return sel;
}

void GeneralQQHiggs::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  for(unsigned int ix=0;ix<3;++ix) hard.push_back(sub->outgoing()[ix]);
  // identify the process and calculate the matrix element
  if(hard[0]->id()<0)             swap(hard[0],hard[1]);
  if(hard[2]->id()==higgs_->id()) swap(hard[2],hard[4]);
  if(hard[3]->id()==higgs_->id()) swap(hard[3],hard[4]);
  if(hard[2]->id()<0)             swap(hard[2],hard[3]);
  if(hard[2]->id()<0)             swap(hard[2],hard[3]);
  if(hard[0]->id()==ParticleID::g) {
    vector<VectorWaveFunction> g1,g2;
    vector<SpinorBarWaveFunction> q;
    vector<SpinorWaveFunction> qbar;
    // off-shell wavefunctions for the spin correlations
    VectorWaveFunction(     g1,hard[0],incoming,false,true,true);
    VectorWaveFunction(     g2,hard[1],incoming,false,true,true);
    SpinorBarWaveFunction(q   ,hard[2],outgoing,true ,true);
    SpinorWaveFunction(   qbar,hard[3],outgoing,true ,true);
    ScalarWaveFunction hwave(        hard[4],outgoing,true);
    g1[1]=g1[2];g2[1]=g2[2];
    ggME(g1,g2,q,qbar,hwave,flow_);
  }
  // q qbar -> Q Qbar Higgs
  else {
    vector<SpinorWaveFunction>    q1,q4;
    vector<SpinorBarWaveFunction> q2,q3;
    // off-shell for spin correlations
    SpinorWaveFunction(   q1,hard[0],incoming,false,true);
    SpinorBarWaveFunction(q2,hard[1],incoming,false,true);
    SpinorBarWaveFunction(q3,hard[2],outgoing,true ,true);
    SpinorWaveFunction(   q4,hard[3],outgoing,true ,true);
    ScalarWaveFunction hwave(        hard[4],outgoing,true);
    qqME(q1,q2,q3,q4,hwave,flow_);
  }
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me_);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<5;++ix) 
    hard[ix]->spinInfo()->productionVertex(hardvertex);
}

void GeneralQQHiggs::setProcessInfo(unsigned int quark, PDPtr hin,
				    AbstractFFSVertexPtr vertex,
				    unsigned int shapeOpt,
				    unsigned int proc) {
  quarkFlavour_ = quark;
  higgs_ = hin;
  QQHVertex_ = vertex;
  process_ = proc;
  shapeOpt_ = shapeOpt;
}
