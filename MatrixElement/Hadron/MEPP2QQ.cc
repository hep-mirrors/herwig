// -*- C++ -*-
//
// MEPP2QQ.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2QQ class.
//

#include "MEPP2QQ.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;


MEPP2QQ::MEPP2QQ() : _quarkflavour(6), _process(0),
		     _bottomopt(1), _topopt(1), _maxflavour(5) {
}

void MEPP2QQ::rebind(const TranslationMap & trans) {
  _gggvertex  = trans.translate( _gggvertex);
  _qqgvertex  = trans.translate( _qqgvertex);
  _gluon      = trans.translate( _gluon);
  for(unsigned int ix=0;ix<_quark.size();++ix) 
    _quark[ix]=trans.translate(_quark[ix]);
  for(unsigned int ix=0;ix<_antiquark.size();++ix)
    _antiquark[ix]=trans.translate(_quark[ix]);
  HwMEBase::rebind(trans);
}

IVector MEPP2QQ::getReferences() {
  IVector ret = HwMEBase::getReferences();
  ret.push_back(_gggvertex);
  ret.push_back(_qqgvertex);
  ret.push_back(_gluon);
  for(unsigned int ix=0;ix<_quark.size();++ix)
    ret.push_back(_quark[ix]);
  for(unsigned int ix=0;ix<_antiquark.size();++ix)
    ret.push_back(_antiquark[ix]);
  return ret;
}

void MEPP2QQ::doinit() {
  HwMEBase::doinit();
  // handling of masses
  if(_quarkflavour==6) {
    massOption(vector<unsigned int>(2,_topopt));
  }
  else {
    massOption(vector<unsigned int>(2,_bottomopt));
  }
  // get the vertex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _qqgvertex = hwsm->vertexFFG();
    _gggvertex = hwsm->vertexGGG();
  }
  else throw InitException() 
	 << "Wrong type of StandardModel object in "
	 << "MEPP2QQ::doinit() the Herwig version must be used" 
	 << Exception::runerror;
  // get the particle data objects
  _gluon=getParticleData(ParticleID::g);
  for(int ix=1;ix<=6;++ix) {
    _quark.push_back(    getParticleData( ix));
    _antiquark.push_back(getParticleData(-ix));
  }
}

Energy2 MEPP2QQ::scale() const {
  Energy2 mq2 = max(meMomenta()[2].mass2(),meMomenta()[3].mass2());
  Energy2 s(0.5*sHat()),u(0.5*(uHat()-mq2)),t(0.5*(tHat()-mq2));
  return 4.*s*t*u/(s*s+t*t+u*u);
}

void MEPP2QQ::persistentOutput(PersistentOStream & os) const {
  os << _gggvertex << _qqgvertex << _quarkflavour << _bottomopt << _topopt
     << _process << _gluon << _quark << _antiquark << _maxflavour;
}

void MEPP2QQ::persistentInput(PersistentIStream & is, int) {
  is >> _gggvertex >> _qqgvertex >> _quarkflavour >> _bottomopt >> _topopt
     >> _process >> _gluon >> _quark >> _antiquark >> _maxflavour;
}

unsigned int MEPP2QQ::orderInAlphaS() const {
  return 2;
}

unsigned int MEPP2QQ::orderInAlphaEW() const {
  return 0;
}

ClassDescription<MEPP2QQ> MEPP2QQ::initMEPP2QQ;
// Definition of the static class description member.

void MEPP2QQ::Init() {

  static ClassDocumentation<MEPP2QQ> documentation
    ("The MEPP2QQ class implements the matrix element for"
     " heavy quark pair production");

  static Switch<MEPP2QQ,unsigned int> interfaceQuarkType
    ("QuarkType",
     "The type of quark",
     &MEPP2QQ::_quarkflavour, 6, false, false);
  static SwitchOption interfaceQuarkTypeCharm
    (interfaceQuarkType,
     "Charm",
     "Produce charm-anticharm",
     4);
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

  static Switch<MEPP2QQ,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2QQ::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcessPair
    (interfaceProcess,
     "Pair",
     "Only include quark-antiquark pair production subprocesses",
     1);
  static SwitchOption interfaceProcess1
    (interfaceProcess,
     "gg2qqbar",
     "Include only gg -> q qbar processes",
     2);
  static SwitchOption interfaceProcessqgqg
    (interfaceProcess,
     "qg2qg",
     "Include only q g -> q g processes",
     4);
  static SwitchOption interfaceProcessqbargqbarg
    (interfaceProcess,
     "qbarg2qbarg",
     "Include only qbar g -> qbar g processes",
     5);
  static SwitchOption interfaceProcessqqqq
    (interfaceProcess,
     "qq2qq",
     "Include only q q -> q q processes",
     6);
  static SwitchOption interfaceProcessqbarqbarqbarqbar
    (interfaceProcess,
     "qbarqbar2qbarqbar",
     "Include only qbar qbar -> qbar qbar processes",
     7);
  static SwitchOption interfaceProcessqqbarqqbar
    (interfaceProcess,
     "qqbar2qqbar",
     "Include only q qbar -> q qbar processes",
     8);

  static Switch<MEPP2QQ,unsigned int> interfaceTopMassOption
    ("TopMassOption",
     "Option for the treatment of the top quark mass",
     &MEPP2QQ::_topopt, 1, false, false);
  static SwitchOption interfaceTopMassOptionOnMassShell
    (interfaceTopMassOption,
     "OnMassShell",
     "The top is produced on its mass shell",
     1);
  static SwitchOption interfaceTopMassOption2
    (interfaceTopMassOption,
     "OffShell",
     "The top is generated off-shell using the mass and width generator.",
     2);

  static Switch<MEPP2QQ,unsigned int> interfaceBottomMassOption
    ("CharmBottomMassOption",
     "Option for the treatment of bottom and lighter quark masses",
     &MEPP2QQ::_bottomopt, 1, false, false);
  static SwitchOption interfaceBottomMassOptionOnMassShell
    (interfaceBottomMassOption,
     "OnMassShell",
     "The bottom is produced on its mass shell, this mass used everywhere",
     1);
  static SwitchOption interfaceBottomMassOptionZeroMass
    (interfaceBottomMassOption,
     "ZeroMass",
     "The bottom is generated on mass shell, but zero mass is used in"
     " the matrix element",
     0);

  static Parameter<MEPP2QQ,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the quarks in the process",
     &MEPP2QQ::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);
}

Selector<MEBase::DiagramIndex>
MEPP2QQ::diagrams(const DiagramVector & diags) const {
  // select the diagram, this is easy for us as we have already done it
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if(diags[i]->id()==-int(_diagram)) sel.insert(1.0, i);
    else sel.insert(0., i);
  }
  return sel;
}

double MEPP2QQ::gg2qqbarME(vector<VectorWaveFunction> &g1,
			   vector<VectorWaveFunction> &g2,
			   vector<SpinorBarWaveFunction> & q,
			   vector<SpinorWaveFunction> & qbar,
			   unsigned int iflow) const {
  // scale
  Energy2 mt(scale());
  Energy mass = q[0].mass();
  // matrix element to be stored
  if(iflow!=0) _me.reset(ProductionMatrixElement(PDT::Spin1,PDT::Spin1,
						 PDT::Spin1Half,PDT::Spin1Half));
  // calculate the matrix element
  double output(0.),sumdiag[3]={0.,0.,0.},sumflow[2]={0.,0.};
  Complex diag[3],flow[2];
  VectorWaveFunction interv;
  SpinorWaveFunction inters;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) { 
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      interv=_gggvertex->evaluate(mt,5,_gluon,g1[ihel1],g2[ihel2]);
      for(unsigned int ohel1=0;ohel1<2;++ohel1) { 
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  //first t-channel diagram
	  inters =_qqgvertex->evaluate(mt,3,qbar[ohel2].particle(),
				       qbar[ohel2],g2[ihel2],mass);
	  diag[0]=_qqgvertex->evaluate(mt,inters,q[ohel1],g1[ihel1]);
	  //second t-channel diagram
	  inters =_qqgvertex->evaluate(mt,3,qbar[ohel2].particle(),
				       qbar[ohel2],g1[ihel1],mass);
	  diag[1]=_qqgvertex->evaluate(mt,inters,q[ohel1],g2[ihel2]);
	  // s-channel diagram
	  diag[2]=_qqgvertex->evaluate(mt,qbar[ohel2],q[ohel1],interv);
	  // colour flows
	  flow[0]=diag[0]+diag[2];
	  flow[1]=diag[1]-diag[2];
	  // sums
	  for(unsigned int ix=0;ix<3;++ix) sumdiag[ix] += norm(diag[ix]);
	  for(unsigned int ix=0;ix<2;++ix) sumflow[ix] += norm(flow[ix]);
	  // total
	  output +=real(flow[0]*conj(flow[0])+flow[1]*conj(flow[1])
			-0.25*flow[0]*conj(flow[1]));
	  // store the me if needed
	  if(iflow!=0) _me(2*ihel1,2*ihel2,ohel1,ohel2)=flow[iflow-1];
	}
      }
    }
  }
  // select a colour flow
  _flow=1+UseRandom::rnd2(sumflow[0],sumflow[1]);
  // select a diagram ensuring it is one of those in the selected colour flow
  sumdiag[_flow%2]=0.;
  _diagram=4+UseRandom::rnd3(sumdiag[0],sumdiag[1],sumdiag[2]);
  // final part of colour and spin factors
//    Energy2 mq2=sqr(getParticleData(6)->mass());
//    double tau1=(tHat()-mq2)/sHat(),tau2=(uHat()-mq2)/sHat(),rho=4*mq2/sHat();
//    double test=(1./6./tau1/tau2-3./8.)*sqr(4.*pi*SM().alphaS(mt))*
//      (sqr(tau1)+sqr(tau2)+rho-0.25*sqr(rho)/tau1/tau2);
//    cerr << "testing matrix element " 
//         << output/48./test << "\n";
   return output/48.;
}

double MEPP2QQ::qg2qgME(vector<SpinorWaveFunction> & qin,
			vector<VectorWaveFunction> &g2,
			vector<SpinorBarWaveFunction> & qout,
			vector<VectorWaveFunction> &g4,
			unsigned int iflow) const {
  // scale
  Energy2 mt(scale());
  Energy mass = qout[0].mass();
  // matrix element to be stored
  if(iflow!=0) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
						 PDT::Spin1Half,PDT::Spin1));
  // calculate the matrix element
  double output(0.),sumdiag[3]={0.,0.,0.},sumflow[2]={0.,0.};
  Complex diag[3],flow[2];
  VectorWaveFunction interv;
  SpinorWaveFunction inters,inters2;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) { 
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      inters=_qqgvertex->evaluate(mt,3,qin[ihel1].particle()->CC(),
				  qin[ihel1],g2[ihel2],mass);
      for(unsigned int ohel1=0;ohel1<2;++ohel1) { 
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  // s-channel diagram
	  diag[0]=_qqgvertex->evaluate(mt,inters,qout[ohel1],g4[ohel2]);
	  // first t-channel
	  inters2=_qqgvertex->evaluate(mt,3,qin[ihel1].particle()->CC(),
				       qin[ihel1],g4[ohel2],mass);
	  diag[1]=_qqgvertex->evaluate(mt,inters2,qout[ohel1],g2[ihel2]);
	  // second t-channel
	  interv=_qqgvertex->evaluate(mt,5,_gluon,qin[ihel1],qout[ohel1]);
	  diag[2]=_gggvertex->evaluate(mt,g2[ihel2],g4[ohel2],interv);
	  // colour flows
	  flow[0]=diag[0]-diag[2];
	  flow[1]=diag[1]+diag[2];
	  // sums
	  for(unsigned int ix=0;ix<3;++ix) sumdiag[ix] += norm(diag[ix]);
	  for(unsigned int ix=0;ix<2;++ix) sumflow[ix] += norm(flow[ix]);
	  // total
	  output +=real(flow[0]*conj(flow[0])+flow[1]*conj(flow[1])
			-0.25*flow[0]*conj(flow[1]));
	  // store the me if needed
	  if(iflow!=0) _me(ihel1,2*ihel2,ohel1,2*ohel2)=flow[iflow-1];
	}
      }
    }
  }
  // test code vs me from ESW
  //Energy2 u(uHat()),t(tHat()),s(sHat());
  //double alphas(4.*pi*SM().alphaS(mt));
  //cerr << "testing matrix element "
  //     << 18./output*(-4./9./s/u+1./t/t)*(s*s+u*u)*sqr(alphas) << endl;
  //select a colour flow
  _flow=1+UseRandom::rnd2(sumflow[0],sumflow[1]);
  // select a diagram ensuring it is one of those in the selected colour flow
  sumdiag[_flow%2]=0.;
  _diagram=10+UseRandom::rnd3(sumdiag[0],sumdiag[1],sumdiag[2]);
  // final part of colour and spin factors
  return output/18.;
}

double MEPP2QQ::qbarg2qbargME(vector<SpinorBarWaveFunction> & qin,
			      vector<VectorWaveFunction> &g2,
			      vector<SpinorWaveFunction> & qout,
			      vector<VectorWaveFunction> &g4,
			      unsigned int iflow) const {
  // scale
  Energy2 mt(scale());
  Energy mass = qout[0].mass();
  // matrix element to be stored
  if(iflow!=0) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
						 PDT::Spin1Half,PDT::Spin1));
  // calculate the matrix element
  double output(0.),sumdiag[3]={0.,0.,0.},sumflow[2]={0.,0.};
  Complex diag[3],flow[2];
  VectorWaveFunction interv;
  SpinorBarWaveFunction inters,inters2;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) { 
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      inters=_qqgvertex->evaluate(mt,3,qin[ihel1].particle()->CC(),
				  qin[ihel1],g2[ihel2],mass);
      for(unsigned int ohel1=0;ohel1<2;++ohel1) { 
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  // s-channel diagram
	  diag[0]=_qqgvertex->evaluate(mt,qout[ohel1],inters,g4[ohel2]);
	  // first t-channel
	  inters2=_qqgvertex->evaluate(mt,3,qin[ihel1].particle()->CC(),
				       qin[ihel1],g4[ohel2],mass);
	  diag[1]=_qqgvertex->evaluate(mt,qout[ohel1],inters2,g2[ihel2]);
	  // second t-channel
	  interv=_qqgvertex->evaluate(mt,5,_gluon,qout[ohel1],qin[ihel1]);
	  diag[2]=_gggvertex->evaluate(mt,g2[ihel2],g4[ohel2],interv);
	  // colour flows
	  flow[0]=diag[0]+diag[2];
	  flow[1]=diag[1]-diag[2];
	  // sums
	  for(unsigned int ix=0;ix<3;++ix) sumdiag[ix] += norm(diag[ix]);
	  for(unsigned int ix=0;ix<2;++ix) sumflow[ix] += norm(flow[ix]);
	  // total
	  output +=real(flow[0]*conj(flow[0])+flow[1]*conj(flow[1])
			-0.25*flow[0]*conj(flow[1]));
	  // store the me if needed
	  if(iflow!=0) _me(ihel1,2*ihel2,ohel1,2*ohel2)=flow[iflow-1];
	}
      }
    }
  }
  // test code vs me from ESW
  //Energy2 u(uHat()),t(tHat()),s(sHat());
  //double alphas(4.*pi*SM().alphaS(mt));
  //cerr << "testing matrix element "
  //     << 18./output*(-4./9./s/u+1./t/t)*(s*s+u*u)*sqr(alphas) << endl;
  //select a colour flow
  _flow=1+UseRandom::rnd2(sumflow[0],sumflow[1]);
  // select a diagram ensuring it is one of those in the selected colour flow
  sumdiag[_flow%2]=0.;
  _diagram=13+UseRandom::rnd3(sumdiag[0],sumdiag[1],sumdiag[2]);
  // final part of colour and spin factors
  return output/18.;
}

double MEPP2QQ::qq2qqME(vector<SpinorWaveFunction> & q1,
			  vector<SpinorWaveFunction> & q2,
			  vector<SpinorBarWaveFunction> & q3,
			  vector<SpinorBarWaveFunction> & q4,
			  unsigned int iflow) const {
  // identify special case of identical quarks
  bool identical(q1[0].id()==q2[0].id());
  // scale
  Energy2 mt(scale());
  // matrix element to be stored
  if(iflow!=0) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
						 PDT::Spin1Half,PDT::Spin1Half));
  // calculate the matrix element
  double output(0.),sumdiag[2]={0.,0.};
  Complex diag[2];
  VectorWaveFunction interv;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) { 
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) { 
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  // first diagram
	  interv = _qqgvertex->evaluate(mt,5,_gluon,q1[ihel1],q3[ohel1]);
	  diag[0] = _qqgvertex->evaluate(mt,q2[ihel2],q4[ohel2],interv);
	  // second diagram if identical
	  if(identical) {
	    interv = _qqgvertex->evaluate(mt,5,_gluon,q1[ihel1],q4[ohel2]);
	    diag[1]=_qqgvertex->evaluate(mt,q2[ihel2],q3[ohel1],interv);
	  }
	  else diag[1]=0.;
	  // sum of diagrams
	  for(unsigned int ix=0;ix<2;++ix) sumdiag[ix] += norm(diag[ix]);
	  // total
	  output +=real(diag[0]*conj(diag[0])+diag[1]*conj(diag[1])
			+2./3.*diag[0]*conj(diag[1]));
	  // store the me if needed
	  if(iflow!=0) _me(ihel1,ihel2,ohel1,ohel2)=diag[iflow-1];
	}
      }
    }
  }
  // identical particle symmetry factor if needed
  if(identical) output*=0.5;
  // test code vs me from ESW
  //Energy2 u(uHat()),t(tHat()),s(sHat());
  //double alphas(4.*pi*SM().alphaS(mt));
  //if(identical)
  //  {cerr << "testing matrix element A "
  //   << 18./output*0.5*(4./9.*((s*s+u*u)/t/t+(s*s+t*t)/u/u)
  //		       -8./27.*s*s/u/t)*sqr(alphas) << endl;}
  //else
  //  {cerr << "testing matrix element B "
  //	   << 18./output*(4./9.*(s*s+u*u)/t/t)*sqr(alphas) << endl;}
  //select a colour flow
  _flow=1+UseRandom::rnd2(sumdiag[0],sumdiag[1]);
  // select a diagram ensuring it is one of those in the selected colour flow
  sumdiag[_flow%2]=0.;
  _diagram=16+UseRandom::rnd2(sumdiag[0],sumdiag[1]);
  // final part of colour and spin factors
  return output/18.;
}

double MEPP2QQ::qbarqbar2qbarqbarME(vector<SpinorBarWaveFunction> & q1,
				      vector<SpinorBarWaveFunction> & q2,
				      vector<SpinorWaveFunction> & q3,
				      vector<SpinorWaveFunction> & q4,
				      unsigned int iflow) const {
  // identify special case of identical quarks
  bool identical(q1[0].id()==q2[0].id());
  // scale
  Energy2 mt(scale());
  // matrix element to be stored
  if(iflow!=0)
    {_me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
				       PDT::Spin1Half,PDT::Spin1Half));}
  // calculate the matrix element
  double output(0.),sumdiag[2]={0.,0.};
  Complex diag[2];
  VectorWaveFunction interv;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) { 
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) { 
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  // first diagram
	  interv = _qqgvertex->evaluate(mt,5,_gluon,q3[ohel1],q1[ihel1]);
	  diag[0] = _qqgvertex->evaluate(mt,q4[ohel2],q2[ihel2],interv);
	  // second diagram if identical
	  if(identical) {
	    interv = _qqgvertex->evaluate(mt,5,_gluon,q4[ohel2],q1[ihel1]);
	    diag[1]=_qqgvertex->evaluate(mt,q3[ohel1],q2[ihel2],interv);
	  }
	  else diag[1]=0.;
	  // sum of diagrams
	  for(unsigned int ix=0;ix<2;++ix) sumdiag[ix] += norm(diag[ix]);
	  // total
	  output +=real(diag[0]*conj(diag[0])+diag[1]*conj(diag[1])
			+2./3.*diag[0]*conj(diag[1]));
	  // store the me if needed
	  if(iflow!=0) _me(ihel1,ihel2,ohel1,ohel2)=diag[iflow-1];
	}
      }
    }
  }
  // identical particle symmetry factor if needed
  if(identical){output*=0.5;}
  // test code vs me from ESW
  //    Energy2 u(uHat()),t(tHat()),s(sHat());
  //    double alphas(4.*pi*SM().alphaS(mt));
  //    if(identical)
  //      {cerr << "testing matrix element A "
  //       << 18./output*0.5*(4./9.*((s*s+u*u)/t/t+(s*s+t*t)/u/u)
  //    		       -8./27.*s*s/u/t)*sqr(alphas) << endl;}
  //    else
  //      {cerr << "testing matrix element B "
  //    	   << 18./output*(4./9.*(s*s+u*u)/t/t)*sqr(alphas) << endl;}
  //select a colour flow
  _flow=1+UseRandom::rnd2(sumdiag[0],sumdiag[1]);
  // select a diagram ensuring it is one of those in the selected colour flow
  sumdiag[_flow%2]=0.;
  _diagram=18+UseRandom::rnd2(sumdiag[0],sumdiag[1]);
  // final part of colour and spin factors
  return output/18.;
}

double MEPP2QQ::qqbar2qqbarME(vector<SpinorWaveFunction>    & q1,
				vector<SpinorBarWaveFunction> & q2,
				vector<SpinorBarWaveFunction> & q3,
				vector<SpinorWaveFunction>    & q4,
				unsigned int iflow) const {
  // type of process
  bool diagon[2]={q1[0].id()== -q2[0].id(),
		  q1[0].id()== -q3[0].id()};
  // scale
  Energy2 mt(scale());
  // matrix element to be stored
  if(iflow!=0) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
						 PDT::Spin1Half,PDT::Spin1Half));
  // calculate the matrix element
  double output(0.),sumdiag[2]={0.,0.};
  Complex diag[2];
  VectorWaveFunction interv;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) { 
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) { 
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  // first diagram
	  if(diagon[0]) {
	    interv = _qqgvertex->evaluate(mt,5,_gluon,q1[ihel1],q2[ihel2]);
	    diag[0] = _qqgvertex->evaluate(mt,q4[ohel2],q3[ohel1],interv);
	  }
	  else diag[0]=0.;
	  // second diagram
	  if(diagon[1]) {
	    interv = _qqgvertex->evaluate(mt,5,_gluon,q1[ihel1],q3[ohel1]);
	    diag[1]=_qqgvertex->evaluate(mt,q4[ohel2],q2[ihel2],interv);
	  }
	  else diag[1]=0.;
	  // sum of diagrams
	  for(unsigned int ix=0;ix<2;++ix) sumdiag[ix] += norm(diag[ix]);
	  // total
	  output +=real(diag[0]*conj(diag[0])+diag[1]*conj(diag[1])
			+2./3.*diag[0]*conj(diag[1]));
	  // store the me if needed
	  if(iflow!=0) _me(ihel1,ihel2,ohel1,ohel2)=diag[iflow-1];
	}
      }
    }
  }
  //select a colour flow
  _flow=1+UseRandom::rnd2(sumdiag[0],sumdiag[1]);
  // select a diagram ensuring it is one of those in the selected colour flow
  sumdiag[_flow%2]=0.;
  _diagram=20+UseRandom::rnd2(sumdiag[0],sumdiag[1]);
  // final part of colour and spin factors
//   Energy2 mq2=sqr(getParticleData(6)->mass());
//   double tau1=(tHat()-mq2)/sHat(),tau2=(uHat()-mq2)/sHat(),rho=4*mq2/sHat();
//   cerr << "testing matrix element " 
//        << output/18./(4./9.*sqr(4.*pi*SM().alphaS(mt))*(sqr(tau1)+sqr(tau2)+0.5*rho))
//        << "\n";
  return output/18.;
}

void MEPP2QQ::getDiagrams() const {
  // gg -> q qbar subprocesses
  if(_process==0||_process==1||_process==2) {
    // first t-channel
    add(new_ptr((Tree2toNDiagram(3),_gluon,_antiquark[_quarkflavour-1],_gluon,
		 1,_quark[_quarkflavour-1], 2,_antiquark[_quarkflavour-1],-4)));
    // interchange
    add(new_ptr((Tree2toNDiagram(3),_gluon,_antiquark[_quarkflavour-1],_gluon,
		 2,_quark[_quarkflavour-1], 1,_antiquark[_quarkflavour-1],-5)));
    // s-channel
    add(new_ptr((Tree2toNDiagram(2),_gluon,_gluon, 1, _gluon,
		 3,_quark[_quarkflavour-1], 3, _antiquark[_quarkflavour-1], -6)));
  }
  // q g -> q g subprocesses
  if(_process==0||_process==4) {
    // s-channel
    add(new_ptr((Tree2toNDiagram(2),_quark[_quarkflavour-1], _gluon,
		 1, _quark[_quarkflavour-1], 
		 3, _quark[_quarkflavour-1], 3, _gluon,-10)));
    // quark t-channel
    add(new_ptr((Tree2toNDiagram(3),_quark[_quarkflavour-1],
		 _quark[_quarkflavour-1],_gluon,
		 2,_quark[_quarkflavour-1],1,_gluon,-11)));
    // gluon t-channel
    add(new_ptr((Tree2toNDiagram(3),_quark[_quarkflavour-1],_gluon,_gluon,
		 1,_quark[_quarkflavour-1],2,_gluon,-12)));
  }
  // qbar g -> qbar g subprocesses
  if(_process==0||_process==5) {
    // s-channel
    add(new_ptr((Tree2toNDiagram(2),_antiquark[_quarkflavour-1], _gluon,
		 1, _antiquark[_quarkflavour-1],
		 3, _antiquark[_quarkflavour-1], 3, _gluon,-13)));
    // quark t-channel
    add(new_ptr((Tree2toNDiagram(3),_antiquark[_quarkflavour-1],
		 _antiquark[_quarkflavour-1],_gluon,
		 2,_antiquark[_quarkflavour-1],1,_gluon,-14)));
    // gluon t-channel
    add(new_ptr((Tree2toNDiagram(3),_antiquark[_quarkflavour-1],_gluon,_gluon,
		 1,_antiquark[_quarkflavour-1],2,_gluon,-15)));
  }
  // processes involving two quark lines
  for(unsigned int iy=0;iy<_maxflavour;++iy) {
    // q q -> q q subprocesses
    if(_process==0||_process==6) {
      // gluon t-channel
      add(new_ptr((Tree2toNDiagram(3),_quark[_quarkflavour-1],_gluon,_quark[iy],
		   1,_quark[_quarkflavour-1],2,_quark[iy],-16)));
      // exchange for identical quarks
      if(_quarkflavour-1==iy)
	add(new_ptr((Tree2toNDiagram(3),_quark[_quarkflavour-1],_gluon,_quark[iy],
		     2,_quark[_quarkflavour-1],1,_quark[iy],-17)));
    }
    // qbar qbar -> qbar qbar subprocesses
    if(_process==0||_process==7) {
      // gluon t-channel
      add(new_ptr((Tree2toNDiagram(3),_antiquark[_quarkflavour-1],_gluon,_antiquark[iy],
		   1,_antiquark[_quarkflavour-1],2,_antiquark[iy],-18)));
      // exchange for identical quarks
      if(_quarkflavour-1==iy)
	add(new_ptr((Tree2toNDiagram(3),_antiquark[_quarkflavour-1],_gluon,_antiquark[iy],
		     2,_antiquark[_quarkflavour-1],1,_antiquark[iy],-19)));
    }
    // q qbar -> q qbar
    if(_process==0||_process==1||_process==8) {
      // gluon s-channel
      add(new_ptr((Tree2toNDiagram(2),_quark[iy], _antiquark[iy],
		   1, _gluon, 3, _quark[_quarkflavour-1], 
		   3, _antiquark[_quarkflavour-1],-20)));
      if(_process!=1) {
	// gluon t-channel
	add(new_ptr((Tree2toNDiagram(3),_quark[iy],_gluon,_antiquark[_quarkflavour-1],
		     1,_quark[iy],2,_antiquark[_quarkflavour-1],-21)));
	// gluon t-channel
	add(new_ptr((Tree2toNDiagram(3),_quark[_quarkflavour-1],_gluon,_antiquark[iy],
		     1,_quark[_quarkflavour-1],2,_antiquark[iy],-21)));
      }
    }
  }
}

Selector<const ColourLines *>
MEPP2QQ::colourGeometries(tcDiagPtr diag) const {
  // colour lines for gg to q qbar
  static const ColourLines cggqq[4]={ColourLines("1  4, -1 -2 3, -3 -5"),
				     ColourLines("3  4, -3 -2 1, -1 -5"),
				     ColourLines("2 -1,  1  3 4, -2 -3 -5"),
				     ColourLines("1 -2, -1 -3 -5, 2 3 4")};
  // colour lines for q g to q g
  static const ColourLines cqgqg[4]={ColourLines("1 -2, 2 3 5, 4 -5"),
				     ColourLines("1 5, 3 4,-3 2 -5 "),
				     ColourLines("1 2 -3, 3 5, -5 -2 4"),
				     ColourLines("1 -2 5,3 2 4,-3 -5")};
  // colour lines for qbar g -> qbar g
  static const ColourLines cqbgqbg[4]={ColourLines("-1 2, -2 -3 -5, -4 5"),
				       ColourLines("-1 -5, -3 -4, 3 -2 5"),
				       ColourLines("-1 -2 3, -3 -5, 5 2 -4"),
				       ColourLines("-1 2 -5,-3 -2 -4, 3 5")};
  // colour lines for q q -> q q 
  static const ColourLines cqqqq[2]={ColourLines("1 2 5,3 -2 4"),
				     ColourLines("1 2 4,3 -2 5")};
  // colour lines for qbar qbar -> qbar qbar
  static const ColourLines cqbqbqbqb[2]={ColourLines("-1 -2 -5,-3 2 -4"),
					 ColourLines("-1 -2 -4,-3 2 -5")};
  // colour lines for q qbar -> q qbar
  static const ColourLines cqqbqqb[2]={ColourLines("1 3 4,-2 -3 -5"),
				       ColourLines("1 2 -3,4 -2 -5")};
  // select the colour flow (as all ready picked just insert answer)
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
    // gg -> q qbar subprocess
  case 4: case 5:
    sel.insert(1.0, &cggqq[abs(diag->id())-4]);
    break;
  case 6:
    sel.insert(1.0, &cggqq[1+_flow]);
    break;
    // q g -> q g subprocess
  case 10: case 11:
    sel.insert(1.0, &cqgqg[abs(diag->id())-10]);
    break;
  case 12:
    sel.insert(1.0, &cqgqg[1+_flow]);
    break;
    // q g -> q g subprocess
  case 13: case 14:
    sel.insert(1.0, &cqbgqbg[abs(diag->id())-13]);
    break;
  case 15:
    sel.insert(1.0, &cqbgqbg[1+_flow]);
    break;
    // q q -> q q subprocess
  case 16: case 17:
    sel.insert(1.0, &cqqqq[abs(diag->id())-16]);
    break;
    // qbar qbar -> qbar qbar subprocess
  case 18: case 19:
    sel.insert(1.0, &cqbqbqbqb[abs(diag->id())-18]);
    break;
    // q qbar -> q qbar subprocess
  case 20: case 21:
    sel.insert(1.0, &cqqbqqb[abs(diag->id())-20]);
    break;
  }
  return sel;
}

double MEPP2QQ::me2() const {
  // total matrix element
  double me(0.);
  // gg initiated processes
  if(mePartonData()[0]->id()==ParticleID::g&&mePartonData()[1]->id()==ParticleID::g) {
    VectorWaveFunction      g1w(rescaledMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction      g2w(rescaledMomenta()[1],mePartonData()[1],incoming);
    SpinorBarWaveFunction    qw(rescaledMomenta()[2],mePartonData()[2],outgoing);
    SpinorWaveFunction    qbarw(rescaledMomenta()[3],mePartonData()[3],outgoing);
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
    me=gg2qqbarME(g1,g2,q,qbar,0);
  }
  // quark first processes
  else if(mePartonData()[0]->id()>0) {
    // q g -> q g
    if(mePartonData()[1]->id()==ParticleID::g) {
      SpinorWaveFunction     qinw(rescaledMomenta()[0],mePartonData()[0],incoming);
      VectorWaveFunction      g2w(rescaledMomenta()[1],mePartonData()[1],incoming);
      SpinorBarWaveFunction qoutw(rescaledMomenta()[2],mePartonData()[2],outgoing);
      VectorWaveFunction      g4w(rescaledMomenta()[3],mePartonData()[3],outgoing);
      vector<VectorWaveFunction> g2,g4;
      vector<SpinorWaveFunction> qin;
      vector<SpinorBarWaveFunction> qout;
      for(unsigned int ix=0;ix<2;++ix) {
	qinw.reset(ix);qin.push_back(qinw);
	g2w.reset(2*ix);g2.push_back(g2w);
	qoutw.reset(ix);qout.push_back(qoutw);
	g4w.reset(2*ix);g4.push_back(g4w);
      }
      // calculate the matrix element
      me = qg2qgME(qin,g2,qout,g4,0);
    }
    // q qbar to q qbar 
    else if(mePartonData()[1]->id()<0) {
      SpinorWaveFunction    q1w(rescaledMomenta()[0],mePartonData()[0],incoming);
      SpinorBarWaveFunction q2w(rescaledMomenta()[1],mePartonData()[1],incoming);
      SpinorBarWaveFunction q3w(rescaledMomenta()[2],mePartonData()[2],outgoing);
      SpinorWaveFunction    q4w(rescaledMomenta()[3],mePartonData()[3],outgoing);
      vector<SpinorWaveFunction>    q1,q4;
      vector<SpinorBarWaveFunction> q2,q3;
      for(unsigned int ix=0;ix<2;++ix) {
	q1w.reset(ix);q1.push_back(q1w);
	q2w.reset(ix);q2.push_back(q2w);
	q3w.reset(ix);q3.push_back(q3w);
	q4w.reset(ix);q4.push_back(q4w);
      }
      // calculate the matrix element
      me = qqbar2qqbarME(q1,q2,q3,q4,0);
    }
    // q q -> q q 
    else if(mePartonData()[1]->id()>0) {
      SpinorWaveFunction    q1w(rescaledMomenta()[0],mePartonData()[0],incoming);
      SpinorWaveFunction    q2w(rescaledMomenta()[1],mePartonData()[1],incoming);
      SpinorBarWaveFunction q3w(rescaledMomenta()[2],mePartonData()[2],outgoing);
      SpinorBarWaveFunction q4w(rescaledMomenta()[3],mePartonData()[3],outgoing);
      vector<SpinorWaveFunction>    q1,q2;
      vector<SpinorBarWaveFunction> q3,q4;
      for(unsigned int ix=0;ix<2;++ix) {
	q1w.reset(ix);q1.push_back(q1w);
	q2w.reset(ix);q2.push_back(q2w);
	q3w.reset(ix);q3.push_back(q3w);
	q4w.reset(ix);q4.push_back(q4w);
      }
      // calculate the matrix element
      me = qq2qqME(q1,q2,q3,q4,0);
    }
  }
  // antiquark first processes
  else if(mePartonData()[0]->id()<0) {
    // qbar g -> qbar g
    if(mePartonData()[1]->id()==ParticleID::g) {
      SpinorBarWaveFunction  qinw(rescaledMomenta()[0],mePartonData()[0],incoming);
      VectorWaveFunction      g2w(rescaledMomenta()[1],mePartonData()[1],incoming);
      SpinorWaveFunction    qoutw(rescaledMomenta()[2],mePartonData()[2],outgoing);
      VectorWaveFunction      g4w(rescaledMomenta()[3],mePartonData()[3],outgoing);
      vector<VectorWaveFunction> g2,g4;
      vector<SpinorBarWaveFunction> qin;
      vector<SpinorWaveFunction> qout;
      for(unsigned int ix=0;ix<2;++ix) {
	qinw.reset(ix);qin.push_back(qinw);
	g2w.reset(2*ix);g2.push_back(g2w);
	qoutw.reset(ix);qout.push_back(qoutw);
	g4w.reset(2*ix);g4.push_back(g4w);	      
      }
      // calculate the matrix element
      me = qbarg2qbargME(qin,g2,qout,g4,0);
    }
    // qbar qbar -> qbar qbar
    else if(mePartonData()[1]->id()<0) {
      SpinorBarWaveFunction q1w(rescaledMomenta()[0],mePartonData()[0],incoming);
      SpinorBarWaveFunction q2w(rescaledMomenta()[1],mePartonData()[1],incoming);
      SpinorWaveFunction    q3w(rescaledMomenta()[2],mePartonData()[2],outgoing);
      SpinorWaveFunction    q4w(rescaledMomenta()[3],mePartonData()[3],outgoing);
      vector<SpinorBarWaveFunction> q1,q2;
      vector<SpinorWaveFunction>    q3,q4;
      for(unsigned int ix=0;ix<2;++ix) {
	q1w.reset(ix);q1.push_back(q1w);
	q2w.reset(ix);q2.push_back(q2w);
	q3w.reset(ix);q3.push_back(q3w);
	q4w.reset(ix);q4.push_back(q4w);
      }
      // calculate the matrix element
      me = qbarqbar2qbarqbarME(q1,q2,q3,q4,0);
    }
  }
  else throw Exception() << "Unknown process in MEPP2QQ::me2()" 
			 << Exception::abortnow;
  // return the answer  
  return me;
}

void MEPP2QQ::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // identify the process and calculate the matrix element
  if(hard[0]->id()==ParticleID::g&&hard[1]->id()==ParticleID::g) {
    if(hard[2]->id()<0) swap(hard[2],hard[3]);
    vector<VectorWaveFunction> g1,g2;
    vector<SpinorBarWaveFunction> q;
    vector<SpinorWaveFunction> qbar;
    // off-shell wavefunctions for the spin correlations
    VectorWaveFunction(     g1,hard[0],incoming,false,true,true);
    VectorWaveFunction(     g2,hard[1],incoming,false,true,true);
    SpinorBarWaveFunction(q   ,hard[2],outgoing,true ,true);
    SpinorWaveFunction(   qbar,hard[3],outgoing,true ,true);
    g1[1]=g1[2];g2[1]=g2[2];
    // on-shell for matrix element
    vector<Lorentz5Momentum> momenta;
    cPDVector data;
    for(unsigned int ix=0;ix<4;++ix) {
      momenta.push_back(hard[ix]->momentum());
      data   .push_back(hard[ix]->dataPtr());
    }
    rescaleMomenta(momenta,data);
    VectorWaveFunction      g1w(rescaledMomenta()[0],data[0],incoming);
    VectorWaveFunction      g2w(rescaledMomenta()[1],data[1],incoming);
    SpinorBarWaveFunction    qw(rescaledMomenta()[2],data[2],outgoing);
    SpinorWaveFunction    qbarw(rescaledMomenta()[3],data[3],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      g1w.reset(2*ix); g1  [ix] = g1w  ;
      g2w.reset(2*ix); g2  [ix] = g2w  ;
      qw.reset(ix);    q   [ix] = qw   ;
      qbarw.reset(ix); qbar[ix] = qbarw;
    }
    gg2qqbarME(g1,g2,q,qbar,_flow);
  }
  else if(hard[0]->id()==ParticleID::g||hard[1]->id()==ParticleID::g) {
    if(hard[0]->id()==ParticleID::g) swap(hard[0],hard[1]);
    if(hard[2]->id()==ParticleID::g) swap(hard[2],hard[3]);
    // on-shell for matrix element
    vector<Lorentz5Momentum> momenta;
    cPDVector data;
    for(unsigned int ix=0;ix<4;++ix) {
      momenta.push_back(hard[ix]->momentum());
      data   .push_back(hard[ix]->dataPtr());
    }
    rescaleMomenta(momenta,data);
    vector<VectorWaveFunction> g2,g4;
    VectorWaveFunction(g2,hard[1],incoming,false,true,true);
    VectorWaveFunction(g4,hard[3],outgoing,true ,true,true);
    VectorWaveFunction gin (rescaledMomenta()[1],data[1],incoming);
    VectorWaveFunction gout(rescaledMomenta()[3],data[3],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      gin .reset(2*ix); g2  [ix] = gin ;
      gout.reset(2*ix); g4  [ix] = gout;
    }
    // q g -> q g 
    if(hard[0]->id()>0) {
      vector<SpinorWaveFunction>    q1;
      vector<SpinorBarWaveFunction> q3;
      SpinorWaveFunction(   q1,hard[0],incoming,false,true);
      SpinorBarWaveFunction(q3,hard[2],outgoing,true ,true);
      SpinorWaveFunction     qin (rescaledMomenta()[0],data[0],incoming);
      SpinorBarWaveFunction  qout(rescaledMomenta()[2],data[2],outgoing);
      for(unsigned int ix=0;ix<2;++ix) {
	qin .reset(ix); q1[ix] = qin ;
	qout.reset(ix); q3[ix] = qout;
      }
      qg2qgME(q1,g2,q3,g4,_flow);
    }
    // qbar g -> qbar g
    else {
      vector<SpinorBarWaveFunction> q1;
      vector<SpinorWaveFunction>    q3;
      SpinorBarWaveFunction(q1,hard[0],incoming,false,true);
      SpinorWaveFunction(   q3,hard[2],outgoing,true ,true);
      SpinorBarWaveFunction qin (rescaledMomenta()[0],data[0],incoming);
      SpinorWaveFunction    qout(rescaledMomenta()[2],data[2],outgoing);
      for(unsigned int ix=0;ix<2;++ix) {
	qin .reset(ix); q1[ix] = qin ;
	qout.reset(ix); q3[ix] = qout;
      }
      qbarg2qbargME(q1,g2,q3,g4,_flow);
    }
  }
  else if(hard[0]->id()>0||hard[1]->id()>0) {
    // q q -> q q 
    if(hard[0]->id()>0&&hard[1]->id()>0) {
      if(hard[2]->id()!=hard[0]->id()) swap(hard[2],hard[3]);
      vector<SpinorWaveFunction> q1,q2;
      vector<SpinorBarWaveFunction> q3,q4;
      SpinorWaveFunction(   q1,hard[0],incoming,false,true);
      SpinorWaveFunction(   q2,hard[1],incoming,false,true);
      SpinorBarWaveFunction(q3,hard[2],outgoing,true ,true);
      SpinorBarWaveFunction(q4,hard[3],outgoing,true ,true);
      // on-shell for matrix element
      vector<Lorentz5Momentum> momenta;
      cPDVector data;
      for(unsigned int ix=0;ix<4;++ix) {
	momenta.push_back(hard[ix]->momentum());
	data   .push_back(hard[ix]->dataPtr());
      }
      rescaleMomenta(momenta,data);
      SpinorWaveFunction    q1w(rescaledMomenta()[0],data[0],incoming);
      SpinorWaveFunction    q2w(rescaledMomenta()[1],data[1],incoming);
      SpinorBarWaveFunction q3w(rescaledMomenta()[2],data[2],outgoing);
      SpinorBarWaveFunction q4w(rescaledMomenta()[3],data[3],outgoing);
      for(unsigned int ix=0;ix<2;++ix) {
	q1w.reset(ix);q1[ix] = q1w;
	q2w.reset(ix);q2[ix] = q2w;
	q3w.reset(ix);q3[ix] = q3w;
	q4w.reset(ix);q4[ix] = q4w;
      }
      qq2qqME(q1,q2,q3,q4,_flow);
    }
    // q qbar -> q qbar
    else {
      if(hard[0]->id()<0) swap(hard[0],hard[1]);
      if(hard[2]->id()<0) swap(hard[2],hard[3]);
      vector<SpinorWaveFunction>    q1,q4;
      vector<SpinorBarWaveFunction> q2,q3;
      // off-shell for spin correlations
      SpinorWaveFunction(   q1,hard[0],incoming,false,true);
      SpinorBarWaveFunction(q2,hard[1],incoming,false,true);
      SpinorBarWaveFunction(q3,hard[2],outgoing,true ,true);
      SpinorWaveFunction(   q4,hard[3],outgoing,true ,true);
      // on-shell for matrix element
      vector<Lorentz5Momentum> momenta;
      cPDVector data;
      for(unsigned int ix=0;ix<4;++ix) {
	momenta.push_back(hard[ix]->momentum());
	data   .push_back(hard[ix]->dataPtr());
      }
      rescaleMomenta(momenta,data);
      SpinorWaveFunction    q1w(rescaledMomenta()[0],data[0],incoming);
      SpinorBarWaveFunction q2w(rescaledMomenta()[1],data[1],incoming);
      SpinorBarWaveFunction q3w(rescaledMomenta()[2],data[2],outgoing);
      SpinorWaveFunction    q4w(rescaledMomenta()[3],data[3],outgoing);
      for(unsigned int ix=0;ix<2;++ix) {
	q1w.reset(ix);q1[ix] = q1w;
	q2w.reset(ix);q2[ix] = q2w;
	q3w.reset(ix);q3[ix] = q3w;
	q4w.reset(ix);q4[ix] = q4w;
      }
      qqbar2qqbarME(q1,q2,q3,q4,_flow);
    }
  }
  else if (hard[0]->id()<0&&hard[1]->id()<0) {
    if(hard[2]->id()!=hard[0]->id()) swap(hard[2],hard[3]);
    vector<SpinorBarWaveFunction> q1,q2;
    vector<SpinorWaveFunction> q3,q4;
    SpinorBarWaveFunction(q1,hard[0],incoming,false,true);
    SpinorBarWaveFunction(q2,hard[1],incoming,false,true);
    SpinorWaveFunction(   q3,hard[2],outgoing,true ,true);
    SpinorWaveFunction(   q4,hard[3],outgoing,true ,true);
    // on-shell for matrix element
    vector<Lorentz5Momentum> momenta;
    cPDVector data;
    for(unsigned int ix=0;ix<4;++ix) {
      momenta.push_back(hard[ix]->momentum());
      data   .push_back(hard[ix]->dataPtr());
    }
    rescaleMomenta(momenta,data);
    SpinorBarWaveFunction q1w(rescaledMomenta()[0],data[0],incoming);
    SpinorBarWaveFunction q2w(rescaledMomenta()[1],data[1],incoming);
    SpinorWaveFunction    q3w(rescaledMomenta()[2],data[2],outgoing);
    SpinorWaveFunction    q4w(rescaledMomenta()[3],data[3],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      q1w.reset(ix);q1[ix] = q1w;
      q2w.reset(ix);q2[ix] = q2w;
      q3w.reset(ix);q3[ix] = q3w;
      q4w.reset(ix);q4[ix] = q4w;
    }
    qbarqbar2qbarqbarME(q1,q2,q3,q4,_flow);
  }
  else throw Exception() << "Unknown process in MEPP2QQ::constructVertex()"
			 << Exception::runerror;
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) 
    hard[ix]->spinInfo()->productionVertex(hardvertex);
}
