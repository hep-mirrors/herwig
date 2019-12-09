// -*- C++ -*-
//
// SMHiggsGGHiggsPPDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsGGHiggsPPDecayer class.
//

#include "SMHiggsGGHiggsPPDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "Herwig/Utilities/HiggsLoopFunctions.h"

using namespace Herwig;
using namespace Herwig::HiggsLoopFunctions;
using namespace ThePEG::Helicity;

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<SMHiggsGGHiggsPPDecayer,PerturbativeDecayer>
describeHerwigSMHiggsGGHiggsPPDecayer("Herwig::SMHiggsGGHiggsPPDecayer",
				      "HwPerturbativeHiggsDecay.so");

bool SMHiggsGGHiggsPPDecayer::accept(tcPDPtr parent,
				       const tPDVector & children) const {
  int idp = parent->id();
  int id0 = children[0]->id();
  int id1 = children[1]->id();
  if((idp == ParticleID::h0 && 
      id0 == ParticleID::g     && id1 == ParticleID::g) ||
     (idp == ParticleID::h0 && 
      id0 == ParticleID::gamma && id1 == ParticleID::gamma)||
     (idp == ParticleID::h0 && 
      id0 == ParticleID::Z0 && id1 == ParticleID::gamma)||
     (idp == ParticleID::h0 && 
      id0 == ParticleID::gamma && id1 == ParticleID::Z0)) 
    return true;
  else
    return false;
}

ParticleVector SMHiggsGGHiggsPPDecayer::decay(const Particle & parent,
					      const tPDVector & children) const {
  int imode(2);
  if(children[0]->id() == ParticleID::gamma && 
     children[1]->id() == ParticleID::gamma)
    imode = 1;
  else if(children[0]->id() ==ParticleID::g)
    imode = 0;
  ParticleVector out(generate(true,false,imode,parent));
  //colour flow
  if(children[0]->id() == ParticleID::g &&
     children[1]->id() == ParticleID::g) {
    out[0]->colourNeighbour(out[1]);
    out[0]->antiColourNeighbour(out[1]);
  }
  return out;
}

void SMHiggsGGHiggsPPDecayer::persistentOutput(PersistentOStream & os) const {
  os << _hggvertex << _hppvertex << _hzpvertex  << _h0wgt
     << _minloop << _maxloop << _massopt;
}

void SMHiggsGGHiggsPPDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _hggvertex >> _hppvertex >> _hzpvertex >> _h0wgt
     >> _minloop >> _maxloop >> _massopt;
}

void SMHiggsGGHiggsPPDecayer::Init() {

  static ClassDocumentation<SMHiggsGGHiggsPPDecayer> documentation
    ("This is an implentation of h0->gg or h0->gamma,gamma "
     "decayer using the SMHGGVertex.");
  
  static Reference<SMHiggsGGHiggsPPDecayer,AbstractVVSVertex> 
    interfaceSMHGGVertex
    ("SMHGGVertex",
     "Pointer to SMHGGVertex",
     &SMHiggsGGHiggsPPDecayer::_hggvertex, false, false, true, 
     false, false);
  
  static Reference<SMHiggsGGHiggsPPDecayer,AbstractVVSVertex> 
    interfaceSMHPPVertex
    ("SMHPPVertex",
     "Pointer to SMHPPVertex",
     &SMHiggsGGHiggsPPDecayer::_hppvertex, false, false, true, 
     false, false);
  
  static Reference<SMHiggsGGHiggsPPDecayer,AbstractVVSVertex> 
    interfaceSMHZPVertex
    ("SMHZPVertex",
     "Pointer to SMHZPVertex",
     &SMHiggsGGHiggsPPDecayer::_hzpvertex, false, false, true, 
     false, false);
  
  static ParVector<SMHiggsGGHiggsPPDecayer,double> interfaceMaxWeights
    ("MaxWeights",
     "Maximum weights for the various decays",
     &SMHiggsGGHiggsPPDecayer::_h0wgt, 3, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
  static Parameter<SMHiggsGGHiggsPPDecayer,int> interfaceMinimumInLoop
    ("MinimumInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &SMHiggsGGHiggsPPDecayer::_minloop, 6, 4, 6,
     false, false, Interface::limited);

  static Parameter<SMHiggsGGHiggsPPDecayer,int> interfaceMaximumInLoop
    ("MaximumInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &SMHiggsGGHiggsPPDecayer::_maxloop, 6, 4, 6,
     false, false, Interface::limited);

  static Switch<SMHiggsGGHiggsPPDecayer,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the treatment of the masses in the loop diagrams",
     &SMHiggsGGHiggsPPDecayer::_massopt, 0, false, false);
  static SwitchOption interfaceMassOptionFull
    (interfaceMassOption,
     "Full",
     "Include the full mass dependence",
     0);
  static SwitchOption interfaceMassOptionLarge
    (interfaceMassOption,
     "Large",
     "Use the heavy mass limit",
     1);

}

double SMHiggsGGHiggsPPDecayer::me2(const int, 
				    const Particle & part,
				    const ParticleVector & decay,
				    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1)));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    _swave = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(_rho);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					  incoming,true);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::constructSpinInfo(_vwave[ix],decay[ix],
					    outgoing,true,
					    decay[ix]->id()!=ParticleID::Z0);
    return 0.;
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(_vwave[ix],decay[ix],outgoing,decay[ix]->id()!=ParticleID::Z0);
  //Set up decay matrix
  Energy2 scale(sqr(part.mass()));
  unsigned int v1hel,v2hel;
  AbstractVVSVertexPtr vertex;
  unsigned int vstep1(2),vstep2(2);
  double sym(1.);
  if(decay[0]->id() == ParticleID::g &&
     decay[1]->id() == ParticleID::g) {
    vertex = _hggvertex;
    sym = 2.;
  }
  else if(decay[0]->id() == ParticleID::gamma &&
	  decay[1]->id() == ParticleID::gamma) {
    vertex = _hppvertex;
    sym = 2.;
  }
  else if(decay[0]->id() == ParticleID::Z0 &&
	  decay[1]->id() == ParticleID::gamma) {
    vertex = _hzpvertex;
    vstep1 = 1;
  }
  else if(decay[1]->id() == ParticleID::Z0 &&
	  decay[0]->id() == ParticleID::gamma) {
    vertex = _hzpvertex;
    vstep2 = 1;
  }
  else
    assert(false);
  // loop over the helicities of the outgoing bosons
  for(v1hel = 0;v1hel < 3;v1hel+=vstep1) {
    for(v2hel = 0;v2hel < 3;v2hel+=vstep2) {
      (*ME())(0,v1hel,v2hel) = vertex->evaluate(scale,_vwave[0][v1hel],
						_vwave[1][v2hel],_swave);
    }
  }
  //store matrix element
  double output = ME()->contract(_rho).real()*UnitRemoval::E2/scale;
  //colour factor (N^2 - 1)/4
  if(decay[0]->id() == ParticleID::g) output *= 8.;
  //symmetric final states
  output /= sym;
  // return the answer
  return output;
}

void SMHiggsGGHiggsPPDecayer::doinit() {
  PerturbativeDecayer::doinit();
  // get the width generator for the higgs
  tPDPtr higgs = getParticleData(ParticleID::h0);
  if(_hggvertex) {
    _hggvertex->init();
  }
  else {
    throw InitException() << "SMHiggsGGHiggsPPDecayer::doinit() - " 
			  << "_hggvertex is null";
  }
  if(_hppvertex) {
    _hppvertex->init();
  }
  else {
    throw InitException() << "SMHiggsGGHiggsPPDecayer::doinit() - " 
			  << "_hppvertex is null";
  }
  if(_hzpvertex) {
    _hzpvertex->init();
  }
  else {
    throw InitException() << "SMHiggsGGHiggsZPDecayer::doinit() - " 
			  << "_hzpvertex is null";
  }
  //set up decay modes
  DecayPhaseSpaceModePtr mode;
  tPDVector extpart(3);
  vector<double> wgt(0);
  // glu,glu mode
  extpart[0] = getParticleData(ParticleID::h0);
  extpart[1] = getParticleData(ParticleID::g);
  extpart[2] = getParticleData(ParticleID::g);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  addMode(mode,_h0wgt[0],wgt);
  // gamma,gamma mode
  extpart[1] = getParticleData(ParticleID::gamma);
  extpart[2] = getParticleData(ParticleID::gamma);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  addMode(mode,_h0wgt[1],wgt);
  // Z0,gamma mode
  extpart[1] = getParticleData(ParticleID::Z0);
  extpart[2] = getParticleData(ParticleID::gamma);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  addMode(mode,_h0wgt[2],wgt);
}

void SMHiggsGGHiggsPPDecayer::doinitrun() {
  _hggvertex->initrun();
  _hppvertex->initrun();
  _hzpvertex->initrun();
  PerturbativeDecayer::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      _h0wgt[ix] = mode(ix)->maxWeight();
    }
  }
}

void SMHiggsGGHiggsPPDecayer::dataBaseOutput(ofstream & os,bool header) const {
  if(header) os << "update decayers set parameters=\"";
  // parameters for the PerturbativeDecayer base class
  for(unsigned int ix=0;ix<_h0wgt.size();++ix) {
    os << "newdef " << name() << ":MaxWeights " << ix << " "
	   << _h0wgt[ix] << "\n";
  }
  os << "newdef " << name() << ":SMHGGVertex " << _hggvertex->fullName() << "\n";
  os << "newdef " << name() << ":SMHPPVertex " << _hppvertex->fullName() << "\n";
  os << "newdef " << name() << ":SMHZPVertex " << _hzpvertex->fullName() << "\n";
  PerturbativeDecayer::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" 
		<< fullName() << "\";" << endl;
}

double SMHiggsGGHiggsPPDecayer::matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
						   const ParticleVector & decay3, MEOption,
#ifndef NDEBUG
						   ShowerInteraction inter) {
#else
						   ShowerInteraction ) {
#endif
  assert(inter==ShowerInteraction::QCD);
  // extract partons and LO momentas
  vector<cPDPtr> partons(1,inpart.dataPtr());
  vector<Lorentz5Momentum> lomom(1,inpart.momentum());
  for(unsigned int ix=0;ix<2;++ix) {
    partons.push_back(decay2[ix]->dataPtr());
    lomom.push_back(decay2[ix]->momentum());
  }
  vector<Lorentz5Momentum> realmom(1,inpart.momentum());
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==2) partons.push_back(decay3[ix]->dataPtr());
    realmom.push_back(decay3[ix]->momentum());
  }
  Energy2 scale = sqr(inpart.mass());
  Energy2 lome = loME(inpart.mass());
  double reme = realME(partons,realmom);
  double ratio = reme/lome*scale;
  // // analytic value for mt -> infinity
  // double x1 = 2.*decay3[0]->momentum().t()/inpart.mass();
  // double x2 = 2.*decay3[1]->momentum().t()/inpart.mass();
  // double x3 = 2.*decay3[2]->momentum().t()/inpart.mass();
  // double test = 8.*Constants::pi*3.*(1.+pow(1-x1,4)+pow(1-x2,4)+pow(1-x3,4))
  //   /(1.-x1)/(1.-x2)/(1.-x3);
  // generator()->log() << "TESTING RATIO " << test << " " << ratio << " " << ratio/test << "\n";
  // remember the symmetry factor
  return ratio/3.;
}

double SMHiggsGGHiggsPPDecayer::realME(//const vector<cPDPtr> & partons, 
				       const vector<cPDPtr> &, 
				       const vector<Lorentz5Momentum> & momenta) const {
  // using std::norm;
  // ScalarWaveFunction  hout(momenta[0],partons[0],outgoing);
  // LorentzPolarizationVector g[3][2];
  // // calculate the polarization vectors for the gluons
  // for(unsigned int iw=0;iw<3;++iw) {
  //   VectorWaveFunction gwave(momenta[iw+1],partons[iw+1],outgoing);
  //   for(unsigned int ix=0;ix<2;++ix) {
  //     //if(iw==2) gwave.reset(10);
  //     //else
  //     gwave.reset(2*ix);
  //     g[iw][ix] = gwave.wave();
  //   }
  // }
  Energy2 mh2 = momenta[0].mass2();
  Energy2 s = (momenta[1]+momenta[2]).m2();
  Energy2 t = (momenta[1]+momenta[3]).m2();
  Energy2 u = (momenta[2]+momenta[3]).m2();
  // calculate the loop functions
  Complex A4stu(0.),A2stu(0.),A2tsu(0.),A2ust(0.);
  for(int ix=_minloop;ix<=_maxloop;++ix) {
    // loop functions
    if(_massopt==0) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      A4stu+=A4(s,t,u,mf2);
      A2stu+=A2(s,t,u,mf2);
      A2tsu+=A2(t,s,u,mf2);
      A2ust+=A2(u,s,t,mf2);
    }
    else {
      A4stu=-1./3.;
      A2stu=-sqr(s/mh2)/3.;
      A2tsu=-sqr(t/mh2)/3.;
      A2ust=-sqr(u/mh2)/3.;
    }
  }
  // Complex A3stu=0.5*(A2stu+A2ust+A2tsu-A4stu);
  // // compute the dot products for the matrix element
  // // and polarization vector * momenta
  // Energy2 pdot[3][3];
  // complex<InvEnergy> eps[3][3][2]; 
  // for(unsigned int ig=0;ig<3;++ig) {
  //   for(unsigned int ip=0;ip<3;++ip) {
  //     pdot[ig][ip]=momenta[ig+1]*momenta[ip+1];
  //     for(unsigned int ih=0;ih<2;++ih) {
  // 	if(ig!=ip)
  // 	  eps[ig][ip][ih]=g[ig][ih].dot(momenta[ip+1])/pdot[ig][ip];
  // 	else
  // 	  eps[ig][ip][ih]=ZERO;
  //     }
  //   }
  // }
  // prefactors
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  // Energy3 pre=sqr(mh2)/mw;
  // // compute the matrix element
  // double output(0.);
  // complex<InvEnergy2> wdot[3][3];
  //  for(unsigned int ghel1=0;ghel1<2;++ghel1) {
  //    for(unsigned int ghel2=0;ghel2<2;++ghel2) {
  //      for(unsigned int ghel3=0;ghel3<2;++ghel3) {
  // 	 wdot[0][1]=g[0][ghel1].dot(g[1][ghel2])/pdot[0][1];
  // 	 wdot[0][2]=g[0][ghel1].dot(g[2][ghel3])/pdot[0][2];
  // 	 wdot[1][0]=wdot[0][1];
  // 	 wdot[1][2]=g[1][ghel2].dot(g[2][ghel3])/pdot[1][2];
  // 	 wdot[2][0]=wdot[0][2];
  // 	 wdot[2][1]=wdot[1][2];
  // 	 // last piece
  // 	 Complex diag=pre*A3stu*(eps[0][2][ghel1]*eps[1][0][ghel2]*eps[2][1][ghel3]-
  // 				 eps[0][1][ghel1]*eps[1][2][ghel2]*eps[2][0][ghel3]+
  // 				 (eps[2][0][ghel3]-eps[2][1][ghel3])*wdot[0][1]+
  // 				 (eps[1][2][ghel2]-eps[1][0][ghel2])*wdot[0][2]+
  // 				 (eps[0][1][ghel1]-eps[0][2][ghel1])*wdot[1][2]);
  // 	 // first piece
  // 	 diag+=pre*(+A2stu*(eps[0][1][ghel1]*eps[1][0][ghel2]-wdot[0][1])*
  // 		    (eps[2][0][ghel3]-eps[2][1][ghel3])
  // 		    +A2tsu*(eps[0][2][ghel1]*eps[2][0][ghel3]-wdot[0][2])*
  // 		    (eps[1][2][ghel2]-eps[1][0][ghel2])
  // 		    +A2ust*(eps[1][2][ghel2]*eps[2][1][ghel3]-wdot[1][2])*
  // 		    (eps[0][1][ghel1]-eps[0][2][ghel1]));
  // 	 output+=norm(diag);
  //      }
  //    }
  //  }
  //  // colour factor and what's left of the prefactor
  //  output *= 6.;
   double me=4.*24./s/t/u*pow<4,1>(mh2)/sqr(mw)*
     (norm(A2stu)+norm(A2ust)+norm(A2tsu)+norm(A4stu));
   return me;
}

Energy2 SMHiggsGGHiggsPPDecayer::loME(Energy mh) const {
  Complex loop(0.);
  Energy2 mh2(sqr(mh));
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  for(int ix=_minloop;ix<=_maxloop;++ix) {
    // loop functions
    if(_massopt==0) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      loop += A1(mh2,mf2);
    }
    else {
      loop += 2./3.;
    }
  }
  return 1./Constants::pi*sqr(mh2)/sqr(mw)*norm(loop);
}
