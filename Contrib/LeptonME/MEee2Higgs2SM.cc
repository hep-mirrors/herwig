// -*- C++ -*-
//
// MEee2Higgs2SM.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2Higgs2SM class.
//

#include "MEee2Higgs2SM.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

void MEee2Higgs2SM::doinit() {
  ME2to2Base::doinit();
  _h0      = getParticleData(ThePEG::ParticleID::h0);
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _theFFHVertex = hwsm->vertexFFH();
  }
  else {
    throw InitException() << "Must have Herwig++ StandardModel object in "
			  << "MEee2Higgs2SM::doinit()" << Exception::runerror;
  }
}

Energy2 MEee2Higgs2SM::scale() const {
  return sHat();
}

unsigned int MEee2Higgs2SM::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2Higgs2SM::orderInAlphaEW() const {
  return 2;
}

void MEee2Higgs2SM::getDiagrams() const {
  // specific the diagrams
  tcPDPtr ep       = getParticleData(ParticleID::eplus );
  tcPDPtr em       = getParticleData(ParticleID::eminus);
  // outgoing quarks
  for(int i=1;i<=5;++i) {
    if(_allowed==0 || _allowed==1 || _allowed==i+2) {
      tcPDPtr lm = getParticleData(i);
      tcPDPtr lp = lm->CC();
      add(new_ptr((Tree2toNDiagram(2), em, ep, 1, _h0, 3, lm, 3, lp, -1)));
    }
  }
  // outgoing leptons
  for( int i =11;i<=16;i+=2) {
    if(_allowed==0 || _allowed==2 || _allowed==(i+7)/2) {
      tcPDPtr lm = getParticleData(i);
      tcPDPtr lp = lm->CC();
      add(new_ptr((Tree2toNDiagram(2), em, ep, 1, _h0, 3, lm, 3, lp, -1)));
    }
  }
}

double MEee2Higgs2SM::me2() const {
  double aver=0.;
  // get the order right
  int ielectron(0),ipositron(1),ilp(2),ilm(3);
  if(mePartonData()[0]->id()!=11) swap(ielectron,ipositron);
  if(mePartonData()[2]->id()>mePartonData()[3]->id()) swap(ilp,ilm);
  // the arrays for the wavefunction to be passed to the matrix element
  vector<SpinorWaveFunction> fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  for(unsigned int ihel=0;ihel<2;++ihel) {
    fin .push_back(SpinorWaveFunction(meMomenta()[ielectron],
				      mePartonData()[ielectron],ihel,incoming));
    ain .push_back(SpinorBarWaveFunction(meMomenta()[ipositron],
					 mePartonData()[ipositron],ihel,incoming));
    fout.push_back(SpinorBarWaveFunction(meMomenta()[ilm],
					 mePartonData()[ilm],ihel,outgoing));
    aout.push_back(SpinorWaveFunction(meMomenta()[ilp],
				      mePartonData()[ilp],ihel,outgoing));
  }
  HelicityME(fin,ain,fout,aout,aver);
  if(mePartonData()[ilm]->id()<=6) aver*=3.;
  return aver;
}

// the helicity amplitude matrix element
ProductionMatrixElement MEee2Higgs2SM::HelicityME(vector<SpinorWaveFunction> fin,
						  vector<SpinorBarWaveFunction> ain,
						  vector<SpinorBarWaveFunction> fout,
						  vector<SpinorWaveFunction> aout,
						  double & aver) const {
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (u    spinor)
  // 1 incoming antifermion (vbar spinor)
  // for the outgoing       
  // 0 outgoing fermion     (ubar spinor)
  // 1 outgoing antifermion (v    spinor)
  // me to be returned
  ProductionMatrixElement output(PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  // wavefunctions for the intermediate particles
  ScalarWaveFunction interh;
  // temporary storage of the different diagrams
  Complex diag;
  aver=0.;
  // sum over helicities to get the matrix element
  unsigned int inhel1,inhel2,outhel1,outhel2;
  for(inhel1=0;inhel1<2;++inhel1) {	  
    for(inhel2=0;inhel2<2;++inhel2) {
      interh = _theFFHVertex->evaluate(sHat(),1,_h0,fin[inhel1],ain[inhel2]);
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {
	  diag  = _theFFHVertex->evaluate(sHat(),aout[outhel2],
					  fout[outhel1],interh);
	  output(inhel1,inhel2,outhel1,outhel2)=diag;
	  aver +=real(diag*conj(diag));
	}
      }
    }
  }
  return output;
}

Selector<MEBase::DiagramIndex>
MEee2Higgs2SM::diagrams(const DiagramVector &) const {
  Selector<DiagramIndex> sel;sel.insert(1.0, 0);
  return sel;
}
Selector<const ColourLines *>
MEee2Higgs2SM::colourGeometries(tcDiagPtr diag) const {
  static ColourLines neutral ( " " );
  static ColourLines quarks  ( "-5 4");
  Selector<const ColourLines *> sel;
  int id = abs((diag->partons()[2])->id());
  if (id<=6 )
    sel.insert(1.0, &quarks);
  else 
    sel.insert(1.0, &neutral);
  return sel;
}

void MEee2Higgs2SM::persistentOutput(PersistentOStream & os) const {
  os << _theFFHVertex << _h0 << _allowed;
}

void MEee2Higgs2SM::persistentInput(PersistentIStream & is, int) {
  is >> _theFFHVertex >> _h0 >>_allowed;
}

ClassDescription<MEee2Higgs2SM> MEee2Higgs2SM::initMEee2Higgs2SM;
// Definition of the static class description member.

void MEee2Higgs2SM::Init() {

  static ClassDocumentation<MEee2Higgs2SM> documentation
    ("The MEee2Higgs2SM class implements the matrix element for e+e- to"
     " SM particle via Higgs exchnage and is designed to be a process to test"
     " things involving scalar particles. ");

  static Switch<MEee2Higgs2SM,int> interfaceallowed
    ("Allowed",
     "Allowed outgoing particles",
     &MEee2Higgs2SM::_allowed, 0, false, false);
  static SwitchOption interfaceallowedAll
    (interfaceallowed,
     "All",
     "All SM particles allowed",
     0);
  static SwitchOption interfaceallowed1
    (interfaceallowed,
     "Quarks",
     "Only the quarks allowed",
     1);
  static SwitchOption interfaceallowed2
    (interfaceallowed,
     "Leptons",
     "Only the leptons allowed",
     2);
  static SwitchOption interfacealloweddown
    (interfaceallowed,
     "Down",
     "Only d dbar allowed",
     3);
  static SwitchOption interfaceallowedup
    (interfaceallowed,
     "Up",
     "Only u ubar allowed",
     4);
  static SwitchOption interfaceallowedstrange
    (interfaceallowed,
     "Strange",
     "Only s sbar allowed",
     5);
  static SwitchOption interfaceallowedcharm
    (interfaceallowed,
     "Charm",
     "Only c cbar allowed",
     6);
  static SwitchOption interfaceallowedbottom
    (interfaceallowed,
     "Bottom",
     "Only b bbar allowed",
     7);
  static SwitchOption interfaceallowedtop
    (interfaceallowed,
     "Top",
     "Only t tbar allowed",
     8);
  static SwitchOption interfaceallowedelectron
    (interfaceallowed,
     "Electron",
     "Only e+e- allowed",
     9);
  static SwitchOption interfaceallowedMuon
    (interfaceallowed,
     "Muon",
     "Only mu+mu- allowed",
     10);
  static SwitchOption interfaceallowedTau
    (interfaceallowed,
     "Tau",
     "Only tau+tau- allowed",
     11);
}

void MEee2Higgs2SM::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);hard.push_back(sub->outgoing()[1]);
  if(hard[0]->id()<hard[1]->id()) swap(hard[0],hard[1]);
  if(hard[2]->id()<hard[3]->id()) swap(hard[2],hard[3]);
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction(   fin ,hard[0],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
  SpinorBarWaveFunction(fout,hard[2],outgoing,true ,true);
  SpinorWaveFunction(   aout,hard[3],outgoing,true ,true);
  // calculate the matrix element
  double me;
  ProductionMatrixElement prodme=HelicityME(fin,ain,fout,aout,me);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(prodme);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) {
    (hard[ix]->spinInfo())->productionVertex(hardvertex);
  }
}

void MEee2Higgs2SM::rebind(const TranslationMap & trans)
  {
  // dummy = trans.translate(dummy);
  _theFFHVertex = trans.translate(_theFFHVertex);
  _h0        = trans.translate(_h0);
  ME2to2Base::rebind(trans);
}

IVector MEee2Higgs2SM::getReferences() {
  IVector ret = ME2to2Base::getReferences();
  ret.push_back(_theFFHVertex);
  ret.push_back(_h0);
  return ret;
}
