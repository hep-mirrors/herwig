// -*- C++ -*-
//
// MEee2gZ2qq.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2gZ2qq class.
//

#include "MEee2gZ2qq.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/MatrixElement/General/HardVertex.h"

using namespace Herwig;

void MEee2gZ2qq::doinit() throw(InitException) {
  HwME2to2Base::doinit();
  massOption(true ,_massopt);
  massOption(false,_massopt);
  rescalingOption(3);
  if(_minflav>_maxflav)
    throw InitException() << "The minimum flavour " << _minflav  
			  << "must be lower the than maximum flavour " << _maxflav
			  << " in MEee2gZ2qq::doinit() " 
			  << Exception::runerror;
  // set the particle data objects
  _Z0=getParticleData(ThePEG::ParticleID::Z0);
  _gamma=getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _theFFZVertex = hwsm->vertexFFZ();
    _theFFPVertex = hwsm->vertexFFP();
  }
  else throw InitException() << "Wrong type of StandardModel object in "
			     << "MEee2gZ2qq::doinit() the Herwig++ version must be used" 
			     << Exception::runerror;
}

void MEee2gZ2qq::getDiagrams() const {
  // specific the diagrams
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  tcPDPtr em = getParticleData(ParticleID::eminus);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  // setup the processes
  for(unsigned int i =_minflav;i<=_maxflav;++i) {
    tcPDPtr qk = getParticleData(i);
    tcPDPtr qb = qk->CC();
    add(new_ptr((Tree2toNDiagram(2), em, ep, 1, gamma, 3, qk, 3, qb, -1)));
    add(new_ptr((Tree2toNDiagram(2), em, ep, 1, Z0   , 3, qk, 3, qb, -2)));
  }
}

Energy2 MEee2gZ2qq::scale() const {
  return sHat();
}

unsigned int MEee2gZ2qq::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2gZ2qq::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEee2gZ2qq::diagrams(const DiagramVector & diags) const {
  double lastCont(0.5),lastBW(0.5);
  if ( lastXCombPtr() ) {
    lastCont = meInfo()[0];
    lastBW = meInfo()[1];
  }
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 ) sel.insert(lastCont, i);
    else if ( diags[i]->id() == -2 ) sel.insert(lastBW, i);
  }
  return sel;
}

Selector<const ColourLines *>
MEee2gZ2qq::colourGeometries(tcDiagPtr ) const {
  static const ColourLines c("-5 4");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c);
  return sel;
}

void MEee2gZ2qq::persistentOutput(PersistentOStream & os) const {
  os << _theFFZVertex << _theFFPVertex << _Z0 << _gamma << _minflav 
     << _maxflav << _massopt;
}

void MEee2gZ2qq::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZVertex >> _theFFPVertex >> _Z0 >> _gamma >> _minflav 
     >> _maxflav >> _massopt;
}

ClassDescription<MEee2gZ2qq> MEee2gZ2qq::initMEee2gZ2qq;
// Definition of the static class description member.
void MEee2gZ2qq::Init() {

  static ClassDocumentation<MEee2gZ2qq> documentation
    ("The MEee2gZ2qq class implements the matrix element for e+e- -> q qbar");

  static Parameter<MEee2gZ2qq,unsigned int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The PDG code of the quark with the lowest PDG code to produce.",
     &MEee2gZ2qq::_minflav, 1, 1, 6,
     false, false, Interface::limited);

  static Parameter<MEee2gZ2qq,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The PDG code of the quark with the highest PDG code to produce",
     &MEee2gZ2qq::_maxflav, 5, 1, 6,
     false, false, Interface::limited);

  static Switch<MEee2gZ2qq,unsigned int> interfaceTopMassOption
    ("TopMassOption",
     "Option for the treatment of the top quark mass",
     &MEee2gZ2qq::_massopt, 1, false, false);
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

}

double MEee2gZ2qq::me2() const {
  int ie(0),ip(1),iqk(2),iqb(3);
  // get the order right
  if(mePartonData()[0]->id()!=11) swap(ie,ip);
  if(mePartonData()[2]->id()<0)   swap(iqk,iqb);
  // compute the spinors
  vector<SpinorWaveFunction> fin,aout;
  vector<SpinorBarWaveFunction>  ain,fout;
  SpinorWaveFunction    ein  (rescaledMomenta()[ie ],mePartonData()[ie ],incoming);
  SpinorBarWaveFunction pin  (rescaledMomenta()[ip ],mePartonData()[ip ],incoming);
  SpinorBarWaveFunction qkout(rescaledMomenta()[iqk],mePartonData()[iqk],outgoing);
  SpinorWaveFunction    qbout(rescaledMomenta()[iqb],mePartonData()[iqb],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    ein.reset(ix)  ;fin.push_back( ein  );
    pin.reset(ix)  ;ain.push_back( pin  );
    qkout.reset(ix);fout.push_back(qkout);
    qbout.reset(ix);aout.push_back(qbout);
  }
  // compute the matrix element
  double me,lastCont,lastBW;
  HelicityME(fin,ain,fout,aout,me,lastCont,lastBW);
  // save the components
  DVector save;
  save.push_back(lastCont);
  save.push_back(lastBW);
  meInfo(save);
  // add the QCD K-factor
  int Nf = SM().Nf(scale());
  me *= (1.0 + alphaS()/Constants::pi 
	 + (1.986-0.115*Nf)*sqr(alphaS()/Constants::pi));
  // return the answer
  return me;
}

ProductionMatrixElement MEee2gZ2qq::HelicityME(vector<SpinorWaveFunction>    & fin,
					       vector<SpinorBarWaveFunction> & ain,
					       vector<SpinorBarWaveFunction> & fout,
					       vector<SpinorWaveFunction>    & aout,
					       double & me,
					       double & cont,
					       double & BW ) const {
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
  VectorWaveFunction interZ,interG;
  // temporary storage of the different diagrams
  Complex diag1,diag2;
  // sum over helicities to get the matrix element
  unsigned int inhel1,inhel2,outhel1,outhel2;
  double total[3]={0.,0.};
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      // intermediate Z
      interZ = _theFFZVertex->evaluate(sHat(),1,_Z0,fin[inhel1],ain[inhel2]);
      // intermediate photon
      interG = _theFFPVertex->evaluate(sHat(),1,_gamma,fin[inhel1],ain[inhel2]);
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {		
	  // first the Z exchange diagram
	  diag1 = _theFFZVertex->evaluate(sHat(),aout[outhel2],fout[outhel1],
					  interZ);
	  // then the photon exchange diagram
	  diag2 = _theFFPVertex->evaluate(sHat(),aout[outhel2],fout[outhel1],
					  interG);
	  // add up squares of individual terms
	  total[1] += real(diag1*conj(diag1));
	  total[2] += real(diag2*conj(diag2));
	  // the full thing including interference
	  diag1+=diag2;
	  total[0] += real(diag1*conj(diag1));
	  output(inhel1,inhel2,outhel1,outhel2)=diag1;
	}
      }
    }
  }
  // results
  for(int ix=0;ix<3;++ix){total[ix]*=0.75;}
  cont = total[2];
  BW   = total[1];
  me   = total[0];
  return output;
}

void MEee2gZ2qq::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  if(hard[0]->id()<hard[1]->id()) swap(hard[0],hard[1]);
  if(hard[2]->id()<hard[3]->id()) swap(hard[2],hard[3]);
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  // get wave functions for off-shell momenta for later on
  SpinorWaveFunction(   fin ,hard[0],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
  SpinorBarWaveFunction(fout,hard[2],outgoing,true ,true);
  SpinorWaveFunction(   aout,hard[3],outgoing,true ,true);
  // now rescale the momenta and compute the matrix element with the
  // rescaled momenta for correlations
  vector<Lorentz5Momentum> momenta;
  cPDVector data;
  for(unsigned int ix=0;ix<4;++ix) {
    momenta.push_back(hard[ix]->momentum());
    data   .push_back(hard[ix]->dataPtr());
  }
  rescaleMomenta(momenta,data);
  SpinorWaveFunction    ein  (rescaledMomenta()[0],data[0],incoming);
  SpinorBarWaveFunction pin  (rescaledMomenta()[1],data[1],incoming);
  SpinorBarWaveFunction qkout(rescaledMomenta()[2],data[2],outgoing);
  SpinorWaveFunction    qbout(rescaledMomenta()[3],data[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    ein.reset(ix)  ; fin [ix] = ein  ;
    pin.reset(ix)  ; ain [ix] = pin  ;
    qkout.reset(ix); fout[ix] = qkout;
    qbout.reset(ix); aout[ix] = qbout;
  }
  // calculate the matrix element
  double me,cont,BW;
  ProductionMatrixElement prodme=HelicityME(fin,ain,fout,aout,me,cont,BW);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(prodme);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix)
    dynamic_ptr_cast<SpinfoPtr>(hard[ix]->spinInfo())->setProductionVertex(hardvertex);
}
