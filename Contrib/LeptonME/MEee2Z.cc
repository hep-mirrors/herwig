// -*- C++ -*-
//
// MEee2Z.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2Z class.
//

#include "MEee2Z.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Cuts/Cuts.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

void MEee2Z::getDiagrams() const {
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  tcPDPtr em = getParticleData(ParticleID::eminus);
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  add(new_ptr((Tree2toNDiagram(2), em, ep, 1, Z0,-1)));
}

Energy2 MEee2Z::scale() const {
  return sHat();
}

int MEee2Z::nDim() const {
  return 1;
}

bool MEee2Z::generateKinematics(const double *) {
  // Here you can use nDim() random numbers in the vector provided
  // to generate the internal kinematics. Note that sHat() has
  // already been given from the outside.
  meMomenta()[2]=meMomenta()[0]+meMomenta()[1];
  meMomenta()[2].rescaleMass();
  jacobian(1.0);
  // check passes all the cuts
  vector<LorentzMomentum> out(1,meMomenta()[2]);
  tcPDVector tout(1,mePartonData()[2]);
  // return true if passes the cuts
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}

double MEee2Z::me2() const {
  double aver=0.;
  // the arrays for the wavefunction to be passed to the matrix element
  vector<SpinorWaveFunction> fin;
  vector<SpinorBarWaveFunction> ain;
  vector<VectorWaveFunction> vin;
  SpinorWaveFunction    fwave(meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction awave(meMomenta()[1],mePartonData()[1],incoming);
  for(unsigned int ihel=0;ihel<2;++ihel) {
    fwave.reset(ihel);fin.push_back(fwave);
    awave.reset(ihel);ain.push_back(awave);
  }
  VectorWaveFunction vwave(meMomenta()[2],mePartonData()[2],outgoing);
  for(unsigned int ihel=0;ihel<3;++ihel) {
    vwave.reset(ihel); vin.push_back(vwave);
  }
  ProductionMatrixElement temp=HelicityME(fin,ain,vin,aver);
  // add the Breit-Wigner factors
  Energy width=mePartonData()[2]->width();
  Energy mass =mePartonData()[2]->mass();
  InvEnergy2 fact = width*mass/(sqr(sHat()-mass*mass)+sqr(mass*width));
  return aver*fact*sHat();
}

CrossSection MEee2Z::dSigHatDR() const {
  return (me2()*jacobian()/sHat())*sqr(hbarc);
}

unsigned int MEee2Z::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2Z::orderInAlphaEW() const {
  return 1;
}

Selector<MEBase::DiagramIndex>
MEee2Z::diagrams(const DiagramVector &) const {
  Selector<DiagramIndex> sel;
  sel.insert(1.0, 0);
  return sel;
}

Selector<const ColourLines *>
MEee2Z::colourGeometries(tcDiagPtr) const {
  static const ColourLines neutral ( " " );
  Selector<const ColourLines *> sel;
  sel.insert(1.,&neutral);
  return sel;
}


void MEee2Z::persistentOutput(PersistentOStream & os) const {
  os << _theFFZVertex;
}

void MEee2Z::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZVertex;
}

ClassDescription<MEee2Z> MEee2Z::initMEee2Z;
// Definition of the static class description member.

void MEee2Z::Init() {

  static ClassDocumentation<MEee2Z> documentation
    ("The MEee2Z class implements the e+e- -> Z as a 1->2 process for testing"
     " of spin correlations etc..");

}

void MEee2Z::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  if(hard[0]->id()<hard[1]->id()) swap(hard[0],hard[1]);
  vector<SpinorWaveFunction>    fin;
  vector<SpinorBarWaveFunction> ain;
  vector<VectorWaveFunction>    vin;
  SpinorWaveFunction(   fin,hard[0],incoming,false,true);
  SpinorBarWaveFunction(ain,hard[1],incoming,false,true);
  VectorWaveFunction   (vin,hard[2],outgoing,true,false,true);
  double dummy;
  ProductionMatrixElement prodme=HelicityME(fin,ain,vin,dummy);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(prodme);
  // set the pointers to and from the vertex
  for(unsigned int ix=0;ix<3;++ix) {
    (hard[ix]->spinInfo())->
      productionVertex(hardvertex);
  }
}

// the helicity amplitude matrix element
ProductionMatrixElement MEee2Z::HelicityME(vector<SpinorWaveFunction> fin,
					   vector<SpinorBarWaveFunction> ain,
					   vector<VectorWaveFunction> vout,
					   double & aver) const
{
  ProductionMatrixElement output(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1);
  Complex product;
  // sum over helicities to get the matrix element
  unsigned int inhel1,inhel2,outhel1;
  double me(0.);
  LorentzPolarizationVector vec;
  Complex ii(0.,1.);
  for(inhel1=0;inhel1<2;++inhel1) {	  
    for(inhel2=0;inhel2<2;++inhel2) {
      for(outhel1=0;outhel1<3;++outhel1) {
	product=_theFFZVertex->evaluate(sHat(),fin[inhel1],ain[inhel2],
					vout[outhel1]);
	output(inhel1,inhel2,outhel1)=product;
	me+=real(product*conj(product));
      }
    }
  }
  aver=me/4.;
  return output;
}

void MEee2Z::doinit() {
  MEBase::doinit();
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm)
    { _theFFZVertex = hwsm->vertexFFZ();}
  else
    {throw InitException();}
}
