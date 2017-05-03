// -*- C++ -*-
//
// MEee2VectorMeson.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2VectorMeson class.
//

#include "MEee2VectorMeson.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;


void MEee2VectorMeson::getDiagrams() const {
  tcPDPtr em = getParticleData(ParticleID::eminus);
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  add(new_ptr((Tree2toNDiagram(2), em, ep, 1, _vector,-1)));
}

Energy2 MEee2VectorMeson::scale() const {
  return sHat();
}

int MEee2VectorMeson::nDim() const {
  return 1;
}

void MEee2VectorMeson::setKinematics() {
  MEBase::setKinematics(); // Always call the base class method first.
}

bool MEee2VectorMeson::generateKinematics(const double *r) {
  Lorentz5Momentum pout=meMomenta()[0]+meMomenta()[1];
  pout.rescaleMass();
  meMomenta()[2] = pout;
  jacobian(1.0);
  // check passes all the cuts
  vector<LorentzMomentum> out(1,meMomenta()[2]);
  tcPDVector tout(1,mePartonData()[2]);
  jacobian(1.+0.002*(0.5-r[0]));
  // return true if passes the cuts
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}

unsigned int MEee2VectorMeson::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2VectorMeson::orderInAlphaEW() const {
  return 0;
}

Selector<const ColourLines *>
MEee2VectorMeson::colourGeometries(tcDiagPtr) const {
  static ColourLines neutral ( " " );
  Selector<const ColourLines *> sel;sel.insert(1.,&neutral);
  return sel;
}

void MEee2VectorMeson::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _vector << _massgen << _lineshape;
}

void MEee2VectorMeson::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _vector >> _massgen >> _lineshape;
}

ClassDescription<MEee2VectorMeson> MEee2VectorMeson::initMEee2VectorMeson;
// Definition of the static class description member.

void MEee2VectorMeson::Init() {

  static ClassDocumentation<MEee2VectorMeson> documentation
    ("The MEee2VectorMeson class implements the production of a vector meson"
     " in e+e- collisions and is primilarly intended to test the hadron decay package");

  static Switch<MEee2VectorMeson,bool> interfaceLineShape
    ("LineShape",
     "Option for the vector meson lineshape",
     &MEee2VectorMeson::_lineshape, false, false, false);
  static SwitchOption interfaceLineShapeMassGenerator
    (interfaceLineShape,
     "MassGenerator",
     "Use the mass generator if available",
     true);
  static SwitchOption interfaceLineShapeBreitWigner
    (interfaceLineShape,
     "BreitWigner",
     "Use a Breit-Wigner with the naive running width",
     false);

  static Reference<MEee2VectorMeson,ParticleData> interfaceVectorMeson
    ("VectorMeson",
     "The vector meson produced",
     &MEee2VectorMeson::_vector, false, false, true, false, false);

  static Parameter<MEee2VectorMeson,double> interfaceCoupling
    ("Coupling",
     "The leptonic coupling of the vector meson",
     &MEee2VectorMeson::_coupling, 0.0012, 0.0, 10.0,
     false, false, Interface::limited);

}

Selector<MEBase::DiagramIndex>
MEee2VectorMeson::diagrams(const DiagramVector &) const {
  Selector<DiagramIndex> sel;sel.insert(1.0, 0);
  return sel;
}

CrossSection MEee2VectorMeson::dSigHatDR() const {
  InvEnergy2 wgt;
  Energy  M(_vector->mass()),G(_vector->width());
  Energy2 M2(sqr(M)),GM(G*M);
  if(_massgen&&_lineshape) {
    wgt =Constants::pi*_massgen->weight(sqrt(sHat()))/(sqr(sHat()-M2)+sqr(GM))*UnitRemoval::E2;
  }
  else {
    wgt = sHat()*G/M/(sqr(sHat()-M2)+sqr(sHat()*G/M));
  }
  return me2()*jacobian()*wgt*sqr(hbarc);
}

void MEee2VectorMeson::doinit() {
  MEBase::doinit();
  // mass generator
  tMassGenPtr mass=_vector->massGenerator();
  if(mass) {
    _massgen=dynamic_ptr_cast<GenericMassGeneratorPtr>(mass);
  }
}

double MEee2VectorMeson::me2() const {
  double aver=0.;
  // get the order right
  int ielectron(0),ipositron(1);
  if(mePartonData()[0]->id()!=11) swap(ielectron,ipositron);
  // the vectors for the wavefunction to be passed to the matrix element
  vector<SpinorWaveFunction> fin;
  vector<SpinorBarWaveFunction> ain;
  vector<VectorWaveFunction> vout;
  for(unsigned int ihel=0;ihel<2;++ihel) {
    fin.push_back(SpinorWaveFunction(meMomenta()[ielectron],
				     mePartonData()[ielectron],ihel,incoming));
    ain.push_back(SpinorBarWaveFunction(meMomenta()[ipositron],
					mePartonData()[ipositron],ihel,incoming));
  }
  for(unsigned int ihel=0;ihel<3;++ihel) {
    vout.push_back(VectorWaveFunction(meMomenta()[2],mePartonData()[2],ihel,outgoing));
  }
  ProductionMatrixElement temp=HelicityME(fin,ain,vout,aver);
  return aver;
}

// the helicity amplitude matrix element
ProductionMatrixElement MEee2VectorMeson::HelicityME(vector<SpinorWaveFunction> fin,
						     vector<SpinorBarWaveFunction> ain,
						     vector<VectorWaveFunction> vout,
						     double & aver) const {
  ProductionMatrixElement output(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1);
  Complex product;
  // sum over helicities to get the matrix element
  unsigned int inhel1,inhel2,outhel1;
  double me(0.);
  LorentzPolarizationVector vec;
  Complex ii(0.,1.);
  for(inhel1=0;inhel1<2;++inhel1) {	  
    for(inhel2=0;inhel2<2;++inhel2) {
      vec =  fin[inhel1].wave().vectorCurrent(ain[inhel2].wave());
      vec*=_coupling;
      for(outhel1=0;outhel1<3;++outhel1) {
	product = vec.dot(vout[outhel1].wave());
	output(inhel1,inhel2,outhel1)=product;
	me+=real(product*conj(product));
      }
    }
  }
  aver=(me*UnitRemoval::E2)/sHat();
  return output;
}

void MEee2VectorMeson::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  if(hard[0]->id()<hard[1]->id()) swap(hard[0],hard[1]);
  vector<SpinorWaveFunction>    fin;
  vector<SpinorBarWaveFunction> ain;
  vector<VectorWaveFunction> vout;
  SpinorWaveFunction(   fin ,hard[0],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
  VectorWaveFunction(vout   ,hard[2],outgoing,true,false,true);
  double dummy;
  ProductionMatrixElement prodme=HelicityME(fin,ain,vout,dummy);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(prodme);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<3;++ix) {
    (hard[ix]->spinInfo())->productionVertex(hardvertex);
  }
}





