// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGQBCQBase class.
//

#include "MEGQBCQBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

using namespace Herwig;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<MEGQBCQBase,HwMEBase>
describeHerwigMEGQBCQBase("Herwig::MEGQBCQBase", "HwOniumParameters.so HwMEHadronOnium.so");

void MEGQBCQBase::Init() {
  
  static ClassDocumentation<MEGQBCQBase> documentation
    ("The MEGQBCQBase class is the base class fpr g c -> B_c b processes");
  
  static Parameter<MEGQBCQBase,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &MEGQBCQBase::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Reference<MEGQBCQBase,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &MEGQBCQBase::params_, false, false, true, false, false);

}

void MEGQBCQBase::doinit() {
  HwMEBase::doinit();
  long pid = 100000*(n_-1)+id_;
  state_ = getParticleData(pid);
  if(!state_)
    throw Exception() << "No B_c state with pid = " << pid << "in MEGQBCQBase::doinit()" << Exception::runerror;
  setMassGenerator(dynamic_ptr_cast<GenericMassGeneratorPtr>(state_->massGenerator()));
}

void MEGQBCQBase::persistentOutput(PersistentOStream & os) const {
  os << n_ << params_ << state_;
}

void MEGQBCQBase::persistentInput(PersistentIStream & is, int) {
  is >> n_ >> params_ >> state_;
}

void MEGQBCQBase::getDiagrams() const {
  tcPDPtr g = getParticleData(ParticleID::g);
  tcPDPtr c = getParticleData(4);
  tcPDPtr b = getParticleData(5);
  // c initiated
  add(new_ptr((Tree2toNDiagram(3), g, c, c, 2, state_, 1, b,     -1)));
  add(new_ptr((Tree2toNDiagram(2), g, c, 1, c , 3, state_, 3, b, -2)));
  // cbar initiated
  add(new_ptr((Tree2toNDiagram(3), g, c->CC(), c->CC(), 2, state_->CC(), 1, b->CC(),     -3)));
  add(new_ptr((Tree2toNDiagram(2), g, c->CC(), 1, c->CC() , 3, state_->CC(), 3, b->CC(), -4)));
}

Selector<const ColourLines *>
MEGQBCQBase::colourGeometries(tcDiagPtr diag) const {
  static ColourLines c[4] = {ColourLines("1 5, 3 2 -1"),
			     ColourLines("1 3 5, 2 -1"),
			     ColourLines("-1 -5, -3 -2 1"),
			     ColourLines("-1 -3 -5, -2 1")};
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c[abs(diag->id())-1]);
  return sel;
}

Selector<MEBase::DiagramIndex>
MEGQBCQBase::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 || diags[i]->id()==-3)
      sel.insert(meInfo()[0], i);
    else
      sel.insert(meInfo()[1], i);
  return sel;
}

Energy2 MEGQBCQBase::scale() const {
  return sHat();
}

void MEGQBCQBase::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.reserve(4);
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // get them in the right order
  bool swapped(false);
  if(hard[0]->id()!=ParticleID::g) {
    swapped=true;
  }
  if(swapped) {
    swap(hard[0],hard[1]);
    swap(hard[2],hard[3]);
  }
  // boost to partonic CMS
  Lorentz5Momentum pcms = hard[0]->momentum()+hard[1]->momentum();
  LorentzRotation boost(-pcms.boostVector());
  for(PPtr part : hard) part->transform(boost);
  // set basis states and compute the matrix element
  // B_c state
  if(hard[2]->dataPtr()->iSpin()==PDT::Spin0) {
    ScalarWaveFunction(      hard[2],outgoing,true);
  }
  else if(hard[2]->dataPtr()->iSpin()==PDT::Spin1) {
    vector<VectorWaveFunction> v3;
    VectorWaveFunction(v3,hard[2],outgoing,true ,false,true,vector_phase);
  }
  else
    assert(false);
  // gluon
  vector<VectorWaveFunction> g1;
  VectorWaveFunction(   g1,hard[0],incoming,false,true,true,vector_phase);
  g1[1]=g1[2];
  // other particles
  if(hard[1]->id()>0) {
    vector<SpinorWaveFunction> q2;
    vector<SpinorBarWaveFunction> q4;
    SpinorWaveFunction(   q2,hard[1],incoming,false,true);
    SpinorBarWaveFunction(q4,hard[3],outgoing,true ,true);
  }
  else {
    vector<SpinorBarWaveFunction> q2;
    vector<SpinorWaveFunction>    q4;
    SpinorBarWaveFunction(q2,hard[1],incoming,false,true);
    SpinorWaveFunction(   q4,hard[3],outgoing,true ,true);
  }
  // rescale the momenta
  double rr = hard[1]->dataPtr()->mass()/hard[3]->dataPtr()->mass();
  Energy mBc = hard[2]->mass();
  Energy ecm=sqrt(sHat());
  vector<Lorentz5Momentum> rescaled(4);
  // masses
  rescaled[0].setMass(          ZERO);
  rescaled[1].setMass(rr/(1.+rr)*mBc);
  rescaled[2].setMass(           mBc);
  rescaled[3].setMass(1./(1.+rr)*mBc);
  // incoming
  Energy pin = 0.5*(sHat()-sqr(rescaled[1].mass()))/ecm;
  if(hard[0]->momentum().z()>ZERO) {
    rescaled[0].setZ(pin); rescaled[1].setZ(-pin);
  }
  else {
    rescaled[0].setZ(-pin); rescaled[1].setZ(pin);
  }
  rescaled[0].setT(pin); rescaled[1].setT(0.5*(sHat()+sqr(rescaled[1].mass()))/ecm);
  // outgoing
  Energy q;
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), rescaled[2].mass(), rescaled[3].mass());
  } 
  catch ( ImpossibleKinematics & e) {
    assert(false);
  }
  rescaled[2].setVect(hard[2]->momentum().vect().unit()*q);
  rescaled[3].setVect(hard[3]->momentum().vect().unit()*q);
  rescaled[2].rescaleEnergy();
  rescaled[3].rescaleEnergy();
  rescaledMomenta(rescaled);
  // calculate the matrix element
  me2();
  // construct the vertex
  HardVertexPtr hardvertex = new_ptr(HardVertex());
  // // set the matrix element for the vertex
  hardvertex->ME(me_);
  // set the pointers and to and from the vertex
  for(unsigned int i = 0; i < hard.size(); ++i)
    hard[i]->spinInfo()->productionVertex(hardvertex);
  // boost back to lab
  boost = LorentzRotation(pcms.boostVector());
  for(PPtr part : hard)
    part->transform(boost);
}
