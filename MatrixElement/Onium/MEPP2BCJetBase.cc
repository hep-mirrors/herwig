// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2BCJetBase class.
//

#include "MEPP2BCJetBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
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
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/Rank3TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

using namespace Herwig;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<MEPP2BCJetBase,HwMEBase>
describeHerwigMEPP2BCJetBase("Herwig::MEPP2BCJetBase", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP2BCJetBase::Init() {
  
  static ClassDocumentation<MEPP2BCJetBase> documentation
    ("The MEPP2BCJetBase class is the base class for g q -> B_c q and q qbar -> Bc g processes");
  
  static Parameter<MEPP2BCJetBase,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &MEPP2BCJetBase::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Reference<MEPP2BCJetBase,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &MEPP2BCJetBase::params_, false, false, true, false, false);

  static Switch<MEPP2BCJetBase,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2BCJetBase::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcessGQ
    (interfaceProcess,
     "GQ",
     "Only include the gq/gqbar initiated processes",
     1);
  static SwitchOption interfaceProcessQQbar
    (interfaceProcess,
     "QQbar",
     "Only include the q qbar initiated processes",
     2);

}

void MEPP2BCJetBase::doinit() {
  HwMEBase::doinit();
  long pid = 100000*(n_-1)+id_;
  state_ = getParticleData(pid);
  if(!state_)
    throw Exception() << "No B_c state with pid = " << pid << "in MEPP2BCJetBase::doinit()" << Exception::runerror;
  setMassGenerator(dynamic_ptr_cast<GenericMassGeneratorPtr>(state_->massGenerator()));
}

void MEPP2BCJetBase::persistentOutput(PersistentOStream & os) const {
  os << n_ << process_ << params_ << state_;
}

void MEPP2BCJetBase::persistentInput(PersistentIStream & is, int) {
  is >> n_ >> process_ >> params_ >> state_;
}

void MEPP2BCJetBase::getDiagrams() const {
  tcPDPtr g = getParticleData(ParticleID::g);
  tcPDPtr c = getParticleData(4);
  tcPDPtr b = getParticleData(5);
  if(process_!=2) {
    // c initiated
    add(new_ptr((Tree2toNDiagram(3), g, c, c, 2, state_, 1, b,     -1)));
    add(new_ptr((Tree2toNDiagram(2), g, c, 1, c , 3, state_, 3, b, -2)));
    // b initiated
    add(new_ptr((Tree2toNDiagram(3), g, b, b, 2, state_->CC(), 1, c,     -1)));
    add(new_ptr((Tree2toNDiagram(2), g, b, 1, b , 3, state_->CC(), 3, c, -2)));
    // cbar initiated
    add(new_ptr((Tree2toNDiagram(3), g, c->CC(), c->CC(), 2, state_->CC(), 1, b->CC(),     -3)));
    add(new_ptr((Tree2toNDiagram(2), g, c->CC(), 1, c->CC() , 3, state_->CC(), 3, b->CC(), -4)));
    // bbar initiated
    add(new_ptr((Tree2toNDiagram(3), g, b->CC(), b->CC(), 2, state_, 1, c->CC(),     -3)));
    add(new_ptr((Tree2toNDiagram(2), g, b->CC(), 1, b->CC() , 3, state_, 3, c->CC(), -4)));
  }
  if(process_!=1) {
    // b cbar -> Bc g
    add(new_ptr((Tree2toNDiagram(3), b, c, c->CC(), 1, state_->CC(), 2, g, -5)));
    add(new_ptr((Tree2toNDiagram(3), b, b->CC(), c->CC(), 2, state_->CC(), 1, g, -6)));
    // c bbar -> Bc g
    add(new_ptr((Tree2toNDiagram(3), c, b, b->CC(), 1, state_, 2, g, -5)));
    add(new_ptr((Tree2toNDiagram(3), c, c->CC(), b->CC(), 2, state_, 1, g, -6)));
  }
}

Selector<const ColourLines *>
MEPP2BCJetBase::colourGeometries(tcDiagPtr diag) const {
  static ColourLines c[6] = {ColourLines("1 5, 3 2 -1"),
			     ColourLines("1 3 5, 2 -1"),
			     ColourLines("-1 -5, -3 -2 1"),
			     ColourLines("-1 -3 -5, -2 1"),
			     ColourLines("1 2 5, -3 -5"),
			     ColourLines("1 5, -3 -2 -5")};
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c[abs(diag->id())-1]);
  return sel;
}

Selector<MEBase::DiagramIndex>
MEPP2BCJetBase::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 || diags[i]->id()==-3)
      sel.insert(meInfo()[0], i);
    else
      sel.insert(meInfo()[1], i);
  return sel;
}

Energy2 MEPP2BCJetBase::scale() const {
  return sHat();
}

void MEPP2BCJetBase::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.reserve(4);
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // get them in the right order
  bool swapped(false);
  if(hard[1]->id()==ParticleID::g || hard[0]->id()<0) {
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
  else if(hard[2]->dataPtr()->iSpin()==PDT::Spin2) {
    vector<TensorWaveFunction> t3;
    TensorWaveFunction(t3,hard[2],outgoing,true ,false,true,tensor_phase);
  }
  else if(hard[2]->dataPtr()->iSpin()==PDT::Spin3) {
    vector<Rank3TensorWaveFunction> t3;
    Rank3TensorWaveFunction(t3,hard[2],outgoing,true ,false,true);
  }
  else
    assert(false);
  // masses and momenta
  Energy mBc = hard[2]->mass();
  Energy ecm=sqrt(sHat());
  vector<Lorentz5Momentum> rescaled(4);
  // gluon initiate
  if(hard[0]->id()==ParticleID::g) {
    vector<VectorWaveFunction> g1;
    VectorWaveFunction(   g1,hard[0],incoming,false,true,true,vector_phase);
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
    // masses for rescaling
    double rr = hard[1]->dataPtr()->mass()/hard[3]->dataPtr()->mass();
    // masses
    rescaled[0].setMass(          ZERO);
    rescaled[1].setMass(rr/(1.+rr)*mBc);
    rescaled[2].setMass(           mBc);
    rescaled[3].setMass(1./(1.+rr)*mBc);
  }
  // q qbar
  else {
    vector<SpinorWaveFunction> q1;
    vector<SpinorBarWaveFunction> q2;
    SpinorWaveFunction(   q1,hard[0],incoming,false,true);
    SpinorBarWaveFunction(q2,hard[1],incoming,false,true);
    vector<VectorWaveFunction> g4;
    VectorWaveFunction(   g4,hard[3],outgoing,false,true,true,vector_phase);
    // masses for rescaling
    double rr = mePartonData()[1]->mass()/mePartonData()[0]->mass();
    // masses
    rescaled[0].setMass(1./(1.+rr)*mBc);
    rescaled[1].setMass(rr/(1.+rr)*mBc);
    rescaled[2].setMass(           mBc);
    rescaled[3].setMass(     ZERO     );
  }
  // rescale the momenta
  // incoming
  Energy pin = SimplePhaseSpace::getMagnitude(sHat(), rescaled[0].mass(), rescaled[1].mass());
  rescaled[0].setZ(pin); rescaled[1].setZ(-pin);
  rescaled[0].setT(0.5*(sHat()+sqr(rescaled[0].mass())-sqr(rescaled[1].mass()))/ecm);
  rescaled[1].setT(0.5*(sHat()-sqr(rescaled[0].mass())+sqr(rescaled[1].mass()))/ecm);
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
