// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGGto3P0 class.
//

#include "MEGGto3P0.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

void MEGGto3P0::doinit() {
  HwMEBase::doinit();
  // get the non-perturbative ME
  O1_ = params_->singletMEProduction<1>(state_,n_,1,0);
  // get the mass generator of the onium state
  unsigned int iq = 4+state_;
  long id = iq*110+10001 + (n_-1)*100000;
  tcPDPtr ps = getParticleData(id);
  if(!ps)
    throw Exception() << "No onium particle with id " << id
		      << "in " << fullName();
  if(ps->massGenerator())
    massGen_=dynamic_ptr_cast<GenericMassGeneratorPtr>(ps->massGenerator());
  if(!massGen_)
    throw Exception() << "Must have mass generator for " << ps->PDGName()
		      << "in " << fullName();
}

IBPtr MEGGto3P0::clone() const {
  return new_ptr(*this);
}

IBPtr MEGGto3P0::fullclone() const {
  return new_ptr(*this);
}

void MEGGto3P0::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*GeV2*GeV2) << oenum(state_) << n_ << massGen_;
}

void MEGGto3P0::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*GeV2*GeV2) >> ienum(state_) >> n_ >> massGen_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEGGto3P0,HwMEBase>
describeHerwigMEGGto3P0("Herwig::MEGGto3P0",
			"HwOniumParameters.so HwMEHadronOnium.so");

void MEGGto3P0::Init() {

  static ClassDocumentation<MEGGto3P0> documentation
    ("The MEGGto3P0 class implements the colour singlet matrix element for g g -> 3P0");

  static Reference<MEGGto3P0,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &MEGGto3P0::params_, false, false, true, false, false);
  
  static Switch<MEGGto3P0,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &MEGGto3P0::state_, ccbar, false, false);
  static SwitchOption interfaceStateccbar
    (interfaceState,
     "ccbar",
     "Charmonium state",
     ccbar);
  static SwitchOption interfaceStatebbbar
    (interfaceState,
     "bbbar",
     "Bottomonium state",
     bbbar);
  
  static Parameter<MEGGto3P0,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &MEGGto3P0::n_, 1, 1, 10,
     false, false, Interface::limited);
}

void MEGGto3P0::getDiagrams() const {
  // construct the meson PDG code from quark ids
  unsigned int iq = 4+state_;
  tcPDPtr ps = getParticleData(long(iq*110+10001 + (n_-1)*100000));
  tcPDPtr g = getParticleData(ParticleID::g);
  add(new_ptr((Tree2toNDiagram(2), g, g, 1, ps, -1)));
}

Energy2 MEGGto3P0::scale() const {
  return sHat();
}

int MEGGto3P0::nDim() const {
  return 0;
}

Selector<MEBase::DiagramIndex>
MEGGto3P0::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) sel.insert(1.0, i);
  return sel;
}

Selector<const ColourLines *>
MEGGto3P0::colourGeometries(tcDiagPtr) const {
  static ColourLines cFlow("1 -2, 2 -1");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &cFlow);
  return sel;
}

bool MEGGto3P0::generateKinematics(const double * ) {
  Lorentz5Momentum pout = meMomenta()[0] + meMomenta()[1];
  pout.rescaleMass();
  meMomenta()[2].setMass(pout.mass());
  meMomenta()[2] = LorentzMomentum(pout.x(),pout.y(),pout.z(),pout.t());
  jacobian(1.0);
  // check whether it passes all the cuts: returns true if it does
  vector<LorentzMomentum> out(1,meMomenta()[2]);
  tcPDVector tout(1,mePartonData()[2]);
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}

CrossSection MEGGto3P0::dSigHatDR() const {
  return Constants::pi*sqr(hbarc)*me2()*jacobian()*massGen_->BreitWignerWeight(sqrt(sHat()));
}

double MEGGto3P0::me2() const {
  return 8.*sqr(Constants::pi)*O1_/3./sqr(sHat())/sqrt(sHat())*sqr(standardModel()->alphaS(scale()));
}

void MEGGto3P0::constructVertex(tSubProPtr sub) {
  using namespace ThePEG::Helicity;
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  // gluon wave functions (remove zero longitudinal polarization
  vector<VectorWaveFunction> g1,g2;
  VectorWaveFunction (g1,hard[0],incoming,false,true,true);
  g1[1] = g1[2];
  VectorWaveFunction (g2,hard[1],incoming,false,true,true);
  g2[1] = g2[2];
  // 1S0 wavefunction
  ScalarWaveFunction out(hard[2],outgoing,true);
  // matrix element
  ProductionMatrixElement me(PDT::Spin1,PDT::Spin1,PDT::Spin0);
  // for(unsigned int ih1=0;ih1<2;++ih1) {
  //   complex<Energy> d1 = g1[ih1].wave()*hard[1]->momentum();
  //   for(unsigned int ih2=0;ih2<2;++ih2) {
  //     me(2*ih1,2*ih2,0) = g1[ih1].wave()*g2[ih2].wave()  -2./sqr(hard[2]->mass())*d1*(g2[ih2].wave()*hard[0]->momentum());
  //   }
  // }
  me(0,0,0) =  1.;
  me(2,2,0) = -1.;
  // construct the vertex
  HardVertexPtr hardvertex = new_ptr(HardVertex());
  // // set the matrix element for the vertex
  hardvertex->ME(me);
  // set the pointers and to and from the vertex
  for(unsigned int i = 0; i < 3; ++i)
    hard[i]->spinInfo()->productionVertex(hardvertex);
}
