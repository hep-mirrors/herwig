// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGGto1D2 class.
//

#include "MEGGto1D2.h"
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
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"

using namespace Herwig;

void MEGGto1D2::doinit() {
  HwMEBase::doinit();
  // get the non-perturbative ME
  O1_ = params_->singletMEProduction<2>(state_,n_,0,2);
  // get the mass generator of the onium state
  unsigned int iq = 4+state_;
  long id = iq*110+10005 + (n_-1)*100000;
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

IBPtr MEGGto1D2::clone() const {
  return new_ptr(*this);
}

IBPtr MEGGto1D2::fullclone() const {
  return new_ptr(*this);
}

void MEGGto1D2::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*GeV2*GeV2*GeV2) << oenum(state_) << n_ << massGen_;
}

void MEGGto1D2::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*GeV2*GeV2*GeV2) >> ienum(state_) >> n_ >> massGen_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEGGto1D2,HwMEBase>
describeHerwigMEGGto1D2("Herwig::MEGGto1D2",
			"HwOniumParameters.so HwMEHadronOnium.so");

void MEGGto1D2::Init() {

  static ClassDocumentation<MEGGto1D2> documentation
    ("The MEGGto1D2 class implements the colour singlet matrix element for g g -> 1D2");

  static Reference<MEGGto1D2,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &MEGGto1D2::params_, false, false, true, false, false);
  
  static Switch<MEGGto1D2,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &MEGGto1D2::state_, ccbar, false, false);
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
  
  static Parameter<MEGGto1D2,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &MEGGto1D2::n_, 1, 1, 10,
     false, false, Interface::limited);
}

void MEGGto1D2::getDiagrams() const {
  // construct the meson PDG code from quark ids
  unsigned int iq = 4+state_;
  tcPDPtr ps = getParticleData(long(iq*110+10005 + (n_-1)*100000));
  tcPDPtr g = getParticleData(ParticleID::g);
  add(new_ptr((Tree2toNDiagram(2), g, g, 1, ps, -1)));
}

Energy2 MEGGto1D2::scale() const {
  return sHat();
}

int MEGGto1D2::nDim() const {
  return 0;
}

Selector<MEBase::DiagramIndex>
MEGGto1D2::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) sel.insert(1.0, i);
  return sel;
}

Selector<const ColourLines *>
MEGGto1D2::colourGeometries(tcDiagPtr) const {
  static ColourLines cFlow("1 -2, 2 -1");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &cFlow);
  return sel;
}

bool MEGGto1D2::generateKinematics(const double * ) {
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

CrossSection MEGGto1D2::dSigHatDR() const {
  return Constants::pi*sqr(hbarc)*me2()*jacobian()*massGen_->BreitWignerWeight(sqrt(sHat()));
}

double MEGGto1D2::me2() const {
  return 64./135.*sqr(Constants::pi)*O1_/pow<7,2>(sHat())*sqr(standardModel()->alphaS(scale()));
}

void MEGGto1D2::constructVertex(tSubProPtr sub) {
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
  // 1D2 wavefunction
  vector<TensorWaveFunction> twave;
  TensorWaveFunction(twave,hard[2],outgoing,true,false);
  // matrix element
  ProductionMatrixElement me(PDT::Spin1,PDT::Spin1,PDT::Spin2);

  Lorentz5Momentum pDiff = hard[0]->momentum()-hard[1]->momentum();
  Energy M = hard[2]->mass();
  vector<Complex> tamp(5);
  for(unsigned int i=0;i<5;++i) {
    tamp[i] = twave[i].wave().trace()+2./sqr(M)*(twave[i].wave().preDot(pDiff)*pDiff);
  }
  // for(unsigned int ih1=0;ih1<2;++ih1) {
  //   auto vOff = epsilon(g1[ih1].wave(),hard[0]->momentum(),hard[1]->momentum());
  //   for(unsigned int ih2=0;ih2<2;++ih2) {
  //     Complex vamp = (vOff*g2[ih2].wave())/sqr(M);
  //     for(unsigned int ih3=0;ih3<5;++ih3) {
  // 	me(2*ih1,2*ih2,ih3) = tamp[ih3]*vamp*Complex(0,1)*sqrt(1.5);
  //     }
  //   }
  // }
  me(0,0,2) = -1.;
  me(2,2,2) =  1.;
  // construct the vertex
  HardVertexPtr hardvertex = new_ptr(HardVertex());
  // // set the matrix element for the vertex
  hardvertex->ME(me);
  // set the pointers and to and from the vertex
  for(unsigned int i = 0; i < 3; ++i)
    hard[i]->spinInfo()->productionVertex(hardvertex);
}
