// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGammaGamma2PseudoScalar class.
//

#include "MEGammaGamma2PseudoScalar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

void MEGammaGamma2PseudoScalar::getDiagrams() const {
  tcPDPtr g=getParticleData(ParticleID::gamma);
  for(long pid=331; pid<340; pid+=110) {
    tcPDPtr ps = getParticleData(pid);
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, ps, -1)));
    
  }
}

Selector<MEBase::DiagramIndex>
MEGammaGamma2PseudoScalar::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for (DiagramIndex i = 0; i < diags.size(); ++i)
    sel.insert(1.0, i);
  return sel;
}

void MEGammaGamma2PseudoScalar::setKinematics() {
  HwMEBase::setKinematics(); // Always call the base class method first.
}

bool MEGammaGamma2PseudoScalar::generateKinematics(const double * ) {
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

unsigned int MEGammaGamma2PseudoScalar::orderInAlphaS() const {
  return 0;
}

unsigned int MEGammaGamma2PseudoScalar::orderInAlphaEW() const {
  return 2;
}

Selector<const ColourLines *>
MEGammaGamma2PseudoScalar::colourGeometries(tcDiagPtr) const {
  static const ColourLines cl("");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &cl);
  return sel;
}

IBPtr MEGammaGamma2PseudoScalar::clone() const {
  return new_ptr(*this);
}

IBPtr MEGammaGamma2PseudoScalar::fullclone() const {
  return new_ptr(*this);
}

void MEGammaGamma2PseudoScalar::persistentOutput(PersistentOStream & os) const {
  os << ounit(F00_,1./GeV) << ounit(LambdaP2_,GeV2);
}

void MEGammaGamma2PseudoScalar::persistentInput(PersistentIStream & is, int) {
  is >> iunit(F00_,1./GeV) >> iunit(LambdaP2_,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEGammaGamma2PseudoScalar,HwMEBase>
describeHerwigMEGammaGamma2PseudoScalar("Herwig::MEGammaGamma2PseudoScalar",
					"HwMEGammaGamma.so");

void MEGammaGamma2PseudoScalar::Init() {

  static ClassDocumentation<MEGammaGamma2PseudoScalar> documentation
    ("The MEGammaGamma2PseudoScalar class implements the production of"
     "pi0, eta eta' mesons for gamma gamma collisions");

  static ParVector<MEGammaGamma2PseudoScalar,InvEnergy> interfaceF00
    ("F00",
     "The form factor at zero momentum transfer",
     &MEGammaGamma2PseudoScalar::F00_, 1./GeV, 3, 1./GeV, 0./GeV, 100./GeV,
     false, false, Interface::limited);

  static ParVector<MEGammaGamma2PseudoScalar,Energy2> interfaceLambdaP2
    ("LambdaP2",
     "The square of the pole mass for the form factor",
     &MEGammaGamma2PseudoScalar::LambdaP2_, GeV2, 3, 1.0*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);
}

CrossSection MEGammaGamma2PseudoScalar::dSigHatDR() const {
  Energy width = mePartonData()[2]->width();
  Energy mass  = mePartonData()[2]->mass();
  Energy mhat = sqrt(sHat());
  if(mhat>mePartonData()[2]->massMax()||mhat<mePartonData()[2]->massMin())
    return ZERO;
  using Constants::pi;
  InvEnergy2 bwfact = width*sqrt(sHat())/(sqr(sHat()-sqr(mass))+sqr(mass*width));
  return me2() * jacobian() * bwfact * sqr(hbarc);
}

double MEGammaGamma2PseudoScalar::helicityME(vector<VectorWaveFunction> &p1,
					     vector<VectorWaveFunction> &p2,
					     const vector<Lorentz5Momentum> & momenta,
					     int iloc,
					     bool calc) const {
  // matrix element to be stored
  if(calc) me_.reset(ProductionMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0));
  // form factor
  // calculate the matrix element
  double output(0.);
  Energy rs = sqrt(sHat());
  InvEnergy form = F00_[iloc]/(1.-momenta[0].m2()/LambdaP2_[iloc])/(1.-momenta[1].m2()/LambdaP2_[iloc]);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    auto v1 = Helicity::epsilon(p1[ihel1].wave(),momenta[0],momenta[1]); 
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      Complex amp = form*(v1*p2[ihel2].wave())/rs;
      if(calc) me_(2*ihel1,2*ihel2,0) = amp;
      output += norm(amp);
    }
  }
  // coupling factors
  return 0.25*output*sqr(SM().alphaEM()*4.*Constants::pi);
}

double MEGammaGamma2PseudoScalar::me2() const {
  using ThePEG::Helicity::incoming;
  int iloc = mePartonData()[2]->id()/100 -1;
  VectorWaveFunction      p1w(meMomenta()[0],mePartonData()[0],incoming);
  VectorWaveFunction      p2w(meMomenta()[1],mePartonData()[1],incoming);
  vector<VectorWaveFunction> p1,p2;
  for(unsigned int ix=0;ix<2;++ix) {
    p1w.reset(2*ix);p1.push_back(p1w);
    p2w.reset(2*ix);p2.push_back(p2w);
  }
  // calculate the matrix element
  return helicityME(p1,p2,meMomenta(),iloc,false);
}

void MEGammaGamma2PseudoScalar::constructVertex(tSubProPtr sub) {
  using ThePEG::Helicity::incoming;
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  vector<VectorWaveFunction> p1,p2;
  VectorWaveFunction   (p1  ,hard[0],incoming,false,true);
  VectorWaveFunction   (p2  ,hard[1],incoming,false,true);
  p1[1]=p1[2];
  p2[1]=p2[2];
  vector<Lorentz5Momentum> momenta;
  for(unsigned int ix=0;ix<3;++ix)
    momenta.push_back(hard[ix]->momentum());
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me_);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<3;++ix)
    if(hard[ix]->spinInfo()) hard[ix]->spinInfo()->productionVertex(hardvertex);
}
