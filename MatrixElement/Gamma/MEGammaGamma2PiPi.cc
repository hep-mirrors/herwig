// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGammaGamma2PiPi class.
//

#include "MEGammaGamma2PiPi.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;
using ThePEG::Helicity::incoming;

MEGammaGamma2PiPi::MEGammaGamma2PiPi() {
  massOption(vector<unsigned int>(2,1));
}

void MEGammaGamma2PiPi::getDiagrams() const {
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr pip = getParticleData(ParticleID::piplus);
  tcPDPtr pim = pip->CC();
  // first t-channel
  add(new_ptr((Tree2toNDiagram(3),gamma,pip,gamma,1,pip, 2,pim,-1)));
  // interchange
  add(new_ptr((Tree2toNDiagram(3),gamma,pip,gamma,2,pip, 1,pim,-2)));
}

Energy2 MEGammaGamma2PiPi::scale() const {
  return 2.*sHat()*tHat()*uHat()/(sqr(sHat())+sqr(tHat())+sqr(uHat()));
}

double MEGammaGamma2PiPi::me2() const {
  VectorWaveFunction      p1w(rescaledMomenta()[0],mePartonData()[0],incoming);
  VectorWaveFunction      p2w(rescaledMomenta()[1],mePartonData()[1],incoming);
  vector<VectorWaveFunction> p1,p2;
  for(unsigned int ix=0;ix<2;++ix) {
    p1w.reset(2*ix);p1.push_back(p1w);
    p2w.reset(2*ix);p2.push_back(p2w);
  }
  // calculate the matrix element
  return helicityME(p1,p2,rescaledMomenta(),false);
}

unsigned int MEGammaGamma2PiPi::orderInAlphaS() const {
  return 0;
}

unsigned int MEGammaGamma2PiPi::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEGammaGamma2PiPi::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 )      sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
  return sel;
}

Selector<const ColourLines *>
MEGammaGamma2PiPi::colourGeometries(tcDiagPtr ) const {
  static ColourLines c1("");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c1);
  return sel;
}

IBPtr MEGammaGamma2PiPi::clone() const {
  return new_ptr(*this);
}

IBPtr MEGammaGamma2PiPi::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<MEGammaGamma2PiPi,HwMEBase>
describeHerwigMEGammaGamma2PiPi("Herwig::MEGammaGamma2PiPi", "HwMEGammaGamma.so");

void MEGammaGamma2PiPi::Init() {

  static ClassDocumentation<MEGammaGamma2PiPi> documentation
    ("The MEGammaGamma2PiPi class provides a simple matrix element for gamma gamma -> pi+pi-");

}

void MEGammaGamma2PiPi::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // order of particles
  unsigned int order[4]={0,1,2,3};
  // identify the process and calculate matrix element
  vector<VectorWaveFunction> p1,p2;
  if(hard[order[2]]->id()<0) swap(order[2],order[3]);
  VectorWaveFunction   (p1  ,hard[order[0]],incoming,false,true);
  VectorWaveFunction   (p2  ,hard[order[1]],incoming,false,true);
  p1[1]=p1[2];
  p2[1]=p2[2];
  vector<Lorentz5Momentum> momenta;
  for(unsigned int ix=0;ix<4;++ix)
    momenta.push_back(hard[order[ix]]->momentum());
  helicityME(p1,p2,momenta,true);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me_);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) {
    if(hard[order[ix]]->spinInfo())
      hard[order[ix]]->spinInfo()->productionVertex(hardvertex);
  }
}

double MEGammaGamma2PiPi::helicityME(vector<VectorWaveFunction> &p1,
				     vector<VectorWaveFunction> &p2,
				     const vector<Lorentz5Momentum> & momenta,
				     bool calc) const {
  // for(unsigned int ix=0;ix<4;++ix)
  //   cerr << mePartonData()[ix]->PDGName() <<" " << meMomenta()[ix]/GeV << " " << meMomenta()[ix].mass()/GeV << "\n";
  // double beta = meMomenta()[2].vect().mag()/meMomenta()[2].t();
  // double ct   = meMomenta()[2].cosTheta();
  // double phi  = meMomenta()[2].phi();
  // matrix element to be stored
  if(calc) me_.reset(ProductionMatrixElement(PDT::Spin1,PDT::Spin1,
  					     PDT::Spin0,PDT::Spin0));
  // calculate the matrix element
  double output(0.),sumdiag[2]={0.,0.};
  Complex diag[2];
  Energy2 d1 = momenta[0]*momenta[2], d2 = momenta[1]*momenta[2];
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    complex<Energy> a1 = momenta[2]*p1[ihel1].wave();
    complex<Energy> a2 = momenta[3]*p1[ihel1].wave();
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      complex<Energy> b1 = momenta[2]*p2[ihel2].wave();
      complex<Energy> b2 = momenta[3]*p2[ihel2].wave();
      // t-channel diagram
      diag[0] = a1*b2/d1;
      // u-channel diagram
      diag[1] = a2*b1/d2;
      sumdiag[0] += norm(diag[0]);
      sumdiag[1] += norm(diag[1]);
      Complex amp = p1[ihel1].wave()*p2[ihel2].wave()-diag[0]-diag[1];
      output += norm(amp);
      // store the me if needed
      if(calc) me_(2*ihel1,2*ihel2,0,0) = amp;
    }
  }
  // diagrams
  DVector save;
  save.push_back(sumdiag[0]);
  save.push_back(sumdiag[1]);
  meInfo(save);
  // code to test vs the analytic result
  // Energy2 m2 = sqr(momenta[2].mass());
  // double test = 2. + ((d1 + d2)*m2*(-2*d1*d2 + (d1 + d2)*m2))/(sqr(d1)*sqr(d2));
  // cerr << "testing ME " << (output-test)/(output+test) << " " << output << " " << test << " " << output/test << "\n"; 
  // coupling factors
  return output*sqr(SM().alphaEM()*4.*Constants::pi);
}
