// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGammaGamma2X class.
//

#include "MEGammaGamma2X.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

void MEGammaGamma2X::getDiagrams() const {
  vector<DiagPtr> diags = amp_->getDiagrams(0);
  for(DiagPtr diag : diags) add(diag);
}

void MEGammaGamma2X::setKinematics() {
  HwMEBase::setKinematics(); // Always call the base class method first.
}

IBPtr MEGammaGamma2X::clone() const {
  return new_ptr(*this);
}

IBPtr MEGammaGamma2X::fullclone() const {
  return new_ptr(*this);
}

void MEGammaGamma2X::persistentOutput(PersistentOStream & os) const {
  os << amp_;
}

void MEGammaGamma2X::persistentInput(PersistentIStream & is, int) {
  is >> amp_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEGammaGamma2X,HwMEBase>
describeHerwigMEGammaGamma2X("Herwig::MEGammaGamma2X", "HwMEGammaGamma.so");

void MEGammaGamma2X::Init() {

  static ClassDocumentation<MEGammaGamma2X> documentation
    ("The MEGammaGamma2X class implements the matrix element for gamma gamma ->X");
  
  static Reference<MEGammaGamma2X,GammaGammaAmplitude> interfaceAmplitude
    ("Amplitude",
     "The gamma gamma -> X amplitude",
     &MEGammaGamma2X::amp_, false, false, true, false, false);

}

Selector<const ColourLines *>
MEGammaGamma2X::colourGeometries(tcDiagPtr diag) const {
  return amp_->colourGeometries(0,mePartonData(),diag);
}

bool MEGammaGamma2X::generateKinematics(const double * r) {
  vector<Lorentz5Momentum> pout(meMomenta().size()-2);
  tcPDVector tout(mePartonData().begin()+2,mePartonData().end());
  double jac = amp_->generateKinematics(r,sHat(),pout,tout);
  if(jac<0.) return false;
  jacobian(jac);
  for(unsigned int ix=0;ix<pout.size();++ix)
    meMomenta()[ix+2] = pout[ix];
  vector<LorentzMomentum> p2out(meMomenta().size()-2);
  for(unsigned int ix=0;ix<pout.size();++ix)
    p2out[ix] = pout[ix];
  return lastCuts().passCuts(tout, p2out, mePartonData()[0], mePartonData()[1]);
}

Energy2 MEGammaGamma2X::scale() const {
  return sHat();
}

CrossSection MEGammaGamma2X::dSigHatDR() const {
  return 0.5*sqr(hbarc)*me2()*jacobian()*pow(Constants::twopi,4)/sHat();
}

double MEGammaGamma2X::me2() const {
  using ThePEG::Helicity::incoming;
  VectorWaveFunction      p1w(meMomenta()[0],mePartonData()[0],incoming);
  VectorWaveFunction      p2w(meMomenta()[1],mePartonData()[1],incoming);
  vector<VectorWaveFunction> p1,p2;
  for(unsigned int ix=0;ix<2;++ix) {
    p1w.reset(2*ix);p1.push_back(p1w);
    p2w.reset(2*ix);p2.push_back(p2w);
    //p1w.reset(10);p1.push_back(p1w);
    //p2w.reset(10);p2.push_back(p2w);
  }
  DVector save;
  double output =  amp_->me2(p1,p2,meMomenta()[0].m2(),meMomenta()[1].m2(),sHat(),
			     vector<Lorentz5Momentum>(meMomenta().begin()+2,meMomenta().end()),
			     cPDVector(mePartonData().begin()+2,mePartonData().end()),save);
  meInfo(save);
  return output;
}

Selector<MEBase::DiagramIndex>
MEGammaGamma2X::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    unsigned int id = abs(diags[i]->id())-1;
    if(id<meInfo().size())
      sel.insert(meInfo()[id], i);
    else if ( meInfo().empty() ) sel.insert(1., i);
    else
      assert(false);
  }
  return sel;
}

void MEGammaGamma2X::constructVertex(tSubProPtr sub) {
  using namespace ThePEG::Helicity;
  // extract the particles in the hard process
  tParticleVector hard;
  hard.reserve(sub->outgoing().size()+2);
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  for(unsigned int ix=0;ix<sub->outgoing().size();++ix)
    hard.push_back(sub->outgoing()[ix]);
  // calculate the photon currents (removing null 0 state)
  vector<VectorWaveFunction> p1,p2;
  VectorWaveFunction   (p1,hard[0],incoming,false,true );
  p1[1]=p1[2];p1.resize(2);
  VectorWaveFunction   (p2,hard[1],incoming,false,true );
  p2[1]=p2[2];p2.resize(2);
  tParticleVector pTemp(hard.begin()+2,hard.end());
  ProductionMatrixElement me = amp_->me(p1,p2,pTemp);
  // // for(auto part : hard ) cerr << *part << "\n";
  // // for(unsigned int ih1=0;ih1<2;++ih1)
  // //   for(unsigned int ih2=0;ih2<2;++ih2)
  // // 	  for(unsigned int ih5=0;ih5<hard.back()->dataPtr()->iSpin();++ih5)
  // // 	    cerr << ih1 << " " << ih2 << " " << ih5 << " " << me(ih1,ih2,ih5) << "\n";
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<hard.size();++ix) {
    hard[ix]->spinInfo()->productionVertex(hardvertex);
  }
}



















