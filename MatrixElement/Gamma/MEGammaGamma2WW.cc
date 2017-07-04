// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGammaGamma2WW class.
//

#include "MEGammaGamma2WW.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

MEGammaGamma2WW::MEGammaGamma2WW() : massOption_(2) 
{}

void MEGammaGamma2WW::doinit() {
  HwMEBase::doinit();
  massOption(vector<unsigned int>(2,massOption_));
  rescalingOption(2);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!hwsm)
    throw InitException() << "Must be the Herwig StandardModel class in "
			  << "MEGammaGamma2WW::doinit" << Exception::abortnow;
  WWWVertex_  = hwsm->vertexWWW();
  WWWWVertex_ = hwsm->vertexWWWW();
}

void MEGammaGamma2WW::getDiagrams() const {
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Wp    = getParticleData(ParticleID::Wplus);
  tcPDPtr Wm    = Wp->CC();
  // first t-channel
  add(new_ptr((Tree2toNDiagram(3),gamma,Wp,gamma,1,Wm, 2,Wp,-1)));
  // interchange
  add(new_ptr((Tree2toNDiagram(3),gamma,Wp,gamma,2,Wm, 1,Wp,-2)));
}

Energy2 MEGammaGamma2WW::scale() const {
  return 2.*sHat()*tHat()*uHat()/(sqr(sHat())+sqr(tHat())+sqr(uHat()));
}

unsigned int MEGammaGamma2WW::orderInAlphaS() const {
  return 0;
}

unsigned int MEGammaGamma2WW::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEGammaGamma2WW::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 )      sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
  return sel;
}

Selector<const ColourLines *>
MEGammaGamma2WW::colourGeometries(tcDiagPtr) const {
  static ColourLines c1(""); 
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c1);
  return sel;
}

IBPtr MEGammaGamma2WW::clone() const {
  return new_ptr(*this);
}

IBPtr MEGammaGamma2WW::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEGammaGamma2WW,HwMEBase>
describeHerwigMEGammaGamma2WW("Herwig::MEGammaGamma2WW", "HwMEGammaGamma.so");

void MEGammaGamma2WW::Init() {

  static ClassDocumentation<MEGammaGamma2WW> documentation
    ("The MEGammaGamma2WW class implements the matrix elements"
     " for e+e- -> W+W-");

  static Switch<MEGammaGamma2WW,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the treatment of the W mass",
     &MEGammaGamma2WW::massOption_, 1, false, false);
  static SwitchOption interfaceMassOptionOnMassShell
    (interfaceMassOption,
     "OnMassShell",
     "The W is produced on its mass shell",
     1);
  static SwitchOption interfaceMassOption2
    (interfaceMassOption,
     "OffShell",
     "The W is generated off-shell using the mass and width generator.",
     2);

}

double MEGammaGamma2WW::me2() const {
  VectorWaveFunction    p1w(rescaledMomenta()[0],mePartonData()[0],incoming);
  VectorWaveFunction    p2w(rescaledMomenta()[1],mePartonData()[1],incoming);
  VectorWaveFunction    w1w(rescaledMomenta()[2],mePartonData()[2],outgoing);
  VectorWaveFunction    w2w(rescaledMomenta()[3],mePartonData()[3],outgoing);
  vector<VectorWaveFunction> p1,p2,w1,w2;
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix!=1) {
      p1w.reset(ix);
      p1.push_back(p1w);
      p2w.reset(ix);
      p2.push_back(p2w);
    }
    w1w.reset(ix);
    w1.push_back(w1w);
    w2w.reset(ix);
    w2.push_back(w2w);
  }
  // calculate the matrix element
  return helicityME(p1,p2,w1,w2,false);
}

void MEGammaGamma2WW::persistentOutput(PersistentOStream & os) const {
  os << massOption_ << WWWVertex_ << WWWWVertex_;
}

void MEGammaGamma2WW::persistentInput(PersistentIStream & is, int) {
  is >> massOption_ >> WWWVertex_ >> WWWWVertex_;
}

double MEGammaGamma2WW::helicityME(vector<VectorWaveFunction> & p1,
				   vector<VectorWaveFunction> & p2,
				   vector<VectorWaveFunction> & w1,
				   vector<VectorWaveFunction> & w2, bool calc) const {
  // scale (external photons so scale in couplings is 0)
  Energy2 mt(0.*GeV2);
  // matrix element to be stored
  if(calc) me_.reset(ProductionMatrixElement(PDT::Spin1,PDT::Spin1,
					     PDT::Spin1,PDT::Spin1));
  // calculate the matrix element
  double output(0.),sumdiag[2]={0.,0.};
  Complex diag[3];
  VectorWaveFunction inter;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) { 
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<3;++ohel1) { 
	for(unsigned int ohel2=0;ohel2<3;++ohel2) {
	  // first t-channel diagram
	  inter = WWWVertex_->evaluate(mt,1,w1[ohel1].particle()->CC(),
				       w1[ohel1],p1[ihel1]);
	  diag[0]= WWWVertex_->evaluate(mt,inter,w2[ohel2],p2[ihel2]);
	  //second t-channel diagram
	  inter  = WWWVertex_->evaluate(mt,1,w1[ohel1].particle()->CC(),
					w1[ohel1],p2[ihel2]);
	  diag[1] = WWWVertex_->evaluate(mt,inter,w2[ohel2],p1[ihel1]);
	  // four point diagrams
	  diag[2] = WWWWVertex_->evaluate(mt,0,p1[ihel1],p2[ihel2],
					  w1[ohel1],w2[ohel2]);
	  sumdiag[0] += norm(diag[0]);
	  sumdiag[1] += norm(diag[1]);
	  diag[0] += diag[1] + diag[2];
	  output += norm(diag[0]);
	  // store the me if needed
	  if(calc) me_(2*ihel1,2*ihel2,ohel1,ohel2) = diag[0];
	}
      }
    }
  }
  // diagrams
  DVector save;
  save.push_back(sumdiag[0]);
  save.push_back(sumdiag[1]);
  meInfo(save);
  // colour factors if needed
  if(mePartonData()[2]->coloured()) output *= 3.;
  // spin factors
  return 0.25*output;
}

void MEGammaGamma2WW::constructVertex(tSubProPtr sub) {
   // extract the particles in the hard process
   ParticleVector hard;
   hard.push_back(sub->incoming().first);
   hard.push_back(sub->incoming().second);
   hard.push_back(sub->outgoing()[0]);
   hard.push_back(sub->outgoing()[1]);
   // order of particles
   unsigned int order[4]={0,1,2,3};
   // identify the process and calculate matrix element
   vector<VectorWaveFunction> p1,p2,w1,w2;
   if(hard[order[2]]->id()>0) swap(order[2],order[3]);
   VectorWaveFunction   (p1,hard[order[0]],incoming,false,true );
   VectorWaveFunction   (p2,hard[order[1]],incoming,false,true );
   VectorWaveFunction   (w1,hard[order[2]],outgoing,true ,false);
   VectorWaveFunction   (w2,hard[order[3]],outgoing,true ,false);
   p1[1] = p1[2];
   p2[1] = p2[2];
   helicityME(p1,p2,w1,w2,true);
   // construct the vertex
   HardVertexPtr hardvertex=new_ptr(HardVertex());
   // set the matrix element for the vertex
   hardvertex->ME(me_);
   // set the pointers and to and from the vertex
   for(unsigned int ix=0;ix<4;++ix)
     hard[order[ix]]->spinInfo()->productionVertex(hardvertex);
}
