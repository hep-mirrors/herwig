// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGamma2ffAmplitude class.
//

#include "GammaGamma2ffAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

using namespace Herwig;

IBPtr GammaGamma2ffAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGamma2ffAmplitude::fullclone() const {
  return new_ptr(*this);
}

void GammaGamma2ffAmplitude::persistentOutput(PersistentOStream & os) const {
  os << process_ << vertex_;
}

void GammaGamma2ffAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> process_ >> vertex_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GammaGamma2ffAmplitude,GammaGammaAmplitude>
describeHerwigGammaGamma2ffAmplitude("Herwig::GammaGamma2ffAmplitude", "HwMEGammaGamma.so");

void GammaGamma2ffAmplitude::Init() {

  static ClassDocumentation<GammaGamma2ffAmplitude> documentation
    ("The GammaGamma2ffAmplitude class implements the amplitude for gamma gamma -> f fbar");

  static Switch<GammaGamma2ffAmplitude,int> interfaceProcess
    ("Process",
     "Which process to included",
     &GammaGamma2ffAmplitude::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all SM fermions as outgoing particles",
     0);
  static SwitchOption interfaceProcessQuarks
    (interfaceProcess,
     "Quarks",
     "All include the quarks as outgoing particles",
     1);
  static SwitchOption interfaceProcessLeptons
    (interfaceProcess,
     "Leptons",
     "Only include the leptons as outgoing particles",
     2);
  static SwitchOption interfaceProcessElectron
    (interfaceProcess,
     "Electron",
     "Only include e+e- as outgoing particles",
     5);
  static SwitchOption interfaceProcessMuon
    (interfaceProcess,
     "Muon",
     "Only include mu+mu- as outgoing particles",
     6);
  static SwitchOption interfaceProcessTau
    (interfaceProcess,
     "Tau",
     "Only include tau+tau- as outgoing particles",
     7);
  static SwitchOption interfaceProcessDown
    (interfaceProcess,
     "Down",
     "Only include d dbar as outgoing particles",
     11);
  static SwitchOption interfaceProcessUp
    (interfaceProcess,
     "Up",
     "Only include u ubar as outgoing particles",
     12);
  static SwitchOption interfaceProcessStrange
    (interfaceProcess,
     "Strange",
     "Only include s sbar as outgoing particles",
     13);
  static SwitchOption interfaceProcessCharm
    (interfaceProcess,
     "Charm",
     "Only include c cbar as outgoing particles",
     14);
  static SwitchOption interfaceProcessBottom
    (interfaceProcess,
     "Bottom",
     "Only include b bbar as outgoing particles",
     15);
  static SwitchOption interfaceProcessTop
    (interfaceProcess,
     "Top",
     "Only include t tbar as outgoing particles",
     16);

}

vector<DiagPtr> GammaGamma2ffAmplitude::getDiagrams(unsigned int iopt) const {
  vector<DiagPtr> output; output.reserve(24);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr ep = getParticleData(ParticleID::eplus );
  tcPDPtr em = getParticleData(ParticleID::eminus);
  for(int ix=1;ix<17;++ix) {
    // increment counter to switch between quarks and leptons
    if(ix==7) ix+=4;
    if(ix>11&&ix%2==0) ++ix;
    // is it a valid quark process
    bool quark = ix<=6 && (process_==0 || process_==1 || process_-10==ix);
    // is it a valid lepton process
    bool lepton= ix>=11 && ix<=16 && 
      (process_==0 || process_==2 || ( (ix-9)/2 ==process_-4) );
    // only add diagram if allowed
    if(!quark && !lepton) continue;
    tcPDPtr lm = getParticleData(ix);
    tcPDPtr lp = lm->CC();
    if(iopt==0) {
      // first t-channel
      output.push_back(new_ptr((Tree2toNDiagram(3),gamma,lp,gamma,1,lm, 2,lp,-1)));
      // interchange
      output.push_back(new_ptr((Tree2toNDiagram(3),gamma,lp,gamma,2,lm, 1,lp,-2)));
    }
    else {
      // first t-channel
      output.push_back(new_ptr((Tree2toNDiagram(5),em,gamma,lp,gamma,ep, 1, em, 4, ep, 2,lm, 3,lp,-1)));
      // interchange							              
      output.push_back(new_ptr((Tree2toNDiagram(5),em,gamma,lp,gamma,ep, 1, em, 4, ep, 3,lm, 2,lp,-2)));
    }
  }
  return output;
}

double GammaGamma2ffAmplitude::me2(const vector<VectorWaveFunction> & v1,
				   const vector<VectorWaveFunction> & v2,
				   const Energy2 & , const Energy2 & ,
				   const Energy2 & , 
				   const vector<Lorentz5Momentum> & momenta,
				   const cPDVector & partons,
				   DVector & sumdiag) const {
  SpinorBarWaveFunction    fw(momenta[0],partons[0],outgoing);
  SpinorWaveFunction    fbarw(momenta[1],partons[1],outgoing);
  vector<SpinorBarWaveFunction> f;
  vector<SpinorWaveFunction> fbar;
  for(unsigned int ix=0;ix<2;++ix) {
    fw.reset(ix);f.push_back(fw);
    fbarw.reset(ix);fbar.push_back(fbarw);
  }
  // calculate the matrix element
  double output(0.);
  sumdiag={0.,0.};
  Complex diag[2];
  Energy2 mt(ZERO);
  SpinorWaveFunction inters;
  for(unsigned int ihel1=0;ihel1<v1.size();++ihel1) { 
    for(unsigned int ihel2=0;ihel2<v2.size();++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) { 
      	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  // first t-channel diagram
	  inters =vertex_->evaluate(mt,1,fbar[ohel2].particle()->CC(),
				    fbar[ohel2],v2[ihel2]);
	  diag[0]=vertex_->evaluate(mt,inters,f[ohel1],v1[ihel1]);
	  //second t-channel diagram
	  inters =vertex_->evaluate(mt,1,fbar[ohel2].particle()->CC(),
				    fbar[ohel2],v1[ihel1]);
	  diag[1]=vertex_->evaluate(mt,inters,f[ohel1],v2[ihel2]);
	  sumdiag[0] += norm(diag[0]);
	  sumdiag[1] += norm(diag[1]);
	  diag[0] += diag[1];
	  output += norm(diag[0]);
      	}
      }
    }
  }
  // colour factors if needed
  if(partons[0]->coloured()) output *= 3.;
  // // code to test vs the analytic result
  // Energy2 m2 = sqr(f[0].particle()->mass());
  // Energy2 tm = (v1[0].momentum()+f   [0].momentum()).m2()-m2;
  // Energy2 um = (v1[0].momentum()+fbar[0].momentum()).m2()-m2;
  // double test = 8.*um/tm+8.*tm/um- 32*m2/tm - 32*m2/um
  //   -32*sqr(double(m2/tm)) - 64*sqr(m2)/tm/um - 32*sqr(double(m2/um));
  // test *= sqr(4.*Constants::pi*SM().alphaEM());
  // if(mePartonData()[2]->coloured()) test *= 3.;
  // cerr << "testing ME " << (output-test)/(output+test) << "\n"; 
  // spin factors
  return 0.25*output;
}

void GammaGamma2ffAmplitude::doinit() {
  GammaGammaAmplitude::doinit();
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  // do the initialisation
  if(!hwsm)
    throw InitException() << "Must be the Herwig StandardModel class in "
 			  << "MEGammaGamma2ff::doinit" << Exception::abortnow;
  vertex_ = hwsm->vertexFFP();
}

Selector<const ColourLines *>
GammaGamma2ffAmplitude::colourGeometries(unsigned int iopt, const cPDVector & partons, tcDiagPtr ) const {
  static ColourLines c1("");
  static ColourLines c2("4 -2 -5"); 
  static ColourLines c3("8 -3 -9"); 
  Selector<const ColourLines *> sel;
  if(iopt==0) {
    if(partons[2]->coloured()) sel.insert(1.0, &c2);
    else                       sel.insert(1.0, &c1);
  }
  else {
    if(partons[4]->coloured()) sel.insert(1.0, &c3);
    else                       sel.insert(1.0, &c1);
  }
  return sel;
}


// #include "MEGammaGamma2ff.h"
// #include "ThePEG/Utilities/DescribeClass.h"
// #include "ThePEG/Interface/ClassDocumentation.h"
// #include "ThePEG/Persistency/PersistentOStream.h"
// #include "ThePEG/Persistency/PersistentIStream.h"
// #include "ThePEG/PDT/EnumParticles.h"
// #include "ThePEG/MatrixElement/Tree2toNDiagram.h"
// #include "ThePEG/Interface/Switch.h"
// #include "Herwig/Models/StandardModel/StandardModel.h"
// #include "ThePEG/Handlers/StandardXComb.h"
// #include "Herwig/MatrixElement/HardVertex.h"

// using namespace Herwig;

// MEGammaGamma2ff::MEGammaGamma2ff()  {
//   massOption(vector<unsigned int>(2,1));
// }

// Energy2 MEGammaGamma2ff::scale() const {
//   return 2.*sHat()*tHat()*uHat()/(sqr(sHat())+sqr(tHat())+sqr(uHat()));
// }




// void MEGammaGamma2ff::constructVertex(tSubProPtr sub) {
//   // extract the particles in the hard process
//   ParticleVector hard;
//   hard.push_back(sub->incoming().first);
//   hard.push_back(sub->incoming().second);
//   hard.push_back(sub->outgoing()[0]);
//   hard.push_back(sub->outgoing()[1]);
//   // order of particles
//   unsigned int order[4]={0,1,2,3};
//   // identify the process and calculate matrix element
//   vector<VectorWaveFunction> v1,v2;
//   if(hard[order[2]]->id()<0) swap(order[2],order[3]);
//   vector<SpinorWaveFunction> fbar;
//   vector<SpinorBarWaveFunction>  f;
//   VectorWaveFunction   (v1  ,hard[order[0]],incoming,false,true);
//   VectorWaveFunction   (v2  ,hard[order[1]],incoming,false,true);
//   SpinorBarWaveFunction(f   ,hard[order[2]],outgoing,true );
//   SpinorWaveFunction   (fbar,hard[order[3]],outgoing,true );
//   v1[1]=v1[2];
//   v2[1]=v2[2];
//   helicityME(v1,v2,f,fbar,true);
//   // construct the vertex
//   HardVertexPtr hardvertex=new_ptr(HardVertex());
//   // set the matrix element for the vertex
//   hardvertex->ME(me_);
//   // set the pointers and to and from the vertex
//   for(unsigned int ix=0;ix<4;++ix)
//     hard[order[ix]]->spinInfo()->productionVertex(hardvertex);
// }

Energy GammaGamma2ffAmplitude::generateW(double r, const tcPDVector & partons,Energy Wmax,Energy2 & jacW) {
  Energy Wmin = 2.*partons[0]->constituentMass();
  Energy W = Wmin+r*(Wmax-Wmin);
  jacW = 2.*W*(Wmax-Wmin);
  return W;
}
