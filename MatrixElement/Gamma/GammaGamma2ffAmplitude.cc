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
  tcPDPtr pp = getParticleData(ParticleID::pplus ); 
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
      // first t-channel
      output.push_back(new_ptr((Tree2toNDiagram(5),pp,gamma,lp,gamma,pp, 1, pp, 4, pp, 2,lm, 3,lp,-1)));
      // interchange							              
      output.push_back(new_ptr((Tree2toNDiagram(5),pp,gamma,lp,gamma,pp, 1, pp, 4, pp, 3,lm, 2,lp,-2)));
    }
  }
  return output;
}

ProductionMatrixElement
GammaGamma2ffAmplitude::helicityAmplitude(const vector<VectorWaveFunction> & v1,
					  const vector<VectorWaveFunction> & v2,
					  const vector<SpinorBarWaveFunction> & f,
					  const vector<SpinorWaveFunction>    & fbar,
					  double & output, DVector & dweight) const {
  dweight = {0.,0.};
  ProductionMatrixElement me;
  if(v1.size()==4&&v2.size()==4) {
    me = ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
  				 PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1Half);
  }
  else if(v1.size()==2&&v2.size()==2) {
    me = ProductionMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half);
  }
  else
    assert(false);
  Complex diag[2];
  Energy2 mt(ZERO);
  SpinorWaveFunction inters;
  for(unsigned int ih1A=0;ih1A<v1.size()/2;++ih1A) {
    for(unsigned int ih1B=0;ih1B<2;++ih1B) {
      unsigned int ih1 = 2*ih1A+ih1B;
      for(unsigned int ih2A=0;ih2A<v1.size()/2;++ih2A) {
  	for(unsigned int ih2B=0;ih2B<2;++ih2B) {
	  unsigned int ih2 = 2*ih2A+ih2B;
	  for(unsigned int ohel1=0;ohel1<2;++ohel1) { 
	    for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	      // first t-channel diagram
	      inters =vertex_->evaluate(mt,1,fbar[ohel2].particle()->CC(),
					fbar[ohel2],v2[ih2]);
	      diag[0]=vertex_->evaluate(mt,inters,f[ohel1],v1[ih1]);
	      //second t-channel diagram
	      inters =vertex_->evaluate(mt,1,fbar[ohel2].particle()->CC(),
					fbar[ohel2],v1[ih1]);
	      diag[1]=vertex_->evaluate(mt,inters,f[ohel1],v2[ih2]);
	      dweight[0] += norm(diag[0]);
	      dweight[1] += norm(diag[1]);
	      Complex amp = diag[0]+diag[1];
	      output += norm(amp);
	      if(v1.size()==4) {
		me(ih1A,ih1B,ih2A,ih2B,ohel1,ohel2) = amp;
	      }
	      else {
		me(2*ih1B,2*ih2B,ohel1,ohel2) = amp;
	      }
	    }
	  }
	}
      }
    }
  }
  return me;
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
  helicityAmplitude(v1,v2,f,fbar,output,sumdiag);
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

// MEGammaGamma2ff::MEGammaGamma2ff()  {
//   massOption(vector<unsigned int>(2,1));
// }

// Energy2 MEGammaGamma2ff::scale() const {
//   return 2.*sHat()*tHat()*uHat()/(sqr(sHat())+sqr(tHat())+sqr(uHat()));
// }

ProductionMatrixElement GammaGamma2ffAmplitude::me(const vector<VectorWaveFunction> & v1,
						   const vector<VectorWaveFunction> & v2,
						   tParticleVector & particles) const {
  vector<SpinorWaveFunction> fbar;
  vector<SpinorBarWaveFunction>  f;
  SpinorBarWaveFunction(f   ,particles[0],outgoing,true );
  SpinorWaveFunction   (fbar,particles[1],outgoing,true );
  DVector dwgt;
  double output(0);
  return helicityAmplitude(v1,v2,f,fbar,output,dwgt);
}

Energy GammaGamma2ffAmplitude::generateW(double r, const tcPDVector & partons,Energy Wmax,Energy2 & jacW, Energy2) {
  Energy Wmin = 2.*partons[0]->constituentMass();
  //Energy Wmin=3.*GeV;
  Energy W = Wmin*pow(Wmax/Wmin,r);
  jacW = 2.*sqr(W)*log(Wmax/Wmin);
  return W;
}
