// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGammaGamma2ff class.
//

#include "MEGammaGamma2ff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"

using namespace Herwig;

void MEGammaGamma2ff::doinit() throw(InitException) {
  ME2to2Base::doinit();
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _theFFPVertex = hwsm->vertexFFP();
  }
  else
    throw InitException() << "Must be the Herwig++ StandardModel class in "
			  << "MEGammaGamma2ff::doinit" << Exception::abortnow;
}

void MEGammaGamma2ff::getDiagrams() const {
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  for(unsigned int ix=1;ix<17;++ix) {
    // increment counter to switch between quarks and leptons
    if(ix==7) ix+=4;
    if(ix>11&&ix%2==0) ++ix;
    // is it a valid quark process
    bool quark = ix<=6 && (_process==0 || _process==1 || _process-10==ix);
    // is it a valid lepton process
    bool lepton= ix>=11 && ix<=16 && 
      (_process==0 || _process==2 || ( (ix-9)/2 ==_process-4) );
    // only add diagram if allowed
    if(!quark && !lepton) continue;
    tcPDPtr lm = getParticleData(ix);
    tcPDPtr lp = lm->CC();
    // first t-channel
    add(new_ptr((Tree2toNDiagram(3),gamma,lp,gamma,1,lm, 2,lp,-1)));
    // interchange
    add(new_ptr((Tree2toNDiagram(3),gamma,lp,gamma,2,lm, 1,lp,-2)));
  }
}

Energy2 MEGammaGamma2ff::scale() const {
  return sHat();
}

double MEGammaGamma2ff::me2() const {
  VectorWaveFunction      p1w(meMomenta()[0],mePartonData()[0],incoming);
  VectorWaveFunction      p2w(meMomenta()[1],mePartonData()[1],incoming);
  SpinorBarWaveFunction    fw(meMomenta()[2],mePartonData()[2],outgoing);
  SpinorWaveFunction    fbarw(meMomenta()[3],mePartonData()[3],outgoing);
  vector<VectorWaveFunction> p1,p2;
  vector<SpinorBarWaveFunction> f;
  vector<SpinorWaveFunction> fbar;
  for(unsigned int ix=0;ix<2;++ix) {
      p1w.reset(2*ix);p1.push_back(p1w);
      p2w.reset(2*ix);p2.push_back(p2w);
      fw.reset(ix);f.push_back(fw);
      fbarw.reset(ix);fbar.push_back(fbarw);
    }
  // calculate the matrix element
  return helicityME(p1,p2,f,fbar,false);
}

unsigned int MEGammaGamma2ff::orderInAlphaS() const {
  return 0;
}

unsigned int MEGammaGamma2ff::orderInAlphaEW() const {
  return 2;
}

double MEGammaGamma2ff::helicityME(vector<VectorWaveFunction> &p1,
				   vector<VectorWaveFunction> &p2,
				   vector<SpinorBarWaveFunction> & f,
				   vector<SpinorWaveFunction> & fbar, bool calc) const {
  // scale
  Energy2 mt(scale());
  // matrix element to be stored
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1,PDT::Spin1,
					     PDT::Spin1Half,PDT::Spin1Half));
  // calculate the matrix element
  double output(0.),sumdiag[2]={0.,0.};
  Complex diag[2];
  SpinorWaveFunction inters;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) { 
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) { 
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  //first t-channel diagram
	  inters =_theFFPVertex->evaluate(mt,5,fbar[ohel2].getParticle(),
					  fbar[ohel2],p2[ihel2]);
	  diag[0]=_theFFPVertex->evaluate(mt,inters,f[ohel1],p1[ihel1]);
	  //second t-channel diagram
	  inters =_theFFPVertex->evaluate(mt,5,fbar[ohel2].getParticle(),
					  fbar[ohel2],p1[ihel1]);
	  diag[1]=_theFFPVertex->evaluate(mt,inters,f[ohel1],p2[ihel2]);
	  sumdiag[0] += norm(diag[0]);
	  sumdiag[1] += norm(diag[1]);
	  diag[0] += diag[1];
	  output += norm(diag[0]);
	  // store the me if needed
	  if(calc) _me(2*ihel1,2*ihel2,ohel1,ohel2) = diag[0];
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
  return output/4.;
}

Selector<MEBase::DiagramIndex>
MEGammaGamma2ff::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 )      sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
  return sel;
}

Selector<const ColourLines *>
MEGammaGamma2ff::colourGeometries(tcDiagPtr ) const {
  static ColourLines c1("");
  static ColourLines c2("4 -2 -5"); 
  Selector<const ColourLines *> sel;
  if(mePartonData()[2]->coloured()) sel.insert(1.0, &c2);
  else                              sel.insert(1.0, &c1);
  return sel;
}

void MEGammaGamma2ff::persistentOutput(PersistentOStream & os) const {
  os << _process << _theFFPVertex;
}

void MEGammaGamma2ff::persistentInput(PersistentIStream & is, int) {
  is >> _process >> _theFFPVertex;
}

ClassDescription<MEGammaGamma2ff> MEGammaGamma2ff::initMEGammaGamma2ff;
// Definition of the static class description member.

void MEGammaGamma2ff::Init() {

  static ClassDocumentation<MEGammaGamma2ff> documentation
    ("The MEGammaGamma2ff class implements the matrix elemments for"
     " direct fermion-antifermion production in gamma-gamma collisions.");

  static Switch<MEGammaGamma2ff,unsigned int> interfaceProcess
    ("Process",
     "Which process to included",
     &MEGammaGamma2ff::_process, 0, false, false);
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
