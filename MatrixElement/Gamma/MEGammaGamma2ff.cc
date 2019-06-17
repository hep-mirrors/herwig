// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGammaGamma2ff class.
//

#include "MEGammaGamma2ff.h"
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

MEGammaGamma2ff::MEGammaGamma2ff() : process_(0) {
  massOption(vector<unsigned int>(2,1));
}

void MEGammaGamma2ff::doinit() {
  HwMEBase::doinit();
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!hwsm)
    throw InitException() << "Must be the Herwig StandardModel class in "
			  << "MEGammaGamma2ff::doinit" << Exception::abortnow;
  vertex_ = hwsm->vertexFFP();
}

void MEGammaGamma2ff::getDiagrams() const {
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
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
    // first t-channel
    add(new_ptr((Tree2toNDiagram(3),gamma,lp,gamma,1,lm, 2,lp,-1)));
    // interchange
    add(new_ptr((Tree2toNDiagram(3),gamma,lp,gamma,2,lm, 1,lp,-2)));
  }
}

Energy2 MEGammaGamma2ff::scale() const {
  return 2.*sHat()*tHat()*uHat()/(sqr(sHat())+sqr(tHat())+sqr(uHat()));
}

double MEGammaGamma2ff::me2() const {
  VectorWaveFunction      p1w(rescaledMomenta()[0],mePartonData()[0],incoming);
  VectorWaveFunction      p2w(rescaledMomenta()[1],mePartonData()[1],incoming);
  SpinorBarWaveFunction    fw(rescaledMomenta()[2],mePartonData()[2],outgoing);
  SpinorWaveFunction    fbarw(rescaledMomenta()[3],mePartonData()[3],outgoing);
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
  // scale (external photons so scale in couplings is 0)
  Energy2 mt(0.*GeV2);
  // matrix element to be stored
  if(calc) me_.reset(ProductionMatrixElement(PDT::Spin1,PDT::Spin1,
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
	  inters =vertex_->evaluate(mt,1,fbar[ohel2].particle()->CC(),
				    fbar[ohel2],p2[ihel2]);
	  diag[0]=vertex_->evaluate(mt,inters,f[ohel1],p1[ihel1]);
	  //second t-channel diagram
	  inters =vertex_->evaluate(mt,1,fbar[ohel2].particle()->CC(),
				    fbar[ohel2],p1[ihel1]);
	  diag[1]=vertex_->evaluate(mt,inters,f[ohel1],p2[ihel2]);
	  sumdiag[0] += norm(diag[0]);
	  sumdiag[1] += norm(diag[1]);
	  diag[0] += diag[1];
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
  // // code to test vs the analytic result
  // Energy2 m2 = sqr(f[0].particle()->mass());
  // Energy2 tm = (p1[0].momentum()+f   [0].momentum()).m2()-m2;
  // Energy2 um = (p1[0].momentum()+fbar[0].momentum()).m2()-m2;
  // double test = 8.*um/tm+8.*tm/um- 32*m2/tm - 32*m2/um
  //   -32*sqr(double(m2/tm)) - 64*sqr(m2)/tm/um - 32*sqr(double(m2/um));
  // test *= sqr(4.*Constants::pi*SM().alphaEM());
  // if(mePartonData()[2]->coloured()) test *= 3.;
  // cerr << "testing ME " << (output-test)/(output+test) << "\n"; 
  // spin factors
  return 0.25*output;
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
  os << process_ << vertex_;
}

void MEGammaGamma2ff::persistentInput(PersistentIStream & is, int) {
  is >> process_ >> vertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEGammaGamma2ff,HwMEBase>
describeHerwigMEGammaGamma2ff("Herwig::MEGammaGamma2ff", "HwMEGammaGamma.so");

void MEGammaGamma2ff::Init() {

  static ClassDocumentation<MEGammaGamma2ff> documentation
    ("The MEGammaGamma2ff class implements the matrix elemments for"
     " direct fermion-antifermion production in gamma-gamma collisions.");

  static Switch<MEGammaGamma2ff,int> interfaceProcess
    ("Process",
     "Which process to included",
     &MEGammaGamma2ff::process_, 0, false, false);
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

void MEGammaGamma2ff::constructVertex(tSubProPtr sub) {
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
  vector<SpinorWaveFunction> fbar;
  vector<SpinorBarWaveFunction>  f;
  VectorWaveFunction   (p1  ,hard[order[0]],incoming,false,true);
  VectorWaveFunction   (p2  ,hard[order[1]],incoming,false,true);
  SpinorBarWaveFunction(f   ,hard[order[2]],outgoing,true );
  SpinorWaveFunction   (fbar,hard[order[3]],outgoing,true );
  p1[1]=p1[2];
  p2[1]=p2[2];
  helicityME(p1,p2,f,fbar,true);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me_);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix)
    hard[order[ix]]->spinInfo()->productionVertex(hardvertex);
}
