// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2Higgs class.
//

#include "MEPP2Higgs.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/PDT/GenericMassGenerator.h"
#include "Herwig++/Helicity/Correlations/HardVertex.h"
#include "ThePEG/Cuts/Cuts.h"

using namespace Herwig;

ClassDescription<MEPP2Higgs> MEPP2Higgs::initMEPP2Higgs;
// Definition of the static class description member.

void MEPP2Higgs::Init() {

  static ClassDocumentation<MEPP2Higgs> documentation
    ("The MEPP2Higgs class implements the matrix elements for q qbar/ gg -> H");

  static Switch<MEPP2Higgs,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2Higgs::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcessGluon
    (interfaceProcess,
     "Gluon",
     "Only include gg -> H",
     1);
  static SwitchOption interfaceProcessQuark
    (interfaceProcess,
     "Quark",
     "Only include q qbar -> H",
     2);

  static Parameter<MEPP2Higgs,unsigned int> interfaceMinimumInLoop
    ("MinimumInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &MEPP2Higgs::_minloop, 6, 4, 6,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,unsigned int> interfaceMaximumInLoop
    ("MaximumInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &MEPP2Higgs::_maxloop, 6, 4, 6,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,unsigned int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The minimum flavour of the incoming quarks in the hard process",
     &MEPP2Higgs::_minFlavour, 3, 3, 5,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the incoming quarks in the hard process",
     &MEPP2Higgs::_maxFlavour, 5, 3, 5,
     false, false, Interface::limited);

  static Switch<MEPP2Higgs,bool> interfaceRunningMass
    ("RunningMass",
     "Use the running quark mass in the loop diagrams",
     &MEPP2Higgs::_runLoop, false, false, false);
  static SwitchOption interfaceRunningMassRunning
    (interfaceRunningMass,
     "Running",
     "Use the running mass",
     true);
  static SwitchOption interfaceRunningMassPole
    (interfaceRunningMass,
     "Pole",
     "Use the pole mass",
     false);

  static Switch<MEPP2Higgs,bool> interfaceLineShape
    ("LineShape",
     "Option for the Higgs lineshape",
     &MEPP2Higgs::_lineshape, false, false, false);
  static SwitchOption interfaceLineShapeMassGenerator
    (interfaceLineShape,
     "MassGenerator",
     "Use the mass generator if available",
     true);
  static SwitchOption interfaceLineShapeBreitWigner
    (interfaceLineShape,
     "BreitWigner",
     "Use a Breit-Wigner with the naive running width",
     true);
}

Energy2 MEPP2Higgs::scale() const {
  return sHat();
}

int MEPP2Higgs::nDim() const {
  return 0;
}

unsigned int MEPP2Higgs::orderInAlphaS() const {
  return 2;
}

unsigned int MEPP2Higgs::orderInAlphaEW() const {
  return 1;
}

void MEPP2Higgs::persistentOutput(PersistentOStream & os) const {
  os << _process << _minloop << _maxloop << _minFlavour << _maxFlavour
     << _runLoop << _lineshape << _massgen << _theFFHVertex << _runningmass;  
}

void MEPP2Higgs::persistentInput(PersistentIStream & is, int) {
  is >> _process >> _minloop >> _maxloop >> _minFlavour >> _maxFlavour
     >> _runLoop >> _lineshape >> _massgen >> _theFFHVertex >> _runningmass;
}

void MEPP2Higgs::getDiagrams() const {
  tcPDPtr h0=getParticleData(ParticleID::h0);
  // gg -> H process
  if(_process==0||_process==1) {
    tcPDPtr g=getParticleData(ParticleID::g);
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, h0, -1)));
  }
  // q qbar -> H processes
  if(_process==0||_process==2) {
    for (unsigned int i = _minFlavour; i <= _maxFlavour; ++i ) {
      tcPDPtr q = getParticleData(i);
      tcPDPtr qb = q->CC();
    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, h0, -2)));
    }
  }
}

Selector<MEBase::DiagramIndex>
MEPP2Higgs::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    sel.insert(1.0, i);
  return sel;
}

Selector<const ColourLines *>
MEPP2Higgs::colourGeometries(tcDiagPtr diag) const {
  // This example corresponds to the diagrams specified in the example
  // in the getDiagrams() function.

  static const ColourLines gluon("1 -2,2 -1");
  static const ColourLines quark("1 -2");

  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 ) sel.insert(1.0, &gluon);
  else                    sel.insert(1.0, &quark);
  return sel;
}

void MEPP2Higgs::doinit() throw(InitException) {
  MEBase::doinit();
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!hwsm) InitException() << "Must be the Herwig++ StandardModel class in "
			    << "MEqq2gZ2ll::doinit" << Exception::abortnow;
  // get the vertex
  _theFFHVertex = hwsm->vertexFFH();
  // running mass option
  _runningmass = hwsm->massPtr();
  // mass generator
  tMassGenPtr mass=getParticleData(ParticleID::h0)->massGenerator();
  if(mass) {
    _massgen=dynamic_ptr_cast<GenericMassGeneratorPtr>(mass);
  }
}

bool MEPP2Higgs::generateKinematics(const double *) {
  Lorentz5Momentum pout=meMomenta()[0]+meMomenta()[1];
  pout.rescaleMass();
  meMomenta()[2].setMass(pout.mass());
  meMomenta()[2] = LorentzMomentum(pout.x(),pout.y(),pout.z(),pout.t());
  jacobian(1.0);
  // check passes all the cuts
  vector<LorentzMomentum> out(1,meMomenta()[2]);
  tcPDVector tout(1,mePartonData()[2]);
  // return true if passes the cuts
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}

CrossSection MEPP2Higgs::dSigHatDR() const {
  // compute the weight for the Breit-Wigner piece and final-state
  tcPDPtr h0=getParticleData(ParticleID::h0);
  InvEnergy2 wgt;
  Energy  M(h0->mass()),G(h0->width());
  Energy2 M2(sqr(M)),GM(G*M);
  if(_massgen&&_lineshape) {
    // DGRELL units?
    wgt =Constants::pi*UnitRemoval::E2*_massgen->weight(sqrt(sHat()))/(sqr(sHat()-M2)+sqr(GM));
  }
  else {
    wgt = sHat()*G/M/(sqr(sHat()-M2)+sqr(sHat()*G/M));
  }
  return me2()*jacobian()*wgt*sqr(hbarc);                         
}

double MEPP2Higgs::me2() const {
  useMe();
  double output(0.);
  // g g to g
  if(mePartonData()[0]->id()==ParticleID::g&&
     mePartonData()[1]->id()==ParticleID::g) {
    VectorWaveFunction glin1(meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction glin2(meMomenta()[1],mePartonData()[1],incoming);
    ScalarWaveFunction  hout(meMomenta()[2],mePartonData()[2],outgoing);
    vector<VectorWaveFunction> g1,g2;
    for(unsigned int ix=0;ix<2;++ix) {
      glin1.reset(2*ix);g1.push_back(glin1);
      glin2.reset(2*ix);g2.push_back(glin2);
    }
    // calculate the matrix element
    output = ggME(g1,g2,hout,false); 
  }
  // q qbar to H g
  else if(mePartonData()[0]->id()==-mePartonData()[1]->id()) {
    // calculate the spinors and polarization vectors
    vector<SpinorWaveFunction> fin;
    vector<SpinorBarWaveFunction>  ain;
    SpinorWaveFunction    qin (meMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction qbin(meMomenta()[1],mePartonData()[1],incoming);
    ScalarWaveFunction    hout(meMomenta()[2],mePartonData()[2],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)    ; fin.push_back(  qin);
      qbin.reset(ix)   ; ain.push_back( qbin);
    }
    // calculate the matrix element
    output = qqbarME(fin,ain,hout,false); 
  }
  else
    throw Exception() << "Unknown subprocess in MEPP2Higgs::me2()" 
		      << Exception::runerror;
  // DGRELL return the answer with factor to make dimensionless
  return UnitRemoval::E2*output/sHat();
}

double MEPP2Higgs::ggME(vector<VectorWaveFunction> g1, 
			vector<VectorWaveFunction> g2,
			ScalarWaveFunction &, bool calc) const {
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1,PDT::Spin0);
  // get the kinematic invariants
  Energy2 s(sHat()),et(scale());
  // calculate the loop functions
  Complex A1s(0.);
  for(unsigned int ix=_minloop;ix<=_maxloop;++ix) {
    Energy mf;
    if(_runningmass&&_runLoop) mf = _runningmass->value(s,getParticleData(ix));
    else                       mf = getParticleData(ix)->mass();
    mf= getParticleData(ix)->mass();
    A1s+=A1(s,sqr(mf));
  }
  // prefactors
  using Constants::pi;
  double g(sqrt(4.*pi*SM().alphaEM(et)/SM().sin2ThetaW()));
  double gs(sqrt(4.*pi*SM().alphaS(et)));
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  double pre=UnitRemoval::InvE*(g*s*sqr(gs)/mw/32/sqr(Constants::pi));
  // compute the matrix element
  double output(0.);
  complex<Energy> w1p2[2],w2p1[2];
  for(unsigned int ihel=0;ihel<2;++ihel) {
    w1p2[ihel]=g1[ihel].wave().dot(g2[0].getMomentum());
    w2p1[ihel]=g2[ihel].wave().dot(g1[0].getMomentum());
  }
  Complex ii(0.,1.);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      Complex wdot=g1[ihel1].wave().dot(g2[ihel2].wave());
      Complex me=-ii*A1s*pre*(wdot-2.*w1p2[ihel1]*w2p1[ihel2]/s);
      output+=real(me*conj(me));
      // matrix element if needed
      if(calc) newme(2*ihel1,2*ihel2,0)=me;
    }
  }
  // analytic form for testing
//    double test = sqr(SM().alphaS(et))*SM().alphaEM(et)/SM().sin2ThetaW()/pi*
//      sqr(sHat()/mw)*real(A1s*conj(A1s));
//    test/=4*64;
//    cerr << "testing matrix element " 
// 	<< test << " " 
// 	<< 8.*output/4./64. << " " 
// 	<< test/(8.*output/4./64.)
// 	<< endl;
  // final colour and spin factors (8/64) colour 1/4 spin
  if(calc) _me.reset(newme);
  return output/32.;
}

double MEPP2Higgs::qqbarME(vector<SpinorWaveFunction>    & fin,
			   vector<SpinorBarWaveFunction> & ain,
			   ScalarWaveFunction & hout,bool calc) const {
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0);
  // get the kinematic invariants
  Energy2 et(scale());
  // calculate the loop function
  double output(0.);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      Complex me=_theFFHVertex->evaluate(et,fin[ihel1],ain[ihel2],hout);
      output+=real(me*conj(me));
      if(calc) newme(ihel1,ihel2,0)=me;
    }
  }
  // final colour/spin factors
  if(calc) _me.reset(newme);
  return 3.*output/9./4.;
}

void MEPP2Higgs::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  if(hard[0]->id()<hard[1]->id()) swap(hard[0],hard[1]);
  // identify the process and calculate the matrix element
  if(hard[0]->id()==ParticleID::g&&hard[1]->id()==ParticleID::g)
    {
      vector<VectorWaveFunction> g1,g2;
      vector<SpinorBarWaveFunction> q;
      vector<SpinorWaveFunction> qbar;
      VectorWaveFunction(     g1,hard[0],incoming,false,true,true);
      VectorWaveFunction(     g2,hard[1],incoming,false,true,true);
      ScalarWaveFunction    hout(hard[2],outgoing,true,true);
      g1[1]=g1[2];g2[1]=g2[2];
      ggME(g1,g2,hout,true);
    }
  // q qbar -> q qbar
  else
    {
      vector<SpinorWaveFunction>    q1;
      vector<SpinorBarWaveFunction> q2;
      SpinorWaveFunction(   q1,hard[0],incoming,false,true);
      SpinorBarWaveFunction(q2,hard[1],incoming,false,true);
      ScalarWaveFunction  hout(hard[2],outgoing,true,true);
      qqbarME(q1,q2,hout,true);
    }
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());  
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<3;++ix) {
    dynamic_ptr_cast<SpinfoPtr>(hard[ix]->spinInfo())->
      setProductionVertex(hardvertex);
  }
}





























