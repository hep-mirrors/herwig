// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2ZPrime2ff class.
//

#include "MEqq2ZPrime2ff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "RadiativeZPrimeModel.h"

using namespace RadiativeZPrime;

MEqq2ZPrime2ff::MEqq2ZPrime2ff() : _maxflavour(5), _gammaZ(0), _process(0) {
  massOption(vector<unsigned int>(2,1));
}

void MEqq2ZPrime2ff::doinit() {
  HwMEBase::doinit();
  _zPrime = getParticleData(32); 
  _z0     = getParticleData(ThePEG::ParticleID::Z0);
  _gamma  = getParticleData(ThePEG::ParticleID::gamma);
  tcSMPtr sm = generator()->standardModel();
  tcRadiativeZPrimeModelPtr model = 
    dynamic_ptr_cast<tcRadiativeZPrimeModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the RadiativeZPrimeModel in "
			       << "MEqq2ZPrime2ff::doinit()" << Exception::abortnow;
  _theFFZVertex      = model->vertexFFZ();
  _theFFPVertex      = model->vertexFFP();
  _theFFZPrimeVertex = model->vertexFFZPrime();
}

void MEqq2ZPrime2ff::getDiagrams() const {
  // which intermediates to include
  bool gamma  = _gammaZ==0 || _gammaZ==1;
  bool Z0     = _gammaZ==0 || _gammaZ==2;
  bool Zprime = _gammaZ==0 || _gammaZ==3;
  // loop over the processes we need
  for(unsigned int ix=1;ix<17;++ix) {
    // increment counter to switch between quarks and leptons
    if(ix==7) ix+=4;
    // is it a valid quark process
    bool quark = ix<=6 && (_process==0 || _process==1 || _process-10==ix);
    // is it a valid lepton process
    bool lepton= ix>=11 && ix<=16 
      && (_process==0
	  || _process==2
	  || (_process==3 && ix%2==1)
	  || (_process==4 && ix%2==0)
	  || (ix%2==0 && (ix-10)/2==_process-7)
	  || (ix%2==1 && (ix-9)/2 ==_process-4)
	  );
    // if not a valid process continue
    if(!(quark||lepton)) continue;
    tcPDPtr lm = getParticleData(long(ix));
    tcPDPtr lp = lm->CC();
    for(unsigned int i = 1; i <= _maxflavour; ++i) {
      tcPDPtr q  = getParticleData(long(i));
      tcPDPtr qb = q->CC();
      if(Z0)     add(new_ptr((Tree2toNDiagram(2), q, qb, 1, _z0    , 3, lm, 3, lp, -1)));
      if(gamma)  add(new_ptr((Tree2toNDiagram(2), q, qb, 1, _gamma , 3, lm, 3, lp, -2)));
      if(Zprime) add(new_ptr((Tree2toNDiagram(2), q, qb, 1, _zPrime, 3, lm, 3, lp, -3)));
    }
  }
}

Energy2 MEqq2ZPrime2ff::scale() const {
  return sHat();
}

double MEqq2ZPrime2ff::me2() const {
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction       q(meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction qbar(meMomenta()[1],mePartonData()[1],incoming);
  SpinorBarWaveFunction    f(meMomenta()[2],mePartonData()[2],outgoing);
  SpinorWaveFunction    fbar(meMomenta()[3],mePartonData()[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q.reset(ix)   ; fin.push_back(q);
    qbar.reset(ix); ain.push_back(qbar);
    f.reset(ix)   ;fout.push_back(f);
    fbar.reset(ix);aout.push_back(fbar);
  }
  return qqME(fin,ain,fout,aout,false);
}

unsigned int MEqq2ZPrime2ff::orderInAlphaS() const {
  return 0;
}

unsigned int MEqq2ZPrime2ff::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEqq2ZPrime2ff::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
    else if ( diags[i]->id() == -3 ) sel.insert(meInfo()[2], i);
  }
  return sel;
}

Selector<const ColourLines *>
MEqq2ZPrime2ff::colourGeometries(tcDiagPtr diag) const {
  static const ColourLines c1("1 -2");
  static const ColourLines c2("1 -2,4 -5");
  Selector<const ColourLines *> sel;
  if(abs(mePartonData()[2]->id())<=6) sel.insert(1.0, &c2);
  else                                sel.insert(1.0, &c1);
  return sel;
}


void MEqq2ZPrime2ff::persistentOutput(PersistentOStream & os) const {
  os << _theFFZPrimeVertex << _theFFZVertex << _theFFPVertex << _zPrime << _z0
     << _gamma << _maxflavour << _gammaZ << _process;
}

void MEqq2ZPrime2ff::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZPrimeVertex >> _theFFZVertex >> _theFFPVertex >> _zPrime >> _z0
     >> _gamma >> _maxflavour >> _gammaZ >> _process;
}

ClassDescription<MEqq2ZPrime2ff> MEqq2ZPrime2ff::initMEqq2ZPrime2ff;
// Definition of the static class description member.

void MEqq2ZPrime2ff::Init() {

  static ClassDocumentation<MEqq2ZPrime2ff> documentation
    ("The MEqq2ZPrime2ff class provides the matrix elements for "
     "q qbar -> f fbar in the RadiativeZPrimeModel including the option"
     " of interference with the photon and Z");

  static Parameter<MEqq2ZPrime2ff,unsigned int> interfaceMaxFlavour
    ("MaxFlavour",
     "The heaviest incoming quark flavour this matrix element is allowed to handle",
     &MEqq2ZPrime2ff::_maxflavour, 5, 1, 6,
     false, false, Interface::limited);

  static Switch<MEqq2ZPrime2ff,unsigned int> interfaceGammaZ
    ("GammaZ",
     "Which terms to include",
     &MEqq2ZPrime2ff::_gammaZ, 3, false, false);
  static SwitchOption interfaceGammaZAll
    (interfaceGammaZ,
     "All",
     "Include both gamma, Z and Z' terms",
     0);
  static SwitchOption interfaceGammaZGamma
    (interfaceGammaZ,
     "Gamma",
     "Only include the photon",
     1);
  static SwitchOption interfaceGammaZZ
    (interfaceGammaZ,
     "Z",
     "Only include the Z",
     2);
  static SwitchOption interfaceGammaZZPrime
    (interfaceGammaZ,
     "ZPrime",
     "Only include the Z'",
     3);

  static Switch<MEqq2ZPrime2ff,unsigned int> interfaceProcess
    ("Process",
     "Which process to included",
     &MEqq2ZPrime2ff::_process, 0, false, false);
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
  static SwitchOption interfaceProcessChargedLeptons
    (interfaceProcess,
     "ChargedLeptons",
     "Only include the charged leptons as outgoing particles",
     3);
  static SwitchOption interfaceProcessNeutrinos
    (interfaceProcess,
     "Neutrinos",
     "Only include the neutrinos as outgoing particles",
     4);
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
  static SwitchOption interfaceProcessNu_e
    (interfaceProcess,
     "Nu_e",
     "Only include nu_e ne_ebar as outgoing particles",
     8);
  static SwitchOption interfaceProcessnu_mu
    (interfaceProcess,
     "Nu_mu",
     "Only include nu_mu nu_mubar as outgoing particles",
     9);
  static SwitchOption interfaceProcessnu_tau
    (interfaceProcess,
     "Nu_tau",
     "Only include nu_tau nu_taubar as outgoing particles",
     10);
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

double MEqq2ZPrime2ff::qqME(vector<SpinorWaveFunction>    & fin ,
			    vector<SpinorBarWaveFunction> & ain ,
			    vector<SpinorBarWaveFunction> & fout,
			    vector<SpinorWaveFunction>    & aout,
			    bool calc) const {
  // scale
  Energy2 mb2(scale());
  // matrix element to be stored
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  // which intermediates to include
  bool gamma  = _gammaZ==0 || _gammaZ==1;
  bool Z0     = _gammaZ==0 || _gammaZ==2;
  bool Zprime = _gammaZ==0 || _gammaZ==3;
  // declare the variables we need
  unsigned int ihel1,ihel2,ohel1,ohel2;
  VectorWaveFunction inter[3];
  double me[4]={0.,0.,0.,0.};
  Complex diag[3];
  // sum over helicities to get the matrix element
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      // intermediate for Z
      if(Z0)     inter[0]=_theFFZVertex     ->
	evaluate(mb2,1,_z0    ,fin[ihel1],ain[ihel2]);
      // intermediate for photon
      if(gamma)  inter[1]=_theFFPVertex     ->
	evaluate(mb2,1,_gamma ,fin[ihel1],ain[ihel2]);
      // intermediate for Z'
      if(Zprime) inter[2]=_theFFZPrimeVertex-> 
	evaluate(mb2,1,_zPrime,fin[ihel1],ain[ihel2]);
      for(ohel1=0;ohel1<2;++ohel1) {
	for(ohel2=0;ohel2<2;++ohel2) {
	  // first the Z exchange diagram
	  diag[0] = Z0 ?
	    _theFFZVertex->evaluate(mb2,aout[ohel2],fout[ohel1],inter[0]) : 0.;
	  // first the photon exchange diagram
	  diag[1] = gamma ?
	    _theFFPVertex->evaluate(mb2,aout[ohel2],fout[ohel1],inter[1]) : 0.;
	  // then the Z' diagram
	  diag[2] = Zprime ?
	    _theFFZPrimeVertex->evaluate(mb2,aout[ohel2],fout[ohel1],inter[2]) : 0.;
	  // add up squares of individual terms
	  for(unsigned int ix=0;ix<3;++ix) me[ix+1] = norm(diag[ix]);
	  // the full thing including interference
	  diag[0]+=diag[1]+diag[2];
	  me[0] += norm(diag[0]);
	  if(calc) menew(ihel1,ihel2,ohel1,ohel2) = diag[0];
	}
      }
    }
  }
  // spin and colour factor
  double colspin=1./12.;
  if(abs(fout[0].id())<=6) colspin*=3.;
  // results
  // factor 12 from 4 helicity and 3 colour
  for(int ix=0;ix<3;++ix) me[ix]*=colspin;
  DVector save;
  save.push_back(me[1]);
  save.push_back(me[2]);
  save.push_back(me[3]);
  meInfo(save);
  if(calc) _me.reset(menew);
  return me[0];
}

void MEqq2ZPrime2ff::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // order of particles
  unsigned int order[4]={0,1,2,3};
  if(hard[0]->id()<0) swap(order[0],order[1]);
  if(hard[2]->id()<0) swap(order[2],order[3]);
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction(   fin ,hard[order[0]],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[order[1]],incoming,false,true);
  SpinorBarWaveFunction(fout,hard[order[2]],outgoing,true ,true);
  SpinorWaveFunction(   aout,hard[order[3]],outgoing,true ,true);
  qqME(fin,ain,fout,aout,true);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) 
    hard[order[ix]]->spinInfo()->productionVertex(hardvertex);
}
