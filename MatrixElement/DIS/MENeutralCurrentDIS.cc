// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MENeutralCurrentDIS class.
//

#include "MENeutralCurrentDIS.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"

using namespace Herwig;

MENeutralCurrentDIS::MENeutralCurrentDIS() 
  : _minflavour(1), _maxflavour(5), _gammaZ(0),
    _sinW(0.), _cosW(0.), _mz2(ZERO) {
  vector<unsigned int> mopt(2,0);
  mopt[0] = 1;
  massOption(mopt);
}

void MENeutralCurrentDIS::doinit() {
  DISBase::doinit();
  _z0    = getParticleData(ThePEG::ParticleID::Z0);
  _gamma = getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() 
    << "Must be the Herwig StandardModel class in "
    << "MENeutralCurrentDIS::doinit" << Exception::abortnow;
  // vertices
  _theFFZVertex = hwsm->vertexFFZ();
  _theFFPVertex = hwsm->vertexFFP();
  // electroweak parameters
  _sinW = generator()->standardModel()->sin2ThetaW();
  _cosW = sqrt(1.-_sinW);
  _sinW = sqrt(_sinW);
  _mz2 = sqr(_z0->mass());
}

void MENeutralCurrentDIS::getDiagrams() const {
  // which intermediates to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // create the diagrams
  for(int ix=11;ix<=14;++ix) {
    for(unsigned int iz=0;iz<2;++iz) {
      tPDPtr lep = getParticleData(ix);
      if(iz==1) lep = lep->CC();
      for(int iy=_minflavour;iy<=_maxflavour;++iy) {
	tPDPtr quark = getParticleData(iy);
	// lepton quark scattering via gamma and Z
	if(gamma) add(new_ptr((Tree2toNDiagram(3), lep, _gamma, quark,
			       1, lep, 2, quark, -1)));
	if(Z0)    add(new_ptr((Tree2toNDiagram(3), lep, _z0   , quark,
			       1, lep, 2, quark, -2)));
	// lepton antiquark scattering via gamma and Z
	quark = quark->CC();
	if(gamma) add(new_ptr((Tree2toNDiagram(3), lep, _gamma, quark,
			       1, lep, 2, quark, -3)));
	if(Z0)    add(new_ptr((Tree2toNDiagram(3), lep, _z0   , quark,
			       1, lep, 2, quark, -4)));
      }
    }
  }
}

unsigned int MENeutralCurrentDIS::orderInAlphaS() const {
  return 0;
}

unsigned int MENeutralCurrentDIS::orderInAlphaEW() const {
  return 2;
}

Selector<const ColourLines *>
MENeutralCurrentDIS::colourGeometries(tcDiagPtr diag) const {
  static ColourLines c1("3 5");
  static ColourLines c2("-3 -5");
  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 || diag->id() == -2 )
    sel.insert(1.0, &c1);
  else
    sel.insert(1.0, &c2);
  return sel;
}

void MENeutralCurrentDIS::persistentOutput(PersistentOStream & os) const {
  os << _minflavour << _maxflavour << _gammaZ << _theFFZVertex << _theFFPVertex 
     << _gamma << _z0 << _sinW << _cosW << ounit(_mz2,GeV2);
}

void MENeutralCurrentDIS::persistentInput(PersistentIStream & is, int) {
  is >> _minflavour >> _maxflavour >> _gammaZ >> _theFFZVertex >> _theFFPVertex 
     >> _gamma >> _z0 >> _sinW >> _cosW >> iunit(_mz2,GeV2) ;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MENeutralCurrentDIS,DISBase>
describeHerwigMENeutralCurrentDIS("Herwig::MENeutralCurrentDIS", "HwMEDIS.so");

void MENeutralCurrentDIS::Init() {

  static ClassDocumentation<MENeutralCurrentDIS> documentation
    ("The MENeutralCurrentDIS class implements the matrix elements for leading-order "
     "neutral current deep inelastic scattering.");

  static Parameter<MENeutralCurrentDIS,int> interfaceMaxFlavour
    ("MaxFlavour",
     "The highest incoming quark flavour this matrix element is allowed to handle",
     &MENeutralCurrentDIS::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);

  static Parameter<MENeutralCurrentDIS,int> interfaceMinFlavour
    ("MinFlavour",
     "The lightest incoming quark flavour this matrix element is allowed to handle",
     &MENeutralCurrentDIS::_minflavour, 1, 1, 5,
     false, false, Interface::limited);

  static Switch<MENeutralCurrentDIS,unsigned int> interfaceGammaZ
    ("GammaZ",
     "Which terms to include",
     &MENeutralCurrentDIS::_gammaZ, 0, false, false);
  static SwitchOption interfaceGammaZAll
    (interfaceGammaZ,
     "All",
     "Include both gamma and Z terms",
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
}

Selector<MEBase::DiagramIndex>
MENeutralCurrentDIS::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if      ( diags[i]->id() == -1 || diags[i]->id() == -3 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 || diags[i]->id() == -4 ) sel.insert(meInfo()[1], i);
  }
  return sel;
}

double MENeutralCurrentDIS::helicityME(vector<SpinorWaveFunction>    & f1,
				       vector<SpinorWaveFunction>    & f2,
				       vector<SpinorBarWaveFunction> & a1,
				       vector<SpinorBarWaveFunction> & a2,
				       bool lorder, bool qorder,
				       bool calc) const {
  // scale
  Energy2 mb2(scale());
  // matrix element to be stored
  ProductionMatrixElement menew (PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  ProductionMatrixElement gamma (PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  ProductionMatrixElement Zboson(PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  // which intermediates to include
  bool gam = _gammaZ==0 || _gammaZ==1;
  bool Z0  = _gammaZ==0 || _gammaZ==2;
  // declare the variables we need
  VectorWaveFunction inter[2];
  double me[3]={0.,0.,0.};
  Complex diag1,diag2;
  // sum over helicities to get the matrix element
  unsigned int hel[4];
  unsigned int lhel1,lhel2,qhel1,qhel2;
  for(lhel1=0;lhel1<2;++lhel1) {
    for(lhel2=0;lhel2<2;++lhel2) {
      // intermediate for photon
      if(gam) inter[0]=_theFFPVertex->evaluate(mb2,1,_gamma,f1[lhel1],a1[lhel2]);
      // intermediate for Z
      if(Z0)    inter[1]=_theFFZVertex->evaluate(mb2,1,_z0   ,f1[lhel1],a1[lhel2]);
      for(qhel1=0;qhel1<2;++qhel1) {
	for(qhel2=0;qhel2<2;++qhel2) {
	  hel[0] = lhel1;
	  hel[1] = qhel1;
	  hel[2] = lhel2;
	  hel[3] = qhel2;
	  if(!lorder) swap(hel[0],hel[2]);
	  if(!qorder) swap(hel[1],hel[3]);
	  // first the photon exchange diagram
	  diag1 = gam ?
	    _theFFPVertex->evaluate(mb2,f2[qhel1],a2[qhel2],inter[0]) : 0.;
	  // then the Z exchange diagram
	  diag2 = Z0 ?
	    _theFFZVertex->evaluate(mb2,f2[qhel1],a2[qhel2],inter[1]) : 0.;
	  // add up squares of individual terms
	  me[1] += norm(diag1);
	  gamma (hel[0],hel[1],hel[2],hel[3]) = diag1;
	  me[2] += norm(diag2);
	  Zboson(hel[0],hel[1],hel[2],hel[3]) = diag2;
	  // the full thing including interference
	  diag1 += diag2;
	  me[0] += norm(diag1);
	  menew(hel[0],hel[1],hel[2],hel[3]) = diag1;
	}
      }
    }
  }
  // spin and colour factor
  double colspin = 0.25;
  // results
  for(int ix=0;ix<3;++ix) me[ix] *= colspin;
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
			 beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    me[0] = menew .average(rho[0],rho[1]);
    me[1] = gamma .average(rho[0],rho[1]);
    me[2] = Zboson.average(rho[0],rho[1]);
  }
  DVector save;
  save.push_back(me[1]);
  save.push_back(me[2]);
  meInfo(save);
  if(calc) _me.reset(menew);
  // analytic expression for testing
//   double test = 8.*sqr(4.*Constants::pi*generator()->standardModel()->alphaEM(mb2))*
//     sqr(double(mePartonData()[1]->iCharge())/3.)/sqr(tHat())
//     *(sqr(sHat())+sqr(uHat())+4.*sqr(mePartonData()[0]->mass())*tHat())/4.;
//   cerr << "testing me " << me[0]/test << "\n";
  return me[0];
}

double MENeutralCurrentDIS::me2() const {
  vector<SpinorWaveFunction>    f1,f2;
  vector<SpinorBarWaveFunction> a1,a2;
  bool lorder,qorder;
  SpinorWaveFunction    l1,q1;
  SpinorBarWaveFunction l2,q2;
  // lepton wave functions
  if(mePartonData()[0]->id()>0) {
    lorder=true;
    l1 = SpinorWaveFunction   (meMomenta()[0],mePartonData()[0],incoming);
    l2 = SpinorBarWaveFunction(meMomenta()[2],mePartonData()[2],outgoing);
  }
  else {
    lorder=false;
    l1 = SpinorWaveFunction   (meMomenta()[2],mePartonData()[2],outgoing);
    l2 = SpinorBarWaveFunction(meMomenta()[0],mePartonData()[0],incoming);
  }
  // quark wave functions
  if(mePartonData()[1]->id()>0) {
    qorder = true;
    q1 = SpinorWaveFunction   (meMomenta()[1],mePartonData()[1],incoming);
    q2 = SpinorBarWaveFunction(meMomenta()[3],mePartonData()[3],outgoing);
  }
  else {
    qorder = false;
    q1 = SpinorWaveFunction   (meMomenta()[3],mePartonData()[3],outgoing);
    q2 = SpinorBarWaveFunction(meMomenta()[1],mePartonData()[1],incoming);
  }
  // wavefunctions for various helicities
  for(unsigned int ix=0;ix<2;++ix) {
    l1.reset(ix); f1.push_back(l1);
    l2.reset(ix); a1.push_back(l2);
    q1.reset(ix); f2.push_back(q1);
    q2.reset(ix); a2.push_back(q2);
  }
  return helicityME(f1,f2,a1,a2,lorder,qorder,false);
}

void MENeutralCurrentDIS::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // sort out the ordering
  unsigned int order[4]={0,1,2,3};
  bool lorder(true),qorder(true);
  if(abs(hard[0]->id())<6) swap(hard[0],hard[1]);
  if(abs(hard[2]->id())<6) swap(hard[2],hard[3]);
  if(hard[0]->id()<0) {
    swap(order[0],order[2]);
    lorder = false;
  }
  if(hard[1]->id()<0) {
    swap(order[1],order[3]);
    qorder = false;
  }
  vector<SpinorWaveFunction>    f1,f2;
  vector<SpinorBarWaveFunction> a1,a2;
  SpinorWaveFunction   (f1,hard[order[0]], lorder ? incoming : outgoing, !lorder,true);
  SpinorWaveFunction   (f2,hard[order[1]], qorder ? incoming : outgoing, !qorder,true);
  SpinorBarWaveFunction(a1,hard[order[2]], lorder ? outgoing : incoming,  lorder,true);
  SpinorBarWaveFunction(a2,hard[order[3]], qorder ? outgoing : incoming,  qorder,true);
  helicityME(f1,f2,a1,a2,lorder,qorder,true);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) {
    tSpinPtr spin = hard[ix]->spinInfo();
    if(ix<2) {
      tcPolarizedBeamPDPtr beam = 
	dynamic_ptr_cast<tcPolarizedBeamPDPtr>(hard[ix]->dataPtr());
      if(beam) spin->rhoMatrix() = beam->rhoMatrix();
    }
    spin->productionVertex(hardvertex);
  }
}

double MENeutralCurrentDIS::A(tcPDPtr lin, tcPDPtr,
			      tcPDPtr qin, tcPDPtr, Energy2 q2) const {
  // photon only 
  if(_gammaZ==1) return 0.;
  // everything else
  double r = _gammaZ==0 || _gammaZ==2 ? double(q2/(q2+_mz2)) : 0.;
  double eL,eQ,cvL,caL,cvQ,caQ;;
  if(abs(lin->id())%2==0) {
    eL = _gammaZ==0 ? generator()->standardModel()->enu() : 0.;
    cvL = 0.25*generator()->standardModel()->vnu();
    caL = 0.25*generator()->standardModel()->anu();
  }
  else {
    eL = _gammaZ==0 ? generator()->standardModel()->ee()  : 0.;
    cvL = 0.25*generator()->standardModel()->ve();
    caL = 0.25*generator()->standardModel()->ae();
  }
  if(abs(qin->id())%2==0) {
    eQ = _gammaZ==0 ? generator()->standardModel()->eu() : 0.;
    cvQ = 0.25*generator()->standardModel()->vu();
    caQ = 0.25*generator()->standardModel()->au();
  }
  else {
    eQ = _gammaZ==0 ? generator()->standardModel()->ed() : 0.;
    cvQ = 0.25*generator()->standardModel()->vd();
    caQ = 0.25*generator()->standardModel()->ad();
  }
  double output = 4.*r*caL*caQ*(eL*eQ+2.*r*cvL*cvQ/sqr(_sinW*_cosW))
    /sqr(_sinW*_cosW)/(sqr(eL*eQ)+2.*eL*eQ*r/sqr(_cosW*_sinW)*cvL*cvQ
      +sqr(r/sqr(_cosW*_sinW))*(sqr(cvL)+sqr(caL))*(sqr(cvQ)+sqr(caQ)));
  if(qin->id()<0) output *= -1.;
  if(lin->id()<0) output *= -1.;
  return output;
}
