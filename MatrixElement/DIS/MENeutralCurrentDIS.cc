// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MENeutralCurrentDIS class.
//

#include "MENeutralCurrentDIS.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/MatrixElement/General/HardVertex.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Handlers/StandardXComb.h"

using namespace Herwig;

void MENeutralCurrentDIS::doinit() throw(InitException) {
  ME2to2Base::doinit();
  _z0=getParticleData(ThePEG::ParticleID::Z0);
  _gamma=getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _theFFZVertex = hwsm->vertexFFZ();
    _theFFPVertex = hwsm->vertexFFP();
  }
  else
    throw InitException() << "Must be the Herwig++ StandardModel class in "
			  << "MENeutralCurrentDIS::doinit" << Exception::abortnow;
}

void MENeutralCurrentDIS::getDiagrams() const {
  // which intermediates to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // create the diagrams
  for(unsigned int ix=11;ix<=14;++ix) {
    for(unsigned int iz=0;iz<2;++iz) {
      tPDPtr lep = getParticleData(ix);
      if(iz==1) lep = lep->CC();
      for(unsigned int iy=_minflavour;iy<=_maxflavour;++iy) {
	tPDPtr quark = getParticleData(iy);
	// lepton quark scattering via gamma and Z
	if(gamma) add(new_ptr((Tree2toNDiagram(3), lep, _gamma, quark, 1, lep, 2, quark, -1)));
	if(Z0)    add(new_ptr((Tree2toNDiagram(3), lep, _z0   , quark, 1, lep, 2, quark, -2)));
	// lepton antiquark scattering via gamma and Z
	quark = quark->CC();
	if(gamma) add(new_ptr((Tree2toNDiagram(3), lep, _gamma, quark, 1, lep, 2, quark, -3)));
	if(Z0)    add(new_ptr((Tree2toNDiagram(3), lep, _z0   , quark, 1, lep, 2, quark, -4)));
      }
    }
  }
}

Energy2 MENeutralCurrentDIS::scale() const {
  return sHat();
}

bool MENeutralCurrentDIS::generateKinematics(const double * r) {
  double ctmin = -1.0;
  double ctmax =  1.0;
  meMomenta()[3].setMass(0.*GeV);

  Energy q = 0.0*GeV;
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), meMomenta()[2].mass(), meMomenta()[3].mass());
  } catch ( ImpossibleKinematics ) {
    return false;
  }

  Energy e = sqrt(sHat())/2.0;
		    
  Energy2 m22 = meMomenta()[2].mass2();
  Energy2 m32 = meMomenta()[3].mass2();
  Energy2 e0e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e1e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e0e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 e1e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 pq = 2.0*e*q;

  Energy2 thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[2]);
  if ( thmin > 0.0*GeV2 ) ctmax = min(ctmax, (e0e2 - m22 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[2]);
  if ( thmin > 0.0*GeV2 ) ctmin = max(ctmin, (thmin + m22 - e1e2)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[3]);
  if ( thmin > 0.0*GeV2 ) ctmax = min(ctmax, (e1e3 - m32 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[3]);
  if ( thmin > 0.0*GeV2 ) ctmin = max(ctmin, (thmin + m32 - e0e3)/pq);

  Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
   		     lastCuts().minKT(mePartonData()[3]));
  if ( ptmin > 0.0*GeV ) {
    double ctm = 1.0 - sqr(ptmin/q);
    if ( ctm <= 0.0 ) return false;
    ctmin = max(ctmin, -sqrt(ctm));
    ctmax = min(ctmax, sqrt(ctm));
  }

  if ( ctmin >= ctmax ) return false;
    
  double cth = getCosTheta(ctmin, ctmax, r);

  Energy pt = q*sqrt(1.0-sqr(cth));
  phi(rnd(2.0*Constants::pi));
  meMomenta()[2].setVect(Momentum3(pt*sin(phi()), pt*cos(phi()), q*cth));

  meMomenta()[3].setVect(Momentum3(-pt*sin(phi()),-pt*cos(phi()), -q*cth));

  meMomenta()[2].rescaleEnergy();
  meMomenta()[3].rescaleEnergy();

  vector<LorentzMomentum> out(2);
  out[0] = meMomenta()[2];
  out[1] = meMomenta()[3];
  tcPDVector tout(2);
  tout[0] = mePartonData()[2];
  tout[1] = mePartonData()[3];
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
    return false;

  tHat(pq*cth + m22 - e0e2);
  uHat(m22 + m32 - sHat() - tHat());
  jacobian((pq/sHat())*Constants::pi*jacobian());
  return true;
}

unsigned int MENeutralCurrentDIS::orderInAlphaS() const {
  return 1;
}

unsigned int MENeutralCurrentDIS::orderInAlphaEW() const {
  return 1;
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
     << _gamma << _z0;
}

void MENeutralCurrentDIS::persistentInput(PersistentIStream & is, int) {
  is >> _minflavour >> _maxflavour >> _gammaZ >> _theFFZVertex >> _theFFPVertex 
     >> _gamma >> _z0;
}

ClassDescription<MENeutralCurrentDIS> MENeutralCurrentDIS::initMENeutralCurrentDIS;
// Definition of the static class description member.

void MENeutralCurrentDIS::Init() {

  static ClassDocumentation<MENeutralCurrentDIS> documentation
    ("The MENeutralCurrentDIS class implements the matrix elements for leading-order "
     "neutral current deep inelastic scattering.");

  static Parameter<MENeutralCurrentDIS,unsigned int> interfaceMaxFlavour
    ("MaxFlavour",
     "The highest incoming quark flavour this matrix element is allowed to handle",
     &MENeutralCurrentDIS::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);

  static Parameter<MENeutralCurrentDIS,unsigned int> interfaceMinFlavour
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
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  // which intermediates to include
  bool gamma=_gammaZ==0||_gammaZ==1;
  bool Z0   =_gammaZ==0||_gammaZ==2;
  // declare the variables we need
  unsigned int ihel1,ihel2,ohel1,ohel2;
  VectorWaveFunction inter[2],test;
  double me[3]={0.,0.,0.};
  Complex diag1,diag2;
  // sum over helicities to get the matrix element
  unsigned int hel[4];
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      // intermediate for photon
      if(gamma) inter[0]=_theFFPVertex->evaluate(mb2,1,_gamma,f1[ihel1],a1[ihel2]);
      // intermediate for Z
      if(Z0)    inter[1]=_theFFZVertex->evaluate(mb2,1,_z0   ,f1[ihel1],a1[ihel2]);
      for(ohel1=0;ohel1<2;++ohel1) {
	for(ohel2=0;ohel2<2;++ohel2) {
	  hel[0] = ihel1;
	  hel[1] = ohel1;
	  hel[2] = ihel2;
	  hel[3] = ohel2;
	  if(!lorder) swap(hel[0],hel[2]);
	  if(!qorder) swap(hel[1],hel[3]);
	  // first the photon exchange diagram
	  diag1 = gamma ?
	    _theFFPVertex->evaluate(mb2,f2[ohel1],a2[ohel2],inter[0]) : 0.;
	  // then the Z exchange diagram
	  diag2 = Z0 ?
	    _theFFZVertex->evaluate(mb2,f2[ohel1],a2[ohel2],inter[1]) : 0.;
	  // add up squares of individual terms
	  me[1] += norm(diag1);
	  me[2] += norm(diag2);
	  // the full thing including interference
	  diag1 +=diag2;
	  me[0] += norm(diag1);
	  if(calc) menew(hel[0],hel[1],hel[2],hel[3]) = diag1;
	}
      }
    }
  }
  // spin and colour factor
  double colspin=0.25;
  // results
  for(int ix=0;ix<3;++ix) me[ix]*=colspin;
  DVector save;
  save.push_back(me[1]);
  save.push_back(me[2]);
  meInfo(save);
  if(calc) _me.reset(menew);
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
  unsigned int order[4]={0,1,2,3};
  bool lorder(true),qorder(true);
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
  SpinorWaveFunction   (f1,hard[order[0]],incoming,!lorder,true);
  SpinorWaveFunction   (f2,hard[order[1]],incoming,!qorder,true);
  SpinorBarWaveFunction(a1,hard[order[2]],outgoing, lorder,true);
  SpinorBarWaveFunction(a2,hard[order[3]],outgoing, qorder,true);
  helicityME(f1,f2,a1,a2,lorder,qorder,false);
  // get the spin info objects
  SpinfoPtr spin[4];
  for(unsigned int ix=0;ix<4;++ix) {
    spin[ix]=dynamic_ptr_cast<SpinfoPtr>(hard[ix]->spinInfo());
  }
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) {
    spin[ix]->setProductionVertex(hardvertex);
  }
}



