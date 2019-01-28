// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2ZJet class.
//

#include "MEPP2ZJet.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

MEPP2ZJet::MEPP2ZJet() : _process(0), _maxflavour(5), _zdec(0),
			 _gammaZ(0), _widthopt(1), _pprob(0.5)
{}

void MEPP2ZJet::doinit() {
  HwMEBase::doinit();
  _z0    = getParticleData(ThePEG::ParticleID::Z0   );
  _gamma = getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  ThePEG::Ptr<Herwig::StandardModel>::transient_const_pointer 
    hwsm=ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardModel>
    ::transient_const_pointer>(standardModel());
  // do the initialisation
  if(!hwsm) 
    throw InitException() << "Must be Herwig::StandardModel in MEPP2ZJet::doinit()"
			  << Exception::runerror;
  // set the vertex pointers
  _theFFZVertex = hwsm->vertexFFZ();
  _theFFPVertex = hwsm->vertexFFP();
  _theQQGVertex = hwsm->vertexFFG();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2ZJet,HwMEBase>
describeHerwigMEPP2ZJet("Herwig::MEPP2ZJet", "HwMEHadron.so");

void MEPP2ZJet::Init() {

  static ClassDocumentation<MEPP2ZJet> documentation
    ("The MEPP2ZJet class implements the matrix element for Z/gamma+ jet production");

  static Parameter<MEPP2ZJet,int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEPP2ZJet::_maxflavour, 5, 0, 8, false, false, true);

  static Switch<MEPP2ZJet,int> interfaceZDecay
    ("ZDecay",
     "Which process to included",
     &MEPP2ZJet::_zdec, 0, false, false);
  static SwitchOption interfaceZDecayAll
    (interfaceZDecay,
     "All",
     "Include all SM fermions as outgoing particles",
     0);
  static SwitchOption interfaceZDecayQuarks
    (interfaceZDecay,
     "Quarks",
     "All include the quarks as outgoing particles",
     1);
  static SwitchOption interfaceZDecayLeptons
    (interfaceZDecay,
     "Leptons",
     "Only include the leptons as outgoing particles",
     2);
  static SwitchOption interfaceZDecayChargedLeptons
    (interfaceZDecay,
     "ChargedLeptons",
     "Only include the charged leptons as outgoing particles",
     3);
  static SwitchOption interfaceZDecayNeutrinos
    (interfaceZDecay,
     "Neutrinos",
     "Only include the neutrinos as outgoing particles",
     4);
  static SwitchOption interfaceZDecayElectron
    (interfaceZDecay,
     "Electron",
     "Only include e+e- as outgoing particles",
     5);
  static SwitchOption interfaceZDecayMuon
    (interfaceZDecay,
     "Muon",
     "Only include mu+mu- as outgoing particles",
     6);
  static SwitchOption interfaceZDecayTau
    (interfaceZDecay,
     "Tau",
     "Only include tau+tau- as outgoing particles",
     7);
  static SwitchOption interfaceZDecayNu_e
    (interfaceZDecay,
     "Nu_e",
     "Only include nu_e ne_ebar as outgoing particles",
     8);
  static SwitchOption interfaceZDecaynu_mu
    (interfaceZDecay,
     "Nu_mu",
     "Only include nu_mu nu_mubar as outgoing particles",
     9);
  static SwitchOption interfaceZDecaynu_tau
    (interfaceZDecay,
     "Nu_tau",
     "Only include nu_tau nu_taubar as outgoing particles",
     10);
  static SwitchOption interfaceZDecayDown
    (interfaceZDecay,
     "Down",
     "Only include d dbar as outgoing particles",
     11);
  static SwitchOption interfaceZDecayUp
    (interfaceZDecay,
     "Up",
     "Only include u ubar as outgoing particles",
     12);
  static SwitchOption interfaceZDecayStrange
    (interfaceZDecay,
     "Strange",
     "Only include s sbar as outgoing particles",
     13);
  static SwitchOption interfaceZDecayCharm
    (interfaceZDecay,
     "Charm",
     "Only include c cbar as outgoing particles",
     14);
  static SwitchOption interfaceZDecayBottom
    (interfaceZDecay,
     "Bottom",
     "Only include b bbar as outgoing particles",
     15);
  static SwitchOption interfaceZDecayTop
    (interfaceZDecay,
     "Top",
     "Only include t tbar as outgoing particles",
     16);

  static Switch<MEPP2ZJet,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2ZJet::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcessqqbar
    (interfaceProcess,
     "qqbar",
     "Only include q qbar -> Z/gamma g process",
     1);
  static SwitchOption interfaceProcessqg
    (interfaceProcess,
     "qg",
     "Only include the q g -> Z/gamma q process",
     2);
  static SwitchOption interfaceProcessqbarg
    (interfaceProcess,
     "qbarg",
     "Only include the qbar g -> Z/gamma qbar process",
     3);

  static Parameter<MEPP2ZJet,double> interfacePhotonProbablity
    ("PhotonProbablity",
     "Probability for using the \\f$1/s^2\\f$ piece for the"
     " generation of the gauge boson mass",
     &MEPP2ZJet::_pprob, 0.5, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<MEPP2ZJet,unsigned int> interfaceGammaZ
    ("GammaZ",
     "Which terms to include",
     &MEPP2ZJet::_gammaZ, 0, false, false);
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

  static Switch<MEPP2ZJet,unsigned int> interfaceWidthOption
    ("WidthOption",
     "The option for handling the width of the off-shell W boson",
     &MEPP2ZJet::_widthopt, 1, false, false);
  static SwitchOption interfaceWidthOptionFixedDenominator
    (interfaceWidthOption,
     "FixedDenominator",
     "Use a fxied with in the W propagator but the full matrix element"
     " in the numerator",
     1);
  static SwitchOption interfaceWidthOptionAllRunning
    (interfaceWidthOption,
     "AllRunning",
     "Use a running width in the W propagator and the full matrix "
     "element in the numerator",
     2);

}

void MEPP2ZJet::getDiagrams() const {
  // which intermediates to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // pointer for gluon
  tcPDPtr g = getParticleData(ParticleID::g);
  bool quark,lepton;
  for ( int ix=1; ix<17; ++ix ) {
    // increment counter to switch between quarks and leptons
    if(ix==7) ix+=4;
    // is it a valid quark process
    quark=ix<=6&&(_zdec==0||_zdec==1||_zdec-10==ix);
    // is it a valid lepton process
    lepton=ix>=11&&ix<=16&&
      (_zdec==0||_zdec==2||(_zdec==3&&ix%2==1)||
       (_zdec==4&&ix%2==0)||(ix%2==0&&(ix-10)/2==_zdec-7)||
       (ix%2==1&&(ix-9)/2==_zdec-4));
    // if not a validf process continue
    if(!(quark||lepton)) continue;
    // pointer for Z decay products
    tcPDPtr lm = getParticleData(ix);
    tcPDPtr lp = lm->CC();
    for (int i = 1; i <= _maxflavour; ++i ) {
      tcPDPtr q = getParticleData(i);
      tcPDPtr qb = q->CC();
      // q qbar -> Z g -> l+l- g
      if(_process==0||_process==1) {
	if(gamma) add(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, _gamma,
			       2, g,  4, lm, 4, lp, -1)));
	if(Z0)    add(new_ptr((Tree2toNDiagram(3), q, q, qb, 1,    _z0,
			       2, g,  4, lm, 4, lp, -2)));
	if(gamma) add(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, _gamma,
			       1, g,  4, lm, 4, lp, -3)));
	if(Z0)    add(new_ptr((Tree2toNDiagram(3), q, q, qb, 2,    _z0,
			       1, g,  4, lm, 4, lp, -4)));
      }
      // q g   -> Z q -> l+l- qbar
      if(_process==0||_process==2) {
	if(gamma) add(new_ptr((Tree2toNDiagram(3), q, q, g,    1, _gamma,
			       2, q,  4, lm, 4, lp, -5)));
	if(Z0)    add(new_ptr((Tree2toNDiagram(3), q, q, g,    1,    _z0,
			       2, q,  4, lm, 4, lp, -6)));
	if(gamma) add(new_ptr((Tree2toNDiagram(2), q, g, 1, q, 3, _gamma,
			       3, q,  4, lm, 4, lp, -7)));
	if(Z0)    add(new_ptr((Tree2toNDiagram(2), q, g, 1, q, 3,    _z0,
			       3, q,  4, lm, 4, lp, -8)));
      }
      // qbar g   -> Z qbar -> l+l- qbar
      if(_process==0||_process==3) {
	if(gamma) add(new_ptr((Tree2toNDiagram(3), qb, qb, g,     1, _gamma,
			       2, qb,  4, lm, 4, lp, -9 )));
	if(Z0)    add(new_ptr((Tree2toNDiagram(3), qb, qb, g,     1,    _z0,
			       2, qb,  4, lm, 4, lp, -10)));
	if(gamma) add(new_ptr((Tree2toNDiagram(2), qb,  g, 1, qb, 3, _gamma,
			       3, qb,  4, lm, 4, lp, -11)));
	if(Z0)    add(new_ptr((Tree2toNDiagram(2), qb,  g, 1, qb, 3,    _z0,
			       3, qb,  4, lm, 4, lp, -12)));
      }
    }
  }
}

unsigned int MEPP2ZJet::orderInAlphaS() const {
  return 1;
}

unsigned int MEPP2ZJet::orderInAlphaEW() const {
  return 2;
}

void MEPP2ZJet::persistentOutput(PersistentOStream & os) const {
  os << _theFFZVertex << _theFFPVertex << _theQQGVertex << _z0 << _widthopt
     << _gamma << _process << _maxflavour << _zdec << _pprob << _gammaZ;
}

void MEPP2ZJet::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZVertex >> _theFFPVertex >> _theQQGVertex >> _z0 >> _widthopt
     >> _gamma >> _process >> _maxflavour >> _zdec >> _pprob >> _gammaZ;
}

int MEPP2ZJet::nDim() const {
  return 5;
}

Selector<const ColourLines *>
MEPP2ZJet::colourGeometries(tcDiagPtr diag) const {
  // colour lines for q qbar -> Z g
  static const ColourLines cqqbar[4]={ColourLines("1 2 5,-3 -5"),
				      ColourLines("1 5,-5 2 -3"),
				      ColourLines("1 2 5,-3 -5, 6 -7"),
				      ColourLines("1 5,-5 2 -3, 6 -7")};
  // colour lines for q g -> Z q
  static const ColourLines cqg   [4]={ColourLines("1 2 -3,3 5"),
				      ColourLines("1 -2,2 3 5"),
				      ColourLines("1 2 -3,3 5, 6 -7"),
				      ColourLines("1 -2,2 3 5, 6 -7")};
  // colour lines for qbar q -> Z qbar
  static const ColourLines cqbarg[4]={ColourLines("-1 -2 3,-3 -5"),
				      ColourLines("-1 2,-2 -3 -5"),
				      ColourLines("-1 -2 3,-3 -5, 6 -7"),
				      ColourLines("-1 2,-2 -3 -5, 6 -7")};
  // select the correct line
  unsigned int icol = mePartonData()[3]->coloured() ? 2 : 0;
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
  case  1 : case  2:
    sel.insert(1.0, &cqqbar[icol]);
    break;
  case  3 : case  4:
    sel.insert(1.0, &cqqbar[icol+1]);
    break;
  case  5 : case  6:
    sel.insert(1.0, &cqg[icol]);
    break;
  case  7 : case  8:
    sel.insert(1.0, &cqg[icol+1]);
    break;
  case  9 : case 10:
    sel.insert(1.0, &cqbarg[icol]);
    break;
  case 11 : case 12:
    sel.insert(1.0, &cqbarg[icol+1]);
    break;
  }
  return sel;
}

Selector<MEBase::DiagramIndex>
MEPP2ZJet::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    int id=abs(diags[i]->id());
    if     (id <= 4 ) sel.insert(meInfo()[id-1],i);
    else if(id <= 8 ) sel.insert(meInfo()[id-5],i);
    else if(id <= 12) sel.insert(meInfo()[id-9],i);
  }
  return sel;
}

Energy2 MEPP2ZJet::scale() const {
  return _scale;
}

CrossSection MEPP2ZJet::dSigHatDR() const {
  return me2()*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
}

bool MEPP2ZJet::generateKinematics(const double * r) {
  // initialize jacobian
  jacobian(1.);
  // cms energy
  Energy ecm=sqrt(sHat());
  // first generate the mass of the off-shell gauge boson
  // minimum mass of the 
  tcPDVector ptemp;
  ptemp.push_back(mePartonData()[3]);
  ptemp.push_back(mePartonData()[4]);
  Energy2 minMass2 = max(lastCuts().minSij(mePartonData()[3],mePartonData()[4]),
			 lastCuts().minS(ptemp));
  // minimum pt of the jet
  Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
		     lastCuts().minKT(_z0));
  // maximum mass of the gauge boson so pt is possible
  Energy2 maxMass2 = min(ecm*(ecm-2.*ptmin),lastCuts().maxS(ptemp));
  if(maxMass2<=ZERO||minMass2<ZERO) return false;
  // also impose the limits from the ParticleData object
  minMass2 = max(minMass2,sqr(_z0->massMin()));
  maxMass2 = min(maxMass2,sqr(_z0->massMax()));
  // also impose the limits from the ParticleData object
  if(maxMass2<minMass2) return false;
  // generation of the mass
  Energy  M(_z0->mass()),Gamma(_z0->width());
  Energy2 M2(sqr(M)),MG(M*Gamma);
  double rhomin = atan2((minMass2-M2),MG);
  double rhomax = atan2((maxMass2-M2),MG);
  if(r[1]<_pprob) {
    double rand=r[1]/_pprob;
    _mz2=minMass2*maxMass2/(minMass2+rand*(maxMass2-minMass2));
  }
  else {
    double rand=(r[1]-_pprob)/(1.-_pprob);
    _mz2=M2+MG*tan(rhomin+rand*(rhomax-rhomin));
  }
  Energy mz=sqrt(_mz2);
  InvEnergy2 emjac1 = _pprob*minMass2*maxMass2/(maxMass2-minMass2)/sqr(_mz2);
  InvEnergy2 emjac2 = (1.-_pprob)*MG/(rhomax-rhomin)/(sqr(_mz2-M2)+sqr(MG));
  // jacobian
  jacobian(jacobian()/sHat()/(emjac1+emjac2));
  // set the masses of the outgoing particles to 2-2 scattering
  meMomenta()[2].setMass(ZERO);
  Lorentz5Momentum pz(mz);
  // generate the polar angle of the hard scattering
  double ctmin(-1.0), ctmax(1.0);
  Energy q(ZERO);
  try {
    q = SimplePhaseSpace::getMagnitude(sHat(), meMomenta()[2].mass(),mz);
  } 
  catch ( ImpossibleKinematics & e ) {
    return false;
  }	    
  Energy2 pq = sqrt(sHat())*q;
  if ( ptmin > ZERO ) {
    double ctm = 1.0 - sqr(ptmin/q);
    if ( ctm <= 0.0 ) return false;
    ctmin = max(ctmin, -sqrt(ctm));
    ctmax = min(ctmax,  sqrt(ctm));
  }
  if ( ctmin >= ctmax ) return false;
  double cth = getCosTheta(ctmin, ctmax, r[0]);
  Energy pt  = q*sqrt(1.0-sqr(cth));
  double phi = 2.0*Constants::pi*r[2];
  meMomenta()[2].setVect(Momentum3( pt*sin(phi), pt*cos(phi), q*cth));
  pz.setVect(            Momentum3(-pt*sin(phi),-pt*cos(phi),-q*cth));
  meMomenta()[2].rescaleEnergy();
  pz.rescaleEnergy();
  // set the scale
  _scale = _mz2+sqr(pt);
  // generate the momenta of the Z decay products
  meMomenta()[3].setMass(mePartonData()[3]->mass());
  meMomenta()[4].setMass(mePartonData()[4]->mass());
  Energy q2 = ZERO;
  try {
    q2 = SimplePhaseSpace::getMagnitude(_mz2, meMomenta()[3].mass(),
					meMomenta()[4].mass());
  } catch ( ImpossibleKinematics & e ) {
    return false;
  }
  double cth2 =-1.+2.*r[3];
  double phi2=Constants::twopi*r[4];
  Energy pt2 =q2*sqrt(1.-sqr(cth2));
  Lorentz5Momentum pl[2]={Lorentz5Momentum( pt2*cos(phi2), pt2*sin(phi2), q2*cth2,ZERO,
					    meMomenta()[3].mass()),
			  Lorentz5Momentum(-pt2*cos(phi2),-pt2*sin(phi2),-q2*cth2,ZERO,
					   meMomenta()[4].mass())};
  pl[0].rescaleEnergy();
  pl[1].rescaleEnergy();
  Boost boostv(pz.boostVector());
  pl[0].boost(boostv);
  pl[1].boost(boostv);
  meMomenta()[3] = pl[0];
  meMomenta()[4] = pl[1];
  // check passes all the cuts
  vector<LorentzMomentum> out(3);
  tcPDVector tout(3);
  for(unsigned int ix=0;ix<3;++ix) {
    out[ ix] = meMomenta()[ix+2];
    tout[ix] = mePartonData()[ix+2];
  }
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
    return false;
  // jacobian
  jacobian((pq/sHat())*Constants::pi*jacobian()/8./sqr(Constants::pi)*q2/mz);
  return true;
}

double MEPP2ZJet::me2() const {
  InvEnergy2 output(ZERO);
  // construct spinors for the leptons (always the same)
  vector<SpinorBarWaveFunction> lm;
  vector<SpinorWaveFunction>    lp;
  SpinorBarWaveFunction lmout(meMomenta()[3],mePartonData()[3],outgoing);
  SpinorWaveFunction    lpout(meMomenta()[4],mePartonData()[4],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    lmout.reset(ix);lm.push_back(lmout);
    lpout.reset(ix);lp.push_back(lpout);
  }
  // q g to q Z
  if(mePartonData()[0]->id()<=6&&mePartonData()[0]->id()>0&&
     mePartonData()[1]->id()==ParticleID::g) {
    // polarization states for the particles
    vector<SpinorWaveFunction> fin;
    vector<VectorWaveFunction> gin;
    vector<SpinorBarWaveFunction> fout;
    SpinorWaveFunction    qin (meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction    glin(meMomenta()[1],mePartonData()[1],incoming);
    SpinorBarWaveFunction qout(meMomenta()[2],mePartonData()[2],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix) ; fin.push_back(qin);
      glin.reset(2*ix); gin.push_back(glin);
      qout.reset(ix);fout.push_back(qout);
    }
    output=qgME(fin,gin,fout,lm,lp);
  }
  // qbar g to qbar Z
  else if(mePartonData()[0]->id()>=-6&&mePartonData()[0]->id()<0&&
	  mePartonData()[1]->id()==ParticleID::g) {
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gin;
    vector<SpinorWaveFunction> aout;
    SpinorBarWaveFunction qbin (meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction    glin (meMomenta()[1],mePartonData()[1],incoming);
    SpinorWaveFunction    qbout(meMomenta()[2],mePartonData()[2],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qbin .reset(ix  ); ain .push_back(qbin );
      glin .reset(2*ix); gin .push_back(glin );
      qbout.reset(ix  ); aout.push_back(qbout);
    }
    output=qbargME(ain,gin,aout,lm,lp);
  }
  // q qbar to g Z
  else {
    vector<SpinorWaveFunction>     fin;
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gout;
    SpinorWaveFunction    qin (meMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction qbin(meMomenta()[1],mePartonData()[1],incoming);
    VectorWaveFunction   glout(meMomenta()[2],mePartonData()[2],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)  ;  fin.push_back(qin);
      qbin.reset(ix) ;  ain.push_back(qbin);
      glout.reset(2*ix); gout.push_back(glout);
    }
    output=qqbarME(fin,ain,gout,lm,lp);
  }
  return output*sHat();
}

InvEnergy2 MEPP2ZJet::qqbarME(vector<SpinorWaveFunction> & fin,
			      vector<SpinorBarWaveFunction> & ain,
			      vector<VectorWaveFunction> & gout,
			      vector<SpinorBarWaveFunction> & lm,
			      vector<SpinorWaveFunction> & lp,
			      bool calc) const {
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
					     PDT::Spin1,PDT::Spin1Half,
					     PDT::Spin1Half));
  // diagrams to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // some integers
  unsigned int ihel1,ihel2,ohel1,ohel2,ohel3;
  // compute the leptonic photon and Z currents for speed
  VectorWaveFunction bcurr[2][2][2];
  for(ohel2=0;ohel2<2;++ohel2) {
    for(ohel3=0;ohel3<2;++ohel3) {
      // photon current
      if(gamma) bcurr[0][ohel2][ohel3]=
	_theFFPVertex->evaluate(_mz2,1,_gamma,lp[ohel3],lm[ohel2]);
      // Z current
      if(Z0)    bcurr[1][ohel2][ohel3]=
	_theFFZVertex->evaluate(_mz2,_widthopt,_z0,lp[ohel3],lm[ohel2]);
    }
  }
  // compute the matrix elements
  double me[5]={0.,0.,0.,0.,0.};
  Complex diag[4];
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      for(ohel1=0;ohel1<2;++ohel1) {
	// intermediates for the diagrams
 	inters=_theQQGVertex->evaluate(_scale,5,mePartonData()[0],
				       fin[ihel1],gout[ohel1]);
	interb=_theQQGVertex->evaluate(_scale,5,mePartonData()[1],
				       ain[ihel2],gout[ohel1]);
	for(ohel2=0;ohel2<2;++ohel2) {
	  for(ohel3=0;ohel3<2;++ohel3) {
	    diag[0] = gamma ? 
	      _theFFPVertex->evaluate(_mz2,fin[ihel1],interb,
				      bcurr[0][ohel2][ohel3]) : 0.;
	    diag[1]= Z0 ? 
	      _theFFZVertex->evaluate(_mz2,fin[ihel1],interb,
				      bcurr[1][ohel2][ohel3]) : 0.;
	    diag[2]= gamma ? 
	      _theFFPVertex->evaluate(_mz2,inters,ain[ihel2],
				      bcurr[0][ohel2][ohel3]) : 0.;
	    diag[3]= Z0 ? 
	      _theFFZVertex->evaluate(_mz2,inters,ain[ihel2],
				      bcurr[1][ohel2][ohel3]) : 0.;
	    // diagram contributions
	    me[1] += norm(diag[0]);
	    me[2] += norm(diag[1]);
	    me[3] += norm(diag[2]);
	    me[4] += norm(diag[3]);
	    // total
	    diag[0] += diag[1] + diag[2] + diag[3];
	    me[0]   += norm(diag[0]);
	    if(calc) _me(ihel1,ihel2,2*ohel1,ohel2,ohel3) = diag[0];
	  }
	}
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin = 1./9./4.;
  // and C_F N_c from matrix element
  colspin *= 4.;
  // and for Z decay products
  if(mePartonData()[3]->coloured()) colspin *= 3.;
  DVector save;
  for(unsigned int ix=0;ix<5;++ix) {
    me[ix] *= colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0]*UnitRemoval::InvE2;
}

InvEnergy2 MEPP2ZJet::qgME(vector<SpinorWaveFunction> & fin,
			   vector<VectorWaveFunction> & gin,
			   vector<SpinorBarWaveFunction> & fout,
			   vector<SpinorBarWaveFunction> & lm,
			   vector<SpinorWaveFunction> & lp,
			   bool calc) const {
  // diagrams to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
					     PDT::Spin1Half,PDT::Spin1Half,
					     PDT::Spin1Half));
  // some integers
  unsigned int ihel1,ihel2,ohel1,ohel2,ohel3;
  // compute the leptonic photon and Z currents for speed
  VectorWaveFunction bcurr[2][2][2];
  for(ohel2=0;ohel2<2;++ohel2) {
    for(ohel3=0;ohel3<2;++ohel3) {
      // photon current
      if(gamma) bcurr[0][ohel2][ohel3]=
	_theFFPVertex->evaluate(_mz2,1,_gamma,lp[ohel3],lm[ohel2]);
      // Z current
      if(Z0)    bcurr[1][ohel2][ohel3]=
	_theFFZVertex->evaluate(_mz2,_widthopt,_z0,lp[ohel3],lm[ohel2]);
    }
  }
  // compute the matrix elements
  double me[5]={0.,0.,0.,0.,0.};
  Complex diag[4];
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;
  Energy2 _scale=scale();
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      for(ohel1=0;ohel1<2;++ohel1) {
	// intermediates for the diagrams
	interb=_theQQGVertex->evaluate(_scale,5,mePartonData()[2]->CC(),
				       fout[ohel1],gin[ihel2]);
	inters=_theQQGVertex->evaluate(_scale,5,mePartonData()[0],
				       fin[ihel1],gin[ihel2]);
	for(ohel2=0;ohel2<2;++ohel2) {
	  for(ohel3=0;ohel3<2;++ohel3) {
	    diag[0]=gamma ?
	      _theFFPVertex->evaluate(_mz2,fin[ihel1],interb,
	      	               	      bcurr[0][ohel2][ohel3]) : 0.;
	    diag[1]=Z0    ?
	      _theFFZVertex->evaluate(_mz2,fin[ihel1],interb,
				      bcurr[1][ohel2][ohel3]) : 0.;
	    diag[2]=gamma ?
	      _theFFPVertex->evaluate(_mz2,inters,fout[ohel1],
	      	               	      bcurr[0][ohel2][ohel3]) : 0.;
	    diag[3]=Z0    ?
	      _theFFZVertex->evaluate(_mz2,inters,fout[ohel1],
				      bcurr[1][ohel2][ohel3]) : 0.;
	    // diagram contributions
	    me[1] += norm(diag[0]);
	    me[2] += norm(diag[1]);
	    me[3] += norm(diag[2]);
	    me[4] += norm(diag[3]);
	    // total
	    diag[0] += diag[1] + diag[2] + diag[3];
	    me[0]   += norm(diag[0]);
	    if(calc) _me(ihel1,2*ihel2,ohel1,ohel2,ohel3) = diag[0];
	  }
	}
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin = 1./24./4.;
  // and C_F N_c from matrix element
  colspin *= 4.;
  // and for Z decay products
  if(mePartonData()[3]->coloured()) colspin *= 3.;
  DVector save;
  for(unsigned int ix=0;ix<5;++ix) {
    me[ix] *= colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0]*UnitRemoval::InvE2;
}

InvEnergy2 MEPP2ZJet::qbargME(vector<SpinorBarWaveFunction> & fin,
			  vector<VectorWaveFunction> & gin,
			  vector<SpinorWaveFunction> & fout,
			  vector<SpinorBarWaveFunction> & lm,
			  vector<SpinorWaveFunction> & lp,
			  bool calc) const {
  // diagrams to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
					     PDT::Spin1Half,PDT::Spin1Half,
					     PDT::Spin1Half));
  // some integers
  unsigned int ihel1,ihel2,ohel1,ohel2,ohel3;
  // compute the leptonic photon and Z currents for speed
  VectorWaveFunction bcurr[2][2][2];
  for(ohel2=0;ohel2<2;++ohel2) {
    for(ohel3=0;ohel3<2;++ohel3) {
      // photon current
      if(gamma) bcurr[0][ohel2][ohel3]=
	_theFFPVertex->evaluate(_mz2,1,_gamma,lp[ohel3],lm[ohel2]);
      // Z current
      if(Z0)    bcurr[1][ohel2][ohel3]=
	_theFFZVertex->evaluate(_mz2,_widthopt,_z0,lp[ohel3],lm[ohel2]);
    }
  }
  // compute the matrix elements
  double me[5]={0.,0.,0.,0.,0.};
  Complex diag[4];
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;
  Energy2 _scale=scale();
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      for(ohel1=0;ohel1<2;++ohel1) {
	// intermediates for the diagrams
	inters=_theQQGVertex->evaluate(_scale,5,mePartonData()[2]->CC(),
				       fout[ohel1],gin[ihel2]);
	interb=_theQQGVertex->evaluate(_scale,5,mePartonData()[0],
				       fin[ihel1],gin[ihel2]);
	for(ohel2=0;ohel2<2;++ohel2) {
	  for(ohel3=0;ohel3<2;++ohel3) {
	    diag[0]= gamma ?
	      _theFFPVertex->evaluate(_mz2,inters,fin[ihel1],
	      	               	      bcurr[0][ohel2][ohel3]) : 0.;
	    diag[1]= Z0    ?
	      _theFFZVertex->evaluate(_mz2,inters,fin[ihel1],
	      	               	      bcurr[1][ohel2][ohel3]) : 0.;
	    diag[2]= gamma ?
	      _theFFPVertex->evaluate(_mz2,fout[ohel1],interb,
	      	               	      bcurr[0][ohel2][ohel3]) : 0.;
	    diag[3]= Z0    ?
	      _theFFZVertex->evaluate(_mz2,fout[ohel1],interb,
				      bcurr[1][ohel2][ohel3]) : 0.;
	    // diagram contributions
	    me[1] += norm(diag[0]);
	    me[2] += norm(diag[1]);
	    me[3] += norm(diag[2]);
	    me[4] += norm(diag[3]);
	    // total
	    diag[0] += diag[1] + diag[2] + diag[3];
	    me[0]   += norm(diag[0]);
	    if(calc) _me(ihel1,2*ihel2,ohel1,ohel2,ohel3) = diag[0];
	  }
	}
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin = 1./24./4.;
  // and C_F N_c from matrix element
  colspin *= 4.;
  // and for Z decay products
  if(mePartonData()[3]->coloured()) colspin*=3.;
  DVector save;
  for(unsigned int ix=0;ix<5;++ix) {
    me[ix] *= colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0]*UnitRemoval::InvE2;
}

void MEPP2ZJet::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard(5);
  // incoming
  hard[0]=sub->incoming().first;
  hard[1]=sub->incoming().second;
  if((hard[0]->id()<0&&hard[1]->id()<=6)||
     hard[0]->id()==ParticleID::g) swap(hard[0],hard[1]);
  // outgoing
  for(unsigned int ix=0;ix<3;++ix) {
    unsigned int iloc;
    PPtr mother=sub->outgoing()[ix]->parents()[0];
    if(mother&&(mother->id()==ParticleID::gamma||mother->id()==ParticleID::Z0)) {
      if(sub->outgoing()[ix]->id()>0) iloc=3;
      else                            iloc=4;
    }
    else iloc=2;
    hard[iloc]=sub->outgoing()[ix];
  }
  // wavefunctions for the Z decay products
  vector<SpinorBarWaveFunction> lm;
  vector<SpinorWaveFunction>    lp;
  SpinorBarWaveFunction(lm,hard[3],outgoing,true,true);
  SpinorWaveFunction   (lp,hard[4],outgoing,true,true);
  // identify hard process and calculate matrix element
  // q g to q Z
  if(hard[0]->id()<=6&&hard[0]->id()>0&&hard[1]->id()==ParticleID::g) {
    vector<SpinorWaveFunction> fin;
    vector<VectorWaveFunction> gin;
    vector<SpinorBarWaveFunction> fout;
    SpinorWaveFunction    (fin ,hard[0],incoming,false,true);
    VectorWaveFunction    (gin ,hard[1],incoming,false,true,true);
    SpinorBarWaveFunction (fout,hard[2],outgoing,true ,true);
    gin[1]=gin[2];
    qgME(fin,gin,fout,lm,lp,true);
  }
  // qbar g to qbar Z
  else if(hard[0]->id()>=-6&&hard[0]->id()<0&&hard[1]->id()==ParticleID::g) {
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gin;
    vector<SpinorWaveFunction> aout;
    SpinorBarWaveFunction(ain ,hard[0],incoming,false,true);
    VectorWaveFunction   (gin ,hard[1],incoming,false,true,true);
    SpinorWaveFunction   (aout,hard[2],outgoing,true ,true);
    gin[1]=gin[2];
    qbargME(ain,gin,aout,lm,lp,true);
  }
  // q qbar to g Z
  else {
    vector<SpinorWaveFunction>     fin;
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gout;
    SpinorWaveFunction   (fin ,hard[0],incoming,false,true);
    SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
    VectorWaveFunction   (gout,hard[2],outgoing,true ,true,true);
    gout[1]=gout[2];
    qqbarME(fin,ain,gout,lm,lp,true);
  }
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<5;++ix)
    (hard[ix]->spinInfo())->productionVertex(hardvertex);
}
