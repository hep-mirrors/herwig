// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2WJet class.
//

#include "MEPP2WJet.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
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

MEPP2WJet::MEPP2WJet() : _process(0), _maxflavour(5), _plusminus(0),
			 _wdec(0), _widthopt(1)
{}

void MEPP2WJet::doinit() {
  HwMEBase::doinit();
  _wplus  = getParticleData(ThePEG::ParticleID::Wplus );
  _wminus = getParticleData(ThePEG::ParticleID::Wminus);
  // cast the SM pointer to the Herwig SM pointer
  ThePEG::Ptr<Herwig::StandardModel>::transient_const_pointer 
    hwsm=ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardModel>
    ::transient_const_pointer>(standardModel());
  // do the initialisation
  if(!hwsm) 
    throw InitException() << "Must be Herwig::StandardModel in MEPP2WJet::doinit()"
			  << Exception::runerror;
  // set the vertex pointers
  _theFFWVertex = hwsm->vertexFFW();
  _theQQGVertex = hwsm->vertexFFG();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2WJet,HwMEBase>
describeHerwigMEPP2WJet("Herwig::MEPP2WJet", "HwMEHadron.so");

void MEPP2WJet::Init() {

  static ClassDocumentation<MEPP2WJet> documentation
    ("The MEPP2WJet class implements the matrix element for W + jet production");

  static Parameter<MEPP2WJet,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEPP2WJet::_maxflavour, 5, 0, 8, false, false, true);

  static Switch<MEPP2WJet,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2WJet::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcessqqbar
    (interfaceProcess,
     "qqbar",
     "Only include q qbar -> W g process",
     1);
  static SwitchOption interfaceProcessqg
    (interfaceProcess,
     "qg",
     "Only include the q g -> W q process",
     2);
  static SwitchOption interfaceProcessqbarg
    (interfaceProcess,
     "qbarg",
     "Only include the qbar g -> W qbar process",
     3);

  static Switch<MEPP2WJet,unsigned int> interfacePlusMinus
    ("Wcharge",
     "Which intermediate W bosons to include",
     &MEPP2WJet::_plusminus, 0, false, false);
  static SwitchOption interfacePlusMinusAll
    (interfacePlusMinus,
     "Both",
     "Include W+ and W-",
     0);
  static SwitchOption interfacePlusMinusPlus
    (interfacePlusMinus,
     "Plus",
     "Only include W+",
     1);
  static SwitchOption interfacePlusMinusMinus
    (interfacePlusMinus,
     "Minus",
     "Only include W-",
     2);

  static Switch<MEPP2WJet,unsigned int> interfaceWDecay
    ("WDecay",
     "Which processes to include",
     &MEPP2WJet::_wdec, 0, false, false);
  static SwitchOption interfaceWDecayAll
    (interfaceWDecay,
     "All",
     "Include all SM fermions as outgoing particles",
     0);
  static SwitchOption interfaceWDecayQuarks
    (interfaceWDecay,
     "Quarks",
     "Only include outgoing quarks",
     1);
  static SwitchOption interfaceWDecayLeptons
    (interfaceWDecay,
     "Leptons",
     "All include outgoing leptons",
     2);
  static SwitchOption interfaceWDecayElectron
    (interfaceWDecay,
     "Electron",
     "Only include outgoing e nu_e",
     3);
  static SwitchOption interfaceWDecayMuon
    (interfaceWDecay,
     "Muon",
     "Only include outgoing mu nu_mu",
     4);
  static SwitchOption interfaceWDecayTau
    (interfaceWDecay,
     "Tau",
     "Only include outgoing tauu nu_tau",
     5);
  static SwitchOption interfaceWDecayUpDown
    (interfaceWDecay,
     "UpDown",
     "Only include outgoing u dbar/ d ubar",
     6);
  static SwitchOption interfaceWDecayUpStrange
    (interfaceWDecay,
     "UpStrange",
     "Only include outgoing u sbar/ s ubar",
     7);
  static SwitchOption interfaceWDecayUpBottom
    (interfaceWDecay,
     "UpBottom",
     "Only include outgoing u bbar/ b ubar",
     8);
  static SwitchOption interfaceWDecayCharmDown
    (interfaceWDecay,
     "CharmDown",
     "Only include outgoing c dbar/ d cbar",
     9);
  static SwitchOption interfaceWDecayCharmStrange
    (interfaceWDecay,
     "CharmStrange",
     "Only include outgoing c sbar/ s cbar",
     10);
  static SwitchOption interfaceWDecayCharmBottom
    (interfaceWDecay,
     "CharmBottom",
     "Only include outgoing c bbar/ b cbar",
     11);

  static Switch<MEPP2WJet,unsigned int> interfaceWidthOption
    ("WidthOption",
     "The option for handling the width of the off-shell W boson",
     &MEPP2WJet::_widthopt, 1, false, false);
  static SwitchOption interfaceWidthOptionFixedDenominator
    (interfaceWidthOption,
     "FixedDenominator",
     "Use a fixed with in the W propagator but the full matrix element"
     " in the numerator",
     1);
  static SwitchOption interfaceWidthOptionAllRunning
    (interfaceWidthOption,
     "AllRunning",
     "Use a running width in the W propagator and the full matrix "
     "element in the numerator",
     2);

}

void MEPP2WJet::getDiagrams() const {
  // which intgermediates to include
  bool wplus  = _plusminus==0 || _plusminus==1;
  bool wminus = _plusminus==0 || _plusminus==2;
  // possible incoming and outgoing particles
  typedef std::vector<pair<long,long> > Pairvector;
  // possible parents
  Pairvector parentpair;
  parentpair.reserve(6);
  // don't even think of putting 'break' in here!
  switch(_maxflavour) {
  case 5:
    parentpair.push_back(make_pair(ParticleID::b, ParticleID::cbar));
    parentpair.push_back(make_pair(ParticleID::b, ParticleID::ubar));
    [[fallthrough]];
  case 4:
    parentpair.push_back(make_pair(ParticleID::s, ParticleID::cbar));
    parentpair.push_back(make_pair(ParticleID::d, ParticleID::cbar));
    [[fallthrough]];
  case 3:
    parentpair.push_back(make_pair(ParticleID::s, ParticleID::ubar));
    [[fallthrough]];
  case 2:
    parentpair.push_back(make_pair(ParticleID::d, ParticleID::ubar));
    [[fallthrough]];
  default:
    ;
  }
  // possible children
  Pairvector childpair;
  childpair.reserve(9);
  childpair.push_back(make_pair(ParticleID::eminus,   ParticleID::nu_ebar));
  childpair.push_back(make_pair(ParticleID::muminus,  ParticleID::nu_mubar));
  childpair.push_back(make_pair(ParticleID::tauminus, ParticleID::nu_taubar));
  childpair.push_back(make_pair(ParticleID::d, ParticleID::ubar));
  childpair.push_back(make_pair(ParticleID::s, ParticleID::ubar));
  childpair.push_back(make_pair(ParticleID::b, ParticleID::ubar));
  childpair.push_back(make_pair(ParticleID::d, ParticleID::cbar));
  childpair.push_back(make_pair(ParticleID::s, ParticleID::cbar));
  childpair.push_back(make_pair(ParticleID::b, ParticleID::cbar));
  // gluon for diagrams
  tcPDPtr g = getParticleData(ParticleID::g);
  // loop over the children
  bool lepton,quark;
  Pairvector::const_iterator child = childpair.begin();
  for (; child != childpair.end(); ++child) {
    // allowed leptonic decay
    lepton=child->first>10&&
      (_wdec==0||_wdec==2||
       (abs(child->first)-5)/2==int(_wdec));
    // allowed quark decay
    quark =abs(child->second)<10&&
      (_wdec==0||_wdec==1||
       (abs(child->second)==2&&(abs(child->first)+11)/2==int(_wdec))||
       (abs(child->second)==4&&(abs(child->first)+17)/2==int(_wdec)));
    // if decay not allowed skip
    if(!(quark||lepton)) continue;
    // decay products
    tcPDPtr lNeg1 = getParticleData(child->first);
    tcPDPtr lNeg2 = getParticleData(child->second);
    tcPDPtr lPos1 = lNeg2->CC();
    tcPDPtr lPos2 = lNeg1->CC();
    Pairvector::const_iterator parent = parentpair.begin();
    for (; parent != parentpair.end(); ++parent) {
      // parents
      tcPDPtr qNeg1 = getParticleData(parent->first);
      tcPDPtr qNeg2 = getParticleData(parent->second);
      tcPDPtr qPos1 = qNeg2->CC();
      tcPDPtr qPos2 = qNeg1->CC();
      // diagrams
      // q qbar annhilation processes
      if(_process==0||_process==1) {
	// q qbar -> W- g
	if(wminus) {
	  add(new_ptr((Tree2toNDiagram(3), qNeg1, qNeg2, qNeg2, 1, _wminus,
		       2, g,  4, lNeg1, 4, lNeg2, -1)));
	  add(new_ptr((Tree2toNDiagram(3), qNeg1, qNeg1, qNeg2, 2, _wminus,
		       1, g,  4, lNeg1, 4, lNeg2, -2)));
	}
	// q qbar -> W+ g
	if(wplus) {
	  add(new_ptr((Tree2toNDiagram(3), qPos1, qPos2, qPos2, 1, _wplus,
		       2, g,  4, lPos1, 4, lPos2, -3)));
	  add(new_ptr((Tree2toNDiagram(3), qPos1, qPos1, qPos2, 2, _wplus,
		       1, g,  4, lPos1, 4, lPos2, -4)));
	}
      }
      // q g compton
      if(_process==0||_process==2) {
	if(wminus) {
	  add(new_ptr((Tree2toNDiagram(3), qNeg1, qPos1, g    , 1, _wminus,
		       2, qPos1,  4, lNeg1, 4, lNeg2, -5)));
	  add(new_ptr((Tree2toNDiagram(2), qNeg1, g, 1, qNeg1,  3, _wminus,
		       3, qPos1,  4, lNeg1, 4, lNeg2, -6)));
	}
	if(wplus) {
	  add(new_ptr((Tree2toNDiagram(3), qPos1, qNeg1, g,     1, _wplus,
		       2, qNeg1,  4, lPos1, 4, lPos2, -7)));
	  add(new_ptr((Tree2toNDiagram(2), qPos1, g, 1, qNeg1,  3, _wplus,
		       3, qNeg1,  4, lPos1, 4, lPos2, -8)));
	}
      }
      // qbar g compton
      if(_process==0||_process==3) {
	if(wminus) {
	  add(new_ptr((Tree2toNDiagram(3), qNeg2, qPos2, g,     1, _wminus,
		       2, qPos2,  4, lNeg1, 4, lNeg2, -9 )));
	  add(new_ptr((Tree2toNDiagram(2), qNeg2, g, 1, qNeg2,  3, _wminus,
		       3, qPos2,  4, lNeg1, 4, lNeg2, -10)));
	}
	if(wplus) {
	  add(new_ptr((Tree2toNDiagram(3), qPos2, qNeg2, g,     1, _wplus,
		       2, qNeg2,  4, lPos1, 4, lPos2, -11)));
	  add(new_ptr((Tree2toNDiagram(2), qPos2,  g, 1, qPos2, 3, _wplus,
		       3, qNeg2,  4, lPos1, 4, lPos2, -12)));
	}
      }
    }
  }
}

unsigned int MEPP2WJet::orderInAlphaS() const {
  return 1;
}

unsigned int MEPP2WJet::orderInAlphaEW() const {
  return 2;
}

void MEPP2WJet::persistentOutput(PersistentOStream & os) const {
  os << _theFFWVertex << _theQQGVertex << _wplus << _widthopt
     << _wminus << _process << _maxflavour << _plusminus << _wdec;
}

void MEPP2WJet::persistentInput(PersistentIStream & is, int) {
  is >> _theFFWVertex >> _theQQGVertex >> _wplus >> _widthopt
     >> _wminus >> _process >> _maxflavour >> _plusminus >> _wdec;
}

int MEPP2WJet::nDim() const {
  return 5;
}

Selector<const ColourLines *>
MEPP2WJet::colourGeometries(tcDiagPtr diag) const {
  // colour lines for q qbar -> W g
  static const ColourLines cqqbar[4]={ColourLines("1 -2 5,-3 -5"),
				      ColourLines("1 5, -5 2 -3"),
				      ColourLines("1 -2 5,-3 -5,6 -7"),
				      ColourLines("1 5, -5 2 -3,6 -7")};
  // colour lines for q g -> W q
  static const ColourLines cqg   [4]={ColourLines("1 2 -3,3 5"),
				      ColourLines("1 -2,2 3 5"),
				      ColourLines("1 2 -3,3 5,6 -7"),
				      ColourLines("1 -2,2 3 5,6 -7")};
  // colour lines for qbar q -> W qbar
  static const ColourLines cqbarg[4]={ColourLines("-1 -2 3,-3 -5"),
				      ColourLines("-1 2,-2 -3 -5"),
				      ColourLines("-1 -2 3,-3 -5,6 -7"),
				      ColourLines("-1 2,-2 -3 -5,6 -7")};
  // select the correct line
  unsigned int icol = mePartonData()[3]->coloured() ? 2 : 0;
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
  case  1 : case  3:
    sel.insert(1.0, &cqqbar[icol]);
    break;
  case  2 : case  4:
    sel.insert(1.0, &cqqbar[icol+1]);
    break;
  case  5 : case  7:
    sel.insert(1.0, &cqg[icol]);
    break;
  case  6 : case  8:
    sel.insert(1.0, &cqg[icol+1]);
    break;
  case  9 : case 11:
    sel.insert(1.0, &cqbarg[icol]);
    break;
  case 10 : case 12:
    sel.insert(1.0, &cqbarg[icol+1]);
    break;
  }
  return sel;
}

Selector<MEBase::DiagramIndex>
MEPP2WJet::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    int id=abs(diags[i]->id());
    if     (id <= 2 ) sel.insert(meInfo()[id- 1],i);
    else if(id <= 4 ) sel.insert(meInfo()[id- 3],i);
    else if(id <= 6 ) sel.insert(meInfo()[id- 5],i);
    else if(id <= 8 ) sel.insert(meInfo()[id- 7],i);
    else if(id <= 10) sel.insert(meInfo()[id- 9],i);
    else if(id <= 12) sel.insert(meInfo()[id-11],i);
  }
  return sel;
}

Energy2 MEPP2WJet::scale() const {
  return _scale;
}

CrossSection MEPP2WJet::dSigHatDR() const {
  return me2()*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
}

bool MEPP2WJet::generateKinematics(const double * r) {
  // initialize jacobian
  jacobian(1.);
  // cms energy
  Energy ecm=sqrt(sHat());
  // find the right W pointer
  tcPDPtr wdata = mePartonData()[3]->iCharge()+mePartonData()[4]->iCharge() > 0 ? 
    _wplus :_wminus; 
  // first generate the mass of the off-shell gauge boson
  // minimum mass of the 
  tcPDVector ptemp;
  ptemp.push_back(mePartonData()[3]);
  ptemp.push_back(mePartonData()[4]);
  Energy2 minMass2 = max(lastCuts().minSij(mePartonData()[3],mePartonData()[4]),
			 lastCuts().minS(ptemp));
  // minimum pt of the jet
  Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
		     lastCuts().minKT(wdata));
  // maximum mass of the gauge boson so pt is possible
  Energy2 maxMass2 = min(ecm*(ecm-2.*ptmin),lastCuts().maxS(ptemp));
  if(maxMass2<=ZERO||minMass2<ZERO) return false;
  // also impose the limits from the ParticleData object
  minMass2 = max(minMass2,sqr(wdata->massMin()));
  maxMass2 = min(maxMass2,sqr(wdata->massMax()));
  // return if not kinematically possible
  if(minMass2>maxMass2) return false;
  // generation of the mass
  Energy  M(wdata->mass()),Gamma(wdata->width());
  Energy2 M2(sqr(M)),MG(M*Gamma);
  double rhomin = atan2((minMass2-M2),MG);
  double rhomax = atan2((maxMass2-M2),MG);
  _mw2=M2+MG*tan(rhomin+r[1]*(rhomax-rhomin));
  Energy mw=sqrt(_mw2);
  // jacobian
  jacobian(jacobian()*(sqr(_mw2-M2)+sqr(MG))/MG*(rhomax-rhomin)/sHat());
  // set the masses of the outgoing particles in the 2-2 scattering
  meMomenta()[2].setMass(ZERO);
  Lorentz5Momentum pw(mw);
  // generate the polar angle of the hard scattering
  double ctmin(-1.0), ctmax(1.0);
  Energy q(ZERO);
  try {
    q = SimplePhaseSpace::getMagnitude(sHat(), meMomenta()[2].mass(),mw);
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
  // momenta of particle in hard scattering
  Energy pt = q*sqrt(1.0-sqr(cth));
  double phi=2.0*Constants::pi*r[2];
  meMomenta()[2].setVect(Momentum3( pt*sin(phi), pt*cos(phi), q*cth));
  pw.setVect(            Momentum3(-pt*sin(phi),-pt*cos(phi),-q*cth));
  meMomenta()[2].rescaleEnergy();
  pw.rescaleEnergy();
  // set the scale
  _scale = _mw2+sqr(pt);
  // generate the momenta of the W decay products
  meMomenta()[3].setMass(mePartonData()[3]->mass());
  meMomenta()[4].setMass(mePartonData()[4]->mass());
  Energy q2 = ZERO;
  try {
    q2 = SimplePhaseSpace::getMagnitude(_mw2, meMomenta()[3].mass(),
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
  Boost boostv(pw.boostVector());
  pl[0].boost(boostv);
  pl[1].boost(boostv);
  meMomenta()[3] = pl[0];
  meMomenta()[4] = pl[1]; 
  // check passes all the cuts
  vector<LorentzMomentum> out(3);
  out[0] = meMomenta()[2];
  out[1] = meMomenta()[3];
  out[2] = meMomenta()[4];
  tcPDVector tout(3);
  tout[0] = mePartonData()[2];
  tout[1] = mePartonData()[3];
  tout[2] = mePartonData()[4];
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
    return false;
  // jacobian
  jacobian((pq/sHat())*Constants::pi*jacobian()/8./sqr(Constants::pi)*q2/mw);
  return true;	    
}

double MEPP2WJet::me2() const {
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
  // q g to q W
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
  // qbar g to qbar W
  else if(mePartonData()[0]->id()>=-6&&mePartonData()[0]->id()<0&&
	  mePartonData()[1]->id()==ParticleID::g) {
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gin;
    vector<SpinorWaveFunction> aout;
    SpinorBarWaveFunction qbin (meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction    glin (meMomenta()[1],mePartonData()[1],incoming);
    SpinorWaveFunction    qbout(meMomenta()[2],mePartonData()[2],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qbin.reset(ix) ; ain.push_back(qbin);
      glin.reset(2*ix) ; gin.push_back(glin);
      qbout.reset(ix);aout.push_back(qbout);
    }
    output=qbargME(ain,gin,aout,lm,lp);
  }
  // q qbar to g W
  else {
    vector<SpinorWaveFunction>     fin;
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gout;
    SpinorWaveFunction    qin (meMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction qbin(meMomenta()[1],mePartonData()[1],incoming);
    VectorWaveFunction   glout(meMomenta()[2],mePartonData()[2],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)    ;  fin.push_back(qin);
      qbin.reset(ix)   ;  ain.push_back(qbin);
      glout.reset(2*ix); gout.push_back(glout);
    }
    output=qqbarME(fin,ain,gout,lm,lp);
  }
  return output*sHat();
}

InvEnergy2 MEPP2WJet::qqbarME(vector<SpinorWaveFunction> & fin,
			      vector<SpinorBarWaveFunction> & ain,
			      vector<VectorWaveFunction> & gout,
			      vector<SpinorBarWaveFunction> & lm,
			      vector<SpinorWaveFunction> & lp,
			      bool calc) const {
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
					     PDT::Spin1,PDT::Spin1Half,
					     PDT::Spin1Half));
  // some integers
  unsigned int ihel1,ihel2,ohel1,ohel2,ohel3;
  // find the right W pointer
  tcPDPtr wdata = mePartonData()[3]->iCharge()+mePartonData()[4]->iCharge() > 0
    ? _wplus :_wminus; 
  // compute the W current for speed
  VectorWaveFunction bcurr[2][2];
  for(ohel2=0;ohel2<2;++ohel2) {
    for(ohel3=0;ohel3<2;++ohel3) {
      bcurr[ohel2][ohel3] = _theFFWVertex->evaluate(_mw2,_widthopt,wdata,
						    lp[ohel3],lm[ohel2]);
    }
  }
  double me[3]={0.,0.,0.};
  Complex diag[2];
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
	    diag[0] = _theFFWVertex->evaluate(_mw2,fin[ihel1],interb,
					      bcurr[ohel2][ohel3]);
	    diag[1] = _theFFWVertex->evaluate(_mw2,inters,ain[ihel2],
					      bcurr[ohel2][ohel3]);
	    // diagram contributions
	    me[1] += norm(diag[0]);
	    me[2] += norm(diag[1]);
	    // total
	    diag[0] += diag[1];
	    me[0]   += norm(diag[0]);
	    if(calc) _me(ihel1,ihel2,2*ohel1,ohel2,ohel3) = diag[0];
	  }
	}
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin=1./9./4.;
  // and C_F N_c from matrix element
  colspin *= 4.;
  // colour factor for the W decay
  if(mePartonData()[3]->coloured()) colspin*=3.;
  DVector save;
  for(unsigned int ix=0;ix<3;++ix) {
    me[ix]*=colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0] * UnitRemoval::InvE2;
}

InvEnergy2 MEPP2WJet::qgME(vector<SpinorWaveFunction> & fin,
			   vector<VectorWaveFunction> & gin,
			   vector<SpinorBarWaveFunction> & fout,
			   vector<SpinorBarWaveFunction> & lm,
			   vector<SpinorWaveFunction> & lp,
			   bool calc) const {
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
					     PDT::Spin1Half,PDT::Spin1Half,
					     PDT::Spin1Half));
  // find the right W pointer
  tcPDPtr wdata = mePartonData()[3]->iCharge()+mePartonData()[4]->iCharge() > 0 ? 
    _wplus :_wminus; 
  // some integers
  unsigned int ihel1,ihel2,ohel1,ohel2,ohel3;
  // compute the leptonic W current for speed
  VectorWaveFunction bcurr[2][2];
  for(ohel2=0;ohel2<2;++ohel2) {
    for(ohel3=0;ohel3<2;++ohel3) {
      bcurr[ohel2][ohel3] = _theFFWVertex->evaluate(_mw2,_widthopt,wdata,
						    lp[ohel3],lm[ohel2]);
    }
  }
  // compute the matrix elements
  double me[3]={0.,0.,0.};
  Complex diag[2];
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;
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
	    diag[0]=_theFFWVertex->evaluate(_mw2,fin[ihel1],interb,
					    bcurr[ohel2][ohel3]);
	    diag[1]=_theFFWVertex->evaluate(_mw2,inters,fout[ohel1],
					    bcurr[ohel2][ohel3]);
	    // diagram contributions
	    me[1] += norm(diag[0]);
	    me[2] += norm(diag[1]);
	    // total
	    diag[0] += diag[1];
	    me[0]   += norm(diag[0]);
	    if(calc) _me(ihel1,2*ihel2,ohel1,ohel2,ohel3) = diag[0];
	  }
	}
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin=1./24./4.;
  // and C_F N_c from matrix element
  colspin *=4.;
  // colour factor for the W decay
  if(mePartonData()[3]->coloured()) colspin*=3.;
  DVector save;
  for(unsigned int ix=0;ix<3;++ix) {
    me[ix]*=colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0] * UnitRemoval::InvE2;
}

InvEnergy2 MEPP2WJet::qbargME(vector<SpinorBarWaveFunction> & fin,
			     vector<VectorWaveFunction> & gin,
			     vector<SpinorWaveFunction> & fout,
			     vector<SpinorBarWaveFunction> & lm,
			     vector<SpinorWaveFunction> & lp,
			     bool calc) const {
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
					     PDT::Spin1Half,PDT::Spin1Half,
					     PDT::Spin1Half));
  // find the right W pointer
  tcPDPtr wdata = mePartonData()[3]->iCharge()+mePartonData()[4]->iCharge() > 0 ?
    _wplus :_wminus; 
  // some integers
  unsigned int ihel1,ihel2,ohel1,ohel2,ohel3;
  // compute the leptonic W current for speed
  VectorWaveFunction bcurr[2][2];
  for(ohel2=0;ohel2<2;++ohel2) {
    for(ohel3=0;ohel3<2;++ohel3) {
      bcurr[ohel2][ohel3] = _theFFWVertex->evaluate(_mw2,_widthopt,wdata,
						    lp[ohel3],lm[ohel2]);
    }
  }
  // compute the matrix elements
  double me[3]={0.,0.,0.};
  Complex diag[2];
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;
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
	    diag[0]= _theFFWVertex->evaluate(_mw2,inters,fin[ihel1],
					     bcurr[ohel2][ohel3]);
	    diag[1]= _theFFWVertex->evaluate(_mw2,fout[ohel1],interb,
					     bcurr[ohel2][ohel3]);
	    // diagram contributions
	    me[1] += norm(diag[0]);
	    me[2] += norm(diag[1]);
	    // total
	    diag[0] += diag[1];
	    me[0]   += norm(diag[0]);
	    if(calc) _me(ihel1,2*ihel2,ohel1,ohel2,ohel3) = diag[0];
	  }
	}
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin=1./24./4.;
  // and C_F N_c from matrix element
  colspin *= 4.;
  // colour factor for the W decay
  if(mePartonData()[3]->coloured()) colspin*=3.;
  DVector save;
  for(unsigned int ix=0;ix<3;++ix) {
    me[ix]*=colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0] * UnitRemoval::InvE2;
}

void MEPP2WJet::constructVertex(tSubProPtr sub) {
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
    if(mother&&abs(mother->id())==ParticleID::Wplus) {
      if(sub->outgoing()[ix]->id()>0) iloc=3;
      else                            iloc=4;
    }
    else iloc=2;
    hard[iloc]=sub->outgoing()[ix];
  }
  // wavefunctions for the W decay products
  vector<SpinorBarWaveFunction> lm;
  vector<SpinorWaveFunction>    lp;
  SpinorBarWaveFunction(lm,hard[3],outgoing,true,true);
  SpinorWaveFunction   (lp,hard[4],outgoing,true,true);
  // identify hard process and calculate matrix element
  // q g to q W
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
  // qbar g to qbar W
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
  // q qbar to g W
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
