// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2gZ2SleptonsJet class.
//

#include "MEPP2gZ2SleptonsJet.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Models/Susy/SusyBase.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig++/MatrixElement/HardVertex.h"

using namespace Herwig;

MEPP2gZ2SleptonsJet::MEPP2gZ2SleptonsJet() 
  : _process(0), _maxflavour(5), _zdec(0),
    _gammaZ(0), _widthopt(1), _pprob(0.5)
{}

void MEPP2gZ2SleptonsJet::doinit() {
  MEBase::doinit();
  _z0    = getParticleData(ThePEG::ParticleID::Z0   );
  _gamma = getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  ThePEG::Ptr<Herwig::SusyBase>::transient_const_pointer 
    hwsm=ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::SusyBase>
    ::transient_const_pointer>(standardModel());
  // do the initialisation
  if(!hwsm) 
    throw InitException() << "Must be Herwig::StandardModel in"
			  << " MEPP2gZ2SleptonsJet::doinit()"
			  << Exception::runerror;
  // set the vertex pointers
  _theFFZVertex = hwsm->vertexFFZ();
  _theFFPVertex = hwsm->vertexFFP();
  _theQQGVertex = hwsm->vertexFFG();
  _theVSSVertex = hwsm->vertexWSFSF();
}

ClassDescription<MEPP2gZ2SleptonsJet> MEPP2gZ2SleptonsJet::initMEPP2gZ2SleptonsJet;
// Definition of the static class description member.

void MEPP2gZ2SleptonsJet::Init() {

  static ClassDocumentation<MEPP2gZ2SleptonsJet> documentation
    ("The MEPP2gZ2SleptonsJet class implements the matrix element"
     " for Z/gamma+ jet production");

  static Parameter<MEPP2gZ2SleptonsJet,int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEPP2gZ2SleptonsJet::_maxflavour, 5, 0, 5, false, false, true);

  static Switch<MEPP2gZ2SleptonsJet,int> interfaceZDec
    ("ZDecay",
     "Which processes to generate",
     &MEPP2gZ2SleptonsJet::_zdec, 0, false, false);
  static SwitchOption interfaceZDecAll
    (interfaceZDec,
     "All",
     "Generate all the processes",
     0);
  static SwitchOption interfaceZDece_L
    (interfaceZDec,
     "e_L",
     "Only produce ~e_L",
     1);
  static SwitchOption interfaceZDecmu_L
    (interfaceZDec,
     "mu_L",
     "Onle produce ~mu_L",
     2);
  static SwitchOption interfaceZDectau_1
    (interfaceZDec,
     "tau_1",
     "Only produce tau_1 pairs",
     3);
  static SwitchOption interfaceZDece_R
    (interfaceZDec,
     "e_R",
     "Only produce e_R",
     4);
  static SwitchOption interfaceZDecmu_R
    (interfaceZDec,
     "mu_R",
     "Only produce ~mu_R",
     5);
  static SwitchOption interfaceZDectau_2
    (interfaceZDec,
     "tau_2",
     "Only produce tau_2 pairs",
     6);
  static SwitchOption interfaceZDecnu_e
    (interfaceZDec,
     "nu_e",
     "Only product ~nu_e",
     7);
  static SwitchOption interfaceZDecnu_mu
    (interfaceZDec,
     "nu_mu",
     "Only produce ~nu_mu",
     8);
  static SwitchOption interfaceZDecnu_tau
    (interfaceZDec,
     "nu_tau",
     "Only produce ~nu_tau",
     9);
  static SwitchOption interfaceZDecMixedTau
    (interfaceZDec,
     "MixedTau",
     "Only produce mixing tau_1 tau_2 pairs",
     10);

  static Switch<MEPP2gZ2SleptonsJet,unsigned int> interfaceProcess
    ("Process",
     "Which processes to generate",
     &MEPP2gZ2SleptonsJet::_process, 0, false, false);
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

  static Parameter<MEPP2gZ2SleptonsJet,double> interfacePhotonProbablity
    ("PhotonProbablity",
     "Probability for using the \\f$1/s^2\\f$ piece for the"
     " generation of the gauge boson mass",
     &MEPP2gZ2SleptonsJet::_pprob, 0.5, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<MEPP2gZ2SleptonsJet,unsigned int> interfaceGammaZ
    ("GammaZ",
     "Which terms to include",
     &MEPP2gZ2SleptonsJet::_gammaZ, 0, false, false);
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

  static Switch<MEPP2gZ2SleptonsJet,unsigned int> interfaceWidthOption
    ("WidthOption",
     "The option for handling the width of the off-shell W boson",
     &MEPP2gZ2SleptonsJet::_widthopt, 1, false, false);
  static SwitchOption interfaceWidthOptionFixedDenominator
    (interfaceWidthOption,
     "FixedDenominator",
     "Use a fxied with in the W propagator but the full matrix element"
     " in the numerator",
     1);
  static SwitchOption interfaceWidthOptionAllRunning
    (interfaceWidthOption,
     "AllRunning",
     "Use a running width in the Z propagator and the full matrix "
     "element in the numerator",
     2);

}

void MEPP2gZ2SleptonsJet::getDiagrams() const {
  // which intermediates to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // pointer for gluon
  tcPDPtr g = getParticleData(ParticleID::g);
  for ( int ix=11; ix<17; ++ix ) {
    unsigned int ymax = 1;
    if(ix==11 || ix == 13 ) ymax = 2;
    else if(ix==15)         ymax = 4;
    for(unsigned int iy=0;iy<ymax;++iy) {
      tcPDPtr lm,lp;
      if(iy==0) {
	lm = getParticleData(1000000+ix);
	lp = lm->CC();
      }
      else if(iy==1) {
	lm = getParticleData(2000000+ix);
	lp = lm->CC();
      }
      else if(iy==2) {
	lm = getParticleData( 1000000+ix);
	lp = getParticleData(-2000000-ix);
      }
      else if(iy==3) {
	lm = getParticleData( 2000000+ix);
	lp = getParticleData(-1000000-ix);
      }
      // if not a valid process continue
      if((_zdec==1 && lm->id()!=1000011) ||
	 (_zdec==2 && lm->id()!=1000013) ||
	 (_zdec==3 && (lm->id()!=1000015 || lp->id()!=-1000015)) ||
	 (_zdec==4 && lm->id()!=2000011) ||
	 (_zdec==5 && lm->id()!=2000013) ||
	 (_zdec==6 && (lm->id()!=2000015 || lp->id()!=-2000015)) ||
	 (_zdec==7 && lm->id()!=1000012) ||
	 (_zdec==8 && lm->id()!=1000014) ||
	 (_zdec==9 && lm->id()!=1000016) ||
	 (_zdec==10 && (ix!=15 || lm->id()==-lp->id()))) continue;
      // pointer for Z decay products
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
}

unsigned int MEPP2gZ2SleptonsJet::orderInAlphaS() const {
  return 1;
}

unsigned int MEPP2gZ2SleptonsJet::orderInAlphaEW() const {
  return 2;
}

void MEPP2gZ2SleptonsJet::persistentOutput(PersistentOStream & os) const {
  os << _theFFZVertex << _theFFPVertex << _theQQGVertex << _theVSSVertex
     << _z0 << _widthopt
     << _gamma << _process << _maxflavour << _zdec << _pprob << _gammaZ;
}

void MEPP2gZ2SleptonsJet::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZVertex >> _theFFPVertex >> _theQQGVertex >> _theVSSVertex
     >> _z0 >> _widthopt
     >> _gamma >> _process >> _maxflavour >> _zdec >> _pprob >> _gammaZ;
}

int MEPP2gZ2SleptonsJet::nDim() const {
  return 5;
}

Selector<const ColourLines *>
MEPP2gZ2SleptonsJet::colourGeometries(tcDiagPtr diag) const {
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
MEPP2gZ2SleptonsJet::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    int id=abs(diags[i]->id());
    if     (id <= 4 ) sel.insert(meInfo()[id-1],i);
    else if(id <= 8 ) sel.insert(meInfo()[id-5],i);
    else if(id <= 12) sel.insert(meInfo()[id-9],i);
  }
  return sel;
}

Energy2 MEPP2gZ2SleptonsJet::scale() const {
  return _scale;
}

CrossSection MEPP2gZ2SleptonsJet::dSigHatDR() const {
  return me2()*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
}

bool MEPP2gZ2SleptonsJet::generateKinematics(const double * r) {
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
  if(maxMass2<minMass2) return false;
  // generation of the mass
  Energy  M(_z0->mass()),Gamma(_z0->width());
  Energy2 M2(sqr(M)),MG(M*Gamma);
  double rhomin = atan((minMass2-M2)/MG);
  double rhomax = atan((maxMass2-M2)/MG);
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
  catch ( ImpossibleKinematics ) {
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
  } catch ( ImpossibleKinematics ) {
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

double MEPP2gZ2SleptonsJet::getCosTheta(double ctmin, double ctmax, const double r) {
  double cth = 0.0;
  double zmin = 0.5*(1.0 - ctmax);
  double zmax = 0.5*(1.0 - ctmin);
  if ( zmin <= 0.0 || zmax >= 1.0 ) {
    jacobian((ctmax - ctmin)*jacobian());
    cth = ctmin + r*(ctmax-ctmin);
  } else {
    double A1 = (2.0*zmax - 1.0)/(zmax*(1.0-zmax));
    double A0 = (2.0*zmin - 1.0)/(zmin*(1.0-zmin));
    double A = r*(A1 - A0) + A0;
    double z = A < 2.0? 2.0/(sqrt(sqr(A) + 4.0) + 2 - A):
      0.5*(A - 2.0 + sqrt(sqr(A) + 4.0))/A;
    cth = 1.0 - 2.0*z;
    jacobian((2.0*(A1 - A0)*sqr(z)*sqr(1.0 - z)/(sqr(z) + sqr(1.0 - z)))*jacobian());
  }
  return cth;
}  

double MEPP2gZ2SleptonsJet::me2() const {
  InvEnergy2 output(ZERO);
  // construct spinors for the leptons (always the same)
  ScalarWaveFunction lm(meMomenta()[3],mePartonData()[3],outgoing);
  ScalarWaveFunction lp(meMomenta()[4],mePartonData()[4],outgoing);
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

InvEnergy2 MEPP2gZ2SleptonsJet::qqbarME(vector<SpinorWaveFunction> & fin,
			      vector<SpinorBarWaveFunction> & ain,
			      vector<VectorWaveFunction> & gout,
			      ScalarWaveFunction & lm,
			      ScalarWaveFunction & lp,
			      bool calc) const {
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
					     PDT::Spin1,PDT::Spin0,
					     PDT::Spin0));
  // diagrams to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // some integers
  unsigned int ihel1,ihel2,ohel1;
  // compute the leptonic photon and Z currents for speed
  VectorWaveFunction bcurr[2];
  // photon current
  if(gamma) bcurr[0]= _theVSSVertex->evaluate(_mz2,1,_gamma,lp,lm);
  // Z current
  if(Z0)    bcurr[1]= _theVSSVertex->evaluate(_mz2,_widthopt,_z0,lp,lm);
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
	diag[0] = gamma ? 
	  _theFFPVertex->evaluate(_mz2,fin[ihel1],interb,bcurr[0]) : 0.;
	diag[1]= Z0 ? 
	  _theFFZVertex->evaluate(_mz2,fin[ihel1],interb,bcurr[1]) : 0.;
	diag[2]= gamma ? 
	  _theFFPVertex->evaluate(_mz2,inters,ain[ihel2],bcurr[0]) : 0.;
	diag[3]= Z0 ? 
	  _theFFZVertex->evaluate(_mz2,inters,ain[ihel2],bcurr[1]) : 0.;
	// diagram contributions
	me[1] += norm(diag[0]);
	me[2] += norm(diag[1]);
	me[3] += norm(diag[2]);
	me[4] += norm(diag[3]);
	// total
	diag[0] += diag[1] + diag[2] + diag[3];
	me[0]   += norm(diag[0]);
	if(calc) _me(ihel1,ihel2,2*ohel1,0,0) = diag[0];
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin = 1./9./4.;
  // and C_F N_c from matrix element
  colspin *= 4.;
  DVector save;
  for(unsigned int ix=0;ix<5;++ix) {
    me[ix] *= colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0]*UnitRemoval::InvE2;
}

InvEnergy2 MEPP2gZ2SleptonsJet::qgME(vector<SpinorWaveFunction> & fin,
			   vector<VectorWaveFunction> & gin,
			   vector<SpinorBarWaveFunction> & fout,
			   ScalarWaveFunction & lm,
			   ScalarWaveFunction & lp,
			   bool calc) const {
  // diagrams to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
					     PDT::Spin1Half,PDT::Spin0,
					     PDT::Spin0));
  // some integers
  unsigned int ihel1,ihel2,ohel1;
  // compute the leptonic photon and Z currents for speed
  VectorWaveFunction bcurr[2];
  // photon current
  if(gamma) bcurr[0]= _theVSSVertex->evaluate(_mz2,1,_gamma,lp,lm);
  // Z current
  if(Z0)    bcurr[1]= _theVSSVertex->evaluate(_mz2,_widthopt,_z0,lp,lm);
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
	interb=_theQQGVertex->evaluate(_scale,5,mePartonData()[2],
				       fout[ohel1],gin[ihel2]);
	inters=_theQQGVertex->evaluate(_scale,5,mePartonData()[0],
				       fin[ihel1],gin[ihel2]);
	diag[0]=gamma ?
	  _theFFPVertex->evaluate(_mz2,fin[ihel1],interb,bcurr[0]) : 0.;
	diag[1]=Z0    ?
	  _theFFZVertex->evaluate(_mz2,fin[ihel1],interb,bcurr[1]) : 0.;
	diag[2]=gamma ?
	  _theFFPVertex->evaluate(_mz2,inters,fout[ohel1],bcurr[0]) : 0.;
	diag[3]=Z0    ?
	  _theFFZVertex->evaluate(_mz2,inters,fout[ohel1],bcurr[1]) : 0.;
	// diagram contributions
	me[1] += norm(diag[0]);
	me[2] += norm(diag[1]);
	me[3] += norm(diag[2]);
	me[4] += norm(diag[3]);
	// total
	diag[0] += diag[1] + diag[2] + diag[3];
	me[0]   += norm(diag[0]);
	if(calc) _me(ihel1,2*ihel2,ohel1) = diag[0];
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin = 1./24./4.;
  // and C_F N_c from matrix element
  colspin *= 4.;
  DVector save;
  for(unsigned int ix=0;ix<5;++ix) {
    me[ix] *= colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0]*UnitRemoval::InvE2;
}

InvEnergy2 MEPP2gZ2SleptonsJet::qbargME(vector<SpinorBarWaveFunction> & fin,
					vector<VectorWaveFunction> & gin,
					vector<SpinorWaveFunction> & fout,
					ScalarWaveFunction & lm,
					ScalarWaveFunction & lp,
					bool calc) const {
  // diagrams to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
					     PDT::Spin1Half,PDT::Spin0,
					     PDT::Spin0));
  // some integers
  unsigned int ihel1,ihel2,ohel1;
  // compute the leptonic photon and Z currents for speed
  VectorWaveFunction bcurr[2];
  // photon current
  if(gamma) bcurr[0]= _theVSSVertex->evaluate(_mz2,1,_gamma,lp,lm);
  // Z current
  if(Z0)    bcurr[1]= _theVSSVertex->evaluate(_mz2,_widthopt,_z0,lp,lm);
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
	inters=_theQQGVertex->evaluate(_scale,5,mePartonData()[2],
				       fout[ohel1],gin[ihel2]);
	interb=_theQQGVertex->evaluate(_scale,5,mePartonData()[0],
				       fin[ihel1],gin[ihel2]);
	diag[0]= gamma ?
	  _theFFPVertex->evaluate(_mz2,inters,fin[ihel1],bcurr[0]) : 0.;
	diag[1]= Z0    ?
	  _theFFZVertex->evaluate(_mz2,inters,fin[ihel1],bcurr[1]) : 0.;
	diag[2]= gamma ?
	  _theFFPVertex->evaluate(_mz2,fout[ohel1],interb,bcurr[0]) : 0.;
	diag[3]= Z0    ?
	  _theFFZVertex->evaluate(_mz2,fout[ohel1],interb,bcurr[1]) : 0.;
	// diagram contributions
	me[1] += norm(diag[0]);
	me[2] += norm(diag[1]);
	me[3] += norm(diag[2]);
	me[4] += norm(diag[3]);
	// total
	diag[0] += diag[1] + diag[2] + diag[3];
	me[0]   += norm(diag[0]);
	if(calc) _me(ihel1,2*ihel2,ohel1) = diag[0];
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin = 1./24./4.;
  // and C_F N_c from matrix element
  colspin *= 4.;
  DVector save;
  for(unsigned int ix=0;ix<5;++ix) {
    me[ix] *= colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0]*UnitRemoval::InvE2;
}

void MEPP2gZ2SleptonsJet::constructVertex(tSubProPtr sub) {
//   // extract the particles in the hard process
//   ParticleVector hard(5);
//   // incoming
//   hard[0]=sub->incoming().first;
//   hard[1]=sub->incoming().second;
//   if((hard[0]->id()<0&&hard[1]->id()<=6)||
//      hard[0]->id()==ParticleID::g) swap(hard[0],hard[1]);
//   // outgoing
//   for(unsigned int ix=0;ix<3;++ix) {
//     unsigned int iloc;
//     PPtr mother=sub->outgoing()[ix]->parents()[0];
//     if(mother&&(mother->id()==ParticleID::gamma||mother->id()==ParticleID::Z0)) {
//       if(sub->outgoing()[ix]->id()>0) iloc=3;
//       else                            iloc=4;
//     }
//     else iloc=2;
//     hard[iloc]=sub->outgoing()[ix];
//   }
//   // wavefunctions for the Z decay products
//   vector<SpinorBarWaveFunction> lm;
//   vector<SpinorWaveFunction>    lp;
//   SpinorBarWaveFunction(lm,hard[3],outgoing,true,true);
//   SpinorWaveFunction   (lp,hard[4],outgoing,true,true);
//   // identify hard process and calculate matrix element
//   // q g to q Z
//   if(hard[0]->id()<=6&&hard[0]->id()>0&&hard[1]->id()==ParticleID::g) {
//     vector<SpinorWaveFunction> fin;
//     vector<VectorWaveFunction> gin;
//     vector<SpinorBarWaveFunction> fout;
//     SpinorWaveFunction    (fin ,hard[0],incoming,false,true);
//     VectorWaveFunction    (gin ,hard[1],incoming,false,true,true);
//     SpinorBarWaveFunction (fout,hard[2],outgoing,true ,true);
//     gin[1]=gin[2];
//     qgME(fin,gin,fout,lm,lp,true);
//   }
//   // qbar g to qbar Z
//   else if(hard[0]->id()>=-6&&hard[0]->id()<0&&hard[1]->id()==ParticleID::g) {
//     vector<SpinorBarWaveFunction>  ain;
//     vector<VectorWaveFunction> gin;
//     vector<SpinorWaveFunction> aout;
//     SpinorBarWaveFunction(ain ,hard[0],incoming,false,true);
//     VectorWaveFunction   (gin ,hard[1],incoming,false,true,true);
//     SpinorWaveFunction   (aout,hard[2],outgoing,true ,true);
//     gin[1]=gin[2];
//     qbargME(ain,gin,aout,lm,lp,true);
//   }
//   // q qbar to g Z
//   else {
//     vector<SpinorWaveFunction>     fin;
//     vector<SpinorBarWaveFunction>  ain;
//     vector<VectorWaveFunction> gout;
//     SpinorWaveFunction   (fin ,hard[0],incoming,false,true);
//     SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
//     VectorWaveFunction   (gout,hard[2],outgoing,true ,true,true);
//     gout[1]=gout[2];
//     qqbarME(fin,ain,gout,lm,lp,true);
//   }
//   // construct the vertex
//   HardVertexPtr hardvertex=new_ptr(HardVertex());
//   // set the matrix element for the vertex
//   hardvertex->ME(_me);
//   // set the pointers and to and from the vertex
//   for(unsigned int ix=0;ix<5;++ix)
//     dynamic_ptr_cast<SpinfoPtr>(hard[ix]->spinInfo())->setProductionVertex(hardvertex);
}
