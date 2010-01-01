// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2GammaGammaPowhegNEW class.
//

#include "MEPP2GammaGammaPowhegNEW.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/HardTree.h"
#include <numeric>

using namespace Herwig;

MEPP2GammaGammaPowhegNEW::MEPP2GammaGammaPowhegNEW() 
  : process_(0), maxFlavour_(5), QCDOption_(2),
    QCDScale_(100.*GeV2),QEDOption_(1), QEDScale_(100.*GeV2) {
  massOption(true ,0);
  massOption(false,0);
}

ClassDescription<MEPP2GammaGammaPowhegNEW> 
MEPP2GammaGammaPowhegNEW::initMEPP2GammaGammaPowhegNEW;
// Definition of the static class description member.

void MEPP2GammaGammaPowhegNEW::Init() {

  static ClassDocumentation<MEPP2GammaGammaPowhegNEW> documentation
    ("The MEPP2GammaGammaPowhegNEW class implements the next-to-leading order"
     "matrix elements for photon pair production in the POWHEG scheme");

  static Parameter<MEPP2GammaGammaPowhegNEW,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the quarks in the process",
     &MEPP2GammaGammaPowhegNEW::maxFlavour_, 5, 1, 5,
     false, false, Interface::limited);


  static Switch<MEPP2GammaGammaPowhegNEW,unsigned int> interfaceProcess
    ("Process",
     "Which types of processes to include",
     &MEPP2GammaGammaPowhegNEW::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all leading and next-to-leading order processes"
     " and gluon gluon -> gamma gamma via a quark box.",
     0);
  static SwitchOption interfaceProcessNLO
    (interfaceProcess,
     "NLO",
     "Include all the leading  and next-to-leading order processes",
     1);
  static SwitchOption interfaceProcessGammaGamma
    (interfaceProcess,
     "GammaGamma",
     "Include the leading-order photon pair production processes",
     2);
  static SwitchOption interfaceProcessGammaJet
    (interfaceProcess,
     "GammaJet",
     "Include the leading order gamma+jet production processes",
     3);
  static SwitchOption interfaceProcessNLOGammaGamma
    (interfaceProcess,
     "NLOGammaGamma",
     "Include the NLO photon pair production processes",
     4);
  static SwitchOption interfaceProcess3Body
    (interfaceProcess,
     "3Body",
     "Include the three-body production processes",
     5);
  static SwitchOption interfaceProcessBox
    (interfaceProcess,
     "Box",
     "Include gluon gluon -> gamma gamma via a quark box.",
     6);

  static Switch<MEPP2GammaGammaPowhegNEW,unsigned int> interfaceQCDOption
    ("QCDOption",
     "Option for the QCD supression factor",
     &MEPP2GammaGammaPowhegNEW::QCDOption_, 2, false, false);
  static SwitchOption interfaceQCDOptionNone
    (interfaceQCDOption,
     "None",
     "Don't include the factor",
     0);
  static SwitchOption interfaceQCDOptionFixedScale
    (interfaceQCDOption,
     "FixedScale",
     "Include the factor with a fixed scale",
     1);
  static SwitchOption interfaceQCDOptionVariableScale
    (interfaceQCDOption,
     "VariableScale",
     "Include the factor using the mass of the photon pair as a variable scale",
     2);

  static Parameter<MEPP2GammaGammaPowhegNEW,Energy2> interfaceQCDScale
    ("QCDScale",
     "Fxied scale for the QCD supression factor",
     &MEPP2GammaGammaPowhegNEW::QCDScale_, GeV2, 100.0*GeV2, 10.0*GeV2, 14000.0*GeV2,
     false, false, Interface::limited);

  static Switch<MEPP2GammaGammaPowhegNEW,unsigned int> interfaceQEDOption
    ("QEDOption",
     "Option for the QED supression factor",
     &MEPP2GammaGammaPowhegNEW::QEDOption_, 1, false, false);
  static SwitchOption interfaceQEDOptionNone
    (interfaceQEDOption,
     "None",
     "Don't include the factor",
     0);
  static SwitchOption interfaceQEDOptionFixedScale
    (interfaceQEDOption,
     "FixedScale",
     "Include the factor with a fixed scale",
     1);
  static SwitchOption interfaceQEDOptionVariableScale
    (interfaceQEDOption,
     "VariableScale",
     "Include the factor using the mass of the photon pair as a variable scale",
     2);

  static Parameter<MEPP2GammaGammaPowhegNEW,Energy2> interfaceQEDScale
    ("QEDScale",
     "Fxied scale for the QED supression factor",
     &MEPP2GammaGammaPowhegNEW::QEDScale_, GeV2, 100.0*GeV2, 10.0*GeV2, 14000.0*GeV2,
     false, false, Interface::limited);

}

void MEPP2GammaGammaPowhegNEW::doinit() {
  // get the vedrtex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!hwsm) InitException() << "Wrong type of StandardModel object in "
			    << "MEPP2GammaGammaPowhegNEW::doinit() the Herwig++"
			    << " version must be used" 
			    << Exception::runerror;
  gluonVertex_  = hwsm->vertexFFG();
  photonVertex_ = hwsm->vertexFFP();
  // call the base class
  HwME2to2Base::doinit();
}

IBPtr MEPP2GammaGammaPowhegNEW::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2GammaGammaPowhegNEW::fullclone() const {
  return new_ptr(*this);
}

void MEPP2GammaGammaPowhegNEW::persistentOutput(PersistentOStream & os) const {
  os << process_ << maxFlavour_ << photonVertex_ << gluonVertex_
     << QCDOption_ << ounit(QCDScale_,GeV2) << QEDOption_ << ounit(QEDScale_,GeV2);
}

void MEPP2GammaGammaPowhegNEW::persistentInput(PersistentIStream & is, int) {
  is >> process_ >> maxFlavour_ >> photonVertex_ >> gluonVertex_
     >> QCDOption_ >> iunit(QCDScale_,GeV2) >> QEDOption_ >> iunit(QEDScale_,GeV2);
}

int MEPP2GammaGammaPowhegNEW::nDim() const {
  return HwME2to2Base::nDim() + 3;
}

Energy2 MEPP2GammaGammaPowhegNEW::scale() const {
  Energy2 s(sHat()),u(uHat()),t(tHat());
  return 2.*s*t*u/(s*s+t*t+u*u);
}

void MEPP2GammaGammaPowhegNEW::getDiagrams() const {
  // extract the particle data objects
  tcPDPtr p = getParticleData(ParticleID::gamma);
  tcPDPtr g = getParticleData(ParticleID::g);
  vector<tcPDPtr> q(6),qb(6);
  for ( int i = 1; i <= int(maxFlavour_); ++i ) {
    q [i] = getParticleData(i);
    qb[i] = q[i]->CC();
  }
  /////  DIAGRAMS FOR PHOTON PAIRS AT LO ////////////////
  // leading-order diagrams for q qbar -> gamma gamma
  if( process_ < 5 && process_ !=3 ) {
    for ( unsigned int i = 1; i <= maxFlavour_; ++i ) {
      // t channel
      add(new_ptr((Tree2toNDiagram(3), 
		   q[i], qb[i], qb[i], 1, p, 2, p, -1)));
      // u channel
      add(new_ptr((Tree2toNDiagram(3), 
		   q[i], qb[i], qb[i], 2, p, 1, p, -2)));
    }
  }
  /////  DIAGRAMS FOR PHOTON+JET AT LO   ////////////////
  // for each quark species there are three subprocesses
  if( process_ != 2 && process_ < 4 ) {
    for ( unsigned int i=1; i<= maxFlavour_; ++i ) {
      // q qbar to gamma gluon (two diagrams)
      add(new_ptr((Tree2toNDiagram(3),
		   q[i], qb[i], qb[i], 1, p, 2, g, -11)));
      add(new_ptr((Tree2toNDiagram(3),
		   q[i],  q[i], qb[i], 2, p, 1, g, -12)));
      // q gluon to gamma q (two diagrams)
      add(new_ptr((Tree2toNDiagram(3),
		   q[i], q[i], g, 1, p, 2, q[i], -13)));
      add(new_ptr((Tree2toNDiagram(2),
		   q[i], g, 1, q[i] , 3, p, 3, q[i], -14)));
      // qbar gluon to gamma qbar (two diagrams)
      add(new_ptr((Tree2toNDiagram(3),
		   g, qb[i], qb[i], 2, p, 1, qb[i], -15)));
      add(new_ptr((Tree2toNDiagram(2),
		   g, qb[i], 1, qb[i] , 3, p, 3, qb[i], -16)));
    }
  }
  /////  DIAGRAMS FOR PHOTON PAIRS VIA BOX ////////////////
  // diagrams for g g to gamma gamma
  if( process_ == 0 || process_ == 6 ) {
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, p, 1, p, -3)));
  }
  /////  THREE-BODY DIAGRAMS               ////////////////
  if( process_ == 0 || process_ == 5 ) {
    for(unsigned int i=1;i<=maxFlavour_;++i) {
      // diagrams for q qbar -> gamma gamma gluon
      add(new_ptr((Tree2toNDiagram(4), 
		   q[i], qb[i], qb[i], qb[i], 1, p, 2, p, 3, g, -21)));
      add(new_ptr((Tree2toNDiagram(4), 
		   q[i], qb[i], qb[i], qb[i], 2, p, 1, p, 3, g, -22)));
      add(new_ptr((Tree2toNDiagram(4), 
		   q[i], qb[i], qb[i], qb[i], 1, p, 3, p, 2, g, -23)));
      add(new_ptr((Tree2toNDiagram(4), 
		   q[i], qb[i], qb[i], qb[i], 3, p, 1, p, 2, g, -24)));
      add(new_ptr((Tree2toNDiagram(4), 
		   q[i], qb[i], qb[i], qb[i], 2, p, 3, p, 1, g, -25)));
      add(new_ptr((Tree2toNDiagram(4), 
		   q[i], qb[i], qb[i], qb[i], 3, p, 2, p, 1, g, -26)));
      // diagrams for q g -> gamma gamma q
      add(new_ptr((Tree2toNDiagram(4), 
		   q[i], q[i], q[i], g, 1, p, 2, p, 3, q[i], -31)));
      add(new_ptr((Tree2toNDiagram(4), 
		   q[i], q[i], q[i], g, 2, p, 1, p, 3, q[i], -32)));
      add(new_ptr((Tree2toNDiagram(3), 
		   q[i], q[i], g, 2, q[i], 1, p, 4, p, 5, q[i], -33)));
      add(new_ptr((Tree2toNDiagram(3), 
		   q[i], q[i], g, 2, q[i], 4, p, 1, p, 5, q[i], -34)));
      add(new_ptr((Tree2toNDiagram(2), 
		   q[i], g, 1, q[i], 3, q[i], 3, p, 4, p, 4, q[i], -35)));
      add(new_ptr((Tree2toNDiagram(2), 
		   q[i], g, 1, q[i], 3, q[i], 4, p, 3, p, 4, q[i], -36)));
      // diagrams for g qbar -> gamma gamma q
      add(new_ptr((Tree2toNDiagram(4), 
		   g, qb[i], qb[i], qb[i], 2, p, 3, p, 1, qb[i], -41)));
      add(new_ptr((Tree2toNDiagram(4), 
		   g, qb[i], qb[i], qb[i], 3, p, 2, p, 1, qb[i], -42)));
      add(new_ptr((Tree2toNDiagram(3), 
		   g, qb[i], qb[i], 1, qb[i], 2, p, 4, p, 5, qb[i], -43)));
      add(new_ptr((Tree2toNDiagram(3), 
		   g, qb[i], qb[i], 1, qb[i], 4, p, 2, p, 5, qb[i], -44)));
      add(new_ptr((Tree2toNDiagram(2), 
		   g, qb[i], 1, qb[i], 3, qb[i], 3, p, 4, p, 4, qb[i], -45)));
      add(new_ptr((Tree2toNDiagram(2), 
		   g, qb[i], 1, qb[i], 3, qb[i], 4, p, 3, p, 4, qb[i], -46)));
    }
  }
}

unsigned int MEPP2GammaGammaPowhegNEW::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2GammaGammaPowhegNEW::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEPP2GammaGammaPowhegNEW::diagrams(const DiagramVector & diags) const {
  double diag1 = meInfo()[0];
  double diag2 = meInfo()[1];
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if( diags[i]->id() > -10 ) {
      if ( diags[i]->id() == -1 ) sel.insert(diag1, i);
      else if ( diags[i]->id() == -2 )  sel.insert(diag2, i);
      else if ( diags[i]->id() == -3 )  sel.insert(1.0, i);
    }
    else if( diags[i]->id() < -10 && diags[i]->id() > -20) {
      if      ( abs(diags[i]->id())%2 == 1 ) sel.insert(diag1, i);
      else                                   sel.insert(diag2, i);
    }
  }
  return sel;
}

Selector<const ColourLines *>
MEPP2GammaGammaPowhegNEW::colourGeometries(tcDiagPtr diag) const {
  // colour lines for gamma gamma pairs
  // q qbar colour lines
  static const ColourLines cqqbar("1 -2 -3");
  // g g colour lines
  static const ColourLines cgluon("1 -2,-1 2");
  // q qbar to gamma gluon colour lines
  static const ColourLines qqbar1("1 5, -5 -2 -3");
  static const ColourLines qqbar2("1 2 5, -5 -3");
  // q gluon to gamma q colour lines
  static const ColourLines qg1("1 2 -3, 3 5");
  static const ColourLines qg2("1 -2, 2 3 5");
  // qbar gluon to gamma qbar lines
  static const ColourLines qbarg1("1 -2 -3,-1 -5");
  static const ColourLines qbarg2("1 -2, -1 -3 -5");
  // selector
  Selector<const ColourLines *> sel;
  switch (diag->id()) {
  case -1: case -2:
    sel.insert(1.0, &cqqbar);
    break;
  case -3:
    sel.insert(1.0, &cgluon);
    break;
  case -11 :
    sel.insert(1.0, &qqbar1);
    break;
  case -12 :
    sel.insert(1.0, &qqbar2);
    break;
  case -13 :
    sel.insert(1.0, &qg1);
    break;
  case -14 :
    sel.insert(1.0, &qg2);
    break;
  case -15 :
    sel.insert(1.0, &qbarg1);
    break;
  case -16 :
    sel.insert(1.0, &qbarg2);
  }
  return sel;
}

bool MEPP2GammaGammaPowhegNEW::generateKinematics(const double * r) {
  if(mePartonData().size()==4) {
    // generate the ztilde variable
    int ndim=nDim();
    ztilde_ = r[ndim-1];
    // vtilde
    vtilde_ = r[ndim-2];
    // generate the azimuth
    phi_ = r[ndim-3]*Constants::twopi;
    // call base class and return
    return HwME2to2Base::generateKinematics(r);
  }
  else {
    assert(false);
  }
}

double MEPP2GammaGammaPowhegNEW::me2() const {
  // 2-2 processes
  if(mePartonData().size()==4) {
    // gamma gamma
    if(mePartonData()[2]->id()==ParticleID::gamma &&
       mePartonData()[3]->id()==ParticleID::gamma) {
      if(mePartonData()[0]->id()==-mePartonData()[1]->id()) {
	return qqbarGammaGammaME(mePartonData(),meMomenta(),true);
      }
      else if(mePartonData()[0]->id()==ParticleID::g&&
	      mePartonData()[1]->id()==ParticleID::g) {
	return ggGammaGammaME(mePartonData(),meMomenta(),true);
      }
      else assert(false);
    }
    // gamma +jet
    else {
      if(mePartonData()[0]->id()==-mePartonData()[1]->id()) {
	return qqbarGammaJetME(mePartonData(),meMomenta(),true);
      }
      else if(mePartonData()[0]->id()>0&&
	      mePartonData()[1]->id()==ParticleID::g) {
	return qgGammaJetME(mePartonData(),meMomenta(),true);
      }
      else if(mePartonData()[0]->id()==ParticleID::g&&
	      mePartonData()[1]->id()<0) {
	return gqbarGammaJetME(mePartonData(),meMomenta(),true);
      }
      else assert(false);
    }
  }
  else {
    assert(false);
  }
}

double MEPP2GammaGammaPowhegNEW::
qqbarGammaGammaME(const cPDVector & particles,
		  const vector<Lorentz5Momentum> & momenta,
		  bool first) const {
  assert(particles[0]->id()>0 && particles[0]->id()<= 5);
  assert(particles[1]->id()<0 && particles[1]->id()>=-5);
  assert(particles[2]->id()==ParticleID::gamma);
  assert(particles[3]->id()==ParticleID::gamma);
  SpinorWaveFunction    qin (momenta[0],particles[0],incoming);
  SpinorBarWaveFunction qbin(momenta[1],particles[1],incoming);
  VectorWaveFunction   p1out(momenta[2],particles[2],outgoing);
  VectorWaveFunction   p2out(momenta[3],particles[3],outgoing);
  vector<SpinorWaveFunction> fin;
  vector<SpinorBarWaveFunction>  ain;
  vector<VectorWaveFunction> p1,p2;
  for(unsigned int ix=0;ix<2;++ix) {
    qin.reset(ix)    ;fin.push_back( qin );
    qbin.reset(ix)   ;ain.push_back( qbin);
    p1out.reset(2*ix); p1.push_back(p1out);
    p2out.reset(2*ix); p2.push_back(p2out);
  }
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1,PDT::Spin1);
  // wavefunction for the intermediate particles
  SpinorWaveFunction inter;
  unsigned int inhel1,inhel2,outhel1,outhel2;
  Complex diag[3];
  double me(0.),diag1(0.),diag2(0.);
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {
	  // first diagram
	  inter = photonVertex_->evaluate(ZERO,5,fin[inhel1].particle()->CC(),
					  fin[inhel1],p1[outhel1]);
	  diag[0] = photonVertex_->evaluate(ZERO,inter,ain[inhel2],p2[outhel2]);
	  // second diagram
	  inter = photonVertex_->evaluate(ZERO,5,fin[inhel1].particle()->CC(),
					  fin[inhel1],p2[outhel2]);
	  diag[1] = photonVertex_->evaluate(ZERO,inter,ain[inhel2],p1[outhel1]);
	  // compute the running totals
	  diag[2]=diag[0]+diag[1];
	  diag1 += norm(diag[0]);
	  diag2 += norm(diag[1]);
	  me    += norm(diag[2]);
	  // matrix element
	  if(first) newme(inhel1,inhel2,2*outhel1,2*outhel2)=diag[2];
	}
      }		
    }
  }
  // save the info on the diagrams
  if(first) {
    DVector save;
    save.push_back(diag1);
    save.push_back(diag2);
    meInfo(save);
    _me.reset(newme);
  }
  // return the answer (including colour and spin factor)
  return me/12.;
}
  
double MEPP2GammaGammaPowhegNEW::
ggGammaGammaME(const cPDVector & particles,
	       const vector<Lorentz5Momentum> & momenta,
	       bool first) const {
  assert(particles[0]->id()==ParticleID::g);
  assert(particles[1]->id()==ParticleID::g);
  assert(particles[2]->id()==ParticleID::gamma);
  assert(particles[3]->id()==ParticleID::gamma);
  // get the scales
  Energy2 s((momenta[0]+momenta[1]).m2());
  Energy2 t((momenta[0]-momenta[2]).m2());
  Energy2 u((momenta[0]-momenta[3]).m2());
  Complex me[2][2][2][2];
  double charge(11./9.);
  // ++++
  me[1][1][1][1] = charge*ggme(s,t,u);
  // +++-
  me[1][1][1][0] =-charge;
  // ++-+
  me[1][1][0][1] =-charge;
  // ++--
  me[1][1][0][0] =-charge;
  // +-++
  me[1][0][1][1] =-charge;
  // +-+-
  me[1][0][1][0] = charge*ggme(u,t,s);
  // +--+
  me[1][0][0][1] = charge*ggme(t,s,u);
  // +---
  me[1][0][0][0] = charge;
  // -+++
  me[0][1][1][1] =-charge;
  // -++-
  me[0][1][1][0] =-me[1][0][0][1];
  // -+-+
  me[0][1][0][1] =-me[1][0][1][0];
  // -+--
  me[0][1][0][0] = charge;
  // --++
  me[0][0][1][1] = charge;
  // --+-
  me[0][0][1][0] = charge;
  // ---+
  me[0][0][0][1] = charge;
  // ----
  me[0][0][0][0] =-me[1][1][1][1];
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1,
				PDT::Spin1,PDT::Spin1);
  unsigned int inhel1,inhel2,outhel1,outhel2;
  double sum(0.);
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {
	  sum+=real(     me[inhel1][inhel2][outhel1][outhel2]*
			 conj(me[inhel1][inhel2][outhel1][outhel2]));
	  // matrix element
	  if(first) newme(2*inhel1,2*inhel2,
			  2*outhel1,2*outhel2)=me[inhel1][inhel2][outhel1][outhel2];
	}
      }		
    }
  }
  // final factors
  if(first) _me.reset(newme);
  return 0.5*sum*sqr(SM().alphaS(scale())*SM().alphaEM(ZERO));
}

double MEPP2GammaGammaPowhegNEW::
qqbarGammaJetME(const cPDVector & particles,
		const vector<Lorentz5Momentum> & momenta,
		bool first) const {
  assert(particles[0]->id()>0 && particles[0]->id()<= 5);
  assert(particles[1]->id()<0 && particles[1]->id()>=-5);
  assert(particles[2]->id()==ParticleID::gamma);
  assert(particles[3]->id()==ParticleID::g);
  // calculate the spinors and polarization vectors
  vector<SpinorWaveFunction> fin;
  vector<SpinorBarWaveFunction>  ain;
  vector<VectorWaveFunction> pout,gout;
  SpinorWaveFunction    qin (momenta[0],particles[0],incoming);
  SpinorBarWaveFunction qbin(momenta[1],particles[1],incoming);
  VectorWaveFunction   phout(momenta[2],particles[2],outgoing);
  VectorWaveFunction   glout(momenta[3],particles[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    qin.reset(ix)    ; fin.push_back( qin );
    qbin.reset(ix)   ; ain.push_back( qbin);
    glout.reset(2*ix);gout.push_back(glout);
    phout.reset(2*ix);pout.push_back(phout);
  }
  // calculate the matrix element
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1,PDT::Spin1);
  // wavefunction for the intermediate particles
  SpinorWaveFunction inter;
  unsigned int inhel1,inhel2,outhel1,outhel2;
  Energy2 mt(scale());
  Complex diag[3];
  double me(0.),diag1(0.),diag2(0.);
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {
	  // first diagram
	  inter = gluonVertex_->evaluate(mt,5,fin[inhel1].particle()->CC(),
					 fin[inhel1],gout[outhel1]);
	  diag[0] = photonVertex_->evaluate(ZERO,inter,ain[inhel2],pout[outhel2]);
	  // second diagram
	  inter = photonVertex_->evaluate(ZERO,5,fin[inhel1].particle()->CC(),
					  fin[inhel1],pout[outhel2]);
	  diag[1] = gluonVertex_->evaluate(mt,inter,ain[inhel2],gout[outhel1]);
	  // compute the running totals
	  diag[2]=diag[0]+diag[1];
	  diag1 +=norm(diag[0]);
	  diag2 +=norm(diag[1]);
	  me    +=norm(diag[2]);
	  // matrix element
	  if(first) newme(inhel1,inhel2,2*outhel1,2*outhel2)=diag[2];
	}
      }		
    }
  }
  // save the info on the diagrams
  if(first) {
    DVector save;
    save.push_back(diag1);
    save.push_back(diag2);
    meInfo(save);
    _me.reset(newme);
  }
  return me/9.;
}

double MEPP2GammaGammaPowhegNEW::
qgGammaJetME(const cPDVector & particles,
	     const vector<Lorentz5Momentum> & momenta,
	     bool first) const {
  assert(particles[0]->id()>0 && particles[0]->id()<= 5);
  assert(particles[1]->id()==ParticleID::g);
  assert(particles[2]->id()==ParticleID::gamma);
  assert(particles[3]->id()>0 && particles[3]->id()<= 5);
  // calculate the spinors and polarization vectors
  vector<SpinorWaveFunction> fin;
  vector<SpinorBarWaveFunction>  fout;
  vector<VectorWaveFunction> pout,gin;
  SpinorWaveFunction     qin (momenta[0],particles[0],incoming);
  VectorWaveFunction     glin(momenta[1],particles[1],incoming);  
  VectorWaveFunction    phout(momenta[2],particles[2],outgoing);
  SpinorBarWaveFunction  qout(momenta[3],particles[3],outgoing); 
  for(unsigned int ix=0;ix<2;++ix) {
    qin.reset(ix)    ;fin.push_back(  qin );
    qout.reset(ix)   ;fout.push_back( qout);
    glin.reset(2*ix) ;gin.push_back(  glin);
    phout.reset(2*ix);pout.push_back(phout);
  }
  // calculate the matrix element
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1,
				PDT::Spin1,PDT::Spin1Half);
  // wavefunction for the intermediate particles
  SpinorWaveFunction inter;
  unsigned int inhel1,inhel2,outhel1,outhel2;
  Energy2 mt(scale());
  Complex diag[3];
  double me(0.),diag1(0.),diag2(0.);
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {
	  // first diagram
	  inter = photonVertex_->evaluate(ZERO,5,fin[inhel1].particle()->CC(),
					  fin[inhel1],pout[outhel1]);
	  diag[0]=gluonVertex_->evaluate(mt,inter,fout[outhel2],gin[inhel2]);
	  // second diagram
	  inter = gluonVertex_->evaluate(mt,5,fin[inhel1].particle()->CC(),
					 fin[inhel1],gin[inhel2]);
	  diag[1]=photonVertex_->evaluate(ZERO,inter,fout[outhel2],pout[outhel1]);
	  // compute the running totals
	  diag[2]=diag[0]+diag[1];
	  diag1 +=norm(diag[0]);
	  diag2 +=norm(diag[1]);
	  me    +=norm(diag[2]);
	  // matrix element
	  if(first) newme(inhel1,2*inhel2,2*outhel1,outhel2)=diag[2];
	}
      }		
    }
  }
  // save the info on the diagrams
  if(first) {
    DVector save;
    save.push_back(diag1);
    save.push_back(diag2);
    meInfo(save); 
    _me.reset(newme);
  }
  return me/24.;
}

double MEPP2GammaGammaPowhegNEW::
gqbarGammaJetME(const cPDVector & particles,
		const vector<Lorentz5Momentum> & momenta,
		bool first) const {
  assert(particles[0]->id()==ParticleID::g);
  assert(particles[1]->id()<0 && particles[0]->id()>=-5);
  assert(particles[2]->id()==ParticleID::gamma);
  assert(particles[3]->id()<0 && particles[3]->id()>=-5);
  // calculate the spinors and polarization vectors
  vector<SpinorBarWaveFunction> ain;
  vector<SpinorWaveFunction>  aout;
  vector<VectorWaveFunction> pout,gin;
  VectorWaveFunction     glin(momenta[0],particles[0],incoming);  
  SpinorBarWaveFunction  qin (momenta[1],particles[1],incoming);
  VectorWaveFunction    phout(momenta[2],particles[2],outgoing); 
  SpinorWaveFunction     qout(momenta[3],particles[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    qin.reset(ix)    ;ain.push_back(  qin );
    qout.reset(ix)   ;aout.push_back( qout);
    glin.reset(2*ix) ;gin.push_back(  glin);
    phout.reset(2*ix);pout.push_back(phout);
  }
  // calculate the matrix element
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1,
				PDT::Spin1,PDT::Spin1Half);
  // wavefunction for the intermediate particles
  SpinorBarWaveFunction inter;
  SpinorWaveFunction interb;
  unsigned int inhel1,inhel2,outhel1,outhel2;
  Energy2 mt(scale());
  Complex diag[3];
  double me(0.),diag1(0.),diag2(0.);
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {
	  // first diagram
	  inter = photonVertex_->evaluate(ZERO,5,ain[inhel1].particle()->CC(),
					  ain[inhel1],pout[outhel1]);
	  diag[0]=gluonVertex_->evaluate(mt,aout[outhel2],inter,gin[inhel2]);
	  // second diagram
	  inter = gluonVertex_->evaluate(mt,5,ain[inhel1].particle()->CC(),
					 ain[inhel1],gin[inhel2]);
	  diag[1]=photonVertex_->evaluate(ZERO,aout[outhel2],inter,pout[outhel1]);
	  // compute the running totals
	  diag[2]=diag[0]+diag[1];
	  diag1 +=norm(diag[0]);
	  diag2 +=norm(diag[1]);
	  me    +=norm(diag[2]);
	  // matrix element
	  if(first) newme(inhel1,2*inhel2,2*outhel1,outhel2)=diag[2];
	}
      }		
    }
  }
  // save the info on the diagrams
  if(first) {
    DVector save;
    save.push_back(diag1);
    save.push_back(diag2);
    meInfo(save);
    _me.reset(newme);
  }
  return me/24.;
}

InvEnergy2 MEPP2GammaGammaPowhegNEW::
qqbargME(const vector<tcPDPtr> & particles,
	 const vector<Lorentz5Momentum> & p) const {
  double sum(0.);
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qbarin;
  vector<VectorWaveFunction> pout1,pout2,gout;
  SpinorWaveFunction    q_in   (p[0],particles[0],incoming);
  SpinorBarWaveFunction qbar_in(p[1],particles[1],incoming);
  VectorWaveFunction     p_out1(p[2],particles[2],outgoing);
  VectorWaveFunction     p_out2(p[3],particles[3],outgoing);
  VectorWaveFunction     g_out (p[4],particles[4],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q_in.reset(ix);               qin.push_back(q_in   );
    qbar_in.reset(ix);         qbarin.push_back(qbar_in);
    g_out.reset(2*ix);           gout.push_back(g_out  );
    p_out1.reset(2*ix);         pout1.push_back(p_out1 );
    p_out2.reset(2*ix);         pout2.push_back(p_out2 );
  }
  vector<Complex> diag(6);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int phel1=0;phel1<2;++phel1) {
	for(unsigned int phel2=0;phel2<2;++phel2) {
	  for(unsigned int ghel=0;ghel<2;++ghel) {
 	    // first diagram
	    SpinorWaveFunction inters1 = 
	      photonVertex_->evaluate(ZERO,5,mePartonData()[0],
				      qin[ihel1],pout1[phel1]);
	    SpinorBarWaveFunction inters2 = 
	      photonVertex_->evaluate(ZERO,5,mePartonData()[1],
				      qbarin[ihel2],pout2[phel2]);
	    diag[0] = gluonVertex_->evaluate(scale(),inters1,inters2,gout[ghel]);
	    // second diagram
	    SpinorWaveFunction inters3 = 
	      gluonVertex_->evaluate(scale(),5,mePartonData()[0],qin[ihel1],gout[ghel]);
	    SpinorBarWaveFunction inters4 = 
	      photonVertex_->evaluate(ZERO,5,mePartonData()[1],
				      qbarin[ihel2],pout1[phel1]);
	    diag[1] = photonVertex_->evaluate(ZERO,inters3,inters4,pout2[phel2]);
	    // fourth diagram
	    diag[2] = photonVertex_->evaluate(ZERO,inters3,inters2,pout1[phel1]);
	    // fifth diagram
	    SpinorBarWaveFunction inters5 = 
	      gluonVertex_->evaluate(scale(),5,mePartonData()[1],
				     qbarin[ihel2],gout[ghel]);
	    diag[3] = 
	      photonVertex_->evaluate(ZERO,inters1,inters5,pout2[phel2]);
	    // sixth diagram
	    SpinorWaveFunction inters6 = 
	      photonVertex_->evaluate(ZERO,5,mePartonData()[0],qin[ihel1],pout2[phel2]);
	    diag[4] = gluonVertex_->evaluate(scale(),inters6,inters4,gout[ghel]);
	    // eighth diagram
	    diag[5] = photonVertex_->evaluate(ZERO,inters6,inters5,pout1[phel1]);
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  return sum/9.*UnitRemoval::InvE2;
}

InvEnergy2 MEPP2GammaGammaPowhegNEW::
qgqME(const vector<tcPDPtr> & particles,
      const vector<Lorentz5Momentum> & p) const {
  double sum(0.);
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qout;
  vector<VectorWaveFunction> pout1,pout2,gin;
  SpinorWaveFunction q_in    (p[0],particles[0],incoming);
  VectorWaveFunction g_in    (p[1],particles[1],incoming);
  VectorWaveFunction p_out1  (p[2],particles[2],outgoing);
  VectorWaveFunction p_out2  (p[3],particles[3],outgoing);
  SpinorBarWaveFunction q_out(p[4],particles[4],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q_in .reset(ix);         qin  .push_back(q_in  );
    q_out.reset(ix);         qout .push_back(q_out );
    g_in .reset(2*ix);       gin  .push_back(g_in  );
    p_out1.reset(2*ix);      pout1.push_back(p_out1);
    p_out2.reset(2*ix);      pout2.push_back(p_out2);
  }
  vector<Complex> diag(6);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int phel1=0;phel1<2;++phel1) {
	for(unsigned int phel2=0;phel2<2;++phel2) {
	  for(unsigned int ohel=0;ohel<2;++ohel) {
 	    // first diagram
	    SpinorWaveFunction inters1 = 
	      photonVertex_->evaluate(ZERO,5,qin[ihel1].particle()->CC(),
				      qin[ihel1],pout1[phel1]);
	    SpinorBarWaveFunction inters2 = 
	      photonVertex_->evaluate(ZERO,5,qout[ohel].particle(),
				      qout[ohel],pout2[phel2]);
	    diag[0] = gluonVertex_->evaluate(scale(),inters1,inters2,gin[ihel2]);
	    // second diagram
	    SpinorWaveFunction inters3 = 
	      gluonVertex_->evaluate(scale(),5,qin[ihel1].particle()->CC(),
				     qin[ihel1],gin[ihel2]);
	    SpinorBarWaveFunction inters4 = 
	      photonVertex_->evaluate(ZERO,5,qout[ohel].particle(),
				      qout[ohel],pout1[phel1]);
	    diag[1] = photonVertex_->evaluate(ZERO,inters3,inters4,pout2[phel2]);
	    // fourth diagram
	    diag[2] = photonVertex_->evaluate(ZERO,inters3,inters2,pout1[phel1]);
	    // fifth diagram
	    SpinorBarWaveFunction inters5 = 
	      gluonVertex_->evaluate(scale(),5,qout[ohel].particle(),
				     qout[ohel],gin[ihel2]);
	    diag[3] = 
	      photonVertex_->evaluate(ZERO,inters1,inters5,pout2[phel2]);
	    // sixth diagram
	    SpinorWaveFunction inters6 = 
	      photonVertex_->evaluate(ZERO,5,qin[ihel1].particle()->CC(),
				      qin[ihel1],pout2[phel2]);
	    diag[4] = gluonVertex_->evaluate(scale(),inters6,inters4,gin[ihel2]);
	    // eighth diagram
	    diag[5] = photonVertex_->evaluate(ZERO,inters6,inters5,pout1[phel1]);
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  return sum/24.*UnitRemoval::InvE2;
}

InvEnergy2 MEPP2GammaGammaPowhegNEW::
qbargqbarME(const vector<tcPDPtr> & particles,
	    const vector<Lorentz5Momentum> & p) const {
  double sum(0.);
  vector<SpinorBarWaveFunction> qin;
  vector<SpinorWaveFunction> qout;
  vector<VectorWaveFunction> pout1,pout2,gin;
  VectorWaveFunction g_in   (p[0],particles[0],incoming);
  SpinorBarWaveFunction q_in(p[1],particles[1],incoming);
  VectorWaveFunction p_out1 (p[2],particles[2],outgoing);
  VectorWaveFunction p_out2 (p[3],particles[3],outgoing);
  SpinorWaveFunction q_out  (p[4],particles[4],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q_in .reset(ix);         qin  .push_back(q_in  );
    q_out.reset(ix);         qout .push_back(q_out );
    g_in .reset(2*ix);       gin  .push_back(g_in  );
    p_out1.reset(2*ix);      pout1.push_back(p_out1);
    p_out2.reset(2*ix);      pout2.push_back(p_out2);
  }
  vector<Complex> diag(6);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int phel1=0;phel1<2;++phel1) {
	for(unsigned int phel2=0;phel2<2;++phel2) {
	  for(unsigned int ohel=0;ohel<2;++ohel) {
 	    // first diagram
	    SpinorBarWaveFunction inters1 = 
	      photonVertex_->evaluate(ZERO,5,qin[ihel1].particle()->CC(),
				      qin[ihel1],pout1[phel1]);
	    SpinorWaveFunction inters2 = 
	      photonVertex_->evaluate(ZERO,5,qout[ohel].particle(),
				      qout[ohel],pout2[phel2]);
	    diag[0] = gluonVertex_->evaluate(scale(),inters2,inters1,gin[ihel2]);
	    // second diagram
	    SpinorBarWaveFunction inters3 = 
	      gluonVertex_->evaluate(scale(),5,qin[ihel1].particle()->CC(),
				     qin[ihel1],gin[ihel2]);
	    SpinorWaveFunction inters4 = 
	      photonVertex_->evaluate(ZERO,5,qout[ohel].particle(),
				      qout[ohel],pout1[phel1]);
	    diag[1] = photonVertex_->evaluate(ZERO,inters4,inters3,pout2[phel2]);
	    // fourth diagram
	    diag[2] = photonVertex_->evaluate(ZERO,inters2,inters3,pout1[phel1]);
	    // fifth diagram
	    SpinorWaveFunction inters5 = 
	      gluonVertex_->evaluate(scale(),5,qout[ohel].particle(),
				     qout[ohel],gin[ihel2]);
	    diag[3] = 
	      photonVertex_->evaluate(ZERO,inters5,inters1,pout2[phel2]);
	    // sixth diagram
	    SpinorBarWaveFunction inters6 = 
	      photonVertex_->evaluate(ZERO,5,qin[ihel1].particle()->CC(),
				      qin[ihel1],pout2[phel2]);
	    diag[4] = gluonVertex_->evaluate(scale(),inters4,inters6,gin[ihel2]);
	    // eighth diagram
	    diag[5] = photonVertex_->evaluate(ZERO,inters5,inters6,pout1[phel1]);
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  return sum/24.*UnitRemoval::InvE2;
}

HardTreePtr MEPP2GammaGammaPowhegNEW::generateHardest(ShowerTreePtr) const {
  cerr << "testing got to the ME class\n";
  exit(0);
}
