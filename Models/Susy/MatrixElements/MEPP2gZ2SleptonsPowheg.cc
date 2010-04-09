// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2gZ2SleptonsPowheg class.
//

#include "MEPP2gZ2SleptonsPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/Susy/SusyBase.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;


MEPP2gZ2SleptonsPowheg::MEPP2gZ2SleptonsPowheg() : process_(0) {
  vector<unsigned int> mopt(2,1);
  massOption(mopt);
}

void MEPP2gZ2SleptonsPowheg::doinit() {
  NLODrellYanBase::doinit();
  // get the photon and Z ParticleData objects
  Z0_    = getParticleData(ThePEG::ParticleID::Z0);
  gamma_ = getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcSusyBasePtr hwsm=ThePEG::dynamic_ptr_cast<tcSusyBasePtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    FFZVertex_ = hwsm->vertexFFZ();
    FFPVertex_ = hwsm->vertexFFP();
    FFGVertex_ = hwsm->vertexFFG();
    WSSVertex_ = hwsm->vertexWSFSF();
  }
  else
    throw InitException() << "Must be the Herwig++ SusyBase class in "
			  << "MEPP2gZ2SleptonsPowheg::doinit" 
			  << Exception::abortnow;
}

Selector<const ColourLines *>
MEPP2gZ2SleptonsPowheg::colourGeometries(tcDiagPtr) const {
  static const ColourLines c1("1 -2");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c1);
  return sel;
}

void MEPP2gZ2SleptonsPowheg::getDiagrams() const {
  // loop over the processes we need
  for(int i = 1; i <= maxFlavour_; ++i) {
    tcPDPtr q  = getParticleData(i);
    tcPDPtr qb = q->CC();
    for(int ix=11;ix<17;++ix) {
      // production of left-handed sleptons
      tcPDPtr lm = getParticleData(1000000+ix);
      tcPDPtr lp = lm->CC();
      // always Z
      if(process_ == 0 || (ix%2==1 && process_ == (ix-9)/2) ||
	 (ix%2==0 && process_ == 6 +(ix-9)/2) ) 
	add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Z0_   , 3, lm, 3, lp, -1)));
      // if sneutrinos that's all
      if(ix%2==0) continue; 
      // photon 
      if(process_ == 0 || process_ == (ix-9)/2 )
	add(new_ptr((Tree2toNDiagram(2), q, qb, 1, gamma_, 3, lm, 3, lp, -2)));
      // production of right-handed sleptons
      tcPDPtr rm = getParticleData(2000000+ix);
      tcPDPtr rp = rm->CC();
      if(process_ == 0 || process_ == 3+(ix-9)/2 ) {
	add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Z0_   , 3, rm, 3, rp, -1)));
	add(new_ptr((Tree2toNDiagram(2), q, qb, 1, gamma_, 3, rm, 3, rp, -2)));
      }
      // production of left-right for stau only
      if(ix==15 && (process_==0 || process_==10)) {
	add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Z0_   , 3, rm, 3, lp, -1)));
	add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Z0_   , 3, lm, 3, rp, -1)));
      }
    }
  }
}

unsigned int MEPP2gZ2SleptonsPowheg::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2gZ2SleptonsPowheg::orderInAlphaEW() const {
  return 2;
}

IBPtr MEPP2gZ2SleptonsPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2gZ2SleptonsPowheg::fullclone() const {
  return new_ptr(*this);
}

void MEPP2gZ2SleptonsPowheg::persistentOutput(PersistentOStream & os) const {
  os << FFZVertex_ << FFPVertex_ << WSSVertex_ << FFGVertex_ 
     << Z0_ << gamma_ << process_ << maxFlavour_;
}

void MEPP2gZ2SleptonsPowheg::persistentInput(PersistentIStream & is, int) {
  is >> FFZVertex_ >> FFPVertex_ >> WSSVertex_ >> FFGVertex_ 
     >> Z0_ >> gamma_ >> process_ >> maxFlavour_;
}

ClassDescription<MEPP2gZ2SleptonsPowheg> 
MEPP2gZ2SleptonsPowheg::initMEPP2gZ2SleptonsPowheg;
// Definition of the static class description member.

void MEPP2gZ2SleptonsPowheg::Init() {

  static ClassDocumentation<MEPP2gZ2SleptonsPowheg> documentation
    ("MEPP2gZ2SleptonsPowheg implements the ME calculation"
     " of the fermion-antifermion to sfermion-sfermion hard process.");


  static Switch<MEPP2gZ2SleptonsPowheg,int> interfaceProcess
    ("Process",
     "Which processes to generate",
     &MEPP2gZ2SleptonsPowheg::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Generate all the processes",
     0);
  static SwitchOption interfaceProcesse_L
    (interfaceProcess,
     "e_L",
     "Only produce ~e_L",
     1);
  static SwitchOption interfaceProcessmu_L
    (interfaceProcess,
     "mu_L",
     "Onle produce ~mu_L",
     2);
  static SwitchOption interfaceProcesstau_1
    (interfaceProcess,
     "tau_1",
     "Only produce tau_1 pairs",
     3);
  static SwitchOption interfaceProcesse_R
    (interfaceProcess,
     "e_R",
     "Only produce e_R",
     4);
  static SwitchOption interfaceProcessmu_R
    (interfaceProcess,
     "mu_R",
     "Only produce ~mu_R",
     5);
  static SwitchOption interfaceProcesstau_2
    (interfaceProcess,
     "tau_2",
     "Only produce tau_2 pairs",
     6);
  static SwitchOption interfaceProcessnu_e
    (interfaceProcess,
     "nu_e",
     "Only product ~nu_e",
     7);
  static SwitchOption interfaceProcessnu_mu
    (interfaceProcess,
     "nu_mu",
     "Only produce ~nu_mu",
     8);
  static SwitchOption interfaceProcessnu_tau
    (interfaceProcess,
     "nu_tau",
     "Only produce ~nu_tau",
     9);
  static SwitchOption interfaceProcessMixedTau
    (interfaceProcess,
     "MixedTau",
     "Only produce mixing tau_1 tau_2 pairs",
     10);

  static Parameter<MEPP2gZ2SleptonsPowheg,int> interfaceMaxFlavour
    ("MaxFlavour",
     "The maximum flavour of the incoming quarks",
     &MEPP2gZ2SleptonsPowheg::maxFlavour_, 5, 1, 5,
     false, false, Interface::limited);

}

Selector<MEBase::DiagramIndex>
MEPP2gZ2SleptonsPowheg::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
  }
  return sel;
}

NLODrellYanBase::Singular MEPP2gZ2SleptonsPowheg::virtualME() const {
  Singular output;
  output.eps2 = -2;
  output.eps1 = -3;
  output.finite =-8.+sqr(Constants::pi);
  return output;
}

double MEPP2gZ2SleptonsPowheg::
qqbarME(vector<SpinorWaveFunction>    & sp ,
	vector<SpinorBarWaveFunction> & sbar ,
	ScalarWaveFunction & sca1,ScalarWaveFunction &sca2,
	bool first) const {
  // scale for the process
  const Energy2 q2(scale());
  // storage of the matrix elements for specific diagrams
  vector<double> me(2, 0.);
  double me2(0.);
  // storgage of the individual diagrams
  vector<Complex> diag(2, Complex(0.));
  ProductionMatrixElement pme(PDT::Spin1Half, PDT::Spin1Half, 
			      PDT::Spin0, PDT::Spin0);
  // loop over the helicities and calculate the matrix elements
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      VectorWaveFunction interV = FFZVertex_->evaluate(q2, 1, Z0_, sp[if1], 
						       sbar[if2]);
      diag[0] = WSSVertex_->evaluate(q2, interV, sca2, sca1);
      // if needed photon diagram
      if(mePartonData()[2]->charged() && 
	 mePartonData()[2]->id() == -mePartonData()[3]->id()) {
	interV = FFPVertex_->evaluate(q2, 1, gamma_, sp[if1], 
				      sbar[if2]);
	diag[1] = WSSVertex_->evaluate(q2, interV, sca2, sca1);
      }
      // sum up the matrix elements
      me2 += norm(diag[0]+diag[1]);
      me[0] += norm(diag[0]);
      me[1] += norm(diag[1]);
      pme(if1, if2, 0, 0) = diag[0]+diag[1];
    }
  }
  if(first) {
    DVector save(2);
    for(DVector::size_type ix = 0; ix < 2; ++ix)
      save[ix] = me[ix]/12.;
    meInfo(save);
    me_.reset(pme);
  }
  return me2/12.;
}

double MEPP2gZ2SleptonsPowheg::loME(const cPDVector & particles,
				    const vector<Lorentz5Momentum> & momenta,
				    bool first) const {
  // wavefunctions for the incoming fermions
  vector<SpinorWaveFunction> sp(2);
  vector<SpinorBarWaveFunction> sbar(2);
  for( unsigned int i = 0; i < 2; ++i ) {
    sp[i] = SpinorWaveFunction(momenta[0], particles[0], i,
			       incoming);
    sbar[i] = SpinorBarWaveFunction(momenta[1], particles[1], i,
				    incoming);
  }
  // outgoing scalar wavefunctions
  ScalarWaveFunction sca1(momenta[2], particles[2],
			  Complex(1.), outgoing);
  ScalarWaveFunction sca2(momenta[3], particles[3],
			  Complex(1.), outgoing);
  return qqbarME(sp,sbar,sca1,sca2,first);
}

double MEPP2gZ2SleptonsPowheg::realME(const cPDVector & particles,
				      const vector<Lorentz5Momentum> & momenta) const {
  vector<SpinorWaveFunction> sp(2);
  vector<SpinorBarWaveFunction> sbar(2);
  vector<VectorWaveFunction> gluon(2);
  // wavefunctions for the q qbar -> sf sf g process
  if(particles[0]->id()==-particles[1]->id()) {
    for( unsigned int i = 0; i < 2; ++i ) {
      sp[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      sbar[i] = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
      gluon[i]= VectorWaveFunction   (momenta[4],particles[4],2*i,outgoing);
    }
  }
  else if(particles[0]->id()==ParticleID::g &&
	  particles[1]->id()<0) {
    for( unsigned int i = 0; i < 2; ++i ) {
      sp[i]   = SpinorWaveFunction   (momenta[4],particles[4],  i,outgoing);
      sbar[i] = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
      gluon[i]= VectorWaveFunction   (momenta[0],particles[0],2*i,incoming);
    }
  }
  else if(particles[0]->id()>0 &&
	  particles[1]->id()==ParticleID::g) {
    for( unsigned int i = 0; i < 2; ++i ) {
      sp[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      sbar[i] = SpinorBarWaveFunction(momenta[4],particles[4],  i,outgoing);
      gluon[i]= VectorWaveFunction   (momenta[1],particles[1],2*i,incoming);
    }
  }
  else {
    for(unsigned int ix=0;ix<particles.size();++ix) {
      cerr << particles[ix]->PDGName() << " " << momenta[ix]/GeV << "\n";
    }
    assert(false);
  }
  // wavefunctions for the scalars 
  ScalarWaveFunction sca1(momenta[2], particles[2],Complex(1.), outgoing);
  ScalarWaveFunction sca2(momenta[3], particles[3],Complex(1.), outgoing);
  double output(0.);
  Complex diag[4]={0.,0.,0.,0.};
  Energy2 shat = scale();
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	// first Z diagram
 	SpinorWaveFunction inters = FFGVertex_->evaluate(shat,5,sp[ihel1].particle(),
							 sp[ihel1],gluon[ohel1]);
	VectorWaveFunction interV = FFZVertex_->evaluate(shat, 1, Z0_, inters, 
							 sbar[ihel2]);
	diag[0] = WSSVertex_->evaluate(shat, interV, sca2, sca1);
	// second Z diagram
	SpinorBarWaveFunction interb = FFGVertex_->evaluate(shat,5,sbar[ihel1].particle(),
							    sbar[ihel2],gluon[ohel1]);
	interV = FFZVertex_->evaluate(shat, 1, Z0_, sp[ihel1], 
				      interb);
	diag[1] = WSSVertex_->evaluate(shat, interV, sca2, sca1);
	if(particles[2]->id()==-particles[3]->id()&&particles[2]->charged()) {
	  // first photon diagram
	  SpinorWaveFunction inters = FFGVertex_->evaluate(shat,5,sp[ihel1].particle(),
							   sp[ihel1],gluon[ohel1]);
	  VectorWaveFunction interV = FFPVertex_->evaluate(shat, 1, gamma_, inters, 
							   sbar[ihel2]);
	  diag[2] = WSSVertex_->evaluate(shat, interV, sca2, sca1);
	  // second photon diagram
	  SpinorBarWaveFunction interb = FFGVertex_->evaluate(shat,5,sbar[ihel1].particle(),
							      sbar[ihel2],gluon[ohel1]);
	  interV = FFPVertex_->evaluate(shat, 1, gamma_, sp[ihel1], 
					interb);
	  diag[3] = WSSVertex_->evaluate(shat, interV, sca2, sca1);
	}
	// add them up
	output += norm(diag[0]+diag[1]+diag[2]+diag[3]);
      }
    }
  }
  // colour and spin factors
  if(particles[0]->id()==-particles[1]->id()) {
    output *= 1./9.;
  }
  else  {
    output *= 1./24.;
  }
  // divided by 2 g_S^2
  return 0.5*output/norm(FFGVertex_->norm());
}

void MEPP2gZ2SleptonsPowheg::constructVertex(tSubProPtr sub) {
  //get particles
  ParticleVector ext(4);
  ext[0] = sub->incoming().first;
  ext[1] = sub->incoming().second;
  ext[2] = sub->outgoing()[0];
  ext[3] = sub->outgoing()[1];
  if( ext[0]->id() != mePartonData()[0]->id() ) swap(ext[0], ext[1]);
  if( ext[2]->id() != mePartonData()[2]->id() ) swap(ext[2], ext[3]);
  //First calculate wave functions with off-shell momenta
  //to calculate correct spin information
  vector<SpinorWaveFunction> sp;
  SpinorWaveFunction(sp, ext[0], incoming, false);
  vector<SpinorBarWaveFunction> sbar;
  SpinorBarWaveFunction(sbar, ext[1], incoming, false);
  ScalarWaveFunction sca1(ext[2], outgoing, true);
  ScalarWaveFunction sca2(ext[3], outgoing, true);
  //Need to use rescale momenta to calculate matrix element
  cPDVector data(4);
  vector<Lorentz5Momentum> momenta(4);
  for( size_t i = 0; i < 4; ++i ) {
    data[i] = ext[i]->dataPtr();
    momenta[i] = ext[i]->momentum();
  }
  rescaleMomenta(momenta, data);
  SpinorWaveFunction spr(rescaledMomenta()[0], data[0], incoming);
  SpinorBarWaveFunction sbr(rescaledMomenta()[1], data[1],incoming);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
    spr.reset(ihel);
    sp[ihel] = spr;
    sbr.reset(ihel);
    sbar[ihel] = sbr;
  }
  sca1 = ScalarWaveFunction(rescaledMomenta()[2], data[2], outgoing);
  sca2 = ScalarWaveFunction(rescaledMomenta()[3], data[3], outgoing);
  qqbarME(sp, sbar, sca1, sca2,true);
  HardVertexPtr hv = new_ptr(HardVertex());
  hv->ME(me_);
  for(unsigned int i = 0; i < 4; ++i )
    dynamic_ptr_cast<SpinfoPtr>(ext[i]->spinInfo())->setProductionVertex(hv);  
}





