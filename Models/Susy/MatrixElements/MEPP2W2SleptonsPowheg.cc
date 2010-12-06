// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2W2SleptonsPowheg class.
//

#include "MEPP2W2SleptonsPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/Susy/SusyBase.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include "SusyLoopIntegral.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;


MEPP2W2SleptonsPowheg::MEPP2W2SleptonsPowheg() : process_(0), maxFlavour_(4) {
  vector<unsigned int> mopt(2,1);
  massOption(mopt);
}

void MEPP2W2SleptonsPowheg::doinit() {
  NLODrellYanBase::doinit();
  // get the photon and Z ParticleData objects
  Wplus_    = getParticleData(ThePEG::ParticleID::Wplus);
  Wminus_   = getParticleData(ThePEG::ParticleID::Wminus);
  // cast the SM pointer to the Herwig SM pointer
  tcSusyBasePtr hwsm=ThePEG::dynamic_ptr_cast<tcSusyBasePtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    FFWVertex_ = hwsm->vertexFFW();
    FFGVertex_ = hwsm->vertexFFG();
    WSSVertex_ = hwsm->vertexWSFSF();
  }
  else
    throw InitException() << "Must be the Herwig++ SusyBase class in "
			  << "MEPP2W2SleptonsPowheg::doinit" 
			  << Exception::abortnow;
}

Selector<const ColourLines *>
MEPP2W2SleptonsPowheg::colourGeometries(tcDiagPtr) const {
  static const ColourLines c1("1 -2");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c1);
  return sel;
}

void MEPP2W2SleptonsPowheg::getDiagrams() const {
  // loop over the Wplus processes we need
  for(int i = 2; i <= maxFlavour_; i+=2 ) {
    tcPDPtr qi  = getParticleData(i);
    tcPDPtr qib = qi->CC();
    for(int j = 1; j <= maxFlavour_; j+=2 ) {
      tcPDPtr qj  = getParticleData(j);
      tcPDPtr qjb = qj->CC();
      
      tcPDPtr q  = qi ;
      tcPDPtr qb = qjb;
      
      for(int ix=11; ix<17; ix+=2 ) {
	// production of left-handed sleptons
	tcPDPtr lp = getParticleData(-(1000000+ix));
	tcPDPtr lnu = getParticleData(1000000+ix+1);
	if(process_==0|| int(process_)==4+(ix-9)/2 )
	  add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Wplus_,
		       3, lnu, 3, lp, -1)));
	// production of stau_2
	if(ix==15) {
	  lp = getParticleData(-(2000000+ix));
	  if(process_==0||process_== 8 )
	    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Wplus_,
			 3, lnu, 3, lp, -1)));
	}
      }
    }
  }
  // loop over the Wminus processes we need
  for(int i = 2; i <= 5; i+=2) {
    tcPDPtr qi  = getParticleData(i);
    tcPDPtr qib = qi->CC();
    
    for(int j = 1; j <= 5; j+=2) {
      tcPDPtr qj  = getParticleData(j);
      tcPDPtr qjb = qj->CC();
    
      tcPDPtr q=qj;
      tcPDPtr qb=qib;
    
      for(int ix=11; ix<17; ix+=2) {
	// production of left-handed sleptons
	tcPDPtr lp = getParticleData(1000000+ix);
	tcPDPtr lnu = getParticleData(-(1000000+ix+1));
	if(process_==0|| process_==(ix-9)/2 )
	  add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Wminus_,
		       3, lnu, 3, lp, -2)));
	if(ix==15) {
	  lp = getParticleData(2000000+ix);
	  if(process_==0||process_==4 )
	    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Wminus_,
			 3, lnu, 3, lp, -2)));
	}  
      }
    }
  } 
}

unsigned int MEPP2W2SleptonsPowheg::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2W2SleptonsPowheg::orderInAlphaEW() const {
  return 2;
}

IBPtr MEPP2W2SleptonsPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2W2SleptonsPowheg::fullclone() const {
  return new_ptr(*this);
}

void MEPP2W2SleptonsPowheg::persistentOutput(PersistentOStream & os) const {
  os << FFWVertex_ << WSSVertex_ << FFGVertex_ 
     << Wplus_ << Wminus_ << process_ << maxFlavour_;
}

void MEPP2W2SleptonsPowheg::persistentInput(PersistentIStream & is, int) {
  is >> FFWVertex_ >> WSSVertex_ >> FFGVertex_ 
     >> Wplus_ >> Wminus_ >> process_ >> maxFlavour_;
}

ClassDescription<MEPP2W2SleptonsPowheg> 
MEPP2W2SleptonsPowheg::initMEPP2W2SleptonsPowheg;
// Definition of the static class description member.

void MEPP2W2SleptonsPowheg::Init() {

  static ClassDocumentation<MEPP2W2SleptonsPowheg> documentation
    ("MEPP2W2SleptonsPowheg implements the ME calculation"
     " of the fermion-antifermion to sfermion-sfermion hard process.");


  static Switch<MEPP2W2SleptonsPowheg,int> interfaceProcess
    ("Process",
     "Which processes to generate",
     &MEPP2W2SleptonsPowheg::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Generate all processes",
     0);
  static SwitchOption interfaceProcesse_Lminus
    (interfaceProcess,
     "e_Lminus",
     "Produce e_L-",
     1);
  static SwitchOption interfaceProcessmu_Lminus
    (interfaceProcess,
     "mu_Lminus",
     "Produce mu_L-",
     2);
  static SwitchOption interfaceProcesstau_1minus
    (interfaceProcess,
     "tau_1minus",
     "Produce tau_1-",
     3);
  static SwitchOption interfaceProcesstau_2minus
    (interfaceProcess,
     "tau_2minus",
     "Produce tau_2-",
     4);
  static SwitchOption interfaceProcesse_Lplus
    (interfaceProcess,
     "e_Lplus",
     "Produce e_L-",
     5);
  static SwitchOption interfaceProcessmu_Lplus
    (interfaceProcess,
     "mu_Lplus",
     "Produce mu_L-",
     6);
  static SwitchOption interfaceProcesstau_1plus
    (interfaceProcess,
     "tau_1plus",
     "Produce tau_1-",
     7);
  static SwitchOption interfaceProcesstau_2plus
    (interfaceProcess,
     "tau_2plus",
     "Produce tau_2-",
     8);

  static Parameter<MEPP2W2SleptonsPowheg,int> interfaceMaxFlavour
    ("MaxFlavour",
     "The maximum flavour of the incoming quarks",
     &MEPP2W2SleptonsPowheg::maxFlavour_, 5, 1, 5,
     false, false, Interface::limited);
}

Selector<MEBase::DiagramIndex>
MEPP2W2SleptonsPowheg::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    sel.insert(1., i);
  }
  return sel;
}

NLODrellYanBase::Singular MEPP2W2SleptonsPowheg::virtualME() const {
  // singular pieces of the output
  Singular output;
  output.eps2 = -2;
  output.eps1 = -3;
  // cut-off parameter
  Energy2 eps = 0.1*MeV2;
  // average of outgoing masses
  Energy2 mav2 = 0.25*sqr(meMomenta()[2].mass()+meMomenta()[3].mass());
  // renormalisation/factorization scale
  Energy2 scale2 = scale();
  // gluino mass
  Energy mg  = getParticleData(ParticleID::SUSY_g)->mass();
  Energy2 mg2(sqr(mg));
  // average of left/right squark masses
  Energy  ms = 0.5*(getParticleData(1000000+mePartonData()[0]->id())->mass()+
		    getParticleData(2000000+mePartonData()[0]->id())->mass());
  Energy2 ms2(sqr(ms));
  // the finite piece
  output.finite = -3.*log(sHat()/mav2) + sqr(log(sHat()/mav2))
    - sqr(log(sHat()/scale2)) + sqr(Constants::pi) -1.
    +(1.+2.*(mg2-ms2)/sHat())*(SusyLoopIntegral::B0 (sHat(),ms    ,ms    ,scale2)-
			       SusyLoopIntegral::B0 (ZERO  ,mg    ,ms    ,scale2))
    -3.*SusyLoopIntegral::B0 (sHat(),ZERO  ,ZERO  ,scale2)
    -   SusyLoopIntegral::B0P(ZERO  ,mg    ,ms    ,scale2)*(mg2-ms2)
    +2.*real(SusyLoopIntegral::C0 (eps,eps,sHat(),ms,mg,ms))*
    (sqr(ms2-mg2)/sHat() + mg2);
  // QCD only piece for testing
  //output.finite =-8.+sqr(Constants::pi);
  output.finite *= loWeight();
  return output;
}

double MEPP2W2SleptonsPowheg::
qqbarME(vector<SpinorWaveFunction>    & sp ,
	vector<SpinorBarWaveFunction> & sbar ,
	ScalarWaveFunction & sca1,ScalarWaveFunction &sca2,
	bool first) const {
  // scale for the process
  const Energy2 q2(scale());
  // type of gauge boson
  int icharge = sca1.particle()->iCharge()+sca2.particle()->iCharge();
  tcPDPtr boson = icharge > 0 ? Wplus_ : Wminus_;
  // storage of the answer
  double me2(0.);
  // storage of the individual diagrams
  Complex diag(0.);
  ProductionMatrixElement pme(PDT::Spin1Half, PDT::Spin1Half, 
			      PDT::Spin0, PDT::Spin0);
  // loop over the helicities and calculate the matrix elements
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      VectorWaveFunction interV = FFWVertex_->evaluate(q2, 1, boson, sp[if1], 
						       sbar[if2]);
      diag = WSSVertex_->evaluate(q2, interV, sca2, sca1);
      // sum up the matrix elements
      me2 += norm(diag);
      pme(if1, if2, 0, 0) = diag;
    }
  }
  if(first) me_.reset(pme);
  return me2/12.;
}

double MEPP2W2SleptonsPowheg::loME(const cPDVector & particles,
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

double MEPP2W2SleptonsPowheg::
realME(const cPDVector & particles,
       const vector<Lorentz5Momentum> & momenta) const {
  vector<SpinorWaveFunction> sp(2);
  vector<SpinorBarWaveFunction> sbar(2);
  vector<VectorWaveFunction> gluon(2);
  // wavefunctions for the q qbar -> sf sf g process
  if(particles[0]->id()<=6&&particles[1]->id()<0) {
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
  // type of gauge boson
  int icharge = sca1.particle()->iCharge()+sca2.particle()->iCharge();
  tcPDPtr boson = icharge > 0 ? Wplus_ : Wminus_;
  // matrix element
  double output(0.);
  Complex diag[4]={0.,0.};
  Energy2 shat = scale();
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	// First Wplus diagram
 	SpinorWaveFunction inters = FFGVertex_->evaluate(shat,5,sp[ihel1].particle(),
							 sp[ihel1],gluon[ohel1]);
	VectorWaveFunction interV = FFWVertex_->evaluate(shat, 1, boson, inters, 
							 sbar[ihel2]);
	diag[0] = WSSVertex_->evaluate(shat, interV, sca2, sca1);
	// Second Wplus diagram
	SpinorBarWaveFunction interb = FFGVertex_->evaluate(shat,5,sbar[ihel1].particle(),
							    sbar[ihel2],gluon[ohel1]);
	interV = FFWVertex_->evaluate(shat, 1, boson, sp[ihel1], 
				      interb);
	diag[1] = WSSVertex_->evaluate(shat, interV, sca2, sca1);
	// add them up
	output += norm(diag[0]+diag[1]);
      }
    }
  }
  // colour and spin factors
  // q g or g qbar
  if(particles[0]->id()==ParticleID::g||
     particles[1]->id()==ParticleID::g) {
    output *= 1./24.;
  }
  // q qbar
  else  {
    output *= 1./9.;
  }
  // divided by 2 g_S^2
  return 0.5*output/norm(FFGVertex_->norm());
}

void MEPP2W2SleptonsPowheg::constructVertex(tSubProPtr sub) {
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
    ext[i]->spinInfo()->productionVertex(hv);  
}





