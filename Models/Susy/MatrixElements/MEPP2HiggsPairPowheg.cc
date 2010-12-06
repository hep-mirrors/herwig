// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2HiggsPairPowheg class.
//

#include "MEPP2HiggsPairPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG;

MEPP2HiggsPairPowheg::MEPP2HiggsPairPowheg() : maxFlavour_(5) {
  vector<unsigned int> mopt(2,1);
  massOption(mopt);
}

IBPtr MEPP2HiggsPairPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2HiggsPairPowheg::fullclone() const {
  return new_ptr(*this);
}

void MEPP2HiggsPairPowheg::persistentOutput(PersistentOStream & os) const {
  os << FFZVertex_ << FFPVertex_ << FFWVertex_ << ZSSVertex_ << PSSVertex_
     << WSSVertex_ << FFGVertex_ << Z0_ << gamma_ << Wplus_ << Wminus_
     << outgoing1_ << outgoing2_ << maxFlavour_;
}

void MEPP2HiggsPairPowheg::persistentInput(PersistentIStream & is, int) {
  is >> FFZVertex_ >> FFPVertex_ >> FFWVertex_ >> ZSSVertex_ >> PSSVertex_
     >> WSSVertex_ >> FFGVertex_ >> Z0_ >> gamma_ >> Wplus_ >> Wminus_
     >> outgoing1_ >> outgoing2_ >> maxFlavour_;
}

ClassDescription<MEPP2HiggsPairPowheg> MEPP2HiggsPairPowheg::initMEPP2HiggsPairPowheg;
// Definition of the static class description member.

void MEPP2HiggsPairPowheg::Init() {

  static ClassDocumentation<MEPP2HiggsPairPowheg> documentation
    ("The MEPP2HiggsPairPowheg class simulates the production of a pair of Higgs bosons,"
     "produced via a virtual s-channel vector boson in models with extended Higgs"
     " sectors.");

  static Reference<MEPP2HiggsPairPowheg,ParticleData> interfaceFirstOutgoing
    ("FirstOutgoing",
     "The first outgoing Higgs boson",
     &MEPP2HiggsPairPowheg::outgoing1_, false, false, true, false, false);

  static Reference<MEPP2HiggsPairPowheg,ParticleData> interfaceSecondOutgoing
    ("SecondOutgoing",
     "The second outgoing Higgs boson",
     &MEPP2HiggsPairPowheg::outgoing2_, false, false, true, false, false);

  static Parameter<MEPP2HiggsPairPowheg,int> interfaceMaxFlavour
    ("MaxFlavour",
     "The maximum flavour of the incoming quarks",
     &MEPP2HiggsPairPowheg::maxFlavour_, 5, 1, 5,
     false, false, Interface::limited);
}

unsigned int MEPP2HiggsPairPowheg::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2HiggsPairPowheg::orderInAlphaEW() const {
  return 2;
}

void MEPP2HiggsPairPowheg::doinit() {
  NLODrellYanBase::doinit();
  // get the ParticleData objects
  Wplus_  = getParticleData(ThePEG::ParticleID::Wplus);
  Wminus_ = getParticleData(ThePEG::ParticleID::Wminus);
  Z0_     = getParticleData(ThePEG::ParticleID::Z0);
  gamma_  = getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!hwsm)
    throw InitException() << "Must be the Herwig++ StandardModel class in "
 			  << "MEPP2HiggsPairPowheg::doinit" 
 			  << Exception::abortnow;
  FFZVertex_ = hwsm->vertexFFZ();
  FFPVertex_ = hwsm->vertexFFP();
  FFWVertex_ = hwsm->vertexFFW();
  FFGVertex_ = hwsm->vertexFFG();
  // work out which bosons we need
  int icharge = outgoing1_->iCharge()+outgoing2_->iCharge();
  if(icharge!=-3 && icharge!= +3 && icharge !=0)
    throw InitException() << "Can't produce selected particles via either "
			  << "photon/Z or W exchange in "
			  << "MEPP2HiggsPairPowheg::doinit" 
 			  << Exception::abortnow;
  // search the vertices and find the ones we want
  for(unsigned int ix = 0; ix < hwsm->numberOfVertices(); ++ix ) {
    VertexBasePtr vertex = hwsm->vertex(ix);
    vertex->init();
    // neutral gauge bosons
    if(icharge==0) {
      if(vertex->allowed(ParticleID::gamma,outgoing1_->id(),outgoing2_->id())) {
	if(PSSVertex_)
	  throw InitException() << "More than one possible vertex for "
				<< "photon exchange in "
				<< "MEPP2HiggsPairPowheg::doinit" 
				<< Exception::abortnow;
	PSSVertex_ = dynamic_ptr_cast<AbstractVSSVertexPtr>(vertex);
      }
      if(vertex->allowed(ParticleID::Z0,outgoing1_->id(),outgoing2_->id())) {
	if(ZSSVertex_)
	  throw InitException() << "More than one possible vertex for "
				<< "Z0 exchange in "
				<< "MEPP2HiggsPairPowheg::doinit" 
				<< Exception::abortnow;
	ZSSVertex_ = dynamic_ptr_cast<AbstractVSSVertexPtr>(vertex);
      }
    }
    // W-
    else if(icharge== -3) {
      if(vertex->allowed(ParticleID::Wplus,outgoing1_->id(),outgoing2_->id())) {
	if(WSSVertex_)
	  throw InitException() << "More than one possible vertex for "
				<< "W exchange in "
				<< "MEPP2HiggsPairPowheg::doinit" 
				<< Exception::abortnow;
	WSSVertex_ = dynamic_ptr_cast<AbstractVSSVertexPtr>(vertex);
      }
    }
    // W+
    else if(icharge== +3) {
      if(vertex->allowed(ParticleID::Wminus,outgoing1_->id(),outgoing2_->id())) {
	if(WSSVertex_)
	  throw InitException() << "More than one possible vertex for "
				<< "W exchange in "
				<< "MEPP2HiggsPairPowheg::doinit" 
				<< Exception::abortnow;
	WSSVertex_ = dynamic_ptr_cast<AbstractVSSVertexPtr>(vertex);
      }
    }
  }
  if( (    icharge  == 0 && (!PSSVertex_ && !ZSSVertex_) ) ||
      (abs(icharge) == 3 && !WSSVertex_) )
    throw InitException() << "Failed to find the vertices in "
			  << "MEPP2HiggsPairPowheg::doinit" 
			  << Exception::abortnow;
}

Selector<const ColourLines *> MEPP2HiggsPairPowheg::colourGeometries(tcDiagPtr) const {
  static const ColourLines c1("1 -2");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c1);
  return sel;
}

NLODrellYanBase::Singular MEPP2HiggsPairPowheg::virtualME() const {
  Singular output;
  output.eps2 = -2;
  output.eps1 = -3;
  output.finite =-8.+sqr(Constants::pi);
  output.finite *= loWeight();
  return output;
}

void MEPP2HiggsPairPowheg::getDiagrams() const {
  int icharge = outgoing1_->iCharge()+outgoing2_->iCharge();
  assert( icharge==0 || abs(icharge) == 3 );
  // neutral final-state
  if(icharge==0) {
    bool charged = outgoing1_->charged();
    // loop over the processes we need
    for(int i = 1; i <= maxFlavour_; ++i) {
      tcPDPtr q  = getParticleData(i);
      tcPDPtr qb = q->CC();
      // always Z 
      add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Z0_   , 3,
		   outgoing1_, 3, outgoing2_, -1)));
      if(charged)
	add(new_ptr((Tree2toNDiagram(2), q, qb, 1, gamma_, 3,
		     outgoing1_, 3, outgoing2_, -2)));
    }
  }
  // charged final-state
  else if(icharge==3) {
    // charge conjugates
    tcPDPtr outgoing1c = outgoing1_->CC() ?
      tcPDPtr(outgoing1_->CC()) : tcPDPtr(outgoing1_);
    tcPDPtr outgoing2c = outgoing2_->CC() ?
      tcPDPtr(outgoing2_->CC()) : tcPDPtr(outgoing2_);
    // loop over the Wplus processes we need
    for(int i = 2; i <= maxFlavour_; i+=2 ) {
      tcPDPtr qi  = getParticleData(i);
      tcPDPtr qib = qi->CC();
      for(int j = 1; j <= maxFlavour_; j+=2 ) {
	tcPDPtr qj  = getParticleData(j);
	tcPDPtr qjb = qj->CC();
	if(icharge==3) {
	  add(new_ptr((Tree2toNDiagram(2), qi, qjb, 1, Wplus_,
		       3, outgoing1_, 3, outgoing2_, -1)));
	  add(new_ptr((Tree2toNDiagram(2), qj, qib, 1, Wminus_,
		       3, outgoing1c, 3, outgoing2c, -1)));
	}
	else {
	  add(new_ptr((Tree2toNDiagram(2), qi, qjb, 1, Wplus_,
		       3, outgoing1c, 3, outgoing2c, -1)));
	  add(new_ptr((Tree2toNDiagram(2), qj, qib, 1, Wminus_,
		       3, outgoing1_, 3, outgoing2_, -1)));
	}
      }
    } 
  }
}

Selector<MEBase::DiagramIndex> 
MEPP2HiggsPairPowheg::diagrams(const DiagramVector & diags ) const {
  int icharge = outgoing1_->iCharge()+outgoing2_->iCharge();
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if(icharge==0) {
      if ( diags[i]->id() == -1 ) sel.insert(meInfo()[0], i);
      else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
    }
    else {
      sel.insert(1., i);
    }
  }
  return sel;
}

double MEPP2HiggsPairPowheg::loME(const cPDVector & particles,
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

double MEPP2HiggsPairPowheg::
qqbarME(vector<SpinorWaveFunction>    & sp ,
	vector<SpinorBarWaveFunction> & sbar ,
	ScalarWaveFunction & sca1,ScalarWaveFunction &sca2,
	bool first) const {
  // scale for the process
  const Energy2 q2(scale());
  // type of gauge boson
  int icharge = sca1.particle()->iCharge()+sca2.particle()->iCharge();
  double me2(0.);
  ProductionMatrixElement pme(PDT::Spin1Half, PDT::Spin1Half, 
 			      PDT::Spin0, PDT::Spin0);
  // neutral gauge boson
  if(icharge==0) {
    // storage of the matrix elements for specific diagrams
    vector<double> me(2, 0.);
    // storage of the individual diagrams
    vector<Complex> diag(2, Complex(0.));
    // loop over the helicities and calculate the matrix elements
    for(unsigned int if1 = 0; if1 < 2; ++if1) {
      for(unsigned int if2 = 0; if2 < 2; ++if2) {
	// Z vertex
	if(ZSSVertex_) {
	  VectorWaveFunction interV = FFZVertex_->
	    evaluate(q2, 1, Z0_, sp[if1], sbar[if2]);
	  diag[0] = ZSSVertex_->evaluate(q2, interV, sca2, sca1);
	}
	// photon vertex
	if(PSSVertex_) {
	  VectorWaveFunction interV = FFPVertex_->
	    evaluate(q2, 1, gamma_, sp[if1], sbar[if2]);
	  diag[1] = PSSVertex_->evaluate(q2, interV, sca2, sca1);
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
    }
  }
  // charged gauge boson
  else if(abs(icharge)==3) {
    // s-channel boson
    tcPDPtr boson = icharge > 0 ? Wplus_ : Wminus_;
    // storage of the individual diagrams
    Complex diag(0.);
    // loop over the helicities and calculate the matrix elements
    for(unsigned int if1 = 0; if1 < 2; ++if1) {
      for(unsigned int if2 = 0; if2 < 2; ++if2) {
	VectorWaveFunction interV = FFWVertex_->
	  evaluate(q2, 1, boson, sp[if1], sbar[if2]);
	diag = WSSVertex_->evaluate(q2, interV, sca2, sca1);
	// sum up the matrix elements
	me2 += norm(diag);
	pme(if1, if2, 0, 0) = diag;
      }
    }
  }
  else
    assert(false);
  // return the answer
  if(first) me_.reset(pme);
  return me2/12.;
}

double MEPP2HiggsPairPowheg::realME(const cPDVector & particles,
				    const vector<Lorentz5Momentum> & momenta) const {
  vector<SpinorWaveFunction> sp(2);
  vector<SpinorBarWaveFunction> sbar(2);
  vector<VectorWaveFunction> gluon(2);
  // wavefunctions for the q qbar -> H H g process
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
  // matrix element
  double output(0.);
  Complex diag[4]={0.,0.,0.,0.};
  Energy2 shat = scale();
  // neutral gauge bosons
  if(icharge==0) {
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	  SpinorWaveFunction inters = FFGVertex_->evaluate(shat,5,sp[ihel1].particle(),
							   sp[ihel1],gluon[ohel1]);
	  SpinorBarWaveFunction interb = FFGVertex_->evaluate(shat,5,sbar[ihel1].particle(),
							      sbar[ihel2],gluon[ohel1]);
	  if(ZSSVertex_) {
	    // first Z diagram
	    VectorWaveFunction interV = FFZVertex_->evaluate(shat, 1, Z0_, inters, 
							     sbar[ihel2]);
	    diag[0] = ZSSVertex_->evaluate(shat, interV, sca2, sca1);
	    // second Z diagram
	    interV = FFZVertex_->evaluate(shat, 1, Z0_, sp[ihel1], 
					  interb);
	    diag[1] = ZSSVertex_->evaluate(shat, interV, sca2, sca1);
	  }
	  if(PSSVertex_) {
	    // first photon diagram
	    VectorWaveFunction interV = FFPVertex_->evaluate(shat, 1, gamma_, inters, 
							     sbar[ihel2]);
	    diag[2] = PSSVertex_->evaluate(shat, interV, sca2, sca1);
	    // second photon diagram
	    interV = FFPVertex_->evaluate(shat, 1, gamma_, sp[ihel1], 
					  interb);
	    diag[3] = PSSVertex_->evaluate(shat, interV, sca2, sca1);
	  }
	  // add them up
	  output += norm(diag[0]+diag[1]+diag[2]+diag[3]);
	}
      }
    }
  }
  // charged gauge bosons
  else if(abs(icharge)==3) {
    tcPDPtr boson = icharge > 0 ? Wplus_ : Wminus_;
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
  }
  else
    assert(false);
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

void MEPP2HiggsPairPowheg::constructVertex(tSubProPtr sub) {
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

Energy2 MEPP2HiggsPairPowheg::scale() const {
//   return sHat();
  return sqr(0.5*(mePartonData()[2]->mass()+mePartonData()[3]->mass()));
}
