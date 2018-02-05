// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericHGGVertex class.
//

#include "GenericHGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;

GenericHGGVertex::GenericHGGVertex() : setup_(false), q2Last_(ZERO), coupLast_(0.), idLast_(0) {
  orderInGs(2);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

IBPtr GenericHGGVertex::clone() const {
  return new_ptr(*this);
}

IBPtr GenericHGGVertex::fullclone() const {
  return new_ptr(*this);
}

void GenericHGGVertex::persistentOutput(PersistentOStream & os) const {
  os << bosons_ << setup_ << vertices_ << model_;
}

void GenericHGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> bosons_ >> setup_ >> vertices_ >> model_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GenericHGGVertex,VVSLoopVertex>
describeHerwigGenericHGGVertex("Herwig::GenericHGGVertex", "Herwig.so");

void GenericHGGVertex::Init() {

  static ClassDocumentation<GenericHGGVertex> documentation
    ("The GenericHGGVertex class implements the coupling of"
     " the Higgs bosons to gluons in a generic model.");

  static RefVector<GenericHGGVertex,ParticleData> interfaceBosons
    ("Bosons",
     "The Higgs bosons in the model.",
     &GenericHGGVertex::bosons_, -1, false, false, true, false, false);

}

void GenericHGGVertex::doinit() {
  for(unsigned int ix=0;ix<bosons_.size();++ix) {
    addToList(21,21,bosons_[ix]->id());
  }
  VVSLoopVertex::doinit();
  if(loopToolsInitialized()) Looptools::ltexi();
}

void GenericHGGVertex::setCoupling (Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				    tcPDPtr part3) {
  if(!setup_) initializeVertex();
  assert(part1->id()==ParticleID::g && part2->id()==ParticleID::g);
  // find the particles in the loop
  map<cPDPtr,vector<Interaction> >::iterator it = vertices_.find(part3);
  // check there are some
  if(it==vertices_.end()) {
    norm(0.);
    return;
  }
  Looptools::clearcache();
  // overall factor
  if( q2 != q2Last_ || coupLast_ == 0. || part3->id() != idLast_ ) {
    q2Last_ = q2;
    idLast_ = part3->id();
    coupLast_ = sqr(strongCoupling(q2));
    // loop over the loop particles
    masses.clear();
    type.clear();
    couplings.clear();
    setNParticles(it->second.size());
    for(unsigned int ix=0;ix<it->second.size();++ix) {
      masses.push_back(model_->mass(q2,it->second[ix].particle));
      if(it->second[ix].particle->iSpin()==PDT::Spin0) {
	type.push_back(PDT::Spin0);
	it->second[ix].scalar->setCoupling(q2,part3,it->second[ix].particle,it->second[ix].particle->CC());
	couplings.push_back(make_pair(0.5*it->second[ix].scalar->norm(),0.5*it->second[ix].scalar->norm()));
      }
      else if(it->second[ix].particle->iSpin()==PDT::Spin1Half) {
	type.push_back(PDT::Spin1Half);
	assert(it->second[ix].fermion);
	it->second[ix].fermion->setCoupling(q2,it->second[ix].particle,it->second[ix].particle->CC(),part3);
	Complex coupling = it->second[ix].fermion->norm();
	Complex lc = it->second[ix].fermion->left ();
	Complex rc = it->second[ix].fermion->right();
	couplings.push_back(make_pair(0.5*coupling*lc,0.5*coupling*rc));
      }
      else
	assert(false);
    }
    VVSLoopVertex::setCoupling(q2, part1, part2, part3);
  }
  norm(coupLast_);
}

void GenericHGGVertex::initializeVertex() {
  // get the model
  model_ = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if(!model_) throw InitException();
  // loop over all the vertices
  unsigned int nv(model_->numberOfVertices());
  for(unsigned int ib=0;ib<bosons_.size();++ib) {
    for(unsigned int iv = 0; iv < nv; ++iv) {
      // 3-point vertex with boson as incoming
      if( model_->vertex(iv)->getNpoint()>3) continue;
      if( !model_->vertex(iv)->isIncoming(bosons_[ib])) continue;
      for(unsigned int il = 0; il < 3; ++il) { 
	tPDVector decaylist = model_->vertex(iv)->search(il, bosons_[ib]);
	tPDVector::size_type nd = decaylist.size();
	for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
	  tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
	  if( pb->id() == bosons_[ib]->id() ) swap(pa, pb);
	  if( pc->id() == bosons_[ib]->id() ) swap(pa, pc);
	  // check coloured and particle antiparticle
	  if( pb->CC() != pc || !pb->coloured() || !pc->coloured()) 
	    continue;
	  // only know how to deal with triplets
	  assert(pb->iColour()==PDT::Colour3 || pb->iColour()==PDT::Colour3bar);
	  // scalar loop
	  if(pb->iSpin()==PDT::Spin0) {
	    SSSVertexPtr vertex = dynamic_ptr_cast<SSSVertexPtr>(model_->vertex(iv));
	    if(!vertex) continue;
	    map<cPDPtr,vector<Interaction> >::iterator it = vertices_.find(bosons_[ib]);
	    if(it!=vertices_.end()) {
	      it->second.push_back(Interaction(pb->id()>0?pb:pc,vertex,FFSVertexPtr()));
	    }
	    else {
	      vertices_.insert(make_pair(bosons_[ib],vector<Interaction>(1,Interaction(pb->id()>0?pb:pc,vertex,FFSVertexPtr()))));
	    }
	  }
	  // fermion loop
	  else if(pb->iSpin()==PDT::Spin1Half) {
	    FFSVertexPtr vertex = dynamic_ptr_cast<FFSVertexPtr>(model_->vertex(iv));
	    if(!vertex) continue;
	    map<cPDPtr,vector<Interaction> >::iterator it = vertices_.find(bosons_[ib]);
	    if(it!=vertices_.end()) {
	      it->second.push_back(Interaction(pb->id()>0?pb:pc,SSSVertexPtr(),vertex));
	    }
	    else {
	     vertices_.insert(make_pair(bosons_[ib],vector<Interaction>(1,Interaction(pb->id()>0?pb:pc,SSSVertexPtr(),vertex))));
	    }
	  }
	  else
	    assert(false);
	}
      }
    }
  }
  // set up now
  setup_ = true;
}
