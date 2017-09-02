// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVSDecayer class.
//

#include "VVSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr VVSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VVSDecayer::fullclone() const {
  return new_ptr(*this);
}
void VVSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      VertexBasePtr vertex,
			      map<ShowerInteraction,VertexBasePtr> &,
			      const vector<map<ShowerInteraction,VertexBasePtr> > &,
			      map<ShowerInteraction,VertexBasePtr>) {
  decayInfo(incoming,outgoing);
  vertex_             = dynamic_ptr_cast<AbstractVVSVertexPtr>(vertex);
  perturbativeVertex_ = dynamic_ptr_cast<VVSVertexPtr>        (vertex);
}

void VVSDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << perturbativeVertex_;
}

void VVSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> perturbativeVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VVSDecayer,GeneralTwoBodyDecayer>
describeHerwigVVSDecayer("Herwig::VVSDecayer", "Herwig.so");

void VVSDecayer::Init() {

  static ClassDocumentation<VVSDecayer> documentation
    ("The VVSDecayer class implements the decay of a vector"
     " to a vector and a scalar");

}

double VVSDecayer::me2(const int , const Particle & inpart,
 		       const ParticleVector & decay,
		       MEOption meopt) const {
  bool massless = ( decay[0]->id()==ParticleID::gamma || 
		    decay[0]->id()==ParticleID::g );
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    VectorWaveFunction::
      constructSpinInfo(vectors_[1],decay[0],outgoing,true,massless);
    ScalarWaveFunction::
      constructSpinInfo(decay[1],outgoing,true);
    return 0.;
  }
  VectorWaveFunction::
    calculateWaveFunctions(vectors_[1],decay[0],outgoing,massless);
  ScalarWaveFunction sca(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int in=0;in<3;++in) {
    for(unsigned int out=0;out<3;++out) {
      if(massless&&out==1) ++out;
      (*ME())(in,out,0) = 
	vertex_->evaluate(scale,vectors_[0][in],vectors_[1][out],sca);
    }
  }
  double output=(ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy VVSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_) {
    Energy2 scale(sqr(inpart.second));
    double mu1sq = sqr(outa.second/inpart.second);
    double mu2sq = sqr(outb.second/inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if( outb.first->iSpin() == PDT::Spin0 )
      perturbativeVertex_->setCoupling(sqr(inpart.second), in, 
				       outa.first, outb.first);
    else {
      perturbativeVertex_->setCoupling(sqr(inpart.second), in, 
				       outb.first, outa.first);
      swap(mu1sq, mu2sq);
    }
    double vn = norm(perturbativeVertex_->norm());
    if(vn == ZERO || mu1sq == ZERO) return ZERO;
    double me2 = 2. + 0.25*sqr(1. + mu1sq - mu2sq)/mu1sq;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = vn*me2*pcm/(24.*Constants::pi)/scale*UnitRemoval::E2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

