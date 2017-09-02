// -*- C++ -*-
//
// SSVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSVDecayer class.
//

#include "SSVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr SSVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SSVDecayer::fullclone() const {
  return new_ptr(*this);
}

void SSVDecayer::doinit() {
  perturbativeVertex_      = dynamic_ptr_cast<VSSVertexPtr>         (vertex());
  abstractVertex_          = dynamic_ptr_cast<AbstractVSSVertexPtr> (vertex());
  abstractIncomingVertex_  = dynamic_ptr_cast<AbstractVSSVertexPtr> (incomingVertex());
  abstractFourPointVertex_ = dynamic_ptr_cast<AbstractVVSSVertexPtr>(getFourPointVertex());

  if (outgoingVertices()[0]){
    if (outgoingVertices()[0]->getName()==VertexType::VSS){
      abstractOutgoingVertexS_   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outgoingVertices()[0]);
      abstractOutgoingVertexV_   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outgoingVertices()[1]);
    }
    else {
      abstractOutgoingVertexS_   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outgoingVertices()[1]);
      abstractOutgoingVertexV_   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outgoingVertices()[0]);
    }
  }
  else if (outgoingVertices()[1]){
    if (outgoingVertices()[1]->getName()==VertexType::VSS){
      abstractOutgoingVertexS_   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outgoingVertices()[1]);
      abstractOutgoingVertexV_   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outgoingVertices()[0]);
    }
    else {
      abstractOutgoingVertexS_   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outgoingVertices()[0]);
      abstractOutgoingVertexV_   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outgoingVertices()[1]);
    }
  }
  GeneralTwoBodyDecayer::doinit();
}

void SSVDecayer::persistentOutput(PersistentOStream & os) const {
  os << abstractVertex_           << perturbativeVertex_
     << abstractIncomingVertex_   << abstractOutgoingVertexS_
     << abstractOutgoingVertexV_  << abstractFourPointVertex_;
}

void SSVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> abstractVertex_           >> perturbativeVertex_
     >> abstractIncomingVertex_   >> abstractOutgoingVertexS_
     >> abstractOutgoingVertexV_  >> abstractFourPointVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSVDecayer,GeneralTwoBodyDecayer>
describeHerwigSSVDecayer("Herwig::SSVDecayer", "Herwig.so");

void SSVDecayer::Init() {

  static ClassDocumentation<SSVDecayer> documentation
    ("This implements the decay of a scalar to a vector and a scalar");

}

double SSVDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  unsigned int isc(0),ivec(1);
  if(decay[0]->dataPtr()->iSpin() != PDT::Spin0) swap(isc,ivec);
  if(!ME()) {
    if(ivec==1)
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1)));
    else
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin0)));
  }
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&inpart),incoming);
    swave_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    ScalarWaveFunction::
      constructSpinInfo(decay[isc],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(vector_,decay[ivec],outgoing,true,false);
  }
  VectorWaveFunction::
    calculateWaveFunctions(vector_,decay[ivec],outgoing,false);
  ScalarWaveFunction sca(decay[isc]->momentum(),decay[isc]->dataPtr(),outgoing);
  Energy2 scale(sqr(inpart.mass()));
  //make sure decay matrix element is in the correct order
  double output(0.);
  if(ivec == 0) {
    for(unsigned int ix = 0; ix < 3; ++ix)
      (*ME())(0, ix, 0) = abstractVertex_->evaluate(scale,vector_[ix],sca, swave_);
  }
  else {
    for(unsigned int ix = 0; ix < 3; ++ix)
      (*ME())(0, 0, ix) = abstractVertex_->evaluate(scale,vector_[ix],sca,swave_);
  }
  output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy SSVDecayer:: partialWidth(PMPair inpart, PMPair outa, 
				 PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_) {
    double mu1sq(sqr(outa.second/inpart.second)),
      mu2sq(sqr(outb.second/inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if(outa.first->iSpin() == PDT::Spin0) {
      perturbativeVertex_->setCoupling(sqr(inpart.second), outb.first, outa.first,in);
    }
    else {
      swap(mu1sq,mu2sq);
      perturbativeVertex_->setCoupling(sqr(inpart.second), outa.first, outb.first,in);
    }
    double me2(0.);
    if(mu2sq == 0.) 
      me2 = -2.*mu1sq - 2.;
    else
      me2 = ( sqr(mu2sq - mu1sq) - 2.*(mu2sq + mu1sq) + 1. )/mu2sq;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					       outb.second);
    Energy output = pcm*me2*norm(perturbativeVertex_->norm())/8./Constants::pi;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}


double  SSVDecayer::threeBodyME(const int , const Particle & inpart,
				const ParticleVector & decay, MEOption meopt) {

  int iscal (0), ivect (1), iglu (2);
  // get location of outgoing scalar/vector
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin0) swap(iscal,ivect);

  // no emissions from massive vectors
  if (abstractOutgoingVertexV_ && decay[ivect]->dataPtr()->mass()!=ZERO)
    throw Exception()
      << "No dipoles available for massive vectors in SSVDecayer::threeBodyME"
      << Exception::runerror;

  if(meopt==Initialize) {
    // create scalar wavefunction for decaying particle
    ScalarWaveFunction::calculateWaveFunctions(rho3_,const_ptr_cast<tPPtr>(&inpart),incoming);
    swave3_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    ScalarWaveFunction::
      constructSpinInfo(decay[iscal],outgoing,true);
     VectorWaveFunction::
      constructSpinInfo(vector3_,decay[ivect],outgoing,true,false);
    VectorWaveFunction::
      constructSpinInfo(gluon_,  decay[iglu ],outgoing,true,false);
    return 0.;
  }
  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);  
  if(nflow==2) cfactors[0][1]=cfactors[1][0];
  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin0, PDT::Spin0,
								       PDT::Spin1, PDT::Spin1)));

  // create wavefunctions
  ScalarWaveFunction scal_(decay[iscal]->momentum(),  decay[iscal]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(vector3_,decay[ivect],outgoing,false);
  VectorWaveFunction::calculateWaveFunctions(gluon_,  decay[iglu ],outgoing,true );

  // // gauge invariance test
  // gluon_.clear();
  // for(unsigned int ix=0;ix<3;++ix) {
  //   if(ix==1) gluon_.push_back(VectorWaveFunction());
  //   else {
  //     gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  // 					  decay[iglu ]->dataPtr(),10,
  // 					  outgoing));
  //   }
  // }


  if (! ((abstractIncomingVertex_  && (abstractOutgoingVertexS_ || abstractOutgoingVertexV_)) ||
	 (abstractOutgoingVertexS_ &&  abstractOutgoingVertexV_)))
    throw Exception()
      << "Invalid vertices for QCD radiation in SSV decay in SSVDecayer::threeBodyME"
      << Exception::runerror;


  // sort out colour flows
  int S(1), V(2);
  if (decay[iscal]->dataPtr()->iColour()==PDT::Colour3bar && 
      decay[ivect]->dataPtr()->iColour()==PDT::Colour8)
    swap(S,V);
  else if (decay[ivect]->dataPtr()->iColour()==PDT::Colour3 && 
	   decay[iscal]->dataPtr()->iColour()==PDT::Colour8)
    swap(S,V);

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);

  for(unsigned int iv = 0; iv < 3; ++iv) {
    for(unsigned int ig = 0; ig < 2; ++ig) {
      // radiation from the incoming scalar
      if(inpart.dataPtr()->coloured()) {
	assert(abstractIncomingVertex_);
	ScalarWaveFunction scalarInter = 
	  abstractIncomingVertex_->evaluate(scale,3,inpart.dataPtr(),
					    gluon_[2*ig],swave3_,inpart.mass());
	
	if (swave3_.particle()->PDGName()!=scalarInter.particle()->PDGName())
	  throw Exception()
	    << swave3_    .particle()->PDGName() << " was changed to " 
	    << scalarInter.particle()->PDGName() << " in SSVDecayer::threeBodyME"
	    << Exception::runerror;
	
	double gs    = abstractIncomingVertex_->strongCoupling(scale);
	double sign  = 1.;//inpart.dataPtr()->id()>0 ? 1:-1;	
	Complex diag = sign * abstractVertex_->evaluate(scale,vector3_[iv],scal_,scalarInter)/gs;
	for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	  (*ME[colourFlow[0][ix].first])(0, 0, iv, ig) += 
	    colourFlow[0][ix].second*diag; 
	}
      }
      // radiation from the outgoing scalar
      if(decay[iscal]->dataPtr()->coloured()) {
	assert(abstractOutgoingVertexS_);
	// ensure you get correct outgoing particle from first vertex
	tcPDPtr off = decay[iscal]->dataPtr();
	if(off->CC()) off = off->CC();
	ScalarWaveFunction scalarInter = 
	  abstractOutgoingVertexS_->evaluate(scale,3,off,gluon_[2*ig],scal_,decay[iscal]->mass());
	
	if (scal_.particle()->PDGName()!=scalarInter.particle()->PDGName())
	  throw Exception()
	    << scal_      .particle()->PDGName() << " was changed to " 
	    << scalarInter.particle()->PDGName() << " in SSVDecayer::threeBodyME"
	    << Exception::runerror;

	double gs    = abstractOutgoingVertexS_->strongCoupling(scale);
	double sign  = 1.;//decay[iscal]->dataPtr()->id()>0 ? -1:1;
	Complex diag = sign*abstractVertex_->evaluate(scale,vector3_[iv],scalarInter,swave3_)/gs;
	for(unsigned int ix=0;ix<colourFlow[S].size();++ix) {
	  (*ME[colourFlow[S][ix].first])(0, 0, iv, ig) += 
	    colourFlow[S][ix].second*diag;
	}
      }

      // radiation from outgoing vector
      if(decay[ivect]->dataPtr()->coloured()) {
	assert(abstractOutgoingVertexV_);
	// ensure you get correct outgoing particle from first vertex
	tcPDPtr off = decay[ivect]->dataPtr();
	if(off->CC()) off = off->CC();
	VectorWaveFunction  vectorInter = 
	  abstractOutgoingVertexV_->evaluate(scale,3,off,gluon_[2*ig],
					     vector3_[iv],decay[ivect]->mass());
	    
	if(vector3_[iv].particle()->PDGName()!=vectorInter.particle()->PDGName())
	  throw Exception()
	    << vector3_[iv].particle()->PDGName() << " was changed to " 
	    << vectorInter. particle()->PDGName() << " in SSVDecayer::threeBodyME"
	    << Exception::runerror; 

	double sign  =  1.;//decay[iscal]->id()>0 ? -1:1;
	double gs    = abstractOutgoingVertexV_->strongCoupling(scale);	
	Complex diag =  sign*abstractVertex_->evaluate(scale,vectorInter,scal_,swave3_)/gs;
	for(unsigned int ix=0;ix<colourFlow[V].size();++ix) {
	  (*ME[colourFlow[V][ix].first])(0, 0, iv, ig) += 
	    colourFlow[V][ix].second*diag;
	}
      }
      // radiation from 4 point vertex
      if (abstractFourPointVertex_){
	double gs    = abstractFourPointVertex_->strongCoupling(scale);
	double sign  =  decay[iscal]->id()>0 ? -1:-1;
	Complex diag =  sign*abstractFourPointVertex_->evaluate(scale, gluon_[2*ig], vector3_[iv],
							       scal_, swave3_)/gs;
	for(unsigned int ix=0;ix<colourFlow[3].size();++ix) {
	  (*ME[colourFlow[3][ix].first])(0, 0, iv, ig) += 
	     colourFlow[3][ix].second*diag;
	}
      }
    }
  }
  
  // contract matrices
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
    }
  }
  output*=(4.*Constants::pi);
  // return the answer
  return output;
}
