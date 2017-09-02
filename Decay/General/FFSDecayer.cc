// -*- C++ -*-
//
// FFSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFSDecayer class.
//

#include "FFSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr FFSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FFSDecayer::fullclone() const {
  return new_ptr(*this);
}

void FFSDecayer::doinit() {
  perturbativeVertex_        = dynamic_ptr_cast<FFSVertexPtr>        (vertex());
  abstractVertex_            = dynamic_ptr_cast<AbstractFFSVertexPtr>(vertex());
  abstractIncomingVertex_    = dynamic_ptr_cast<AbstractFFVVertexPtr>(incomingVertex());

  if (outgoingVertices()[0]){
    if (outgoingVertices()[0]->getName()==VertexType::FFV){
      abstractOutgoingVertexF_   = dynamic_ptr_cast<AbstractFFVVertexPtr>(outgoingVertices()[0]);
      abstractOutgoingVertexS_   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outgoingVertices()[1]);
    }
    else {
      abstractOutgoingVertexF_   = dynamic_ptr_cast<AbstractFFVVertexPtr>(outgoingVertices()[1]);
      abstractOutgoingVertexS_   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outgoingVertices()[0]);
    }
  }
  else if (outgoingVertices()[1]){
    if (outgoingVertices()[1]->getName()==VertexType::FFV){
      abstractOutgoingVertexF_   = dynamic_ptr_cast<AbstractFFVVertexPtr>(outgoingVertices()[1]);
      abstractOutgoingVertexS_   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outgoingVertices()[0]);
    }
    else {
      abstractOutgoingVertexF_   = dynamic_ptr_cast<AbstractFFVVertexPtr>(outgoingVertices()[0]);
      abstractOutgoingVertexS_   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outgoingVertices()[1]);
    }
  }
  GeneralTwoBodyDecayer::doinit();
}

void FFSDecayer::persistentOutput(PersistentOStream & os) const {
  os << perturbativeVertex_       << abstractVertex_
     << abstractIncomingVertex_   << abstractOutgoingVertexF_
     << abstractOutgoingVertexS_;
}

void FFSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> perturbativeVertex_       >> abstractVertex_
     >> abstractIncomingVertex_   >> abstractOutgoingVertexF_
     >> abstractOutgoingVertexS_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FFSDecayer,GeneralTwoBodyDecayer>
describeHerwigFFSDecayer("Herwig::FFSDecayer", "Herwig.so");

void FFSDecayer::Init() {

  static ClassDocumentation<FFSDecayer> documentation
    ("The FFSDecayer class implements the decay of a fermion to "
     "a fermion and a scalar.");

}

double FFSDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0)));
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id()    > 0? 0:1;
  else                              itype[0] = 2;
  if(decay[0]->dataPtr()->CC())     itype[1] = decay[0]->id() > 0? 0:1;
  else                              itype[1] = 2;
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));

  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(wave_,rho_,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(wave_[0].wave().Type() != SpinorType::u)
	for(unsigned int ix = 0; ix < 2; ++ix) wave_   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar_,rho_,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(wavebar_[0].wave().Type() != SpinorType::v)
	for(unsigned int ix = 0; ix < 2; ++ix) wavebar_[ix].conjugate();
    }
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(ferm) {
      SpinorWaveFunction::
	constructSpinInfo(wave_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(wavebar_,decay[0],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(wavebar_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(wave_,decay[0],outgoing,true);
    }
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
  }
  if(ferm)
    SpinorBarWaveFunction::
      calculateWaveFunctions(wavebar_,decay[0],outgoing);
  else
    SpinorWaveFunction::
      calculateWaveFunctions(wave_   ,decay[0],outgoing);
  ScalarWaveFunction scal(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      if(ferm) (*ME())(if1, if2, 0) = 
	abstractVertex_->evaluate(scale,wave_[if1],wavebar_[if2],scal);
      else     (*ME())(if2, if1, 0) = 
	abstractVertex_->evaluate(scale,wave_[if1],wavebar_[if2],scal);
    }
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy FFSDecayer::partialWidth(PMPair inpart, PMPair outa,
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_) {
    double mu1(0.),mu2(0.);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if(outa.first->iSpin() == PDT::Spin1Half) {
      mu1 = outa.second/inpart.second;
      mu2 = outb.second/inpart.second;
      perturbativeVertex_->setCoupling(sqr(inpart.second), in, outa.first, outb.first);
    }
    else {
      mu1 = outb.second/inpart.second;
      mu2 = outa.second/inpart.second;
      perturbativeVertex_->setCoupling(sqr(inpart.second), in, outb.first, outa.first);
      
    }
    double c2 = norm(perturbativeVertex_->norm());
    Complex cl = perturbativeVertex_->left();
    Complex cr = perturbativeVertex_->right();
    double me2 = c2*( (norm(cl) + norm(cr))*(1. + sqr(mu1) - sqr(mu2))
		      + 2.*mu1*(conj(cl)*cr + conj(cr)*cl).real() );
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					outb.second);
    Energy output = me2*pcm/16./Constants::pi;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}


double FFSDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay, MEOption meopt) {

  int iscal (0), iferm (1), iglu (2);
  // get location of outgoing fermion/scalar
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin0) swap(iscal,iferm);
  // work out whether inpart is a fermion or antifermion
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[iferm]->dataPtr()->CC()) itype[1] = decay[iferm]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;

  bool ferm(itype[0] == 0 || itype[1] == 0 || 
	   (itype[0] == 2 && itype[1] == 2 && decay[iscal]->id() < 0));  

  if(meopt==Initialize) {
    // create spinor (bar) for decaying particle
    if(ferm) {
      SpinorWaveFunction::calculateWaveFunctions(wave3_, rho3_, const_ptr_cast<tPPtr>(&inpart), 
						 incoming);
      if(wave3_[0].wave().Type() != SpinorType::u)
   	for(unsigned int ix = 0; ix < 2; ++ix) wave3_[ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar3_,rho3_, const_ptr_cast<tPPtr>(&inpart), 
						    incoming);
      if(wavebar3_[0].wave().Type() != SpinorType::v)
   	for(unsigned int ix = 0; ix < 2; ++ix) wavebar3_[ix].conjugate();
    }
  }
  // setup spin information when needed 
  if(meopt==Terminate) {
    if(ferm) {
      SpinorWaveFunction::
	constructSpinInfo(wave3_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(wavebar3_,decay[iferm],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(wavebar3_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(wave3_,decay[iferm],outgoing,true);
    }
    ScalarWaveFunction::constructSpinInfo(        decay[iscal],outgoing,true);
    VectorWaveFunction::constructSpinInfo(gluon_, decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calulate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  if(nflow==2) cfactors[0][1] = cfactors[1][0];

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half, PDT::Spin0,
								       PDT::Spin1Half, PDT::Spin1)));
  // create wavefunctions
  if (ferm)  SpinorBarWaveFunction::calculateWaveFunctions(wavebar3_, decay[iferm],outgoing);
  else       SpinorWaveFunction::   calculateWaveFunctions(wave3_   , decay[iferm],outgoing);
  
  ScalarWaveFunction swave3_(decay[iscal]->momentum(), decay[iscal]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(gluon_,   decay[iglu ],outgoing,true);

  // // gauge invariance test
  //   gluon_.clear();
  // for(unsigned int ix=0;ix<3;++ix) {
  //   if(ix==1) gluon_.push_back(VectorWaveFunction());
  //   else {
  //     gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  // 					  decay[iglu ]->dataPtr(),10,
  // 					  outgoing));
  //   }
  // }

  if (! ((abstractIncomingVertex_  && (abstractOutgoingVertexF_ || abstractOutgoingVertexS_)) ||
	 (abstractOutgoingVertexF_ &&  abstractOutgoingVertexS_)))
    throw Exception()
      << "Invalid vertices for QCD radiation in FFS decay in FFSDecayer::threeBodyME"
      << Exception::runerror;


  // sort out colour flows
  int F(1), S(2);
  if (decay[iscal]->dataPtr()->iColour()==PDT::Colour3 && 
      decay[iferm]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,S);
  else if (decay[iferm]->dataPtr()->iColour()==PDT::Colour3bar && 
	   decay[iscal]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,S);


  Complex diag;
  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);

  for(unsigned int ifi = 0; ifi < 2; ++ifi) {
    for(unsigned int ifo = 0; ifo < 2; ++ifo) {
      for(unsigned int ig = 0; ig < 2; ++ig) {
   	// radiation from the incoming fermion
   	if(inpart.dataPtr()->coloured()) {
   	  assert(abstractIncomingVertex_);
	  double gs = abstractIncomingVertex_->strongCoupling(scale);	  
	  if (ferm){
	    SpinorWaveFunction spinorInter =
	      abstractIncomingVertex_->evaluate(scale,3,inpart.dataPtr(),wave3_[ifi],
					       gluon_[2*ig],inpart.mass());

	    if (wave3_[ifi].particle()->PDGName()!=spinorInter.particle()->PDGName())
	      throw Exception()
		<< wave3_[ifi].particle()->PDGName()  << " was changed to " 
		<< spinorInter.particle()->PDGName()  << " in FFSDecayer::threeBodyME"
		<< Exception::runerror;
	    diag = abstractVertex_->evaluate(scale,spinorInter,wavebar3_[ifo],swave3_)/gs;
	  }
	  else {
	    SpinorBarWaveFunction spinorBarInter = 
	      abstractIncomingVertex_->evaluate(scale,3,inpart.dataPtr(),wavebar3_[ifi],
					       gluon_[2*ig],inpart.mass());

	    if (wavebar3_[ifi].particle()->PDGName()!=spinorBarInter.particle()->PDGName())
	      throw Exception()
		<< wavebar3_[ifi].particle()->PDGName()  << " was changed to " 
		<< spinorBarInter.particle()->PDGName()  << " in FFSDecayer::threeBodyME"
		<< Exception::runerror;
	    diag = abstractVertex_->evaluate(scale,wave3_[ifo], spinorBarInter,swave3_)/gs;
	  }
	  for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	    (*ME[colourFlow[0][ix].first])(ifi, 0, ifo, ig) += 
	       colourFlow[0][ix].second*diag;
	  }
	}
	  
  	// radiation from outgoing fermion
  	if(decay[iferm]->dataPtr()->coloured()) {
  	  assert(abstractOutgoingVertexF_);
	  // ensure you get correct outgoing particle from first vertex
	  tcPDPtr off = decay[iferm]->dataPtr();
	  if(off->CC()) off = off->CC();

	  double gs   = abstractOutgoingVertexF_->strongCoupling(scale);	  	  
	  if (ferm) {	    
	    SpinorBarWaveFunction spinorBarInter = 
	      abstractOutgoingVertexF_->evaluate(scale,3,off,wavebar3_[ifo],
						gluon_[2*ig],decay[iferm]->mass());
	    
	    if(wavebar3_[ifo].particle()->PDGName()!=spinorBarInter.particle()->PDGName())
	      throw Exception()
		<< wavebar3_[ifo].particle()->PDGName() << " was changed to " 
		<< spinorBarInter.particle()->PDGName() << " in FFSDecayer::threeBodyME"
		<< Exception::runerror;
	    diag = abstractVertex_->evaluate(scale,wave3_[ifi],spinorBarInter,swave3_)/gs;
	  }
	  else {
	    SpinorWaveFunction spinorInter = 
	      abstractOutgoingVertexF_->evaluate(scale,3,off,wave3_[ifo],
						gluon_[2*ig],decay[iferm]->mass());
	      
	    if(wave3_[ifo].particle()->PDGName()!=spinorInter.particle()->PDGName())
	      throw Exception()
		<< wave3_[ifo].particle()->PDGName() << " was changed to " 
		<< spinorInter.particle()->PDGName() << " in FFSDecayer::threeBodyME"
		<< Exception::runerror;
	    diag = abstractVertex_->evaluate(scale,spinorInter,wavebar3_[ifi],swave3_)/gs;
	  }
	  for(unsigned int ix=0;ix<colourFlow[F].size();++ix) {
	    (*ME[colourFlow[F][ix].first])(ifi, 0, ifo, ig) += 
	      colourFlow[F][ix].second*diag;
	  }
  	}

  	// radiation from outgoing scalar
  	if(decay[iscal]->dataPtr()->coloured()) {
  	  assert(abstractOutgoingVertexS_);
	  // ensure you get correct ougoing particle from first vertex
	  tcPDPtr off = decay[iscal]->dataPtr();
	  if(off->CC()) off = off->CC();
	  
	  double gs = abstractOutgoingVertexS_->strongCoupling(scale);
	  ScalarWaveFunction  scalarInter = 
	    abstractOutgoingVertexS_->evaluate(scale,3,off,gluon_[2*ig],
					      swave3_,decay[iscal]->mass());
	    
	  if(swave3_.particle()->PDGName()!=scalarInter.particle()->PDGName())
	    throw Exception()
	      << swave3_    .particle()->PDGName() << " was changed to " 
	      << scalarInter.particle()->PDGName() << " in FFSDecayer::threeBodyME"
	      << Exception::runerror; 
	  if (ferm){
	    diag = abstractVertex_->evaluate(scale,wave3_[ifi],wavebar3_[ifo],scalarInter)/gs;
	  }
	  else {
	    diag = abstractVertex_->evaluate(scale,wave3_[ifo],wavebar3_[ifi],scalarInter)/gs;
	  }
	  for(unsigned int ix=0;ix<colourFlow[S].size();++ix) {
  	    (*ME[colourFlow[S][ix].first])(ifi, 0, ifo, ig) += 
	      colourFlow[S][ix].second*diag;
	  }
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
