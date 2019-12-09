// -*- C++ -*-
//
// FFVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFVDecayer class.
//

#include "FFVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr FFVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FFVDecayer::fullclone() const {
  return new_ptr(*this);
}

void FFVDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_             .push_back(dynamic_ptr_cast<AbstractFFVVertexPtr>(vert));
    perturbativeVertex_ .push_back(dynamic_ptr_cast<FFVVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(inV.at(inter));
    if(outV[0].at(inter)) {
      if (outV[0].at(inter)->getName()==VertexType::FFV)
	outgoingVertexF_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[0].at(inter));
      else
	outgoingVertexV_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[0].at(inter));
    }
    if(outV[1].at(inter)) {
      if (outV[1].at(inter)->getName()==VertexType::FFV)
	outgoingVertexF_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[1].at(inter));
      else
	outgoingVertexV_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[1].at(inter));
    }
  }
}

void FFVDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_           << perturbativeVertex_
     << incomingVertex_   << outgoingVertexF_
     << outgoingVertexV_;
}

void FFVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_           >> perturbativeVertex_
     >> incomingVertex_   >> outgoingVertexF_
     >> outgoingVertexV_;
}

double FFVDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay, 
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  // type of process
  int itype[2];
  if(inpart.dataPtr()->CC())    itype[0] = inpart.id() > 0 ? 0 : 1;
  else                          itype[0] = 2;
  if(decay[0]->dataPtr()->CC()) itype[1] = decay[0]->id() > 0 ? 0 : 1;
  else                          itype[1] = 2;  
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
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
    // fix rho if no correlations
    fixRho(rho_);
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
    VectorWaveFunction::
      constructSpinInfo(vector_,decay[1],outgoing,true,false);
  }
  Energy2 scale(sqr(inpart.mass()));
  if(ferm)
    SpinorBarWaveFunction::
      calculateWaveFunctions(wavebar_,decay[0],outgoing);
  else
    SpinorWaveFunction::
      calculateWaveFunctions(wave_   ,decay[0],outgoing);
  bool massless = decay[1]->dataPtr()->mass()==ZERO;
  VectorWaveFunction::
    calculateWaveFunctions(vector_,decay[1],outgoing,massless);
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      for(unsigned int vhel = 0; vhel < 3; ++vhel) {
	if(massless && vhel == 1) ++vhel;
	if(ferm)
	  (*ME())(if1, if2,vhel) = 0.;
	else
	  (*ME())(if2, if1, vhel) = 0.;
	for(auto vertex : vertex_) {
	  if(ferm)
	    (*ME())(if1, if2,vhel) += 
	      vertex->evaluate(scale,wave_[if1],wavebar_[if2],vector_[vhel]);
	  else
	    (*ME())(if2, if1, vhel) += 
	      vertex->evaluate(scale,wave_[if1],wavebar_[if2],vector_[vhel]);
	}
      }
    }
  }
  double output=(ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy FFVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    double mu1(outa.second/inpart.second),mu2(outb.second/inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if( outa.first->iSpin() == PDT::Spin1Half)
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in,
				       outa.first, outb.first);
    else {
      swap(mu1,mu2);
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second),in,
				       outb.first,outa.first);
    }
    Complex cl(perturbativeVertex_[0]->left()),cr(perturbativeVertex_[0]->right());
    double me2(0.);
    if( mu2 > 0. ) {
      me2 = (norm(cl) + norm(cr))*(1. + sqr(mu1*mu2) + sqr(mu2) 
				   - 2.*sqr(mu1) - 2.*sqr(mu2*mu2) 
				   +  sqr(mu1*mu1))
	- 6.*mu1*sqr(mu2)*(conj(cl)*cr + conj(cr)*cl).real();
      me2 /= sqr(mu2);
    }
    else
      me2 = 2.*( (norm(cl) + norm(cr))*(sqr(mu1) + 1.) 
		 - 4.*mu1*(conj(cl)*cr + conj(cr)*cl).real() );
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					outb.second);
    Energy output = norm(perturbativeVertex_[0]->norm())*me2*pcm/16./Constants::pi; 
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer 
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FFVDecayer,GeneralTwoBodyDecayer>
describeHerwigFFVDecayer("Herwig::FFVDecayer", "Herwig.so");

void FFVDecayer::Init() {

  static ClassDocumentation<FFVDecayer> documentation
    ("The FFVDecayer class implements the decay of a fermion to a fermion and a vector boson");

}

double  FFVDecayer::threeBodyME(const int , const Particle & inpart,
				const ParticleVector & decay,
				ShowerInteraction inter, MEOption meopt) {

  int iferm (0), ivect (1), iglu (2);
  // get location of outgoing lepton/vector
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin1Half) swap(iferm,ivect);
  // work out whether inpart is a fermion or antifermion
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[iferm]->dataPtr()->CC()) itype[1] = decay[iferm]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;

  bool ferm(itype[0] == 0 || itype[1] == 0 || 
	   (itype[0] == 2 && itype[1] == 2 && decay[ivect]->id() < 0));  

  // no emissions from massive vectors
  bool massless = decay[ivect]->dataPtr()->mass()==ZERO;

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
    VectorWaveFunction::constructSpinInfo(vector3_, decay[ivect],outgoing,true,massless);
    VectorWaveFunction::constructSpinInfo(gluon_,   decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calulate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half, PDT::Spin1Half,
								       PDT::Spin1,     PDT::Spin1)));

  // create wavefunctions
  if (ferm)  SpinorBarWaveFunction::calculateWaveFunctions(wavebar3_, decay[iferm],outgoing);
  else       SpinorWaveFunction::   calculateWaveFunctions(wave3_   , decay[iferm],outgoing);
  
  VectorWaveFunction::calculateWaveFunctions(vector3_, decay[ivect],outgoing,massless);
  VectorWaveFunction::calculateWaveFunctions(gluon_,   decay[iglu ],outgoing,true );

  // gauge invariance test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  					  decay[iglu ]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif

  if (! ((incomingVertex_[inter]  && (outgoingVertexF_[inter]  || outgoingVertexV_[inter])) ||
	 (outgoingVertexF_[inter] &&  outgoingVertexV_[inter])))
    throw Exception()
      << "Invalid vertices for QCD radiation in FFV decay in FFVDecayer::threeBodyME"
      << Exception::runerror;


  // sort out colour flows
  int F(1), V(2);
  if (decay[iferm]->dataPtr()->iColour()==PDT::Colour3bar && 
      decay[ivect]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,V);
  else if (decay[ivect]->dataPtr()->iColour()==PDT::Colour3 && 
	   decay[iferm]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,V);

  Complex diag;
  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int ifi = 0; ifi < 2; ++ifi) {
    for(unsigned int ifo = 0; ifo < 2; ++ifo) {
      for(unsigned int iv = 0; iv < 3; ++iv) {
	for(unsigned int ig = 0; ig < 2; ++ig) {
	  // radiation from the incoming fermion
	  if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(incomingVertex_[inter]);
	    if (ferm){
	      SpinorWaveFunction spinorInter =
		incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),wave3_[ifi],
						  gluon_[2*ig],inpart.mass());
	      
	      assert(wave3_[ifi].particle()->id()==spinorInter.particle()->id());
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,spinorInter,wavebar3_[ifo],vector3_[iv]);
	    }
	    else {
	      SpinorBarWaveFunction spinorBarInter = 
		incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),wavebar3_[ifi],
						  gluon_[2*ig],inpart.mass());
	      
	      assert(wavebar3_[ifi].particle()->id()==spinorBarInter.particle()->id());
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,wave3_[ifo], spinorBarInter,vector3_[iv]);
	    }
	    if(!couplingSet) {
	      gs = abs(incomingVertex_[inter]->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	      (*ME[colourFlow[0][ix].first])(ifi, ifo, iv, ig) += 
		colourFlow[0][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }
	  // radiation from outgoing fermion
	  if((decay[iferm]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[iferm]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertexF_[inter]);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[iferm]->dataPtr();
	    if(off->CC()) off = off->CC(); 
	    if (ferm) {	    
	      SpinorBarWaveFunction spinorBarInter = 
		outgoingVertexF_[inter]->evaluate(scale,3,off,wavebar3_[ifo],
						  gluon_[2*ig],decay[iferm]->mass());
	      
	      assert(wavebar3_[ifo].particle()->id()==spinorBarInter.particle()->id());
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,wave3_[ifi],spinorBarInter,vector3_[iv]);
	    }
	    else {
	      SpinorWaveFunction spinorInter = 
		outgoingVertexF_[inter]->evaluate(scale,3,off,wave3_[ifo],
						  gluon_[2*ig],decay[iferm]->mass());
		
	      assert(wave3_[ifo].particle()->id()==spinorInter.particle()->id());
	      
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,spinorInter,wavebar3_[ifi],vector3_[iv]);
	    }
	    if(!couplingSet) {
	      gs = abs(outgoingVertexF_[inter]->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[F].size();++ix) {
	      (*ME[colourFlow[F][ix].first])(ifi, ifo, iv, ig) += 
		 colourFlow[F][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }
	  
	  // radiation from outgoing vector
	  if((decay[ivect]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[ivect]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertexV_[inter]);
	    // ensure you get correct ougoing particle from first vertex
	    tcPDPtr off = decay[ivect]->dataPtr();
	    if(off->CC()) off = off->CC();
	    VectorWaveFunction  vectorInter = 
	      outgoingVertexV_[inter]->evaluate(scale,3,off,gluon_[2*ig],
						vector3_[iv],decay[ivect]->mass());
	    
	    assert(vector3_[iv].particle()->id()==vectorInter.particle()->id());
	    if (ferm) {
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,wave3_[ifi],wavebar3_[ifo],vectorInter);
	    }
	    else {
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,wave3_[ifo],wavebar3_[ifi],vectorInter);
	    }
	    if(!couplingSet) {
	      gs = abs(outgoingVertexV_[inter]->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[V].size();++ix) {
	      (*ME[colourFlow[V][ix].first])(ifi, ifo, iv, ig) += 
		colourFlow[V][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }
	}
	if(massless) ++iv;
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
  // divide by alpha(S,eM)
  output *= (4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
  double ratio = output/total;
  if(abs(ratio)>1e-20) {
    generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
    for(unsigned int ix=0;ix<decay.size();++ix)
      generator()->log() << *decay[ix] << "\n";
    generator()->log() << "Test of gauge invariance " << ratio << "\n";
  }
#endif
  // return the answer
  return output;
}
