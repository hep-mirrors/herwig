// -*- C++ -*-
//
// SVVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SVVDecayer class.
//

#include "SVVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr SVVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SVVDecayer::fullclone() const {
  return new_ptr(*this);
}

void SVVDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractVVSVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<VVSVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(inV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[1].at(inter));
  }
}

void SVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << perturbativeVertex_
     << incomingVertex_   << outgoingVertex1_
     << outgoingVertex2_;
}

void SVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> perturbativeVertex_
     >> incomingVertex_   >> outgoingVertex1_
     >> outgoingVertex2_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SVVDecayer,GeneralTwoBodyDecayer>
describeHerwigSVVDecayer("Herwig::SVVDecayer", "Herwig.so");

void SVVDecayer::Init() {

  static ClassDocumentation<SVVDecayer> documentation
    ("This implements the decay of a scalar to 2 vector bosons.");

}

double SVVDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector& decay, 
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1)));
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = decay[ix]->mass()==ZERO;
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&inpart),incoming);
    swave_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(rho_);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::
	constructSpinInfo(vectors_[ix],decay[ix],outgoing,true,photon[ix]);
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(vectors_[ix],decay[ix],outgoing,photon[ix]);
  
  
  Energy2 scale(sqr(inpart.mass()));
  unsigned int iv1,iv2;
  for(iv2 = 0; iv2 < 3; ++iv2) {
    if( photon[1] && iv2 == 1 ) ++iv2;
    for(iv1=0;iv1<3;++iv1) {
      if( photon[0] && iv1 == 1) ++iv1;
      (*ME())(0, iv1, iv2) = 0.;
      for(auto vert : vertex_)
	(*ME())(0, iv1, iv2) += vert->evaluate(scale,vectors_[0][iv1],
					       vectors_[1][iv2],swave_);
    }
  }
  double output = ME()->contract(rho_).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy SVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy2 scale(sqr(inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(scale, outa.first , 
				    outb.first, in);
    double mu1sq = sqr(outa.second/inpart.second);
    double mu2sq = sqr(outb.second/inpart.second);
    double m1pm2 = mu1sq + mu2sq;
    double me2(0.); 
    if( mu1sq > 0. && mu2sq > 0.)
      me2 = ( m1pm2*(m1pm2 - 2.) + 8.*mu1sq*mu2sq + 1.)/4./mu1sq/mu2sq;
    else if( mu1sq == 0. || mu2sq == 0. )
      me2 = 3.;
    else 
      me2 = 4.;
    
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = norm(perturbativeVertex_[0]->norm())*
      me2*pcm/(8*Constants::pi)/scale*UnitRemoval::E2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double SVVDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  if(meopt==Initialize) {
    // create scalar wavefunction for decaying particle
    ScalarWaveFunction::
      calculateWaveFunctions(rho3_,const_ptr_cast<tPPtr>(&inpart),incoming);
    swave3_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    VectorWaveFunction::
      constructSpinInfo(vectors3_[0],decay[0],outgoing,true,false);
    VectorWaveFunction::
      constructSpinInfo(vectors3_[1],decay[1],outgoing,true,false);
    VectorWaveFunction::
      constructSpinInfo(gluon_      ,decay[2],outgoing,true,false);
    return 0.;
  }
  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin0, PDT::Spin1,
								       PDT::Spin1, PDT::Spin1)));
  bool massless[2];
  for(unsigned int ix=0;ix<2;++ix)
    massless[ix] = decay[ix]->mass()!=ZERO;
  // create wavefunctions
  VectorWaveFunction::calculateWaveFunctions(vectors3_[0],decay[0],outgoing,massless[0]);
  VectorWaveFunction::calculateWaveFunctions(vectors3_[1],decay[1],outgoing,massless[1]);
  VectorWaveFunction::calculateWaveFunctions(gluon_      ,decay[2],outgoing,true);

  // gauge test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[2]->momentum(),
  					  decay[2]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif

  // get the outgoing vertices
  AbstractVVVVertexPtr outgoingVertex1;
  AbstractVVVVertexPtr outgoingVertex2;
  identifyVertices(inpart,decay, outgoingVertex1, outgoingVertex2,inter);

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif

   for(unsigned int iv1 = 0; iv1 < 3; ++iv1) {
     if(massless[0] && iv1==1) continue;
     for(unsigned int iv2 = 0; iv2 < 3; ++iv2) {
       if(massless[1] && iv2==1) continue;
       for(unsigned int ig = 0; ig < 2; ++ig) {
	 // radiation from the incoming vector
	 if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(incomingVertex_[inter]);	
	   ScalarWaveFunction scalarInter = 
	     incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),gluon_[2*ig],
					      swave3_,inpart.mass());
	
	   assert(swave3_.particle()->id()==scalarInter.particle()->id());
	
	   Complex diag = 0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectors3_[0][iv1],
				      vectors3_[1][iv2],scalarInter);
	   if(!couplingSet) {
	     gs = abs(incomingVertex_[inter]->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	     (*ME[colourFlow[0][ix].first])(0, iv1, iv2, ig) += 
	       colourFlow[0][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
	 }
	 // radiation from the 1st outgoing vector
	 if((decay[0]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (decay[0]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(outgoingVertex1);
	// ensure you get correct outgoing particle from first vertex
	   tcPDPtr off = decay[0]->dataPtr();
	   if(off->CC()) off = off->CC();
	   VectorWaveFunction vectorInter = 
	     outgoingVertex1->evaluate(scale,3,off,gluon_[2*ig],vectors3_[0][iv1],decay[0]->mass());
	   
	   assert(vectors3_[0][iv1].particle()->id()==vectorInter.particle()->id());
	 
	   Complex diag =0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectorInter,vectors3_[1][iv2],swave3_);
	   if(!couplingSet) {
	     gs = abs(outgoingVertex1->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	     (*ME[colourFlow[1][ix].first])(0, iv1, iv2, ig) += 
	       colourFlow[1][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
	 }
	 
	 if((decay[1]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (decay[1]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(outgoingVertex2);
	   // ensure you get correct outgoing particle from first vertex
	   tcPDPtr off = decay[1]->dataPtr();
	   if(off->CC()) off = off->CC();
	   VectorWaveFunction vectorInter = 
	     outgoingVertex2->evaluate(scale,3,off, gluon_[2*ig],vectors3_[1][iv2],decay[1]->mass());
	
	   assert(vectors3_[1][iv2].particle()->id()==vectorInter.particle()->id());
	
	   Complex diag = 0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectors3_[0][iv1],vectorInter,swave3_);
	   if(!couplingSet) {
	     gs = abs(outgoingVertex2->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	     (*ME[colourFlow[2][ix].first])(0, iv1, iv2, ig) += 
	       colourFlow[2][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
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
  // divide by alpha_(S,EM)
  output*=(4.*Constants::pi)/sqr(gs);
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

void SVVDecayer::identifyVertices(const Particle & inpart, const ParticleVector & decay, 
				  AbstractVVVVertexPtr & outgoingVertex1, 
				  AbstractVVVVertexPtr & outgoingVertex2,
				  ShowerInteraction inter) {
  if(inter==ShowerInteraction::QCD) {
    // work out which scalar each outgoing vertex corresponds to 
    // two outgoing vertices
    if( inpart.dataPtr()       ->iColour()==PDT::Colour0     &&
	((decay[0]->dataPtr()->iColour()==PDT::Colour3     &&
	  decay[1]->dataPtr()->iColour()==PDT::Colour3bar) ||
	 (decay[0]->dataPtr()->iColour()==PDT::Colour8     &&
	  decay[1]->dataPtr()->iColour()==PDT::Colour8))){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex2_[inter];
	outgoingVertex2 = outgoingVertex1_[inter];
      }
    }
    else if(inpart.dataPtr()       ->iColour()==PDT::Colour8 &&
	    decay[0]->dataPtr()->iColour()==PDT::Colour3 &&
	    decay[1]->dataPtr()->iColour()==PDT::Colour3bar){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex2_[inter];
	outgoingVertex2 = outgoingVertex1_[inter];
      }
    }
    
    // one outgoing vertex
    else if(inpart.dataPtr()->iColour()==PDT::Colour3){
      if(decay[0]->dataPtr()->iColour()==PDT::Colour3 &&  
	 decay[1]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertex1 = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertex1 = outgoingVertex2_[inter];
      }
      else if (decay[0]->dataPtr()->iColour()==PDT::Colour3 &&
	       decay[1]->dataPtr()->iColour()==PDT::Colour8){
	if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[1]->dataPtr()->id()))){
	  outgoingVertex1 = outgoingVertex2_[inter];
	  outgoingVertex2 = outgoingVertex1_[inter];
	}
	else {
	  outgoingVertex1 = outgoingVertex1_[inter];
	  outgoingVertex2 = outgoingVertex2_[inter];
	}
      }
    }
    else if(inpart.dataPtr()->iColour()==PDT::Colour3bar){
      if(decay[1]->dataPtr()->iColour()==PDT::Colour3bar &&  
	 decay[0]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertex2 = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (decay[0]->dataPtr()->iColour()==PDT::Colour8 &&
	       decay[1]->dataPtr()->iColour()==PDT::Colour3bar){
	if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[0]->dataPtr()->id()))){
	  outgoingVertex1 = outgoingVertex1_[inter];
	  outgoingVertex2 = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertex1 = outgoingVertex2_[inter];
	  outgoingVertex2 = outgoingVertex1_[inter];
	}
      }
    }
    
    if (! ((incomingVertex_[inter]  && (outgoingVertex1  || outgoingVertex2)) ||
	   ( outgoingVertex1 &&  outgoingVertex2)))
      throw Exception()
	<< "Invalid vertices for QCD radiation in SVV decay in SVVDecayer::identifyVertices"
	<< Exception::runerror;
  }
  else {
    if(decay[0]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[0]->dataPtr())))
	outgoingVertex1 = outgoingVertex1_[inter];
      else
	outgoingVertex1 = outgoingVertex2_[inter];
    }
    if(decay[1]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[1]->dataPtr())))
	outgoingVertex2 = outgoingVertex1_[inter];
      else
	outgoingVertex2 = outgoingVertex2_[inter];
    }
  }
}
