// -*- C++ -*-
//
// VVVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVVDecayer class.
//

#include "VVVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr VVVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VVVDecayer::fullclone() const {
  return new_ptr(*this);
}

void VVVDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> fourV) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_             .push_back(dynamic_ptr_cast<AbstractVVVVertexPtr>(vert));
    perturbativeVertex_ .push_back(dynamic_ptr_cast<VVVVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter]  = dynamic_ptr_cast<AbstractVVVVertexPtr>(inV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[1].at(inter));
    fourPointVertex_[inter] = dynamic_ptr_cast<AbstractVVVVVertexPtr>(fourV.at(inter));
  }
}

void VVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << perturbativeVertex_
     << incomingVertex_  << outgoingVertex1_
     << outgoingVertex2_ << fourPointVertex_;
}

void VVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> perturbativeVertex_
     >> incomingVertex_   >> outgoingVertex1_
     >> outgoingVertex2_ >> fourPointVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VVVDecayer,GeneralTwoBodyDecayer>
describeHerwigVVVDecayer("Herwig::VVVDecayer", "Herwig.so");

void VVVDecayer::Init() {

  static ClassDocumentation<VVVDecayer> documentation
    ("The VVVDecayer class implements the decay of a vector boson "
     "into 2 vector bosons");

}

double VVVDecayer::me2(const int , const Particle & inpart,
                       const ParticleVector & decay,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  bool massless[2];
  for(unsigned int ix=0;ix<2;++ix) 
    massless[ix] = (decay[ix]->id()==ParticleID::gamma ||
		    decay[ix]->id()==ParticleID::g);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::
	constructSpinInfo(vectors_[ix+1],decay[ix],outgoing,true,massless[ix]);
    return 0.;
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(vectors_[ix+1],decay[ix],outgoing,massless[ix]);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int iv3=0;iv3<3;++iv3) {
    for(unsigned int iv2=0;iv2<3;++iv2) {
      for(unsigned int iv1=0;iv1<3;++iv1) {
	(*ME())(iv1,iv2,iv3) = 0.;
	for(auto vert : vertex_) {
	  (*ME())(iv1,iv2,iv3) += vert->
	    evaluate(scale,vectors_[1][iv2],vectors_[2][iv3],vectors_[0][iv1]);
	}
      }
    }
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy VVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in,
					outa.first, outb.first);
    double mu1(outa.second/inpart.second), mu1sq(sqr(mu1)),
      mu2(outb.second/inpart.second), mu2sq(sqr(mu2));
    double vn = norm(perturbativeVertex_[0]->norm());
    if(vn == ZERO || mu1sq == ZERO || mu2sq == ZERO) return ZERO;
    double me2 = 
      (mu1 - mu2 - 1.)*(mu1 - mu2 + 1.)*(mu1 + mu2 - 1.)*(mu1 + mu2 + 1.)
      * (sqr(mu1sq) + sqr(mu2sq) + 10.*(mu1sq*mu2sq + mu1sq + mu2sq) + 1.)
      /4./mu1sq/mu2sq;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy pWidth = vn*me2*pcm/24./Constants::pi;
    // colour factor
    pWidth *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return pWidth;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double VVVDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  if(meopt==Initialize) {
    // create vector wavefunction for decaying particle
    VectorWaveFunction::calculateWaveFunctions(vector3_, rho3_, const_ptr_cast<tPPtr>(&inpart), 
					       incoming, false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::
      constructSpinInfo(vector3_ ,const_ptr_cast<tPPtr>(&inpart),outgoing,true,false);
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
  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1, PDT::Spin1,
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

   for(unsigned int iv0 = 0; iv0 < 3; ++iv0) {
     for(unsigned int iv1 = 0; iv1 < 3; ++iv1) {
       if(massless[0] && iv1==1) continue;
       for(unsigned int iv2 = 0; iv2 < 3; ++iv2) {
	 if(massless[1] && iv2==1) continue;
	 for(unsigned int ig = 0; ig < 2; ++ig) {
	   // radiation from the incoming vector
	   if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	      (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	     assert(incomingVertex_[inter]);
	     VectorWaveFunction vectorInter = 
	       incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),vector3_[iv0],
						gluon_[2*ig],inpart.mass());
	     
	     assert(vector3_[iv0].particle()->id()==vectorInter.particle()->id());
	     
	     Complex diag = 0.;
	     for(auto vertex : vertex_)
	       diag += vertex->evaluate(scale,vectorInter,vectors3_[0][iv1],
					vectors3_[1][iv2]);
	     if(!couplingSet) {
	       gs = abs(incomingVertex_[inter]->norm());
	       couplingSet = true;
	     }
	     for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	       (*ME[colourFlow[0][ix].first])(iv0, iv1, iv2, ig) += 
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
	     
	     Complex diag = 0.;
	     for(auto vertex : vertex_)
	       diag += vertex->evaluate(scale,vector3_[iv0],vectorInter,vectors3_[1][iv2]);
	     if(!couplingSet) {
	       gs = abs(outgoingVertex1->norm());
	       couplingSet = true;
	     }
	     for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	       (*ME[colourFlow[1][ix].first])(iv0, iv1, iv2, ig) += 
		 colourFlow[1][ix].second*diag;
	     }
#ifdef GAUGE_CHECK
	     total+=norm(diag);
#endif
	   }
	   // radiation from the second outgoing vector
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
	       diag += vertex->evaluate(scale,vector3_[iv0],vectors3_[0][iv1],vectorInter);
	     if(!couplingSet) {
	       gs = abs(outgoingVertex2->norm());
	       couplingSet = true;
	     }
	     for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	       (*ME[colourFlow[2][ix].first])(iv0, iv1, iv2, ig) += 
		 colourFlow[2][ix].second*diag;
	     }
#ifdef GAUGE_CHECK
	     total+=norm(diag);
#endif
	   }
	   // 4-point vertex
	   if (fourPointVertex_[inter]) {
	     Complex diag = fourPointVertex_[inter]->evaluate(scale,0, vector3_[iv0],
							      vectors3_[0][iv1],
							      vectors3_[1][iv2],gluon_[2*ig]);
	     for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	       (*ME[colourFlow[3][ix].first])(iv0, iv1, iv2, ig) += 
		 colourFlow[3][ix].second*diag;
	     }
#ifdef GAUGE_CHECK
	     total+=norm(diag);
#endif
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

void VVVDecayer::identifyVertices(const Particle & inpart, const ParticleVector & decay, 
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
	<< "Invalid vertices for QCD radiation in VVV decay in VVVDecayer::identifyVertices"
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
