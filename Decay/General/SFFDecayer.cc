// -*- C++ -*-
//
// SFFDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SFFDecayer class.
//

#include "SFFDecayer.h"
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

IBPtr SFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void SFFDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_             .push_back(dynamic_ptr_cast<AbstractFFSVertexPtr>(vert));
    perturbativeVertex_ .push_back(dynamic_ptr_cast<FFSVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(inV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[1].at(inter));
  }
}

void SFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_           << perturbativeVertex_ 
     << incomingVertex_   << outgoingVertex1_
     << outgoingVertex2_;
}

void SFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_           >> perturbativeVertex_ 
     >> incomingVertex_   >> outgoingVertex1_
     >> outgoingVertex2_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SFFDecayer,GeneralTwoBodyDecayer>
describeHerwigSFFDecayer("Herwig::SFFDecayer", "Herwig.so");

void SFFDecayer::Init() {

  static ClassDocumentation<SFFDecayer> documentation
    ("This class implements to decay of a scalar to 2 fermions");

}

double SFFDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
  // work out which is the fermion and antifermion
  int iferm(1),ianti(0);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0||itype[1]==1||(itype[0]==2&&itype[1]==2)) swap(iferm,ianti);

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
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
    return 0.;
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(wavebar_,decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(wave_   ,decay[ianti],outgoing);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int ifm = 0; ifm < 2; ++ifm){
    for(unsigned int ia = 0; ia < 2; ++ia) {
      if(iferm > ianti) (*ME())(0, ia, ifm) = 0.;
      else              (*ME())(0, ifm, ia) = 0.;
      for(auto vert : vertex_) {
	if(iferm > ianti){
	  (*ME())(0, ia, ifm) += vert->evaluate(scale,wave_[ia],
						wavebar_[ifm],swave_);
	}
	else {
	  (*ME())(0, ifm, ia) += vert->evaluate(scale,wave_[ia],
						wavebar_[ifm],swave_);
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

Energy SFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(sqr(inpart.second), outb.first, outa.first,
				     in);
    double mu1(outa.second/inpart.second),mu2(outb.second/inpart.second);
    double c2 = norm(perturbativeVertex_[0]->norm());
    Complex al(perturbativeVertex_[0]->left()), ar(perturbativeVertex_[0]->right());
    double me2 = -c2*( (norm(al) + norm(ar))*( sqr(mu1) + sqr(mu2) - 1.)
		       + 2.*(ar*conj(al) + al*conj(ar)).real()*mu1*mu2 );
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					outb.second);
    Energy output = me2*pcm/(8*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double SFFDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  // work out which is the fermion and antifermion
  int ianti(0), iferm(1), iglu(2);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0 && itype[1]!=0) swap(iferm, ianti);
  if(itype[0]==2 && itype[1]==1) swap(iferm, ianti);
  if(itype[0]==0 && itype[1]==0 && decay[0]->dataPtr()->id()<decay[1]->dataPtr()->id()) 
    swap(iferm, ianti);
  if(itype[0]==1 && itype[1]==1 && decay[0]->dataPtr()->id()<decay[1]->dataPtr()->id()) 
    swap(iferm, ianti);
  
  if(meopt==Initialize) {
    // create scalar wavefunction for decaying particle
    ScalarWaveFunction::
      calculateWaveFunctions(rho3_,const_ptr_cast<tPPtr>(&inpart),incoming);
    swave3_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar3_ ,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(wave3_    ,decay[ianti],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(gluon_    ,decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,     PDT::Spin1Half,
								       PDT::Spin1Half, PDT::Spin1)));
  // create wavefunctions
  SpinorBarWaveFunction::
    calculateWaveFunctions(wavebar3_, decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(wave3_   , decay[ianti],outgoing);
  VectorWaveFunction::
    calculateWaveFunctions(gluon_   , decay[iglu ],outgoing,true);

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

  // identify fermion and/or anti-fermion vertex
  AbstractFFVVertexPtr outgoingVertexF;
  AbstractFFVVertexPtr outgoingVertexA;
  identifyVertices(iferm, ianti, inpart, decay, outgoingVertexF, outgoingVertexA,
		   inter);

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);

  Energy2 scale(sqr(inpart.mass()));
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int ifm = 0; ifm < 2; ++ifm) {
    for(unsigned int ia = 0; ia < 2; ++ia) {
      for(unsigned int ig = 0; ig < 2; ++ig) {
	// radiation from the incoming scalar
	if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	   (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	  assert(incomingVertex_[inter]);

	  ScalarWaveFunction scalarInter = 
	    incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),
					      gluon_[2*ig],swave3_,inpart.mass());

	  assert(swave3_.particle()->id()==scalarInter.particle()->id());

	  if(!couplingSet) {
	    gs = abs(incomingVertex_[inter]->norm());
	    couplingSet = true;
	  }
	  Complex diag = 0.;
	  for(auto vertex : vertex_)
	    diag += vertex->evaluate(scale,wave3_[ia],
				     wavebar3_[ifm],scalarInter);
	  for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	    (*ME[colourFlow[0][ix].first])(0, ia, ifm, ig) += 
	      colourFlow[0][ix].second*diag;
	  }
#ifdef GAUGE_CHECK
	  total+=norm(diag);
#endif
	}
	// radiation from outgoing fermion
	if((decay[iferm]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	   (decay[iferm]->dataPtr()->charged()  && inter==ShowerInteraction::QED)) {
	  assert(outgoingVertexF);
	  // ensure you get correct outgoing particle from first vertex
	  tcPDPtr off = decay[iferm]->dataPtr();
	  if(off->CC()) off = off->CC();
	  SpinorBarWaveFunction interS = 
	    outgoingVertexF->evaluate(scale,3,off,wavebar3_[ifm],
				      gluon_[2*ig],decay[iferm]->mass());
	  
	  assert(wavebar3_[ifm].particle()->id()==interS.particle()->id());
	  
	  if(!couplingSet) {
	    gs = abs(outgoingVertexF->norm());
	    couplingSet = true;
	  }
	  Complex diag = 0.;
	  for(auto vertex : vertex_)
	    diag += vertex->evaluate(scale,wave3_[ia], interS,swave3_);
	  for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	    (*ME[colourFlow[1][ix].first])(0, ia, ifm, ig) += 
	      colourFlow[1][ix].second*diag;
	  }
#ifdef GAUGE_CHECK
	  total+=norm(diag);
#endif
	}
	
	// radiation from outgoing antifermion
	if((decay[ianti]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	   (decay[ianti]->dataPtr()->charged()  && inter==ShowerInteraction::QED)) {
	  assert(outgoingVertexA);
	  // ensure you get correct outgoing particle from first vertex
	  tcPDPtr off = decay[ianti]->dataPtr();
	  if(off->CC()) off = off->CC();
	  SpinorWaveFunction  interS = 
	    outgoingVertexA->evaluate(scale,3,off,wave3_[ia],
				      gluon_[2*ig],decay[ianti]->mass());
	  
	  assert(wave3_[ia].particle()->id()==interS.particle()->id());

	  if(!couplingSet) {
	    gs = abs(outgoingVertexA->norm());
	    couplingSet = true;
	  }
	  Complex diag = 0.;
	  for(auto vertex : vertex_)
	    diag += vertex->evaluate(scale,interS,wavebar3_[ifm],swave3_);
	  for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	    (*ME[colourFlow[2][ix].first])(0, ia, ifm, ig) += 
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
  // divide by alpha(S,EM)
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

void SFFDecayer::identifyVertices(const int iferm, const int ianti,
				  const Particle & inpart, const ParticleVector & decay, 
				  AbstractFFVVertexPtr & outgoingVertexF, 
				  AbstractFFVVertexPtr & outgoingVertexA,
				  ShowerInteraction inter) {
  // QCD
  if(inter==ShowerInteraction::QCD) {
    // work out which fermion each outgoing vertex corresponds to 
    // two outgoing vertices
    if( inpart.dataPtr()       ->iColour()==PDT::Colour0     &&
	((decay[iferm]->dataPtr()->iColour()==PDT::Colour3     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) ||
	 (decay[iferm]->dataPtr()->iColour()==PDT::Colour8     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour8))) {
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]) {
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	outgoingVertexF = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    else if(inpart.dataPtr()       ->iColour()==PDT::Colour8 &&
	    decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&
	    decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) {
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]) {
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	outgoingVertexF = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    
    // one outgoing vertex
    else if(inpart.dataPtr()->iColour()==PDT::Colour3){
      if(decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&  
	 decay[ianti]->dataPtr()->iColour()==PDT::Colour0){
      if     (outgoingVertex1_[inter]) outgoingVertexF = outgoingVertex1_[inter];
      else if(outgoingVertex2_[inter]) outgoingVertexF = outgoingVertex2_[inter];
      }
      else if (decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&
	       decay[ianti]->dataPtr()->iColour()==PDT::Colour8) {
	if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[ianti]->dataPtr()))) {
	  outgoingVertexF = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
	else {
	  outgoingVertexF = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
      }
      else if(decay[iferm]->dataPtr()->iColour()==PDT::Colour3bar &&
	      decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) {
	if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	  outgoingVertexF = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertexF = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
      }
    }
    else if(inpart.dataPtr()->iColour()==PDT::Colour3bar){
      if(decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
	 decay[iferm]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertexA = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (decay[iferm]->dataPtr()->iColour()==PDT::Colour8 &&
	       decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
	if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	  outgoingVertexF = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertexF = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
      }
      else if(decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&
	      decay[ianti]->dataPtr()->iColour()==PDT::Colour3) {
	if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	  outgoingVertexF = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertexF = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
      }
    }
    else if(inpart.dataPtr()->iColour()==PDT::Colour6 ||
	    inpart.dataPtr()->iColour()==PDT::Colour6bar) {
      if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else {
	outgoingVertexF = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
  
    if (! ((incomingVertex_[inter]  && (outgoingVertexF  || outgoingVertexA)) ||
	   ( outgoingVertexF &&  outgoingVertexA))) {
      throw Exception()
	<< "Invalid vertices for QCD radiation in SFF decay in SFFDecayer::identifyVertices"
	<< Exception::runerror;
    }
  }
  // QED
  else {
    if(decay[iferm]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr())))
	outgoingVertexF = outgoingVertex1_[inter];
      else
	outgoingVertexF = outgoingVertex2_[inter];
    }
    if(decay[ianti]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[ianti]->dataPtr())))
	outgoingVertexA = outgoingVertex1_[inter];
      else
	outgoingVertexA = outgoingVertex2_[inter];
    }
  }
}
