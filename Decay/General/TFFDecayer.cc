// -*- C++ -*-
//
// TFFDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TFFDecayer class.
//

#include "TFFDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
IBPtr TFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr TFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void TFFDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> &,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> fourV) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractFFTVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<FFTVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    fourPointVertex_[inter] = dynamic_ptr_cast<AbstractFFVTVertexPtr>(fourV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr> (outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr> (outV[1].at(inter));
  }
}

void TFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_          << perturbativeVertex_
     << outgoingVertex1_ << outgoingVertex2_
     << fourPointVertex_;
}

void TFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_          >> perturbativeVertex_
     >> outgoingVertex1_ >> outgoingVertex2_
     >> fourPointVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TFFDecayer,GeneralTwoBodyDecayer>
describeHerwigTFFDecayer("Herwig::TFFDecayer", "Herwig.so");

void TFFDecayer::Init() {

  static ClassDocumentation<TFFDecayer> documentation
    ("The TFFDecayer class implements the decay of a tensor particle "
     "to 2 fermions ");
  
}

double TFFDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  unsigned int iferm(0),ianti(1);
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin2,PDT::Spin1Half,PDT::Spin1Half)));
  if(decay[0]->id()>=0) swap(iferm,ianti);
  if(meopt==Initialize) {
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&inpart),
			     incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  if(meopt==Terminate) {
    TensorWaveFunction::
      constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&inpart),
			incoming,true,false);
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
  unsigned int thel,fhel,ahel;
  for(thel=0;thel<5;++thel) {
    for(fhel=0;fhel<2;++fhel) {
      for(ahel=0;ahel<2;++ahel) {
	if(iferm > ianti) {
	  (*ME())(thel,fhel,ahel) = 0.;
	  for(auto vert : vertex_)
	    (*ME())(thel,fhel,ahel) += 
	      vert->evaluate(scale,wave_[ahel],
			     wavebar_[fhel],tensors_[thel]);
	}
	else {
	  (*ME())(thel,ahel,fhel) = 0.;
	  for(auto vert : vertex_)
	    (*ME())(thel,ahel,fhel) += 
	      vert->evaluate(scale,wave_[ahel],
			     wavebar_[fhel],tensors_[thel]);
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

Energy TFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy2 scale = sqr(inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(scale, in, outa.first, outb.first);
    double musq = sqr(outa.second/inpart.second);
    double b = sqrt(1- 4.*musq);
    double me2 = b*b*(5-2*b*b)*scale/120.*UnitRemoval::InvE2;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = norm(perturbativeVertex_[0]->norm())*me2*pcm/(8.*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}


double TFFDecayer::threeBodyME(const int , const Particle & inpart,
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
    // create tensor wavefunction for decaying particle
    TensorWaveFunction::
      calculateWaveFunctions(tensors3_, rho3_, const_ptr_cast<tPPtr>(&inpart), incoming, false);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    TensorWaveFunction::
      constructSpinInfo(tensors3_, const_ptr_cast<tPPtr>(&inpart),incoming,true, false);
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

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin2,     PDT::Spin1Half,
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
  
  if (! (outgoingVertex1_[inter] && outgoingVertex2_[inter]))
    throw Exception()
      << "Invalid vertices for QCD radiation in TFF decay in TFFDecayer::threeBodyME"
      << Exception::runerror;

  // identify fermion and/or anti-fermion vertex
  AbstractFFVVertexPtr outgoingVertexF = outgoingVertex1_[inter];
  AbstractFFVVertexPtr outgoingVertexA = outgoingVertex2_[inter];

  if(outgoingVertex1_[inter]!=outgoingVertex2_[inter] &&
     outgoingVertex1_[inter]->isIncoming(getParticleData(decay[ianti]->id())))
    swap (outgoingVertexF, outgoingVertexA);  
  
  if(! (inpart.dataPtr()->iColour()==PDT::Colour0)){
    throw Exception()
      << "Invalid vertices for QCD radiation in TFF decay in TFFDecayer::threeBodyME"
      << Exception::runerror;
  }

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int it = 0; it < 5; ++it) {  
    for(unsigned int ifm = 0; ifm < 2; ++ifm) {
      for(unsigned int ia = 0; ia < 2; ++ia) {
	for(unsigned int ig = 0; ig < 2; ++ig) {

	  // radiation from outgoing fermion
	  if((decay[iferm]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[iferm]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertexF);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[iferm]->dataPtr();
	    if(off->CC()) off = off->CC();
	    SpinorBarWaveFunction interS = 
	      outgoingVertexF->evaluate(scale,3,off,wavebar3_[ifm],
					gluon_[2*ig],decay[iferm]->mass());
	  
	    assert(wavebar3_[ifm].particle()->id()==interS.particle()->id());

	    Complex diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,wave3_[ia], interS,tensors3_[it]);
	    if(!couplingSet) {
	      gs = abs(outgoingVertexF->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	      (*ME[colourFlow[1][ix].first])(it, ifm, ia, ig) += 
		colourFlow[1][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }

	  // radiation from outgoing antifermion
	  if((decay[ianti]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[ianti]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertexA);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[ianti]->dataPtr();
	    if(off->CC()) off = off->CC();
	    SpinorWaveFunction  interS = 
	      outgoingVertexA->evaluate(scale,3,off,wave3_[ia],
					gluon_[2*ig],decay[ianti]->mass());
	    
	    assert(wave3_[ia].particle()->id()==interS.particle()->id());

	    Complex diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,interS,wavebar3_[ifm],tensors3_[it]);
	    if(!couplingSet) {
	      gs = abs(outgoingVertexA->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	      (*ME[colourFlow[2][ix].first])(it, ifm, ia, ig) += 
		colourFlow[2][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }

	  // radiation from 4 point vertex
	  if (fourPointVertex_[inter]) {
	    Complex diag = fourPointVertex_[inter]->evaluate(scale, wave3_[ia], wavebar3_[ifm],
							     gluon_[2*ig], tensors3_[it]);
	    for(unsigned int ix=0;ix<colourFlow[3].size();++ix) {
	      (*ME[colourFlow[3][ix].first])(it, ifm, ia, ig) += 
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
  // divide by alpha_(s,em)
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

