// -*- C++ -*-
//
// TVVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TVVDecayer class.
//

#include "TVVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr TVVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr TVVDecayer::fullclone() const {
  return new_ptr(*this);
}

void TVVDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> &,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> fourV) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractVVTVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<VVTVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    fourPointVertex_[inter] = dynamic_ptr_cast<AbstractVVVTVertexPtr>(fourV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr> (outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr> (outV[1].at(inter));
  }
}

void TVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_          << perturbativeVertex_
     << outgoingVertex1_ << outgoingVertex2_
     << fourPointVertex_;
}

void TVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_          >> perturbativeVertex_
     >> outgoingVertex1_ >> outgoingVertex2_
     >> fourPointVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TVVDecayer,GeneralTwoBodyDecayer>
describeHerwigTVVDecayer("Herwig::TVVDecayer", "Herwig.so");

void TVVDecayer::Init() {

  static ClassDocumentation<TVVDecayer> documentation
    ("This class implements the decay of a tensor to 2 vector bosons");

}

double TVVDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin1)));
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = decay[ix]->mass()==ZERO;
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
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::
	constructSpinInfo(vectors_[ix],decay[ix],outgoing,true,photon[ix]);
    return 0.;
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(vectors_[ix],decay[ix],outgoing,photon[ix]);
  Energy2 scale(sqr(inpart.mass()));
  unsigned int thel,v1hel,v2hel;
  for(thel=0;thel<5;++thel) {
    for(v1hel=0;v1hel<3;++v1hel) {
      for(v2hel=0;v2hel<3;++v2hel) {
	(*ME())(thel,v1hel,v2hel) = 0.;
	  for(auto vert : vertex_)
	    (*ME())(thel,v1hel,v2hel) += vert->evaluate(scale,
							vectors_[0][v1hel],
							vectors_[1][v2hel],
							tensors_[thel]);
	if(photon[1]) ++v2hel;
      }
      if(photon[0]) ++v1hel;
    }
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}
  
Energy TVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy2 scale(sqr(inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(scale, outa.first, outb.first, in);
    double mu2 = sqr(outa.second/inpart.second);
    double b = sqrt(1 - 4.*mu2);
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy2 me2;
    if(outa.second > ZERO && outb.second > ZERO)
      me2 = scale*(30 - 20.*b*b + 3.*pow(b,4))/120.; 
    else 
      me2 = scale/10.;
    
    Energy output = norm(perturbativeVertex_[0]->norm())*me2*pcm
      /(8.*Constants::pi)*UnitRemoval::InvE2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double TVVDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  bool massless[2];
  for(unsigned int ix=0;ix<2;++ix)
    massless[ix] = decay[ix]->mass()==ZERO;
  int iglu(2);  
  if(meopt==Initialize) {
    // create tensor wavefunction for decaying particle
    TensorWaveFunction::
      calculateWaveFunctions(tensors3_, rho3_, const_ptr_cast<tPPtr>(&inpart), incoming, false);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    TensorWaveFunction::
      constructSpinInfo(tensors3_, const_ptr_cast<tPPtr>(&inpart),incoming,true, false);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::
	constructSpinInfo(vectors3_[ix],decay[ix   ],outgoing,true, massless[ix]);
    VectorWaveFunction::
        constructSpinInfo(gluon_       ,decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin2, PDT::Spin1,
								       PDT::Spin1, PDT::Spin1)));
  // create wavefunctions
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(vectors3_[ix],decay[ix   ],outgoing,massless[ix]);
  VectorWaveFunction::
      calculateWaveFunctions(gluon_       ,decay[iglu ],outgoing,true);

  // gauge test
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
  
  // work out which vector each outgoing vertex corresponds to 
  if(outgoingVertex1_[inter]!=outgoingVertex2_[inter] &&
     outgoingVertex1_[inter]->isIncoming(getParticleData(decay[1]->id())))
    swap(outgoingVertex1_[inter], outgoingVertex2_[inter]);
  
  if (! (outgoingVertex1_[inter] && outgoingVertex2_[inter]))
    throw Exception()
      << "Invalid vertices for radiation in TVV decay in TVVDecayer::threeBodyME"
      << Exception::runerror;

  if( !(!inpart.dataPtr()->coloured() && inter ==ShowerInteraction::QCD) &&
      !(!inpart.dataPtr()->charged()  && inter ==ShowerInteraction::QED))
    throw Exception()
      << "Invalid vertices for radiation in TVV decay in TVVDecayer::threeBodyME"
      << Exception::runerror;


  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int it = 0; it < 5; ++it) {  
    for(unsigned int iv0 = 0; iv0 < 3; ++iv0) {
      for(unsigned int iv1 = 0; iv1 < 3; ++iv1) {
	for(unsigned int ig = 0; ig < 2; ++ig) {

	  // radiation from first outgoing vector
	  if((decay[0]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[0]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertex1_[inter]);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[0]->dataPtr();
	    if(off->CC()) off = off->CC();
	    VectorWaveFunction vectInter = 
	      outgoingVertex1_[inter]->evaluate(scale,3,off,gluon_[2*ig],
						 vectors3_[0][iv0],decay[0]->mass());
	  
	    assert(vectors3_[0][iv0].particle()->PDGName()==vectInter.particle()->PDGName());

	    Complex diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,vectors3_[1][iv1], 
				       vectInter,tensors3_[it]);
	    if(!couplingSet) {
	      gs = abs(outgoingVertex1_[inter]->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	      (*ME[colourFlow[1][ix].first])(it, iv0, iv1, ig) += 
		colourFlow[1][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }

	  // radiation from second outgoing vector
	  if((decay[1]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[1]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertex2_[inter]);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[1]->dataPtr();
	    if(off->CC()) off = off->CC();
	    VectorWaveFunction  vectInter = 
	      outgoingVertex2_[inter]->evaluate(scale,3,off,vectors3_[1][iv1],
						gluon_[2*ig],decay[1]->mass());
	    
	    assert(vectors3_[1][iv1].particle()->PDGName()==vectInter.particle()->PDGName());
	    
	    Complex diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,vectInter,vectors3_[0][iv0],
					tensors3_[it]);
	    if(!couplingSet) {
	      gs = abs(outgoingVertex2_[inter]->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	      (*ME[colourFlow[2][ix].first])(it, iv0, iv1, ig) += 
		colourFlow[2][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }

	  // radiation from 4 point vertex
	  if (fourPointVertex_[inter]) {
	    Complex diag = fourPointVertex_[inter]->evaluate(scale, vectors3_[0][iv0],
							     vectors3_[1][iv1],gluon_[2*ig], 
							     tensors3_[it]);
	    for(unsigned int ix=0;ix<colourFlow[3].size();++ix) {
	      (*ME[colourFlow[3][ix].first])(it, iv0, iv1, ig) += 
		colourFlow[3][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }
	}
	if(massless[1]) ++iv1;
      }
      if(massless[0]) ++iv0;
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
