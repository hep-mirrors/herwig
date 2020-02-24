// -*- C++ -*-
//
// SSVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
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

void SSVDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> fourV) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex)
    vertex_            .push_back(dynamic_ptr_cast<AbstractVSSVertexPtr>(vert));
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter]  = dynamic_ptr_cast<AbstractVSSVertexPtr>(inV.at(inter));
    fourPointVertex_[inter] = dynamic_ptr_cast<AbstractVVSSVertexPtr>(fourV.at(inter));
    outgoingVertexS_[inter] = AbstractVSSVertexPtr();
    outgoingVertexV_[inter] = AbstractVVVVertexPtr();  
    if(outV[0].at(inter)) {
      if (outV[0].at(inter)->getName()==VertexType::VSS)
	outgoingVertexS_[inter]   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[0].at(inter));
      else
	outgoingVertexV_[inter]   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[0].at(inter));
    }
    if(outV[1].at(inter)) {
      if (outV[1].at(inter)->getName()==VertexType::VSS)
	outgoingVertexS_[inter]   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[1].at(inter));
      else 
	outgoingVertexV_[inter]   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[1].at(inter));
    }
  }
}

void SSVDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_
     << incomingVertex_   << outgoingVertexS_
     << outgoingVertexV_  << fourPointVertex_;
}

void SSVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_
     >> incomingVertex_   >> outgoingVertexS_
     >> outgoingVertexV_  >> fourPointVertex_;
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
    // fix rho if no correlations
    fixRho(rho_);
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
    for(unsigned int ix = 0; ix < 3; ++ix) {
      (*ME())(0, ix, 0) = 0.;
      for(auto vert : vertex_)
	(*ME())(0, ix, 0) += vert->evaluate(scale,vector_[ix],sca, swave_);
    }
  }
  else {
    for(unsigned int ix = 0; ix < 3; ++ix) {
      (*ME())(0, 0, ix) = 0.;
      for(auto vert : vertex_)
	(*ME())(0, 0, ix) += vert->evaluate(scale,vector_[ix],sca,swave_);
    }
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
  return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
}


double  SSVDecayer::threeBodyME(const int , const Particle & inpart,
				const ParticleVector & decay,
				ShowerInteraction inter, MEOption meopt) {
  int iscal (0), ivect (1), iglu (2);
  // get location of outgoing scalar/vector
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin0) swap(iscal,ivect);

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
  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin0, PDT::Spin0,
								       PDT::Spin1, PDT::Spin1)));

  // create wavefunctions
  ScalarWaveFunction scal(decay[iscal]->momentum(),  decay[iscal]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(vector3_,decay[ivect],outgoing,false);
  VectorWaveFunction::calculateWaveFunctions(gluon_,  decay[iglu ],outgoing,true );

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

  if (! ((incomingVertex_[inter]  && (outgoingVertexS_[inter] || outgoingVertexV_[inter])) ||
	 (outgoingVertexS_[inter] &&  outgoingVertexV_[inter])))
    throw Exception()
      << "Invalid vertices for radiation in SSV decay in SSVDecayer::threeBodyME"
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
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int iv = 0; iv < 3; ++iv) {
    for(unsigned int ig = 0; ig < 2; ++ig) {
      // radiation from the incoming scalar
      if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	 (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	assert(incomingVertex_[inter]);
	ScalarWaveFunction scalarInter = 
	  incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),
					   gluon_[2*ig],swave3_,inpart.mass());
	
	assert(swave3_.particle()->id()==scalarInter.particle()->id());
	Complex diag = 0.;
	for(auto vertex : vertex_)
	  diag += vertex->evaluate(scale,vector3_[iv],scal,scalarInter);
	if(!couplingSet) {
	  gs = abs(incomingVertex_[inter]->norm());
	  couplingSet = true;
	}
	for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	  (*ME[colourFlow[0][ix].first])(0, 0, iv, ig) += 
	    colourFlow[0][ix].second*diag; 
	}
#ifdef GAUGE_CHECK
	total+=norm(diag);
#endif
      }
      // radiation from the outgoing scalar
      if((decay[iscal]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	 (decay[iscal]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	assert(outgoingVertexS_[inter]);
	// ensure you get correct outgoing particle from first vertex
	tcPDPtr off = decay[iscal]->dataPtr();
	if(off->CC()) off = off->CC();
	ScalarWaveFunction scalarInter = 
	  outgoingVertexS_[inter]->evaluate(scale,3,off,gluon_[2*ig],scal,decay[iscal]->mass());
	
	assert(scal.particle()->id()==scalarInter.particle()->id());

	if(!couplingSet) {
	  gs = abs(outgoingVertexS_[inter]->norm());
	  couplingSet = true;
	}
	Complex diag = 0.;
	for(auto vertex : vertex_)
	  diag += vertex->evaluate(scale,vector3_[iv],scalarInter,swave3_);
	for(unsigned int ix=0;ix<colourFlow[S].size();++ix) {
	  (*ME[colourFlow[S][ix].first])(0, 0, iv, ig) += 
	    colourFlow[S][ix].second*diag;
	}
#ifdef GAUGE_CHECK
	total+=norm(diag);
#endif
      }

      // radiation from outgoing vector
      if((decay[ivect]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	 (decay[ivect]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	assert(outgoingVertexV_[inter]);
	// ensure you get correct outgoing particle from first vertex
	tcPDPtr off = decay[ivect]->dataPtr();
	if(off->CC()) off = off->CC();
	VectorWaveFunction  vectorInter = 
	  outgoingVertexV_[inter]->evaluate(scale,3,off,gluon_[2*ig],
					    vector3_[iv],decay[ivect]->mass());
	
	assert(vector3_[iv].particle()->id()==vectorInter.particle()->id());
	
	if(!couplingSet) {
	  gs = abs(outgoingVertexV_[inter]->norm());
	  couplingSet = true;
	}	
	Complex diag = 0.;
	for(auto vertex : vertex_)
	  diag += vertex->evaluate(scale,vectorInter,scal,swave3_);
	for(unsigned int ix=0;ix<colourFlow[V].size();++ix) {
	  (*ME[colourFlow[V][ix].first])(0, 0, iv, ig) += 
	    colourFlow[V][ix].second*diag;
	}
#ifdef GAUGE_CHECK
	total+=norm(diag);
#endif
      }
      // radiation from 4 point vertex
      if (fourPointVertex_[inter]) {
	Complex diag =  fourPointVertex_[inter]->evaluate(scale, gluon_[2*ig], vector3_[iv],
							  scal, swave3_);
	for(unsigned int ix=0;ix<colourFlow[3].size();++ix) {
	  (*ME[colourFlow[3][ix].first])(0, 0, iv, ig) += 
	     colourFlow[3][ix].second*diag;
	}
#ifdef GAUGE_CHECK
	total+=norm(diag);
#endif
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
