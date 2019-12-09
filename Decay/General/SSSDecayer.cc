// -*- C++ -*-
//
// SSSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSSDecayer class.
//

#include "SSSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr SSSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SSSDecayer::fullclone() const {
  return new_ptr(*this);
}

void SSSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractSSSVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<SSSVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(inV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[1].at(inter));
  }
}

void SSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_           << perturbativeVertex_
     << incomingVertex_   << outgoingVertex1_
     << outgoingVertex2_;
}

void SSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_           >> perturbativeVertex_
     >> incomingVertex_   >> outgoingVertex1_
     >> outgoingVertex2_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSSDecayer,GeneralTwoBodyDecayer>
describeHerwigSSSDecayer("Herwig::SSSDecayer", "Herwig.so");

void SSSDecayer::Init() {

  static ClassDocumentation<SSSDecayer> documentation
    ("This class implements the decay of a scalar to 2 scalars.");

}

double SSSDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0)));
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
      ScalarWaveFunction::
	constructSpinInfo(decay[ix],outgoing,true);
  }
  ScalarWaveFunction s1(decay[0]->momentum(),decay[0]->dataPtr(),outgoing);
  ScalarWaveFunction s2(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  Energy2 scale(sqr(inpart.mass()));
  (*ME())(0,0,0) = 0.;
  for(auto vert : vertex_) {
    (*ME())(0,0,0) += vert->evaluate(scale,s1,s2,swave_);
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy SSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0] && !perturbativeVertex_[0]->kinematics()) {
    Energy2 scale(sqr(inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(scale, in, outa.first, outb.first);
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					       outb.second);
    double c2 = norm(perturbativeVertex_[0]->norm());
    Energy pWidth = c2*pcm/8./Constants::pi/scale*UnitRemoval::E2;
    // colour factor
    pWidth *= colourFactor(inpart.first,outa.first,outb.first);
    return pWidth;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double SSSDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  // work out which is the scalar and anti scalar
  int ianti(0), iscal(1), iglu(2);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0 && itype[1]!=0) swap(ianti, iscal);
  if(itype[0]==2 && itype[1]==1) swap(ianti, iscal);
  if(itype[0]==0 && itype[1]==0 && abs(decay[0]->dataPtr()->id())>abs(decay[1]->dataPtr()->id())) 
    swap(iscal, ianti);
  if(itype[0]==1 && itype[1]==1 && abs(decay[0]->dataPtr()->id())<abs(decay[1]->dataPtr()->id())) 
    swap(iscal, ianti);

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
    ScalarWaveFunction::
      constructSpinInfo(decay[ianti],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(gluon_,decay[iglu ],outgoing,true,false);
    return 0.;
  }
  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin0, PDT::Spin0,
								       PDT::Spin0, PDT::Spin1)));

  // create wavefunctions
  ScalarWaveFunction scal(decay[iscal]->momentum(), decay[iscal]->dataPtr(),outgoing);
  ScalarWaveFunction anti(decay[ianti]->momentum(), decay[ianti]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(gluon_,decay[iglu ],outgoing,true);

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

  AbstractVSSVertexPtr outgoingVertexS;
  AbstractVSSVertexPtr outgoingVertexA;
  identifyVertices(iscal, ianti, inpart, decay, outgoingVertexS, outgoingVertexA,inter);

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
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
	diag += vertex->evaluate(scale,scal,anti,scalarInter);
      if(!couplingSet) {
	gs = abs(incomingVertex_[inter]->norm());
	couplingSet = true;
      }
      for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	(*ME[colourFlow[0][ix].first])(0, 0, 0, ig) += 
	   colourFlow[0][ix].second*diag; 
      }
#ifdef GAUGE_CHECK
      total+=norm(diag);
#endif
    }
    // radiation from the outgoing scalar
    if((decay[iscal]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
       (decay[iscal]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
      assert(outgoingVertexS);
      // ensure you get correct outgoing particle from first vertex
      tcPDPtr off = decay[iscal]->dataPtr();
      if(off->CC()) off = off->CC();
      ScalarWaveFunction scalarInter = 
	outgoingVertexS->evaluate(scale,3,off,gluon_[2*ig],scal,decay[iscal]->mass());

      assert(scal.particle()->id()==scalarInter.particle()->id());

      Complex diag = 0.;
      for(auto vertex : vertex_)
	diag += vertex->evaluate(scale,swave3_,anti,scalarInter);
      if(!couplingSet) {
	gs = abs(outgoingVertexS->norm());
	couplingSet = true;
      }
      for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	(*ME[colourFlow[1][ix].first])(0, 0, 0, ig) += 
	   colourFlow[1][ix].second*diag;
      }
#ifdef GAUGE_CHECK
      total+=norm(diag);
#endif
    }
    // radiation from the outgoing anti scalar
    if((decay[ianti]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
       (decay[ianti]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
      assert(outgoingVertexA);
      // ensure you get correct outgoing particle from first vertex
      tcPDPtr off = decay[ianti]->dataPtr();
      if(off->CC()) off = off->CC();
      ScalarWaveFunction scalarInter = 
	outgoingVertexA->evaluate(scale,3,off, gluon_[2*ig],anti,decay[ianti]->mass());

      assert(anti.particle()->id()==scalarInter.particle()->id());

      Complex diag = 0.;
      for(auto vertex : vertex_)
	diag += vertex->evaluate(scale,swave3_,scal,scalarInter);
      if(!couplingSet) {
	gs = abs(outgoingVertexA->norm());
	couplingSet = true;
      }
      for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	(*ME[colourFlow[2][ix].first])(0, 0, 0, ig) += 
	  colourFlow[2][ix].second*diag;
      }
#ifdef GAUGE_CHECK
      total+=norm(diag);
#endif
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


void SSSDecayer::identifyVertices(const int iscal, const int ianti,
				  const Particle & inpart, const ParticleVector & decay, 
				  AbstractVSSVertexPtr & outgoingVertexS, 
				  AbstractVSSVertexPtr & outgoingVertexA,
				  ShowerInteraction inter){
  // QCD
  if(inter==ShowerInteraction::QCD) {
    // work out which scalar each outgoing vertex corresponds to 
    // two outgoing vertices
    if( inpart.dataPtr()       ->iColour()==PDT::Colour0     &&
	((decay[iscal]->dataPtr()->iColour()==PDT::Colour3     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) ||
	 (decay[iscal]->dataPtr()->iColour()==PDT::Colour8     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour8))){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    else if(inpart.dataPtr()       ->iColour()==PDT::Colour8 &&
	    decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&
	    decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    // one outgoing vertex
    else if(inpart.dataPtr()    ->iColour()==PDT::Colour3){ 
      if(decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&  
	 decay[ianti]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertexS = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertexS = outgoingVertex2_[inter]; 
      }
      else if (decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&  
	       decay[ianti]->dataPtr()->iColour()==PDT::Colour8){
	if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[ianti]->dataPtr()->id()))){
	  outgoingVertexS = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
	else {
	  outgoingVertexS = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
      }
    }
    else if(inpart.dataPtr()    ->iColour()==PDT::Colour3bar){
      if(decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
	 decay[iscal]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertexA = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
	       decay[iscal]->dataPtr()->iColour()==PDT::Colour8){
	if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[iscal]->dataPtr()->id()))){
	  outgoingVertexS = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertexS = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
      }
    }
    
    if (! ((incomingVertex_[inter]  && (outgoingVertexS || outgoingVertexA)) ||
	   ( outgoingVertexS &&  outgoingVertexA)))
      throw Exception()
	<< "Invalid vertices for QCD radiation in SSS decay in SSSDecayer::identifyVertices"
	<< Exception::runerror;
  }
  // QED
  else {
    if(decay[iscal]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iscal]->dataPtr())))
	outgoingVertexS = outgoingVertex1_[inter];
      else
	outgoingVertexS = outgoingVertex2_[inter];
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
