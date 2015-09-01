// -*- C++ -*-
//
// VSSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VSSDecayer class.
//

#include "VSSDecayer.h"
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

IBPtr VSSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VSSDecayer::fullclone() const {
  return new_ptr(*this);
}

void VSSDecayer::doinit() {
  _perturbativeVertex      = dynamic_ptr_cast<VSSVertexPtr>        (getVertex());
  _abstractVertex          = dynamic_ptr_cast<AbstractVSSVertexPtr>(getVertex());
  _abstractIncomingVertex  = dynamic_ptr_cast<AbstractVVVVertexPtr>(getIncomingVertex());
  _abstractOutgoingVertex1 = dynamic_ptr_cast<AbstractVSSVertexPtr>(getOutgoingVertices()[0]);
  _abstractOutgoingVertex2 = dynamic_ptr_cast<AbstractVSSVertexPtr>(getOutgoingVertices()[1]);
  GeneralTwoBodyDecayer::doinit();
}

void VSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex           << _perturbativeVertex
     << _abstractIncomingVertex   << _abstractOutgoingVertex1
     << _abstractOutgoingVertex2;
}

void VSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex           >> _perturbativeVertex
     >> _abstractIncomingVertex   >> _abstractOutgoingVertex1
     >> _abstractOutgoingVertex2;
}

ClassDescription<VSSDecayer> VSSDecayer::initVSSDecayer;
// Definition of the static class description member.

void VSSDecayer::Init() {

  static ClassDocumentation<VSSDecayer> documentation
    ("This implements the decay of a vector to 2 scalars");

}

double VSSDecayer::me2(const int , const Particle & inpart,
 		       const ParticleVector & decay, 
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    for(unsigned int ix=0;ix<2;++ix)
      ScalarWaveFunction::
	constructSpinInfo(decay[ix],outgoing,true);
    return 0.;
  }
  ScalarWaveFunction sca1(decay[0]->momentum(),decay[0]->dataPtr(),outgoing);
  ScalarWaveFunction sca2(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int ix=0;ix<3;++ix) {
    (*ME())(ix,0,0) = _abstractVertex->evaluate(scale,_vectors[ix],sca1,sca2);
  }
  double output=(ME()->contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy VSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    _perturbativeVertex->setCoupling(sqr(inpart.second), in, outa.first,
				     outb.first);
    double mu1sq = sqr(outa.second/inpart.second);
    double mu2sq = sqr(outb.second/inpart.second);
    double me2 = sqr(mu1sq - mu2sq) - 2.*(mu1sq + mu2sq);
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = -norm(_perturbativeVertex->norm())*me2*pcm /
      (24.*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double VSSDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay, MEOption meopt) {

  bool massless = inpart.mass()==ZERO;
 
  // work out which is the scalar and anti-scalar
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
    // create vector wavefunction for decaying particle
    VectorWaveFunction::calculateWaveFunctions(_vector3, _rho3, const_ptr_cast<tPPtr>(&inpart), 
					       incoming, massless);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    VectorWaveFunction::
      constructSpinInfo(_vector3 ,const_ptr_cast<tPPtr>(&inpart),outgoing,true,massless);
    ScalarWaveFunction::constructSpinInfo(       decay[iscal],outgoing,true);
    ScalarWaveFunction::constructSpinInfo(       decay[ianti],outgoing,true);
    VectorWaveFunction::constructSpinInfo(_gluon,decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  if(nflow==2) cfactors[0][1]=cfactors[1][0];
  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1, PDT::Spin0,
								       PDT::Spin0, PDT::Spin1)));

  // create wavefunctions
  ScalarWaveFunction scal(decay[iscal]->momentum(), decay[iscal]->dataPtr(),outgoing);
  ScalarWaveFunction anti(decay[ianti]->momentum(), decay[ianti]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(_gluon,decay[iglu ],outgoing,true);

  // gauge test
  // _gluon.clear();
  // for(unsigned int ix=0;ix<3;++ix) {
  //   if(ix==1) _gluon.push_back(VectorWaveFunction());
  //   else {
  //     _gluon.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  // 					  decay[iglu ]->dataPtr(),10,
  // 					  outgoing));
  //   }
  // }

  // identify scalar and/or anti-scalar vertex
  AbstractVSSVertexPtr abstractOutgoingVertexS;
  AbstractVSSVertexPtr abstractOutgoingVertexA;
  identifyVertices(iscal, ianti, inpart, decay, abstractOutgoingVertexS, abstractOutgoingVertexA);

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);

  for(unsigned int iv = 0; iv < 3; ++iv) {
    for(unsigned int ig = 0; ig < 2; ++ig) {
      // radiation from the incoming vector
      if(inpart.dataPtr()->coloured()) {
	assert(_abstractIncomingVertex);	
	VectorWaveFunction vectorInter = 
	  _abstractIncomingVertex->evaluate(scale,3,inpart.dataPtr(),_vector3[iv],
					    _gluon[2*ig],inpart.mass());
	
	if (_vector3[iv].particle()->PDGName()!=vectorInter.particle()->PDGName())
	  throw Exception()
	    << _vector3[iv].particle()->PDGName() << " was changed to " 
	    << vectorInter .particle()->PDGName() << " in VSSDecayer::threeBodyME"
	    << Exception::runerror;
	
	double gs    = _abstractIncomingVertex->strongCoupling(scale);
	Complex diag = _abstractVertex->evaluate(scale,vectorInter,scal,anti)/gs;
	for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	  (*ME[colourFlow[0][ix].first])(iv, 0, 0, ig) += 
	    colourFlow[0][ix].second*diag;
	}
      }
      // radiation from the outgoing scalar
      if(decay[iscal]->dataPtr()->coloured()) {
	assert(abstractOutgoingVertexS);
	// ensure you get correct outgoing particle from first vertex
	tcPDPtr off = decay[iscal]->dataPtr();
	if(off->CC()) off = off->CC();
	ScalarWaveFunction scalarInter = 
	  abstractOutgoingVertexS->evaluate(scale,3,off,_gluon[2*ig],scal,decay[iscal]->mass());
	
	if (scal.particle()->PDGName()!=scalarInter.particle()->PDGName())
	  throw Exception()
	    << scal       .particle()->PDGName() << " was changed to " 
	    << scalarInter.particle()->PDGName() << " in VSSDecayer::threeBodyME"
	    << Exception::runerror;
	
	double gs    = abstractOutgoingVertexS->strongCoupling(scale);
	Complex diag =_abstractVertex->evaluate(scale,_vector3[iv],anti,scalarInter)/gs;
	for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	  (*ME[colourFlow[1][ix].first])(iv, 0, 0, ig) += 
	    colourFlow[1][ix].second*diag;
	}
      }
      
      if(decay[ianti]->dataPtr()->coloured()) {
	assert(abstractOutgoingVertexA);
	// ensure you get correct outgoing particle from first vertex
	tcPDPtr off = decay[ianti]->dataPtr();
	if(off->CC()) off = off->CC();
	ScalarWaveFunction scalarInter = 
	  abstractOutgoingVertexA->evaluate(scale,3,off, _gluon[2*ig],anti,decay[ianti]->mass());
	
	if (anti.particle()->PDGName()!=scalarInter.particle()->PDGName())
	  throw Exception()
	    << anti       .particle()->PDGName() << " was changed to " 
	    << scalarInter.particle()->PDGName() << " in VSSDecayer::threeBodyME"
	    << Exception::runerror;
	
	double gs    = abstractOutgoingVertexA->strongCoupling(scale);
	Complex diag =_abstractVertex->evaluate(scale,_vector3[iv],scal,scalarInter)/gs;
	
	for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	  (*ME[colourFlow[2][ix].first])(iv, 0, 0, ig) += 
	    colourFlow[2][ix].second*diag;
	}
      }
    }
    if(massless) ++iv;
  }

  // contract matrices 
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],_rho3)).real();
    }
  }
  output*=(4.*Constants::pi);

  // return the answer
  return output;
}


void VSSDecayer::identifyVertices(const int iscal, const int ianti,
				  const Particle & inpart, const ParticleVector & decay, 
				  AbstractVSSVertexPtr & abstractOutgoingVertexS, 
				  AbstractVSSVertexPtr & abstractOutgoingVertexA){

  // work out which scalar each outgoing vertex corresponds to 
  // two outgoing vertices
  if( inpart.dataPtr()       ->iColour()==PDT::Colour0     &&
    ((decay[iscal]->dataPtr()->iColour()==PDT::Colour3     &&
      decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) ||
     (decay[iscal]->dataPtr()->iColour()==PDT::Colour8     &&
      decay[ianti]->dataPtr()->iColour()==PDT::Colour8))){
    if(_abstractOutgoingVertex1==_abstractOutgoingVertex2){
      abstractOutgoingVertexS = _abstractOutgoingVertex1;
      abstractOutgoingVertexA = _abstractOutgoingVertex2;
    }
    else if (_abstractOutgoingVertex1->isIncoming(getParticleData(decay[iscal]->id()))){
      abstractOutgoingVertexS = _abstractOutgoingVertex1;
      abstractOutgoingVertexA = _abstractOutgoingVertex2;
    }
    else if (_abstractOutgoingVertex2->isIncoming(getParticleData(decay[iscal]->id()))){
      abstractOutgoingVertexS = _abstractOutgoingVertex2;
      abstractOutgoingVertexA = _abstractOutgoingVertex1;
    }
  }
  else if(inpart.dataPtr()       ->iColour()==PDT::Colour8 &&
	  decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
    if(_abstractOutgoingVertex1==_abstractOutgoingVertex2){
      abstractOutgoingVertexS = _abstractOutgoingVertex1;
      abstractOutgoingVertexA = _abstractOutgoingVertex2;
    }
    else if (_abstractOutgoingVertex1->isIncoming(getParticleData(decay[iscal]->id()))){
      abstractOutgoingVertexS = _abstractOutgoingVertex1;
      abstractOutgoingVertexA = _abstractOutgoingVertex2;
    }
    else if (_abstractOutgoingVertex2->isIncoming(getParticleData(decay[iscal]->id()))){
      abstractOutgoingVertexS = _abstractOutgoingVertex2;
      abstractOutgoingVertexA = _abstractOutgoingVertex1;
    }
  }

  // one outgoing vertex
  else if(inpart.dataPtr()->iColour()==PDT::Colour3){
    if(decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&  
       decay[ianti]->dataPtr()->iColour()==PDT::Colour0){
      if     (_abstractOutgoingVertex1) abstractOutgoingVertexS = _abstractOutgoingVertex1;
      else if(_abstractOutgoingVertex2) abstractOutgoingVertexS = _abstractOutgoingVertex2;
    }
    else if (decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&
	     decay[ianti]->dataPtr()->iColour()==PDT::Colour8){
      if (_abstractOutgoingVertex1->isIncoming(getParticleData(decay[ianti]->dataPtr()->id()))){
	abstractOutgoingVertexS = _abstractOutgoingVertex2;
	abstractOutgoingVertexA = _abstractOutgoingVertex1;
      }
      else {
	abstractOutgoingVertexS = _abstractOutgoingVertex1;
	abstractOutgoingVertexA = _abstractOutgoingVertex2;
      }
    }
  }
  else if(inpart.dataPtr()->iColour()==PDT::Colour3bar){
    if(decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
       decay[iscal]->dataPtr()->iColour()==PDT::Colour0){
      if     (_abstractOutgoingVertex1) abstractOutgoingVertexA = _abstractOutgoingVertex1;
      else if(_abstractOutgoingVertex2) abstractOutgoingVertexA = _abstractOutgoingVertex2;
    }
    else if (decay[iscal]->dataPtr()->iColour()==PDT::Colour8 &&
	     decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
      if (_abstractOutgoingVertex1->isIncoming(getParticleData(decay[iscal]->dataPtr()->id()))){
	abstractOutgoingVertexS = _abstractOutgoingVertex1;
	abstractOutgoingVertexA = _abstractOutgoingVertex2;
      }
      else {
	abstractOutgoingVertexS = _abstractOutgoingVertex2;
	abstractOutgoingVertexA = _abstractOutgoingVertex1;
      }
    }
  }

  if (! ((_abstractIncomingVertex  && (abstractOutgoingVertexS  || abstractOutgoingVertexA)) ||
	 ( abstractOutgoingVertexS &&  abstractOutgoingVertexA)))
    throw Exception()
      << "Invalid vertices for QCD radiation in VSS decay in VSSDecayer::identifyVertices"
      << Exception::runerror;

  // // prohibit all for now since all unchecked
  // if (true)
  //   throw Exception()
  //     << "Invalid vertices for QCD radiation in VSS decay in VSSDecayer::identifyVertices"
  //     << Exception::runerror;

}

