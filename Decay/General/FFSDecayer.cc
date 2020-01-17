// -*- C++ -*-
//
// FFSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
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

void FFSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
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
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(inV.at(inter));
    outgoingVertexF_[inter] = AbstractFFVVertexPtr();
    outgoingVertexS_[inter] = AbstractVSSVertexPtr();
    if(outV[0].at(inter)) {
      if (outV[0].at(inter)->getName()==VertexType::FFV)
	outgoingVertexF_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[0].at(inter));
      else
	outgoingVertexS_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[0].at(inter));
    }
    if(outV[1].at(inter)) {
      if (outV[1].at(inter)->getName()==VertexType::FFV)
	outgoingVertexF_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[1].at(inter));
      else
	outgoingVertexS_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[1].at(inter));
    }
  }
}

void FFSDecayer::persistentOutput(PersistentOStream & os) const {
  os << perturbativeVertex_       << vertex_
     << incomingVertex_   << outgoingVertexF_
     << outgoingVertexS_;
}

void FFSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> perturbativeVertex_       >> vertex_
     >> incomingVertex_   >> outgoingVertexF_
     >> outgoingVertexS_;
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
      if(ferm) (*ME())(if1, if2, 0) = 0.;
      else     (*ME())(if2, if1, 0) = 0.;
      for(auto vert : vertex_) {
	if(ferm) (*ME())(if1, if2, 0) += 
		   vert->evaluate(scale,wave_[if1],wavebar_[if2],scal);
	else     (*ME())(if2, if1, 0) += 
		   vert->evaluate(scale,wave_[if1],wavebar_[if2],scal);
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

Energy FFSDecayer::partialWidth(PMPair inpart, PMPair outa,
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    double mu1(0.),mu2(0.);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if(outa.first->iSpin() == PDT::Spin1Half) {
      mu1 = outa.second/inpart.second;
      mu2 = outb.second/inpart.second;
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in, outa.first, outb.first);
    }
    else {
      mu1 = outb.second/inpart.second;
      mu2 = outa.second/inpart.second;
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in, outb.first, outa.first);
      
    }
    double c2 = norm(perturbativeVertex_[0]->norm());
    Complex cl = perturbativeVertex_[0]->left();
    Complex cr = perturbativeVertex_[0]->right();
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
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  int iscal (0), iferm (1), iglu (2);
  // get location of outgoing fermion/scalar
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin0) swap(iscal,iferm);
  // work out whether inpart is a fermion or antifermion
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[iferm]->dataPtr()->CC()) itype[1] = decay[iferm]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;

  bool ferm(false);
  if(itype[0] == itype[1] ) {
    ferm = itype[0]==0 || (itype[0]==2 && decay[iscal]->id() < 0);
  }
  else if(itype[0] == 2) {
    ferm = itype[1]==0;
  }
  else if(itype[1] == 2) {
    ferm = itype[0]==0;
  }
  else if((itype[0] == 1 && itype[1] == 0) ||
	  (itype[0] == 0 && itype[1] == 1)) {
    if(abs(inpart.id())<=16) {
      ferm = itype[0]==0;
    }
    else if(abs(decay[iferm]->id())<=16) {
      ferm = itype[1]==0;
    }
    else {
      ferm = true;
    }
  }
  else
    assert(false);
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

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half, PDT::Spin0,
								       PDT::Spin1Half, PDT::Spin1)));
  // create wavefunctions
  if (ferm)  SpinorBarWaveFunction::calculateWaveFunctions(wavebar3_, decay[iferm],outgoing);
  else       SpinorWaveFunction::   calculateWaveFunctions(wave3_   , decay[iferm],outgoing);
  
  ScalarWaveFunction swave3_(decay[iscal]->momentum(), decay[iscal]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(gluon_,   decay[iglu ],outgoing,true);

  // gauge invariance test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),decay[iglu ]->dataPtr(),10,
					  outgoing));
    }
  }
#endif
  
  if (! ((incomingVertex_[inter]  && (outgoingVertexF_[inter] || outgoingVertexS_[inter])) ||
	 (outgoingVertexF_[inter] &&  outgoingVertexS_[inter])))
    throw Exception()
      << "Invalid vertices for radiation in FFS decay in FFSDecayer::threeBodyME"
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
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int ifi = 0; ifi < 2; ++ifi) {
    for(unsigned int ifo = 0; ifo < 2; ++ifo) {
      for(unsigned int ig = 0; ig < 2; ++ig) {
   	// radiation from the incoming fermion
   	if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	   (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
   	  assert(incomingVertex_[inter]);
	  if (ferm) {
	    SpinorWaveFunction spinorInter =
	      incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),wave3_[ifi],
					       gluon_[2*ig],inpart.mass());

	    assert(wave3_[ifi].particle()->id()==spinorInter.particle()->id());
	    diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,spinorInter,wavebar3_[ifo],swave3_);
	  }
	  else {
	    SpinorBarWaveFunction spinorBarInter = 
	      incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),wavebar3_[ifi],
					       gluon_[2*ig],inpart.mass());

	    assert(wavebar3_[ifi].particle()->id()==spinorBarInter.particle()->id());
	    diag = 0.;
	    for(auto vertex :vertex_)
	      diag+= vertex->evaluate(scale,wave3_[ifo], spinorBarInter,swave3_);
	  }
	  if(!couplingSet) {
	    gs = abs(incomingVertex_[inter]->norm());
	    couplingSet = true;
	  }
	  for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	    (*ME[colourFlow[0][ix].first])(ifi, 0, ifo, ig) += 
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
	    for(auto vertex :vertex_)
	      diag+= vertex->evaluate(scale,wave3_[ifi],spinorBarInter,swave3_);
	  }
	  else {
	    SpinorWaveFunction spinorInter = 
	      outgoingVertexF_[inter]->evaluate(scale,3,off,wave3_[ifo],
						gluon_[2*ig],decay[iferm]->mass());
	      
	    assert(wave3_[ifo].particle()->id()==spinorInter.particle()->id());
	    diag = 0.;
	    for(auto vertex :vertex_)
	      diag+= vertex->evaluate(scale,spinorInter,wavebar3_[ifi],swave3_);
	  }
	  if(!couplingSet) {
	    gs = abs(outgoingVertexF_[inter]->norm());
	    couplingSet = true;
	  }
	  for(unsigned int ix=0;ix<colourFlow[F].size();++ix) {
	    (*ME[colourFlow[F][ix].first])(ifi, 0, ifo, ig) += 
	      colourFlow[F][ix].second*diag;
	  }
#ifdef GAUGE_CHECK
	  total+=norm(diag);
#endif
  	}

  	// radiation from outgoing scalar
   	if((decay[iscal]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	   (decay[iscal]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
  	  assert(outgoingVertexS_[inter]);
	  // ensure you get correct ougoing particle from first vertex
	  tcPDPtr off = decay[iscal]->dataPtr();
	  if(off->CC()) off = off->CC();
	  ScalarWaveFunction  scalarInter = 
	    outgoingVertexS_[inter]->evaluate(scale,3,off,gluon_[2*ig],
					      swave3_,decay[iscal]->mass());
	    
	  assert(swave3_.particle()->id()==scalarInter.particle()->id());
	  if (ferm){
	    diag = 0.;
	    for(auto vertex :vertex_)
	      diag += vertex->evaluate(scale,wave3_[ifi],wavebar3_[ifo],scalarInter);
	  }
	  else {
	    diag = 0.;
	    for(auto vertex :vertex_)
	      diag += vertex->evaluate(scale,wave3_[ifo],wavebar3_[ifi],scalarInter);
	  }
	  if(!couplingSet) {
	    gs = abs(outgoingVertexS_[inter]->norm());
	    couplingSet = true;
	  }
	  for(unsigned int ix=0;ix<colourFlow[S].size();++ix) {
  	    (*ME[colourFlow[S][ix].first])(ifi, 0, ifo, ig) += 
	      colourFlow[S][ix].second*diag;
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
