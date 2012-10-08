// -*- C++ -*-
//
// FFSDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFSDecayer class.
//

#include "FFSDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr FFSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FFSDecayer::fullclone() const {
  return new_ptr(*this);
}

void FFSDecayer::doinit() {
  _perturbativeVertex = dynamic_ptr_cast<FFSVertexPtr>        (getVertex());
  _abstractVertex     = dynamic_ptr_cast<AbstractFFSVertexPtr>(getVertex());
  _abstractIncomingVertex    = dynamic_ptr_cast<AbstractFFVVertexPtr>(getIncomingVertex());
  _abstractOutgoingVertex1   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[0]);
  _abstractOutgoingVertex2   = dynamic_ptr_cast<AbstractVSSVertexPtr>(getOutgoingVertices()[1]);
  GeneralTwoBodyDecayer::doinit();
}

void FFSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _perturbativeVertex       << _abstractVertex
     << _abstractIncomingVertex   << _abstractOutgoingVertex1
     << _abstractOutgoingVertex2;
}

void FFSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _perturbativeVertex       >> _abstractVertex
     >> _abstractIncomingVertex   >> _abstractOutgoingVertex1
     >> _abstractOutgoingVertex2;
}

ClassDescription<FFSDecayer> FFSDecayer::initFFSDecayer;
// Definition of the static class description member.

void FFSDecayer::Init() {

  static ClassDocumentation<FFSDecayer> documentation
    ("The FFSDecayer class implements the decay of a fermion to "
     "a fermion and a scalar.");

}

double FFSDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[0]->dataPtr()->CC()) itype[1] = decay[0]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));
  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(_wave,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(_wave[0].wave().Type() != u_spinortype)
	for(unsigned int ix = 0; ix < 2; ++ix) _wave   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(_wavebar,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(_wavebar[0].wave().Type() != v_spinortype)
	for(unsigned int ix = 0; ix < 2; ++ix) _wavebar[ix].conjugate();
    }
    ME(DecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0));
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(ferm) {
      SpinorWaveFunction::
	constructSpinInfo(_wave,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(_wavebar,decay[0],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_wavebar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(_wave,decay[0],outgoing,true);
    }
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
  }
  if(ferm)
    SpinorBarWaveFunction::
      calculateWaveFunctions(_wavebar,decay[0],outgoing);
  else
    SpinorWaveFunction::
      calculateWaveFunctions(_wave   ,decay[0],outgoing);
  ScalarWaveFunction scal(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      if(ferm) ME()(if1, if2, 0) = 
	_abstractVertex->evaluate(scale,_wave[if1],_wavebar[if2],scal);
      else     ME()(if2, if1, 0) = 
	_abstractVertex->evaluate(scale,_wave[if1],_wavebar[if2],scal);
    }
  }
  double output = (ME().contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy FFSDecayer::partialWidth(PMPair inpart, PMPair outa,
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    double mu1(0.),mu2(0.);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if(outa.first->iSpin() == PDT::Spin1Half) {
      mu1 = outa.second/inpart.second;
      mu2 = outb.second/inpart.second;
      _perturbativeVertex->setCoupling(sqr(inpart.second), in, outa.first, outb.first);
    }
    else {
      mu1 = outb.second/inpart.second;
      mu2 = outa.second/inpart.second;
      _perturbativeVertex->setCoupling(sqr(inpart.second), in, outb.first, outa.first);
      
    }
    double c2 = norm(_perturbativeVertex->norm());
    Complex cl = _perturbativeVertex->left();
    Complex cr = _perturbativeVertex->right();
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
			       const ParticleVector & decay, MEOption meopt) {

  ofstream output("output.top", ios::app);
  output << inpart.PDGName() << " -> " << decay[0]->PDGName() 
	 << " \t" << decay[1]->PDGName() << "\n";

  //fermion = 1, antifermion =-1, scalar = 0, majorana fermion = 2
  int iscal (0), ilept (1), iglu (2);
  // work out whether inpart is a fermion or antifermion
  int in;
  if(inpart.dataPtr()->CC()) in = inpart.id()>0 ? 1:-1;
  else                        in = 2;
  int out[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if     (decay[ix]->dataPtr()->iSpin()==PDT::Spin0) out[ix] = 0;
    else if(decay[ix]->dataPtr()->CC())                out[ix] = decay[ix]->id()>0 ? 1:-1;
    else                                               out[ix] = 2;
  }

  if(out[1]==0) swap(iscal, ilept);
  if(in==2){
    if     (out[ilept]== 1)       in =  1;
    else if(out[ilept]==-1)       in = -1;
    else if(decay[iscal]->id()>0) in = -1;
    else if(decay[iscal]->id()<0) in =  1;
  }
  
  output << "inparticle " << in << endl;
  output << "out lepton " << decay[ilept]->PDGName() << "\t" << out[ilept] << endl;
  output << "out scalar " << decay[iscal]->PDGName() << "\t" << out[iscal] << endl;

  // if(meopt==Initialize) {
  //   ScalarWaveFunction::
  //     calculateWaveFunctions(_rho3,const_ptr_cast<tPPtr>(&inpart),incoming);
  //   _swave3 = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  // }
  // if(meopt==Terminate) {
  //   ScalarWaveFunction::
  //     constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
  //   SpinorBarWaveFunction::
  //     constructSpinInfo(_wavebar3 ,decay[iferm],outgoing,true);
  //   SpinorWaveFunction::
  //     constructSpinInfo(_wave3    ,decay[ianti],outgoing,true);
  //   VectorWaveFunction::
  //     constructSpinInfo(_vwave3   ,decay[iglu ],outgoing,true,false);
  //   return 0.;
  // }

  // unsigned int nflow;
  // vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  // vector<DecayMatrixElement> ME(nflow,DecayMatrixElement(PDT::Spin0,     PDT::Spin1Half,
  // 							 PDT::Spin1Half, PDT::Spin1));

  // SpinorBarWaveFunction::
  //   calculateWaveFunctions(_wavebar3, decay[iferm],outgoing);
  // SpinorWaveFunction::
  //   calculateWaveFunctions(_wave3   , decay[ianti],outgoing);
  // VectorWaveFunction::
  //   calculateWaveFunctions(_vwave3  , decay[iglu ],outgoing,true);

  // // gauge test
  // //_vwave3.clear();
  // //for(unsigned int ix=0;ix<3;++ix) {
  // //if(ix==1) _vwave3.push_back(VectorWaveFunction());
  // //else {
  // //  _vwave3.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  // //					   decay[iglu ]->dataPtr(),10,
  // //					   outgoing));
  // //}
  // //}

  // AbstractFFVVertexPtr abstractOutgoingVertexF;
  // AbstractFFVVertexPtr abstractOutgoingVertexA;
  
  // if(decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&  
  //    decay[ianti]->dataPtr()->iColour()==PDT::Colour0){
  //   if     (_abstractOutgoingVertex1) abstractOutgoingVertexF = _abstractOutgoingVertex1;
  //   else if(_abstractOutgoingVertex2) abstractOutgoingVertexF = _abstractOutgoingVertex2;
  // }
  // else if(decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
  // 	  decay[iferm]->dataPtr()->iColour()==PDT::Colour0){
  //   if     (_abstractOutgoingVertex1) abstractOutgoingVertexA = _abstractOutgoingVertex1;
  //   else if(_abstractOutgoingVertex2) abstractOutgoingVertexA = _abstractOutgoingVertex2;
  // }
  // else if(inpart.       dataPtr()->iColour()==PDT::Colour0     &&
  // 	  decay[iferm]->dataPtr()->iColour()==PDT::Colour3     && 
  // 	  decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
  //   abstractOutgoingVertexF = _abstractOutgoingVertex1;
  //   abstractOutgoingVertexA = _abstractOutgoingVertex2;
  // }

  // if (not ((_abstractIncomingVertex  && (abstractOutgoingVertexF || abstractOutgoingVertexA)) ||
  // 	   ( abstractOutgoingVertexF &&  abstractOutgoingVertexA)))
  //   throw Exception()
  //   << "Invalid vertices for QCD radiation in SFF decay in SFFDecayer::threeBodyME"
  //   << Exception::runerror;

  // tcPDPtr Scalar=getParticleData(inpart.id());
  // Energy2 scale(sqr(inpart.mass()));
  // double gs=0;
  // for(unsigned int ifm = 0; ifm < 2; ++ifm) {
  //   for(unsigned int ia = 0; ia < 2; ++ia) {
  //     for(unsigned int iv = 0; iv < 2; ++iv) {
  // 	// radiation from the incoming scalar
  // 	if(inpart.dataPtr()->coloured()) {
  // 	  assert(_abstractIncomingVertex);

  // 	  ScalarWaveFunction scalarInter = 
  // 	    _abstractIncomingVertex->evaluate(scale,3,inpart.dataPtr(),
  // 					      _vwave3[2*iv],_swave3,inpart.mass());

  // 	  if (_swave3.particle()->PDGName()!=scalarInter.particle()->PDGName())
  // 	    throw Exception()
  // 	      << _swave3    .particle()->PDGName() << " was changed to " 
  // 	      << scalarInter.particle()->PDGName() << " in SFFDecayer::threeBodyME"
  // 	      << Exception::runerror;

  // 	  gs = _abstractIncomingVertex->strongCoupling(scale);
  // 	  double sign  = inpart.dataPtr()->id()>0 ? 1:-1;
  // 	  Complex diag = sign * _abstractVertex->evaluate(scale,_wave3[ia],
  // 							  _wavebar3[ifm],scalarInter)/gs;
  // 	  for(unsigned int ix=0;ix<colourFlows()[0].size();++ix) {
  // 	    ME[colourFlows()[0][ix].first](0, ia, ifm, iv) += colourFlows()[0][ix].second*diag;
  // 	  }
  // 	}	

  // 	//radiation from outgoing fermion
  // 	if(decay[iferm]->dataPtr()->coloured()) {
  // 	  assert(abstractOutgoingVertexF);
  // 	  tcPDPtr off = decay[iferm]->dataPtr();
  // 	  if(off->CC()) off = off->CC();
  // 	  SpinorBarWaveFunction interS = 
  // 	    abstractOutgoingVertexF->evaluate(scale,3,off,_wavebar3[ifm],
  // 					       _vwave3[2*iv],decay[iferm]->mass());
	  
  // 	  if(_wavebar3[ifm].particle()->PDGName()!=interS.particle()->PDGName())
  // 	    throw Exception()
  // 	      << _wavebar3[ifm].particle()->PDGName() << " was changed to " 
  // 	      << interS        .particle()->PDGName() << " in SFFDecayer::threeBodyME"
  // 	      << Exception::runerror;

  // 	  gs = abstractOutgoingVertexF->strongCoupling(scale);
  // 	  Complex diag = _abstractVertex->evaluate(scale,_wave3[ia],
  // 						   interS,_swave3)/gs;
  // 	  for(unsigned int ix=0;ix<colourFlows()[1].size();++ix) {
  // 	    ME[colourFlows()[1][ix].first](0, ia, ifm, iv) += colourFlows()[1][ix].second*diag;
  // 	  }
  // 	}

  // 	//radiation from outgoing antifermion
  // 	if(decay[ianti]->dataPtr()->coloured()) {
  // 	  assert(abstractOutgoingVertexA);
  // 	  tcPDPtr off = decay[ianti]->dataPtr();
  // 	  if(off->CC()) off = off->CC();
  // 	  SpinorWaveFunction  interS = 
  // 	    abstractOutgoingVertexA->evaluate(scale,3,off,_wave3[ia],
  // 					      _vwave3[2*iv],decay[ianti]->mass());

  // 	  if(_wave3[ia].particle()->PDGName()!=interS.particle()->PDGName())
  // 	    throw Exception()
  // 	      << _wave3[ia].particle()->PDGName() << " was changed to " 
  // 	      << interS    .particle()->PDGName() << " in SFFDecayer::threeBodyME"
  // 	      << Exception::runerror;

  // 	  gs = abstractOutgoingVertexA->strongCoupling(scale);
  // 	  Complex diag =  _abstractVertex->evaluate(scale,interS,
  // 						    _wavebar3[ifm],_swave3)/gs;
  // 	  for(unsigned int ix=0;ix<colourFlows()[2].size();++ix) {
  // 	    ME[colourFlows()[2][ix].first](0, ia, ifm, iv) += colourFlows()[2][ix].second*diag; 
  // 	  }
  // 	}
  //     }
  //   }
  // }

  // //colour and identical particle factors
  // double output=0.;
  // for(unsigned int ix=0; ix<nflow; ++ix){
  //   for(unsigned int iy=0; iy<nflow; ++iy){
  //     output+=cfactors[ix][iy]*(ME[ix].contract(ME[iy],_rho3)).real();
  //   }
  // }
  // output*=(4.*Constants::pi);
  // // return the answer
  // return output;

  return 0.;
}
