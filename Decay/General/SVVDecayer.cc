// -*- C++ -*-
//
// SVVDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SVVDecayer class.
//

#include "SVVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/GeneralVVSVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr SVVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SVVDecayer::fullclone() const {
  return new_ptr(*this);
}

void SVVDecayer::doinit() {
  GeneralTwoBodyDecayer::doinit();
  _abstractVertex     = dynamic_ptr_cast<AbstractVVSVertexPtr>(getVertex());
  _generalVertex      = dynamic_ptr_cast<GeneralVVSVertexPtr >(getVertex());
  _perturbativeVertex = dynamic_ptr_cast<VVSVertexPtr        >(getVertex());
  GeneralTwoBodyDecayer::doinit();
}

void SVVDecayer::doinitrun() {
  getVertex()->initrun();
  GeneralTwoBodyDecayer::doinitrun();
}

void SVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _generalVertex << _perturbativeVertex;
}

void SVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _generalVertex >> _perturbativeVertex;
}

ClassDescription<SVVDecayer> SVVDecayer::initSVVDecayer;
// Definition of the static class description member.

void SVVDecayer::Init() {

  static ClassDocumentation<SVVDecayer> documentation
    ("This implements the decay of a scalar to 2 vector bosons.");

}

double SVVDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector& decay, 
		       MEOption meopt) const {
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = decay[ix]->mass()==ZERO;
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
    _swave = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
    ME(DecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1));
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::
	constructSpinInfo(_vectors[ix],decay[ix],outgoing,true,photon[ix]);
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(_vectors[ix],decay[ix],outgoing,photon[ix]);
  
  
  Energy2 scale(sqr(inpart.mass()));
  unsigned int iv1,iv2;
  for(iv2 = 0; iv2 < 3; ++iv2) {
    if( photon[1] && iv2 == 1 ) ++iv2;
    for(iv1=0;iv1<3;++iv1) {
      if( photon[0] && iv1 == 1) ++iv1;
      ME()(0, iv1, iv2) = _abstractVertex->evaluate(scale,_vectors[0][iv1],
						    _vectors[1][iv2],_swave);
    }
  }
  double output = ME().contract(_rho).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy SVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    Energy2 scale(sqr(inpart.second));
    _perturbativeVertex->setCoupling(scale, outa.first , 
				    outb.first, inpart.first);
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
    Energy output = norm(_perturbativeVertex->getNorm())*
      me2*pcm/(8*Constants::pi)/scale*UnitRemoval::E2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else if(_generalVertex) {
    Lorentz5Momentum in(ZERO,ZERO,ZERO,inpart.second),out1,out2;
    Kinematics::twoBodyDecay(in,outa.second,outb.second,
			     Axis(0.,0.,1.),out1,out2);
    _generalVertex->calculateKinematics(in,out1,out2);
    Energy2 scale(sqr(inpart.second));
    _generalVertex->setCoupling(scale, outa.first, outb.first, inpart.first);
    //get loop coefficients
    Complex a00(_generalVertex->a00());Complex a11(_generalVertex->a11());
    Complex a12(_generalVertex->a12());Complex a21(_generalVertex->a21());
    Complex a22(_generalVertex->a22());Complex aEp(_generalVertex->aEp());
    double mu1(outa.second/inpart.second),mu1sq(mu1*mu1);
    double mu2(outb.second/inpart.second),mu2sq(mu2*mu2);
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Complex me2 = (aEp*aEp + a11*(2.*a11-2.*a21+a22))*mu1sq*mu1sq
      +2.*a21*a21*mu2sq*mu1sq + 2.*mu1sq*(a11*(a21-a22) - aEp*aEp)
      + mu2sq*mu2sq*(aEp*aEp + a22*(a11-2.*a21+2.*a22)) + aEp*aEp
      + 2.*mu2sq*(aEp*aEp +(aEp*aEp +(a11-a21)*(a21-a22))*mu1sq +a11*a22-a21*a22)
      + a11*a22 + a12*(mu1sq + mu2sq -1)*(-2.*a11*mu1sq +(a21-2.*a22)*mu2sq 
					  +a21*(mu1sq-1));
    me2 *= scale*scale*UnitRemoval::InvE4;
    me2 += 2.*a00*double(scale*UnitRemoval::InvE2)*((2.*a11 - a12 - a21)*mu1sq
						    - mu2sq*(a12+a21-2.*a22)
						    + a21 + a12) + 8.*a00*a00;
    me2 /= 2;
    Energy output = norm(_generalVertex->getNorm())*me2.real()*pcm
      /(8.*Constants::pi)/scale*UnitRemoval::E2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}
