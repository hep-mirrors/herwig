// -*- C++ -*-
//
// SFFDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SFFDecayer class.
//

#include "SFFDecayer.h"
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

IBPtr SFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void SFFDecayer::doinit() {
  _perturbativeVertex        = dynamic_ptr_cast<FFSVertexPtr>        (getVertex());
  _abstractVertex            = dynamic_ptr_cast<AbstractFFSVertexPtr>(getVertex());
  _abstractIncomingVertex    = dynamic_ptr_cast<AbstractVSSVertexPtr>(getIncomingVertex());
  _abstractOutgoingVertex1   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[0]);
  _abstractOutgoingVertex2   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[1]);
  GeneralTwoBodyDecayer::doinit();
}

void SFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex           << _perturbativeVertex 
     << _abstractIncomingVertex   << _abstractOutgoingVertex1
     << _abstractOutgoingVertex2;
}

void SFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex           >> _perturbativeVertex 
     >> _abstractIncomingVertex   >> _abstractOutgoingVertex1
     >> _abstractOutgoingVertex2;
}

ClassDescription<SFFDecayer> SFFDecayer::initSFFDecayer;
// Definition of the static class description member.

void SFFDecayer::Init() {

  static ClassDocumentation<SFFDecayer> documentation
    ("This class implements to decay of a scalar to 2 fermions");

}

double SFFDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,MEOption meopt) const {
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
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
    _swave = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
    ME(DecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half));
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    SpinorBarWaveFunction::
      constructSpinInfo(_wavebar,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(_wave   ,decay[ianti],outgoing,true);
    return 0.;
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(_wavebar,decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(_wave   ,decay[ianti],outgoing);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int ifm = 0; ifm < 2; ++ifm){
    for(unsigned int ia = 0; ia < 2; ++ia) {
      if(iferm > ianti){
	ME()(0, ia, ifm) = _abstractVertex->evaluate(scale,_wave[ia],
						     _wavebar[ifm],_swave);
      }
      else {
	ME()(0, ifm, ia) = _abstractVertex->evaluate(scale,_wave[ia],
						     _wavebar[ifm],_swave);
	
      }
    }
  }
  double output = (ME().contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy SFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    _perturbativeVertex->setCoupling(sqr(inpart.second), outb.first, outa.first,
				     in);
    double mu1(outa.second/inpart.second),mu2(outb.second/inpart.second);
    double c2 = norm(_perturbativeVertex->norm());
    Complex al(_perturbativeVertex->left()), ar(_perturbativeVertex->right());
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
			       const ParticleVector & decay, MEOption meopt) {

  // work out which is the fermion and antifermion
  int ianti(0), iferm(1), iglu(2);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0||itype[1]==1||(itype[0]==2&&itype[1]==2)) swap(iferm,ianti);
   
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho3,const_ptr_cast<tPPtr>(&inpart),incoming);
    _swave3 = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    SpinorBarWaveFunction::
      constructSpinInfo(_wavebar3 ,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(_wave3    ,decay[ianti],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(_vwave3   ,decay[iglu ],outgoing,true,false);

    return 0.;
  }


  //unsigned int nflow = getColourFactors().size();

  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors2(inpart, decay, nflow);

  vector<DecayMatrixElement> ME(nflow,DecayMatrixElement(PDT::Spin0,     PDT::Spin1Half,
							 PDT::Spin1Half, PDT::Spin1));

  SpinorBarWaveFunction::
    calculateWaveFunctions(_wavebar3, decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(_wave3   , decay[ianti],outgoing);
  VectorWaveFunction::
    calculateWaveFunctions(_vwave3  , decay[iglu ],outgoing,true);

  // gauge test
  _vwave3.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) _vwave3.push_back(VectorWaveFunction());
    else {
      _vwave3.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
					   decay[iglu ]->dataPtr(),10,
					   outgoing));
      //const complex<double> xcomp (-decay[iglu]->momentum().x()/GeV, 0);
      //const complex<double> ycomp (-decay[iglu]->momentum().y()/GeV, 0);
      //const complex<double> zcomp (-decay[iglu]->momentum().z()/GeV, 0);
      //const complex<double> tcomp (-decay[iglu]->momentum().t()/GeV, 0);
    
      //_vwave3.push_back(VectorWaveFunction(-decay[iglu ]->momentum(),
      //				   decay[iglu ]->dataPtr(),
      //				   xcomp, ycomp, zcomp, tcomp));
    }
  }
   cerr << inpart.PDGName() << " -> " << decay[iferm]->dataPtr()->PDGName() 
	<< " " << decay[ianti]->dataPtr()->PDGName() << " " 
	<< decay[iglu]->dataPtr()->PDGName() << endl;

   cerr << "vector wavefunction[0] (" << _vwave3[0].x() << ", " << _vwave3[0].y() << ", " 
        << _vwave3[0].z() << ", " << _vwave3[0].t() << ")" << endl;
   cerr << "vector wavefunction[1] (" << _vwave3[1].x() << ", " << _vwave3[1].y() << ", " 
        << _vwave3[1].z() << ", " << _vwave3[1].t() << ")" << endl;
   cerr << "vector wavefunction[2] (" << _vwave3[2].x() << ", " << _vwave3[2].y() << ", " 
        << _vwave3[2].z() << ", " << _vwave3[2].t() << ")" << endl;
   cerr << "gluon momentum " << decay[iglu]->momentum()/GeV << endl;

  // if(_abstractIncomingVertex){
  //   set<tPDPtr>::const_iterator it;
  //   for (it=_abstractIncomingVertex->incoming().begin();
  // 	 it!=_abstractIncomingVertex->incoming().end();++it){
  //     cerr << "testing In incoming " << (**it).PDGName() << endl;
  //   }
  // }

  // if (_abstractIncomingVertex){
  //   tPDPtr inparticle = getParticleData(inpart.id());
  //   tPDPtr gluon = getParticleData(ParticleID::g);
  //   cerr << "Incoming Vertex? 1" << endl;
  //   cerr << "allowed? " << _abstractIncomingVertex->isIncoming(inparticle)
  // 	 << _abstractIncomingVertex->isOutgoing(inparticle) 
  // 	 << _abstractIncomingVertex->isOutgoing(gluon) << endl;
  //   cerr << "norm " << _abstractIncomingVertex->norm() << endl;
  //   cerr << "number of points " << _abstractIncomingVertex->getNpoint() << endl;
  // }


  AbstractFFVVertexPtr abstractOutgoingVertexF;
  AbstractFFVVertexPtr abstractOutgoingVertexA;
  
  if(decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&  
     decay[ianti]->dataPtr()->iColour()==PDT::Colour0){
    if     (_abstractOutgoingVertex1) abstractOutgoingVertexF = _abstractOutgoingVertex1;
    else if(_abstractOutgoingVertex2) abstractOutgoingVertexF = _abstractOutgoingVertex2;
  }
  else if(decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
	  decay[iferm]->dataPtr()->iColour()==PDT::Colour0){
    if     (_abstractOutgoingVertex1) abstractOutgoingVertexA = _abstractOutgoingVertex1;
    else if(_abstractOutgoingVertex2) abstractOutgoingVertexA = _abstractOutgoingVertex2;
  }
  else if((inpart.       dataPtr()->iColour()==PDT::Colour0     &&
	   decay[iferm]->dataPtr()->iColour()==PDT::Colour3     && 
	   decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) &&
	  (-decay[ianti]->dataPtr()->id() == decay[iferm]->dataPtr()->id())){
    abstractOutgoingVertexF = _abstractOutgoingVertex1;
    abstractOutgoingVertexA = _abstractOutgoingVertex2;
  }

  // if(abstractOutgoingVertexF) {
  //   tPDPtr particle = getParticleData(decay[iferm]->dataPtr()->id());
  //   tPDPtr gluon = getParticleData(ParticleID::g);
  //   cerr << abstractOutgoingVertexF->isIncoming(particle)
  // 	 << abstractOutgoingVertexF->isOutgoing(particle) 
  // 	 << abstractOutgoingVertexF->isOutgoing(gluon) << endl;
  // }
  // if(abstractOutgoingVertexA) {
  //   tPDPtr particle = getParticleData(decay[ianti]->dataPtr()->id());
  //   tPDPtr gluon = getParticleData(ParticleID::g);
  //   cerr << abstractOutgoingVertexA->isIncoming(particle)
  // 	 << abstractOutgoingVertexA->isOutgoing(particle) 
  // 	 << abstractOutgoingVertexA->isOutgoing(gluon) << endl;
  // }
       
  cerr << "in momentum " << inpart.momentum()/GeV 
       << inpart.momentum().mass()/GeV << endl;
  cerr << "fermion momentum " << decay[iferm]->momentum()/GeV 
       << decay[iferm]->momentum().mass()/GeV << endl;
  cerr << "antifermion momentum " << decay[ianti]->momentum()/GeV 
       << decay[ianti]->momentum().mass()/GeV << endl;


  tcPDPtr Scalar=getParticleData(inpart.id());
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int ifm = 0; ifm < 2; ++ifm) {
    for(unsigned int ia = 0; ia < 2; ++ia) {
      for(unsigned int iv = 0; iv < 2; ++iv) {
	// radiation from the incoming scalar
	if(inpart.dataPtr()->coloured()) {
	  assert(_abstractIncomingVertex);
	  ScalarWaveFunction scalarInter = 
	    _abstractIncomingVertex->evaluate(scale,3,inpart.dataPtr(),
					      _vwave3[2*iv],_swave3,inpart.mass());
	  Complex diag = _abstractVertex->evaluate(scale,_wave3[ia],
						   _wavebar3[ifm],scalarInter);
	  for(unsigned int ix=0;ix<colourFlows()[0].size();++ix) {
	    if (iferm>ianti){
	      ME[colourFlows()[0][ix].first](0, ia, ifm, iv) += colourFlows()[0][ix].second*diag;
	      cerr << "incoming diag ME(0," << ia << "," << ifm << "," << iv << ") " << colourFlows()[0][ix].second*diag << endl;
	    }
	    else{
	      ME[colourFlows()[0][ix].first](0, ifm, ia, iv) += colourFlows()[0][ix].second*diag;
	      cerr << "incoming diag ME(0," << ifm << "," << ia << "," << iv << ") " << colourFlows()[0][ix].second*diag << endl;
	    }
	    
	  }
	}	
	//radiation from outgoing fermion
	if(decay[iferm]->dataPtr()->coloured()) {
	  assert(abstractOutgoingVertexF);
	  SpinorBarWaveFunction interS = 
	    abstractOutgoingVertexF->evaluate(scale,3,decay[iferm]->dataPtr(),_wavebar3[ifm],
					       _vwave3[2*iv],decay[iferm]->mass());
	  Complex diag = _abstractVertex->evaluate(scale,_wave3[ia],
						   interS,_swave3);
	  for(unsigned int ix=0;ix<colourFlows()[1].size();++ix) {
	    if (iferm>ianti){
	      ME[colourFlows()[1][ix].first](0, ia, ifm, iv) += colourFlows()[1][ix].second*diag;
	      cerr << "fermion diag ME(0," << ia << "," << ifm << "," << iv << ") " << colourFlows()[1][ix].second*diag << endl;
	    }
	    else{
	      ME[colourFlows()[1][ix].first](0, ifm, ia, iv) += colourFlows()[1][ix].second*diag;
	      cerr << "fermion diag ME(0," << ifm << "," << ia << "," << iv << ") " << colourFlows()[1][ix].second*diag << endl;
	    }
	    
	  }
	}
	//radiation from outgoing antifermion
	if(decay[ianti]->dataPtr()->coloured()) {
	  assert(abstractOutgoingVertexA);
	  SpinorWaveFunction  interS = 
	    abstractOutgoingVertexA->evaluate(scale,3,decay[ianti]->dataPtr(),_wave3[ia],
					       _vwave3[2*iv],decay[ianti]->mass());
	  Complex diag = _abstractVertex->evaluate(scale,interS,
						   _wavebar3[ifm],_swave3);
	  for(unsigned int ix=0;ix<colourFlows()[2].size();++ix) {
	    if (iferm>ianti){
	      ME[colourFlows()[2][ix].first](0, ia, ifm, iv) += colourFlows()[2][ix].second*diag;
	      cerr << "antifermion diag ME(0," << ia << "," << ifm << "," << iv << ") " << colourFlows()[2][ix].second*diag << endl;
	    }
	    else{
	      ME[colourFlows()[2][ix].first](0, ifm, ia, iv) += colourFlows()[2][ix].second*diag;
	      cerr << "antifermion diag ME(0," << ifm << "," << ia << "," << iv << ") " << colourFlows()[2][ix].second*diag << endl;
	    }
	    
	  }
	}
      }
    }
  }
  
  // colour and identical particle factors
  //cerr << "colour factors " << getColourFactors()[0][0] << endl;
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix].contract(ME[iy],_rho3)).real();   
    }
  }
  // return the answer
  return output;
}
