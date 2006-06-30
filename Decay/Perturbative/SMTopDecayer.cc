// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMTopDecayer class.
//

#include "SMTopDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// // #include "SMTopDecayer.tcc"
#endif
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/Correlations/DecayVertex.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"


namespace Herwig {
  using namespace ThePEG;
  using ThePEG::Helicity::RhoDMatrix;
  using Helicity::VectorWaveFunction;
  using Helicity::Direction;
  using Helicity::incoming;
  using Helicity::outgoing;

  bool SMTopDecayer::accept(const DecayMode & dm) const {
    if(abs(dm.parent()->id()) == ParticleID::t)
      return true;
    else
      return false;
  }
  
  
  ParticleVector SMTopDecayer::decay(const DecayMode & dm,
				     const Particle & parent) const {
    int id0(0),id1(0),id2(0);
    for(ParticleMSet::const_iterator it = dm.products().begin();
	it != dm.products().end();++it) {
      if(abs((*it)->id()) == ParticleID::b){id0 = (*it)->id();}
      else {
	//leptons
	if(abs((*it)->id()) == ParticleID::nu_e ||
	   abs((*it)->id()) == ParticleID::nu_mu ||
	   abs((*it)->id()) == ParticleID::nu_tau)
	  {id1 = (*it)->id();}
	
	if(abs((*it)->id()) == ParticleID::eminus ||
	   abs((*it)->id()) == ParticleID::muminus ||
	   abs((*it)->id()) == ParticleID::tauminus)
	  {id2 = (*it)->id();}
	//quarks
	if(abs((*it)->id()) == ParticleID::d||
	   abs((*it)->id()) == ParticleID::s)
	  {id1=(*it)->id();}
	if(abs((*it)->id()) == ParticleID::u||
	   abs((*it)->id()) == ParticleID::c)
	  {id2=(*it)->id();}
      }
    }
    unsigned int imode(0);
    if(abs(id2) >=11 && abs(id2)<=16) { imode = (abs(id1)%3)/2;}
    if(abs(id1) == ParticleID::d && abs(id2) == ParticleID::u){imode = 3;}
    if(abs(id1) == ParticleID::d && (id2) == ParticleID::c){imode = 4;}
    if(abs(id1) == ParticleID::s) {imode = 4 + abs(id2)/2;}
    if(abs(id1) == ParticleID::b) {imode = 6 + abs(id2)/2;}
    bool cc = parent.id() == ParticleID::tbar;
    ParticleVector out(generate(true,cc,imode,parent));
    //arrange colour flow
    PPtr pparent=const_ptr_cast<PPtr>(&parent);
    out[1]->incomingColour(pparent,out[1]->id()<0);
    ParticleVector products = out[0]->children();
    if(products[0]->hasColour())
      products[0]->colourNeighbour(products[1],true);
    else if(products[0]->hasAntiColour())
      products[0]->colourNeighbour(products[1],false);
    return out;
  }    
  void SMTopDecayer::persistentOutput(PersistentOStream & os) const {
    os << _wvertex << _wquarkwgt << _wleptonwgt; 
  }
  
  void SMTopDecayer::persistentInput(PersistentIStream & is, int) {
    is >> _wvertex >> _wquarkwgt >> _wleptonwgt;
  }
  
  ClassDescription<SMTopDecayer> SMTopDecayer::initSMTopDecayer;
  // Definition of the static class description member.
  
  void SMTopDecayer::Init() {
    
    static ClassDocumentation<SMTopDecayer> documentation
      ("This is the implementation of the SMTopDecayer which "
       "decays top quarks into bottom quarks and either leptons  "
       "or quark-antiquark pairs.");
 
  }

  double SMTopDecayer::me2(bool vertex, const int ichan, 
			   const Particle & inpart,
			   const ParticleVector & decay) const
  {
    RhoDMatrix rhot(PDT::Spin1Half);
    rhot.average();   //diagonalise  
    DecayMatrixElement topMe(PDT::Spin1Half,PDT::Spin1Half,
			     PDT::Spin1Half,PDT::Spin1Half);
    Energy2 scale(inpart.mass()*inpart.mass());
    if(inpart.id() == ParticleID::t) {
      //Vectors to hold all heliticies of spinors
      vector<SpinorWaveFunction> twave,awave;
      vector<SpinorBarWaveFunction> bwave,fwave;
      //Set-up spinors for external particles
      SpinorWaveFunction(twave,rhot,const_ptr_cast<tPPtr>(&inpart),
			 incoming,true,vertex);
      SpinorBarWaveFunction(bwave,decay[0],outgoing,true,vertex);
      SpinorWaveFunction(awave,decay[1],outgoing,true,vertex);
      SpinorBarWaveFunction(fwave,decay[2],outgoing,true,vertex);
      //Define intermediate vector wave-function for Wplus 
      tcPDPtr Wplus(getParticleData(ParticleID::Wplus));
      VectorWaveFunction inter;
      unsigned int thel,bhel,fhel,afhel;
      for(thel = 0;thel<2;++thel){
	for(bhel = 0;bhel<2;++bhel){	  
	  inter = _wvertex->evaluate(scale,1,Wplus,twave[thel],bwave[bhel]);
	  for(afhel=0;afhel<2;++afhel){
	    for(fhel=0;fhel<2;++fhel){
	      topMe(thel,bhel,afhel,fhel) = 
		_wvertex->evaluate(scale,awave[afhel],
				   fwave[fhel],inter);
	    }
	  }
	}
      }
    }
    if(inpart.id() == ParticleID::tbar)
      {
	//Vectors to hold all heliticies of spinors
	vector<SpinorWaveFunction> bbarWave,awave;
	vector<SpinorBarWaveFunction> tbarWave,fwave;
	//Set-up spinors for externl particles
	SpinorBarWaveFunction(tbarWave,rhot,const_ptr_cast<tPPtr>(&inpart),
			      incoming,true,vertex);
	SpinorWaveFunction(bbarWave,decay[0],outgoing,true,vertex);
	SpinorBarWaveFunction(fwave,decay[1],outgoing,true,vertex);
	SpinorWaveFunction(awave,decay[2],outgoing,true,vertex);
	VectorWaveFunction inter;
	tcPDPtr Wminus(getParticleData(ParticleID::Wminus));
	unsigned int tbhel,bbhel,afhel,fhel;
	for(tbhel = 0;tbhel<2;++tbhel){
	  for(bbhel = 0;bbhel<2;++bbhel){
	    inter = _wvertex->
	      evaluate(scale,1,Wminus,bbarWave[bbhel],tbarWave[tbhel]);
	    for(afhel=0;afhel<2;++afhel){
	      for(fhel=0;fhel<2;++fhel){
		topMe(tbhel,bbhel,afhel,fhel) = 
		  _wvertex->evaluate(scale,awave[afhel],
				   fwave[fhel],inter);
	      }
	    }
	  }
	}
      }
    
    ME(topMe);
  double output = (topMe.contract(rhot)).real();
  if(abs(decay[1]->id())<=4){output *=3.;}
  return output;
  }
}


