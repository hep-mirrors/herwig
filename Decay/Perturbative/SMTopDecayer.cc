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

SMTopDecayer::SMTopDecayer() 
  :_wquarkwgt(6,0.),_wleptonwgt(3,0.) 
{
  _wleptonwgt[0] = 0.434654;
  _wleptonwgt[1] = 0.431726;
  _wleptonwgt[2] = 0.426912;
  _wquarkwgt[0]  = 1.21429;
  _wquarkwgt[1]  = 0.0654597;
  _wquarkwgt[2]  = 0.0666334;
  _wquarkwgt[3]  = 1.23453;
  _wquarkwgt[4]  = 5.72811e-06;
  _wquarkwgt[5]  = 0.000705673;
  generateIntermediates(true);
}
  
bool SMTopDecayer::accept(const DecayMode & dm) const {
  if(abs(dm.parent()->id()) != ParticleID::t) return false;
  int id0(0),id1(0),id2(0);
  for(ParticleMSet::const_iterator it = dm.products().begin();
      it != dm.products().end();++it) {
    int id=(**it).id(),absid(abs(id));
    if(absid==ParticleID::b&&double(id)/double(dm.parent()->id())>0) {
      id0=id;
    }
    else {
      switch (absid) {
      case ParticleID::nu_e: 
      case ParticleID::nu_mu:
      case ParticleID::nu_tau:
	id1 = id;
	break;
      case ParticleID::eminus:
      case ParticleID::muminus:
      case ParticleID::tauminus:
	id2 = id;
	break;
      case ParticleID::d:
      case ParticleID::s:
	id1 = id;
	break;
      case ParticleID::u:
      case ParticleID::c:
	id2=id;
	break;
      default :
	break;
      }
    }
  }
  if(id0==0||id1==0||id2==0) return false;
  if(double(id1)/double(id2)>0) return false;
  return true;
}
  
ParticleVector SMTopDecayer::decay(const DecayMode & dm,
				   const Particle & parent) const {
  int id1(0),id2(0);
  for(ParticleMSet::const_iterator it = dm.products().begin();
      it != dm.products().end();++it) {
    int id=(**it).id(),absid=abs(id);
    if(absid == ParticleID::b && double(id)/double(parent.id())>0) continue;
    //leptons
    if(absid > 10 && absid%2==0) id1=absid;
    if(absid > 10 && absid%2==1) id2=absid;
    //quarks
    if(absid < 10 && absid%2==0) id2=absid;
    if(absid < 10 && absid%2==1) id1=absid;
  }
  unsigned int imode(0);
  if(id2 >=11 && id2<=16) imode = (id1-12)/2;
  else imode = id1+1+id2/2;
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

double SMTopDecayer::me2(bool vertex, const int, 
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
	      topMe(tbhel,bbhel,fhel,afhel) = 
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

void SMTopDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  
  //get vertices from SM object
  tcHwSMPtr hwsm = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(hwsm)
    {
      _wvertex = hwsm->vertexFFW();
      //initialise
      _wvertex->init();
    }
  else{throw InitException();}
  //set up decay modes
  tPDPtr Wplus(getParticleData(ParticleID::Wplus));
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr Wchannel;
  PDVector extpart(4);
  vector<double> wgt(1,1.0);
  extpart[0] = getParticleData(ParticleID::t);
  extpart[1] = getParticleData(ParticleID::b);
  
  //lepton modes
  for(int i=11; i<17;i+=2) {
    extpart[2] = getParticleData(-i);
    extpart[3] = getParticleData(i+1);
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    Wchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    Wchannel->addIntermediate(extpart[0],0,0.0,-1,1);
    Wchannel->addIntermediate(Wplus,0,0.0,2,3);
    Wchannel->init();
    mode->addChannel(Wchannel);
    addMode(mode,_wleptonwgt[(i-11)/2],wgt);
  }
  //quark modes
  unsigned int iz=0;
  for(int ix=1;ix<6;ix+=2)
    {
      for(int iy=2;iy<6;iy+=2)
 	{
 	  // check that the combination of particles is allowed
 	  if(_wvertex->allowed(-ix,iy,ParticleID::Wminus))
 	    {
 	      extpart[2] = getParticleData(-ix);
 	      extpart[3] = getParticleData( iy);
	      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
	      Wchannel = new_ptr(DecayPhaseSpaceChannel(mode));
	      Wchannel->addIntermediate(extpart[0],0,0.0,-1,1);
	      Wchannel->addIntermediate(Wplus,0,0.0,2,3);
	      Wchannel->init();
	      mode->addChannel(Wchannel);
 	      addMode(mode,_wquarkwgt[iz],wgt);
 	      ++iz;
 	    }
 	  else
 	    {throw InitException() << "SMTopDecayer::doinit() the W vertex" 
 				   << "cannot handle all the quark modes" 
 				   << Exception::abortnow;}
	}
    }
}

}

