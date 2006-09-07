// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarWaveFunction class.
//
// Author: Peter Richardson
//

#include "ScalarWaveFunction.h" 

namespace Herwig {
namespace Helicity {

// special set-up constructor
ScalarWaveFunction::ScalarWaveFunction(tPPtr part,Direction dir,bool time,bool vertex)
{
  if(vertex)
    {
      tScalarSpinPtr inspin;
      if(part->spinInfo()){inspin=dynamic_ptr_cast<tScalarSpinPtr>(part->spinInfo());}
      if(dir==outgoing)
	{if(!inspin){part->spinInfo(new_ptr(ScalarSpinInfo(part->momentum(),time)));}}
      else
	{
	  if(inspin){inspin->decay();}
	  else
	    {
	      if(part->spinInfo())
		{throw ThePEG::Helicity::HelicityConsistencyError() 
		    << "Wrong type of SpinInfo for the incoming particle in "
		    << "ScalarWaveFunction() " << part->PDGName()
		    << Exception::warning;}
	      SpinPtr newspin(new_ptr(ScalarSpinInfo(part->momentum(),time)));
	      inspin= dynamic_ptr_cast<tScalarSpinPtr>(newspin);
	      if(time) inspin->decayed(true);
	      part->spinInfo(newspin);
	    }
	}
    }
  direction(dir);
  setMomentum(part->momentum());
  checkParticle(part->dataPtr());
  _wf=1.;
}

// special set-up constructor
ScalarWaveFunction::ScalarWaveFunction(tPPtr part,RhoDMatrix& rho,Direction dir,
				       bool time,bool vertex)
{
  tScalarSpinPtr inspin;
  if(part->spinInfo())
    {inspin=dynamic_ptr_cast<tScalarSpinPtr>(part->spinInfo());}
  if(vertex)
    {
      if(dir==outgoing)
	{if(!inspin){part->spinInfo(new_ptr(ScalarSpinInfo(part->momentum(),time)));}}
      else
	{
	  if(inspin){inspin->decay();}
	  else
	    {
	      if(part->spinInfo())
		{throw ThePEG::Helicity::HelicityConsistencyError() 
		    << "Wrong type of SpinInfo for the incoming particle in "
		    << "ScalarWaveFunction() "
		    << Exception::warning;}
	      SpinPtr newspin=new_ptr(ScalarSpinInfo(part->momentum(),time));
	      inspin= dynamic_ptr_cast<tScalarSpinPtr>(newspin);
	      if(time) inspin->decayed(true);
	      part->spinInfo(newspin);
	    }
	}
    }
  rho=RhoDMatrix(PDT::Spin0);rho.average();
  // set up the wavefunction
  direction(dir);
  setMomentum(part->momentum());
  _wf=1.;
  checkParticle(part->dataPtr());
}

}
}
