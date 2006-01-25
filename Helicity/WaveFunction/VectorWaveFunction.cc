// non-inlined functions of VectorWaveFunction class
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorWaveFunction class.
//
// Author: Peter Richardson
//

#include "VectorWaveFunction.h"

namespace Herwig {
namespace Helicity {

// calculate the Wavefunction
void VectorWaveFunction::calculateWaveFunction(unsigned int ihel,VectorPhase vphase)
{
  Direction dir=direction();
  if(dir==intermediate)
    {throw ThePEG::Helicity::HelicityConsistencyError() 
	<< "In VectorWaveFunction::calcluateWaveFunction "
	<< "particle must be incoming or outgoing not intermediate" 
	<< Exception::abortnow;}
  // check a valid helicity combination
  if(ihel==0 || ihel==2||(ihel==1&&mass()>0.))
    {
      int jhel=ihel-1;
      // extract the momentum components
      double fact=-1.; if(dir==incoming){fact=1.;}
      Energy ppx=fact*px(),ppy=fact*py(),ppz=fact*pz(),pee=fact*e(),pmm=mass();
      // calculate some kinematic quantites;
      Energy pt = ppx*ppx+ppy*ppy;
      Energy pabs = sqrt(pt+ppz*ppz);
      pt = sqrt(pt);
      // overall phase of the vector
      Complex phase;
      if(vphase==vector_phase)
	{
	  if(pt==0. || ihel==1){phase = 1.;}
	  else if(ihel==0)     {phase = Complex(ppx,-fact*ppy)/pt;}
	  else                 {phase = Complex(ppx, fact*ppy)/pt;}
	}
      else
	{
	  phase = 1.;
	}
      if(ihel!=1) phase = phase/sqrt(2.);
      // first the +/-1 helicity states
      if(ihel!=1)
	{
	  // first the no pt case
	  if(pt==0.)
	    {
	      double sgnz;
	      if(ppz<0){sgnz=-1.;}
	      else{sgnz=1.;}
	      _wf[0]=-complex<double>(jhel)*phase;
	      _wf[1]= sgnz*phase*complex<double>(0,-fact);
	      _wf[2]=0.;
	      _wf[3]=0.;
	    }
	  else
	    {
	      double opabs=1./pabs;
	      double opt  =1./pt;
	      _wf[0]=phase*opt*Complex(-jhel*ppz*ppx*opabs, fact*ppy);
	      _wf[1]=phase*opt*Complex(-jhel*ppz*ppy*opabs,-fact*ppx);
	      _wf[2]=jhel*pt*opabs*phase;
	      _wf[3]=0.;
	    }
	}
      // 0 component for massive vectors
      else
	{
	  if(pabs==0)
	    {
	      _wf[0] = 0.;
	      _wf[1] = 0.;
	      _wf[2] = 1.;
	      _wf[3] = 0.;
	    }
	  else
	    {
	      double empabs=pee/pmm/pabs;
	      _wf[0] = empabs*ppx;
	      _wf[1] = empabs*ppy;
	      _wf[2] = empabs*ppz;
	      _wf[3] = pabs/pmm;
	    }
	}
    }
  // special return the momentum as a check of gauge invariance
  else if(ihel==10)
    {
      _wf[0] = px();
      _wf[1] = py();
      _wf[2] = pz();
      _wf[3] = e();
    }
  // issue warning and return zero
  else
    {
      ThePEG::Helicity::HelicityConsistencyError() 
	<< "Invalid Helicity = " << ihel << " requested for Vector " 
	<< getParticle()->PDGName() << Exception::abortnow;
      for(int ix=0;ix<4;++ix){_wf[ix]=0.;}
    }
}

}
}
