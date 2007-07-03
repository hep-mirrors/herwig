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
  if(ihel==0 || ihel==2||(ihel==1&&mass()>Energy()))
    {
      int jhel=ihel-1;
      // extract the momentum components
      double fact=-1.; if(dir==incoming){fact=1.;}
      Energy ppx=fact*px(),ppy=fact*py(),ppz=fact*pz(),pee=fact*e(),pmm=mass();
      // calculate some kinematic quantites;
      Energy2 pt2 = ppx*ppx+ppy*ppy;
      Energy pabs = sqrt(pt2+ppz*ppz);
      Energy pt = sqrt(pt2);
      // overall phase of the vector
      Complex phase;
      if(vphase==vector_phase)
	{
	  if(pt==Energy() || ihel==1){phase = 1.;}
	  else if(ihel==0)     {phase = Complex(ppx/pt,-fact*ppy/pt);}
	  else                 {phase = Complex(ppx/pt, fact*ppy/pt);}
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
	  if(pt==Energy())
	    {
	      double sgnz;
	      if(ppz<Energy()){sgnz=-1.;}
	      else{sgnz=1.;}
	      _wf.setX(-complex<double>(jhel)*phase);
	      _wf.setY(sgnz*phase*complex<double>(0,-fact));
	      _wf.setZ(0.);
	      _wf.setT(0.);
	    }
	  else
	    {
	      InvEnergy opabs=1./pabs;
	      InvEnergy opt  =1./pt;
	      _wf.setX(phase*complex<double>(-jhel*ppz*ppx*opabs*opt, fact*ppy*opt));
	      _wf.setY(phase*complex<double>(-jhel*ppz*ppy*opabs*opt,-fact*ppx*opt));
	      _wf.setZ(double(jhel*pt*opabs)*phase);
	      _wf.setT(0.);
	    }
	}
      // 0 component for massive vectors
      else
	{
	  if(pabs==Energy())
	    {
	      _wf.setX(0.);
	      _wf.setY(0.);
	      _wf.setZ(1.);
	      _wf.setT(0.);
	    }
	  else
	    {
	      InvEnergy empabs=pee/pmm/pabs;
	      _wf.setX(double(empabs*ppx));
	      _wf.setY(double(empabs*ppy));
	      _wf.setZ(double(empabs*ppz));
	      _wf.setT(double(pabs/pmm));
	    }
	}
    }
  // special return the momentum as a check of gauge invariance
  else if(ihel==10)
    {
      _wf.setX(double(px()/MeV));
      _wf.setY(double(py()/MeV));
      _wf.setZ(double(pz()/MeV));
      _wf.setT(double(e()/MeV));
    }
  // issue warning and return zero
  else
    {
      ThePEG::Helicity::HelicityConsistencyError() 
	<< "Invalid Helicity = " << ihel << " requested for Vector " 
	<< getParticle()->PDGName() << Exception::abortnow;
      for(int ix=0;ix<4;++ix){_wf.setX(0.);_wf.setY(0.);_wf.setZ(0.);_wf.setT(0.);}
    }
}

}
}
