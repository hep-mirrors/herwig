// non-inlined functions of TensorWaveFunction class
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorWaveFunction class.
//
// Author: Peter Richardson
//
#include "TensorWaveFunction.h"

namespace Herwig {
namespace Helicity {

// calculate the actual wavefunction
void TensorWaveFunction::calculateWaveFunction(unsigned int ihel, TensorPhase tphase)
{
  int jhel=ihel-2;
  Direction dir=direction();
  if(dir==intermediate)
    {ThePEG::Helicity::HelicityConsistencyError() 
	<< "In TensorWaveFunction::calcluateWaveFunction "
	<< "particle must be incoming or outgoing not intermediate" 
	<< Exception::abortnow;}
  // check for a valid helicty combination
  if((jhel<=2 && jhel>=-2   && mass() >0.) || 
     ((jhel==2 || jhel==-2) && mass()==0.)) 
    {
      // extract the momentum components
      double fact=-1.; if(dir==incoming){fact=1.;}
      Energy ppx=fact*px(),ppy=fact*py(),ppz=fact*pz(),pee=fact*e(),pmm=mass();
      // calculate some kinematic quantites;
      Energy pt = ppx*ppx+ppy*ppy;
      Energy pabs = sqrt(pt+ppz*ppz);
      pt = sqrt(pt);
      // polarization vectors
      complex<double> epsp[4],epsm[4],eps0[4];
      // + helicity vector if needed
      if(jhel>=0)
	{
	  // calculate the overall phase
	  complex<double>phase;
	  if(tphase==tensor_phase)
	    {
	      if(pt==0.){phase=1.;}
	      else{phase = complex<double>(ppx,-fact*ppy)/pt;}
	    }
	  else{phase = 1.;} 
	  phase = phase/sqrt(2.);
	  // first the no pt case
	  if(pt==0.)
	    {
	      double sgnz;
	      if(ppz<0){sgnz=-1.;}
	      else{sgnz=1.;}
	      epsp[0]=-phase;
	      epsp[1]= sgnz*phase*complex<double>(0,-fact);
	      epsp[2]=0.;
	      epsp[3]=0.;
	    }
	  else
	    {
	      double opabs=1./pabs;
	      double opt  =1./pt;
	      epsp[0]=phase*complex<double>(-ppz*ppx*opabs*opt,
					  fact*ppy*opt);
	      epsp[1]=phase*complex<double>(-ppz*ppy*opabs*opt,
					  -fact*ppx*opt);
	      epsp[2]=pt*opabs*phase;
	      epsp[3]=0.;
	    }
	}
      // - helicity vector if needed
      if(jhel<=0)
	{
	  // calculate the overall phase
	  complex<double> phase;
	  if(tphase==tensor_phase)
	    {
	      if(pt==0.){phase=1.;}
	      else{phase = complex<double>(ppx,fact*ppy)/pt;}
	    }
	  else{phase = 1.;}
	  phase = phase/sqrt(2.);
	  // first the no pt case
	  if(pt==0.)
	    {
	      double sgnz;
	      if(ppz<0){sgnz=-1.;}
	      else{sgnz=1.;}
	      epsm[0]= phase;
	      epsm[1]= sgnz*phase*complex<double>(0,-fact);
	      epsm[2]=0.;
	      epsm[3]=0.;
	    }
	  else
	    {
	      double opabs=1./pabs;
	      double opt  =1./pt;
	      epsm[0]=phase*complex<double>(ppz*ppx*opabs*opt,
					  fact*ppy*opt);
	      epsm[1]=phase*complex<double>(ppz*ppy*opabs*opt,
					   -fact*ppx*opt);
	      epsm[2]=-pt*opabs*phase;
	      epsm[3]=0.;
	    }
	}
      // 0 helicity vector if needed
      if(jhel<=1 && jhel>=-1)
	{
	  if(pabs==0)
	    {
	      eps0[0] = 0.;
	      eps0[1] = 0.;
	      eps0[2] = 1.;
	      eps0[3] = 0.;
	    }
	  else
	    {
	      double empabs=pee/pmm/pabs;
	      eps0[0] = empabs*ppx;
	      eps0[1] = empabs*ppy;
	      eps0[2] = empabs*ppz;
	      eps0[3] = pabs/pmm;
	    }
	}
      // put the polarization vectors together to get the wavefunction
      double ort;
      switch (jhel)
	{ 
	case 2:
	  for(int ix=0;ix<4;++ix)
	    {for(int iy=0;iy<4;++iy){_wf(ix,iy)=epsp[ix]*epsp[iy];}}
	  break;
	case 1:
	  ort = 1./sqrt(2.);
	  for(int ix=0;ix<4;++ix)
	    {for(int iy=0;iy<4;++iy){_wf(ix,iy)=ort*( epsp[ix]*eps0[iy]
						    +eps0[ix]*epsp[iy]);}}
	  break;
	case 0:
	  ort = 1./sqrt(6.);
	  for(int ix=0;ix<4;++ix)
	    {for(int iy=0;iy<4;++iy){_wf(ix,iy)=ort*(    epsp[ix]*epsm[iy]
						    +   epsm[ix]*epsp[iy]
						    +2.*eps0[ix]*eps0[iy]);}}
	  break;
	case -1:
	  ort = 1./sqrt(2.);
	  for(int ix=0;ix<4;++ix)
	    {for(int iy=0;iy<4;++iy){_wf(ix,iy)=ort*( epsm[ix]*eps0[iy]
						    +eps0[ix]*epsm[iy]);}}
	  break;
	case -2:
	  for(int ix=0;ix<4;++ix)
	    {for(int iy=0;iy<4;++iy){_wf(ix,iy)=epsm[ix]*epsm[iy];}}
	  break;
	default:
	  ThePEG::Helicity::HelicityConsistencyError() 
	    << "Invalid Helicity = " << jhel << " requested for Tensor" 
	    << Exception::abortnow;
	  for(int ix=0;ix<4;++ix){for(int iy=0;iy<4;++iy){_wf(ix,iy)=0.;}}
	  break;
	}
    }
  else
    {
      ThePEG::Helicity::HelicityConsistencyError() 
	<< "Invalid Helicity = " << jhel << " requested for Tensor" 
	<< Exception::abortnow;
      for(int ix=0;ix<4;++ix){for(int iy=0;iy<4;++iy){_wf(ix,iy)=0.;}}
    }
}

}
}
