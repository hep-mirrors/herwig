// non-inlined functions of SpinorWaveFunction class
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpinorWaveFunction class.
//
// Author: Peter Richardson
//

#include "SpinorWaveFunction.h"

namespace Herwig {
namespace Helicity {

// calculate the Wavefunction
void SpinorWaveFunction::calculateWaveFunction(int idirac,int ihel)
{
  // check ihelicity is O.K.
  int ipart = getDirection();
  if(ihel!=1 && ihel!=-1)
    {
      cerr << "Invalid Helicity = " << ihel << " requested for Spinor" << endl;
      for(int ix=0;ix<4;++ix){_wf[ix]=0.;}
    }
  else
    {
      // extract the momentum components
      Energy ppx=-ipart*px(),ppy=-ipart*py(),ppz=-ipart*pz(),pee=-ipart*e(),pmm=mass();
      // define and calculate some kinematic quantities
      Energy ptran  = ppx*ppx+ppy*ppy;
      Energy pabs   = sqrt(ptran+ppz*ppz);
      ptran  = sqrt(ptran);
      // first need to evalulate the 2-component helicity spinors 
      // this is the same regardless of which definition of the spinors
      // we are using
      complex <double> hel_wf[2];
      // compute the + spinor for + helicty particles and - helicity antiparticles
      if((ipart==-1 && ihel== 1) || (ipart==1 && ihel==-1))
	{
	  // no transverse momentum 
	  if(ptran==0.)
	    {
	      if(ppz>=0)
		{
		  hel_wf[0] = 1;
		  hel_wf[1] = 0;
		}
	      else
		{
		  hel_wf[0] = 0;
		  hel_wf[1] = 1;
		}
	    }
	  else
	    {
	      double denominator = 1./sqrt(2.*pabs);
	      double rtppluspz;
	      if(ppz>=0)
		{rtppluspz = sqrt(pabs+ppz);}
	      else
		{rtppluspz = ptran/sqrt(pabs-ppz);} 
	      hel_wf[0] = denominator*rtppluspz;
	      hel_wf[1] = denominator/rtppluspz*
		complex<double>(ppx,ppy);
	    }
	}
      // compute the - spinor for - helicty particles and + helicity antiparticles
      else
	{
	  // no transverse momentum
	  if(ptran==0.)
	    {
	      if(ppz>=0.)
		{
		  hel_wf[0] = 0;
		  hel_wf[1] = 1;
		}
	  // transverse momentum 
	      else
		{
		  hel_wf[0] = -1;
		  hel_wf[1] =  0;
		}
	    }
	  else
	    {
	      double denominator = 1./sqrt(2.*pabs);
	      double rtppluspz;
	      if(ppz>=0.)
		{rtppluspz = sqrt(pabs+ppz);}
	      else
		{rtppluspz = ptran/sqrt(pabs-ppz);}
	      hel_wf[0] = denominator/rtppluspz*
		complex<double>(-ppx,ppy);
	      hel_wf[1] = denominator*rtppluspz;
	    }
	}
      // decide which definition of the spinors we are using
      Energy eplusm,eminusm,eplusp,eminusp,upper,lower;
      switch(idirac)
	{
	  // Haber lower energy
	case 1:
	  eplusm = sqrt(pee+pmm);
	  eminusm = pabs/eplusm;
	  if(ipart==-1)
	    {
	      upper = eplusm;
	      if(ihel==1) 
		{lower = eminusm;}
	      else
		{lower =-eminusm;}
	    }
	  else
	    {
	      upper = eminusm;
	      if(ihel==1)
		{lower =-eplusm;}
	      else
		{lower = eplusm;}
	    }
	  break;
	case 2:
	  // HELAS
	  eplusp = sqrt(pee+pabs);
	  if(pmm!=0.) 
	    {eminusp=pmm/eplusp;}
	  else
	    {eminusp=0.;}
	  // set up the coefficients for the different cases
	  if(ipart==-1)
	    {
	      if(ihel==1)
		{
		  upper = eminusp;
		  lower = eplusp;
		}
	      else
		{
		  upper = eplusp;
		  lower = eminusp;
		}
	    }
	  else
	    {
	      if(ihel==1)
		{
		  upper = -eplusp;
		  lower = eminusp;
		}
	      else
		{
		  upper = eminusp;
		  lower =-eplusp;
		}
	    }
	  break;
	  // invalid choice
	default:
	  cerr << "Invalid choice of Dirac representation in "
	       << "SpinorWaveFunction::calculateWaveFunction() " << endl; 
	  break;
	}
      // now finally we can construct the spinors
      if(ipart<0)
	{_wf=LorentzSpinor(idirac,1,upper*hel_wf[0],upper*hel_wf[1],lower*hel_wf[0],
			   lower*hel_wf[1]);}
      else
	{_wf=LorentzSpinor(idirac,2,upper*hel_wf[0],upper*hel_wf[1],lower*hel_wf[0],
			   lower*hel_wf[1]);}
    }
}

}
}
