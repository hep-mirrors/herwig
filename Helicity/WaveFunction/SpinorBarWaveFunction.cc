// non-inlined functions of SpinorBarWaveFunction class
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpinorBarWaveFunction class.
//
// Author: Peter Richardson
//

#include "SpinorBarWaveFunction.h"

namespace Herwig {
namespace Helicity {
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::SpinorType;
using ThePEG::Helicity::u_spinortype;
using ThePEG::Helicity::v_spinortype;

// calculate the Wavefunction
void SpinorBarWaveFunction::calculateWaveFunction(int ihel,DiracRep dirac)
{
  Direction dir=direction();
  if(dir==intermediate){cerr << "In SpinorBarWaveFunction::calcluateWaveFunction "
			     << "particle must be incoming or outgoing not intermediate" 
			     << endl;}
  // check ihelicity is O.K.
  if(ihel!=1 && ihel!=-1)
    {
 
      cerr << "Invalid Helicity = " << ihel << " requested for SpinorBar" << endl;
      for(int ix=0;ix<4;++ix){_wf[ix]=0.;}
    }
  else
    {
      // extract the momentum components
      double fact=-1.; if(dir==incoming){fact=1.;}
      Energy ppx=fact*px(),ppy=fact*py(),ppz=fact*pz(),pee=fact*e(),pmm=mass();
      // define and calculate some kinematic quantities
      Energy ptran  = ppx*ppx+ppy*ppy;
      Energy pabs   = sqrt(ptran+ppz*ppz);
      ptran  = sqrt(ptran);
      // first need to evalulate the 2-component helicity spinors 
      // this is the same regardless of which definition of the spinors
      // we are using
      complex <double> hel_wf[2];
      // compute the + spinor for + helicty particles and - helicity antiparticles
      if((dir==outgoing && ihel== 1) || (dir==incoming && ihel==-1))
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
		complex<double>(ppx,-ppy);
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
		complex<double>(-ppx,-ppy);
	      hel_wf[1] = denominator*rtppluspz;
	    }
	}
      // decide which definition of the spinors we are using
      Energy eplusm,eminusm,upper,lower,eplusp,eminusp;
      switch(dirac)
	{
	  // Haber lower energy
	case HaberDRep:
	  eplusm = sqrt(pee+pmm);
	  eminusm = pabs/eplusm;
	  if(dir==outgoing)
	    {
	      upper = eplusm;
	      if(ihel==1) 
		{lower =-eminusm;}
	      else
		{lower = eminusm;}
	    }
	  else
	    {
	      upper = eminusm;
	      if(ihel==1)
		{lower = eplusm;}
	      else
		{lower =-eplusm;}
	    }
	  break;
	case HELASDRep:
	  // HELAS
	  eplusp = sqrt(pee+pabs);
	  if(pmm!=0.) 
	    {eminusp=pmm/eplusp;}
	  else
	    {eminusp=0.;}
	  // set up the coefficients for the different cases
	  if(dir==outgoing)
	    {
	      if(ihel==1)
		{
		  upper = eplusp;
		  lower = eminusp;
		}
	      else
		{
		  upper = eminusp;
		  lower = eplusp;
		}
	    }
	  else
	    {
	      if(ihel==1)
		{
		  upper = eminusp;
		  lower = -eplusp;
		}
	      else
		{
		  upper =-eplusp;
		  lower = eminusp;
		}
	    }
	  break;
	}
      // now finally we can construct the spinors
      if(dir==incoming)
	{_wf=LorentzSpinorBar(upper*hel_wf[0],upper*hel_wf[1],
			      lower*hel_wf[0],lower*hel_wf[1],v_spinortype,dirac);}
      else
	{_wf=LorentzSpinorBar(upper*hel_wf[0],upper*hel_wf[1],
			      lower*hel_wf[0],lower*hel_wf[1],u_spinortype,dirac);}
    }
}

}
}
