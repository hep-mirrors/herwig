// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonOctetOctetScalarDecayer class.
//

#include "SU3BaryonOctetOctetScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SU3BaryonOctetOctetScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


namespace Herwig {
using namespace ThePEG;

SU3BaryonOctetOctetScalarDecayer::~SU3BaryonOctetOctetScalarDecayer() {}

bool SU3BaryonOctetOctetScalarDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed=false;
  if(_incomingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0=dm.parent()->id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  unsigned int ix=0;
  do
    {
      if(id0==_incomingB[ix])
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){allowed=true;}
	}
      else if(id0==-_incomingB[ix])
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){allowed=true;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){allowed=true;}
	}
      ++ix;
    }
  while(ix<_incomingB.size()&&!allowed);
  return allowed;
}

ParticleVector SU3BaryonOctetOctetScalarDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int imode=-1;
  int id=parent.id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  unsigned int ix=0;bool cc;
  do 
    {
      if(id==_incomingB[ix])
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){imode=ix;cc=false;}
	}
      else if(id==-_incomingB[ix])
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){imode=ix;cc=true;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){imode=ix;cc=true;}
	}
      ++ix;
    }
  while(ix<_incomingB.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}


void SU3BaryonOctetOctetScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _sf << _sd << _parity << _fpi << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _eproton << _eneutron << _esigma0 
     << _esigmap << _esigmam << _elambda << _exi0 << _exim << _incomingB << _outgoingB 
     << _outgoingM << _maxweight << _A1 << _A2 << _A3 << _B1 << _B2 << _B3;
}

void SU3BaryonOctetOctetScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _sf >> _sd >> _parity >> _fpi >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _eproton >> _eneutron >> _esigma0 
     >> _esigmap >> _esigmam >> _elambda >> _exi0 >> _exim >> _incomingB >> _outgoingB 
     >> _outgoingM >> _maxweight >> _A1 >> _A2 >> _A3 >> _B1 >> _B2 >> _B3;
}

ClassDescription<SU3BaryonOctetOctetScalarDecayer> SU3BaryonOctetOctetScalarDecayer::initSU3BaryonOctetOctetScalarDecayer;
// Definition of the static class description member.

void SU3BaryonOctetOctetScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonOctetOctetScalarDecayer> documentation
    ("The \\classname{SU3BaryonOctetOctetScalarDecayer} class is designed for the"
     " decay of excited baryon resonances assuming SU(3) symmetry");

  static Parameter<SU3BaryonOctetOctetScalarDecayer,double> interfaceFcoupling
    ("Fcoupling",
     "The F coupling of the baryon resonances",
     &SU3BaryonOctetOctetScalarDecayer::_sf, 0.0, -10.0, 10.0,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,double> interfaceDcoupling
    ("Dcoupling",
     "The D coupling of the baryon resonances",
     &SU3BaryonOctetOctetScalarDecayer::_sd, 0.0, -10.0, 10.0,
     false, false, true);

  static Switch<SU3BaryonOctetOctetScalarDecayer,int> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonOctetOctetScalarDecayer::_parity, 0, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "The multiplets have the same parity.",
     0);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "The multiplets have different parities.",
     1);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonOctetOctetScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the lighter proton-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_proton, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the lighter neutron-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_neutron, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the lighter Sigma+-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_sigmap, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_sigma0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the lighter Sigma--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_sigmam, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_lambda, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the lighter Xi0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_xi0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the lighter Xi--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_xim, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedProton
    ("ExcitedProton",
     "The PDG code for the heavier proton-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_eproton, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedNeutron
    ("ExcitedNeutron",
     "The PDG code for the heavier neutron-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_eneutron, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedSigmap
    ("ExcitedSigma+",
     "The PDG code for the heavier Sigma+-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_esigmap, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedSigma0
    ("ExcitedSigma0",
     "The PDG code for the heavier Sigma0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_esigma0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedSigmam
    ("ExcitedSigma-",
     "The PDG code for the heavier Sigma--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_esigmam, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_elambda, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedXi0
    ("ExcitedXi0",
     "The PDG code for the heavier Xi0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_exi0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedXim
    ("ExcitedXi-",
     "The PDG code for the heavier Xi--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_exim, 0, -100000, 100000,
     false, false, true);

  static ParVector<SU3BaryonOctetOctetScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonOctetOctetScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);
}

// couplings for spin-1/2 to spin-1/2 spin-0
void SU3BaryonOctetOctetScalarDecayer::halfHalfScalarCoupling(int imode,Complex& A,Complex& B) const
 {A=_A1[imode];B=_B1[imode];}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonOctetOctetScalarDecayer::halfThreeHalfScalarCoupling(int imode,Complex& A,Complex& B) const
{A=_A1[imode];B=_B1[imode];}



// set up the decay modes
void SU3BaryonOctetOctetScalarDecayer::setupModes(unsigned int iopt) const
{
  if(_incomingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.resize(0);_incomingB.resize(0);_outgoingM.resize(0);}
  // set up for the various different decay modes
  vector<double> factor;
  vector<int> intemp,outtemp,mestemp;
  double ort(1./sqrt(2.)),ors(1./sqrt(6.)),rt(sqrt(2.));
  // decays of the excited proton
  intemp.push_back(_eproton);outtemp.push_back(_neutron);mestemp.push_back(211);
  factor.push_back((_sd+_sf));
  intemp.push_back(_eproton);outtemp.push_back(_proton);mestemp.push_back(111);
  factor.push_back(-ort*(_sd+_sf));
  intemp.push_back(_eproton);outtemp.push_back(_sigmap);mestemp.push_back(311);
  factor.push_back(_sd-_sf);
  intemp.push_back(_eproton);outtemp.push_back(_sigma0);mestemp.push_back(321);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_eproton);outtemp.push_back(_proton);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_eproton);outtemp.push_back(_lambda);mestemp.push_back(321);
  factor.push_back(-ors*(_sd+3.*_sf));
  // decays of the excited neutron
  intemp.push_back(_eneutron);outtemp.push_back(_proton);mestemp.push_back(-211);
  factor.push_back((_sd+_sf));
  intemp.push_back(_eneutron);outtemp.push_back(_neutron);mestemp.push_back(111);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_eneutron);outtemp.push_back(_sigmam);mestemp.push_back(321);
  factor.push_back(_sd-_sf);
  intemp.push_back(_eneutron);outtemp.push_back(_sigma0);mestemp.push_back(311);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_eneutron);outtemp.push_back(_neutron);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_eneutron);outtemp.push_back(_lambda);mestemp.push_back(311);
  factor.push_back(-ors*(_sd+3.*_sf));
  // decays of the excited lambda
  intemp.push_back(_elambda);outtemp.push_back(_sigma0);mestemp.push_back(111);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_sigmap);mestemp.push_back(-211);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_sigmam);mestemp.push_back(211);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_lambda);mestemp.push_back(221);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_xim);mestemp.push_back(321);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_elambda);outtemp.push_back(_xi0);mestemp.push_back(311);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_elambda);outtemp.push_back(_proton);mestemp.push_back(-321);
  factor.push_back(-ors*(3.*_sf+_sd));
  intemp.push_back(_elambda);outtemp.push_back(_neutron);mestemp.push_back(-311);
  factor.push_back(-ors*(3.*_sf+_sd));
  // decays of the excited sigma+
  intemp.push_back(_esigmap);outtemp.push_back(_sigmap);mestemp.push_back(111);
  factor.push_back(rt*_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_sigma0);mestemp.push_back(211);
  factor.push_back(-rt*_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_xi0);mestemp.push_back(321);
  factor.push_back(_sd+_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_proton);mestemp.push_back(-311);
  factor.push_back(_sd-_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_lambda);mestemp.push_back(211);
  factor.push_back(2.*ors*_sd);
  intemp.push_back(_esigmap);outtemp.push_back(_sigmap);mestemp.push_back(221);
  factor.push_back(2.*ors*_sd);
  // decays of the excited sigma0
  intemp.push_back(_esigma0);outtemp.push_back(_sigmam);mestemp.push_back(211);
  factor.push_back(rt*_sf);
  intemp.push_back(_esigma0);outtemp.push_back(_sigmap);mestemp.push_back(-211);
  factor.push_back(-rt*_sf);
  intemp.push_back(_esigma0);outtemp.push_back(_xim);mestemp.push_back(321);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_xi0);mestemp.push_back(311);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_proton);mestemp.push_back(-321);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_neutron);mestemp.push_back(-311);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_sigma0);mestemp.push_back(221);
  factor.push_back(2.*ors*_sd);
  intemp.push_back(_esigma0);outtemp.push_back(_lambda);mestemp.push_back(111);
  factor.push_back(2.*ors*_sd);
  // decays of the excited simga-
  intemp.push_back(_esigmam);outtemp.push_back(_sigma0);mestemp.push_back(-211);
  factor.push_back(rt*_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_sigmam);mestemp.push_back(111);
  factor.push_back(-rt*_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_xim);mestemp.push_back(311);
  factor.push_back(_sd+_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_neutron);mestemp.push_back(-321);
  factor.push_back(_sd-_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_lambda);mestemp.push_back(-211);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_esigmam);outtemp.push_back(_sigmam);mestemp.push_back(221);
  factor.push_back(2.*_sd*ors);
  // decays of the excited xi-
  intemp.push_back(_exim);outtemp.push_back(_sigmam);mestemp.push_back(-311);
  factor.push_back(_sd+_sf);
  intemp.push_back(_exim);outtemp.push_back(_sigma0);mestemp.push_back(-321);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_exim);outtemp.push_back(_xi0);mestemp.push_back(-211);
  factor.push_back(_sd-_sf);
  intemp.push_back(_exim);outtemp.push_back(_xim);mestemp.push_back(111);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_exim);outtemp.push_back(_lambda);mestemp.push_back(-321);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_exim);outtemp.push_back(_xim);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf+_sd));
  // decays of the excited xi0
  intemp.push_back(_exi0);outtemp.push_back(_sigmap);mestemp.push_back(-321);
  factor.push_back(_sd+_sf);
  intemp.push_back(_exi0);outtemp.push_back(_sigma0);mestemp.push_back(-311);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_exi0);outtemp.push_back(_xim);mestemp.push_back(211);
  factor.push_back(_sd-_sf);
  intemp.push_back(_exi0);outtemp.push_back(_xi0);mestemp.push_back(111);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_exi0);outtemp.push_back(_lambda);mestemp.push_back(-311);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_exi0);outtemp.push_back(_xi0);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf+_sd));
  Energy m0,m1;
  int inspin,outspin;
  for(unsigned int ix=0;ix<intemp.size();++ix)
    {
      if(intemp[ix]!=0&&outtemp[ix]!=0&&mestemp[ix]!=0)
	{
	  _incomingB.push_back(intemp[ix]);
	  _outgoingB.push_back(outtemp[ix]);
	  _outgoingM.push_back(mestemp[ix]);
	  if(iopt==1)
	    {
	      m0 = getParticleData(_incomingB.back())->mass();
	      m1 = getParticleData(_outgoingB.back())->mass();
	      inspin  = getParticleData(_incomingB.back())->iSpin();
	      outspin = getParticleData(_outgoingB.back())->iSpin();
	      if(inspin==2&&outspin==2)
		{
		  if(_parity==0)
		    {
		      _A1.push_back(0.);
		      _B1.push_back(factor[ix]*rt/_fpi*(m0+m1));
		    }
		  else
		    {
		      _A1.push_back(factor[ix]*rt/_fpi*(m0-m1));
		      _B1.push_back(0.);
		    }
		  _A2.push_back(0.);_B2.push_back(0.);
		  _A3.push_back(0.);_B3.push_back(0.);
		}
	      else if(inspin==4&&outspin==2)
		{
		  if(_parity==0)
		    {
		      _A1.push_back(factor[ix]*rt/_fpi*(m0+m1));
		      _B1.push_back(0.);
		    }
		  else
		    {
		      _A1.push_back(0.);
		      _B1.push_back(factor[ix]*rt/_fpi*(m0+m1));
		    }
		  _A2.push_back(0.);_B2.push_back(0.);
		  _A3.push_back(0.);_B3.push_back(0.);
		}
	      else
		{throw DecayIntegratorError() << "Invalid combination of spins in "
					      << "SU3BaryonOctetOctetScalarDecayer::" 
					   << "setupModes()" << Exception::abortnow;}
	    }
	}
    }
}

}
