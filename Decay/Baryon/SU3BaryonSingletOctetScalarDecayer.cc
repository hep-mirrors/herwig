// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonSingletOctetScalarDecayer class.
//

#include "SU3BaryonSingletOctetScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SU3BaryonSingletOctetScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

SU3BaryonSingletOctetScalarDecayer::~SU3BaryonSingletOctetScalarDecayer() {}

int SU3BaryonSingletOctetScalarDecayer::modeNumber(bool & cc,const DecayMode & dm) const
{
  int imode(-1);
  if(_outgoingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(dm.products().size()!=2){return imode;}
  // ids of the particles
  int id0(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);
  cc =false;
  do
    {
      if(id0==_elambda)
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){imode=ix;cc=false;}
	}
      else if(id0==-_elambda)
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
  while(ix<_outgoingB.size()&&imode<0);
  return imode;
}

void SU3BaryonSingletOctetScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _C << _parity << _fpi << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _elambda << _outgoingB 
     << _outgoingM << _maxweight << _prefactor;
}

void SU3BaryonSingletOctetScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _C >> _parity >> _fpi >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _elambda >> _outgoingB 
     >> _outgoingM >> _maxweight >> _prefactor;
}

ClassDescription<SU3BaryonSingletOctetScalarDecayer> SU3BaryonSingletOctetScalarDecayer::initSU3BaryonSingletOctetScalarDecayer;
// Definition of the static class description member.

void SU3BaryonSingletOctetScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonSingletOctetScalarDecayer> documentation
    ("The SU3BaryonSingletOctetScalarDecayer class is designed for"
     "the decay of an excited SU(3) singlet baryon");

  static Parameter<SU3BaryonSingletOctetScalarDecayer,double> interfaceCcoupling
    ("Coupling",
     "The C coupling of the baryon resonances",
     &SU3BaryonSingletOctetScalarDecayer::_C, 0.39, -10.0, 10.0,
     false, false, true);

  static Switch<SU3BaryonSingletOctetScalarDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonSingletOctetScalarDecayer::_parity, true, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "Same parity",
     true);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "Opposite parity",
     false);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonSingletOctetScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the lighter proton-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_proton, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the lighter neutron-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_neutron, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the lighter Sigma+-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigmap, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the lighter Sigma--like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigmam, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the lighter Xi0-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_xi0, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the lighter Xi--like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_xim, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_elambda, 13122, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonSingletOctetScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonSingletOctetScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-0
void SU3BaryonSingletOctetScalarDecayer::halfHalfScalarCoupling(int imode,Energy m0,
								Energy m1,Energy m2,
								Complex& A,
								Complex& B) const
{
  if(_parity)
    {
      A=0.;
      B=_prefactor[imode]*(m0+m1);
    }
  else
    {
      A=_prefactor[imode]*(m0-m1);
      B=0.;
    }
}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonSingletOctetScalarDecayer::threeHalfHalfScalarCoupling(int imode,Energy m0,
								     Energy m1,Energy m2,
								     Complex& A,
								     Complex& B) const
{
  if(_parity)
    {
      A=_prefactor[imode]*(m0+m1);
      B=0.;
    }
  else
    {
      A=0.;
      B=_prefactor[imode]*(m0+m1);
    }
}

// set up the decay modes
void SU3BaryonSingletOctetScalarDecayer::setupModes(unsigned int iopt) const
{
  if(_outgoingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.resize(0);_outgoingM.resize(0);}
  // set up for the various different decay modes
  vector<int> outtemp,mestemp;
  double rt(sqrt(2.));
  if(_elambda==0){throw DecayIntegratorError() << "Invalid incoming particle in "
					       << "SU3BaryonSingletOctetScalarDecayer::" 
					       << "setupModes()" << Exception::abortnow;}
  // decays of the excited lambda
  outtemp.push_back(_sigma0);mestemp.push_back(111);
  outtemp.push_back(_sigmap);mestemp.push_back(-211);
  outtemp.push_back(_sigmam);mestemp.push_back(211);
  outtemp.push_back(_lambda);mestemp.push_back(221);
  outtemp.push_back(_xim);mestemp.push_back(321);
  outtemp.push_back(_xi0);mestemp.push_back(311);
  outtemp.push_back(_proton);mestemp.push_back(-321);
  outtemp.push_back(_neutron);mestemp.push_back(-311);
  PDVector extpart(3);
  extpart[0]=getParticleData(_elambda);
  int inspin(extpart[0]->iSpin()),outspin;
  for(unsigned int ix=0;ix<outtemp.size();++ix)
    {
      if(outtemp[ix]!=0&&mestemp[ix]!=0)
	{
	  extpart[1]=getParticleData(outtemp[ix]);
	  extpart[2]=getParticleData(mestemp[ix]);
	  if(extpart[0]->massMax()>extpart[1]->massMin()+extpart[2]->massMin())
	    {
	      _outgoingB.push_back(outtemp[ix]);
	      _outgoingM.push_back(mestemp[ix]);
	      if(iopt==1)
		{
		  outspin = extpart[1]->iSpin();
		  if(inspin==2&&outspin==2)
		    {_prefactor.push_back(_C*rt/_fpi);}
		  else if(inspin==4&&outspin==2)
		    {_prefactor.push_back(_C*rt/_fpi);}
		  else
		    {throw DecayIntegratorError() 
			<< "Invalid combination of spins in "
			<< "SU3BaryonSingletOctetScalarDecayer::" 
			<< "setupModes()" << Exception::abortnow;}
		}
	    }
	}
    }
}

void SU3BaryonSingletOctetScalarDecayer::dataBaseOutput(ofstream & output,
							bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "set " << fullName() << ":Coupling " << _C << "\n";
  output << "set " << fullName() << ":Parity " << _parity<< "\n";
  output << "set " << fullName() << ":Fpi " << _fpi << "\n";
  output << "set " << fullName() << ":Proton " << _proton << "\n";
  output << "set " << fullName() << ":Neutron " << _neutron << "\n";
  output << "set " << fullName() << ":Sigma+ " << _sigmap << "\n";
  output << "set " << fullName() << ":Sigma0 " << _sigma0 << "\n";
  output << "set " << fullName() << ":Sigma- " << _sigmam << "\n";
  output << "set " << fullName() << ":Lambda " << _lambda << "\n";
  output << "set " << fullName() << ":Xi0 " << _xi0 << "\n";
  output << "set " << fullName() << ":Xi- " << _xim << "\n"; 
  output << "set " << fullName() << ":ExcitedLambda " << _elambda << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix)
    {output << "insert " << fullName() << ":MaxWeight " << ix << " " 
	    << _maxweight[ix] << "\n";}
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

}
