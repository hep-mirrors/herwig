// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonDecupletOctetScalarDecayer class.
//

#include "SU3BaryonDecupletOctetScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SU3BaryonDecupletOctetScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig{
using namespace ThePEG;

SU3BaryonDecupletOctetScalarDecayer::~SU3BaryonDecupletOctetScalarDecayer() {}
int SU3BaryonDecupletOctetScalarDecayer::modeNumber(bool & cc,const DecayMode & dm) const
{
  int imode(-1);
  if(_incomingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(dm.products().size()!=2){return imode;}
  // ids of the particles
  int id0(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);
  cc=false;
  do
    {
      if(id0==_incomingB[ix])
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){imode=ix;cc=false;}
	}
      else if(id0==-_incomingB[ix])
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
  return imode;
}

void SU3BaryonDecupletOctetScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _C << _parity << _fpi << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _deltapp << _deltap << _delta0 << _deltam
     << _sigmasp << _sigmas0 << _sigmasm << _omega << _xism << _xis0 << _incomingB 
     << _outgoingB << _outgoingM << _maxweight << _prefactor;
}

void SU3BaryonDecupletOctetScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _C >> _parity >> _fpi >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _deltapp >> _deltap >> _delta0 >> _deltam
     >> _sigmasp >> _sigmas0 >> _sigmasm >> _omega >> _xism >> _xis0 >> _incomingB 
     >> _outgoingB >> _outgoingM >> _maxweight >> _prefactor;
}

ClassDescription<SU3BaryonDecupletOctetScalarDecayer> SU3BaryonDecupletOctetScalarDecayer::initSU3BaryonDecupletOctetScalarDecayer;
// Definition of the static class description member.

void SU3BaryonDecupletOctetScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonDecupletOctetScalarDecayer> documentation
    ("The SU3BaryonDecupletOctetScalarDecayer class is designed for the"
     " decay of an SU(3) decuplet baryon to an SU(3) octet baryon and a pseudoscalar"
     " meson from the lightest multiplet.");

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,double> interfaceCcoupling
    ("Ccoupling",
     "The C coupling for the decuplet decays.",
     &SU3BaryonDecupletOctetScalarDecayer::_C, 1.5, -10.0, 10.0,
     false, false, true);

  static Switch<SU3BaryonDecupletOctetScalarDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonDecupletOctetScalarDecayer::_parity, true, false, false);
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

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonDecupletOctetScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the proton-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_proton, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the neutron-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_neutron, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the Sigma+-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmap, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the Sigma0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the Sigma--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmam, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the Lambda-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the Xi0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xi0, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the Xi--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xim, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceDeltapp
    ("Delta++",
     "The PDG code for the Delta++ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltapp, 2224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceDeltap
    ("Delta+",
     "The PDG code for the Delta+ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltap, 2214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceDelta0
    ("Delta0",
     "The PDG code for the Delta0 like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_delta0, 2114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceDeltam
    ("Delta-",
     "The PDG code for the Delta- like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltam, 1114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigmasp
    ("Sigma*+",
     "The PDG code for the Sigma*+ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmasp, 3224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigmas0
    ("Sigma*0",
     "The PDG code for the Sigma*0 like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmas0, 3214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigmasm
    ("Sigma*-",
     "The PDG code for the Sigma*- like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmasm, 3114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceOmega
    ("Omega",
     "The PDG code for the Omega like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_omega, 3334, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceXis0
    ("Xi*0",
     "The PDG code for the Xi*0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xis0, 3324, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceXism
    ("Xi*-",
     "The PDG code for the Xi*--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xism, 3314, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonDecupletOctetScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonDecupletOctetScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonDecupletOctetScalarDecayer::
threeHalfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy,
			    Complex& A, Complex& B) const
{
  if(_parity){A=_prefactor[imode]*(m0+m1);B=0.;}
  else{A=0.;B=_prefactor[imode]*(m0+m1);}
}

// set up the decay modes
void SU3BaryonDecupletOctetScalarDecayer::setupModes(unsigned int iopt) const
{
  if(_incomingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.resize(0);_incomingB.resize(0);_outgoingM.resize(0);}
  vector<double> factor;
  vector<int> intemp,outtemp,mestemp;
  double ort(1./sqrt(2.)),ors(1./sqrt(6.)),rt(sqrt(2.)),orr(1./sqrt(3.));
  // decays of the delta++
  intemp.push_back(_deltapp);outtemp.push_back(_proton);mestemp.push_back(211);
  factor.push_back(_C);
  intemp.push_back(_deltapp);outtemp.push_back(_sigmap);mestemp.push_back(321);
  factor.push_back(-_C);
  // decays of the delta+
  intemp.push_back(_deltap);outtemp.push_back(_neutron);mestemp.push_back(211);
  factor.push_back(_C*orr);
  intemp.push_back(_deltap);outtemp.push_back(_proton);mestemp.push_back(111);
  factor.push_back(_C*rt*orr);
  intemp.push_back(_deltap);outtemp.push_back(_sigma0);mestemp.push_back(321);
  factor.push_back(_C*rt*orr);
  intemp.push_back(_deltap);outtemp.push_back(_sigmap);mestemp.push_back(311);
  factor.push_back(_C*orr);
  // decays of the delta0
  intemp.push_back(_delta0);outtemp.push_back(_proton);mestemp.push_back(-211);
  factor.push_back(_C*orr);
  intemp.push_back(_delta0);outtemp.push_back(_neutron);mestemp.push_back(111);
  factor.push_back(_C*rt*orr);
  intemp.push_back(_delta0);outtemp.push_back(_sigma0);mestemp.push_back(311);
  factor.push_back(_C*rt*orr);
  intemp.push_back(_delta0);outtemp.push_back(_sigmam);mestemp.push_back(321);
  factor.push_back(_C*orr);
  // decays of the delta-
  intemp.push_back(_deltam);outtemp.push_back(_neutron);mestemp.push_back(-211);
  factor.push_back(-_C);
  intemp.push_back(_deltam);outtemp.push_back(_sigmam);mestemp.push_back(311);
  factor.push_back(_C);
  // sigma*+
  intemp.push_back(_sigmasp);outtemp.push_back(_sigmap);mestemp.push_back(111);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmasp);outtemp.push_back(_sigma0);mestemp.push_back(211);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmasp);outtemp.push_back(_proton);mestemp.push_back(-311);
  factor.push_back(_C*orr);
  intemp.push_back(_sigmasp);outtemp.push_back(_xi0);mestemp.push_back(321);
  factor.push_back(_C*orr);
  intemp.push_back(_sigmasp);outtemp.push_back(_sigmap);mestemp.push_back(221);
  factor.push_back(_C*ort);
  intemp.push_back(_sigmasp);outtemp.push_back(_lambda);mestemp.push_back(211);
  factor.push_back(_C*ort);
  // sigma*0
  intemp.push_back(_sigmas0);outtemp.push_back(_neutron);mestemp.push_back(-311);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_proton);mestemp.push_back(-321);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_xim);mestemp.push_back(321);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_xi0);mestemp.push_back(311);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigmam);mestemp.push_back(211);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigmap);mestemp.push_back(-211);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_lambda);mestemp.push_back(111);
  factor.push_back(_C*ort);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigma0);mestemp.push_back(211);
  factor.push_back(_C*ort);
  // sigma*-
  intemp.push_back(_sigmasm);outtemp.push_back(_sigmam);mestemp.push_back(111);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmasm);outtemp.push_back(_sigma0);mestemp.push_back(-211);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmasm);outtemp.push_back(_neutron);mestemp.push_back(-321);
  factor.push_back(_C*orr);
  intemp.push_back(_sigmasm);outtemp.push_back(_xim);mestemp.push_back(311);
  factor.push_back(_C*orr);
  intemp.push_back(_sigmasm);outtemp.push_back(_lambda);mestemp.push_back(-211);
  factor.push_back(_C*ort);
  intemp.push_back(_sigmasm);outtemp.push_back(_sigmam);mestemp.push_back(221);
  factor.push_back(_C*ort);
  // xi*0
  intemp.push_back(_xis0);outtemp.push_back(_xim);mestemp.push_back(211);
  factor.push_back(_C*orr);
  intemp.push_back(_xis0);outtemp.push_back(_xi0);mestemp.push_back(111);
  factor.push_back(_C*ors);
  intemp.push_back(_xis0);outtemp.push_back(_sigmap);mestemp.push_back(-321);
  factor.push_back(_C*orr);
  intemp.push_back(_xis0);outtemp.push_back(_sigma0);mestemp.push_back(-311);
  factor.push_back(_C*ors);
  intemp.push_back(_xis0);outtemp.push_back(_xi0);mestemp.push_back(221);
  factor.push_back(_C*ort);
  intemp.push_back(_xis0);outtemp.push_back(_lambda);mestemp.push_back(-311);
  factor.push_back(_C*ort);
  // xi*-
  intemp.push_back(_xism);outtemp.push_back(_xi0);mestemp.push_back(-211);
  factor.push_back(_C*orr);
  intemp.push_back(_xism);outtemp.push_back(_xim);mestemp.push_back(111);
  factor.push_back(_C*ors);
  intemp.push_back(_xism);outtemp.push_back(_sigmam);mestemp.push_back(-311);
  factor.push_back(_C*orr);
  intemp.push_back(_xism);outtemp.push_back(_sigma0);mestemp.push_back(-321);
  factor.push_back(_C*ors);
  intemp.push_back(_xism);outtemp.push_back(_xim);mestemp.push_back(221);
  factor.push_back(_C*ort);
  intemp.push_back(_xism);outtemp.push_back(_lambda);mestemp.push_back(-321);
  factor.push_back(_C*ort);
  // set up the modes
  PDVector extpart(3);
  for(unsigned int ix=0;ix<intemp.size();++ix)
    {
      if(intemp[ix]!=0&&outtemp[ix]!=0&&mestemp[ix]!=0)
	{
	  extpart[0]=getParticleData(intemp[ix]);
	  extpart[1]=getParticleData(outtemp[ix]);
	  extpart[2]=getParticleData(mestemp[ix]);
	  if(extpart[0]->massMax()>extpart[1]->massMin()+extpart[2]->massMin())
	    {
	      _incomingB.push_back(intemp[ix]);
	      _outgoingB.push_back(outtemp[ix]);
	      _outgoingM.push_back(mestemp[ix]);
	      if(iopt==1)
		{_prefactor.push_back(factor[ix]/_fpi);}
	    }
	}
    }
}

void SU3BaryonDecupletOctetScalarDecayer::dataBaseOutput(ofstream & output, 
							 bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "set " << fullName() << ":Ccoupling " << _C<< "\n";
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
  output << "set " << fullName() << ":Delta++ " << _deltapp << "\n";
  output << "set " << fullName() << ":Delta+ " << _deltap << "\n";
  output << "set " << fullName() << ":Delta0 " << _delta0 << "\n";
  output << "set " << fullName() << ":Delta- " << _deltam << "\n";
  output << "set " << fullName() << ":Sigma*+ " << _sigmasp << "\n";
  output << "set " << fullName() << ":Sigma*0 " << _sigmas0 << "\n";
  output << "set " << fullName() << ":Sigma*- " << _sigmasm << "\n";
  output << "set " << fullName() << ":Omega " << _omega << "\n";
  output << "set " << fullName() << ":Xi*0 " << _xis0 << "\n";
  output << "set " << fullName() << ":Xi*- " << _xism << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix)
    {output << "insert " << fullName() << ":MaxWeight " << ix << " " 
	    << _maxweight[ix] << "\n";}
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

}
