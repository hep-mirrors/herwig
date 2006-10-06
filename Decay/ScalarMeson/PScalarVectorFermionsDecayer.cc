// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarVectorFermionsDecayer class.
//

#include "PScalarVectorFermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalarVectorFermionsDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "Herwig++/PDT/ThreeBodyAllOn1IntegralCalculator.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::ScalarWaveFunction;
using Helicity::VectorWaveFunction;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::EpsFunction;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using Helicity::incoming;
using Helicity::outgoing;

PScalarVectorFermionsDecayer::PScalarVectorFermionsDecayer() 
{
  // pi0 -> gamma e+e-
  _incoming.push_back(111);_outgoingV.push_back( 22);
  _outgoingf.push_back(11);_outgoinga.push_back(-11);
  _coupling.push_back(0.00761872/GeV);_maxweight.push_back(0.025);
  _includeVMD.push_back(2);_VMDid.push_back(113);
  _VMDmass.push_back(0.7758*GeV);_VMDwidth.push_back(0.1503*GeV);
  // eta -> gamma e+e-/mu+/mu-
  _incoming.push_back(221);_outgoingV.push_back( 22);
  _incoming.push_back(221);_outgoingV.push_back( 22);
  _outgoingf.push_back(11);_outgoinga.push_back(-11);
  _outgoingf.push_back(13);_outgoinga.push_back(-13);
  _coupling.push_back(0.007554164/GeV);_maxweight.push_back(2.5);
  _coupling.push_back(0.007554164/GeV);_maxweight.push_back(2.);
  _includeVMD.push_back(2);_VMDid.push_back(113);
  _includeVMD.push_back(2);_VMDid.push_back(113);
  _VMDmass.push_back(0.7758*GeV);_VMDwidth.push_back(0.1503*GeV);
  _VMDmass.push_back(0.7758*GeV);_VMDwidth.push_back(0.1503*GeV);
  // eta' -> gamma e+e-/mu+mu-
  _incoming.push_back(331);_outgoingV.push_back( 22);
  _incoming.push_back(331);_outgoingV.push_back( 22);
  _outgoingf.push_back(11);_outgoinga.push_back(-11);
  _outgoingf.push_back(13);_outgoinga.push_back(-13);
  _coupling.push_back(0.0104/GeV);_maxweight.push_back(4.5);
  _coupling.push_back(0.0104/GeV);_maxweight.push_back(3.0);
  _includeVMD.push_back(2);_VMDid.push_back(113);
  _VMDmass.push_back(0.7758*GeV);_VMDwidth.push_back(0.1503*GeV);
  // eta' -> gamma mu+mu-  
  _includeVMD.push_back(2);_VMDid.push_back(113);
  _VMDmass.push_back(0.7758*GeV);_VMDwidth.push_back(0.1503*GeV);
  // initial size of the arrays
  _initsize = _incoming.size();
  // intermediates
  generateIntermediates(false);
}

void PScalarVectorFermionsDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize=_coupling.size();
  if(isize!=_incoming.size()  || isize!=_outgoingV.size()|| isize!=_outgoingf.size()||
     isize!=_outgoinga.size() || isize!=_maxweight.size()|| isize!=_includeVMD.size()||
     isize!=_VMDid.size()     || isize!=_VMDmass.size()  || isize!=_VMDwidth.size())
    {throw InitException() << "Inconsistent parameters in PScalarVectorFermionsDecayer"
			   << Exception::abortnow;}
  // create the integration channel for each mode 
  PDVector extpart(4);
  tPDPtr gamma(getParticleData(ParticleID::gamma));
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode;
  vector<double> wgt(1,1.);
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      extpart[0] = getParticleData(_incoming[ix]);
      extpart[1] = getParticleData(_outgoingV[ix]);
      extpart[2] = getParticleData(_outgoingf[ix]);
      extpart[3] = getParticleData(_outgoinga[ix]);
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
      newchannel->addIntermediate(gamma     ,1,-1.1, 2,3);
      mode->addChannel(newchannel);
      addMode(mode,_maxweight[ix],wgt);
    }
  // set up the values for the VMD factor if needed (copy the default mass and width 
  //                                                 into the array)
  for(unsigned ix=0;ix<isize;++ix)
    {
      if(_includeVMD[ix]==1)
	{
	  _VMDmass[ix]=getParticleData(_VMDid[ix])->mass();
	  _VMDwidth[ix]=getParticleData(_VMDid[ix])->width();
	}
    }
}

PScalarVectorFermionsDecayer::~PScalarVectorFermionsDecayer() {}

int PScalarVectorFermionsDecayer::modeNumber(bool & cc,const DecayMode & dm) const
{
  int imode(-1);
  // must be three outgoing particles
  if(dm.products().size()!=3){return imode;}
  // ids of the particles
  int id0(dm.parent()->id()),idf[2],idv(0);
  unsigned int nf(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit)
    {
      if((**pit).iSpin()==PDT::Spin1){idv=(**pit).id();}
      else{idf[nf]=(**pit).id();++nf;}
    }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id0&&_outgoingV[ix]==idv)
	{if((idf[0]==_outgoingf[ix]&&idf[1]==_outgoinga[ix])||
	    (idf[1]==_outgoingf[ix]&&idf[0]==_outgoinga[ix])){imode=ix;}}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  cc=false;
  return imode;
}

void PScalarVectorFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoingV << _outgoingf << _outgoinga << _maxweight
     << _includeVMD << _VMDid << _VMDmass << _VMDwidth;
}

void PScalarVectorFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoingV >> _outgoingf >> _outgoinga >> _maxweight
     >> _includeVMD >> _VMDid >> _VMDmass >> _VMDwidth;
}

ClassDescription<PScalarVectorFermionsDecayer> PScalarVectorFermionsDecayer::initPScalarVectorFermionsDecayer;
// Definition of the static class description member.

void PScalarVectorFermionsDecayer::Init() {

  static ClassDocumentation<PScalarVectorFermionsDecayer> documentation
    ("The PScalarVectorFermionsDecayer class is designed"
     " for the decay of a pseudoscalar meson to a photon and a"
     "fermion-antifermion pair");

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalarVectorFermionsDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceOutcomingV
    ("OutgoingVector",
     "The PDG code for the outgoing pseudoscalar",
     &PScalarVectorFermionsDecayer::_outgoingV,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceOutcomingF
    ("OutgoingFermion",
     "The PDG code for the outgoing fermion",
     &PScalarVectorFermionsDecayer::_outgoingf,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceOutcomingA
    ("OutgoingAntiFermion",
     "The PDG code for the outgoing antifermion",
     &PScalarVectorFermionsDecayer::_outgoinga,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalarVectorFermionsDecayer::_coupling,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalarVectorFermionsDecayer::_maxweight,
     0, 0, 0, 0.0, 100., false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceIncludeVMD
    ("IncludeVMD",
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &PScalarVectorFermionsDecayer::_includeVMD,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceVMDID
    ("VMDID",
     "The PDG code for the particle to be used for the VMD factor.",
     &PScalarVectorFermionsDecayer::_VMDid,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,Energy> interfaceVMDmass
    ("VMDmass",
     "The mass to use for the particle in the VMD factor",
     &PScalarVectorFermionsDecayer::_VMDmass,
     0, 0, 0, 0., 10000., false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,Energy> interfaceVMDwidth
    ("VMDwidth",
     "The width to use for the particle in the VMD factor",
     &PScalarVectorFermionsDecayer::_VMDwidth,
     0, 0, 0, 0., 10000., false, false, true);

}

double PScalarVectorFermionsDecayer::me2(bool vertex, const int,
					 const Particle & inpart,
					 const ParticleVector & decay) const
{
  // workaround for gcc 3.2.3 bug
  //ALB ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);
  // vectors containing the spinors and polarization vectors
  vector<LorentzSpinor> wave;
  vector<LorentzSpinorBar> wavebar;
  vector<LorentzPolarizationVector> vwave;
  // set up the spin info for the outgoing particles
  VectorWaveFunction(vwave     ,decay[0],outgoing,true,true,vertex);
  SpinorBarWaveFunction(wavebar,decay[1],outgoing,true,vertex);
  SpinorWaveFunction(   wave   ,decay[2],outgoing,true,vertex);
  // now compute the matrix element
  Complex ii(0.,1.);
  Lorentz5Momentum pff(decay[1]->momentum()+decay[2]->momentum());
  pff.rescaleMass();
  Energy2 mff2(pff.mass()*pff.mass());
  // compute the prefactor
  Complex pre(_coupling[imode()]/mff2);
  // the VMD factor
  if(_includeVMD[imode()]>0)
    {
      Energy2 mrho2=_VMDmass[imode()]*_VMDmass[imode()];
      Energy2 mwrho=_VMDmass[imode()]*_VMDwidth[imode()];
      pre*= (-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
    }
  LorentzPolarizationVector eps,fcurrent;
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half);
  vector<unsigned int> ispin(4);ispin[0]=0;
  for(ispin[3]=0;ispin[3]<2;++ispin[3])
    {
      for(ispin[2]=0;ispin[2]<2;++ispin[2])
	{
	  fcurrent = wave[ispin[3]].vectorCurrent(wavebar[ispin[2]]);
	  // compute the current for this part
	  eps = EpsFunction::product(decay[0]->momentum(),pff,fcurrent);
	  for(ispin[1]=0;ispin[1]<3;++ispin[1])
	    {newME(ispin)=pre*(vwave[ispin[1]]*eps);}
	}	  
    }
  ME(newME);
  RhoDMatrix rhoin=RhoDMatrix(PDT::Spin0);rhoin.average();
  /*code to test the matrix element against the analytic result
  Energy  m[4]={inpart.mass(),decay[0]->mass(),decay[1]->mass(),decay[2]->mass()};
  Energy m2[4]={m[0]*m[0],m[1]*m[1],m[2]*m[2],m[3]*m[3]};
  Lorentz5Momentum p12=decay[0]->momentum()+decay[1]->momentum();p12.rescaleMass();
  Energy m122(p12.mass2());
  cout << "testing the matrix element " 
       <<(pre*conj(pre)).real()*(
				 -2*m122*m122*mff2 - mff2*mff2*mff2 + 
				 m2[1]*(2*m2[2]*m2[3] - 2*m2[3]*m2[3] + 
					m2[1]*(m2[2] - 2*m[2]*m[3] - m2[3])) - 
				 2*m[2]*(m[2]*m2[2] - 2*m2[1]*m[3] - m[2]*m2[3])*
				 m2[0] - (m2[2] + 2*m[2]*m[3] - m2[3])*
				 m2[0]*m2[0] +mff2*mff2*
				 (2*m2[1] + (m[2] - m[3])*(m[2] - m[3]) + 2*m2[0]) - 
				 mff2*(m2[1]*m2[1] + 2*m2[1]*m[2]*(m[2] - 2*m[3]) + 
				       2*m2[2]*m2[3] - 2*(2*m[2] - m[3])*m[3]*m2[0] + 
				       m2[0]*m2[0]) + 2*m122*
				 (-mff2*mff2 - (m2[2] - m2[3])*
				  (m2[1] - m2[0]) + 
				  mff2*(m2[1] + m2[2] + m2[3] + 
					m2[0]))) + me << endl;
  */
  return newME.contract(rhoin).real();
}

// method to return an object to calculate the 3 or higher body partial width
WidthCalculatorBasePtr 
PScalarVectorFermionsDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // workout which mode we are doing
  int imode(-1);
  // ids of the particles
  int id0(dm.parent()->id()),idf[2],idv(0);
  unsigned int nf(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit) {
    if((**pit).iSpin()==PDT::Spin1){idv=(**pit).id();}
    else{idf[nf]=(**pit).id();++nf;}
  }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do {
    if(_incoming[ix]==id0&&_outgoingV[ix]==idv) {
      if((idf[0]==_outgoingf[ix]&&idf[1]==_outgoinga[ix])||
	 (idf[1]==_outgoingf[ix]&&idf[0]==_outgoinga[ix])) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  // get the masses we need
  Energy m[3]={getParticleData(_outgoingV[imode])->mass(),
	       getParticleData(_outgoingf[imode])->mass(),
	       getParticleData(_outgoinga[imode])->mass()};
  return 
    new_ptr(ThreeBodyAllOn1IntegralCalculator<PScalarVectorFermionsDecayer>
	    (3,-1000.,-0.9,*this,imode,m[0],m[1],m[2]));
}


double PScalarVectorFermionsDecayer::threeBodydGammads(const int imodeb,
						       const Energy2 q2, const 
						       Energy2 mff2, const  Energy m1,
						       const Energy m2,
						       const  Energy m3) const {
  // the masses of the external particles
  Energy q=sqrt(q2);
  Energy2 m12=m1*m1;
  Energy2 m22=m2*m2;
  Energy2 m32=m3*m3;
  // calculate the prefactor
  Complex pre=_coupling[imodeb],ii(0.,1.);
  pre /= mff2;
  // the VMD factor
  if(_includeVMD[imodeb]>0) {
    Energy2 mrho2=_VMDmass[imodeb]*_VMDmass[imodeb];
    Energy2 mwrho=_VMDmass[imodeb]*_VMDwidth[imodeb];
    pre*= (-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
  }
  double factor=real(pre*conj(pre));
  // compute the pieces from the integration limits
  Energy mff=sqrt(mff2);
  Energy e2star = 0.5*(mff2-m32+m22)/mff;
  Energy e1star = 0.5*(q2-mff2-m12)/mff;
  Energy e1sm = sqrt(e1star*e1star-m12);
  Energy e2sm = sqrt(e2star*e2star-m22);
  Energy2 a = 2*e1star*e2star+m12+m22;
  Energy2 b = 2*e1sm*e2sm;
  // term independent of s3
  double me = 2*b*(2*(m12*(mff2*mff2 + 4*mff2*m2*m3 -(m22 - m32)*(m22 - m32)) + 
		      2*m2*(m12 +m22)*m3*(-mff2 +m22 + q2))
		   +(m12 +m22)*(m12 +m22)*(-mff2 +m22 - 2*m2*m3 - m32)
		   -(mff2 +m22 + 2*m2*m3 - m32)*(-mff2 +m22 + q2)*(-mff2 +m22 + q2));
  // linear term
  me+= 2.*a*b*(2*(-mff2*mff2 - (m22 - m32)*(m12 - q2) + 
		  mff2*(m12 + m22 + m32 + q2)));
  // quadratic term
  me+=-4.*mff2*b*(3.*a*a+b*b)/3.;
  me*=-factor;
  // phase space factors
  return me/256./pi/pi/pi/q2/q;
}

// output the setup information for the particle database
void PScalarVectorFermionsDecayer::dataBaseOutput(ofstream & output,
						  bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":Incoming   " << ix << "  " 
		 << _incoming[ix]   << "\n";
	  output << "set " << fullName() << ":OutgoingVector  " 
		 << ix << "  " << _outgoingV[ix]  << "\n";
	  output << "set " << fullName() << ":OutgoingFermion  " 
		 << ix << "  " << _outgoingf[ix]  << "\n";
	  output << "set " << fullName() << ":OutgoingAntiFermion " 
		 << ix << "  " << _outgoinga[ix]  << "\n";
	  output << "set " << fullName() << ":Coupling   " << ix << "  " 
		 << _coupling[ix]   << "\n";
	  output << "set " << fullName() << ":MaxWeight  " << ix << "  " 
		 << _maxweight[ix]  << "\n";
	  output << "set " << fullName() << ":IncludeVMD " << ix << "  " 
		 << _includeVMD[ix] << "\n";
	  output << "set " << fullName() << ":VMDID      " << ix << "  " 
		 << _VMDid[ix]      << "\n";
	  output << "set " << fullName() << ":VMDmass    " << ix << "  " 
		 << _VMDmass[ix]    << "\n";
	  output << "set " << fullName() << ":VMDwidth   " << ix << "  " 
		 << _VMDwidth[ix]   << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":Incoming   " << ix << "  " 
		 << _incoming[ix]   << "\n";
	  output << "insert " << fullName() << ":OutgoingVector  " 
		 << ix << "  " << _outgoingV[ix]  << "\n";
	  output << "insert " << fullName() << ":OutgoingFermion  " 
		 << ix << "  " << _outgoingf[ix]  << "\n";
	  output << "insert " << fullName() << ":OutgoingAntiFermion " 
		 << ix << "  " << _outgoinga[ix]  << "\n";
	  output << "insert " << fullName() << ":Coupling   " << ix << "  " 
		 << _coupling[ix]   << "\n";
	  output << "insert " << fullName() << ":IncludeVMD " << ix << "  " 
		 << _includeVMD[ix] << "\n";
	  output << "insert " << fullName() << ":VMDID      " << ix << "  " 
		 << _VMDid[ix]      << "\n";
	  output << "insert " << fullName() << ":VMDmass    " << ix << "  " 
		 << _VMDmass[ix]    << "\n";
	  output << "insert " << fullName() << ":VMDwidth   " << ix << "  " 
		 << _VMDwidth[ix]   << "\n";
	}
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}


