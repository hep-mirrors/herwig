// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonPScalarFermionsDecayer class.
//
//  Author: Peter Richardson
//

#include "VectorMesonPScalarFermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonPScalarFermionsDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "Herwig++/PDT/ThreeBodyAllOn1IntegralCalculator.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::ScalarWaveFunction;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::EpsFunction;
using Helicity::outgoing;

VectorMesonPScalarFermionsDecayer::VectorMesonPScalarFermionsDecayer() 
{
  // omega -> pi e+e- /mu+mu-
  _incoming.push_back( 223);_outgoingP.push_back( 111);
  _incoming.push_back( 223);_outgoingP.push_back( 111);
  _outgoingf.push_back(11);_outgoinga.push_back(-11);
  _outgoingf.push_back(13);_outgoinga.push_back(-13);
  _coupling.push_back(0.2096/GeV);_maxweight.push_back(5.);_weight.push_back(0.);
  _coupling.push_back(0.2096/GeV);_maxweight.push_back(3.4);_weight.push_back(0.);
  _includeVMD.push_back(2);_VMDid.push_back(113);
  _includeVMD.push_back(2);_VMDid.push_back(113);
  _VMDmass.push_back(0.7758*GeV);_VMDwidth.push_back(0.1503*GeV);
  _VMDmass.push_back(0.7758*GeV);_VMDwidth.push_back(0.1503*GeV);
  // phi -> eta e+e-/mu+mu-
  _incoming.push_back( 333);_outgoingP.push_back( 221);
  _incoming.push_back( 333);;_outgoingP.push_back( 221);
  _outgoingf.push_back(11);_outgoinga.push_back(-11);
  _outgoingf.push_back(13);_outgoinga.push_back(-13);
  _coupling.push_back(0.0643/GeV);_maxweight.push_back(4.7);_weight.push_back(0.);
  _coupling.push_back(0.0643/GeV);_maxweight.push_back(3.5);_weight.push_back(0.);
  _includeVMD.push_back(2);_VMDid.push_back(113);
  _includeVMD.push_back(2);_VMDid.push_back(113);
  _VMDmass.push_back(0.7758*GeV);_VMDwidth.push_back(0.1503*GeV);
  _VMDmass.push_back(0.7758*GeV);_VMDwidth.push_back(0.1503*GeV);
  // phi -> pi e+e-/mu+mu-
  _incoming.push_back( 333);;_outgoingP.push_back( 111);
  _incoming.push_back( 333);;_outgoingP.push_back( 111);
  _outgoingf.push_back(11);_outgoinga.push_back(-11);
  _outgoingf.push_back(13);_outgoinga.push_back(-13);
  _coupling.push_back(0.0120094/GeV);_maxweight.push_back(6.0);_weight.push_back(0.10);
  _coupling.push_back(0.0120094/GeV);_maxweight.push_back(4.0);_weight.push_back(0.25);
  _includeVMD.push_back(2);_VMDid.push_back(113);
  _includeVMD.push_back(2);_VMDid.push_back(113);
  _VMDmass.push_back(0.7758*GeV);_VMDwidth.push_back(0.1503*GeV);
  _VMDmass.push_back(0.7758*GeV);_VMDwidth.push_back(0.1503*GeV);
  // the initial size of the arrays
  _initsize=_coupling.size();
}

void VectorMesonPScalarFermionsDecayer::doinit() throw(InitException) {
  VectorMesonDecayerBase::doinit();
  // check the parameters are consistent
  unsigned int isize(_coupling.size());
  if(isize!=_incoming.size()  || isize!=_outgoingP.size()|| isize!=_outgoingf.size()||
     isize!=_outgoinga.size() || isize!=_maxweight.size()|| isize!=_includeVMD.size()||
     isize!=_VMDid.size()     || isize!=_VMDmass.size()  || isize!=_VMDwidth.size()||
     isize!=_weight.size())
    {throw InitException() << "Inconsistent parameters in VectorMesonPScalar"
			   << "FermionsDecayer" << Exception::abortnow;}
  // create the integration channel for each mode 
  PDVector extpart(4);
  tPDPtr gamma(getParticleData(ParticleID::gamma)),rho;
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr newmode;
  vector<double> wgt(2,1.);
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      rho=getParticleData(_VMDid[ix]);
      extpart[0] = getParticleData(_incoming[ix]);
      extpart[1] = getParticleData(_outgoingP[ix]);
      extpart[2] = getParticleData(_outgoingf[ix]);
      extpart[3] = getParticleData(_outgoinga[ix]);
      newmode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      // photon channel
      newchannel=new_ptr(DecayPhaseSpaceChannel(newmode));
      newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
      newchannel->addIntermediate(gamma     ,1,-1.1, 2,3);
      newmode->addChannel(newchannel);
      wgt[0]=1.-_weight[ix];
      // vmd channel
      newchannel=new_ptr(DecayPhaseSpaceChannel(newmode));
      newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
      newchannel->addIntermediate(rho       ,0, 0.0, 2,3);
      newmode->addChannel(newchannel);
      wgt[1]=_weight[ix];
      addMode(newmode,_maxweight[ix],wgt);
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

VectorMesonPScalarFermionsDecayer::~VectorMesonPScalarFermionsDecayer() {}

bool VectorMesonPScalarFermionsDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed(false);
  // must be three outgoing particles
  if(dm.products().size()!=3){return allowed;}
  // ids of the particles
  int id0(dm.parent()->id()),idf[2],ids;
  unsigned int nf(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit)
    {
      if((**pit).iSpin()==PDT::Spin0){ids=(**pit).id();}
      else{idf[nf]=(**pit).id();++nf;}
    }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id0&&_outgoingP[ix]==ids)
	{if((idf[0]==_outgoingf[ix]&&idf[1]==_outgoinga[ix])||
	    (idf[1]==_outgoingf[ix]&&idf[0]==_outgoinga[ix])){allowed=true;}}
      ++ix;
    }
  while(!allowed&&ix<_incoming.size());
  return allowed;
}

ParticleVector VectorMesonPScalarFermionsDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  // workout which mode we are doing
  int imode(-1);
  // ids of the particles
  int id0(dm.parent()->id()),idf[2],ids;
  unsigned int nf(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit)
    {
      if((**pit).iSpin()==PDT::Spin0){ids=(**pit).id();}
      else{idf[nf]=(**pit).id();++nf;}
    }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id0&&_outgoingP[ix]==ids)
	{if((idf[0]==_outgoingf[ix]&&idf[1]==_outgoinga[ix])||
	    (idf[1]==_outgoingf[ix]&&idf[0]==_outgoinga[ix])){imode=ix;}}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // perform the decay
  bool cc(false);
  return generate(false,cc,imode,parent); 
}


void VectorMesonPScalarFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoingP << _outgoingf << _outgoinga << _maxweight
     << _weight << _includeVMD << _VMDid << _VMDmass << _VMDwidth;
}

void VectorMesonPScalarFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoingP >> _outgoingf >> _outgoinga >> _maxweight
     >> _weight >> _includeVMD >> _VMDid >> _VMDmass >> _VMDwidth;
}

ClassDescription<VectorMesonPScalarFermionsDecayer> VectorMesonPScalarFermionsDecayer::initVectorMesonPScalarFermionsDecayer;
// Definition of the static class description member.

void VectorMesonPScalarFermionsDecayer::Init() {

  static ClassDocumentation<VectorMesonPScalarFermionsDecayer> documentation
    ("The \\classname{VectorMesonPScalarFermionsDecayer} class is designed to "
     "perform the decay of a vector meson to a pseudoscalar meson and a "
     "fermion-antifermion pair.");

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonPScalarFermionsDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceOutcomingP
    ("OutgoingPseudoScalar",
     "The PDG code for the outgoing pseudoscalar",
     &VectorMesonPScalarFermionsDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceOutcomingF
    ("OutgoingFermion",
     "The PDG code for the outgoing fermion",
     &VectorMesonPScalarFermionsDecayer::_outgoingf,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceOutcomingA
    ("OutgoingAntiFermion",
     "The PDG code for the outgoing antifermion",
     &VectorMesonPScalarFermionsDecayer::_outgoinga,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonPScalarFermionsDecayer::_coupling,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonPScalarFermionsDecayer::_maxweight,
     0, 0, 0, 0., 200.0, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,double> interfaceWeight
    ("Weight",
     "The weight for vector meson phasse space channel",
     &VectorMesonPScalarFermionsDecayer::_weight,
     0, 0, 0, 0., 200.0, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceIncludeVMD
    ("IncludeVMD",
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &VectorMesonPScalarFermionsDecayer::_includeVMD,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceVMDID
    ("VMDID",
     "The PDG code for the particle to be used for the VMD factor.",
     &VectorMesonPScalarFermionsDecayer::_VMDid,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,double> interfaceVMDmass
    ("VMDmass",
     "The mass to use for the particle in the VMD factor",
     &VectorMesonPScalarFermionsDecayer::_VMDmass,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,double> interfaceVMDwidth
    ("VMDwidth",
     "The width to use for the particle in the VMD factor",
     &VectorMesonPScalarFermionsDecayer::_VMDwidth,
     0, 0, 0, -10000000, 10000000, false, false, true);

}

// the hadronic currents    
vector<LorentzPolarizationVector>  
VectorMesonPScalarFermionsDecayer::decayCurrent(const bool vertex,
						const int, const Particle & inpart,
						const ParticleVector & decay) const
{
  // construct the spin information objects for the  decay products
  vector<LorentzSpinor> wave;
  vector<LorentzSpinorBar> wavebar;
  ScalarWaveFunction(decay[0],outgoing,true,vertex);
  SpinorWaveFunction(   wave   ,decay[2],outgoing,true,vertex);
  SpinorBarWaveFunction(wavebar,decay[1],outgoing,true,vertex);
  // now compute the currents
  vector<LorentzPolarizationVector> temp(4);
  Lorentz5Momentum pff(decay[1]->momentum()+decay[2]->momentum());
  pff.rescaleMass();
  // the factor for the off-shell photon
  Energy2 mff2(pff.mass2());
  // prefactor
  Complex pre(_coupling[imode()]/mff2),ii(0.,1.);
  // the VMD factor
  if(_includeVMD[imode()]>0)
    {
      Energy2 mrho2(_VMDmass[imode()]*_VMDmass[imode()]);
      Energy2 mwrho(_VMDmass[imode()]*_VMDwidth[imode()]);
      pre*= (-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
    }
  unsigned int ix,iy;
  for(ix=0;ix<2;++ix)
    {for(iy=0;iy<2;++iy)
	{temp[2*iy+ix]=pre*EpsFunction::product(inpart.momentum(),pff,
						 wave[ix].vectorCurrent(wavebar[iy]));}}
  /* code for the spin averaged me for testing only
  Energy  m[4]={inpart.mass(),decay[0]->mass(),decay[1]->mass(),decay[2]->mass()};
  Energy m2[4]={m[0]*m[0],m[1]*m[1],m[2]*m[2],m[3]*m[3]};
  Lorentz5Momentum p12=decay[0]->momentum()+decay[1]->momentum();p12.rescaleMass();
  Energy m122(p12.mass2());
  cout << "testing the matrix element " 
       << -1./3.*(pre*conj(pre)).real()*(-2*m122*m122*mff2 - mff2*mff2*mff2 + 
		   m2[1]*(2*m2[2]*m2[3] - 2*m2[3]*m2[3] + 
			  m2[1]*(m2[2] - 2*m[2]*m[3] - m2[3])) - 
		   2*m[2]*(m2[2]*m[2] - 2*m2[1]*m[3] - m[2]*m2[3])*
		   m2[0] - (m2[2] + 2*m[2]*m[3] - m2[3])*
		   m2[0]*m2[0] + mff2*mff2*
		   (2*m2[1] + (m[2] - m[3])*(m[2] - m[3]) + 2*m2[0]) - 
		   mff2*(m2[1]*m2[1] + 2*m2[1]*m[2]*(m[2] - 2*m[3]) + 
			 2*m2[2]*m2[3] - 2*(2*m[2] - m[3])*m[3]*m2[0] + 
			 m2[0]*m2[0]) + 2*m122*
		   (-mff2*mff2 - (m[2] - m[3])*(m[2] + m[3])*(m[1] - m[0])*
		    (m[1] + m[0]) + mff2*
		    (m2[1] + m2[2] + m2[3] + m2[0])))
       << endl;
   */
  // return the answer
  return temp;
 }

WidthCalculatorBasePtr 
VectorMesonPScalarFermionsDecayer::threeBodyMEIntegrator(const DecayMode & dm) const
{
  // workout which mode we are doing
  int imode(-1);
  // ids of the particles
  int id0(dm.parent()->id()),idf[2],ids;
  unsigned int nf(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit)
    {
      if((**pit).iSpin()==PDT::Spin0){ids=(**pit).id();}
      else{idf[nf]=(**pit).id();++nf;}
    }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id0&&_outgoingP[ix]==ids)
	{if((idf[0]==_outgoingf[ix]&&idf[1]==_outgoinga[ix])||
	    (idf[1]==_outgoingf[ix]&&idf[0]==_outgoinga[ix])){imode=ix;}}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // get the masses we need
  Energy m[3]={getParticleData(_outgoingP[imode])->mass(),
	       getParticleData(_outgoingf[imode])->mass(),
	       getParticleData(_outgoinga[imode])->mass()};
  return 
    new_ptr(ThreeBodyAllOn1IntegralCalculator(3,-1000.,-0.8,
					      const_ptr_cast<tDecayIntegratorPtr>(this),
					      imode,m[0],m[1],m[2]));
} 

double VectorMesonPScalarFermionsDecayer::threeBodydGammads(int imodeb,Energy2 q2,
							    Energy2 mff2, Energy m1, 
							    Energy m2, Energy m3)
{
  // the masses of the external particles
  Energy q(sqrt(q2)),m12(m1*m1),m22(m2*m2),m32(m3*m3);
  // prefactor
  Complex pre(_coupling[imodeb]/mff2),ii(0.,1.);
  // the VMD factor
  if(_includeVMD[imodeb]>0)
    {
      Energy2 mrho2(_VMDmass[imodeb]*_VMDmass[imodeb]),
	mwrho(_VMDmass[imodeb]*_VMDwidth[imodeb]);
      pre*= (-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
    }
  double factor(real(pre*conj(pre)));
  // compute the pieces from the integration limits
  Energy mff(sqrt(mff2)),e2star(0.5*(mff2-m32+m22)/mff),e1star(0.5*(q2-mff2-m12)/mff),
    e1sm(sqrt(e1star*e1star-m12)),e2sm(sqrt(e2star*e2star-m22));
  Energy2 a(2*e1star*e2star+m12+m22),b(2*e1sm*e2sm);
  // term independent of s3
  double me = 2*b*(-mff2*mff2*mff2 +m12*m12*(m22 - 2*m2*m3 - m32) - 
		   2*m22*(m22 - m32)*q2 -(m22 + 2*m2*m3 - m32)*q2*q2 + 
		   mff2*mff2*(2*m12 +(m2-m3)*(m2-m3)+2*q2) + 2*m12*m3*
		   ((m22-m32)*m3 + 2*m2*q2) - 
		   mff2*(m12*m12 + 2*m12*m2*(m2 - 2*m3) + 2*m22*m32 - 
			 2*(2*m2 - m3)*m3*q2 + q2*q2));
  // linear term
  me+= 2.*a*b*(2*(-mff2*mff2-(m22-m32)*(m12-q2)+mff2*(m12+m22+m32+q2)));
  // quadratic term
  me+=2*b*(3.*a*a+b*b)/3.*(-2*mff2);
  me*=-factor;
  // phase space factors
  return me/768./pi/pi/pi/q2/q;
}

// output the setup information for the particle database
void VectorMesonPScalarFermionsDecayer::dataBaseOutput(ofstream & output)
{
  output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  output << "set " << fullName() << ":Iteration " << _niter << "\n";
  output << "set " << fullName() << ":Ntry " << _ntry << "\n";
  output << "set " << fullName() << ":Points " << _npoint << "\n";
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":Incoming   " << ix << "  " 
		 << _incoming[ix]   << "\n";
	  output << "set " << fullName() << ":OutgoingPseudoScalar  " 
		 << ix << "  " << _outgoingP[ix]  << "\n";
	  output << "set " << fullName() << ":OutgoingFermion  " 
		 << ix << "  " << _outgoingf[ix]  << "\n";
	  output << "set " << fullName() << ":OutgoingAntiFermion " 
		 << ix << "  " << _outgoinga[ix]  << "\n";
	  output << "set " << fullName() << ":Coupling   " << ix << "  " 
		 << _coupling[ix]   << "\n";
	  output << "set " << fullName() << ":MaxWeight  " << ix << "  " 
		 << _maxweight[ix]  << "\n";
	  output << "set " << fullName() << ":Weight  "    << ix << "  " 
		 << _weight[ix]  << "\n";
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
	  output << "insert " << fullName() << ":OutgoingPseudoScalar  " 
		 << ix << "  " << _outgoingP[ix]  << "\n";
	  output << "insert " << fullName() << ":OutgoingFermion  " 
		 << ix << "  " << _outgoingf[ix]  << "\n";
	  output << "insert " << fullName() << ":OutgoingAntiFermion " 
		 << ix << "  " << _outgoinga[ix]  << "\n";
	  output << "insert " << fullName() << ":Coupling   " << ix << "  " 
		 << _coupling[ix]   << "\n";
	  output << "insert " << fullName() << ":Weight  "    << ix << "  " 
		 << _weight[ix]  << "\n";
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
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
}
