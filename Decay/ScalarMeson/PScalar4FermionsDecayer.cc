// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalar4FermionsDecayer class.
//

#include "PScalar4FermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/epsilon.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalar4FermionsDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Helicity::ScalarWaveFunction;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using Helicity::incoming;
using Helicity::outgoing;

void PScalar4FermionsDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize=_coupling.size();
  if(isize!=_incoming.size()  || isize!=_outgoing1.size() || isize!=_outgoing2.size()||
     isize!=_maxweight.size() || isize!=_includeVMD.size()|| isize!=_VMDid.size()    ||
     isize!=_VMDmass.size()  || isize!=_VMDwidth.size())
    {throw InitException() << "Inconsistent parameters in PScalar4FermionsDecayer"
			   << Exception::abortnow;}
  // create the integration channels for each mode 
  PDVector extpart(5);
  tPDPtr gamma=getParticleData(ParticleID::gamma);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode;
  vector<double> wgt;
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      wgt.resize(1);wgt[0]=1.;
      extpart[0] = getParticleData(_incoming[ix]);
      extpart[1] = getParticleData( _outgoing1[ix]);
      extpart[2] = getParticleData(-_outgoing1[ix]);
      extpart[3] = getParticleData( _outgoing2[ix]);
      extpart[4] = getParticleData(-_outgoing2[ix]);
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      // first channel always need this
      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(extpart[0],0, 0.0,-1,-2);
      newchannel->addIntermediate(gamma     ,1,-1.1, 1,2);
      newchannel->addIntermediate(gamma     ,1,-1.1, 3,4);
      mode->addChannel(newchannel);
      if(_outgoing1[ix]==_outgoing2[ix])
	{
	  newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	  newchannel->addIntermediate(extpart[0],0, 0.0,-1,-2);
	  newchannel->addIntermediate(gamma     ,1,-1.1, 3,2);
	  newchannel->addIntermediate(gamma     ,1,-1.1, 1,4);
	  mode->addChannel(newchannel);
	  wgt.resize(2);wgt[0]=0.5;wgt[1]=0.5;
	}
      else{wgt.resize(1);wgt[0]=1.;}
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

PScalar4FermionsDecayer::~PScalar4FermionsDecayer() {}

int PScalar4FermionsDecayer::modeNumber(bool & cc,const DecayMode & dm) const
{
  int imode(-1);
  // must be four outgoing particles
  if(dm.products().size()!=4){return imode;}
  // get the id's of the outgoing particles
  int id[4]; bool done[4]; unsigned int ix(0),iy(0);
  // ids of the particles
  int id0(dm.parent()->id()),idtemp(-1),idl1(-1),idl2(-1),idt[2];
  ParticleMSet::const_iterator pit = dm.products().begin();
  for ( ;pit!=dm.products().end();++pit)
    {id[ix]=(**pit).id();done[ix]=false;++ix;}
  // find the two lepton pairs
  // find the first fermion
  ix=0;
  do{if(id[ix]>0&!done[ix]){done[ix]=true;idtemp=id[ix];}++ix;}
  while(ix<4&&idtemp<0);
  if(idtemp<0){return imode;}
  // find its antiparticle
  ix=0;
  do{if(id[ix]==-idtemp&!done[ix]){done[ix]=true;idl1=idtemp;}++ix;}
  while(ix<4&&idl1<0);
  if(idl1<0){return imode;}
  // find the second particle antiparticle pair
  idtemp=-1;
  for(ix=0;ix<4;++ix){if(!done[ix]){idt[iy]=id[ix];++iy;}}
  if(idt[0]==-idt[1]){idl2=abs(idt[0]);}
  if(idl2<0){return imode;}
  // loop over the modes and see if this is one of them
  ix=0;
  do
    {
      if(_incoming[ix]==id0)
	{
	  if((idl1==_outgoing1[ix]&&idl2==_outgoing2[ix])||
	     (idl2==_outgoing1[ix]&&idl1==_outgoing2[ix])){imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  cc=false;
  return imode;
}

void PScalar4FermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,1/MeV) 
     << _incoming << _outgoing1 << _outgoing2 << _maxweight 
     << _includeVMD << _VMDid 
     << ounit(_VMDmass,MeV) << ounit(_VMDwidth,MeV);
}

void PScalar4FermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,1/MeV) 
     >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight 
     >> _includeVMD >> _VMDid 
     >> iunit(_VMDmass,MeV) >> iunit(_VMDwidth,MeV);
}

ClassDescription<PScalar4FermionsDecayer> PScalar4FermionsDecayer::initPScalar4FermionsDecayer;
// Definition of the static class description member.

void PScalar4FermionsDecayer::Init() {

  static ClassDocumentation<PScalar4FermionsDecayer> documentation
    ("The PScalar4FermionsDecayer class is designed for the decay"
     " of a pseudosclar meson to four fermions. It is intended for the decay of"
     "the pion to two electron-positron pairs.");

  static ParVector<PScalar4FermionsDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalar4FermionsDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,int> interfaceOutcoming1
    ("Outgoing1",
     "The PDG code for the first outgoing fermion",
     &PScalar4FermionsDecayer::_outgoing1,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,int> interfaceOutcoming2
    ("Outgoing2",
     "The PDG code for the second outgoing fermion",
     &PScalar4FermionsDecayer::_outgoing2,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalar4FermionsDecayer::_coupling,
     1/MeV, 0, 0/MeV, -10000/MeV, 10000/MeV, false, false, true);

  static ParVector<PScalar4FermionsDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalar4FermionsDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,int> interfaceIncludeVMD
    ("IncludeVMD",
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &PScalar4FermionsDecayer::_includeVMD,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,int> interfaceVMDID
    ("VMDID",
     "The PDG code for the particle to be used for the VMD factor.",
     &PScalar4FermionsDecayer::_VMDid,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,Energy> interfaceVMDmass
    ("VMDmass",
     "The mass to use for the particle in the VMD factor",
     &PScalar4FermionsDecayer::_VMDmass,
     0*MeV, 0, 0*MeV, -10000*MeV, 10000*MeV, false, false, true);

  static ParVector<PScalar4FermionsDecayer,Energy> interfaceVMDwidth
    ("VMDwidth",
     "The width to use for the particle in the VMD factor",
     &PScalar4FermionsDecayer::_VMDwidth,
     0*MeV, 0, 0*MeV, -10000*MeV, 10000*MeV, false, false, true);

}

double PScalar4FermionsDecayer::me2(bool vertex, const int,
				    const Particle & inpart,
				    const ParticleVector & decay) const
{
  bool identical((_outgoing1[imode()]==_outgoing2[imode()]));
  // workaround for gcc 3.2.3 bug
  //ALB ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);
  // vectors for the spinors
  vector<LorentzSpinor<SqrtEnergy> > wave[2];
  vector<LorentzSpinorBar<SqrtEnergy> > wavebar[2];

  // workaround for gcc 3.2.3 bug
  // set up the spin info for the outgoing particles
  for(unsigned int ix=0;ix<2;++ix)
    {
      //ALB SpinorBarWaveFunction(wavebar[ix],decay[2*ix  ],outgoing,true,vertex);
      //ALB SpinorWaveFunction   (   wave[ix],decay[2*ix+1],outgoing,true,vertex);
      vector<LorentzSpinorBar<SqrtEnergy> > mytempLSbar;
      SpinorBarWaveFunction(mytempLSbar,decay[2*ix],outgoing,true,vertex);
      wavebar[ix]=mytempLSbar;
      vector<LorentzSpinor<SqrtEnergy> > mytempLS;
      SpinorWaveFunction(mytempLS,decay[2*ix+1],outgoing,true,vertex);
      wave[ix]=mytempLS;
    }

  // momenta of the outgoing photons
  Lorentz5Momentum momentum[4];
  momentum[0]=decay[0]->momentum()+decay[1]->momentum();momentum[0].rescaleMass();
  momentum[1]=decay[2]->momentum()+decay[3]->momentum();momentum[1].rescaleMass();
  if(identical)
    {
      momentum[2]=decay[2]->momentum()+decay[1]->momentum();momentum[2].rescaleMass();
      momentum[3]=decay[0]->momentum()+decay[3]->momentum();momentum[3].rescaleMass();
    }
  // compute the currents for the two leptonic decays
  LorentzPolarizationVectorE current[4][2][2];
  unsigned int it,ix,iy,iz;
  for(iz=0;iz<2;++iz)
    {
      if(iz==0){it=1;}
      else{it=0;}
      for(ix=0;ix<2;++ix)
	{
	  for(iy=0;iy<2;++iy)
	    {
	      current[iz][iy][ix] = wave[iz][ix].vectorCurrent(wavebar[iz][iy]);
	      // the second two currents      
	      if(identical)
		{current[iz+2][iy][ix] = wave[it][ix].vectorCurrent(wavebar[iz][iy]);}
	    }
	}
    }
  // invariants
  Energy2 m12(momentum[0].mass()*momentum[0].mass());
  Energy2 m34(momentum[1].mass()*momentum[1].mass());
  Energy2 m14(0.*MeV2), m23(0.*MeV2);
  complex<InvEnergy4> prop1(1./m12/m34),prop2(0./sqr(MeV2));
  Complex ii(0.,1.);
  if(identical)
    {
      m14=momentum[2].mass()*momentum[2].mass();
      m23=momentum[3].mass()*momentum[3].mass();
      prop2=1./m14/m23;
    }
  // the VMD factor if needed
  if(_includeVMD[imode()]>0) {
    Energy2 mrho2(_VMDmass[imode()]*_VMDmass[imode()]);
    Energy2 mwrho(_VMDmass[imode()]*_VMDwidth[imode()]);
    prop1 = prop1*(-mrho2+ii*mwrho)/(m12-mrho2+ii*mwrho)*
                  (-mrho2+ii*mwrho)/(m34-mrho2+ii*mwrho);
    if(identical) {
      prop2 = prop2*(-mrho2+ii*mwrho)/(m14-mrho2+ii*mwrho)*
	            (-mrho2+ii*mwrho)/(m23-mrho2+ii*mwrho);
    }
  }
  // prefactor
  Complex pre(_coupling[imode()]*4.*Constants::pi
	      *SM().alphaEM()*inpart.mass());
  Complex diag;
  // now compute the matrix element
  LorentzVector<complex<Energy3> > eps;
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half,
			   PDT::Spin1Half,PDT::Spin1Half);
  vector<unsigned int> ispin(5,0);
  for(ispin[1]=0; ispin[1]<2;++ispin[1])
    {
      for(ispin[2]=0;ispin[2]<2;++ispin[2])
	{
	  for(ispin[3]=0;ispin[3]<2;++ispin[3])
	    {
	      for(ispin[4]=0;ispin[4]<2;++ispin[4])
		{
		  // the first diagram
		  eps = Helicity::
		    epsilon(current[0][ispin[1]][ispin[2]],
					 momentum[1],
					 current[1][ispin[3]][ispin[4]]);
		  diag = prop1*(eps*momentum[0]);
		  // exchanged diagram if identical particles
		  //  (sign due normal ordering) 
	          if(identical)
		    {
		      eps = Helicity::
			epsilon(current[2][ispin[1]][ispin[4]],
					     momentum[3],
					     current[3][ispin[3]][ispin[2]]);
		      diag-= prop2*(eps*momentum[2]);
		    }
		  newME(ispin)=pre*diag;
		}
	    }
	}
    }
  ME(newME);
  RhoDMatrix rhoin=RhoDMatrix(PDT::Spin0);rhoin.average();
  double me=newME.contract(rhoin).real();
  if(identical){me*=0.25;}
  return me;
}

// output the setup info for the particle database
void PScalar4FermionsDecayer::dataBaseOutput(ofstream & output,
					     bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":Incoming   " << ix << " " 
		 << _incoming[ix]   << "\n";
	  output << "set " << fullName() << ":Outgoing1  " << ix << " " 
		 << _outgoing1[ix]  << "\n";
	  output << "set " << fullName() << ":Outgoing2  " << ix << " " 
		 << _outgoing2[ix]  << "\n";
	  output << "set " << fullName() << ":Coupling   " << ix << " " 
		 << _coupling[ix]*MeV   << "\n";
	  output << "set " << fullName() << ":MaxWeight  " << ix << " " 
		 << _maxweight[ix]  << "\n";
	  output << "set " << fullName() << ":IncludeVMD " << ix << " " 
		 << _includeVMD[ix] << "\n";
	  output << "set " << fullName() << ":VMDID      " << ix << " " 
		 << _VMDid[ix]      << "\n";
	  output << "set " << fullName() << ":VMDmass    " << ix << " " 
		 << _VMDmass[ix]/MeV    << "\n";
	  output << "set " << fullName() << ":VMDwidth   " << ix << " " 
		 << _VMDwidth[ix]/MeV   << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":Incoming   " << ix << " " 
		 << _incoming[ix]   << "\n";
	  output << "insert " << fullName() << ":Outgoing1  " << ix << " " 
		 << _outgoing1[ix]  << "\n";
	  output << "insert " << fullName() << ":Outgoing2  " << ix << " " 
		 << _outgoing2[ix]  << "\n";
	  output << "insert " << fullName() << ":Coupling   " << ix << " " 
		 << _coupling[ix]*MeV   << "\n";
	  output << "insert " << fullName() << ":MaxWeight  " << ix << " " 
		 << _maxweight[ix]  << "\n";
	  output << "insert " << fullName() << ":IncludeVMD " << ix << " " 
		 << _includeVMD[ix] << "\n";
	  output << "insert " << fullName() << ":VMDID      " << ix << " " 
		 << _VMDid[ix]      << "\n";
	  output << "insert " << fullName() << ":VMDmass    " << ix << " " 
		 << _VMDmass[ix]/MeV    << "\n";
	  output << "insert " << fullName() << ":VMDwidth   " << ix << " " 
		 << _VMDwidth[ix]/MeV   << "\n";
	}
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
