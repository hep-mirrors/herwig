// -*- C++ -*-
//
// OniumToOniumPiPiDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OniumToOniumPiPiDecayer class.
//

#include "OniumToOniumPiPiDecayer.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void OniumToOniumPiPiDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  for(unsigned int ix=0;ix<_maxweight.size();++ix) {
    if(initialize()) _maxweight[ix] = mode(ix)->maxWeight();
  }
}

OniumToOniumPiPiDecayer::OniumToOniumPiPiDecayer() {
  // Upsilon(3S)->Upsilon(1S) pi pi
  _incoming.push_back(200553);
  _outgoing.push_back(   553);
  _maxweight.push_back(1.);
  _maxweight.push_back(1.);
  _coupling.push_back(3.92e-6);
  _reA.push_back( 1.   /MeV2);_imA.push_back( ZERO);
  _reB.push_back(-2.523/MeV2);_imB.push_back( 1.189/MeV2);
  _reC.push_back( ZERO);_imC.push_back( ZERO);
  // Upsilon(3S)->Upsilon(2S) pi pi
  _incoming.push_back(200553);
  _outgoing.push_back(100553);
  _maxweight.push_back(1.);
  _maxweight.push_back(1.);
  _coupling.push_back(311e-6);
  _reA.push_back( 1.   /MeV2);_imA.push_back( ZERO);
  _reB.push_back(-0.395/MeV2);_imB.push_back( 0.001/MeV2);
  _reC.push_back( ZERO);_imC.push_back( ZERO);
  // Upsilon(2S)->Upsilon(1S) pi pi
  _incoming.push_back(100553);
  _outgoing.push_back(   553);
  _maxweight.push_back(1.);
  _maxweight.push_back(1.);
  _coupling.push_back(61.4e-6);
  _reA.push_back( 1.   /MeV2);_imA.push_back( ZERO);
  _reB.push_back(-0.753/MeV2);_imB.push_back( ZERO);
  _reC.push_back( ZERO);_imC.push_back( ZERO);
  // Upsilon(4S)->Upsilon(1S) pi pi
  _incoming.push_back(300553);
  _outgoing.push_back(   553);
  _maxweight.push_back(1.);
  _maxweight.push_back(1.);
  _coupling.push_back(1.77e-6);
  _reA.push_back( 1.   /MeV2);_imA.push_back( ZERO);
  _reB.push_back( ZERO);_imB.push_back( ZERO);
  _reC.push_back( ZERO);_imC.push_back( ZERO);
  // Upsilon(4S)->Upsilon(2S) pi pi
  _incoming.push_back(300553);
  _outgoing.push_back(100553);
  _maxweight.push_back(1.);
  _maxweight.push_back(1.);
  _coupling.push_back(68.8e-6);
  _reA.push_back( 1.   /MeV2);_imA.push_back( ZERO);
  _reB.push_back(-2.35   /MeV2);_imB.push_back( 0.55/MeV2);
  _reC.push_back( ZERO);_imC.push_back( ZERO);
  // psi(2s)->psi(1S) pi pi
  _incoming.push_back(100443);
  _outgoing.push_back(   443);
  _maxweight.push_back(1.);
  _maxweight.push_back(1.);
  _coupling.push_back(66.2e-6);
  _reA.push_back( 1.   /MeV2);_imA.push_back( ZERO);
  _reB.push_back(-0.336/MeV2);_imB.push_back( ZERO);
  _reC.push_back( ZERO);_imC.push_back( ZERO);
  // psi(3770)->psi(1S) pi pi
  _incoming.push_back(30443);
  _outgoing.push_back(   443);
  _maxweight.push_back(1.);
  _maxweight.push_back(1.);
  _coupling.push_back(20.6e-6);
  _reA.push_back( 1.   /MeV2);_imA.push_back( ZERO);
  _reB.push_back( ZERO);_imB.push_back( ZERO);
  _reC.push_back( ZERO);_imC.push_back( ZERO);
  // Initial size of the vectors
  _initsize=_incoming.size();
  // don'y generate the intermediates in the phase-space
  generateIntermediates(false);
}

void OniumToOniumPiPiDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistency of the vectors
  unsigned int isize=_incoming.size();
  if(_outgoing.size()!=isize||_maxweight.size()!=2*isize||
     _coupling.size()!=isize||
     _reA     .size()!=isize||_imA.size()      !=isize||
     _reB     .size()!=isize||_imB.size()      !=isize||
     _reC     .size()!=isize||_imC.size()      !=isize)
    throw InitException() << "Inconsistent size of the parameter vectors in "
			  << "OniumToOniumPiPiDecayer"
			  << Exception::runerror;
  // construct the complex couplings
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    _cA.push_back(complex<InvEnergy2>(_reA[ix],_imA[ix]));
    _cB.push_back(complex<InvEnergy2>(_reB[ix],_imB[ix]));
    _cC.push_back(complex<InvEnergy2>(_reC[ix],_imC[ix]));
  }
  // construct the decay channels
  tPDVector extpart(4);
  tPDPtr pip(getParticleData(ParticleID::piplus ));
  tPDPtr pim(getParticleData(ParticleID::piminus));
  tPDPtr pi0(getParticleData(ParticleID::pi0    ));
  tPDPtr rho0(getParticleData(113)); 
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr newchannel;
  vector<double> dummyweights(1,1.);
  for(unsigned int ix=0;ix<isize;++ix) {
    extpart[0]=getParticleData(_incoming[ix]);
    extpart[1]=getParticleData(_outgoing[ix]);
    for(unsigned int iy=0;iy<2;++iy) {
      // pi+ pi-
      if(iy==0) {
	extpart[2]=pip;
	extpart[3]=pim;
      }
      // pi0 pi0
      else {
	extpart[2]=pi0;
	extpart[3]=pi0;
      }
      // construct the phase-space mode
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(extpart[0],0, 0.0,1,-1);
      newchannel->addIntermediate(rho0,1,0.0, 2,3);
      mode->addChannel(newchannel);
      // reset the resonance parameters
      mode->resetIntermediate(rho0,2*extpart[0]->mass(),2*extpart[0]->mass());
      // add the mode
      addMode(mode,_maxweight[2*ix+iy],dummyweights);
    }
  }
}


void OniumToOniumPiPiDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoing << _maxweight << _initsize << ounit(_reA,1./GeV2) 
     << ounit(_imA,1./GeV2) << ounit(_cA,1./GeV2) << ounit(_reB,1./GeV2) 
     << ounit(_imB,1./GeV2) << ounit(_cB,1./GeV2) << ounit(_reC,1./GeV2) 
     << ounit(_imC,1./GeV2) << ounit(_cC,1./GeV2);
}

void OniumToOniumPiPiDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoing >> _maxweight >> _initsize >> iunit(_reA,1./GeV2) 
     >> iunit(_imA,1./GeV2) >> iunit(_cA,1./GeV2) >> iunit(_reB,1./GeV2) 
     >> iunit(_imB,1./GeV2) >> iunit(_cB,1./GeV2) >> iunit(_reC,1./GeV2) 
     >> iunit(_imC,1./GeV2) >> iunit(_cC,1./GeV2);
}

ClassDescription<OniumToOniumPiPiDecayer> OniumToOniumPiPiDecayer::initOniumToOniumPiPiDecayer;
// Definition of the static class description member.

void OniumToOniumPiPiDecayer::Init() {

  static ClassDocumentation<OniumToOniumPiPiDecayer> documentation
    ("The OniumToOniumPiPiDecayer class uses the matrix element of "
     "Brown and Cahn, PRL35, 1 (1975), for"
     " the decay of onium resonaces to lighter states and pion pairs."
     " The results of hep-ex/9909038 are used for psi'->psi and "
     " arXiv:0706.2317 for Upsilon(3S) and Upsilon(2S) decays."
     " The remaining parameters are choosen to approximately reproduce"
     " the distributions from hep-ex/0604031 and hep-ex/0508023.",
     "The decays of onium resonances to lighter states and pion pairs were modelled"
     " using the matrix element of \\cite{Brown:1975dz}. The results of "
     "\\cite{Bai:1999mj} are used for $\\psi'\\to\\psi$ and "
     "\\cite{Cronin-Hennessy:2007sj} for $\\Upsilon(3S)$ and $\\Upsilon(2S)$ decays."
     " The remaining parameters are choosen to approximately reproduce"
     " the distributions from \\cite{Aubert:2006bm} and \\cite{Adam:2005mr}.",
     "\\bibitem{Brown:1975dz} L.~S.~Brown and R.~N.~Cahn,"
     "Phys.\\ Rev.\\ Lett.\\  {\\bf 35} (1975) 1."
     "%%CITATION = PRLTA,35,1;%%\n"
     "\\bibitem{Bai:1999mj} J.~Z.~Bai {\\it et al.}  [BES Collaboration],"
     "Phys.\\ Rev.\\  D {\\bf 62} (2000) 032002 [arXiv:hep-ex/9909038]."
     "%%CITATION = PHRVA,D62,032002;%%\n"
     "\\bibitem{Cronin-Hennessy:2007sj} D.~Cronin-Hennessy{\\it et al.} "
     "[CLEO Collaboration], arXiv:0706.2317 [hep-ex]."
     "%%CITATION = ARXIV:0706.2317;%%\n"
     "\\bibitem{Aubert:2006bm} B.~Aubert {\\it et al.}  [BABAR Collaboration],"
     "Phys.\\ Rev.\\ Lett.\\  {\\bf 96} (2006) 232001 [arXiv:hep-ex/0604031]."
     "%%CITATION = PRLTA,96,232001;%%\n"
     "\\bibitem{Adam:2005mr} N.~E.~Adam {\\it et al.}  [CLEO Collaboration],"
     "Phys.\\ Rev.\\ Lett.\\  {\\bf 96} (2006) 082004 [arXiv:hep-ex/0508023]."
     "%%CITATION = PRLTA,96,082004;%%");

  static ParVector<OniumToOniumPiPiDecayer,long> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming onium state",
     &OniumToOniumPiPiDecayer::_incoming, -1, long(0), -10000000, 10000000,
     false, false, Interface::limited);

  static ParVector<OniumToOniumPiPiDecayer,long> interfaceOutgoing
    ("Outgoing",
     "The PDG code for the outgoing onium state",
     &OniumToOniumPiPiDecayer::_outgoing, -1, long(0), -10000000, 10000000,
     false, false, Interface::limited);

  static ParVector<OniumToOniumPiPiDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode, there should be two "
     "for each mode as we have pi+ pi- and pi0 pi0",
     &OniumToOniumPiPiDecayer::_maxweight, -1, 1.0, 0.0, 10000.0,
     false, false, Interface::limited);

  static ParVector<OniumToOniumPiPiDecayer,double> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay",
     &OniumToOniumPiPiDecayer::_coupling, -1, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static ParVector<OniumToOniumPiPiDecayer,InvEnergy2> interfaceReA
    ("ReA",
     "The real part of the A coupling",
     &OniumToOniumPiPiDecayer::_reA, 1./MeV2, -1, 1.0/MeV2, -1000.0/MeV2, 1000.0/MeV2,
     false, false, Interface::limited);

  static ParVector<OniumToOniumPiPiDecayer,InvEnergy2> interfaceImA
    ("ImA",
     "The imaginary part of the A coupling",
     &OniumToOniumPiPiDecayer::_imA, 1./MeV2, -1, 1.0/MeV2, -1000.0/MeV2, 1000.0/MeV2,
     false, false, Interface::limited);

  static ParVector<OniumToOniumPiPiDecayer,InvEnergy2> interfaceReB
    ("ReB",
     "The real part of the B coupling",
     &OniumToOniumPiPiDecayer::_reB, 1./MeV2, -1, 1.0/MeV2, -1000.0/MeV2, 1000.0/MeV2,
     false, false, Interface::limited);

  static ParVector<OniumToOniumPiPiDecayer,InvEnergy2> interfaceImB
    ("ImB",
     "The imaginary part of the B coupling",
     &OniumToOniumPiPiDecayer::_imB, 1./MeV2, -1, 1.0/MeV2, -1000.0/MeV2, 1000.0/MeV2,
     false, false, Interface::limited);

  static ParVector<OniumToOniumPiPiDecayer,InvEnergy2> interfaceReC
    ("ReC",
     "The real part of the C coupling",
     &OniumToOniumPiPiDecayer::_reC, 1./MeV2, -1, 1.0/MeV2, -1000.0/MeV2, 1000.0/MeV2,
     false, false, Interface::limited);

  static ParVector<OniumToOniumPiPiDecayer,InvEnergy2> interfaceImC
    ("ImC",
     "The imaginary part of the C coupling",
     &OniumToOniumPiPiDecayer::_imC, 1./MeV2, -1, 1.0/MeV2, -1000.0/MeV2, 1000.0/MeV2,
     false, false, Interface::limited);
}

int OniumToOniumPiPiDecayer::modeNumber(bool & cc,tcPDPtr parent,
					const tPDVector & children) const {
  cc=false;
  int imode(-1);
  long idin(parent->id());
  if(children.size()!=3) return -1;
  unsigned int npip(0),npim(0),npi0(0);
  long idother(0),id;
  for(tPDVector::const_iterator pit=children.begin();
      pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)       ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else if(id==ParticleID::pi0)     ++npi0;
    else idother=id;
  }
  // check pi+ pi- or pi0 pi0 and outgoing state
  if(!((npip==1&&npim==1)||npi0==2)||idother==0) return -1;
  unsigned int ix=0;
  do {
    if(idin==_incoming[ix]&&idother==_outgoing[ix]) imode=ix;
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return npi0==2 ? 2*imode+1 : 2*imode;
}

double OniumToOniumPiPiDecayer::me2(const int,
				    const Particle & inpart,
				    const ParticleVector & decay,
				    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors[0],_rho,
						const_ptr_cast<tPPtr>(&inpart),
						incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors[0],const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    VectorWaveFunction::constructSpinInfo(_vectors[1],decay[0],
					  outgoing,true,false);
    for(unsigned int ix=1;ix<3;++ix)
      ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
    return 0.;
  }
  VectorWaveFunction::calculateWaveFunctions(_vectors[1],decay[0],outgoing,false);
  // compute the matrix element
  complex<InvEnergy2> A(_cA[imode()/2]),B(_cB[imode()/2]),C(_cC[imode()/2]);
  Energy2 q2  =(decay[1]->momentum()+decay[2]->momentum()).m2();
  Energy2 mpi2=sqr(decay[1]->mass());
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      Complex dota = _vectors[0][ix].dot(_vectors[1][iy]);
      complex<Energy2> dotb = 
	(_vectors[0][ix]*decay[1]->momentum())*(_vectors[1][iy]*decay[2]->momentum())+
	(_vectors[0][ix]*decay[2]->momentum())*(_vectors[1][iy]*decay[1]->momentum());
      (*ME())(ix,iy,0,0)= _coupling[imode()/2]*
	Complex(A*dota*(q2-2.*mpi2)+B*dota*decay[1]->momentum().e()*decay[2]->momentum().e()
	 +C*dotb);
    }
  }
  // matrix element
  double output=ME()->contract(_rho).real();
  if(imode()%2==1) output*=0.5;
  // test of the matrix element
//   Energy2 s1=(decay[1]->momentum()+decay[2]->momentum()).m2();
//   Energy2 s2=(decay[0]->momentum()+decay[2]->momentum()).m2();
//   Energy2 s3=(decay[1]->momentum()+decay[0]->momentum()).m2();
//   double test=threeBodyMatrixElement(imode(),sqr(inpart.mass()),
// 				     s3,s2,s1,decay[0]->mass(),
// 				     decay[1]->mass(),decay[2]->mass());
  // return the answer
  return output;
}

// output the setup information for the particle database
void OniumToOniumPiPiDecayer::dataBaseOutput(ofstream & output,
					     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "newdef " << name() << ":Outgoing " << ix << " " 
	     << _outgoing[ix] << "\n";
      output << "newdef " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix] << "\n";
      output << "newdef " << name() << ":ReA " << ix << " " 
	     << _reA[ix]*MeV2 << "\n";
      output << "newdef " << name() << ":ImA " << ix << " " 
	     << _imA[ix]*MeV2 << "\n";
      output << "newdef " << name() << ":ReB " << ix << " " 
	     << _reB[ix]*MeV2 << "\n";
      output << "newdef " << name() << ":ImB " << ix << " " 
	     << _imB[ix]*MeV2 << "\n";
      output << "newdef " << name() << ":ReC " << ix << " " 
	     << _reC[ix]*MeV2 << "\n";
      output << "newdef " << name() << ":ImC " << ix << " " 
	     << _imC[ix]*MeV2 << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << name() << ":Outgoing " << ix << " " 
	     << _outgoing[ix] << "\n";
      output << "insert " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix] << "\n";
      output << "insert " << name() << ":ReA " << ix << " " 
	     << _reA[ix]*MeV2 << "\n";
      output << "insert " << name() << ":ImA " << ix << " " 
	     << _imA[ix]*MeV2 << "\n";
      output << "insert " << name() << ":ReB " << ix << " " 
	     << _reB[ix]*MeV2 << "\n";
      output << "insert " << name() << ":ImB " << ix << " " 
	     << _imB[ix]*MeV2 << "\n";
      output << "insert " << name() << ":ReC " << ix << " " 
	     << _reC[ix]*MeV2 << "\n";
      output << "insert " << name() << ":ImC " << ix << " " 
	     << _imC[ix]*MeV2 << "\n";
    }
  }
  for(unsigned int ix=0;ix<_maxweight.size();++ix) {
    if(ix<2*_initsize) {
      output << "newdef " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

WidthCalculatorBasePtr OniumToOniumPiPiDecayer::
threeBodyMEIntegrator(const DecayMode & dm) const {
  int imode(-1);
  long idin(dm.parent()->id());
  unsigned int npip(0),npim(0),npi0(0);
  long idother(0),id;
  for(ParticleMSet::const_iterator pit=dm.products().begin();
      pit!=dm.products().end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)       ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else if(id==ParticleID::pi0)     ++npi0;
    else idother=id;
  }
  unsigned int ix=0;
  do {
    if(idin==_incoming[ix]&&idother==_outgoing[ix]) imode=ix;
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  imode = npi0==2 ? 2*imode+1 : 2*imode;
  // construct the integrator
  vector<double> inweights(1,1.);
  Energy scale=getParticleData(_incoming[ix-1])->mass();
  Energy m1=getParticleData(_outgoing[ix-1])->mass();
  Energy mpi = npi0==2 ? getParticleData(ParticleID::pi0)->mass() :
    getParticleData(ParticleID::piplus)->mass();
  vector<int> intype(1,3);
  vector<Energy> inmass (1,scale);
  vector<Energy> inwidth(1,scale);
  vector<double> inpow(1,0.0);
  return new_ptr(ThreeBodyAllOnCalculator<OniumToOniumPiPiDecayer>
		 (inweights,intype,inmass,inwidth,inpow,
		  *this,imode,m1,mpi,mpi));
}

double OniumToOniumPiPiDecayer::
threeBodyMatrixElement(const int imode, const Energy2 q2,
		       const  Energy2 s3, const Energy2 s2, const Energy2 s1, 
		       const Energy m1, const Energy m2, const Energy m3) const {
  Energy q=sqrt(q2);
  Energy e2 = 0.5*(q2+sqr(m2)-s2)/q;
  Energy e3 = 0.5*(q2+sqr(m3)-s3)/q;
  Complex amp = _cA[imode/2]*(s1-sqr(m2)-sqr(m3))+_cB[imode/2]*e2*e3;
  Energy2 dot = 0.5*(q2+sqr(m1)-s1);
  double output=(2.+sqr(dot/q/m1))*real(amp*conj(amp))*sqr(_coupling[imode/2])/3.;
  if(imode%2==1) output*=0.5;
  return output;
}
