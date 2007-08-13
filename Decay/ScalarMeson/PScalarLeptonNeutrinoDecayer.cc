// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarLeptonNeutrinoDecayer class.
//

#include "PScalarLeptonNeutrinoDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

PScalarLeptonNeutrinoDecayer::PScalarLeptonNeutrinoDecayer() 
  : _incoming(6), _decayconstant(6), _leptons(6), _maxweighte(6), 
    _maxweightmu(6), _maxweighttau(6), _GF(1.16637E-5/GeV2) {
  // pion decay
  _incoming[0] = 211; _decayconstant[0] = 127.3*MeV; _leptons[0] = 2; 
  _maxweighte[0] = 0.000126639; _maxweightmu[0] = 0.986777; 
  _maxweighttau[0] = 0.0; 
  // kaon decay
  _incoming[1] = 321; _decayconstant[1] = 36.31*MeV; _leptons[1] = 2; 
  _maxweighte[1] = 1.6599e-05; _maxweightmu[1] = 0.646136; 
  _maxweighttau[1] = 0.0; 
  // D_s decay
  _incoming[2] = 431; _decayconstant[2] = 286.1*MeV; _leptons[2] = 3; 
  _maxweighte[2] = 8.29008e-08; _maxweightmu[2] = 0.00352387; 
  _maxweighttau[2] = 0.0342748; 
  // D decay
  _incoming[3] = 411; _decayconstant[3] = 50.55*MeV; _leptons[3] = 3; 
  _maxweighte[3] = 1.67007e-07; _maxweightmu[3] = 0.00709453; 
  _maxweighttau[3] = 0.0187692; 
  // B_c decay
  _incoming[4] = 541; _decayconstant[4] = 16.0*MeV; _leptons[4] = 3; 
  _maxweighte[4] = 1.61603e-09; _maxweightmu[4] = 6.90525e-05; 
  _maxweighttau[4] = 0.0166333; 
  // B_u decays
  _incoming[5] = 521; _decayconstant[5] = 0.6970*MeV; _leptons[5] = 3; 
  _maxweighte[5] = 1.07731e-11; _maxweightmu[5] = 4.60215e-07; 
  _maxweighttau[5] = 9.31156e-05; 
  // initial size of the vectors
  _initsize=_incoming.size();
  // intermediates
  generateIntermediates(false);
}

void PScalarLeptonNeutrinoDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize(_incoming.size());
  if(isize!=_decayconstant.size()||isize!=_leptons.size()||isize!=_maxweighte.size()||
     isize!=_maxweightmu.size()||isize!=_maxweighttau.size())
    throw InitException() << "Inconsistent parameters in PScalarLeptonNeutrinoDecayer"
			  << Exception::abortnow;
  // create the integration channels
  PDVector extpart(3);  
  tPDPtr nu[3]={getParticleData(ParticleID::nu_e),
		getParticleData(ParticleID::nu_mu),
		getParticleData(ParticleID::nu_tau)};
  tPDPtr nubar[3]={getParticleData(ParticleID::nu_ebar),
		   getParticleData(ParticleID::nu_mubar),
		   getParticleData(ParticleID::nu_taubar)};
  tPDPtr lep[3]={getParticleData(ParticleID::eminus),
		 getParticleData(ParticleID::muminus),
		 getParticleData(ParticleID::tauminus)};
  tPDPtr lepbar[3]={getParticleData(ParticleID::eplus),
		    getParticleData(ParticleID::muplus),
		    getParticleData(ParticleID::tauplus)};
  int charge;
  vector<double> dummyweights;
  double wgt;
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0]=getParticleData(_incoming[ix]);
    charge=extpart[0]->iCharge();
    for(unsigned int iy=0;iy<_leptons[ix];++iy) {
      if(charge>0) {
	extpart[1]=lepbar[iy];
	extpart[2]=nu[iy];
      }
      else {
	extpart[1]=lep[iy];
	extpart[2]=nubar[iy];
      }
      mode = new DecayPhaseSpaceMode(extpart,this);
      if(iy==0)      wgt = _maxweighte[ix];
      else if(iy==1) wgt = _maxweightmu[ix];
      else           wgt = _maxweighttau[ix];
      addMode(mode,wgt,dummyweights);
    }
  }
}

int PScalarLeptonNeutrinoDecayer::modeNumber(bool & cc,tcPDPtr parent,
					     const PDVector & children) const {
  int imode(-1);
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0(parent->id()),id0bar(id0);
  if(parent->CC()){id0bar=-id0;}
  PDVector::const_iterator pit = children.begin();
  int id;
  unsigned ilep(4);
  for(;pit!=children.end();++pit) {
    id=abs((**pit).id());
    if(id>=11&&id<=16&&id%2==0) ilep=(id-10)/2;
  }
  // find the channel we need
  bool found(false);
  int ichan(-1);
  unsigned int ix(0);
  do {
    if(id0   ==_incoming[ix]||id0bar==_incoming[ix]) {
      found=true;ichan+=ilep;
      cc=id0bar==_incoming[ix];
    }
    else {
      ichan+=_leptons[ix];
    }
    ++ix;
  }
  while (!found&&ix<_incoming.size());
  if(found) imode=ichan;
  return imode;
}


void PScalarLeptonNeutrinoDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << ounit(_decayconstant,GeV)
     << _leptons << _maxweighte << _maxweightmu 
     << _maxweighttau << ounit(_GF,1/GeV2);
}

void PScalarLeptonNeutrinoDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> iunit(_decayconstant,GeV) 
     >> _leptons >> _maxweighte >> _maxweightmu 
     >> _maxweighttau >> iunit(_GF,1/GeV2);
}

ClassDescription<PScalarLeptonNeutrinoDecayer> PScalarLeptonNeutrinoDecayer::initPScalarLeptonNeutrinoDecayer;
// Definition of the static class description member.

void PScalarLeptonNeutrinoDecayer::Init() {

  static ClassDocumentation<PScalarLeptonNeutrinoDecayer> documentation
    ("The PScalarLeptonNeutrinoDecayer class is the base class for"
     " the implementation of leptonic decay of a pseudoscalar meson to a lepton"
     "and a neutrino.");

  static ParVector<PScalarLeptonNeutrinoDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalarLeptonNeutrinoDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarLeptonNeutrinoDecayer,unsigned int> interfaceLepton
    ("Leptons",
     "The allowed outgoing leptons, 1 is only electrons, 2 is electrons and muons and "
     "3 is all lepton flavours",
     &PScalarLeptonNeutrinoDecayer::_leptons,
     0, 0, 0, 1, 3, false, false, true);

  static ParVector<PScalarLeptonNeutrinoDecayer,double> interfaceMaxWeighte
    ("MaxWeightElectron",
     "The maximum weight for the electron decay mode",
     &PScalarLeptonNeutrinoDecayer::_maxweighte,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<PScalarLeptonNeutrinoDecayer,double> interfaceMaxWeightmu
    ("MaxWeightMuon",
     "The maximum weight for the muon decay mode",
     &PScalarLeptonNeutrinoDecayer::_maxweightmu,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<PScalarLeptonNeutrinoDecayer,double> interfaceMaxWeighttau
    ("MaxWeightTau",
     "The maximum weight for the tau decay mode",
     &PScalarLeptonNeutrinoDecayer::_maxweighttau,
     0, 0, 0, 0., 100., false, false, true);

  static Parameter<PScalarLeptonNeutrinoDecayer,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &PScalarLeptonNeutrinoDecayer::_GF, 1./GeV2, 1.16639E-5/GeV2,
     0./GeV2, 1.0e-4*1./GeV2,false, false, false);

  static ParVector<PScalarLeptonNeutrinoDecayer,Energy> interfaceDecayConstant
    ("DecayConstant",
     "The decay constant for the incoming pseudoscaalr meson.",
     &PScalarLeptonNeutrinoDecayer::_decayconstant,
     MeV, 0, 0*MeV, 0*MeV, 1000.*MeV, false, false, true);
}

double PScalarLeptonNeutrinoDecayer::me2(bool vertex, const int,
					 const Particle & inpart,
					 const ParticleVector & decay) const {
  int icoup(0),id(abs(inpart.id())),idferm(0);
  // work out which decay constant to use
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(id==abs(_incoming[ix])) icoup=ix;
  }
  // workaround for gcc 3.2.3 bug
  // check if the decay particle has spin info
  //ALB ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);
  
  // find the particles
  unsigned int iferm(0),ianti(0);
  for(unsigned ix=0;ix<decay.size();++ix) {
    id=decay[ix]->id();
    if(id<=-11&&id>=-16) {
      ianti=ix;
    }
    else if(id>=11&&id<=16) {
      iferm=ix;
      idferm=id;
    }
  }
  // spinors for the lepton and neutrino
  vector<LorentzSpinor<SqrtEnergy> > wave;
  vector<LorentzSpinorBar<SqrtEnergy> > wbar;
  // construct the spininfo's of the outgoing particles
  SpinorWaveFunction(   wave,decay[ianti],outgoing,true,vertex);
  SpinorBarWaveFunction(wbar,decay[iferm],outgoing,true,vertex);
  // the prefactor
  Energy premass;
  if(idferm%2==0) premass=decay[ianti]->mass();
  else            premass=decay[iferm]->mass();
  InvEnergy pre = premass * 2.*_decayconstant[icoup]*_GF/inpart.mass();
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half);
  vector<unsigned int> ispin(3,0);
  for(ispin[ianti+1]=0;ispin[ianti+1]<2;++ispin[ianti+1]) {
    for(ispin[iferm+1]=0;ispin[iferm+1]<2;++ispin[iferm+1]) {
      newME(ispin)= idferm%2==0 ? 
	pre*wave[ispin[ianti+1]].rightScalar(wbar[ispin[iferm+1]]) :
	pre*wave[ispin[ianti+1]].leftScalar( wbar[ispin[iferm+1]]) ;
    }
  }
  ME(newME);
  RhoDMatrix rhoin(PDT::Spin0);rhoin.average();
//   // test of the matrix element
//   double me=newME.contract(rhoin).real();
//   Energy mass = idferm%2==0 ? decay[ianti]->mass() : decay[iferm]->mass();
//   double test = sqr(_decayconstant[icoup]*_GF*2.*mass/inpart.mass())*
//     (sqr(inpart.mass())-sqr(mass));
//   cout << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " " 
//        << me << " " << (me-test)/(me+test) << endl;
  return 0.5*newME.contract(rhoin).real();
}


void PScalarLeptonNeutrinoDecayer::dataBaseOutput(ofstream & output,
						  bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "set " << fullName() << ":GFermi " << _GF*GeV2 << "\n";
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "set " << fullName() << ":Incoming   " << ix << " "
	     << _incoming[ix]   << "\n";
      output << "set " << fullName() << ":Leptons    " << ix << " "
	     << _leptons[ix]   << "\n";
      output << "set " << fullName() << ":MaxWeightElectron " << ix << " "
	     << _maxweighte[ix]   << "\n";
      output << "set " << fullName() << ":MaxWeightMuon "     << ix << " "
	     << _maxweightmu[ix]   << "\n";
      output << "set " << fullName() << ":MaxWeightTau "      << ix << " "
	     << _maxweighttau[ix]   << "\n";
      output << "set " << fullName() << ":DecayConstant "     << ix << " "
	     << _decayconstant[ix]/GeV  << "\n";
    }
    else {
      output << "insert " << fullName() << ":Incoming   " << ix << " "
	     << _incoming[ix]   << "\n";
      output << "insert " << fullName() << ":Leptons    " << ix << " "
	     << _leptons[ix]   << "\n";
      output << "insert " << fullName() << ":MaxWeightElectron " << ix << " "
	     << _maxweighte[ix]   << "\n";
      output << "insert " << fullName() << ":MaxWeightMuon "     << ix << " "
	     << _maxweightmu[ix]   << "\n";
      output << "insert " << fullName() << ":MaxWeightTau "      << ix << " "
	     << _maxweighttau[ix]   << "\n";
      output << "insert " << fullName() << ":DecayConstant "     << ix << " "
	     << _decayconstant[ix]/GeV  << "\n";
    }
  }
  if(header) 
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
