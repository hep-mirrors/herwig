// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2BaryonsDecayer class.
//

#include "VectorMeson2BaryonsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMeson2BaryonsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}

void VectorMeson2BaryonsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=gm_.size();
  if(isize!=incoming_.size()  || isize!=outgoingf_.size()||
     isize!=outgoinga_.size() || isize!= maxweight_.size()||
     isize!= phi_.size() || isize!= ge_.size())
    throw InitException() << "Inconsistent parameters in VectorMeson2"
			   << "BaryonsDecayer::doiin() " << Exception::runerror;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoingf_[ix]),
		     getParticleData(outgoinga_[ix])};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

VectorMeson2BaryonsDecayer::VectorMeson2BaryonsDecayer()
  : gm_       ({0.00163222616377,0.00158881446341,0.00163819673075,0.000944247071286,0.00141162244048 ,0.00150424724997 ,0.00093350341391,0.00107528142563,0.000860759359361,0.00103222902272,0.00098818310509,0.00119542938517 ,0.000941856299908,0.00108249234586 ,0.00102912698738 ,0.000561082932891,0.000526423800011,0.000500391020504,0.000678994938986}),
    ge_       ({0.00135736265293,0.0015117595557 ,0.00136696938558,0.001988067505   ,0.000852659560453,0.000801719253474,0.00150640470181,0.00153402258271,0.0015323974485  ,0.              ,0.00084599445285,0.000621041230885,0.000599393812308,0.000327664954147,0.000664382626732,0.000260326162249,0.000362882855521,0.00025226758188 ,0.000397805035356}),
    phi_      ({0.              ,0.              ,0.740           ,0.               ,0.               ,0.               ,0.              ,0.              ,0.               ,0.              ,0.              ,0.               ,0.               ,0.               ,0.               ,0.               ,0.               ,0.               ,0.               }),
    incoming_ ({443             ,443             ,443             ,443              ,443              ,443              ,443             ,443             ,443              ,100443          ,100443          ,100443           ,100443           ,100443           ,100443           ,100443           ,100443           ,100443           ,100443           }),
    outgoingf_({ 2212           , 2112           , 3122           , 3212            , 3312            , 3322            , 3114           , 3224           , 3214            , 2212           , 2112           , 3122            , 3212            , 3312            , 3322            , 3114            , 3224            , 3214            , 3314            }),
    outgoinga_({-2212           ,-2112           ,-3122           ,-3212            ,-3312            ,-3322            ,-3114           ,-3224           ,-3214            ,-2212           ,-2112           ,-3122            ,-3212            ,-3312            ,-3322            ,-3114            ,-3224            ,-3214            ,-3314            }),
    maxweight_({1.6             ,1.6             ,2.1             ,1.5              ,2.               ,2.1              ,7.              ,10.             ,6.5              ,1.7             ,1.7             ,2.5              ,2.5              ,2.               ,2.               ,30.              ,35.              ,41.              , 0.001           })
{}

int VectorMeson2BaryonsDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(incoming_[ix]==id   ) {
      if((id1   ==outgoingf_[ix]&&id2   ==outgoinga_[ix])||
	 (id2   ==outgoingf_[ix]&&id1   ==outgoinga_[ix])) imode=ix;
    }
    if(incoming_[ix]==idbar) {
      if((id1bar==outgoingf_[ix]&&id2bar==outgoinga_[ix])||
	 (id2bar==outgoingf_[ix]&&id1bar==outgoinga_[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  return imode;
}

IBPtr VectorMeson2BaryonsDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VectorMeson2BaryonsDecayer::fullclone() const {
  return new_ptr(*this);
}

void VectorMeson2BaryonsDecayer::persistentOutput(PersistentOStream & os) const {
  os << gm_ << ge_ << phi_ << incoming_ << outgoingf_ << outgoinga_ << maxweight_;
}

void VectorMeson2BaryonsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> gm_ >> ge_ >> phi_ >> incoming_ >> outgoingf_ >> outgoinga_ >> maxweight_;
}


DescribeClass<VectorMeson2BaryonsDecayer,DecayIntegrator>
describeHerwigVectorMeson2BaryonsDecayer("Herwig::VectorMeson2BaryonsDecayer", "HwVMDecay.so");

void VectorMeson2BaryonsDecayer::Init() {

  static ClassDocumentation<VectorMeson2BaryonsDecayer> documentation
    ("The VectorMeson2BaryonsDecayer class is designed for "
     "the decay of vector mesons to baryons.");

  static ParVector<VectorMeson2BaryonsDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMeson2BaryonsDecayer::incoming_,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2BaryonsDecayer,int> interfaceOutcoming1
    ("OutgoingBaryons",
     "The PDG code for the outgoing fermion",
     &VectorMeson2BaryonsDecayer::outgoingf_,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2BaryonsDecayer,int> interfaceOutcoming2
    ("OutgoingAntiBaryons",
     "The PDG code for the second outgoing anti-fermion",
     &VectorMeson2BaryonsDecayer::outgoinga_,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2BaryonsDecayer,double> interfaceGM
    ("GM",
     "The value of the GM form factor",
     &VectorMeson2BaryonsDecayer::gm_, -1, 0., -1000., 1000.,
     false, false, Interface::limited);
  
  static ParVector<VectorMeson2BaryonsDecayer,double> interfaceGE
    ("GE",
     "The value of the GE form factor",
     &VectorMeson2BaryonsDecayer::ge_, -1, 0., -1000., 1000.,
     false, false, Interface::limited);
  
  static ParVector<VectorMeson2BaryonsDecayer,double> interfacePhi
    ("Phi",
     "The phase of the GE form factor",
     &VectorMeson2BaryonsDecayer::phi_, -1, 0., -Constants::pi, Constants::pi,
     false, false, Interface::limited);

  static ParVector<VectorMeson2BaryonsDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMeson2BaryonsDecayer::maxweight_,
     0, 0, 0, -10000000, 10000000, false, false, true);
}

void VectorMeson2BaryonsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  if(outgoingf_[imode()]!=decay[iferm]->id()) swap(iferm,ianti);
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // outgoing fermion
  if(decay[iferm]->dataPtr()->iSpin()==PDT::Spin1Half)
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  else
    RSSpinorBarWaveFunction::
      constructSpinInfo(wave2bar_,decay[iferm],outgoing,true);
  // outgoing antifermion
  if(decay[ianti]->dataPtr()->iSpin()==PDT::Spin1Half)
    SpinorWaveFunction::
      constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
  else
    RSSpinorWaveFunction::
      constructSpinInfo(wave2_,decay[ianti],outgoing,true);
}

double VectorMeson2BaryonsDecayer::me2(const int,const Particle & part,
				       const tPDVector & outgoing,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  // initialze me
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,outgoing[0]->iSpin(),outgoing[0]->iSpin())));
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(outgoingf_[imode()]!=outgoing[iferm]->id()) swap(iferm,ianti);
  // initialization
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  // spin 1/2
  if(outgoing[0]->iSpin()==PDT::Spin1Half && outgoing[1]->iSpin()==PDT::Spin1Half) {
    wave_.resize(2);
    wavebar_.resize(2);
    for(unsigned int ix=0;ix<2;++ix) {
      wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[iferm],ix,Helicity::outgoing);
      wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[ianti],ix,Helicity::outgoing);
    }
    // coefficients
    Complex GM = gm_[imode()];
    Complex GE = ge_[imode()]*exp(Complex(0.,phi_[imode()]));
    LorentzPolarizationVector c2 = -2.*outgoing[0]->mass()/(4.*sqr(outgoing[0]->mass())-sqr(part.mass()))*
      (GM-GE)*(momenta[iferm]-momenta[ianti]);
    // now compute the currents
    LorentzPolarizationVector temp;
    //double mesum(0.);
    for(unsigned ix=0;ix<2;++ix) {
      for(unsigned iy=0;iy<2;++iy) {
	LorentzPolarizationVector temp = (GM*wave_[ix].vectorCurrent(wavebar_[iy])+c2*wave_[ix].scalar(wavebar_[iy]))/part.mass();
	for(unsigned int iz=0;iz<3;++iz) {
	  if(iferm>ianti) (*ME())(iz,ix,iy)=vectors_[iz].dot(temp);
	  else            (*ME())(iz,iy,ix)=vectors_[iz].dot(temp);
	  //mesum += norm(vectors_[iz].dot(temp));
	}
      }
    }
    double me = ME()->contract(rho_).real();
    // cerr << "testing decay " << part.PDGName() << " -> " << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << "\n";
    // cerr << "testing ME " << mesum/3. << " " << me << " " << 4./3.*(norm(GM)+2.*sqr(outgoing[0]->mass()/part.mass())*norm(GE)) << "\n";
    // cerr << "testing gamma " << mesum/3./8./Constants::pi*sqrt(sqr(part.mass())-4.*sqr(outgoing[0]->mass()))/MeV << "\n";
    // return the answer
    return me;
  }
  // spin 3/2
  else if(outgoing[0]->iSpin()==PDT::Spin3Half && outgoing[0]->iSpin()==PDT::Spin3Half) {
    wave2_.resize(4);
    wave2bar_.resize(4);
    wave_.resize(4);
    wavebar_.resize(4);
    RSSpinorBarWaveFunction swave(momenta[iferm],outgoing[iferm],Helicity::outgoing);
    RSSpinorWaveFunction    awave(momenta[ianti],outgoing[ianti],Helicity::outgoing);
    LorentzPolarizationVector vtemp = part.momentum()/part.mass();
    for(unsigned int ix=0;ix<4;++ix) {
      swave.reset(ix);
      awave.reset(ix);
      wave2bar_[ix] = swave.dimensionedWf();
      wavebar_ [ix] = wave2bar_[ix].dot(vtemp);
      wave2_   [ix] = awave.dimensionedWf();
      wave_    [ix] = wave2_[ix].dot(vtemp);
    }
    // coefficients
    Complex GM = gm_[imode()];
    Complex GE = ge_[imode()]*exp(Complex(0.,phi_[imode()]));
    LorentzPolarizationVector c2 = -2.*outgoing[0]->mass()/(4.*sqr(outgoing[0]->mass())-sqr(part.mass()))*
      (GM-GE)*(momenta[iferm]-momenta[ianti]);
    // now compute the currents
    for(unsigned ix=0;ix<4;++ix) {
      for(unsigned iy=0;iy<4;++iy) {
	// q(al)q(be) piece
	LorentzPolarizationVector temp2 = (GM*wave_[ix].vectorCurrent(wavebar_[iy])+c2*wave_[ix].scalar(wavebar_[iy]))*
	  2.*part.mass()/(4.*sqr(outgoing[0]->mass())-sqr(part.mass()));
	// g(al)g(be) * GM-GE piece
	LorentzPolarizationVector temp3 = wave2_[ix].generalScalar(wave2bar_[iy],1.,1.)*c2/part.mass();
	// g(al)g(be) * gamma^mu
	LorentzPolarizationVector temp1(GM/part.mass()*(wave2bar_[iy](0,3)*wave2_[ix](0,0) + wave2bar_[iy](0,2)*wave2_[ix](0,1) - 
							wave2bar_[iy](0,1)*wave2_[ix](0,2) - wave2bar_[iy](0,0)*wave2_[ix](0,3) + 
							wave2bar_[iy](1,3)*wave2_[ix](1,0) + wave2bar_[iy](1,2)*wave2_[ix](1,1) - 
							wave2bar_[iy](1,1)*wave2_[ix](1,2) - wave2bar_[iy](1,0)*wave2_[ix](1,3) + 
							wave2bar_[iy](2,3)*wave2_[ix](2,0) + wave2bar_[iy](2,2)*wave2_[ix](2,1) - 
							wave2bar_[iy](2,1)*wave2_[ix](2,2) - wave2bar_[iy](2,0)*wave2_[ix](2,3) - 
							wave2bar_[iy](3,3)*wave2_[ix](3,0) - wave2bar_[iy](3,2)*wave2_[ix](3,1) + 
							wave2bar_[iy](3,1)*wave2_[ix](3,2) + wave2bar_[iy](3,0)*wave2_[ix](3,3)),
					Complex(0,1)*GM/part.mass()*(wave2bar_[iy](0,3)*wave2_[ix](0,0) - wave2bar_[iy](0,2)*wave2_[ix](0,1) - 
								     wave2bar_[iy](0,1)*wave2_[ix](0,2) + wave2bar_[iy](0,0)*wave2_[ix](0,3) + 
								     wave2bar_[iy](1,3)*wave2_[ix](1,0) - wave2bar_[iy](1,2)*wave2_[ix](1,1) - 
								     wave2bar_[iy](1,1)*wave2_[ix](1,2) + wave2bar_[iy](1,0)*wave2_[ix](1,3) + 
								     wave2bar_[iy](2,3)*wave2_[ix](2,0) - wave2bar_[iy](2,2)*wave2_[ix](2,1) - 
								     wave2bar_[iy](2,1)*wave2_[ix](2,2) + wave2bar_[iy](2,0)*wave2_[ix](2,3) - 
								     wave2bar_[iy](3,3)*wave2_[ix](3,0) + wave2bar_[iy](3,2)*wave2_[ix](3,1) + 
								     wave2bar_[iy](3,1)*wave2_[ix](3,2) - wave2bar_[iy](3,0)*wave2_[ix](3,3)),
					GM/part.mass()*(wave2bar_[iy](0,2)*wave2_[ix](0,0) - wave2bar_[iy](0,3)*wave2_[ix](0,1) - 
							wave2bar_[iy](0,0)*wave2_[ix](0,2) + wave2bar_[iy](0,1)*wave2_[ix](0,3) + 
							wave2bar_[iy](1,2)*wave2_[ix](1,0) - wave2bar_[iy](1,3)*wave2_[ix](1,1) - 
							wave2bar_[iy](1,0)*wave2_[ix](1,2) + wave2bar_[iy](1,1)*wave2_[ix](1,3) + 
							wave2bar_[iy](2,2)*wave2_[ix](2,0) - wave2bar_[iy](2,3)*wave2_[ix](2,1) - 
							wave2bar_[iy](2,0)*wave2_[ix](2,2) + wave2bar_[iy](2,1)*wave2_[ix](2,3) - 
							wave2bar_[iy](3,2)*wave2_[ix](3,0) + wave2bar_[iy](3,3)*wave2_[ix](3,1) + 
							wave2bar_[iy](3,0)*wave2_[ix](3,2) - wave2bar_[iy](3,1)*wave2_[ix](3,3)),
					GM/part.mass()*(-wave2bar_[iy](0,2)*wave2_[ix](0,0) - wave2bar_[iy](0,3)*wave2_[ix](0,1) - 
							wave2bar_[iy](0,0)*wave2_[ix](0,2) - wave2bar_[iy](0,1)*wave2_[ix](0,3) - 
							wave2bar_[iy](1,2)*wave2_[ix](1,0) - wave2bar_[iy](1,3)*wave2_[ix](1,1) - 
							wave2bar_[iy](1,0)*wave2_[ix](1,2) - wave2bar_[iy](1,1)*wave2_[ix](1,3) - 
							wave2bar_[iy](2,2)*wave2_[ix](2,0) - wave2bar_[iy](2,3)*wave2_[ix](2,1) - 
							wave2bar_[iy](2,0)*wave2_[ix](2,2) - wave2bar_[iy](2,1)*wave2_[ix](2,3) + 
							wave2bar_[iy](3,2)*wave2_[ix](3,0) + wave2bar_[iy](3,3)*wave2_[ix](3,1) + 
							wave2bar_[iy](3,0)*wave2_[ix](3,2) + wave2bar_[iy](3,1)*wave2_[ix](3,3)));
	LorentzPolarizationVector temp = temp1+temp2+temp3;
	for(unsigned int iz=0;iz<3;++iz) {
	  if(iferm>ianti) (*ME())(iz,ix,iy)=vectors_[iz].dot(temp);
	  else            (*ME())(iz,iy,ix)=vectors_[iz].dot(temp);
	}
      }
    }
    // double mesum = ME()->contract(RhoDMatrix(PDT::Spin1)).real();
    // generator()->log() << "testing decay " << part.PDGName() << " -> " << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << "\n";
    // generator()->log() << "testing ME " << mesum  << " " << me << " " << 1./3.*(40.*norm(GM)/9.+16.*sqr(outgoing[0]->mass()/part.mass())*norm(GE)) << "\n";
    // return the answer
    return ME()->contract(rho_).real();
  }
  else
    assert(false);
}

// output the setup information for the particle database
void VectorMeson2BaryonsDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    if(ix<initsize_) {
      output << "newdef " << name() << ":Incoming " << ix << " "
	     << incoming_[ix] << "\n";
      output << "newdef " << name() << ":OutgoingFermion " << ix << " "
	     << outgoingf_[ix] << "\n";
      output << "newdef " << name() << ":OutgoingAntiFermion "  << ix << " "
	     << outgoinga_[ix] << "\n";
      output << "newdef " << name() << ":GM " << ix << " "
	     << gm_[ix] << "\n";
      output << "newdef " << name() << ":GE " << ix << " "
	     << ge_[ix] << "\n";
      output << "newdef " << name() << ":Phi " << ix << " "
	     << phi_[ix] << "\n";
      output << "newdef " << name() << ":MaxWeight " << ix << " "
	     << maxweight_[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming " << ix << " "
	     << incoming_[ix] << "\n";
      output << "insert " << name() << ":OutgoingFermion "  << ix << " "
	     << outgoingf_[ix] << "\n";
      output << "insert " << name() << ":OutgoingAntiFermion "  << ix << " "
	     << outgoinga_[ix] << "\n";
      output << "insert " << name() << ":GM " << ix << " "
	     << gm_[ix] << "\n";
      output << "insert " << name() << ":GE " << ix << " "
	     << ge_[ix] << "\n";
      output << "insert " << name() << ":Phi " << ix << " "
	     << phi_[ix] << "\n";
      output << "insert " << name() << ":MaxWeight " << ix << " "
	     << maxweight_[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
