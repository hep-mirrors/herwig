// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2BaryonsDecayer class.
//

#include "VectorMeson2BaryonsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
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
  if(isize!=incoming_.size()  || isize!=outgoing_.size()  ||
     isize!= maxweight_.size()||
     isize!= phi_.size() || isize!= ge_.size())
    throw InitException() << "Inconsistent parameters in VectorMeson2"
			   << "BaryonsDecayer::doiin() " << Exception::runerror;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

//   : gm_       ({0.00163222616377,0.00158881446341,0.00163819673075,0.000944247071286,0.00141162244048 ,0.00150424724997 ,0.00093350341391,0.00107528142563,0.000860759359361,0.00103222902272,0.00098818310509,0.00119542938517 ,0.000941856299908,0.00108249234586 ,0.00102912698738 ,0.000561082932891,0.000526423800011,0.000500391020504,0.000678994938986,0.00102277134344,0.00288222568159}),
//     ge_       ({0.00135736265293,0.0015117595557 ,0.00136696938558,0.001988067505   ,0.000852659560453,0.000801719253474,0.00150640470181,0.00153402258271,0.0015323974485  ,0.              ,0.00084599445285,0.000621041230885,0.000599393812308,0.000327664954147,0.000664382626732,0.000260326162249,0.000362882855521,0.00025226758188 ,0.000397805035356,0.00233118573611,0.00194199600339}),
//     phi_      ({0.              ,0.              ,0.740           ,0.               ,0.               ,0.               ,0.              ,0.              ,0.               ,0.              ,0.              ,0.               ,0.               ,0.               ,0.               ,0.               ,0.               ,0.               ,0.               ,-0.270          ,0.379           }),
//     incoming_ ({443             ,443             ,443             ,443              ,443              ,443              ,443             ,443             ,443              ,100443          ,100443          ,100443           ,100443           ,100443           ,100443           ,100443           ,100443           ,100443           ,100443           ,443             ,100443          }),
//     outgoingf_({ 2212           , 2112           , 3122           , 3212            , 3312            , 3322            , 3114           , 3224           , 3214            , 2212           , 2112           , 3122            , 3212            , 3312            , 3322            , 3114            , 3224            , 3214            , 3314            , 3222           , 3222           }),
//     outgoinga_({-2212           ,-2112           ,-3122           ,-3212            ,-3312            ,-3322            ,-3114           ,-3224           ,-3214            ,-2212           ,-2112           ,-3122            ,-3212            ,-3312            ,-3322            ,-3114            ,-3224            ,-3214            ,-3314            ,-3222           ,-3222           }),
//     maxweight_({1.6             ,1.6             ,2.1             ,1.5              ,2.               ,2.1              ,5.              ,7.              ,5.5              ,1.7             ,1.7             ,2.5              ,1.6              ,2.               ,2.               ,5.               ,5.               ,4.5              , 0.0012          ,2.0             ,15.             })
// {}

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
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(incoming_[ix]==idbar) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
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
  os << ge_ << gm_ << phi_ << incoming_ << outgoing_ << maxweight_;
}

void VectorMeson2BaryonsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> ge_ >> gm_ >> phi_ >> incoming_ >> outgoing_ >> maxweight_;
}


DescribeClass<VectorMeson2BaryonsDecayer,DecayIntegrator>
describeHerwigVectorMeson2BaryonsDecayer("Herwig::VectorMeson2BaryonsDecayer", "HwVMDecay.so");

void VectorMeson2BaryonsDecayer::Init() {

  static ClassDocumentation<VectorMeson2BaryonsDecayer> documentation
    ("The VectorMeson2BaryonsDecayer class is designed for "
     "the decay of vector mesons to baryons.");
  
  static Command<VectorMeson2BaryonsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, baryon, antibaryon), GE, GM, phase and max weight for a decay",
     &VectorMeson2BaryonsDecayer::setUpDecayMode, false);
}

void VectorMeson2BaryonsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=decay[iferm]->id()) swap(iferm,ianti);
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
  if(outgoing_[imode()].first!=outgoing[iferm]->id()) swap(iferm,ianti);
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
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << gm_[ix] << " " << ge_[ix] << " " << phi_[ix] << " " << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string VectorMeson2BaryonsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  // if(pData->iSpin()!=PDT::Spin1Half)
  //   return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  // if(pData->iSpin()!=PDT::Spin1Half)
  //   return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1/2";
  // get the couplings
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  ge_.push_back(stof(stype));
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  gm_.push_back(stof(stype));
  // and phase
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  phi_.push_back(stof(stype));
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  maxweight_.push_back(wgt);
  // success
  return "";
}
