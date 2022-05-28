// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HQETStrongDecayer class.
//

#include "HQETStrongDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/Rank3TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;

HQETStrongDecayer::HQETStrongDecayer()
  : fPi_(130.2*MeV), g_(0.566), h_(0.544), hp_(0.413), k_(0.407), kp_(0.242), gtilde_(0.283),
    psiL_(0.), psiS_(0.041), deltaEta_(1./43.7), Lambda_(1.*GeV)
{}

void HQETStrongDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size()||isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in HQETStrongDecayer"
    			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) {
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    }
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

void HQETStrongDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

IBPtr HQETStrongDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr HQETStrongDecayer::fullclone() const {
  return new_ptr(*this);
}

void HQETStrongDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(fPi_,MeV) << g_ << h_ << hp_ << k_ << kp_ << gtilde_
     << psiL_ << psiS_ << deltaEta_ << ounit(Lambda_,GeV)
     << incoming_ << outgoing_ << maxWeight_;
}

void HQETStrongDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(fPi_,MeV) >> g_ >> h_ >> hp_ >> k_ >> kp_ >> gtilde_
     >> psiL_ >> psiS_ >> deltaEta_ >> iunit(Lambda_,GeV)
     >> incoming_ >> outgoing_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<HQETStrongDecayer,DecayIntegrator>
describeHerwigHQETStrongDecayer("Herwig::HQETStrongDecayer", "HwHMDecay.so");

void HQETStrongDecayer::Init() {

  static ClassDocumentation<HQETStrongDecayer> documentation
    ("The HQETStrongDecayer class performs the strong decays of excited heavy mesons using HQET results.");
  
  static Parameter<HQETStrongDecayer,Energy> interfacefPi
    ("fPi",
     "The pion decay constant",
     &HQETStrongDecayer::fPi_, MeV, 130.2*MeV, 100.0*MeV, 200.0*MeV,
     false, false, Interface::limited);

  static Parameter<HQETStrongDecayer,double> interfaceg
    ("g",
     "The coupling for 1S (0-,1-) decays",
     &HQETStrongDecayer::g_, 0.566, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<HQETStrongDecayer,double> interfaceh
    ("h",
     "The coupling for 1P (0+,1+) decays",
     &HQETStrongDecayer::h_, 0.544, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<HQETStrongDecayer,double> interfacehp
    ("hp",
     "The coupling for 1P (1+,2+) decays",
     &HQETStrongDecayer::hp_, 0.413, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<HQETStrongDecayer,double> interfacek
    ("k",
     "The coupling for 1D (2-,3-) decays",
     &HQETStrongDecayer::k_, 0.407, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<HQETStrongDecayer,double> interfacekp
    ("kp",
     "The coupling for 1D (1-,2-) decays",
     &HQETStrongDecayer::kp_, 0.242, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<HQETStrongDecayer,double> interfacegtilde
    ("gtilde",
     "The coupling for 2S (0-,1-) decays",
     &HQETStrongDecayer::gtilde_, 0.283, 0.0, 1.0,
     false, false, Interface::limited);

  static ParVector<HQETStrongDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &HQETStrongDecayer::maxWeight_,
     0, 0, 0, 0., 100000., false, false, true);

  static Parameter<HQETStrongDecayer,double> interfaceDeltaEta
    ("DeltaEta",
     "The mixing parameter for eta-pi0 of isospin violating decays",
     &HQETStrongDecayer::deltaEta_, 1./43.7, 0.0, 1.,
     false, false, Interface::limited);

  static Parameter<HQETStrongDecayer,Energy> interfacefLambda
    ("Lambda",
     "Strong decays momentum scale",
     &HQETStrongDecayer::Lambda_, GeV, 1.*GeV, .1*GeV, 2.*GeV,
     false, false, Interface::limited);

  static Parameter<HQETStrongDecayer,double> interfacefpsiL
    ("psiL",
     "D_1 mixing angle for up and down heavy mesons",
     &HQETStrongDecayer::psiL_, 0., -M_PI/2., M_PI/2.,
     false, false, Interface::limited);

  static Parameter<HQETStrongDecayer,double> interfacefpsiS
    ("psiS",
     "D_1 mixing angle for strange heavy mesons",
     &HQETStrongDecayer::psiS_, 0.041, -M_PI/2., M_PI/2.,
     false, false, Interface::limited);

  static Command<HQETStrongDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling(1/GeV^2) and max weight for a decay",
     &HQETStrongDecayer::setUpDecayMode, false);
}

int HQETStrongDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
    	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
    	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
    	imode=ix;
    	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void HQETStrongDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  switch (part.dataPtr()->iSpin()) {
  case PDT::Spin0:
    Helicity::ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
						    Helicity::incoming,true);
    break;
  case PDT::Spin1:
    Helicity::VectorWaveFunction::constructSpinInfo(vecIn_,const_ptr_cast<tPPtr>(&part),
  						    Helicity::incoming,true,false);
    break;
  case PDT::Spin2:
    Helicity::TensorWaveFunction::constructSpinInfo(tensorIn_,const_ptr_cast<tPPtr>(&part),
						    Helicity::incoming,true,false);
    break;
  case PDT::Spin3:
    Helicity::Rank3TensorWaveFunction::constructSpinInfo(spin3In_,const_ptr_cast<tPPtr>(&part),
							 Helicity::incoming,true,false);
    break;
  default:
    assert(false);
  }
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<decay.size();++ix) {
    switch (decay[ix]->dataPtr()->iSpin()) {
    case PDT::Spin1:
      Helicity::VectorWaveFunction::constructSpinInfo(vecOut_,decay[ix],
						      Helicity::outgoing,true,false);
      break;
    case PDT::Spin0:
      Helicity::ScalarWaveFunction::constructSpinInfo(decay[ix],Helicity::outgoing,true);
      break;
    default:
      assert(false);
    }
  }
}

// matrix elememt for the process
double HQETStrongDecayer::me2(const int, const Particle & part,
  			      const tPDVector & outgoing,
  			      const vector<Lorentz5Momentum> & momenta,
  			      MEOption meopt) const {
  if(!ME()) {
    if(part.data().iSpin()!=PDT::Spin3)
      ME(new_ptr(TwoBodyDecayMatrixElement(part.data().iSpin(),
					   outgoing[0]->iSpin(),
					   outgoing[1]->iSpin())));
    else
      ME(new_ptr(GeneralDecayMatrixElement(part.data().iSpin(),
					   outgoing[0]->iSpin(),
					   outgoing[1]->iSpin())));
  }
  // stuff for incoming particle
  if(meopt==Initialize) {
    switch(part.data().iSpin()) {
    case PDT::Spin0:
      rho_ = RhoDMatrix(PDT::Spin0);
      ScalarWaveFunction::
	calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part), Helicity::incoming);
      break;
    case PDT::Spin1:
      rho_ = RhoDMatrix(PDT::Spin1);
      Helicity::VectorWaveFunction::calculateWaveFunctions(vecIn_,rho_,const_ptr_cast<tPPtr>(&part),
  							   Helicity::incoming,false);
      break;
    case PDT::Spin2:
      rho_ = RhoDMatrix(PDT::Spin2);
      Helicity::TensorWaveFunction::calculateWaveFunctions(tensorIn_,rho_,const_ptr_cast<tPPtr>(&part),
							   Helicity::incoming,false);
      break;
    case PDT::Spin3:
      rho_ = RhoDMatrix(PDT::Spin3);
      Rank3TensorWaveFunction::calculateWaveFunctions(spin3In_,rho_,const_ptr_cast<tPPtr>(&part),
						      Helicity::incoming,false);
      break;
    default:
      assert(false);
    }
  }
  // outgoing vector wavefunnctions if needed
  if(outgoing[0]->iSpin()==PDT::Spin1) {
    vecOut_={HelicityFunctions::polarizationVector(-momenta[0],0,Helicity::outgoing),
	     HelicityFunctions::polarizationVector(-momenta[0],1,Helicity::outgoing),
	     HelicityFunctions::polarizationVector(-momenta[0],2,Helicity::outgoing)};
  }
  // decay frame momentum
  Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
					     momenta[1].mass());
  // identify the incoming multiplet
  unsigned int itemp = abs(part.id())-abs(part.id())%1000;
  // calculate the matrix element
  double test(0.);
  // 0- from (0-,1-)
  if     ((itemp==0 ||itemp==100000) && part.data().iSpin()==PDT::Spin0) {
    double coup = itemp==0 ? g_ : gtilde_;
    // 0- to 1- 0-
    if(outgoing[0]->iSpin()==PDT::Spin1) {
      InvEnergy fact = -2.*coup/fPi_*sqrt(momenta[0].mass()/part.mass());
      for(unsigned int ix=0;ix<3;++ix) {
	(*ME())(0,ix,0) = Complex(fact*(vecOut_[ix]*momenta[1]));
      }
      // analytic test of the answer
      test = 4.*sqr(coup)*part.mass()/momenta[0].mass()*sqr(pcm)/sqr(fPi_);
    }
    else
      assert(false);
  }
  // 1- from (0-,1-)
  else if((itemp==0 ||itemp==100000) && part.data().iSpin()==PDT::Spin1) {
    double coup = itemp==0 ? g_ : gtilde_;
    // 1- to 0- 0-
    if(outgoing[0]->iSpin()==PDT::Spin0) {
      InvEnergy fact = -2.*coup/fPi_*sqrt(momenta[0].mass()/part.mass());
      for(unsigned int ix=0;ix<3;++ix) {
	(*ME())(ix,0,0) = Complex(fact*(vecIn_[ix]*momenta[1]));
      }
      // analytic test of the answer
      test = 4.*sqr(coup)*momenta[0].mass()*sqr(pcm)/3./sqr(fPi_)/part.mass();
    }
    // 1- to 1- 0-
    else if(outgoing[0]->iSpin()==PDT::Spin1) {
      InvEnergy2 fact = -2.*coup/fPi_*sqrt(momenta[0].mass()/part.mass())/part.mass();
      for(unsigned int ix=0;ix<3;++ix) {
	LorentzVector<complex<double> > vtemp = fact*epsilon(momenta[1],part.momentum(),vecIn_[ix]);
	for(unsigned int iy=0;iy<3;++iy) {
	  (*ME())(ix,iy,0) = vtemp*vecOut_[iy];
	}
      }
      // analytic test of the answer
      test = 8.*sqr(coup)*momenta[0].mass()*sqr(pcm)/3./sqr(fPi_)/part.mass();
    }
    else
      assert(false);
  }
  // 1+ from (1+,2+) or (0+,1+)
  else if((itemp==10000||itemp==20000)  && part.data().iSpin()==PDT::Spin1) {
    if(outgoing[0]->iSpin()==PDT::Spin1) {
      // mixing
      double psi = (abs(part.id())%100)/10!=3 ? psiL_ : psiS_;
      double m1(cos(psi)),m2(sin(psi));
      if(itemp==20000) {
	swap(m1,m2);
	m1 *=-1.;
      }
      InvEnergy2 A = -2.*sqrt(2.*momenta[0].mass()/3./part.mass())*hp_*m1/fPi_/Lambda_;
      double     B = -h_*m2/1/fPi_*sqrt(momenta[0].mass()/part.mass())/part.mass()*
	(sqr(part.mass())-sqr(momenta[0].mass())+sqr(momenta[1].mass()));
      B += A*sqr(pcm);
      // matrix element
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  (*ME())(ix,iy,0)  = B*vecOut_[iy].dot(vecIn_[ix]) + 3.*A*vecOut_[iy].dot(momenta[1])*vecIn_[ix].dot(momenta[1]);
	}
      }
      // analytic test of the answer
      test = (3.*sqr(A*part.mass()*sqr(pcm)) +sqr(B)/3.*(3.*sqr(momenta[0].mass())+sqr(pcm))
	      -A*B*sqr(pcm)*(sqr(part.mass())+sqr(momenta[0].mass())-sqr(momenta[1].mass())))/sqr(momenta[0].mass());
    }
    else
      assert(false);
  }
  // 2+ from (1+,2+)
  else if(itemp==0      && part.data().iSpin()==PDT::Spin2) {
    // 2+ -> 0- 0-
    if(outgoing[0]->iSpin()==PDT::Spin0) {
      InvEnergy2 fact = -4.*hp_/fPi_*sqrt(momenta[0].mass()/part.mass())/Lambda_;
      for(unsigned int ix=0;ix<5;++ix) {
	(*ME())(ix,0,0) = Complex(fact*((tensorIn_[ix]*momenta[1])*momenta[1]));
      }
      // analytic test of the answer
      test = 32.*sqr(hp_)*momenta[0].mass()*sqr(sqr(pcm))/15./sqr(fPi_)/sqr(Lambda_)/part.mass();
    }
    // 2+ -> 1- 0-
    else if(outgoing[0]->iSpin()==PDT::Spin1) {
      InvEnergy3 fact = -4.*hp_/fPi_*sqrt(momenta[0].mass()/part.mass())/Lambda_/part.mass();
      for(unsigned int iy=0;iy<3;++iy) {
	LorentzVector<complex<InvEnergy> > vtemp = fact*epsilon(momenta[0],vecOut_[iy],momenta[1]);
	for(unsigned int ix=0;ix<5;++ix) {
	  (*ME())(ix,iy,0) = Complex((momenta[1]*tensorIn_[ix]).dot(vtemp));
	}
      }
      // analytic test of the answer
      test = 16.*sqr(hp_)*momenta[0].mass()*sqr(sqr(pcm))/5./sqr(fPi_)/sqr(Lambda_)/part.mass();
    }
    else
      assert(false);
  }
  // 0+ from (0+,1+)
  else if(itemp==10000  && part.data().iSpin()==PDT::Spin0) {
    // 0+ -> 0- 0-
    if(outgoing[0]->iSpin()==PDT::Spin0) {
      (*ME())(0,0,0) = 2.*h_/fPi_*sqrt(momenta[0].mass()/part.mass())*(part.momentum()*momenta[1])/part.mass();
      // analytic test of the answer
      test = sqr(h_)/sqr(fPi_)*momenta[0].mass()/pow<3,1>(part.mass())*sqr(sqr(part.mass())-sqr(momenta[0].mass())+sqr(momenta[1].mass()));
    }
    else
      assert(false);
  }
  // 1- from (1-,2-)
  else if(itemp==30000  && part.data().iSpin()==PDT::Spin1) {
    // 1- -> 0- 0-
    if(outgoing[0]->iSpin()==PDT::Spin0) {
      InvEnergy fact = 4.*kp_*sqrt(2.*momenta[0].mass()/3./part.mass())/part.mass()/fPi_/Lambda_*(part.momentum()*momenta[1]);
      for(unsigned int ix=0;ix<3;++ix) {
	(*ME())(ix,0,0) = Complex(fact*(vecIn_[ix]*momenta[1]));
      }
      // analytic test of the answer
      test = 8.*momenta[0].mass()*sqr(kp_*pcm*(sqr(part.mass())-sqr(momenta[0].mass())+sqr(momenta[1].mass())))
	/9./sqr(fPi_*Lambda_*part.mass())/part.mass();
    }
    // 1- -> 1- 0-
    else if(outgoing[0]->iSpin()==PDT::Spin1) {
      InvEnergy2 fact = -2.*kp_*sqrt(2.*momenta[0].mass()/3./part.mass())/fPi_/Lambda_*(part.momentum()*momenta[1])/sqr(part.mass());	
      for(unsigned int ix=0;ix<3;++ix) {
      	LorentzVector<complex<double> > vtemp = fact*epsilon(momenta[1],part.momentum(),vecIn_[ix]);
      	for(unsigned int iy=0;iy<3;++iy) {
      	  (*ME())(ix,iy,0) = vtemp*vecOut_[iy];
      	}
      }
      // analytic test of the answer
      test = 4.*momenta[0].mass()*sqr(kp_*pcm*(sqr(part.mass())-sqr(momenta[0].mass())+sqr(momenta[1].mass())))
	/9./sqr(fPi_*Lambda_*part.mass())/part.mass();
    }
    else
      assert(false);
  }
  // 2- from (1-,2-)
  else if(itemp==10000 && part.data().iSpin()==PDT::Spin2) {
    // 2- -> 1- 0-
    if(outgoing[0]->iSpin()==PDT::Spin1) {
      InvEnergy fact = 4.*kp_*sqrt(momenta[0].mass()/part.mass())/fPi_/Lambda_*(part.momentum()*momenta[1])/part.mass();
      for(unsigned int ix=0;ix<5;++ix) {
	LorentzVector<complex<double> > vtemp = fact*(momenta[1]*tensorIn_[ix]);
	for(unsigned int iy=0;iy<3;++iy) {
	  (*ME())(ix,iy,0) = vecOut_[iy].dot(vtemp);
	}
      }
      // analytic test of the answer
      test = 2.*sqr(kp_*pcm*(sqr(part.mass())-sqr(momenta[0].mass())+sqr(momenta[1].mass())))/
	15./sqr(fPi_*Lambda_)/momenta[1].mass()/pow<5,1>(part.mass())*
	(sqr(part.mass())*(sqr(part.mass())+8.*sqr(momenta[0].mass())-2.*sqr(momenta[1].mass()))+sqr(sqr(momenta[0].mass())-sqr(momenta[1].mass())));
    }
    else
      assert(false);
  }
  // 2- from (2-,3-)
  else if(itemp==20000  && part.data().iSpin()==PDT::Spin2) {
    // 2- -> 1- 0-
    if(outgoing[0]->iSpin()==PDT::Spin1) {
      InvEnergy3 fact = 4.*k_*sqrt(momenta[0].mass()/part.mass())/sqrt(15.)/fPi_/sqr(Lambda_);
      for(unsigned int ix=0;ix<5;++ix) {
	LorentzVector<complex<Energy> > vtemp = momenta[1]*tensorIn_[ix]+tensorIn_[ix]*momenta[1];
	complex<Energy2> dot = 0.5*(vtemp*momenta[1]);
	for(unsigned int iy=0;iy<3;++iy) {
	  (*ME())(ix,iy,0) = Complex(fact*(sqr(pcm)*vecOut_[iy].dot(vtemp) +5.*dot*(vecOut_[iy]*momenta[1])));
	}
      }
      // analytic test of the answer
      test = 32.*sqr(k_)*pow<6,1>(pcm)/225./sqr(fPi_*sqr(Lambda_))/momenta[0].mass()/pow<3,1>(part.mass())*
	(sqr(part.mass())*(16.*sqr(part.mass())-2.*sqr(momenta[0].mass())+8.*sqr(momenta[1].mass()))+sqr(sqr(momenta[0].mass())-sqr(momenta[1].mass())));
    }
    else
      assert(false);
  }
  // 3- from (2-,3-)
  else if(itemp==0      && part.data().iSpin()==PDT::Spin3) {
    // 3- -> 0- 0-
    if(outgoing[0]->iSpin()==PDT::Spin0) {
      InvEnergy3 fact = -4.*k_/fPi_/sqr(Lambda_)*sqrt(momenta[0].mass()/part.mass());
      for(unsigned int ix=0;ix<7;++ix) {
	(*ME())(ix,0,0) = Complex(fact*((spin3In_[ix].dot(momenta[1],0)*momenta[1])*momenta[1]));
      }
      // analytic test of answer
      test = 32.*sqr(k_)*momenta[0].mass()/35./part.mass()*pow<6,1>(pcm)/sqr(fPi_*sqr(Lambda_));
    }
    // 3- -> 1- 0-
    else if(outgoing[0]->iSpin()==PDT::Spin1) {
      InvEnergy4 fact = 4.*k_/fPi_/sqr(Lambda_)*sqrt(momenta[0].mass()/part.mass())/part.mass();
      for(unsigned int ix=0;ix<7;++ix) {
	LorentzVector<complex<InvEnergy2> > vtemp = fact*(spin3In_[ix].dot(momenta[1],0)*momenta[1]);
	for(unsigned int iy=0;iy<3;++iy) {
	  LorentzVector<complex<Energy2> > vtemp2 =epsilon(momenta[1],part.momentum(),vecOut_[iy]);
	  (*ME())(ix,iy,0) = Complex(vtemp*vtemp2); 
	}
      }
      // analytic test of answer
      test = 128.*sqr(k_)*momenta[0].mass()/part.mass()*pow<6,1>(pcm)/105./sqr(fPi_*sqr(Lambda_));
    }
    else
      assert(false);
  }
  else
    assert(false);
  // spin average
  double output = ME()->contract(rho_).real();
  // testing
  double ratio = abs(output-test)/(output+test);
  if(ratio>1e-14 && generator()->state()!=InterfacedBase::runready)
    generator()->log() << "testing matrix element for " << part.PDGName() << " -> "
		       << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
		       << output << " " << test << " " << ratio << endl;
  // isospin factors
  if(abs(outgoing[1]->id())==ParticleID::pi0) {
    int ispect = (abs(part.id())%100)/10;
    output *= ispect<3 ? 0.5 : 0.125*sqr(deltaEta_);
  }
  // return the answer
  return output;
}

bool HQETStrongDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
				      double & coupling) const {
  // extract ids of the particles
  int id(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  pair<tPDPtr,tPDPtr> out;
  out.first=*pit;
  ++pit;
  out.second=*pit;
  bool order=true;
  if (abs(out.first->id())<abs(out.second->id())) {
    order=false;
    swap(out.first,out.second);
  }
  // identify the incoming multiplet
  unsigned int itemp = abs(id)-abs(id)%1000;
  coupling=1.;
  if     (itemp==0 && dm.parent()->iSpin()==PDT::Spin0) {
    assert(out.first->iSpin()!=PDT::Spin0);
    mecode = 101;
  }
  else if(itemp==0 && dm.parent()->iSpin()==PDT::Spin1) {
    mecode = out.first->iSpin()==PDT::Spin0 ? 102 : 103;
  }
  else if((itemp==10000||itemp==20000)  && dm.parent()->iSpin()==PDT::Spin1) {
    assert(out.first->iSpin()!=PDT::Spin0);
    mecode = 104;
  }
  else if(itemp==0      && dm.parent()->iSpin()==PDT::Spin2) {
    mecode = out.first->iSpin()==PDT::Spin0 ? 105 : 106;
  }
  else if(itemp==10000  && dm.parent()->iSpin()==PDT::Spin0) {
    mecode = 107;
  }
  else if(itemp==30000  && dm.parent()->iSpin()==PDT::Spin1) {
    mecode = out.first->iSpin()==PDT::Spin0 ? 108 : 109;
  }
  else if(itemp==10000 && dm.parent()->iSpin()==PDT::Spin2) {
    assert(out.first->iSpin()!=PDT::Spin0);
    mecode = 110;
  }
  else if(itemp==20000  && dm.parent()->iSpin()==PDT::Spin2) {
    assert(out.first->iSpin()!=PDT::Spin0);
    mecode = 111;
  }
  else if(itemp==0      && dm.parent()->iSpin()==PDT::Spin3) {
    mecode = out.first->iSpin()==PDT::Spin0 ? 112 : 113;
  }
  else if(itemp==100000 && dm.parent()->iSpin()==PDT::Spin0) {
    assert(out.first->iSpin()!=PDT::Spin0);
    mecode = 114;
  }
  else if(itemp==100000 && dm.parent()->iSpin()==PDT::Spin1) {
    mecode = out.first->iSpin()==PDT::Spin0 ? 115 : 116;
  }
  // isospin factors
  if(abs(out.second->id())==ParticleID::pi0) {
    int ispect = (abs(dm.parent()->id())%100)/10;
    coupling *= sqrt(ispect<3 ? 0.5 : 0.125*sqr(deltaEta_));
  }
  // return the order
  return order;
}

void HQETStrongDecayer::dataBaseOutput(ofstream & output,
  				       bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  // couplings
  output << "newdef " << name() << ":fPi    " << fPi_/MeV    << "\n";
  output << "newdef " << name() << ":g      " << g_          << "\n";
  output << "newdef " << name() << ":h      " << h_          << "\n";
  output << "newdef " << name() << ":hp     " << hp_         << "\n";
  output << "newdef " << name() << ":k      " << k_          << "\n";
  output << "newdef " << name() << ":kp     " << kp_         << "\n";
  output << "newdef " << name() << ":gtilde " << gtilde_     << "\n";
  output << "newdef " << name() << ":Lambda " << Lambda_/GeV << "\n";
  output << "newdef " << name() << ":psiL   " << psiL_       << "\n";
  output << "newdef " << name() << ":psiS   " << psiS_       << "\n";
  output << "newdef " << name() << ":DeltaEta " << deltaEta_ << "\n";
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode "
	   << incoming_[ix] << " " << outgoing_[ix].first << " "
	   << outgoing_[ix].second << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\""
  		    << fullName() << "\";" << endl;
}

string HQETStrongDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(int(pData->iSpin())%2!=1)
    return "Incoming particle with id " + std::to_string(in) + "does not integer spin";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0 && pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 0/1";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
