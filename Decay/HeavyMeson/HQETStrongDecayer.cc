// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HQETStrongDecayer class.
//

#include "HQETStrongDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;

HQETStrongDecayer::HQETStrongDecayer()
  : fPi_(130.2*MeV), g_(0.565), h_(0.565), f_(0.565/3.), deltaEta_(1./43.7), Lambda_(1.*GeV),
    incoming_ ({413,413,423,433,                       //D*   decay modes: VtoSS
                415,415,425,425,435,435,               //D*_2 decay modes: TtoSS
                415,415,425,425,435,435,               //D*_2 decay modes: TtoVS
                10413,10413,10423,10423,10433,10433,   //D_1  decay modes: VtoVS
                10411,10411,10421,10421,10431}),       //D*_0 decay modes: StoSS
    outgoingH_({421,411,421,431,
                411,421,411,421,411,421,
                413,423,413,423,413,423,
                413,423,413,423,413,423,
                411,421,411,421,431}),
    outgoingL_({211,111,111,111,
                111,211,-211,111,311,321,
                111,211,-211,111,311,321,
                111,211,-211,111,311,321,
                111,211,-211,111,111}),
    type_     ({1,  1,  1, -1,
                2,  2,  2,  2,  2,  2,
                3,  3,  3,  3,  3,  3,
                4,  4,  4,  4,  4,  4,
                5,  5,  5,  5,  5}),
    maxWeight_({1., 1., 1., 1.,
                1., 1., 1., 1., 1., 1.,
                1., 1., 1., 1., 1., 1.,
                1., 1., 1., 1., 1., 1.,
                1., 1., 1., 1., 1.})
{}

  void HQETStrongDecayer::doinit() {
    DecayIntegrator::doinit();
    // check consistence of the parameters
    unsigned int isize=incoming_.size();
    if(isize!=outgoingH_.size()||isize!=outgoingL_.size()||
       isize!=maxWeight_.size()||isize!=type_     .size())
      throw InitException() << "Inconsistent parameters in HQETStrongDecayer"
    			  << Exception::abortnow;
    // set up the integration channels
    PhaseSpaceModePtr mode;
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
       tPDPtr     in = getParticleData(incoming_[ix]);
       tPDVector out = {getParticleData(outgoingH_[ix]),
    		      getParticleData(outgoingL_[ix])};
      if(in&&out[0]&&out[1])
        mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
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
    os << ounit(fPi_,MeV) << g_ << h_ << f_<< deltaEta_ << ounit(Lambda_,GeV) << maxWeight_;
  }

  void HQETStrongDecayer::persistentInput(PersistentIStream & is, int) {
    is >> iunit(fPi_,MeV) >> g_ >> h_ >> f_ >> deltaEta_ >> iunit(Lambda_,GeV) >> maxWeight_;
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
       "The coupling for D* decays",
       &HQETStrongDecayer::g_, 0.565, 0.0, 1.0,
       false, false, Interface::limited);

    static Parameter<HQETStrongDecayer,double> interfaceh
      ("h",
      "The coupling for D_1, D*_2 and D*_2s decays",
      &HQETStrongDecayer::h_, 0.565, 0.0, 1.0,
      false, false, Interface::limited);

    static Parameter<HQETStrongDecayer,double> interfacef
      ("f",
      "The coupling for D*_0 and D'_1 decays",
      &HQETStrongDecayer::f_, 0.565/3., 0.0, 1.0,
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
        if((id1   ==outgoingH_[ix]&&id2   ==outgoingL_[ix])||
    	 (id2   ==outgoingH_[ix]&&id1   ==outgoingL_[ix])) imode=ix;
      }
      if(idbar==incoming_[ix]) {
        if((id1bar==outgoingH_[ix]&&id2bar==outgoingL_[ix])||
    	 (id2bar==outgoingH_[ix]&&id1bar==outgoingL_[ix])) {
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
      if(abs(type_[imode()])==1) {
        ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0)));
      }
      if(abs(type_[imode()])==2) {
        ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin0,PDT::Spin0)));
      }
      if(abs(type_[imode()])==3) {
        ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin0)));
      }
      if(abs(type_[imode()])==4) {
        ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
      }
      if(abs(type_[imode()])==5) {
        ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0)));
      }
    }
    // stuff for incoming particle
    if(meopt==Initialize) {
      if(abs(type_[imode()])==1 || abs(type_[imode()])==4) {
        rho_ = RhoDMatrix(PDT::Spin1);
        Helicity::VectorWaveFunction::calculateWaveFunctions(vecIn_,rho_,const_ptr_cast<tPPtr>(&part),
  							   Helicity::incoming,false);
      }
      else if(abs(type_[imode()])==2 || abs(type_[imode()])==3) {
        rho_ = RhoDMatrix(PDT::Spin2);
        Helicity::TensorWaveFunction::calculateWaveFunctions(tensorIn_,rho_,const_ptr_cast<tPPtr>(&part),
                   Helicity::incoming,false);
      }
      else if(abs(type_[imode()])==5) {
        rho_ = RhoDMatrix(PDT::Spin0);
        Helicity::ScalarWaveFunction::calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),
                   Helicity::incoming);
      }
      else {
        cerr << "Unknown decay mode: type " << type_[imode()] << " for " << part << "\n";
        assert(false);
      }
    }
    // calculate the matrix element
    Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
					       momenta[1].mass()); //test subject
    
    double test(0.);
    // HeavyVectorMeson to PScalarMeson + PScalarMeson
    if(abs(type_[imode()])==1) {
      InvEnergy fact = -2.*g_/fPi_*sqrt(momenta[0].mass()/part.mass());
      for(unsigned int ix=0;ix<3;++ix) {
        (*ME())(ix,0,0) = fact*(vecIn_[ix]*momenta[1]);
      }
      // analytic test of the answer
      test = 4.*sqr(g_)*momenta[0].mass()*sqr(pcm)/3./sqr(fPi_)/part.mass();
    }
    // HeavyTensorMeson to PScalarMeson + PScalarMeson
    else if(abs(type_[imode()])==2) {
      InvEnergy2 fact = -2.*h_/fPi_*sqrt(momenta[0].mass()/part.mass())/Lambda_;
      for(unsigned int ix=0;ix<5;++ix) {
        (*ME())(ix,0,0) = fact*(tensorIn_[ix]*momenta[1]*momenta[0]);
      }
      // analytic test of the answer
      test = 8.*sqr(h_)*momenta[0].mass()*sqr(sqr(pcm))/15./sqr(fPi_)/sqr(Lambda_)/part.mass();
    }
    // HeavyTensorMeson to VectorMeson + PScalarMeson
    else if(abs(type_[imode()])==3) {
      // get the polarization vectors
      vecOut_={
	HelicityFunctions::polarizationVector(-momenta[0],0,Helicity::outgoing),
	HelicityFunctions::polarizationVector(-momenta[0],1,Helicity::outgoing),
	HelicityFunctions::polarizationVector(-momenta[0],2,Helicity::outgoing)};
      InvEnergy3 fact = -2.*h_/fPi_*sqrt(momenta[0].mass()/part.mass())/Lambda_/part.mass();
      for(unsigned int ix=0;ix<5;++ix) {
        for(unsigned int iy=0;iy<3;++iy) {
          LorentzVector<complex<InvEnergy> > vtemp =
	    fact*epsilon(momenta[0],vecOut_[iy],momenta[1]);
          (*ME())(ix,iy,0) = (momenta[1]*tensorIn_[ix]).dot(vtemp);
        }
      }
      // analytic test of the answer
      test = 4.*sqr(h_)*momenta[0].mass()*sqr(sqr(pcm))/5./sqr(fPi_)/sqr(Lambda_)/part.mass();
    }
    // PVectorMeson to VectorMeson + PScalarMeson
    else if(abs(type_[imode()])==4) {
      // get the polarization vectors
      vecOut_={
        HelicityFunctions::polarizationVector(-momenta[0],0,Helicity::outgoing),
        HelicityFunctions::polarizationVector(-momenta[0],1,Helicity::outgoing),
        HelicityFunctions::polarizationVector(-momenta[0],2,Helicity::outgoing)};
      InvEnergy2 fact = sqrt(2./3.)*(h_/fPi_)*sqrt(momenta[0].mass()/part.mass())/Lambda_*momenta[0].mass()/part.mass();
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  (*ME())(ix,iy,0)=Complex(fact*(vecOut_[iy].dot(vecIn_[ix])
					 *(momenta[1].mass2()-sqr(part.momentum()*momenta[1]/part.mass()))
					 - 3.*(vecIn_[ix]*momenta[1])*(vecOut_[iy]*momenta[1])));
	}
      }
      // analytic test of the answer
      test = 4.*sqr(h_)*momenta[0].mass()*sqr(sqr(pcm))/3./sqr(fPi_)/sqr(Lambda_)/part.mass()*
	(25.-2.*sqr(momenta[0].mass()/part.mass())+10.*sqr(momenta[1].mass()/part.mass())+sqr((sqr(momenta[0].mass())-sqr(momenta[1].mass()))/sqr(part.mass())))/24.;
    }
    // ScalarMeson to ScalarMeson + ScalarMeson
    else if(abs(type_[imode()])==5) {
      InvEnergy fact = f_/fPi_*sqrt(momenta[0].mass()/part.mass());
      (*ME())(0,0,0) = fact*(part.momentum()*momenta[0]/part.mass()
			     - part.momentum()*momenta[1]/momenta[0].mass());
      // analytic test of the answer
      test = sqr(f_)/(4.*sqr(fPi_))*(momenta[0].mass()/part.mass())
	* sqr(part.mass()-momenta[0].mass())
	* sqr(part.mass()+momenta[0].mass()-momenta[1].mass())
	* sqr(part.mass()+momenta[0].mass()+momenta[1].mass())
	/ sqr(part.mass()*momenta[0].mass());
    }
    else {
      assert(false);
    }
    double output = ME()->contract(rho_).real();
    // testing
    double ratio = (output-test)/(output+test);
    generator()->log() << "testing matrix element for " << part.PDGName() << " -> "
		       << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
		       << output << " " << test << " " << ratio << endl;
    // isospin factors
    if(abs(outgoing[1]->id())==ParticleID::pi0) {
      output *= type_[imode()]>0 ? 0.5 : 0.125*sqr(deltaEta_);
    }
    // return the answer
    return output;
  }

  bool HQETStrongDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
  					       double & coupling) const {
    // int imode(-1);
    // int id(dm.parent()->id());
    // int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
    // ParticleMSet::const_iterator pit(dm.products().begin());
    // int id1((**pit).id());
    // int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
    // ++pit;
    // int id2((**pit).id());
    // int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
    // unsigned int ix(0); bool order(false);
    // do {
    //   if(id   ==incoming_[ix]) {
    //     if(id1==outgoingH_[ix]&&id2==outgoingL_[ix]) {
    // 	imode=ix;
    // 	order=true;
    //     }
    //     if(id2==outgoingH_[ix]&&id1==outgoingL_[ix]) {
    // 	imode=ix;
    // 	order=false;
    //     }
    //   }
    //   if(idbar==incoming_[ix]&&imode<0) {
    //     if(id1bar==outgoingH_[ix]&&id2bar==outgoingL_[ix]) {
    // 	imode=ix;
    // 	order=true;
    //     }
    //     if(id2bar==outgoingH_[ix]&&id1bar==outgoingL_[ix]) {
    // 	imode=ix;
    // 	order=false;
    //     }
    //   }
    //   ++ix;
    // }
    // while(ix<incoming_.size()&&imode<0);
    // coupling=_coupling[imode]*dm.parent()->mass();
    // mecode=7;
    // return order;
  }

  void HQETStrongDecayer::dataBaseOutput(ofstream & output,
  				       bool header) const {
    if(header) output << "update decayers set parameters=\"";
    // parameters for the DecayIntegrator base class
    DecayIntegrator::dataBaseOutput(output,false);
    // the rest of the parameters
    // couplings
    output << "newdef " << name() << ":fPi " << fPi_/MeV << "\n";
    output << "newdef " << name() << ":g   " << g_      << "\n";
    if(header) output << "\n\" where BINARY ThePEGName=\""
  		    << fullName() << "\";" << endl;
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      output << "newdef " << name() << ":MaxWeight " << ix << " " << maxWeight_[ix] << "\n";
    }
  }
