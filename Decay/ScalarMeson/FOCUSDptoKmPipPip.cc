// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DtoKPiPiFOCUS class.
//

#include "DtoKPiPiFOCUS.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DtoKPiPiFOCUS::DtoKPiPiFOCUS() : rD0_(5./GeV), rRes_(1.5/GeV),
				 g1_(0.31072*GeV), g2_(-0.02323*GeV),
				 beta_(3.389*GeV), theta_(286), gamma_({304,126,211}),
				 c1_({1.655,0.780,-0.954}), c2_(17.182), c3_(0.734),
				 maxWgt_(1.),weights_({0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1}),
				 mRes_({0.896 *GeV, 1.414*GeV, 1.717*GeV, 1.4324*GeV}),
				 wRes_({0.0503*GeV, 0.232*GeV, 0.322*GeV, 0.109 *GeV})
{}

IBPtr DtoKPiPiFOCUS::clone() const {
  return new_ptr(*this);
}

IBPtr DtoKPiPiFOCUS::fullclone() const {
  return new_ptr(*this);
}

void DtoKPiPiFOCUS::persistentOutput(PersistentOStream & os) const {
  os << ounit(beta_,GeV) << theta_ << gamma_
     << ounit(g1_,GeV) << ounit(g2_,GeV) << c1_ << c2_ << c3_
     << KHalf_ << KThreeHalf_ << maxWgt_ << weights_
     << ounit(rD0_,1./GeV) << ounit(rRes_,1./GeV)
     << ounit(mRes_,GeV) << ounit(wRes_,GeV);
}

void DtoKPiPiFOCUS::persistentInput(PersistentIStream & is, int) {
  is >> iunit(beta_,GeV) >> theta_ >> gamma_
     >> iunit(g1_,GeV) >> iunit(g2_,GeV) >> c1_ >> c2_ >> c3_
     >> KHalf_ >> KThreeHalf_ >> maxWgt_ >> weights_
     >> iunit(rD0_,1./GeV) >> iunit(rRes_,1./GeV)
     >> iunit(mRes_,GeV) >> iunit(wRes_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DtoKPiPiFOCUS,DecayIntegrator>
describeHerwigDtoKPiPiFOCUS("Herwig::DtoKPiPiFOCUS", "HwSMDecay.so");

void DtoKPiPiFOCUS::Init() {

  static ClassDocumentation<DtoKPiPiFOCUS> documentation
    ("There is no documentation for the DtoKPiPiFOCUS class");
  
  static Reference<DtoKPiPiFOCUS,KMatrix> interfaceKHalf
    ("KHalf",
     "The K-matrix for the I=1/2 s-wave component",
     &DtoKPiPiFOCUS::KHalf_, false, false, true, false, false);
  
  static Reference<DtoKPiPiFOCUS,KMatrix> interfaceKThreeHalf
    ("KThreeHalf",
     "The K-matrix for the I=1/2 s-wave component",
     &DtoKPiPiFOCUS::KThreeHalf_, false, false, true, false, false);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceg1
    ("g1",
     "The coupling of the first channel to K*0",
     &DtoKPiPiFOCUS::g1_, 0.31072*GeV, -10.*GeV, 10.*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceg2
    ("g2",
     "The coupling of the first channel to K*0",
     &DtoKPiPiFOCUS::g2_, -0.02323*GeV, -10.*GeV, 10.*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfacebeta
    ("beta",
     "The beta coefficient for the P-vector",
     &DtoKPiPiFOCUS::beta_, 3.389*GeV, -10.*GeV, 10.*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacetheta
    ("theta",
     "The theta phase (in degrees) for the P-vector",
     &DtoKPiPiFOCUS::theta_, 286, 0.0, 360.,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiFOCUS,double> interfacegamma
    ("gamma",
     "The gamma phases in degress",
     &DtoKPiPiFOCUS::gamma_, 3, 0., 0.0, 360.0,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiFOCUS,double> interfacec1
    ("c1",
     "The c1 coefficients for the P-vector",
     &DtoKPiPiFOCUS::c1_, -1, 0., -100., 100.,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiFOCUS,double> interfacec2
    ("c2",
     "The c2 coefficients for the P-vector",
     &DtoKPiPiFOCUS::c2_, -1, 0., -100., 100.,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiFOCUS,double> interfacec3
    ("c3",
     "The c3 coefficients for the P-vector",
     &DtoKPiPiFOCUS::c3_, -1, 0., -100., 100.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,InvEnergy> interfaceDRadius
    ("DRadius",
     "The radius parameter for the Blatt-Weisskopf form-factor for the D",
     &DtoKPiPiFOCUS::rD0_, 1./GeV, 5./GeV, ZERO, 10./GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,InvEnergy> interfaceResonanceRadius
    ("ResonanceRadius",
     "The radius parameter for the Blatt-Weisskopf form-factor for the"
     "intermediate resonances",
     &DtoKPiPiFOCUS::rRes_, 1./GeV, 1.5/GeV, ZERO, 10./GeV,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiFOCUS,Energy> interfaceMasses
    ("Masses",
     "The masses of the resonances",
     &DtoKPiPiFOCUS::mRes_, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static ParVector<DtoKPiPiFOCUS,Energy> interfaceWidths
    ("Widths",
     "The widths of the resonances",
     &DtoKPiPiFOCUS::wRes_, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

int DtoKPiPiFOCUS::modeNumber(bool & cc,tcPDPtr parent,
			      const tPDVector & children) const {
  int id0(parent->id());
  // incoming particle must be D+
  if(abs(id0)!=ParticleID::Dplus) return -1;
  cc = id0<0;
  // must be three decay products
  if(children.size()!=3) return -1;
  unsigned int npip(0),npim(0),nkm(0),nk0(0),npi0(0);
  int id;
  for(tPDPtr child : children) {
    id=child->id();
    if(id          ==ParticleID::piplus)  ++npip;
    else if(id     ==ParticleID::pi0)     ++npi0;
    else if(id     ==ParticleID::piminus) ++npim;
    else if(abs(id)==ParticleID::K0)      ++nk0;
    else if(id     ==ParticleID::K_L0)    ++nk0;
    else if(id     ==ParticleID::K_S0)    ++nk0;
    else if(abs(id)==ParticleID::Kplus)   ++nkm;
  }
  if((id0==ParticleID::Dplus &&(nkm==1&&npip==2))||
     (id0==ParticleID::Dminus&&(nkm==1&&npim==2))) return 0;
  else return -1;
}

void DtoKPiPiFOCUS::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double DtoKPiPiFOCUS::me2(const int ichan, const Particle & part,
			  const tPDVector & ,
			  const vector<Lorentz5Momentum> & momenta,
			  MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  Complex amp(0.);
  for(unsigned int ipi=1;ipi<3;++ipi) {
    unsigned int iother = ipi==1 ? 2 : 1;
    Lorentz5Momentum pres=momenta[0]+momenta[ipi];
    pres.rescaleMass();
    // momentum for the D Blatt-Weisskopf factor
    Energy pD = Kinematics::pstarTwoBodyDecay(part.mass(),pres.mass(),momenta[iother].mass());
    // momentum for the resonance Blatt-Weisskopf factor
    Energy pAB = Kinematics::pstarTwoBodyDecay(pres.mass(),momenta[0].mass(),momenta[ipi].mass());
    // loop over the resonances
    double rD = rD0_*pD, rR = rRes_*pAB;
    // most resonaces vectors so calculate that now
    double Fd = 1./sqrt(1.+sqr(rD));
    double Fr = 1./sqrt(1.+sqr(rR));
    unsigned int ipow=3;
    for(unsigned int ires=0;ires<4;++ires) {
      if(ichan>=0 && ichan+1!=int(ires*2+ipi)) continue;
      Energy pR = Kinematics::pstarTwoBodyDecay(mRes_[ires],momenta[0].mass(),momenta[ipi].mass());
      // last resonance is the K_2
      if(ires==3) {
	Fd = 1./sqrt(9.+3.*sqr(rD)+Math::powi(rD,4));
	Fr = 1./sqrt(9.+3.*sqr(rR)+Math::powi(rR,4));
	ipow = 5;
      }
      double rB = rRes_*pR;
      double Fb = ires<3 ? 1./sqrt(1.+sqr(rB)) : 1./sqrt(9.+3.*sqr(rB)+Math::powi(rB,4));
      Energy gamma = wRes_[ires]*mRes_[ires]/pres.mass()*sqr(Fr/Fb)*Math::powi(pAB/pR,ipow);
      complex<InvEnergy2> BW = Fd*Fr/(sqr(mRes_[ires])-pres.mass2()-Complex(0.,1.)*gamma*mRes_[ires]);
      // // and the decay momenta
      // 
      // Energy pR  = Kinematics::pstarTwoBodyDecay(mres,mA,mB);
      // double Fd(1.),Fr(1.),s(1.);
      // switch(ispin) {
      // case 0:
      //   // default values of parameters are correct
      //   break;
      // case 1:
      //   Fr = sqrt((1.+sqr(_rres*pR ))/(1.+sqr(_rres*pAB )));
      //   Fd = 
      //   s = ((mAC-mBC)*(mAC+mBC)+(mD-mC)*(mD+mC)*(mB-mA)*(mB+mA)/sqr(mres))/GeV2;
      //   break;
      // case 2:
      //   Fr = sqrt((9.+3.*sqr(_rres*pR  )+Math::powi(_rres*pR  ,4))/
      // 	      (9.+3.*sqr(_rres*pAB )+Math::powi(_rres*pAB ,4)));
      //   Fd = 
      //   s = sqr(((mBC-mAC)*(mBC+mAC)+(mD-mC)*(mD+mC)*(mA-mB)*(mA+mB)/sqr(mres))/GeV2)
      //     -(mAB*mAB-2.*mD*mD-2.*mC*mC+sqr((mD-mC)*(mD+mC))/sqr(mres))*
      //      (mAB*mAB-2.*mA*mA-2.*mB*mB+sqr((mA-mB)*(mA+mB))/sqr(mres))/3./GeV2/GeV2;
      //   break;
      // default:
      //   throw Exception() << "D0toKPiPiCLEO::amplitude spin is too high ispin = " 
      // 		      << ispin << Exception::runerror;
      // }
      // }
    }
    if(ichan>=0 && ichan!=int(7+ipi)) continue;
    // finally the scalar piece
    static const Energy2 sc=2.*GeV2;
    ublas::vector<Complex> pVector(2);
    double fact = 1.-pres.mass2()/KHalf_->poles()[0];
    Complex f1=exp(Complex(0.,1.)*gamma_[0]/180.*Constants::pi)*fact;
    double shat = (pres.mass2()-sc)/GeV2;
    pVector[0]=0.;
    for(unsigned int ix=0;ix<c1_.size();++ix) {
      pVector[0] += f1*c1_[ix];
      f1 *= shat;
    }
    pVector[1]=0.;
    Complex f2=exp(Complex(0.,1.)*gamma_[1]/180.*Constants::pi)*fact;
    for(unsigned int ix=0;ix<c2_.size();++ix) {
      pVector[1] += f2*c2_[ix];
      f2 *= shat;
    }
    ublas::vector<Complex> p3Vector(1);
    p3Vector[0]=0.;
    Complex f3=exp(Complex(0.,1.)*gamma_[2]/180.*Constants::pi);
    for(unsigned int ix=0;ix<c3_.size();++ix) {
      p3Vector[0] += f2*c3_[ix];
      f3 *= shat;
    }
    Complex phase = exp(Complex(0.,1.)*theta_/180.*Constants::pi);
    pVector[0] += Complex(beta_*phase*g1_/KHalf_->poles()[0]);
    pVector[1] += beta_*phase*g2_/KHalf_->poles()[0];
    amp +=      KHalf_->amplitudes(pres.mass2(), pVector,true)[0];
    amp += KThreeHalf_->amplitudes(pres.mass2(),p3Vector,true)[0];
  }
  // now compute the matrix element
  (*ME())(0,0,0,0)=amp;
  return norm(amp);
}

void DtoKPiPiFOCUS::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // // parameters
  // output << "newdef " << name() << ":KmPipPipNonResonantMagnitude " 
  // 	 << _a1NR      << "\n";
  // output << "newdef " << name() << ":KmPipPipNonResonantPhase     " 
  // 	 << _phi1NR    << "\n";
  // output << "newdef " << name() << ":KmPipPipK892Magnitude        " 
  // 	 << _a1K892    << "\n";
  // output << "newdef " << name() << ":KmPipPipK892Phase            " 
  // 	 << _phi1K892  << "\n";
  // output << "newdef " << name() << ":KmPipPipK1430Magnitude       " 
  // 	 << _a1K1430   << "\n";
  // output << "newdef " << name() << ":KmPipPipK1430Phase           " 
  // 	 << _phi1K1430 << "\n";
  // output << "newdef " << name() << ":KmPipPipK1680Magnitude       " 
  // 	 << _a1K1680   << "\n";
  // output << "newdef " << name() << ":KmPipPipK1680Phase           " 
  // 	 << _phi1K1680 << "\n";
  // output << "newdef " << name() << ":KmPipPi0NonResonantMagnitude " 
  // 	 << _a2NR      << "\n";
  // output << "newdef " << name() << ":KmPipPi0NonResonantPhase     " 
  // 	 << _phi2NR    << "\n";
  // output << "newdef " << name() << ":KmPipPi0K8920Magnitude       " 
  // 	 << _a2K8920   << "\n";
  // output << "newdef " << name() << ":KmPipPi0K8920Phase           " 
  // 	 << _phi2K8920 << "\n";
  // output << "newdef " << name() << ":KmPipPi0K892mMagnitude       " 
  // 	 << _a2K892m   << "\n";
  // output << "newdef " << name() << ":KmPipPi0K892mPhase           " 
  // 	 << _phi2K892m << "\n";
  // output << "newdef " << name() << ":KmPipPi0RhoMagnitude         " 
  // 	 << _a2rho     << "\n";
  // output << "newdef " << name() << ":KmPipPi0RhoPhase             " 
  // 	 << _phi2rho   << "\n";
  // output << "newdef " << name() << ":K0PipPimNonResonantMagnitude " 
  // 	 << _a3NR      << "\n";
  // output << "newdef " << name() << ":K0PipPimNonResonantPhase     " 
  // 	 << _phi3NR    << "\n";
  // output << "newdef " << name() << ":K0PipPimK892Magnitude        " 
  // 	 << _a3K892    << "\n";
  // output << "newdef " << name() << ":K0PipPimK892Phase            " 
  // 	 << _phi3K892  << "\n";
  // output << "newdef " << name() << ":K0PipPimRhoMagnitude         " 
  // 	 << _a3rho     << "\n";
  // output << "newdef " << name() << ":K0PipPimRhoPhase             " 
  // 	 << _phi3rho   << "\n";
  // output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  // output << "newdef " << name() << ":K8920Mass      " << _mK8920/GeV << "\n";
  // output << "newdef " << name() << ":K8920Width     " << _wK8920/GeV << "\n";
  // output << "newdef " << name() << ":K892MinusMass  " << _mK892m/GeV << "\n";
  // output << "newdef " << name() << ":K892MinusWidth " << _wK892m/GeV << "\n";
  // output << "newdef " << name() << ":K1680Mass      " << _mK1680/GeV << "\n";
  // output << "newdef " << name() << ":K1680Width     " << _wK1680/GeV << "\n";
  // output << "newdef " << name() << ":K1430Mass      " << _mK1430/GeV << "\n";
  // output << "newdef " << name() << ":K1430Width     " << _wK1430/GeV << "\n";
  // output << "newdef " << name() << ":Rho0Mass       " << _mrho0 /GeV << "\n";
  // output << "newdef " << name() << ":Rho0Width      " << _wrho0 /GeV << "\n";
  // output << "newdef " << name() << ":RhoPlusMass    " << _mrhop /GeV << "\n";
  // output << "newdef " << name() << ":RhoPlusWidth   " << _wrhop /GeV << "\n";
  // for(unsigned int ix=0;ix<_maxwgt.size();++ix) {
  //   output << "insert " << name() << ":MaximumWeights " 
  // 	   << ix << " " << _maxwgt[ix] << "\n";
  // }
  // for(unsigned int ix=0;ix<_weights.size();++ix) {
  //   output << "insert " << name() << ":Weights " 
  // 	   << ix << " " << _weights[ix] << "\n";
  // }
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}

void DtoKPiPiFOCUS::doinit() {
  DecayIntegrator::doinit();
  // set up the phase space
  // intermediate resonances
  tPDPtr resonances[5]={getParticleData(-313   ),
			getParticleData(-100313),
			getParticleData(-30313 ),
			getParticleData(-315   ),
			getParticleData(-10311 )};
  // D+ -> K-pi+pi+
  tPDPtr in    =  getParticleData(ParticleID::Dplus);
  tPDVector out= {getParticleData(ParticleID::Kminus),
  		  getParticleData(ParticleID::piplus),
  		  getParticleData(ParticleID::piplus)};
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,maxWgt_));
  vector<PhaseSpaceChannel> channels;
  for(tPDPtr res : resonances) {
    channels.push_back((PhaseSpaceChannel(mode),0,res,0,3,1,1,1,2));
    channels.push_back((PhaseSpaceChannel(mode),0,res,0,2,1,1,1,3));
  }
  // add the mode
  for(unsigned int ix=0;ix<channels.size();++ix) {
    channels[ix].weight(weights_[ix]);
    mode->addChannel(channels[ix]);
  }
  addMode(mode);
}

void DtoKPiPiFOCUS::doinitrun() {
  DecayIntegrator::doinitrun();
}
