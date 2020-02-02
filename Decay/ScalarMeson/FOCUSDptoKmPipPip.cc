// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FOCUSDptoKmPipPip class.
//

#include "FOCUSDptoKmPipPip.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FOCUSDptoKmPipPip::FOCUSDptoKmPipPip() : WeakDalitzDecay(5./GeV,1.5/GeV),
					 g1_(0.31072*GeV), g2_(-0.02323*GeV),
					 beta_(3.389*GeV), theta_(286), gamma_({304,126,211}),
					 c1_({1.655,0.780,-0.954}), c2_(17.182), c3_(0.734)
{}

IBPtr FOCUSDptoKmPipPip::clone() const {
  return new_ptr(*this);
}

IBPtr FOCUSDptoKmPipPip::fullclone() const {
  return new_ptr(*this);
}

void FOCUSDptoKmPipPip::doinit() {
  WeakDalitzDecay::doinit();
  // static const double degtorad = Constants::pi/180.;
  // // create the resonances
  // addResonance(DalitzResonance(getParticleData(-313    ), 896  *MeV, 50.3*MeV,0,1,2,-1.   ,   0.          ));
  // addResonance(DalitzResonance(getParticleData(-313    ), 896  *MeV, 50.3*MeV,0,2,1,-1.   ,   0.          ));
  // addResonance(DalitzResonance(getParticleData(-10311  ),1463  *MeV,163.8*MeV,0,1,2,3.   ,  49.7*degtorad));
  // addResonance(DalitzResonance(getParticleData(-10311  ),1463  *MeV,163.8*MeV,0,2,1,3.   ,  49.7*degtorad));
  // addResonance(DalitzResonance(getParticleData(-315    ),1432.4*MeV,  109*MeV,0,1,2,0.962, -29.9*degtorad));
  // addResonance(DalitzResonance(getParticleData(-315    ),1432.4*MeV,  109*MeV,0,2,1,0.962, -29.9*degtorad));
  // addResonance(DalitzResonance(getParticleData(-30313  ),1717  *MeV,  322*MeV,0,1,2,-6.5  ,  29.0*degtorad));
  // addResonance(DalitzResonance(getParticleData(-30313  ),1717  *MeV,  322*MeV,0,2,1,-6.5  ,  29.0*degtorad));
  // addResonance(DalitzResonance(getParticleData(-9000311), 809  *MeV,  470*MeV,0,1,2,5.01 ,-163.7*degtorad));
  // addResonance(DalitzResonance(getParticleData(-9000311), 809  *MeV,  470*MeV,0,2,1,5.01 ,-163.7*degtorad));
  // // D+ -> K- pi+ pi+
  // createMode(getParticleData(ParticleID::Dplus),
  // 	     {getParticleData(ParticleID::Kminus),
  // 		 getParticleData(ParticleID::piplus),
  // 		 getParticleData(ParticleID::piplus)});
}

void FOCUSDptoKmPipPip::doinitrun() {
  WeakDalitzDecay::doinitrun();
}

void FOCUSDptoKmPipPip::persistentOutput(PersistentOStream & os) const {
  os << ounit(beta_,GeV) << theta_ << gamma_
     << ounit(g1_,GeV) << ounit(g2_,GeV) << c1_ << c2_ << c3_
     << KHalf_ << KThreeHalf_;
}

void FOCUSDptoKmPipPip::persistentInput(PersistentIStream & is, int) {
  is >> iunit(beta_,GeV) >> theta_ >> gamma_
     >> iunit(g1_,GeV) >> iunit(g2_,GeV) >> c1_ >> c2_ >> c3_
     >> KHalf_ >> KThreeHalf_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FOCUSDptoKmPipPip,WeakDalitzDecay>
describeHerwigFOCUSDptoKmPipPip("Herwig::FOCUSDptoKmPipPip", "HwSMDecay.so");

void FOCUSDptoKmPipPip::Init() {

  static ClassDocumentation<FOCUSDptoKmPipPip> documentation
    ("The FOCUSDptoKmPipPip class implements the fits of FOCUS, "
     "Phys.Lett. B653 (2007) 1-11 for the decay D+ -> K-pi+pi-.",
     "The decay $D^+\\to K^-\\pi^+\\pi^+$ was simulated using the"
     " model of FOCUS \\cite{Pennington:2007se}.",
     "\\bibitem{Pennington:2007se}\n"
     "J.~M.~Link {\\it et al.} [FOCUS Collaboration],\n"
     "%``Dalitz plot analysis of the $D^{+} \\to K^{-} \\pi^{+} \\pi^{+}$ decay in the FOCUS experiment,''\n"
     "Phys.\\ Lett.\\ B {\\bf 653} (2007) 1\n"
     "doi:10.1016/j.physletb.2007.06.070\n"
     "[arXiv:0705.2248 [hep-ex]].\n");

  static Reference<FOCUSDptoKmPipPip,KMatrix> interfaceKHalf
    ("KHalf",
     "The K-matrix for the I=1/2 s-wave component",
     &FOCUSDptoKmPipPip::KHalf_, false, false, true, false, false);
  
  static Reference<FOCUSDptoKmPipPip,KMatrix> interfaceKThreeHalf
    ("KThreeHalf",
     "The K-matrix for the I=1/2 s-wave component",
     &FOCUSDptoKmPipPip::KThreeHalf_, false, false, true, false, false);

  static Parameter<FOCUSDptoKmPipPip,Energy> interfaceg1
    ("g1",
     "The coupling of the first channel to K*0",
     &FOCUSDptoKmPipPip::g1_, 0.31072*GeV, -10.*GeV, 10.*GeV,
     false, false, Interface::limited);

  static Parameter<FOCUSDptoKmPipPip,Energy> interfaceg2
    ("g2",
     "The coupling of the first channel to K*0",
     &FOCUSDptoKmPipPip::g2_, -0.02323*GeV, -10.*GeV, 10.*GeV,
     false, false, Interface::limited);

  static Parameter<FOCUSDptoKmPipPip,Energy> interfacebeta
    ("beta",
     "The beta coefficient for the P-vector",
     &FOCUSDptoKmPipPip::beta_, 3.389*GeV, -10.*GeV, 10.*GeV,
     false, false, Interface::limited);

  static Parameter<FOCUSDptoKmPipPip,double> interfacetheta
    ("theta",
     "The theta phase (in degrees) for the P-vector",
     &FOCUSDptoKmPipPip::theta_, 286, 0.0, 360.,
     false, false, Interface::limited);

  static ParVector<FOCUSDptoKmPipPip,double> interfacegamma
    ("gamma",
     "The gamma phases in degress",
     &FOCUSDptoKmPipPip::gamma_, 3, 0., 0.0, 360.0,
     false, false, Interface::limited);

  static ParVector<FOCUSDptoKmPipPip,double> interfacec1
    ("c1",
     "The c1 coefficients for the P-vector",
     &FOCUSDptoKmPipPip::c1_, -1, 0., -100., 100.,
     false, false, Interface::limited);

  static ParVector<FOCUSDptoKmPipPip,double> interfacec2
    ("c2",
     "The c2 coefficients for the P-vector",
     &FOCUSDptoKmPipPip::c2_, -1, 0., -100., 100.,
     false, false, Interface::limited);

  static ParVector<FOCUSDptoKmPipPip,double> interfacec3
    ("c3",
     "The c3 coefficients for the P-vector",
     &FOCUSDptoKmPipPip::c3_, -1, 0., -100., 100.,
     false, false, Interface::limited);

}

int FOCUSDptoKmPipPip::modeNumber(bool & cc,tcPDPtr parent,
				const tPDVector & children) const {
  int id0(parent->id());
  // incoming particle must be D+
  if(abs(id0)!=ParticleID::Dplus) return -1;
  cc = id0<0;
  // must be three decay products
  if(children.size()!=3) return -1;
  unsigned int npip(0),npim(0),nk(0);
  for(tPDPtr child : children) {
    long id= child->id();
    if(id          ==ParticleID::piplus)  ++npip;
    else if(id     ==ParticleID::piminus) ++npim;
    else if(abs(id)==ParticleID::Kplus)   ++nk;
  }
  if((id0==ParticleID::Dplus &&(nk==1&&npip==2))||
     (id0==ParticleID::Dminus&&(nk==1&&npim==2))) return 0;
  else return -1;
}


Complex FOCUSDptoKmPipPip::amplitude(int ichan) const {
  Complex output(0.);
  // unsigned int imin=0, imax(resonances().size());
  // if(ichan>=0) {
  //   imin=ichan;
  //   imax=imin+1;
  // }
  // for(unsigned int ix=imin;ix<imax;++ix) {
  //   if(ix==2 || ix==3 || ix==8 || ix==9) {
  //     output += resAmp(ix,true);
  //   }
  //   else
  //     output += resAmp(ix);
  // }
  // if(ichan<0) {
  //   output += 7.4*Complex(cos(-0.3211405823669566),sin(-0.3211405823669566));
  // }
  return output;
}
