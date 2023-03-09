// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardModelHadronSpectrum class.
//

#include "StandardModelHadronSpectrum.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/Repository.h>


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

namespace {
  bool weightIsLess (pair<long,double> a, pair<long,double> b) {
    return a.second < b.second;
  }

  /**
   * Return true if the particle pointer corresponds to a diquark 
   * or anti-diquark carrying b flavour; false otherwise.
   */
  inline bool isDiquarkWithB(tcPDPtr par1) {
    if (!par1) return false;
    long id1 = par1->id();
    return DiquarkMatcher::Check(id1)  &&  (abs(id1)/1000)%10 == ParticleID::b;
  }
  
  /**
   * Return true if the particle pointer corresponds to a diquark
   *  or anti-diquark carrying c flavour; false otherwise.
   */
  inline bool isDiquarkWithC(tcPDPtr par1) {
    if (!par1) return false;
    long id1 = par1->id();
    return ( DiquarkMatcher::Check(id1)  &&  
       ( (abs(id1)/1000)%10 == ParticleID::c  
         || (abs(id1)/100)%10 == ParticleID::c ) );
  }

}


StandardModelHadronSpectrum::StandardModelHadronSpectrum(unsigned int opt) 
  : HadronSpectrum(),
    _pwtDquark( 1.0 ),_pwtUquark( 1.0 ),_pwtSquark( 1.0 ),_pwtCquark( 0.0 ),
    _pwtBquark( 0.0 ),
    _sngWt( 1.0 ),_decWt( 1.0 ), 
    _weight1S0(Nmax,1.),_weight3S1(Nmax,1.),_weight1P1(Nmax,1.),_weight3P0(Nmax,1.),
    _weight3P1(Nmax,1.),_weight3P2(Nmax,1.),_weight1D2(Nmax,1.),_weight3D1(Nmax,1.),
    _weight3D2(Nmax,1.),_weight3D3(Nmax,1.),
    _topt(opt),_trial(0), 
    _limBottom(), _limCharm(), _limExotic() 
{
  // The mixing angles
  // the ideal mixing angle
  const double idealAngleMix = atan( sqrt(0.5) ) * 180.0 / Constants::pi;
  // \eta-\eta' mixing angle
  _etamix   = -23.0;
  // phi-omega mixing angle
  _phimix   = +36.0;
  // h_1'-h_1 mixing angle
  _h1mix    = idealAngleMix;
  // f_0(1710)-f_0(1370) mixing angle
  _f0mix    = idealAngleMix;
  // f_1(1420)-f_1(1285)\f$ mixing angle
  _f1mix    = idealAngleMix;
  // f'_2-f_2\f$ mixing angle
  _f2mix    = +26.0;
  // eta_2(1870)-eta_2(1645) mixing angle
  _eta2mix  = idealAngleMix;
  // phi(???)-omega(1650) mixing angle
  _omhmix   = idealAngleMix;
  // phi_3-omega_3 mixing angle
  _ph3mix   = +28.0;
  // eta(1475)-eta(1295) mixing angle
  _eta2Smix = idealAngleMix;
  // phi(1680)-omega(1420) mixing angle
  _phi2Smix = idealAngleMix;
}


StandardModelHadronSpectrum::~StandardModelHadronSpectrum() {}


void StandardModelHadronSpectrum::persistentOutput(PersistentOStream & os) const {
  os << _pwtDquark  << _pwtUquark << _pwtSquark 
     << _pwtCquark << _pwtBquark
     << _etamix << _phimix << _h1mix << _f0mix << _f1mix << _f2mix 
     << _eta2mix << _omhmix << _ph3mix << _eta2Smix << _phi2Smix 
     << _weight1S0 << _weight3S1 << _weight1P1 << _weight3P0 << _weight3P1 
     << _weight3P2 << _weight1D2 << _weight3D1 << _weight3D2 << _weight3D3
     << _sngWt << _decWt << _repwt
     << _limBottom << _limCharm << _limExotic;
}

void StandardModelHadronSpectrum::persistentInput(PersistentIStream & is, int) {
  is >> _pwtDquark  >> _pwtUquark >> _pwtSquark 
     >> _pwtCquark >> _pwtBquark 
     >> _etamix >> _phimix >> _h1mix >> _f0mix >> _f1mix >> _f2mix 
     >> _eta2mix >> _omhmix >> _ph3mix >> _eta2Smix >> _phi2Smix 
     >> _weight1S0 >> _weight3S1 >> _weight1P1 >> _weight3P0 >> _weight3P1 
     >> _weight3P2 >> _weight1D2 >> _weight3D1 >> _weight3D2 >> _weight3D3
     >> _sngWt >> _decWt >> _repwt
     >> _limBottom >> _limCharm >> _limExotic;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<StandardModelHadronSpectrum,HadronSpectrum>
describeHerwigStandardModelHadronSpectrum("Herwig::StandardModelHadronSpectrum", "Herwig.so");

void StandardModelHadronSpectrum::Init() {

  static ClassDocumentation<StandardModelHadronSpectrum> documentation
    ("There is no documentation for the StandardModelHadronSpectrum class");

  static Parameter<StandardModelHadronSpectrum,double>
    interfacePwtDquark("PwtDquark","Weight for choosing a quark D",
		       &StandardModelHadronSpectrum::_pwtDquark, 0, 1.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<StandardModelHadronSpectrum,double>
    interfacePwtUquark("PwtUquark","Weight for choosing a quark U",
		       &StandardModelHadronSpectrum::_pwtUquark, 0, 1.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<StandardModelHadronSpectrum,double>
    interfacePwtSquark("PwtSquark","Weight for choosing a quark S",
		       &StandardModelHadronSpectrum::_pwtSquark, 0, 1.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<StandardModelHadronSpectrum,double>
    interfacePwtCquark("PwtCquark","Weight for choosing a quark C",
		       &StandardModelHadronSpectrum::_pwtCquark, 0, 0.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<StandardModelHadronSpectrum,double>
    interfacePwtBquark("PwtBquark","Weight for choosing a quark B",
		       &StandardModelHadronSpectrum::_pwtBquark, 0, 0.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<StandardModelHadronSpectrum,double>
    interfaceSngWt("SngWt","Weight for singlet baryons",
                  &StandardModelHadronSpectrum::_sngWt, 0, 1.0, 0.0, 10.0,
		   false,false,false);

  static Parameter<StandardModelHadronSpectrum,double>
    interfaceDecWt("DecWt","Weight for decuplet baryons",
                  &StandardModelHadronSpectrum::_decWt, 0, 1.0, 0.0, 10.0,
		   false,false,false);

  //
  // mixing angles
  //
  // the ideal mixing angle
  const double idealAngleMix = atan( sqrt(0.5) ) * 180.0 / Constants::pi;

  static Parameter<StandardModelHadronSpectrum,double> interface11S0Mixing
    ("11S0Mixing",
     "The mixing angle for the I=0 mesons from the 1 1S0 multiplet,"
     " i.e. eta and etaprime.",
     &StandardModelHadronSpectrum::_etamix, -23., -180., 180.,
     false, false, Interface::limited);

  static Parameter<StandardModelHadronSpectrum,double> interface13S1Mixing
    ("13S1Mixing",
     "The mixing angle for the I=0 mesons from the 1 3S1 multiplet,"
     " i.e. phi and omega.",
     &StandardModelHadronSpectrum::_phimix, +36., -180., 180.,
     false, false, Interface::limited);

  static Parameter<StandardModelHadronSpectrum,double> interface11P1Mixing
    ("11P1Mixing",
     "The mixing angle for the I=0 mesons from the 1 1P1 multiplet,"
     " i.e. h_1' and h_1.",
     &StandardModelHadronSpectrum::_h1mix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<StandardModelHadronSpectrum,double> interface13P0Mixing
    ("13P0Mixing",
     "The mixing angle for the I=0 mesons from the 1 3P0 multiplet,"
     " i.e. f_0(1710) and f_0(1370).",
     &StandardModelHadronSpectrum::_f0mix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<StandardModelHadronSpectrum,double> interface13P1Mixing
    ("13P1Mixing",
     "The mixing angle for the I=0 mesons from the 1 3P1 multiplet,"
     " i.e. f_1(1420) and f_1(1285).",
     &StandardModelHadronSpectrum::_f1mix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<StandardModelHadronSpectrum,double> interface13P2Mixing
    ("13P2Mixing",
     "The mixing angle for the I=0 mesons from the 1 3P2 multiplet,"
     " i.e. f'_2 and f_2.",
     &StandardModelHadronSpectrum::_f2mix, 26.0, -180., 180.,
     false, false, Interface::limited);

  static Parameter<StandardModelHadronSpectrum,double> interface11D2Mixing
    ("11D2Mixing",
     "The mixing angle for the I=0 mesons from the 1 1D2 multiplet,"
     " i.e. eta_2(1870) and eta_2(1645).",
     &StandardModelHadronSpectrum::_eta2mix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<StandardModelHadronSpectrum,double> interface13D0Mixing
    ("13D0Mixing",
     "The mixing angle for the I=0 mesons from the 1 3D0 multiplet,"
     " i.e. eta_2(1870) phi(?) and omega(1650).",
     &StandardModelHadronSpectrum::_omhmix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<StandardModelHadronSpectrum,double> interface13D1Mixing
    ("13D1Mixing",
     "The mixing angle for the I=0 mesons from the 1 3D1 multiplet,"
     " i.e. phi_3 and omega_3.",
     &StandardModelHadronSpectrum::_ph3mix, 28.0, -180., 180.,
     false, false, Interface::limited);

  static Parameter<StandardModelHadronSpectrum,double> interface21S0Mixing
    ("21S0Mixing",
     "The mixing angle for the I=0 mesons from the 2 1S0 multiplet,"
     " i.e. eta(1475) and eta(1295).",
     &StandardModelHadronSpectrum::_eta2Smix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<StandardModelHadronSpectrum,double> interface23S1Mixing
    ("23S1Mixing",
     "The mixing angle for the I=0 mesons from the 1 3S1 multiplet,"
     " i.e. phi(1680) and omega(1420).",
     &StandardModelHadronSpectrum::_phi2Smix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);
  //
  //  the meson weights
  //
  static ParVector<StandardModelHadronSpectrum,double> interface1S0Weights
    ("1S0Weights",
     "The weights for the 1S0 multiplets start with n=1.",
     &StandardModelHadronSpectrum::_weight1S0, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<StandardModelHadronSpectrum,double> interface3S1Weights
    ("3S1Weights",
     "The weights for the 3S1 multiplets start with n=1.",
     &StandardModelHadronSpectrum::_weight3S1, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<StandardModelHadronSpectrum,double> interface1P1Weights
    ("1P1Weights",
     "The weights for the 1P1 multiplets start with n=1.",
     &StandardModelHadronSpectrum::_weight1P1, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<StandardModelHadronSpectrum,double> interface3P0Weights
    ("3P0Weights",
     "The weights for the 3P0 multiplets start with n=1.",
     &StandardModelHadronSpectrum::_weight3P0, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<StandardModelHadronSpectrum,double> interface3P1Weights
    ("3P1Weights",
     "The weights for the 3P1 multiplets start with n=1.",
     &StandardModelHadronSpectrum::_weight3P1, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<StandardModelHadronSpectrum,double> interface3P2Weights
    ("3P2Weights",
     "The weights for the 3P2 multiplets start with n=1.",
     &StandardModelHadronSpectrum::_weight3P2, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<StandardModelHadronSpectrum,double> interface1D2Weights
    ("1D2Weights",
     "The weights for the 1D2 multiplets start with n=1.",
     &StandardModelHadronSpectrum::_weight1D2, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<StandardModelHadronSpectrum,double> interface3D1Weights
    ("3D1Weights",
     "The weights for the 3D1 multiplets start with n=1.",
     &StandardModelHadronSpectrum::_weight3D1, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<StandardModelHadronSpectrum,double> interface3D2Weights
    ("3D2Weights",
     "The weights for the 3D2 multiplets start with n=1.",
     &StandardModelHadronSpectrum::_weight3D2, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<StandardModelHadronSpectrum,double> interface3D3Weights
    ("3D3Weights",
     "The weights for the 3D3 multiplets start with n=1.",
     &StandardModelHadronSpectrum::_weight3D3, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Switch<StandardModelHadronSpectrum,unsigned int> interfaceTrial
    ("Trial",
     "A Debugging option to only produce certain types of hadrons",
     &StandardModelHadronSpectrum::_trial, 0, false, false);
  static SwitchOption interfaceTrialAll
    (interfaceTrial,
     "All",
     "Produce all the hadrons",
     0);
  static SwitchOption interfaceTrialPions
    (interfaceTrial,
     "Pions",
     "Only produce pions",
     1);
  static SwitchOption interfaceTrialSpin2
    (interfaceTrial,
     "Spin2",
     "Only mesons with spin less than or equal to two are produced",
     2);
  static SwitchOption interfaceTrialSpin3
    (interfaceTrial,
     "Spin3",
     "Only hadrons with spin less than or equal to three are produced",
     3);

  static Parameter<StandardModelHadronSpectrum,double>
    interfaceSingleHadronLimitBottom ("SingleHadronLimitBottom",
				      "Threshold for one-hadron decay of b-cluster",
				      &StandardModelHadronSpectrum::_limBottom,
				      0, 0.0, 0.0, 100.0,false,false,false);

  static Parameter<StandardModelHadronSpectrum,double>
    interfaceSingleHadronLimitCharm ("SingleHadronLimitCharm",
				     "threshold for one-hadron decay of c-cluster",
				     &StandardModelHadronSpectrum::_limCharm,
				     0, 0.0, 0.0, 100.0,false,false,false);

  static Parameter<StandardModelHadronSpectrum,double>
    interfaceSingleHadronLimitExotic ("SingleHadronLimitExotic",
				      "threshold for one-hadron decay of exotic cluster",
				      &StandardModelHadronSpectrum::_limExotic,
				      0, 0.0, 0.0, 100.0,false,false,false);

  static Switch<StandardModelHadronSpectrum,unsigned int> interfaceBelowThreshold
    ("BelowThreshold",
     "Option fo the selection of the hadrons if the cluster is below the pair threshold",
     &StandardModelHadronSpectrum::belowThreshold_, 0, false, false);
  static SwitchOption interfaceBelowThresholdLightest
    (interfaceBelowThreshold,
     "Lightest",
     "Force cluster to decay to the lightest hadron with the appropriate flavours",
     0);
  static SwitchOption interfaceBelowThresholdAll
    (interfaceBelowThreshold,
     "All",
     "Select from all the hadrons below the two hadron threshold according to their spin weights",
     1);

}


PDPtr StandardModelHadronSpectrum::makeDiquark(tcPDPtr par1, tcPDPtr par2) const {
    long id1 = par1->id();
    long id2 = par2->id();
    long pspin = id1==id2 ? 3 : 1;
    long idnew = makeDiquarkID(id1,id2, pspin);
    return getParticleData(idnew);
}

Energy StandardModelHadronSpectrum::hadronPairThreshold(tcPDPtr par1, tcPDPtr par2) const {
  // Determine the sum of the nominal masses of the two lightest hadrons
  // with the right flavour numbers as the cluster under consideration.
  // Notice that we don't need real masses (drawn by a Breit-Wigner 
  // distribution) because the lightest pair of hadrons does not involve
  // any broad resonance.
  Energy threshold = massLightestHadronPair(par1,par2);
  // Special: it allows one-hadron decays also above threshold.
  if (isExotic(par1,par2)) 
    threshold *= (1.0 + UseRandom::rnd()*_limExotic);
  else if (hasBottom(par1,par2)) 
    threshold *= (1.0 + UseRandom::rnd()*_limBottom);
  else if (hasCharm(par1,par2)) 
    threshold *= (1.0 + UseRandom::rnd()*_limCharm);
  return threshold;
}
  
double StandardModelHadronSpectrum::mixingStateWeight(long id) const {
  switch(id) {
  case ParticleID::eta:      return 0.5*probabilityMixing(_etamix  ,1);
  case ParticleID::etaprime: return 0.5*probabilityMixing(_etamix  ,2);
  case ParticleID::phi:      return 0.5*probabilityMixing(_phimix  ,1);
  case ParticleID::omega:    return 0.5*probabilityMixing(_phimix  ,2);
  case ParticleID::hprime_1: return 0.5*probabilityMixing(_h1mix   ,1);
  case ParticleID::h_1:      return 0.5*probabilityMixing(_h1mix   ,2);
  case 10331:                return 0.5*probabilityMixing(_f0mix   ,1);
  case 10221:                return 0.5*probabilityMixing(_f0mix   ,2);
  case ParticleID::fprime_1: return 0.5*probabilityMixing(_f1mix   ,1);
  case ParticleID::f_1:      return 0.5*probabilityMixing(_f1mix   ,2);
  case ParticleID::fprime_2: return 0.5*probabilityMixing(_f2mix   ,1);
  case ParticleID::f_2:      return 0.5*probabilityMixing(_f2mix   ,2);
  case 10335:                return 0.5*probabilityMixing(_eta2mix ,1);
  case 10225:		     return 0.5*probabilityMixing(_eta2mix ,2);
  case 30333:		     return 0.5*probabilityMixing(_omhmix  ,1);
  case 30223:		     return 0.5*probabilityMixing(_omhmix  ,2);
  case 337:                  return 0.5*probabilityMixing(_ph3mix  ,1);
  case 227:		     return 0.5*probabilityMixing(_ph3mix  ,2);
  case 100331:               return 0.5*probabilityMixing(_eta2mix ,1);
  case 100221:		     return 0.5*probabilityMixing(_eta2mix ,2);
  case 100333:               return 0.5*probabilityMixing(_phi2Smix,1);
  case 100223:		     return 0.5*probabilityMixing(_phi2Smix,2);
  default:                   return 1./3.;
  }
}

void StandardModelHadronSpectrum::doinit() {
  // set the weights for the various excited mesons
  // set all to one to start with
  for (int l = 0; l < Lmax; ++l ) {
    for (int j = 0; j < Jmax; ++j) {
      for (int n = 0; n < Nmax; ++n) {
	_repwt[l][j][n] = 1.0;
      }
    }
  }
  // set the others from the relevant vectors
  for( int ix=0;ix<max(int(_weight1S0.size()),int(Nmax));++ix)
    _repwt[0][0][ix]=_weight1S0[ix];
  for( int ix=0;ix<max(int(_weight3S1.size()),int(Nmax));++ix)
    _repwt[0][1][ix]=_weight3S1[ix];
  for( int ix=0;ix<max(int(_weight1P1.size()),int(Nmax));++ix)
    _repwt[1][1][ix]=_weight1P1[ix];
  for( int ix=0;ix<max(int(_weight3P0.size()),int(Nmax));++ix)
    _repwt[1][0][ix]=_weight3P0[ix];
  for( int ix=0;ix<max(int(_weight3P1.size()),int(Nmax));++ix)
    _repwt[1][1][ix]=_weight3P1[ix];
  for( int ix=0;ix<max(int(_weight3P2.size()),int(Nmax));++ix)
    _repwt[1][2][ix]=_weight3P2[ix];
  for( int ix=0;ix<max(int(_weight1D2.size()),int(Nmax));++ix)
    _repwt[2][2][ix]=_weight1D2[ix];
  for( int ix=0;ix<max(int(_weight3D1.size()),int(Nmax));++ix)
    _repwt[2][1][ix]=_weight3D1[ix];
  for( int ix=0;ix<max(int(_weight3D2.size()),int(Nmax));++ix)
    _repwt[2][2][ix]=_weight3D2[ix];
  for( int ix=0;ix<max(int(_weight3D3.size()),int(Nmax));++ix)
    _repwt[2][3][ix]=_weight3D3[ix];

  // find the maximum
  map<long,double>::iterator pit =
    max_element(_pwt.begin(),_pwt.end(),weightIsLess); 
  const double pmax = pit->second;
  for(pit=_pwt.begin(); pit!=_pwt.end(); ++pit) {
    pit->second/=pmax;
  }
  HadronSpectrum::doinit();
}

void StandardModelHadronSpectrum::constructHadronTable() {
  // initialise the table
  _table.clear();
  for(unsigned int ix=0; ix<_partons.size(); ++ix) {
    for(unsigned int iy=0; iy<_partons.size(); ++iy) {
      if (!(DiquarkMatcher::Check(_partons[ix]->id()) 
	    && DiquarkMatcher::Check(_partons[iy]->id())))
      _table[make_pair(_partons[ix]->id(),_partons[iy]->id())] = KupcoData();
    }
  }
  // get the particles from the event generator
  ParticleMap particles = generator()->particles();
  // loop over the particles
  //double maxdd(0.),maxss(0.),maxrest(0.);
  for(ParticleMap::iterator it=particles.begin(); 
      it!=particles.end(); ++it) {
    long pid = it->first;
    tPDPtr particle = it->second;
    int pspin = particle->iSpin();
    // Don't include hadrons which are explicitly forbidden
    if(find(_forbidden.begin(),_forbidden.end(),particle)!=_forbidden.end()) 
      continue;
    // Don't include non-hadrons or antiparticles
    if(pid < 100) continue;
    // remove diffractive particles
    if(pspin == 0) continue;
    // K_0S and K_0L not made make K0 and Kbar0
    if(pid==ParticleID::K_S0||pid==ParticleID::K_L0) continue;
    // Debugging options
    // Only include those with 2J+1 less than...5
    if(_trial==2 && pspin >= 5) continue;
    // Only include those with 2J+1 less than...7
    if(_trial==3 && pspin >= 7) continue;
    // Only include pions
    if(_trial==1 && pid!=111 && pid!=211) continue;
    // shouldn't be coloured
    if(particle->coloured()) continue;
    // Get the flavours
    const int x4 = (pid/1000)%10; 
    const int x3 = (pid/100 )%10;
    const int x2 = (pid/10  )%10;
    const int x7 = (pid/1000000)%10;
    const bool wantSusy = x7 == 1 || x7 == 2;
    // Skip non-hadrons (susy particles, etc...)
    if(x3 == 0 || x2 == 0) continue;
    // Skip particles which are neither SM nor SUSY 
    if(x7 >= 3 && x7 != 9) continue;
    int flav1,flav2;
    // meson
    if(x4 == 0) {
      flav1 = x2;
      flav2 = x3;
    }
    // baryon
    else {
      flav2 = x4;
      // insert the spin 1 diquark, sort out the rest later
      flav1 = makeDiquarkID(x2,x3,3);
    }
    if (wantSusy) flav2 += 1000000 * x7;
    insertToHadronTable(particle,flav1,flav2);
  }
  // normalise the weights


  if(_topt == 0) {
    HadronTable::const_iterator tit;
    KupcoData::iterator it;
    for(tit=_table.begin();tit!=_table.end();++tit) {
      double weight=0;
      for(it = tit->second.begin(); it!=tit->second.end(); ++it)
	weight=max(weight,it->overallWeight);
      weight = 1./weight;
    }

    //   double weight;
    //   if(tit->first.first==tit->first.second) {
    // 	if(tit->first.first==1||tit->first.first==2) weight=1./maxdd;
    // 	else if (tit->first.first==3)                weight=1./maxss;
    // 	else                                         weight=1./maxrest;
    //   }
    //   else                                           weight=1./maxrest;
    //   for(it = tit->second.begin(); it!=tit->second.end(); ++it) {
    // 	it->rescale(weight);
    //   }
    // }
  }
}

double StandardModelHadronSpectrum::strangeWeight(const Energy, tcPDPtr, tcPDPtr) const {
  assert(false);
}

void StandardModelHadronSpectrum::insertMeson(HadronInfo a, int flav1, int flav2) {
  // identical light flavours
  if(flav1 == flav2 && flav1<=3) {
    // ddbar> uubar> admixture states
    if(flav1==1) {
      a.overallWeight *= 0.5;
      _table[make_pair(1,1)].insert(a);
      _table[make_pair(2,2)].insert(a);
    }
    // load up ssbar> uubar> ddbar> admixture states
    else {
      // uubar ddbar pieces
      a.wt = mixingStateWeight(a.id);
      a.overallWeight *= a.wt;
      _table[make_pair(1,1)].insert(a);
      _table[make_pair(2,2)].insert(a);
      a.overallWeight /=a.wt;
      // ssbar piece
      a.wt = 1.- 2.*a.wt;
      if(a.wt > 0) {
        a.overallWeight *= a.wt;
        _table[make_pair(3,3)].insert(a);
      }
    }
  }
  else {
    _table[make_pair(flav1,flav2)].insert(a);
    if(flav1 != flav2) _table[make_pair(flav2,flav1)].insert(a);
  }
}


long StandardModelHadronSpectrum::makeDiquarkID(long id1, long id2, long pspin) const {

  assert( id1 * id2 > 0  
          && QuarkMatcher::Check(id1)  
	  && QuarkMatcher::Check(id2)) ;
  long ida = abs(id1);
  long idb = abs(id2);
  if (ida < idb) swap(ida,idb);
  if (pspin != 1 && pspin != 3) assert(false);
  long idnew = ida*1000 + idb*100 + pspin;
  // Diquarks made of quarks of the same type: uu, dd, ss, cc, bb,
  // have spin 1, and therefore the less significant digit (which
  // corresponds to 2*J+1) is 3 rather than 1 as all other Diquarks.
  if (id1 == id2 && pspin == 1) {
    //cerr<<"WARNING: spin-0 diquiark of the same type cannot exist."
    //    <<" Switching to spin-1 diquark.\n";
    idnew = ida*1000 + idb*100 + 3;
  }

  return id1 > 0 ? idnew : -idnew;
}

bool StandardModelHadronSpectrum::hasBottom(tcPDPtr par1, tcPDPtr par2, tcPDPtr par3) const {
  long id1 = par1 ? par1->id() : 0;
  if ( !par2  &&  !par3 ) {
    return 
      abs(id1) == ThePEG::ParticleID::b    ||
      isDiquarkWithB(par1)                 ||
      ( MesonMatcher::Check(id1)  
	&& (abs(id1)/100)%10  == ThePEG::ParticleID::b ) ||
      ( BaryonMatcher::Check(id1) 
	&& (abs(id1)/1000)%10 == ThePEG::ParticleID::b );
  } 
  else {
    long id2 = par2 ? par2->id() : 0;
    long id3 = par3 ? par3->id() : 0;
    return 
      abs(id1) == ThePEG::ParticleID::b  ||  isDiquarkWithB(par1)  || 
      abs(id2) == ThePEG::ParticleID::b  ||  isDiquarkWithB(par2)  || 
      abs(id3) == ThePEG::ParticleID::b  ||  isDiquarkWithB(par3); 
  }
}


bool StandardModelHadronSpectrum::hasCharm(tcPDPtr par1, tcPDPtr par2, tcPDPtr par3) const {
  long id1 = par1 ? par1->id(): 0;
  if (!par2  &&  !par3) {
    return
      abs(id1) == ThePEG::ParticleID::c     ||
      isDiquarkWithC(par1)                  ||
      ( MesonMatcher::Check(id1) && 
        ((abs(id1)/100)%10 == ThePEG::ParticleID::c ||
	 (abs(id1)/10)%10 == ThePEG::ParticleID::c) ) ||
      ( BaryonMatcher::Check(id1) && 
        ((abs(id1)/1000)%10 == ThePEG::ParticleID::c  ||
	 (abs(id1)/100)%10  == ThePEG::ParticleID::c  ||
	 (abs(id1)/10)%10   == ThePEG::ParticleID::c) );
  } 
  else {
 long id2 = par2 ? par1->id(): 0;
 long id3 = par3 ? par1->id(): 0;
    return 
      abs(id1) == ThePEG::ParticleID::c  ||  isDiquarkWithC(par1)  || 
      abs(id2) == ThePEG::ParticleID::c  ||  isDiquarkWithC(par2)  || 
      abs(id3) == ThePEG::ParticleID::c  ||  isDiquarkWithC(par3); 
  }
}  

bool StandardModelHadronSpectrum::isExotic(tcPDPtr par1, tcPDPtr par2, tcPDPtr par3) const {
  /// \todo make this more general
  long id1 = par1 ? par1->id(): 0;
  long id2 = par2 ? par2->id(): 0;
  long id3 = par3 ? par3->id(): 0;
return 
  ( (id1/1000000)% 10 != 0 && (id1/1000000)% 10 != 9 ) ||
  ( (id2/1000000)% 10 != 0 && (id2/1000000)% 10 != 9 ) ||
  ( (id3/1000000)% 10 != 0 && (id3/1000000)% 10 != 9 ) ||
  abs(id1)==6||abs(id2)==6;
}


bool StandardModelHadronSpectrum::canBeBaryon(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3) const {
  assert(par1 && par2);
  long id1 = par1->id(), id2 = par2->id();
  if (!par3) {
    if( id1*id2 < 0) return false;
    if(DiquarkMatcher::Check(id1))
return abs(int(par2->iColour())) == 3 && !DiquarkMatcher::Check(id2); 
    if(DiquarkMatcher::Check(id2))
return abs(int(par1->iColour())) == 3;
    return false;
  } 
  else {
    // In this case, to be a baryon, all three components must be (anti-)quarks
    // and with the same sign.
    return (par1->iColour() == 3 && par2->iColour() == 3 && par3->iColour() == 3) ||
(par1->iColour() == -3 && par2->iColour() == -3 && par3->iColour() == -3);
  }
}
