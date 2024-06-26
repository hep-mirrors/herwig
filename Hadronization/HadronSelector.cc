// -*- C++ -*-
//
// HadronSelector.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HadronSelector class.
//

#include "HadronSelector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/CurrentGenerator.h>
#include <ThePEG/Repository/Repository.h>
#include "CheckId.h"
#include <ThePEG/Utilities/DescribeClass.h>
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;

DescribeAbstractClass<HadronSelector,Interfaced>
describeHadronSelector("Herwig::HadronSelector","Herwig.so");

namespace {
  // debug helper
  void dumpTable(const HadronTable & tbl) {
    typedef HadronTable::const_iterator TableIter;
    for (TableIter it = tbl.begin(); it != tbl.end(); ++it) {
      cerr << it->first.first << ' '
       	   << it->first.second << '\n';
      for (KupcoData::const_iterator jt = it->second.begin();
      	   jt != it->second.end(); ++jt) {
      	cerr << '\t' << *jt << '\n';
      }
    }
  }

  bool weightIsLess (pair<long,double> a, pair<long,double> b) {
    return a.second < b.second;
  }
}

HadronSelector::HadronSelector(unsigned int opt)
  : _weight1S0(Nmax,1.),_weight3S1(Nmax,1.),_weight1P1(Nmax,1.),_weight3P0(Nmax,1.),
    _weight3P1(Nmax,1.),_weight3P2(Nmax,1.),_weight1D2(Nmax,1.),_weight3D1(Nmax,1.),
    _weight3D2(Nmax,1.),_weight3D3(Nmax,1.),
    _repwt(Lmax,vector<vector<double> >(Jmax,vector<double>(Nmax))),
    _topt(opt),_trial(0),
    _limBottom(), _limCharm(), _limExotic(), belowThreshold_(0)
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

void HadronSelector::persistentOutput(PersistentOStream & os) const {
  os << _partons 
     << _etamix << _phimix << _h1mix << _f0mix << _f1mix << _f2mix
     << _eta2mix << _omhmix << _ph3mix << _eta2Smix << _phi2Smix
     << _weight1S0 << _weight3S1 << _weight1P1 << _weight3P0 << _weight3P1
     << _weight3P2 << _weight1D2 << _weight3D1 << _weight3D2 << _weight3D3
     << _forbidden << _repwt << _pwt
     << _limBottom << _limCharm << _limExotic << belowThreshold_
     << _table << lightestHadrons_;
}

void HadronSelector::persistentInput(PersistentIStream & is, int) {
  is >> _partons 
     >> _etamix >> _phimix >> _h1mix >> _f0mix >> _f1mix >> _f2mix
     >> _eta2mix >> _omhmix >> _ph3mix >> _eta2Smix >> _phi2Smix
     >> _weight1S0 >> _weight3S1 >> _weight1P1 >> _weight3P0 >> _weight3P1
     >> _weight3P2 >> _weight1D2 >> _weight3D1 >> _weight3D2 >> _weight3D3
     >> _forbidden >> _repwt >> _pwt
     >> _limBottom >> _limCharm >> _limExotic >> belowThreshold_
     >> _table >> lightestHadrons_;
}

void HadronSelector::Init() {

  static ClassDocumentation<HadronSelector> documentation
    ("There is no documentation for the HadronSelector class");

  static RefVector<HadronSelector,ParticleData> interfacePartons
    ("Partons",
     "The partons which are to be considered as the consistuents of the hadrons.",
     &HadronSelector::_partons, -1, false, false, true, false, false);

  static RefVector<HadronSelector,ParticleData> interfaceForbidden
    ("Forbidden",
     "The PDG codes of the particles which cannot be produced in the hadronization.",
     &HadronSelector::_forbidden, -1, false, false, true, false, false);

  //
  // mixing angles
  //
  // the ideal mixing angle
  const double idealAngleMix = atan( sqrt(0.5) ) * 180.0 / Constants::pi;

  static Parameter<HadronSelector,double> interface11S0Mixing
    ("11S0Mixing",
     "The mixing angle for the I=0 mesons from the 1 1S0 multiplet,"
     " i.e. eta and etaprime.",
     &HadronSelector::_etamix, -23., -180., 180.,
     false, false, Interface::limited);

  static Parameter<HadronSelector,double> interface13S1Mixing
    ("13S1Mixing",
     "The mixing angle for the I=0 mesons from the 1 3S1 multiplet,"
     " i.e. phi and omega.",
     &HadronSelector::_phimix, +36., -180., 180.,
     false, false, Interface::limited);

  static Parameter<HadronSelector,double> interface11P1Mixing
    ("11P1Mixing",
     "The mixing angle for the I=0 mesons from the 1 1P1 multiplet,"
     " i.e. h_1' and h_1.",
     &HadronSelector::_h1mix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<HadronSelector,double> interface13P0Mixing
    ("13P0Mixing",
     "The mixing angle for the I=0 mesons from the 1 3P0 multiplet,"
     " i.e. f_0(1710) and f_0(1370).",
     &HadronSelector::_f0mix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<HadronSelector,double> interface13P1Mixing
    ("13P1Mixing",
     "The mixing angle for the I=0 mesons from the 1 3P1 multiplet,"
     " i.e. f_1(1420) and f_1(1285).",
     &HadronSelector::_f1mix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<HadronSelector,double> interface13P2Mixing
    ("13P2Mixing",
     "The mixing angle for the I=0 mesons from the 1 3P2 multiplet,"
     " i.e. f'_2 and f_2.",
     &HadronSelector::_f2mix, 26.0, -180., 180.,
     false, false, Interface::limited);

  static Parameter<HadronSelector,double> interface11D2Mixing
    ("11D2Mixing",
     "The mixing angle for the I=0 mesons from the 1 1D2 multiplet,"
     " i.e. eta_2(1870) and eta_2(1645).",
     &HadronSelector::_eta2mix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<HadronSelector,double> interface13D0Mixing
    ("13D0Mixing",
     "The mixing angle for the I=0 mesons from the 1 3D0 multiplet,"
     " i.e. eta_2(1870) phi(?) and omega(1650).",
     &HadronSelector::_omhmix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<HadronSelector,double> interface13D1Mixing
    ("13D1Mixing",
     "The mixing angle for the I=0 mesons from the 1 3D1 multiplet,"
     " i.e. phi_3 and omega_3.",
     &HadronSelector::_ph3mix, 28.0, -180., 180.,
     false, false, Interface::limited);

  static Parameter<HadronSelector,double> interface21S0Mixing
    ("21S0Mixing",
     "The mixing angle for the I=0 mesons from the 2 1S0 multiplet,"
     " i.e. eta(1475) and eta(1295).",
     &HadronSelector::_eta2Smix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);

  static Parameter<HadronSelector,double> interface23S1Mixing
    ("23S1Mixing",
     "The mixing angle for the I=0 mesons from the 1 3S1 multiplet,"
     " i.e. phi(1680) and omega(1420).",
     &HadronSelector::_phi2Smix, idealAngleMix, -180., 180.,
     false, false, Interface::limited);
  //
  //  the meson weights
  //
  static ParVector<HadronSelector,double> interface1S0Weights
    ("1S0Weights",
     "The weights for the 1S0 multiplets start with n=1.",
     &HadronSelector::_weight1S0, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<HadronSelector,double> interface3S1Weights
    ("3S1Weights",
     "The weights for the 3S1 multiplets start with n=1.",
     &HadronSelector::_weight3S1, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<HadronSelector,double> interface1P1Weights
    ("1P1Weights",
     "The weights for the 1P1 multiplets start with n=1.",
     &HadronSelector::_weight1P1, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<HadronSelector,double> interface3P0Weights
    ("3P0Weights",
     "The weights for the 3P0 multiplets start with n=1.",
     &HadronSelector::_weight3P0, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<HadronSelector,double> interface3P1Weights
    ("3P1Weights",
     "The weights for the 3P1 multiplets start with n=1.",
     &HadronSelector::_weight3P1, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<HadronSelector,double> interface3P2Weights
    ("3P2Weights",
     "The weights for the 3P2 multiplets start with n=1.",
     &HadronSelector::_weight3P2, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<HadronSelector,double> interface1D2Weights
    ("1D2Weights",
     "The weights for the 1D2 multiplets start with n=1.",
     &HadronSelector::_weight1D2, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<HadronSelector,double> interface3D1Weights
    ("3D1Weights",
     "The weights for the 3D1 multiplets start with n=1.",
     &HadronSelector::_weight3D1, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<HadronSelector,double> interface3D2Weights
    ("3D2Weights",
     "The weights for the 3D2 multiplets start with n=1.",
     &HadronSelector::_weight3D2, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<HadronSelector,double> interface3D3Weights
    ("3D3Weights",
     "The weights for the 3D3 multiplets start with n=1.",
     &HadronSelector::_weight3D3, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Switch<HadronSelector,unsigned int> interfaceTrial
    ("Trial",
     "A Debugging option to only produce certain types of hadrons",
     &HadronSelector::_trial, 0, false, false);
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

  static Parameter<HadronSelector,double>
    interfaceSingleHadronLimitBottom ("SingleHadronLimitBottom",
				      "Threshold for one-hadron decay of b-cluster",
				      &HadronSelector::_limBottom,
				      0, 0.0, 0.0, 100.0,false,false,false);

  static Parameter<HadronSelector,double>
    interfaceSingleHadronLimitCharm ("SingleHadronLimitCharm",
				     "threshold for one-hadron decay of c-cluster",
				     &HadronSelector::_limCharm,
				     0, 0.0, 0.0, 100.0,false,false,false);

  static Parameter<HadronSelector,double>
    interfaceSingleHadronLimitExotic ("SingleHadronLimitExotic",
				      "threshold for one-hadron decay of exotic cluster",
				      &HadronSelector::_limExotic,
				      0, 0.0, 0.0, 100.0,false,false,false);

  static Switch<HadronSelector,unsigned int> interfaceBelowThreshold
    ("BelowThreshold",
     "Option fo the selection of the hadrons if the cluster is below the pair threshold",
     &HadronSelector::belowThreshold_, 0, false, false);
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

double HadronSelector::mixingStateWeight(long id) const {
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

void HadronSelector::doinit() {
  Interfaced::doinit();
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
  // construct the hadron tables
  constructHadronTable();
  // lightest members (hadrons)
  for(const PDPtr & p1 : partons()) {
    for(const PDPtr & p2 : partons()) {
      tcPDPair lp = lightestHadronPair(p1,p2);
      if(lp.first && lp.second)
	lightestHadrons_[make_pair(p1->id(),p2->id())] = lp;
    }
  }
  // for debugging
  if(Debug::level >= 10 ) 
    dumpTable(table());
}

void HadronSelector::insertToHadronTable(tPDPtr &particle, int flav1, int flav2) {
  // inserting a new Hadron in the hadron table.
  long pid  = particle->id();
  int pspin = particle->iSpin();
  HadronInfo a(pid, particle,specialWeight(pid),particle->mass());
  // set the weight to the number of spin states
  a.overallWeight = pspin*a.swtef;
  // mesons
  if(pspin%2==1)     insertMeson(a,flav1,flav2);
  // spin-1/2 baryons
  else if(pspin==2) insertOneHalf(a,flav1,flav2);
  // spin -3/2 baryons
  else if(pspin==4) insertThreeHalf(a,flav1,flav2);
  // all other cases
  else {
    assert(false);
  }
}

void HadronSelector::constructHadronTable() {
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
    int flav1,flav2;
    // meson
    if(x4 == 0) {
      flav1 = x2;
      flav2 = x3;
      if (wantSusy) flav2 += 1000000 * x7;
    }
    // baryon
    else {
      flav2 = x4;
      if (wantSusy) flav2 += 1000000 * x7;
      // insert the spin 1 diquark, sort out the rest later
      flav1 = CheckId::makeDiquarkID(x2,x3,3);
    }
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

double HadronSelector::mesonWeight(long id) const {
  // Total angular momentum
  int j  = ((id % 10) - 1) / 2;
  // related to Orbital angular momentum l
  int nl = (id/10000 )%10;
  int l  = -999;
  int n  = (id/100000)%10;  // Radial excitation
  if(j == 0) l = nl;
  else if(nl == 0) l = j - 1;
  else if(nl == 1  || nl == 2) l = j;
  else if(nl == 3) l = j + 1;
  // Angular or Radial excited meson
  if((l||j||n) && l>=0  &&  l<Lmax  &&  j<Jmax  &&  n<Nmax) {
    return sqr(_repwt[l][j][n]);
  }
  // rest is not excited or
  // has spin >= 5/2 (ispin >= 6), haven't got those
  else
    return 1.0;
}


int HadronSelector::signHadron(tcPDPtr idQ1, tcPDPtr idQ2,
			       tcPDPtr hadron) const {
  // This method receives in input three PDG ids, whose the
  // first two have proper signs (corresponding to particles, id > 0,
  // or antiparticles, id < 0 ), whereas the third one must
  // be always positive (particle not antiparticle),
  // corresponding to:
  //  --- quark-antiquark, or antiquark-quark, or
  //      quark-diquark, or diquark-quark, or
  //      antiquark-antidiquark, or antidiquark-antiquark
  //      for the first two input (idQ1, idQ2);
  //  --- meson or baryon for the third input (idHad):
  // The method returns:
  //  --- + 1  if the two partons (idQ1, idQ2) are exactly
  //           the constituents for the hadron idHad;
  //  --- - 1  if the two partons (idQ1, idQ2) are exactly
  //           the constituents for the anti-hadron -idHad;
  //  --- + 0  otherwise.
  // The method it is therefore useful to decide the
  // sign of the id of the produced hadron as appeared
  // in the vector _vecHad (where only hadron idHad > 0 are present)
  // given the two constituent partons.
  int sign = 0;
  long idHad = hadron->id();
  assert(idHad > 0);
  int chargeIn  = idQ1->iCharge() + idQ2->iCharge();
  int chargeOut = hadron->iCharge();
  // same charge
  if(     chargeIn ==  chargeOut && chargeIn  !=0 ) sign = +1;
  else if(chargeIn == -chargeOut && chargeIn  !=0 ) sign = -1;
  else if(chargeIn == 0          && chargeOut == 0 ) {
    // In the case of same null charge, there are four cases:
    //  i) K0-like mesons, B0-like mesons, Bs-like mesons
    //     the PDG convention is to consider them "antiparticle" (idHad < 0)
    //     if the "dominant" (heavier) flavour (respectively, s, b)
    //     is a quark (idQ > 0): for instance, B0s = (b, sbar) has id < 0
    //     Remember that there is an important exception for K0L (id=130) and
    //     K0S (id=310): they don't have antiparticles, therefore idHad > 0
    //     always. We use below the fact that K0L and K0S are the unique
    //     hadrons having 0 the first (less significant) digit of their id.
    //  2) D0-like mesons: the PDG convention is to consider them "particle"
    //     (idHad > 0) if the charm flavour is carried by a c: (c,ubar) has id>0
    //  3) the remaining mesons should not have antiparticle, therefore their
    //     sign is always positive.
    //  4) for baryons, that is when one of idQ1 and idQ2 is a (anti-) quark and
    //     the other one is a (anti-) diquark the sign is negative when both
    //     constituents are "anti", that is both with id < 0; positive otherwise.
    // meson
    if(abs(int(idQ1->iColour()))== 3 && abs(int(idQ2->iColour())) == 3 &&
      !DiquarkMatcher::Check(idQ1->id()) && !DiquarkMatcher::Check(idQ2->id()))
    {
      int idQa = abs(idQ1->id());
      int idQb = abs(idQ2->id());
      int dominant = idQ2->id();

      if(idQa > idQb) {
	swap(idQa,idQb);
	dominant = idQ1->id();
      }

      if((idQa==ParticleID::d && idQb==ParticleID::s) ||
	 (idQa==ParticleID::d && idQb==ParticleID::b) ||
	 (idQa==ParticleID::s && idQb==ParticleID::b)) {
	// idHad%10 is zero for K0L,K0S
	if (dominant < 0 || idHad%10 == 0) sign = +1;
	else if(dominant > 0)              sign = -1;
      }
      else if((idQa==ParticleID::u && idQb==ParticleID::c) ||
	      (idQa==ParticleID::u && idQb==ParticleID::t) ||
	      (idQa==ParticleID::c && idQb==ParticleID::t)) {
	if     (dominant > 0) sign = +1;
	else if(dominant < 0) sign = -1;
      }
      else if(idQa==idQb) sign = +1;
      // sets sign for Susy particles
      else sign = (dominant > 0) ? +1 : -1;
    }
    // baryon
    else if(DiquarkMatcher::Check(idQ1->id()) || DiquarkMatcher::Check(idQ2->id())) {
      if     (idQ1->id() > 0 && idQ2->id() > 0) sign = +1;
      else if(idQ1->id() < 0 && idQ2->id() < 0) sign = -1;
    }
  }
  if (sign == 0) {
    cerr << "Could not work out sign for "
	 << idQ1->PDGName() << ' '
	 << idQ2->PDGName() << " => "
	 << hadron->PDGName() << '\n';
    assert(false);
  }
  return sign;
}

tcPDPair HadronSelector::lightestHadronPair(tcPDPtr ptr1, tcPDPtr ptr2) const {
  Energy currentSum = Constants::MaxEnergy;
  tcPDPair output;
  for(unsigned int ix=0; ix<partons().size(); ++ix) {
    HadronTable::const_iterator
      tit1=table().find(make_pair(abs(ptr1->id()),partons()[ix]->id())),
      tit2=table().find(make_pair(partons()[ix]->id(),abs(ptr2->id())));
    if( tit1==table().end() || tit2==table().end()) continue;
    if(tit1->second.empty()||tit2->second.empty()) continue;
    Energy s = tit1->second.begin()->mass + tit2->second.begin()->mass;
    if(currentSum > s) {
      currentSum = s;
      output.first  = tit1->second.begin()->ptrData;
      output.second = tit2->second.begin()->ptrData;
    }
  }
  return output;
}

tcPDPtr HadronSelector::lightestHadron(tcPDPtr ptr1, tcPDPtr ptr2) const {
  assert(ptr1 && ptr2);
  // find entry in the table
  pair<long,long> ids = make_pair(abs(ptr1->id()),abs(ptr2->id()));
  HadronTable::const_iterator tit=_table.find(ids);
  // throw exception if flavours wrong
  if (tit==_table.end())
    throw Exception() << "Could not find "
		      << ids.first << ' ' << ids.second
		      << " in _table. "
		      << "In HadronSelector::lightestHadron()"
		      << Exception::eventerror;
  if(tit->second.empty())
    throw Exception() << "HadronSelector::lightestHadron "
		      << "could not find any hadrons containing "
		      << ptr1->id() << ' ' << ptr2->id() << '\n'
		      << tit->first.first << ' '
		      << tit->first.second << Exception::eventerror;
  // find the lightest hadron
  int sign = signHadron(ptr1,ptr2,tit->second.begin()->ptrData);
  tcPDPtr candidate = sign > 0 ?
    tit->second.begin()->ptrData : tit->second.begin()->ptrData->CC();
  // \todo 20 GeV limit is temporary fudge to let SM particles go through.
  // \todo Use isExotic instead?
  if (candidate->mass() > 20*GeV
      && candidate->mass() < ptr1->constituentMass() + ptr2->constituentMass()) {
    generator()->log() << "HadronSelector::lightestHadron: "
		       << "chosen candidate " << candidate->PDGName()
		       << " is lighter than its constituents "
		       << ptr1->PDGName() << ", " << ptr2->PDGName() << '\n'
		       << candidate->mass()/GeV << " < " << ptr1->constituentMass()/GeV
		       << " + " << ptr2->constituentMass()/GeV << '\n'
		       << "Check your particle data tables.\n";
    assert(false);
  }
  return candidate;
}

vector<pair<tcPDPtr,double> >
HadronSelector::hadronsBelowThreshold(Energy threshold, tcPDPtr ptr1,
				      tcPDPtr ptr2) const {
  // The method assumes ptr3 == 0 rest not implemented
  assert(ptr1 && ptr2);
  // find entry in the table
  pair<long,long> ids = make_pair(abs(ptr1->id()),abs(ptr2->id()));
  HadronTable::const_iterator tit=_table.find(ids);
  // throw exception if flavours wrong
  if (tit==_table.end())
    throw Exception() << "Could not find "
		      << ids.first << ' ' << ids.second
		      << " in _table. "
		      << "In HadronSelector::hadronsBelowThreshold()"
		      << Exception::eventerror;
  if(tit->second.empty())
    throw Exception() << "HadronSelector::hadronsBelowThreshold() "
		      << "could not find any hadrons containing "
		      << ptr1->id() << ' ' << ptr2->id() << '\n'
		      << tit->first.first << ' '
		      << tit->first.second << Exception::eventerror;
  vector<pair<tcPDPtr,double> > candidates;
  KupcoData::const_iterator hit = tit->second.begin();
  // find the hadrons
  while(hit!=tit->second.end()&&hit->mass<threshold) {
    // find the hadron
    int sign = signHadron(ptr1,ptr2,hit->ptrData);
    tcPDPtr candidate = sign > 0 ? hit->ptrData : hit->ptrData->CC();
    // \todo 20 GeV limit is temporary fudge to let SM particles go through.
    // \todo Use isExotic instead?
    if (candidate->mass() > 20*GeV
	&& candidate->mass() < ptr1->constituentMass() + ptr2->constituentMass()) {
      generator()->log() << "HadronSelector::hadronsBelowTheshold: "
			 << "chosen candidate " << candidate->PDGName()
			 << " is lighter than its constituents "
			 << ptr1->PDGName() << ", " << ptr2->PDGName() << '\n'
			 << candidate->mass()/GeV << " < " << ptr1->constituentMass()/GeV
			 << " + " << ptr2->constituentMass()/GeV << '\n'
			 << "Check your particle data tables.\n";
      assert(false);
    }
    candidates.push_back(make_pair(candidate,hit->overallWeight));
    ++hit;
  }
  return candidates;
}


tcPDPtr HadronSelector::chooseSingleHadron(tcPDPtr par1, tcPDPtr par2,
					   Energy mass) const {
  // Determine the sum of the nominal masses of the two lightest hadrons
  // with the right flavour numbers as the cluster under consideration.
  // Notice that we don't need real masses (drawn by a Breit-Wigner
  // distribution) because the lightest pair of hadrons does not involve
  // any broad resonance.
  Energy threshold = massLightestHadronPair(par1,par2);
  // Special: it allows one-hadron decays also above threshold.
  if (CheckId::isExotic(par1,par2))
    threshold *= (1.0 + UseRandom::rnd()*_limExotic);
  else if (CheckId::hasBottom(par1,par2))
    threshold *= (1.0 + UseRandom::rnd()*_limBottom);
  else if (CheckId::hasCharm(par1,par2))
    threshold *= (1.0 + UseRandom::rnd()*_limCharm);

  // only do one hadron decay is mass less than the threshold
  if(mass>=threshold) return tcPDPtr();

  // select the hadron
  tcPDPtr hadron;
  // old option pick the lightest hadron
  if(belowThreshold_ == 0) {
    hadron= lightestHadron(par1,par2);
  }
  // new option select from those available
  else if(belowThreshold_ == 1) {
    vector<pair<tcPDPtr,double> > hadrons =
      hadronsBelowThreshold(threshold,par1,par2);
    if(hadrons.size()==1) {
      hadron = hadrons[0].first;
    }
    else if(hadrons.empty()) {
      hadron= lightestHadron(par1,par2);
    }
    else {
      double totalWeight=0.;
      for(unsigned int ix=0;ix<hadrons.size();++ix) {
	totalWeight += hadrons[ix].second;
      }
      totalWeight *= UseRandom::rnd();
      for(unsigned int ix=0;ix<hadrons.size();++ix) {
	if(totalWeight<=hadrons[ix].second) {
	  hadron = hadrons[ix].first;
	  break;
	}
	else
	  totalWeight -= hadrons[ix].second;
      }
      assert(hadron);
    }
  }
  else
    assert(false);
  return hadron;
}

tcPDPair HadronSelector::chooseHadronPair(const Energy cluMass,
					  tcPDPtr par1, tcPDPtr par2) const {
  useMe();
  // if either of the input partons is a diquark don't allow diquarks to be
  // produced
  bool diquark0 = !(DiquarkMatcher::Check(par1->id()) || DiquarkMatcher::Check(par2->id()));
  bool diquark1 = diquark0;
  bool quark = true;
  // decide is baryon or meson production
  if(diquark0) std::tie(quark,diquark0,diquark1) = selectBaryon(cluMass,par1,par2);
  // weights for the different possibilities
  Energy weight, wgtsum(ZERO);
  // loop over all hadron pairs with the allowed flavours
  static vector<Kupco> hadrons;
  hadrons.clear();
  for(unsigned int ix=0;ix<partons().size();++ix) {
    tcPDPtr quarktopick  = partons()[ix];
    if(!quark  &&  abs(int(quarktopick->iColour())) == 3
       && !DiquarkMatcher::Check(quarktopick->id())) continue;
    if(abs(int(quarktopick->iColour())) == 3
       && DiquarkMatcher::Check(quarktopick->id()) &&
       ((!diquark0 && quarktopick->iSpin()==1) ||
	(!diquark1 && quarktopick->iSpin()==3))) continue;
    HadronTable::const_iterator
      tit1 = table().find(make_pair(abs(par1->id()),quarktopick->id()));
    HadronTable::const_iterator
      tit2 = table().find(make_pair(quarktopick->id(),abs(par2->id())));
    // If not in table skip
    if(tit1 == table().end()||tit2==table().end()) continue;
    // tables empty skip
    const KupcoData & T1 = tit1->second;
    const KupcoData & T2 = tit2->second;
    if(T1.empty()||T2.empty()) continue;
    // if too massive skip
    if(cluMass <= T1.begin()->mass +
                  T2.begin()->mass) continue;
    // quark weight
    double quarkWeight =  pwt(quarktopick->id());
    // special for strange
    if(abs(quarktopick->id()) == 3)
      quarkWeight = strangeWeight(cluMass,par1,par2);
    // loop over the hadrons
    KupcoData::const_iterator H1,H2;
    for(H1 = T1.begin();H1 != T1.end(); ++H1) {
      for(H2 = T2.begin();H2 != T2.end(); ++H2) {
 	// break if cluster too light
 	if(cluMass < H1->mass + H2->mass) break;
	weight = quarkWeight * H1->overallWeight * H2->overallWeight *
	  Kinematics::pstarTwoBodyDecay(cluMass, H1->mass, H2->mass);
	int signQ = 0;
	assert (par1 && quarktopick);
	assert (par2);

	assert(quarktopick->CC());

	if(CheckId::canBeHadron(par1, quarktopick->CC())
	   && CheckId::canBeHadron(quarktopick, par2))
	   signQ = +1;
	else if(CheckId::canBeHadron(par1, quarktopick)
		&& CheckId::canBeHadron(quarktopick->CC(), par2))
	   signQ = -1;
	else {
	  cerr << "Could not make sign for" << par1->id()<< " " << quarktopick->id()
	       << " " << par2->id() << "\n";
	  assert(false);
	}

	if (signQ  == -1)
	  quarktopick = quarktopick->CC();
	// construct the object with the info
	Kupco a(quarktopick, H1->ptrData, H2->ptrData, weight);
	hadrons.push_back(a);
	wgtsum += weight;
      }
    }
  }
  if (hadrons.empty())
    return make_pair(tcPDPtr(),tcPDPtr());
  // select the hadron
  wgtsum *= UseRandom::rnd();
  unsigned int ix=0;
  do {
    wgtsum-= hadrons[ix].weight;
    ++ix;
  }
  while(wgtsum > ZERO && ix < hadrons.size());
  if(ix == hadrons.size() && wgtsum > ZERO)
      return make_pair(tcPDPtr(),tcPDPtr());
  --ix;
  assert(hadrons[ix].idQ);
  int signHad1 = signHadron(par1, hadrons[ix].idQ->CC(), hadrons[ix].hadron1);
  int signHad2 = signHadron(par2, hadrons[ix].idQ, hadrons[ix].hadron2);
  assert( signHad1 != 0 && signHad2 != 0 );
  return make_pair
    ( signHad1 > 0 ? hadrons[ix].hadron1 : tcPDPtr(hadrons[ix].hadron1->CC()),
      signHad2 > 0 ? hadrons[ix].hadron2 : tcPDPtr(hadrons[ix].hadron2->CC()));
}

std::tuple<bool,bool,bool> HadronSelector::selectBaryon(const Energy, tcPDPtr, tcPDPtr )  const {
  assert(false);
}

double HadronSelector::strangeWeight(const Energy, tcPDPtr, tcPDPtr) const {
  assert(false);
}

void HadronSelector::insertMeson(HadronInfo a, int flav1, int flav2) {
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

void HadronSelector::insertOneHalf(HadronInfo a, int flav1, int flav2) {
  assert(DiquarkMatcher::Check(flav1));
  long iq1 = flav1/1000;
  long iq2 = (flav1/100)%10;
  if(iq1!=iq2 && flav1%10==3) flav1-=2;
  if(iq1==iq2) {
    if(iq1==flav2) {
      a.overallWeight *= 1.5;
      _table[make_pair(flav1,flav2)].insert(a);
      _table[make_pair(flav2,flav1)].insert(a);
    }
    else {
      table()[make_pair(flav1,flav2)].insert(a);
      table()[make_pair(flav2,flav1)].insert(a);
      long f3 = CheckId::makeDiquarkID(iq1,flav2,1);
      table()[make_pair(iq1,f3 )].insert(a);
      table()[make_pair(f3 ,iq1)].insert(a);
    }
  }
  else if(iq1==flav2) {
    // ud1 u type
    table()[make_pair(flav1,flav2)].insert(a);
    table()[make_pair(flav2,flav1)].insert(a);
    // and uu1 d type
    long f3 = CheckId::makeDiquarkID(iq1,iq1,3);
    a.overallWeight *= a.wt;
    table()[make_pair(f3 ,iq2)].insert(a);
    table()[make_pair(iq2, f3)].insert(a);
  }
  else if(iq2==flav2) assert(false);
  else {
    table()[make_pair(flav1,flav2)].insert(a);
    table()[make_pair(flav2,flav1)].insert(a);
    long f3 = CheckId::makeDiquarkID(iq1,flav2,1);
    table()[make_pair(iq2,f3)].insert(a);
    table()[make_pair(f3,iq2)].insert(a);
    // 3rd perm
    f3 = CheckId::makeDiquarkID(iq2,flav2,1);
    table()[make_pair(iq1,f3)].insert(a);
    table()[make_pair(f3,iq1)].insert(a);
  }
}

void HadronSelector::insertThreeHalf(HadronInfo a, int flav1, int flav2) {
  assert(DiquarkMatcher::Check(flav1));
  long iq1 = flav1/1000;
  long iq2 = (flav1/100)%10;
  if(iq1!=iq2 && flav1%10==3) flav1-=2;
  if(iq1==iq2) {
    if(iq1==flav2) {
      a.overallWeight *= 1.5;
      _table[make_pair(flav1,flav2)].insert(a);
      _table[make_pair(flav2,flav1)].insert(a);
    }
    else {
      table()[make_pair(flav1,flav2)].insert(a);
      table()[make_pair(flav2,flav1)].insert(a);
      long f3 = CheckId::makeDiquarkID(iq1,flav2,1);
      table()[make_pair(iq1,f3 )].insert(a);
      table()[make_pair(f3 ,iq1)].insert(a);
    }
  }
  else if(iq1==flav2) {
    // ud1 u type
    table()[make_pair(flav1,flav2)].insert(a);
    table()[make_pair(flav2,flav1)].insert(a);
    // and uu1 d type
    long f3 = CheckId::makeDiquarkID(iq1,iq1,3);
    a.overallWeight *= a.wt;
    table()[make_pair(f3 ,iq2)].insert(a);
    table()[make_pair(iq2, f3)].insert(a);
  }
  else {
    table()[make_pair(flav1,flav2)].insert(a);
    table()[make_pair(flav2,flav1)].insert(a);
    long f3 = CheckId::makeDiquarkID(iq1,flav2,1);
    table()[make_pair(iq2,f3)].insert(a);
    table()[make_pair(f3,iq2)].insert(a);
    // 3rd perm
    f3 = CheckId::makeDiquarkID(iq2,flav2,1);
    table()[make_pair(iq1,f3)].insert(a);
    table()[make_pair(f3,iq1)].insert(a);
  }
}

PDPtr HadronSelector::makeDiquark(tcPDPtr par1, tcPDPtr par2) {
  long id1 = par1->id();
  long id2 = par2->id();
  long pspin = id1==id2 ? 3 : 1;
  long idnew = CheckId::makeDiquarkID(id1,id2, pspin);
  return getParticleData(idnew);
}
