// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DarkHadronSpectrum class.
//

#include "DarkHadronSpectrum.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ParMap.h"

#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/CurrentGenerator.h>
#include <ThePEG/Repository/Repository.h>


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

namespace {
  bool weightIsLess (pair<long,double> a, pair<long,double> b) {
    return a.second < b.second;
  }
}


DarkHadronSpectrum::DarkHadronSpectrum(unsigned int opt) 
  : HadronSpectrum(),
    _topt(opt),_trial(0),
    _nlightquarks(2),
    _nheavyquarks(0) {}


DarkHadronSpectrum::~DarkHadronSpectrum() {}


void DarkHadronSpectrum::persistentOutput(PersistentOStream & os) const {
  os << _sngWt << _decWt << _repwt << _pwtDIquark
     << _nlightquarks << _nheavyquarks
     << belowThreshold_;
}

void DarkHadronSpectrum::persistentInput(PersistentIStream & is, int) {
  is >> _sngWt >> _decWt >> _repwt >> _pwtDIquark
     >> _nlightquarks >> _nheavyquarks
     >> belowThreshold_;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<DarkHadronSpectrum,HadronSpectrum>
describeHerwigDarkHadronSpectrum("Herwig::DarkHadronSpectrum", "Herwig.so");

void DarkHadronSpectrum::Init() {

  static ClassDocumentation<DarkHadronSpectrum> documentation
    ("There is no documentation for the DarkHadronSpectrum class");

  static ParMap<DarkHadronSpectrum,double>
    interfacePwt("Pwt","Weights for choosing the quarks",
		       &DarkHadronSpectrum::_pwt, 0, -1, 1.0, 0.0, 10.0, 
     false, false, false);

  static Parameter<DarkHadronSpectrum,double>
    interfacePwtDIquark("PwtDIquark","Weight for choosing a DIquark",
			&DarkHadronSpectrum::_pwtDIquark, 0, 1.0, 0.0, 100.0,
			false,false,false);

  static Parameter<DarkHadronSpectrum,int> interfaceNLightQuarks
    ("NumLightQuarks",
     "The number of light quarks (roughly mass degenerate)",
     &DarkHadronSpectrum::_nlightquarks, 0, 2, 0, 9, false,false,false);

  static Parameter<DarkHadronSpectrum,int> interfaceNHeavyQuarks
    ("NumHeavyQuarks",
     "The number of heavy quarks (non-mass degenerate)",
     &DarkHadronSpectrum::_nheavyquarks, 0, 0, 0, 9, false,false,false);
  //
  // mixing angles
  //
  // the ideal mixing angle

  /*
  //  the meson weights
  //
  static ParVector<DarkHadronSpectrum,double> interface1S0Weights
    ("1S0Weights",
     "The weights for the 1S0 multiplets start with n=1.",
     &DarkHadronSpectrum::_weight1S0, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static ParVector<DarkHadronSpectrum,double> interface3S1Weights
    ("3S1Weights",
     "The weights for the 3S1 multiplets start with n=1.",
     &DarkHadronSpectrum::_weight3S1, Nmax, 1.0, 0.0, 100.0,
     false, false, Interface::limited);
  */

  static Switch<DarkHadronSpectrum,unsigned int> interfaceTrial
    ("Trial",
     "A Debugging option to only produce certain types of hadrons",
     &DarkHadronSpectrum::_trial, 0, false, false);
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

  static Switch<DarkHadronSpectrum,unsigned int> interfaceBelowThreshold
    ("BelowThreshold",
     "Option for the selection of the hadrons if the cluster is below the pair threshold",
     &DarkHadronSpectrum::belowThreshold_, 0, false, false);
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

Energy DarkHadronSpectrum::hadronPairThreshold(tcPDPtr par1, tcPDPtr par2) const {
  // Determine the sum of the nominal masses of the two lightest hadrons
  // with the right flavour numbers as the cluster under consideration.
  // Notice that we don't need real masses (drawn by a Breit-Wigner 
  // distribution) because the lightest pair of hadrons does not involve
  // any broad resonance.
  return massLightestHadronPair(par1,par2);
}
  

double DarkHadronSpectrum::mixingStateWeight(long id) const {
  switch(id) {
  //case ParticleID::eta:      return 0.5*probabilityMixing(_etamix  ,1);
  //case ParticleID::etaprime: return 0.5*probabilityMixing(_etamix  ,2);
  default:                   return 1./3.;
  }
}

void DarkHadronSpectrum::doinit() {
  if (_nlightquarks + _nheavyquarks > 9) {
    throw InitException() << "Can have a maximum of 9 dark quarks!";
  }

  // the default partons allowed
  // the quarks
  for (long ix : hadronizingQuarks()) {
    _partons.push_back(getParticleData(ix));
  }
  // the diquarks
  for(long ix : hadronizingQuarks()) {
    for(long iy : hadronizingQuarks()) {
        // Only add diquarks if they've been defined
      if(ix==iy) {
	if (getParticleData(makeDiquarkID(ix,iy,long(3))))
	  _partons.push_back(getParticleData(makeDiquarkID(ix,iy,long(3))));
      } else {
	if (getParticleData(makeDiquarkID(ix,iy,long(1))))
	  _partons.push_back(getParticleData(makeDiquarkID(ix,iy,long(1))));
      }
    }
  }
  // set the weights for the various excited mesons
  // set all to one to start with
  for (int l = 0; l < Lmax; ++l ) {
    for (int j = 0; j < Jmax; ++j) {
      for (int n = 0; n < Nmax; ++n) {
	_repwt[l][j][n] = 1.0;
      }
    }
  }
  // weights for the different quarks etc

  // NB: Assuming max 9 quarks, first digit used in tagging hadrons and diquarks!
  vector<long> light = lightHadronizingQuarks();
  for(unsigned int ix=0; ix<light.size(); ++ix){
    for(unsigned int iy=0; iy<ix; ++iy){
      _pwt[_DarkHadOffset + (light[ix]%10)*1000 + (light[iy]%10)*100 + 1] =
        0.5 * _pwtDIquark * _pwt[light[ix]] * _pwt[light[iy]];
    }
    _pwt[_DarkHadOffset + (light[ix]%10)*1100 + 3] =
      _pwtDIquark * _pwt[light[ix]] * _pwt[light[ix]];
  }

  // Set any other weights to zero
  // NB: Only producing light diquarks!
  for(unsigned int ix=0; ix<_partons.size(); ++ix) {
    if (_pwt.find(_partons[ix]->id()) == _pwt.end()){
        _pwt[_partons[ix]->id()]=0.;
    }
  }

  // find the maximum
  map<long,double>::iterator pit =
    max_element(_pwt.begin(),_pwt.end(),weightIsLess); 
  const double pmax = pit->second;
  for(pit=_pwt.begin(); pit!=_pwt.end(); ++pit) {
    pit->second/=pmax;
  }
  HadronSpectrum::doinit();
}

void DarkHadronSpectrum::constructHadronTable() {
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
    // Require dark hadrons
    if ((pid/100000)!=49) continue;
    // Get the flavours
    const int x4 = (pid/1000)%10; 
    const int x3 = (pid/100 )%10;
    const int x2 = (pid/10  )%10;
    int flav1;
    int flav2;
    // Skip non-hadrons (susy particles, etc...)
    if(x3 == 0 || x2 == 0) continue;
    else if(x4 == 0) { // meson
      flav1 = x2; 
      flav2 = x3; 
    } 
    else { // baryon
      // insert the spin 1 diquark, sort out the rest later
      flav1 = makeDiquarkID(x2,x3,3);
      flav2 = x4;
    }
    insertToHadronTable(particle,flav1,flav2);
  }

  // normalise weights to one for first option
  if(_topt == 0) {
    HadronTable::const_iterator tit;
    KupcoData::iterator it;
    for(tit=_table.begin();tit!=_table.end();++tit) {
      double weight=0;
      for(it = tit->second.begin(); it!=tit->second.end(); ++it)
	weight=max(weight,it->overallWeight);
      weight = 1./weight;
    }
  }
}

void DarkHadronSpectrum::insertMeson(HadronInfo a, int flav1, int flav2) {
  // identical light flavours
  if(flav1 == flav2 && flav1<=_nlightquarks) {
    vector<long> light = lightHadronizingQuarks();
    // light quark admixture states
    a.overallWeight *= 1. / light.size();
    for(unsigned int ix=0; ix<light.size(); ++ix){
      _table[make_pair(light[ix],light[ix])].insert(a);
	}
  }
  else {
    _table[make_pair(_DarkHadOffset+flav1, _DarkHadOffset+flav2)].insert(a);
    if(flav1 != flav2) _table[make_pair(_DarkHadOffset+flav2, _DarkHadOffset+flav1)].insert(a);
  }
}


double DarkHadronSpectrum::mesonWeight(long) const {
  // Don't currently have radial excitations (practically clashes with dark had offset
  // in pdgId codes; theoretically doesn't make much sense to consider such complex
  // states). For now just return 1
  return 1.0;
}


long DarkHadronSpectrum::makeDiquarkID(long id1, long id2, long pspin) const {

  assert(id1 * id2 > 0);
  // Currently not checking if these are dark quarks to allow giving either ID
  // or 1 digit flav, is this a good idea long term?
  //       && DarkQuarkMatcher::Check(id1)  
  //       && DarkQuarkMatcher::Check(id2)) ;
  long ida = abs(id1)%10;
  long idb = abs(id2)%10;
  if (ida < idb) swap(ida,idb);
  if (pspin != 1 && pspin != 3) assert(false);
  long idnew = ida*1000 + idb*100 + pspin + _DarkHadOffset;
  // Diquarks made of quarks of the same type: uu, dd, ss, cc, bb, 
  // have spin 1, and therefore the less significant digit (which
  // corresponds to 2*J+1) is 3 rather than 1 as all other Diquarks.
  if (id1 == id2 && pspin == 1) {
    //cerr<<"WARNING: spin-0 diquiark of the same type cannot exist."
    //    <<" Switching to spin-1 diquark.\n";
    idnew = ida*1000 + idb*100 + 3 + _DarkHadOffset;
  }

  return id1 > 0 ? idnew : -idnew;
}

bool DarkHadronSpectrum::isExotic(tcPDPtr, tcPDPtr, tcPDPtr) const {
  // Don't list dark particles as exotic to allow them to be treated as either light
  // or heavy in the hadronisation
  return false;
}

bool DarkHadronSpectrum::canBeBaryon(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3) const {
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
