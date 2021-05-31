// -*- C++ -*-
//
// Hw7Selector.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Hw7Selector class.
//

#include "Hw7Selector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Utilities/Selector.h"
#include "ThePEG/Repository/UseRandom.h"
#include "CheckId.h"
#include <cassert>
#include <ThePEG/Utilities/DescribeClass.h>

using namespace Herwig;

DescribeClass<Hw7Selector,HadronSelector>
describeHw7Selector("Herwig::Hw7Selector","Herwig.so");

IBPtr Hw7Selector::clone() const {
  return new_ptr(*this);
}

IBPtr Hw7Selector::fullclone() const {
  return new_ptr(*this);
}

void Hw7Selector::doinit() {
  // the default partons allowed
  // the quarks
  for ( int ix=1; ix<=5; ++ix ) {
    partons().push_back(getParticleData(ix));
  }
  // the diquarks
  for(unsigned int ix=1;ix<=5;++ix) {
    for(unsigned int iy=1; iy<=ix;++iy) {
      partons().push_back(getParticleData(CheckId::makeDiquarkID(ix,iy,long(3))));
      if(ix!=iy)
        partons().push_back(getParticleData(CheckId::makeDiquarkID(ix,iy,long(1))));
    }
  }
  // weights for the different quarks etc
  for(unsigned int ix=0; ix<partons().size(); ++ix) {
    pwt()[partons()[ix]->id()]=0.;
  }
  pwt()[1]  = _pwtDquark;
  pwt()[2]  = _pwtUquark;
  pwt()[3]  = _pwtSquark;
  pwt()[4]  = _pwtCquark;
  pwt()[5]  = _pwtBquark;
  pwt()[1103] = _pwtDIquarkS1 * _pwtDquark * _pwtDquark;
  pwt()[2101] = _pwtDIquarkS0 * _pwtUquark * _pwtDquark;
  pwt()[2103] = _pwtDIquarkS1 * _pwtUquark * _pwtDquark;
  pwt()[2203] = _pwtDIquarkS1 * _pwtUquark * _pwtUquark;
  pwt()[3101] = _pwtDIquarkS0 * _pwtSquark * _pwtDquark;
  pwt()[3103] = _pwtDIquarkS1 * _pwtSquark * _pwtDquark;
  pwt()[3201] = _pwtDIquarkS0 * _pwtSquark * _pwtUquark;
  pwt()[3203] = _pwtDIquarkS1 * _pwtSquark * _pwtUquark;
  pwt()[3303] = _pwtDIquarkS1 * _pwtSquark * _pwtSquark;
  HadronSelector::doinit();
  // lightest members (baryons)
  for(const PDPtr & p1 : partons()) {
    if(DiquarkMatcher::Check(p1->id())) continue;
    for(const PDPtr & p2 : partons()) {
      if(DiquarkMatcher::Check(p2->id())) continue;
      lightestBaryonsS0_[make_pair(p1->id(),p2->id())] = lightestBaryonPair(p1,p2,1);
      lightestBaryonsS1_[make_pair(p1->id(),p2->id())] = lightestBaryonPair(p1,p2,3);
    }
  }
}

void Hw7Selector::persistentOutput(PersistentOStream & os) const {
  os << _pwtDquark  << _pwtUquark << _pwtSquark
     << _pwtCquark << _pwtBquark << _pwtDIquarkS0 << _pwtDIquarkS1
     << _sngWt << _decWt 
     << _mode << _enhanceSProb << ounit(_m0Decay,GeV) << _massMeasure
     << _scHadronWtFactor << _sbHadronWtFactor
     << lightestBaryonsS0_ << lightestBaryonsS1_;
}

void Hw7Selector::persistentInput(PersistentIStream & is, int) {
  is >> _pwtDquark  >> _pwtUquark >> _pwtSquark
     >> _pwtCquark >> _pwtBquark >> _pwtDIquarkS0 >> _pwtDIquarkS1
     >> _sngWt >> _decWt 
     >> _mode >> _enhanceSProb >> iunit(_m0Decay,GeV) >> _massMeasure
     >> _scHadronWtFactor >> _sbHadronWtFactor
     >> lightestBaryonsS0_ >> lightestBaryonsS1_;
}

void Hw7Selector::Init() {

  static ClassDocumentation<Hw7Selector> documentation
    ("The Hw7Selector class implements the Herwig algorithm for selecting"
     " the hadrons",
     "The hadronization used the selection algorithm described in \\cite{Kupco:1998fx}.",
     "%\\cite{Kupco:1998fx}\n"
     "\\bibitem{Kupco:1998fx}\n"
     "  A.~Kupco,\n"
     "  ``Cluster hadronization in HERWIG 5.9,''\n"
     "  arXiv:hep-ph/9906412.\n"
     "  %%CITATION = HEP-PH/9906412;%%\n"
     );
  
  static Parameter<Hw7Selector,double>
    interfacePwtDquark("PwtDquark","Weight for choosing a quark D",
		       &Hw7Selector::_pwtDquark, 0, 1.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<Hw7Selector,double>
    interfacePwtUquark("PwtUquark","Weight for choosing a quark U",
		       &Hw7Selector::_pwtUquark, 0, 1.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<Hw7Selector,double>
    interfacePwtSquark("PwtSquark","Weight for choosing a quark S",
		       &Hw7Selector::_pwtSquark, 0, 1.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<Hw7Selector,double>
    interfacePwtCquark("PwtCquark","Weight for choosing a quark C",
		       &Hw7Selector::_pwtCquark, 0, 0.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<Hw7Selector,double>
    interfacePwtBquark("PwtBquark","Weight for choosing a quark B",
		       &Hw7Selector::_pwtBquark, 0, 0.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<Hw7Selector,double>
    interfacePwtDIquarkS0("PwtDIquarkS0","Weight for choosing a spin-0 DIquark",
			&Hw7Selector::_pwtDIquarkS0, 0, 1.0, 0.0, 100.0,
			false,false,false);

  static Parameter<Hw7Selector,double>
    interfacePwtDIquarkS1("PwtDIquarkS1","Weight for choosing a spin-1 DIquark",
      &Hw7Selector::_pwtDIquarkS1, 0, 1.0, 0.0, 100.0,
    	false,false,false);

  static Parameter<Hw7Selector,double>
    interfaceSngWt("SngWt","Weight for singlet baryons",
                  &Hw7Selector::_sngWt, 0, 1.0, 0.0, 10.0,
		   false,false,false);

  static Parameter<Hw7Selector,double>
    interfaceDecWt("DecWt","Weight for decuplet baryons",
                  &Hw7Selector::_decWt, 0, 1.0, 0.0, 10.0,
		   false,false,false);
  
  static Switch<Hw7Selector,unsigned int> interfaceMode
    ("Mode",
     "Which algorithm to use",
     &Hw7Selector::_mode, 1, false, false);

  static SwitchOption interfaceModeKupco
    (interfaceMode,
     "Kupco",
     "Use the Kupco approach",
     0);

  static SwitchOption interfaceModeHw7
    (interfaceMode,
     "Hw7",
     "Use the Herwig approach",
     1);

  static Switch<Hw7Selector,int> interfaceEnhanceSProb
    ("EnhanceSProb",
     "Option for enhancing strangeness",
     &Hw7Selector::_enhanceSProb, 0, false, false);

  static SwitchOption interfaceEnhanceSProbNo
    (interfaceEnhanceSProb,
     "No",
     "No strangeness enhancement.",
     0);

  static SwitchOption interfaceEnhanceSProbScaled
    (interfaceEnhanceSProb,
     "Scaled",
     "Scaled strangeness enhancement",
     1);

  static SwitchOption interfaceEnhanceSProbExponential
    (interfaceEnhanceSProb,
     "Exponential",
     "Exponential strangeness enhancement",
     2);

   static Switch<Hw7Selector,int> interfaceMassMeasure
     ("MassMeasure",
      "Option to use different mass measures",
      &Hw7Selector::_massMeasure,0,false,false);

   static SwitchOption interfaceMassMeasureMass
     (interfaceMassMeasure,
      "Mass",
      "Mass Measure",
      0);

   static SwitchOption interfaceMassMeasureLambda
     (interfaceMassMeasure,
      "Lambda",
      "Lambda Measure",
      1);

   static Parameter<Hw7Selector,double> interfacescHadronWtFactor
     ("scHadronWtFactor",
     "Wight factor for strenge-charm heavy hadrns",
     &Hw7Selector::_scHadronWtFactor, 1., 0., 10.,
     false, false, Interface::limited);

   static Parameter<Hw7Selector,double> interfacesbHadronWtFactor
     ("sbHadronWtFactor",
     "Wight factor for strenge-bottom heavy hadrns",
     &Hw7Selector::_sbHadronWtFactor, 1., 0., 10.,
     false, false, Interface::limited);

  static Parameter<Hw7Selector,Energy> interfaceDecayMassScale
    ("DecayMassScale",
     "Cluster decay mass scale",
     &Hw7Selector::_m0Decay, GeV, 1.0*GeV, 0.1*GeV, 50.*GeV,
     false, false, Interface::limited);

}

double Hw7Selector::baryonWeight(long id) const {
  const int pspin = id % 10;
  if(pspin == 2) {
    // Singlet (Lambda-like) baryon
    if( (id/100)%10 < (id/10 )%10 ) return sqr(_sngWt);
  }
  // Decuplet baryon
  else if (pspin == 4)              return sqr(_decWt);
  return 1.;
}

namespace {

double kinFunction(Energy m0, tcPDPair pd) {
  Energy m1=pd.first ->mass();
  Energy m2=pd.second->mass();
  return pow((1.-sqr((m1+m2)/m0))*(1.-sqr((m1-m2)/m0)),0.25);
}
  
}

std::tuple<bool,bool,bool> Hw7Selector::selectBaryon(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const {
  useMe();
  std::tuple<bool,bool,bool> output(true,true,true);
  if(_mode ==1) {
    pair<long,long> ids(abs(par1->id()),abs(par2->id()));
    map<pair<long,long>,tcPDPair>::const_iterator lightest = lightestBaryonsS0_.find(ids);
    assert(lightest!=lightestBaryonsS0_.end());
    tcPDPair lightS0 = lightest->second;
    lightest = lightestBaryonsS1_.find(ids);
    assert(lightest!=lightestBaryonsS1_.end());
    tcPDPair lightS1 = lightest->second;
    Energy thresS0 = lightS0.first->mass()+lightS0.second->mass();
    Energy thresS1 = lightS1.first->mass()+lightS1.second->mass();
    // baryons not possible
    if(cluMass<=thresS0&&cluMass<=thresS1) {
      std::get<1>(output) = false;
      std::get<2>(output) = false;
    }
    // all baryons possible
    else if(cluMass>thresS0&&cluMass>thresS1) {
      double p0 = _pwtDIquarkS0*kinFunction(cluMass,lightS0);
      double p1 = _pwtDIquarkS1*kinFunction(cluMass,lightS1);
      double test=UseRandom::rnd()*(1.+p0+p1);
      if(test<=1.) {
	std::get<1>(output) = false;
	std::get<2>(output) = false;
      }
      else if(test<=1.+p0) {
	std::get<0>(output) = false;
	std::get<2>(output) = false;
      }
      else {
	std::get<0>(output) = false;
	std::get<1>(output) = false;
      }
    }
    // just spin-0 diquarks
    else if(cluMass>thresS0) {
      std::get<2>(output) = false;
      double p0 = _pwtDIquarkS0*kinFunction(cluMass,lightS0);
      double test=UseRandom::rnd()*(1.+p0);
      if(test<=1.) {
	std::get<1>(output) = false;
      }
      else {
	std::get<0>(output) = false;
      }
    }
    // just spin-1 diquarks
    else if(cluMass>thresS1) {
      std::get<1>(output) = false;
      double p1 = _pwtDIquarkS1*kinFunction(cluMass,lightS1);
      double test=UseRandom::rnd()*(1.+p1);
      if(test<=1.) {
	std::get<2>(output) = false;
      }
      else {
	std::get<0>(output) = false;
      }
    }
    else
      assert(false);
  }
  return output;
}

double Hw7Selector::strangeWeight(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const {
  // Decoupling the weight of heavy strenge hadrons
  if(_enhanceSProb == 0 && abs(par1->id()) == 4) {
    return pwt(3)*_scHadronWtFactor;
  }
  else if(_enhanceSProb == 0 && abs(par1->id()) == 5) {
    return pwt(3)*_sbHadronWtFactor;
  }
  // Scaling strangeness enhancement
  else if(_enhanceSProb == 1) {
    double scale = double(sqr(_m0Decay/cluMass));
    return (_maxScale < scale) ? 0. : pow(pwt(3),scale);
  }
  // Exponential strangeness enhancement
  else if(_enhanceSProb == 2) {
    Energy2 mass2;
    Energy endpointmass = par1->mass() + par2->mass();
    // Choose to use either the cluster mass
    // or to use the lambda measure
    mass2 = (_massMeasure == 0) ? sqr(cluMass) :
      sqr(cluMass) - sqr(endpointmass);
    double scale = double(sqr(_m0Decay)/mass2);
    return (_maxScale < scale) ? 0. : exp(-scale);
  }
  return pwt(3);
}

void Hw7Selector::insertOneHalf(HadronInfo a, int flav1, int flav2) {
  assert(DiquarkMatcher::Check(flav1));
  long iq1 = flav1/1000;
  long iq2 = (flav1/100)%10;
  if(iq1==iq2) {
    assert(flav2!=iq1);
    // first the uu1 d type piece
    a.wt = 1./3.;
    a.overallWeight *= a.wt;
    table()[make_pair(flav1,flav2)].insert(a);
    table()[make_pair(flav2,flav1)].insert(a);
    // also need ud1 u type
    long f3 = CheckId::makeDiquarkID(iq1,flav2,3);
    a.overallWeight /= a.wt;
    a.wt = 1./6.;
    a.overallWeight *= a.wt;
    table()[make_pair(iq1,f3 )].insert(a);
    table()[make_pair(f3 ,iq1)].insert(a);
    // and       ud0 u type
    f3 = CheckId::makeDiquarkID(iq1,flav2,1);
    a.overallWeight /= a.wt;
    a.wt = 0.5;
    a.overallWeight *= a.wt;
    table()[make_pair(iq1,f3 )].insert(a);
    table()[make_pair(f3 ,iq1)].insert(a);
  }
  else if(iq1==flav2) {
    // ud1 u type
    a.wt = 1./6.;
    a.overallWeight *= a.wt;
    table()[make_pair(flav1,flav2)].insert(a);
    table()[make_pair(flav2,flav1)].insert(a);
    // also need ud0 u type
    long f3 = CheckId::makeDiquarkID(iq1,iq2,1);
    a.overallWeight /= a.wt;
    a.wt = 0.5;
    a.overallWeight *= a.wt;
    table()[make_pair(f3    ,flav2)].insert(a);
    table()[make_pair(flav2 ,f3   )].insert(a);
    // and uu1 d type
    f3 = CheckId::makeDiquarkID(iq1,iq1,3);
    a.overallWeight /= a.wt;
    a.wt = 1./3.;
    a.overallWeight *= a.wt;
    table()[make_pair(f3 ,iq2)].insert(a);
    table()[make_pair(iq2, f3)].insert(a);
  }
  else if(iq2==flav2) assert(false);
  else {
    // determine if light quarks in spin 0 or spin 1
    long it1 = (a.id/100)%10;
    long it2 = (a.id/10 )%10;
    // first perm
    double wgt0(1./4.),wgt1(1./12.);
    // only spin-0
    if(it1<it2) {
      long f3 = CheckId::makeDiquarkID(iq1,iq2,1);
      a.wt = 1./3.;
      a.overallWeight *= a.wt;
      table()[make_pair(f3,flav2)].insert(a);
      table()[make_pair(flav2,f3)].insert(a);
      swap(wgt0,wgt1);
    }
    // only spin-1
    else {
      a.wt = 1./3.;
      a.overallWeight *= a.wt;
      table()[make_pair(flav1,flav2)].insert(a);
      table()[make_pair(flav2,flav1)].insert(a);
    }
    // not sure here
    // spin 0
    a.overallWeight /= a.wt;
    a.wt = wgt0;
    a.overallWeight *= a.wt;
    // second perm
    long f3 = CheckId::makeDiquarkID(iq1,flav2,1);
    table()[make_pair(iq2,f3)].insert(a);
    table()[make_pair(f3,iq2)].insert(a);
    // 3rd perm
    f3 = CheckId::makeDiquarkID(iq2,flav2,1);
    table()[make_pair(iq1,f3)].insert(a);
    table()[make_pair(f3,iq1)].insert(a);
    // spin 1
    a.overallWeight /= a.wt;
    a.wt = wgt1;
    a.overallWeight *= a.wt;
    // second perm
    f3 = CheckId::makeDiquarkID(iq1,flav2,3);
    table()[make_pair(iq2,f3)].insert(a);
    table()[make_pair(f3,iq2)].insert(a);
    // 3rd perm
    f3 = CheckId::makeDiquarkID(iq2,flav2,3);
    table()[make_pair(iq1,f3)].insert(a);
    table()[make_pair(f3,iq1)].insert(a);
  }
}

void Hw7Selector::insertThreeHalf(HadronInfo a, int flav1, int flav2) {
  long iq1 = flav1/1000;
  long iq2 = (flav1/100)%10;
  // all the same
  if(iq1==iq2 && iq1==flav2) {
    a.wt = 1.;
    a.overallWeight *= a.wt;
    table()[make_pair(flav1,flav2)].insert(a);
    table()[make_pair(flav2,flav1)].insert(a);
  }
  else if(iq1==iq2) {
    // first option uu1 d
    a.wt = 1./3.;
    a.overallWeight *= a.wt;
    table()[make_pair(flav1,flav2)].insert(a);
    table()[make_pair(flav2,flav1)].insert(a);
    // also need ud1 u type
    long f3 = CheckId::makeDiquarkID(iq1,flav2,3);
    a.overallWeight /= a.wt;
    a.wt = 2./3.;
    a.overallWeight *= a.wt;
    table()[make_pair(iq1,f3 )].insert(a);
    table()[make_pair(f3 ,iq1)].insert(a);
  }
  else if(iq1==flav2) {
    // also need ud1 u type
    a.wt = 2./3.;
    a.overallWeight *= a.wt;
    table()[make_pair(flav1,flav2)].insert(a);
    table()[make_pair(flav2,flav1)].insert(a);
    // and uu1 d type
    long f3 = CheckId::makeDiquarkID(iq1,iq1,3);
    a.overallWeight /= a.wt;
    a.wt = 1./3.;
    a.overallWeight *= a.wt;
    table()[make_pair(iq2,f3)].insert(a);
    table()[make_pair(f3,iq2)].insert(a);
  }
  else if(iq2==flav2) assert(false);
  else {
    // just need the three different combinations
    // first perm
    a.wt = 1./3.;
    a.overallWeight *= a.wt;
    table()[make_pair(flav1,flav2)].insert(a);
    table()[make_pair(flav2,flav1)].insert(a);
    // 2nd perm
    long f3 = CheckId::makeDiquarkID(iq1,flav2,3);
    table()[make_pair(iq2,f3)].insert(a);
    table()[make_pair(f3,iq2)].insert(a);
    // 3rd perm
    f3 = CheckId::makeDiquarkID(iq2,flav2,3);
    table()[make_pair(iq1,f3)].insert(a);
    table()[make_pair(f3,iq1)].insert(a);
  }
}

PDPtr Hw7Selector::makeDiquark(tcPDPtr par1, tcPDPtr par2) {
  long id1 = par1->id(), id2 = par2->id();
  long pspin = 3;
  if(id1!=id2) {
    if(UseRandom::rnd()<_pwtDIquarkS0/(_pwtDIquarkS0+_pwtDIquarkS1)) pspin = 1;
  }
  long idnew = CheckId::makeDiquarkID(id1,id2, pspin);
  return getParticleData(idnew);
}

tcPDPair Hw7Selector::lightestBaryonPair(tcPDPtr ptr1, tcPDPtr ptr2, int pspin) const {
  // Make sure that we don't have any diquarks as input, return arbitrarily
  // large value if we do
  Energy currentSum = Constants::MaxEnergy;
  tcPDPair output;
  for(unsigned int ix=0; ix<partons().size(); ++ix) {
    if(!DiquarkMatcher::Check(partons()[ix]->id())) continue;
    if(partons()[ix]->iSpin()!=pspin) continue;
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
