// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HadronSpectrum class.
//

#include "HadronSpectrum.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <ThePEG/Repository/CurrentGenerator.h>

#include "ThePEG/Interface/RefVector.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HadronSpectrum::HadronSpectrum() 
  : Interfaced(),
    _repwt(Lmax,vector<vector<double> >(Jmax,vector<double>(Nmax))),
    _sngWt( 1.0 ),_decWt( 1.0 ), belowThreshold_(0) {}

HadronSpectrum::~HadronSpectrum() {}

void HadronSpectrum::doinit() {
  Interfaced::doinit();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void HadronSpectrum::persistentOutput(PersistentOStream & os) const {
  os << _table << _partons << _forbidden;
}

void HadronSpectrum::persistentInput(PersistentIStream & is, int) {
  is >> _table >> _partons >> _forbidden;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<HadronSpectrum,Interfaced>
  describeHerwigHadronSpectrum("Herwig::HadronSpectrum", "Herwig.so");

void HadronSpectrum::Init() {

  static ClassDocumentation<HadronSpectrum> documentation
    ("There is no documentation for the HadronSpectrum class");

  static RefVector<HadronSpectrum,ParticleData> interfacePartons
    ("Partons",
     "The partons which are to be considered as the consistuents of the hadrons.",
     &HadronSpectrum::_partons, -1, false, false, true, false, false);

  static RefVector<HadronSpectrum,ParticleData> interfaceForbidden
    ("Forbidden",
     "The PDG codes of the particles which cannot be produced in the hadronization.",
     &HadronSpectrum::_forbidden, -1, false, false, true, false, false);

}


tcPDPtr HadronSpectrum::chooseSingleHadron(tcPDPtr par1, tcPDPtr par2,
							Energy mass) const {
  Energy threshold = hadronPairThreshold(par1,par2);
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

pair<tcPDPtr,tcPDPtr> HadronSpectrum::lightestHadronPair(tcPDPtr ptr1, tcPDPtr ptr2,
								      tcPDPtr ptr3) const {
  // throw exception if id3!=0 as doesn't work
  if ( ptr3 ) throw Exception() 
    << "ptr3!=0 not yet implemented in HadronSpectrum::lightestHadronPair"
    << Exception::abortnow;

  // charge
  int totalcharge = ptr1->iCharge() + ptr2->iCharge();
  if ( ptr3 ) totalcharge += ptr3->iCharge();

  tcPDPtr vIdHad1 = tcPDPtr(), vIdHad2 = tcPDPtr();
  Energy MinMass = ZERO;
  for (long pid : lightestQuarks()) {
    tcPDPtr idPartner = getParticleData(pid);
    // Change sign to idPartner (transform it into a anti-quark) if it is not
    // possible to form a meson or a baryon.
    assert (ptr1 && idPartner);
    if (!canBeHadron(ptr1, idPartner)) idPartner = idPartner->CC();

    tcPDPtr Had1 = lightestHadron(ptr1, idPartner);
    tcPDPtr Had2 = lightestHadron(ptr2, idPartner->CC());
    if (Had1 && Had2 && Had1->iCharge() + Had2->iCharge() == totalcharge) {
      Energy mass = Had1->mass() + Had2->mass();
      if (MinMass == ZERO || mass < MinMass) {
          MinMass = mass;
          vIdHad1 = Had1;
          vIdHad2 = Had2;
      }
    }  
  }
  // Take the lightest pair compatible with charge conservation.
  return make_pair(vIdHad1, vIdHad2);
}

tcPDPtr HadronSpectrum::lightestHadron(tcPDPtr ptr1, tcPDPtr ptr2,
#ifndef NDEBUG
				      tcPDPtr ptr3) const {
#else
				      tcPDPtr ) const {
#endif
  // The method assumes ptr3 == 0 rest not implemented
  assert(ptr1 && ptr2 && !ptr3);
  // find entry in the table
  pair<long,long> ids = make_pair(abs(ptr1->id()),abs(ptr2->id()));
  HadronTable::const_iterator tit=_table.find(ids);
  // throw exception if flavours wrong
  if (tit==_table.end()) 
    throw Exception() << "Could not find " 
		      << ids.first << ' ' << ids.second 
		      << " in _table. "
		      << "In HadronSpectrum::lightestHadron()"
		      << Exception::eventerror;
  if(tit->second.empty())
    throw Exception() << "HadronSpectrum::lightestHadron "
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
    generator()->log() << "HadronSpectrum::lightestHadron: "
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
HadronSpectrum::hadronsBelowThreshold(Energy threshold, tcPDPtr ptr1,
				      tcPDPtr ptr2,
#ifndef NDEBUG
				      tcPDPtr ptr3) const {
#else
				      tcPDPtr ) const {
#endif
  // The method assumes ptr3 == 0 rest not implemented
  assert(ptr1 && ptr2 && !ptr3);
  // find entry in the table
  pair<long,long> ids = make_pair(abs(ptr1->id()),abs(ptr2->id()));
  HadronTable::const_iterator tit=_table.find(ids);
  // throw exception if flavours wrong
  if (tit==_table.end()) 
    throw Exception() << "Could not find " 
		      << ids.first << ' ' << ids.second 
		      << " in _table. "
		      << "In HadronSpectrum::hadronsBelowThreshold()"
		      << Exception::eventerror;
  if(tit->second.empty())
    throw Exception() << "HadronSpectrum::hadronsBelowThreshold() "
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
      generator()->log() << "HadronSpectrum::hadronsBelowTheshold: "
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

Energy HadronSpectrum::massLightestBaryonPair(tcPDPtr ptr1, tcPDPtr ptr2) const {
  // Make sure that we don't have any diquarks as input, return arbitrarily
  // large value if we do
  Energy currentSum = Constants::MaxEnergy; 
  for(unsigned int ix=0; ix<_partons.size(); ++ix) {
    if(!DiquarkMatcher::Check(_partons[ix]->id())) continue;
    HadronTable::const_iterator 
      tit1=_table.find(make_pair(abs(ptr1->id()),_partons[ix]->id())),
      tit2=_table.find(make_pair(_partons[ix]->id(),abs(ptr2->id())));
    if( tit1==_table.end() || tit2==_table.end()) continue;
    if(tit1->second.empty()||tit2->second.empty()) continue;
    Energy s = tit1->second.begin()->mass + tit2->second.begin()->mass;
    if(currentSum > s) currentSum = s;
  }
  return currentSum;
}


double HadronSpectrum::specialWeight(long id) const {
  const int pspin = id % 10;  
  // Only K0L and K0S have pspin == 0, should
  // not get them until Decay step
  assert( pspin != 0 );
  // Baryon : J = 1/2 or 3/2
  if(pspin == 2) {
    // Singlet (Lambda-like) baryon
    if( (id/100)%10 < (id/10 )%10 ) return sqr(_sngWt);   
    // octet
    else                            return 1.;
  } 
  // Decuplet baryon
  else if (pspin == 4) {
    return sqr(_decWt);
  }    
  // Meson
  else if(pspin % 2 == 1) {
    // Total angular momentum
    int j  = (pspin - 1) / 2;
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
  }
  // rest is not excited or 
  // has spin >= 5/2 (ispin >= 6), haven't got those
  return 1.0;
}

int HadronSpectrum::signHadron(tcPDPtr idQ1, tcPDPtr idQ2, 
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

PDPtr HadronSpectrum::makeDiquark(tcPDPtr par1, tcPDPtr par2) const {
    long id1 = par1->id();
    long id2 = par2->id();
    long idnew = makeDiquarkID(id1,id2);
    assert(!CurrentGenerator::isVoid());
    return CurrentGenerator::current().getParticleData(idnew);
}

bool HadronSpectrum::canBeMeson(tcPDPtr par1,tcPDPtr par2) const {
  assert(par1 && par2);
  long id1 = par1->id();
  long id2 = par2->id();
  // a Meson must not have any diquarks
  if(DiquarkMatcher::Check(id1) || DiquarkMatcher::Check(id2)) return false;
  return ( abs(int(par1->iColour()))== 3  && 
     abs(int(par2->iColour())) == 3 &&  
     id1*id2 < 0);
}

bool HadronSpectrum::canBeBaryon(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3) const {
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
  
