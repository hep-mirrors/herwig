// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HadronSpectrum class.
//

#include "HadronSpectrum.h"
#include "ClusterHadronizationHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <ThePEG/Repository/CurrentGenerator.h>
#include "Herwig/Utilities/Kinematics.h"

#include "ThePEG/Interface/RefVector.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

namespace {
  // debug helper
  void dumpTable(const HadronSpectrum::HadronTable & tbl) {
    typedef HadronSpectrum::HadronTable::const_iterator TableIter;
    for (TableIter it = tbl.begin(); it != tbl.end(); ++it) {
      cerr << it->first.first << ' ' 
  	   << it->first.second << '\n';
      for (HadronSpectrum::KupcoData::const_iterator jt = it->second.begin();
  	   jt != it->second.end(); ++jt) {
  	cerr << '\t' << *jt << '\n';
      }
    }
  }
}


HadronSpectrum::HadronSpectrum() 
  : Interfaced(),
    belowThreshold_(0),
    _repwt(Lmax,vector<vector<double> >(Jmax,vector<double>(Nmax))) {}

HadronSpectrum::~HadronSpectrum() {}

void HadronSpectrum::doinit() {
  Interfaced::doinit();
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
  if (Debug::level >= 10) 
    dumpTable(table());
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void HadronSpectrum::persistentOutput(PersistentOStream & os) const {
  os << _table << _partons << _forbidden
     << belowThreshold_ << _repwt << _pwt << lightestHadrons_;
}

void HadronSpectrum::persistentInput(PersistentIStream & is, int) {
  is >> _table >> _partons >> _forbidden
     >> belowThreshold_ >> _repwt >> _pwt >> lightestHadrons_;
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

void HadronSpectrum::insertToHadronTable(tPDPtr &particle, int flav1, int flav2) {
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

void HadronSpectrum::insertOneHalf(HadronInfo a, int flav1, int flav2) {
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
      _table[make_pair(flav1,flav2)].insert(a);
      _table[make_pair(flav2,flav1)].insert(a);
      long f3 = makeDiquarkID(iq1,flav2,1);
      _table[make_pair(iq1,f3 )].insert(a);
      _table[make_pair(f3 ,iq1)].insert(a);
    }
  }
  else if(iq1==flav2) {
    // ud1 u type
    _table[make_pair(flav1,flav2)].insert(a);
    _table[make_pair(flav2,flav1)].insert(a);
    // and uu1 d type
    long f3 = makeDiquarkID(iq1,iq1,3);
    a.overallWeight *= a.wt;
    _table[make_pair(f3 ,iq2)].insert(a);
    _table[make_pair(iq2, f3)].insert(a);
  }
  else if(iq2==flav2) assert(false);
  else {
    _table[make_pair(flav1,flav2)].insert(a);
    _table[make_pair(flav2,flav1)].insert(a);
    long f3 = makeDiquarkID(iq1,flav2,1);
    _table[make_pair(iq2,f3)].insert(a);
    _table[make_pair(f3,iq2)].insert(a);
    // 3rd perm
    f3 = makeDiquarkID(iq2,flav2,1);
    _table[make_pair(iq1,f3)].insert(a);
    _table[make_pair(f3,iq1)].insert(a);
  }
}

void HadronSpectrum::insertThreeHalf(HadronInfo a, int flav1, int flav2) {
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
      _table[make_pair(flav1,flav2)].insert(a);
      _table[make_pair(flav2,flav1)].insert(a);
      long f3 = makeDiquarkID(iq1,flav2,1);
      _table[make_pair(iq1,f3 )].insert(a);
      _table[make_pair(f3 ,iq1)].insert(a);
    }
  }
  else if(iq1==flav2) {
    // ud1 u type
    _table[make_pair(flav1,flav2)].insert(a);
    _table[make_pair(flav2,flav1)].insert(a);
    // and uu1 d type
    long f3 = makeDiquarkID(iq1,iq1,3);
    a.overallWeight *= a.wt;
    _table[make_pair(f3 ,iq2)].insert(a);
    _table[make_pair(iq2, f3)].insert(a);
  }
  else {
    _table[make_pair(flav1,flav2)].insert(a);
    _table[make_pair(flav2,flav1)].insert(a);
    long f3 = makeDiquarkID(iq1,flav2,1);
    _table[make_pair(iq2,f3)].insert(a);
    _table[make_pair(f3,iq2)].insert(a);
    // 3rd perm
    f3 = makeDiquarkID(iq2,flav2,1);
    _table[make_pair(iq1,f3)].insert(a);
    _table[make_pair(f3,iq1)].insert(a);
  }
}


tcPDPtr HadronSpectrum::chooseSingleHadron(tcPDPtr par1, tcPDPtr par2,
							Energy mass) const {
  Energy threshold = hadronPairThreshold(par1,par2);
  // only do one hadron decay if mass less than the threshold
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

tcPDPair HadronSpectrum::chooseHadronPair(const Energy cluMass,
					  tcPDPtr par1, tcPDPtr par2) const {
  useMe();
  // if either of the input partons is a diquark don't allow diquarks to be
  // produced
  bool isDiquark1 = DiquarkMatcher::Check(par1->id());
  bool isDiquark2 = DiquarkMatcher::Check(par2->id());
  bool noDiquarkInCluster = !(isDiquark1 || isDiquark2);
  bool oneDiquarkInCluster = (isDiquark1 != isDiquark2);
  bool quark = true;
  // decide is baryon or meson production
  if(noDiquarkInCluster) std::tie(quark,noDiquarkInCluster,oneDiquarkInCluster)
  						= selectBaryon(cluMass,par1,par2);
  // weights for the different possibilities
  Energy weight, wgtsum(ZERO);
  // loop over all hadron pairs with the allowed flavours
  static vector<Kupco> hadrons;
  hadrons.clear();
  for(unsigned int ix=0;ix<partons().size();++ix) {
    tcPDPtr quarktopick  = partons()[ix];
    if(!quark && std::find(hadronizingQuarks().begin(), hadronizingQuarks().end(),
        abs(quarktopick->id())) != hadronizingQuarks().end()) continue;
    if(DiquarkMatcher::Check(quarktopick->id()) &&
       ((!noDiquarkInCluster && quarktopick->iSpin()==1) ||
	(!oneDiquarkInCluster && quarktopick->iSpin()==3))) continue;
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
    quarkWeight = specialQuarkWeight(quarkWeight,quarktopick->id(),
            cluMass,par1,par2);
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

	if(canBeHadron(par1, quarktopick->CC())
	   && canBeHadron(quarktopick, par2))
	   signQ = +1;
	else if(canBeHadron(par1, quarktopick)
		&& canBeHadron(quarktopick->CC(), par2))
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

std::tuple<bool,bool,bool> HadronSpectrum::selectBaryon(const Energy, tcPDPtr, tcPDPtr )  const {
  assert(false);
}

tcPDPair HadronSpectrum::lightestHadronPair(tcPDPtr ptr1, tcPDPtr ptr2) const {
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

tcPDPtr HadronSpectrum::lightestHadron(tcPDPtr ptr1, tcPDPtr ptr2) const {
  assert(ptr1 && ptr2);
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
				      tcPDPtr ptr2) const {
  assert(ptr1 && ptr2);
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

double HadronSpectrum::mesonWeight(long id) const {
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
    if(std::find(hadronizingQuarks().begin(), hadronizingQuarks().end(),
                 abs(idQ1->id())) != hadronizingQuarks().end() &&
       std::find(hadronizingQuarks().begin(), hadronizingQuarks().end(),
                 abs(idQ2->id())) != hadronizingQuarks().end())
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
    long pspin = id1==id2 ? 3 : 1;
    long idnew = makeDiquarkID(id1,id2, pspin);
    assert(!CurrentGenerator::isVoid());
    return CurrentGenerator::current().getParticleData(idnew);
}

bool HadronSpectrum::canBeMeson(tcPDPtr par1,tcPDPtr par2) const {
  assert(par1 && par2);
  long id1 = par1->id();
  long id2 = par2->id();
  // a Meson must not have any diquarks
  if(DiquarkMatcher::Check(id1) || DiquarkMatcher::Check(id2)) return false;
  return (std::find(hadronizingQuarks().begin(), hadronizingQuarks().end(),
                    abs(id1)) != hadronizingQuarks().end() &&
          std::find(hadronizingQuarks().begin(), hadronizingQuarks().end(),
                    abs(id2)) != hadronizingQuarks().end() &&
          id1*id2 < 0);
}
  
