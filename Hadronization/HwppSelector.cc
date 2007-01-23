// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwppSelector class.
//

#include "HwppSelector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/Utilities/Selector.h"
#include "ThePEG/Repository/UseRandom.h"
#include "CheckId.h"
#include <cassert>

using namespace Herwig;

void HwppSelector::doinit() throw(InitException) {
  HadronSelector::doinit();
  for(unsigned int ix=0;ix<partons().size();++ix) {
    for(unsigned int iy=0;iy<partons().size();++iy) {
      _baryonmass[make_pair(partons()[ix],partons()[iy])]=
	massLightestBaryonPair(partons()[ix],partons()[iy]);
    }
  }
}

void HwppSelector::persistentOutput(PersistentOStream & os) const {
  os << _mode << _baryonmass;
}

void HwppSelector::persistentInput(PersistentIStream & is, int) {
  is >> _mode >> _baryonmass;
}

ClassDescription<HwppSelector> HwppSelector::initHwppSelector;
// Definition of the static class description member.

void HwppSelector::Init() {

  static ClassDocumentation<HwppSelector> documentation
    ("The HwppSelector class implements the Herwig++ algorithm for selecting"
     " the hadrons");


  static Switch<HwppSelector,unsigned int> interfaceMode
    ("Mode",
     "Which algorithm to use",
     &HwppSelector::_mode, 1, false, false);
  static SwitchOption interfaceModeKupco
    (interfaceMode,
     "Kupco",
     "Use the Kupco approach",
     0);
  static SwitchOption interfaceModeHwpp
    (interfaceMode,
     "Hwpp",
     "Use the Herwig++ approach",
     1);

}

pair<tcPDPtr,tcPDPtr> HwppSelector::chooseHadronPair(const Energy cluMass,
						     const long id1,
						     const long id2, const long)
  throw(Veto, Stop, Exception) {
  // if either of the input partons is a diquark don't allow diquarks to be 
  // produced
  bool diquark = !(DiquarkMatcher::Check(id1) || DiquarkMatcher::Check(id2));
  bool quark=true;
  // if the Herwig++ algorithm 
  if(_mode==1) {
    if(cluMass > _baryonmass[make_pair(abs(id1),abs(id2))] && 
       UseRandom::rnd() > 1./(1.+pwtDIquark())) {
      diquark=true;
      quark=false;
    }
    else {
      diquark=false;
      quark=true;
    }
  }
  // weights for the different possibilities
  vector<Kupco> hadrons;
  Energy weight,wgtsum(0.);
  // loop over all hadron pairs with the allowed flavours
  for(unsigned int ix=0;ix<partons().size();++ix) {
    if(!quark  &&  QuarkMatcher::Check(partons()[ix])) continue;
    if(!diquark&&DiquarkMatcher::Check(partons()[ix])) continue;
    map<pair<long,long>,KupcoData>::const_iterator 
      tit1=table().find(make_pair(abs(id1),partons()[ix]));
    map<pair<long,long>,KupcoData>::const_iterator 
      tit2=table().find(make_pair(partons()[ix],abs(id2)));
    // if not in table skip
    if(tit1==table().end()||tit2==table().end()) continue;
    // tables empty skip
    if(tit1->second.empty()||tit2->second.empty()) continue;
    // if too massive skip
    if(cluMass <= tit1->second.begin()->mass + 
                  tit2->second.begin()->mass) continue; 
    // loop over the hadrons
    KupcoData::iterator H1,H2;
    for(H1 = tit1->second.begin();H1 != tit1->second.end(); ++H1) {
      for(H2 = tit2->second.begin();H2 != tit2->second.end(); ++H2) {
 	// break if cluster too light
 	if(cluMass < H1->mass + H2->mass) break;
 	// calculate the weight
 	weight = pwt()[partons()[ix]]*H1->overallWeight*H2->overallWeight*
 	  Kinematics::pstarTwoBodyDecay(cluMass, H1->mass, H2->mass );
	int signQ = 0;
	long idQ = partons()[ix];
	if((CheckId::canBeMeson(id1,-idQ) || CheckId::canBeBaryon(id1,-idQ)) && 
	   (CheckId::canBeMeson(idQ, id2) || CheckId::canBeBaryon(idQ, id2)))
	  signQ = +1;
	else if((CheckId::canBeMeson( id1,idQ) || CheckId::canBeBaryon( id1,idQ)) && 
		(CheckId::canBeMeson(-idQ,id2) || CheckId::canBeBaryon(-idQ,id2)))
	  signQ = -1;
	// construct the object with the info
	Kupco a(signQ * idQ,H1->ptrData,H2->ptrData,weight);
	hadrons.push_back(a);
	wgtsum+=weight;
      }
    }
  }
  if (hadrons.empty()) 
    return make_pair(tcPDPtr(),tcPDPtr());
  // select the hadron
  wgtsum *=UseRandom::rnd();
  unsigned int ix=0;
  do {
    wgtsum-=hadrons[ix].weight;
    ++ix;
  }
  while(wgtsum>0&&ix<hadrons.size());
  if(ix==hadrons.size()&&wgtsum>0) return make_pair(tcPDPtr(),tcPDPtr());
  --ix;
  int signHad1 = signHadron(id1,-hadrons[ix].idQ, hadrons[ix].hadron1);
  int signHad2 = signHadron(id2, hadrons[ix].idQ, hadrons[ix].hadron2);
  if(signHad1==0||signHad2==0)
    throw Exception() << "HwppSelector::selectPair " 
		      << "***Inconsistent Hadron " 
		      << hadrons[ix].idQ << " " 
		      << hadrons[ix].hadron1->id() << " " 
		      << hadrons[ix].hadron2->id() << " " 
		      << signHad1 << " " << signHad2
		      << Exception::runerror;
  return make_pair
    ( signHad1 > 0 ? hadrons[ix].hadron1 : tcPDPtr(hadrons[ix].hadron1->CC()),
      signHad2 > 0 ? hadrons[ix].hadron2 : tcPDPtr(hadrons[ix].hadron2->CC()));
}
