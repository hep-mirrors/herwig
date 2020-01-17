// -*- C++ -*-
//
// Hw64Selector.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Hw64Selector class.
//

#include "Hw64Selector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/UseRandom.h"
#include "CheckId.h"
#include "Herwig/Utilities/Kinematics.h"
#include <ThePEG/Utilities/DescribeClass.h>

using namespace Herwig;

DescribeNoPIOClass<Hw64Selector,HadronSelector>
describeHw64Selector("Herwig::Hw64Selector","");

IBPtr Hw64Selector::clone() const {
  return new_ptr(*this);
}

IBPtr Hw64Selector::fullclone() const {
  return new_ptr(*this);
}

void Hw64Selector::Init() {

  static ClassDocumentation<Hw64Selector> documentation
    ("There is no documentation for the Hw64Selector class");

}

pair<tcPDPtr,tcPDPtr> Hw64Selector::chooseHadronPair(const Energy cluMass,tcPDPtr par1, 
						     tcPDPtr par2,tcPDPtr) const
  {
  bool diquark = !(DiquarkMatcher::Check(par1->id()) || DiquarkMatcher::Check(par2->id()));



  pair<tcPDPtr,tcPDPtr> lighthad = lightestHadronPair(par1, par2);
  if(!lighthad.first || !lighthad.second)
    throw Exception() << "Hw64Selector::chooseHadronPair "
		      << "We have 0's! First id = " << par1->id() << " second = " 
		      << par2->id() << ". This is probably a problem with either"
		      << " undecayed heavy particles or colour connections" 
		      << Exception::eventerror;
  // calculate maximum momentum
  Energy PCMax = Kinematics::pstarTwoBodyDecay(cluMass,lighthad.first->mass(),
					       lighthad.second->mass());
  tcPDPtr had1 = tcPDPtr();
  tcPDPtr had2 = tcPDPtr();
  int ntry = 0;
  tcPDPtr quark = tcPDPtr();
  const int nmax = 5000;
  Energy p;
  do {
    quark = partons()[UseRandom::irnd(partons().size())];
    if(diquark && DiquarkMatcher::Check(quark->id())) continue;
    if(pwt(quark->id()) <= UseRandom::rnd()) continue;
    pair<long,long> pid(abs(par1->id()),quark->id());
    KupcoData::const_iterator it1,it2;
    const HadronTable::const_iterator tit = table().find(pid);
    assert(tit != table().end());
    const KupcoData & hdata = tit->second;    
    do {
      it1 = hdata.begin();
      advance(it1,int(hdata.size()*UseRandom::rnd()));
    } 
    while(it1 != hdata.end() && it1->overallWeight < UseRandom::rnd());
    
    had1 = it1->ptrData;
    pid = make_pair(quark->id(),abs(par2->id())); 
    do {
      it2 = hdata.begin();
      advance(it2,int(hdata.size()*UseRandom::rnd()));
    }
    while(it2 != hdata.end() && it2->overallWeight < UseRandom::rnd());
    had2 = it2->ptrData;
    if(had1 && had2) {
      p = Kinematics::pstarTwoBodyDecay(cluMass, it1->mass, it2->mass);
      if(p/PCMax < UseRandom::rnd()) { 
	had1 = had2 = tcPDPtr(); 
	ntry++;
      }
    }
  }
  while((!had1|| !had2) && ntry < nmax);
  if(ntry >= nmax) return lighthad;
  int signHad1 = 0;
  int signHad2 = 0;
  if(CheckId::canBeHadron(par1,quark->CC()) && CheckId::canBeHadron(quark,par2)) {
    signHad1 = signHadron(par1, quark->CC(), had1);
    signHad2 = signHadron(par2, quark, had2);
  }
  else if(CheckId::canBeHadron(par1,quark) && CheckId::canBeHadron(quark->CC(),par2)) {
    signHad1 = signHadron(par1, quark, had1);
    signHad2 = signHadron(par2, quark->CC(), had2);
  }
	 
 else throw Exception() << "Hw64Selector::chooseHadronPair()"
			 << Exception::abortnow;
  return make_pair( signHad1 > 0 ? had1 : tcPDPtr(had1->CC()),
		    signHad2 > 0 ? had2 : tcPDPtr(had2->CC()));
}
