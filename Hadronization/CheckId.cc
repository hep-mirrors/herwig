// -*- C++ -*-
//
// CheckId.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CheckId class.
//
#include "CheckId.h"
#include <cassert>

using namespace Herwig;

namespace {
  
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

long CheckId::makeDiquarkID(long id1, long id2) {

  assert( id1 * id2 > 0  
          && QuarkMatcher::Check(id1)  
	  && QuarkMatcher::Check(id2)) ;
  long ida = abs(id1);
  long idb = abs(id2);
  if (ida < idb) swap(ida,idb);

  long idnew = ida*1000 + idb*100 + 1;
  // Diquarks made of quarks of the same type: uu, dd, ss, cc, bb, 
  // have spin 1, and therefore the less significant digit (which
  // corresponds to 2*J+1) is 3 rather than 1 as all other Diquarks.
  if (id1 == id2) idnew += 2;

  return id1 > 0 ? idnew : -idnew;
}

PDPtr CheckId::makeDiquark(tcPDPtr par1, tcPDPtr par2) {
    long id1 = par1->id();
    long id2 = par2->id();
    long idnew = makeDiquarkID(id1,id2);
    assert(!CurrentGenerator::isVoid());
    return CurrentGenerator::current().getParticleData(idnew);
}

bool CheckId::canBeMeson(tcPDPtr par1,tcPDPtr par2) {
  assert(par1 && par2);
  long id1 = par1->id();
  long id2 = par2->id();
  // a Meson must not have any diquarks
  if(DiquarkMatcher::Check(id1) || DiquarkMatcher::Check(id2)) return false;
  return ( abs(int(par1->iColour()))== 3  && 
     abs(int(par2->iColour())) == 3 &&  
     id1*id2 < 0);
}



bool CheckId::canBeBaryon(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3) {
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
  


bool CheckId::hasBottom(tcPDPtr par1, tcPDPtr par2, tcPDPtr par3) {
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


bool CheckId::hasCharm(tcPDPtr par1, tcPDPtr par2, tcPDPtr par3) {
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

bool CheckId::isExotic(tcPDPtr par1, tcPDPtr par2, tcPDPtr par3) {
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
