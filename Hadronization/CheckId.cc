// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CheckId class.
//
#include "CheckId.h"

using namespace Herwig;

long CheckId::diquarkId(const long id1, const long id2) {
  if ( id1 * id2 < 0  ||  
       ! QuarkMatcher::Check(id1)  ||  ! QuarkMatcher::Check(id2) ) return 0;
  long ida = abs(id1);
  long idb = abs(id2);
  if (ida < idb) {
    ida = abs(id2);
    idb = abs(id1);
  }
  long idnew = ida*1000 + idb*100 + 1;
  // Diquarks made of quarks of the same type: uu, dd, ss, cc, bb, 
  // have spin 1, and therefore the less significant digit (which
  // corresponds to 2*J+1) is 3 rather than 1 as all other Diquarks.
  if (id1 == id2) idnew += 2;
  return id1 > 0 ? idnew : -idnew ;
}

bool CheckId::canBeBaryon(const long id1, const long id2, const long id3) {
  bool result = false;
  if (id3 == 0) {
    // In this case, to be a baryon, one component must be a (anti-)quark and
    // the other a (anti-)diquark; furthermore, only:
    //   quark - diquark   or  anti-quark - anti-diquark  
    // (that means same sign for both components) are allowed.
    result = id1*id2 > 0  &&
      ( ( QuarkMatcher::Check(id1)   && DiquarkMatcher::Check(id2) ) || 
	( DiquarkMatcher::Check(id1) &&   QuarkMatcher::Check(id2) ) );
  } 
  else {
    // In this case, to be a baryon, all three components must be (anti-)quarks
    // and with the same sign.
    result = id1*id2 > 0  &&  id2*id3 > 0  && 
      QuarkMatcher::Check(id1) && QuarkMatcher::Check(id2) && QuarkMatcher::Check(id3);
  }
  return result;    
}

bool CheckId::hasBeauty(const long id1, const long id2, const long id3) {
  if ( id2 == 0  &&  id3 == 0 ) {
    return 
      abs(id1) == ThePEG::ParticleID::b    ||
      isDiquarkWithB(id1)                  ||
      MesonMatcher::Check(id1)  && (abs(id1)/100)%10  == ThePEG::ParticleID::b ||
      BaryonMatcher::Check(id1) && (abs(id1)/1000)%10 == ThePEG::ParticleID::b;
  } 
  else {
    return 
      abs(id1) == ThePEG::ParticleID::b  ||  isDiquarkWithB(id1)  || 
      abs(id2) == ThePEG::ParticleID::b  ||  isDiquarkWithB(id2)  || 
      abs(id3) == ThePEG::ParticleID::b  ||  isDiquarkWithB(id3); 
  }
}  

bool CheckId::hasCharm(const long id1, const long id2, const long id3) {
  if ( id2 == 0  &&  id3 == 0 ) {
    return
      abs(id1) == ThePEG::ParticleID::c     ||
      isDiquarkWithC(id1)                   ||
      MesonMatcher::Check(id1) && ( (abs(id1)/100)%10   == ThePEG::ParticleID::c   ||
				    (abs(id1)/10)%10    == ThePEG::ParticleID::c ) ||
      BaryonMatcher::Check(id1) && ( (abs(id1)/1000)%10 == ThePEG::ParticleID::c  ||
				     (abs(id1)/100)%10  == ThePEG::ParticleID::c  ||
				     (abs(id1)/10)%10   == ThePEG::ParticleID::c );
  } 
  else {
    return 
      abs(id1) == ThePEG::ParticleID::c  ||  isDiquarkWithC(id1)  || 
      abs(id2) == ThePEG::ParticleID::c  ||  isDiquarkWithC(id2)  || 
      abs(id3) == ThePEG::ParticleID::c  ||  isDiquarkWithC(id3); 
  }
}  
