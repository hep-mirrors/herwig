// -*- C++ -*-
//
#ifndef HERWIG_CheckId_H
#define HERWIG_CheckId_H
//
// This is the declaration of the CheckId class.

#include "Herwig++/Config/Herwig.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include <ThePEG/PDT/EnumParticles.h>


namespace Herwig {

using namespace ThePEG;

/** \ingroup Hadronization
 *
 *  This is a pure static class which provides some useful methods
 *  for checking the PDG id of particles.
 *  Notice that the name "quark" in the methods below means any of 
 *  the following:
 *  d , u , s , c , b   anti-d , anti-u , anti-s , anti-c , anti-b
 *  that is we do not distinguish between particle or antiparticle 
 *  ( because can be done directly by the user: 
 *        id > 0  for  particles ;  id < 0  for  anti-particles )
 *  and we do not include  t  (and  anti-t ) because we are interested
 *  in components of mesons and baryons.
 *  Similarly for the name  "diquark"  which include all  diquarks (id > 0)
 *  and  anti-diquarks  (id < 0)  not made with  t  ( anti-t ) component. 
 *
 *  NB) For Other useful methods (even some implemented in CheckId class!)
 *      @see StandardMatchers
 */
class CheckId {
  
public:
  
  /**
   * Return the id of the diquark (anti-diquark) made by the two 
   * quarks (antiquarks) of id specified in input (id1, id2).
   * Return 0 if something goes wrong. 
   */
  static long diquarkId(const long id1, const long id2);
  
  /**
   * Return true if the two ids in input can be the components of a meson;
   *false otherwise.
   */
  static inline bool canBeMeson(const long id1, const long id2);
  
  /**
   * Return true if the two or three ids in input can be the components 
   * of a baryon; false otherwise.
   */
  static bool canBeBaryon(const long id1, const long id2, const long id3=0);
  
  /**
   * Return true if any of the possible three input ids has b-flavour; 
   * false otherwise. In the case that only the first id is specified,
   * then the corresponding particle can be: a (anti-)quark, a (anti-)diquark
   * a (anti-)meson, a (anti-)baryon; in the other cases, each id 
   * is assumed to be either (anti-)quark or (anti-)diquark.
   */
  static bool hasBeauty(const long id1, const long id2=0, const long id3=0);
  
  /**
   * Return true if any of the possible three input ids has c-flavour; 
   * false otherwise. In the case that only the first id is specified,
   * then the corresponding particle can be: a (anti-)quark, a (anti-)diquark
   * a (anti-)meson, a (anti-)baryon; in the other cases, each id 
   * is assumed to be either (anti-)quark or (anti-)diquark.
   */
  static bool hasCharm(const long id1, const long id2=0, const long id3=0);
  
private:
  
  /**
   * Return true if the id corresponds to a diquark or anti-diquark 
   * carrying b flavour; false otherwise.
   */
  static inline bool isDiquarkWithB(const long id);
  
  /**
   * Return true if the id corresponds to a diquark or anti-diquark 
   * carrying c flavour; false otherwise.
   */
  static inline bool isDiquarkWithC(const long id);

private:
  
  /**
   * Pure static class so default constructor is private
   */
  CheckId();
  
  /**
   * Pure static class so copy constructor is private
   */
  CheckId(const CheckId & x);
  
  /**
   *  Assignmet is private as static
   */
  CheckId & operator=(const CheckId & x);
  
};
  
}

#include "CheckId.icc"

#endif /* HERWIG_CheckId_H */

