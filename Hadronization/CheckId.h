// -*- C++ -*-
//
#ifndef HERWIG_CheckId_H
#define HERWIG_CheckId_H
//
// This is the declaration of the CheckId class.

#include "Herwig++/Config/Herwig.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/PDT/ParticleData.h"
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
   * Return the particle data of the diquark (anti-diquark) made by the two 
   * quarks (antiquarks) par1, par2.
   * @param par1 (anti-)quark data pointer
   * @param par2 (anti-)quark data pointer
   */
  static PDPtr makeDiquark(tcPDPtr par1, tcPDPtr par2);


  /**
   * Return the id of the diquark (anti-diquark) made by the two 
   * quarks (antiquarks) of id specified in input (id1, id2).
   * Caller must ensure that id1 and id2 are quarks.
   */
  static long makeDiquarkID(long id1, long id2);
  
  /**
   * Return true if the two particles in input can be the components of a meson;
   *false otherwise.
   */
  static bool canBeMeson(tcPDPtr par1,tcPDPtr par2);
  
  /**
   * Return true if the two or three particles in input can be the components 
   * of a baryon; false otherwise.
   */
  static bool canBeBaryon(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3 = PDPtr());

   /**
   * Return true if the two or three particles in input can be the components 
   * of a hadron; false otherwise.
   */
  static bool canBeHadron(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3 = PDPtr());
 
  /**
   * Return true if any of the possible three input particles has
   * b-flavour; 
   * false otherwise. In the case that only the first particle is specified,
   * it can be: an (anti-)quark, an (anti-)diquark
   * an (anti-)meson, an (anti-)baryon; in the other cases, each pointer
   * is assumed to be either (anti-)quark or (anti-)diquark.
   */
  static bool hasBottom(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr());
  /**
   * Return true if any of the possible three input particles has 
   * c-flavour; 
   * false otherwise.In the case that only the first pointer is specified,
   * it can be: a (anti-)quark, a (anti-)diquark
   * a (anti-)meson, a (anti-)baryon; in the other cases, each pointer
   * is assumed to be either (anti-)quark or (anti-)diquark.
   */
  static bool hasCharm(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr());
  /**
   * Return true, if any of the possible input particle pointer is an exotic quark, e.g. Susy quark;
   * false otherwise.   
   */
  static bool isExotic(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr());
private:
  
  /**
   * Return true if the particle pointer corresponds to a diquark 
   * or anti-diquark carrying b flavour; false otherwise.
   */
  static bool isDiquarkWithB(tcPDPtr par1);
  
  /**
   * Return true if the particle pointer corresponds to a diquark
   *  or anti-diquark carrying c flavour; false otherwise.
   */
  static bool isDiquarkWithC(tcPDPtr par1);

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

#endif /* HERWIG_CheckId_H */

