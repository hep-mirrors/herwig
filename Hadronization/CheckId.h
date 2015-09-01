// -*- C++ -*-
//
// CheckId.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef HERWIG_CheckId_H
#define HERWIG_CheckId_H
//
// This is the declaration of the CheckId class.

#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/PDT/ParticleData.h"
#include <ThePEG/PDT/EnumParticles.h>
#include "ThePEG/Repository/CurrentGenerator.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Hadronization
 *
 *  This is a namespace which provides some useful methods
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
 *  NB) For Other useful methods (even some implemented in CheckId)
 *      @see StandardMatchers
 */
namespace CheckId {
  
  /**
   * Return the id of the diquark (anti-diquark) made by the two 
   * quarks (antiquarks) of id specified in input (id1, id2).
   * Caller must ensure that id1 and id2 are quarks.
   */
  long makeDiquarkID(long id1, long id2);
  
  /**
   * Return the particle data of the diquark (anti-diquark) made by the two 
   * quarks (antiquarks) par1, par2.
   * @param par1 (anti-)quark data pointer
   * @param par2 (anti-)quark data pointer
   */
  PDPtr makeDiquark(tcPDPtr par1, tcPDPtr par2);

  /**
   * Return true if the two particles in input can be the components of a meson;
   *false otherwise.
   */
  bool canBeMeson(tcPDPtr par1,tcPDPtr par2);

  /**
   * Return true if the two or three particles in input can be the components 
   * of a baryon; false otherwise.
   */
  bool canBeBaryon(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3 = PDPtr());
  
   /**
   * Return true if the two or three particles in input can be the components 
   * of a hadron; false otherwise.
   */
  inline bool canBeHadron(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3 = PDPtr())  {
    return (!par3 && canBeMeson(par1,par2)) || canBeBaryon(par1,par2,par3);
  }

 
  /**
   * Return true if any of the possible three input particles has
   * b-flavour; 
   * false otherwise. In the case that only the first particle is specified,
   * it can be: an (anti-)quark, an (anti-)diquark
   * an (anti-)meson, an (anti-)baryon; in the other cases, each pointer
   * is assumed to be either (anti-)quark or (anti-)diquark.
   */
  bool hasBottom(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr());
  /**
   * Return true if any of the possible three input particles has 
   * c-flavour; 
   * false otherwise.In the case that only the first pointer is specified,
   * it can be: a (anti-)quark, a (anti-)diquark
   * a (anti-)meson, a (anti-)baryon; in the other cases, each pointer
   * is assumed to be either (anti-)quark or (anti-)diquark.
   */
  bool hasCharm(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr());
  /**
   * Return true, if any of the possible input particle pointer is an exotic quark, e.g. Susy quark;
   * false otherwise.   
   */
  bool isExotic(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr());

}
  
}

#endif /* HERWIG_CheckId_H */
