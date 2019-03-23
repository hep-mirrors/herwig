// -*- C++ -*-
//
// QTildeFinder.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_QTildeFinder_H
#define HERWIG_QTildeFinder_H
//
// This is the declaration of the QTildeFinder class.
//

#include "Herwig/Shower/QTilde/Base/PartnerFinder.h"
#include "Herwig/Shower/Core/ShowerConfig.h"
#include "ThePEG/Interface/Interfaced.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Shower
 *
 *  The QTildeFinder class is responsible for finding the partners and
 *  setting the initial evolution scales for the shower evolution described
 *  in JHEP 0312:045,2003.
 *
 * @see \ref QTildeFinderInterfaces "The interfaces"
 * defined for QTildeFinder.
 */
class QTildeFinder: public PartnerFinder {

public:

  /**
   * The default constructor.
   */
  QTildeFinder() :  _finalFinalConditions(0),
		    _initialFinalDecayConditions(0),
		    _initialInitialConditions(0) {}

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

public:

  /**
   *  Calculate the initial evolution scales given momenta
   */
  pair<Energy,Energy> calculateFinalFinalScales(Lorentz5Momentum p1, Lorentz5Momentum p2,
						bool colouredfirst);

  /**
   *  Calculate the initial evolution scales given momenta
   */
  pair<Energy,Energy> calculateInitialInitialScales(const Lorentz5Momentum& p1, 
						    const Lorentz5Momentum& p2);

  /**
   *  Calculate the initial evolution scales given momenta
   */
  pair<Energy,Energy> calculateInitialFinalScales(const Lorentz5Momentum& pb, const Lorentz5Momentum& pc,
						  const bool isDecayCase);

protected:

  /**
   * Given a pair of particles, supposedly partners w.r.t. an interaction,
   * this method returns their initial evolution scales as a pair.
   * If something wrong happens, it returns the null (ZERO,ZERO) pair. 
   * This method is used by the above setXXXInitialEvolutionScales 
   * methods.
   */
  //@{
  /**
   *  Calculate the initial evolution scales for two final-state particles
   */
  virtual pair<Energy,Energy> calculateFinalFinalScales(const ShowerPPair &);

  /**
   *  Calculate the initial evolution scales for two initial-state particles
   */
  virtual pair<Energy,Energy> calculateInitialInitialScales(const ShowerPPair &);

  /**
   *  Calculate the initial evolution scales for one initial 
   *  and one final-state particles
   */
  virtual pair<Energy,Energy> calculateInitialFinalScales(const ShowerPPair &,
							  const bool isDecayCase);
  //@}

  /**
   * Access function for the initial conditions for the shower
   */
  //@{
  /**
   * Initial conditions for the shower of a final-final colour connection
   * - 0 is the symmetric choice
   * - 1 is maximal emmision from the coloured particle
   * - 2 is maximal emmision from the anticoloured particle
   * - 3 is randomly selected maximal emmision
   */
  unsigned int finalFinalConditions() const 
  {return _finalFinalConditions;}

  /**
   * Initial conditions for the shower of an initial-final decay colour connection
   * - 0 is the symmetric choice
   * - 1 is maximal emission from the decay product
   * - 2 is the smooth choice
   */
  unsigned int initialFinalDecayConditions() const
  {return _initialFinalDecayConditions;}
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeFinder & operator=(const QTildeFinder &) = delete;

private:

  /**
   *  Flags controlling the initial conditions for the shower
   */
  //@{
  /**
   * Initial conditions for the shower with a final-final colour
   * connection
   */
  unsigned int _finalFinalConditions; 

  /**
   * Initial conditions for the shower with an initial-final decay colour
   * connection. This is done according to the top decay colour 
   *  connection calculation in JHEP12(2003)_045. The options act as follows:
   *  0: This is the default 'symmetric' choice which more or less divides
   *     the phase space evenly between the parent and its charged child.
   *  1: This 'maximal' choice maximises the phase space available for 
   *     gluons emitted from the charged child.
   *  2: This (experimental) 'smooth' choice does not suffer from
   *     a discontinuity at the boundary between the region populated by
   *     emissions from the charged child and the region populated by emissions
   *     from the parent. This does, however, mean that the phase space 
   *     available for emissions from the charged child is fairly minimal.
   */
  unsigned int _initialFinalDecayConditions;

  /**
   * Initial conditions for the shower with an initial-initial colour
   * connection. This is done according to the colour connection 
   * calculation in JHEP12(2003)_045. The options act as follows:
   *  0: This is the default 'symmetric' choice which more or less divides
   *     the phase space evenly between the two incoming partons.
   *  1: This increases the phase space for emission from "parton b".
   *  2: This increases the phase space for emission from "parton c".
   */
  unsigned int _initialInitialConditions;
  //@}
};

}

#endif /* HERWIG_QTildeFinder_H */
