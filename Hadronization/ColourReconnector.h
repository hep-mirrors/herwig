// -*- C++ -*-
//
// ColourReconnector.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ColourReconnector_H
#define HERWIG_ColourReconnector_H

#include <ThePEG/Interface/Interfaced.h>
#include "CluHadConfig.h"
#include "ColourReconnector.fh"


namespace Herwig {

using namespace ThePEG;

/** \ingroup Hadronization
 *  \class ColourReconnector
 *  \brief Class for changing colour reconnections of partons.
 *  \author Alberto Ribon
 * 
 *  This class does the nonperturbative colour rearrangement, after the 
 *  nonperturbative gluon splitting and the "normal" cluster formation. 
 *  It uses the list of particles in the event record, and the collections of
 *  "usual" clusters which is passed to the main method. If the colour 
 *  reconnection is actually accepted, then the previous collections of "usual"
 *  clusters is first deleted and then the new one is created.
 *
 *  Note: by default this class does nothing. It can be inherited and overridden
 *  in future hadronization models.
 * * @see \ref ColourReconnectorInterfaces "The interfaces"
 * defined for ColourReconnector.
 */
class ColourReconnector: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline ColourReconnector();
  //@}

  /**
   * Does the colour rearrangment.
   *
   * Does the colour rearrangement, starting from the list of particles
   * in the event record, and the collection of "usual" clusters passed
   * in input. If the actual rearrangement is accepted, the new collection 
   * of clusters is overriden to the intial one.
   */
  void rearrange(EventHandler & ch,
                 ClusterVector & clusters) throw(Veto, Stop, Exception);
    
public:

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ColourReconnector> initColourReconnector;

  /**
   * Private and non-existent assignment operator.
   */
  ColourReconnector & operator=(const ColourReconnector &);

  /**
   * Do we do colour reconnections?
   */
  int _clreco;
};


}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

template <>
/**
 * The following template specialization informs ThePEG about the
 * base class of ColourReconnector.
 */
struct BaseClassTrait<Herwig::ColourReconnector,1> {
  /** Typedef of the base class of ColourReconnector. */
  typedef Interfaced NthBase;
};

template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
struct ClassTraits<Herwig::ColourReconnector>
  : public ClassTraitsBase<Herwig::ColourReconnector> {
  /** Return the class name.*/
  static string className() { return "Herwig::ColourReconnector"; }
};

/** @endcond */

}


#include "ColourReconnector.icc"

#endif /* HERWIG_ColourReconnector_H */
