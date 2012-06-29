// -*- C++ -*-
//
// ShowerModel.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerModel_H
#define HERWIG_ShowerModel_H
//
// This is the declaration of the ShowerModel class.
//
#include "ThePEG/Interface/Interfaced.h"
#include "KinematicsReconstructor.fh"
#include "PartnerFinder.fh"
#include "SudakovFormFactor.fh"
#include "ShowerModel.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  The ShowerModel class is a container for all the objects needed to implement a
 *  specific model of the shower evolution, as opposed to those which are independent
 *  of the evolution.
 *
 *  In general there are four types of object 
 * - The KinematicsReconstructor object which is responsible for reconstruction
 *   of the shower kinematics after the evolution.
 * - The PartnerFinder which is responsible for finding the partner and setting the
 *   initial evolution scale
 * - A vector of SudakovFormFactor objects which will usually all be instances
 *   of a class implementing the SudakovFormFactor for a specific model with
 *   different splitting functions for different branchings
 *
 *  For each model the checkConsistency member must be implemented to check that
 *  the correct objects for the model are used.
 *
 * @see \ref ShowerModelInterfaces "The interfaces"
 * defined for ShowerModel.
 */
class ShowerModel: public Interfaced {

public:  

  /**
   *  Access methods to access the objects
   */
  //@{
  /**
   *  Access to the KinematicsReconstructor object
   */
  tKinematicsReconstructorPtr kinematicsReconstructor() const { return _reconstructor; }

  /**
   *  Access to the PartnerFinder object
   */
  tPartnerFinderPtr partnerFinder() const { return _partnerfinder; }

  /**
   *  Access to the SudakovFormFactor objects
   */
  const vector<SudakovPtr> & sudakovFormFactors() const { return _sudakovs; }
  //@}

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();


protected:

  /**
   *  The checkConsitency member which must be implemented in classes
   *  inheriting from this one.
   */
  virtual void checkConsistency() =0;

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<ShowerModel> initShowerModel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerModel & operator=(const ShowerModel &);

private:

  /**
   *  Pointer to the various objects
   */
  //@{
  /**
   *  Pointer to the KinematicsReconstructor object
   */
  KinematicsReconstructorPtr _reconstructor;

  /**
   *  Pointer to the PartnerFinder object
   */
  PartnerFinderPtr _partnerfinder;

  /**
   *  Pointers to the SudakovFormFactor objects
   */
  vector<SudakovPtr> _sudakovs;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ShowerModel. */
template <>
struct BaseClassTrait<Herwig::ShowerModel,1> {
  /** Typedef of the first base class of ShowerModel. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ShowerModel class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ShowerModel>
  : public ClassTraitsBase<Herwig::ShowerModel> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ShowerModel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ShowerModel is implemented. It may also include several, space-separated,
   * libraries if the class ShowerModel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ShowerModel_H */
