// -*- C++ -*-
//
// HEJFactory.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_HEJFactory_H
#define Herwig_HEJFactory_H
//
// This is the declaration of the HEJFactory class.
//

#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "HEJ/Jets.h"

namespace Herwig {

using namespace ThePEG;

using ::CMultijet;

class HEJMEBase;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief HEJFactory sets up HEJ matrix elements.
 *
 * @see \ref HEJFactoryInterfaces "The interfaces"
 * defined for HEJFactory.
 */
class HEJFactory: public SubProcessHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  HEJFactory();

  /**
   * The destructor.
   */
  virtual ~HEJFactory();
  //@}

public:

  /**
   * Return true if this object needs to be initialized before all
   * other objects (except those for which this function also returns
   * true).  This default version always returns false, but subclasses
   * may override it to return true.
   */
  virtual bool preInitialize() const { return true; }

  /**
   * Ask the factory to create a CMultijet object and pass it to the
   * individual matrix elements.
   */
  void makeCMultijet(tStdXCombPtr);

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

protected:

  /**
   * Clone the HEJ ME and bind it to the given incoming
   * flavours and number of gluons.
   */
  void makeME(cPDPtr ini, cPDPtr inj, unsigned int n);

  /**
   * Setup the generation flags.
   */
  SGeneratorFlags setup(tStdXCombPtr) const;

private:

  /**
   * The HEJ matrix element to be used.
   */
  Ptr<HEJMEBase>::ptr theHEJME;

  /**
   * The set of incoming partons to be considered.
   */
  PDVector theIncomingFlavours;

  /**
   * The maximum number of gluons to take into account.
   */
  unsigned int theNGluons;

  /**
   * The minimum required jet pt
   */
  Energy theJetPTMin;

  /**
   * The minimum required extremal parton jet
   */
  Energy theExtremalPTMin;

  /**
   * The scale choice
   */
  int theScaleChoice;

  /**
   * A fixed scale, if fixed scale is chosen
   */
  Energy theFixedScale;

  /**
   * The jets object to be used.
   */
  CMultijet* theJets;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HEJFactory & operator=(const HEJFactory &);

};

}

#endif /* Herwig_HEJFactory_H */
