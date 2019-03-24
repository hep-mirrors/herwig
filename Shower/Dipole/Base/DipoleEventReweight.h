// -*- C++ -*-
//
// DipoleEventReweight.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_DipoleEventReweight_H
#define Herwig_DipoleEventReweight_H
//
// This is the declaration of the DipoleEventReweight class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/StandardModel/AlphaSBase.h"


namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 * \brief Reweight full final states produced by the shower
 *
 * @see \ref DipoleEventReweightInterfaces "The interfaces"
 * defined for DipoleEventReweight.
 */
class DipoleEventReweight: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleEventReweight();

  /**
   * The destructor.
   */
  virtual ~DipoleEventReweight();
  //@}

public:

  /**
   * Return the weight for the given incoming, outgoing coloured and
   * hard colour neutral particles, after an emission was generated
   */
  virtual double weight(const PPair& in, const PList& out, const PList& hard,
			Ptr<AlphaSBase>::tptr as) const = 0;

  /**
   * Return the weight which is applied to a cascade only once even if there was
   * no emission from the cascade at all. 
   */
  virtual double weightCascade(const PPair& in, const PList& out, const PList& hard,
				  Ptr<AlphaSBase>::tptr as) const = 0;

  /**
   * Return true, if the event reweight should be applied to the hard scattering
   */
  virtual bool firstInteraction() const = 0;

  /**
   * Return true, if the event reweight should be applied to secondary interactions
   */
  virtual bool secondaryInteractions() const = 0;

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleEventReweight & operator=(const DipoleEventReweight &) = delete;

};

}

#endif /* Herwig_DipoleEventReweight_H */
