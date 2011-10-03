// -*- C++ -*-
//
// ColourReconnector.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
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
 *  \author Alberto Ribon, Christian Roehr
 * 
 *  This class does the nonperturbative colour rearrangement, after the 
 *  nonperturbative gluon splitting and the "normal" cluster formation. 
 *  It uses the list of particles in the event record, and the collections of
 *  "usual" clusters which is passed to the main method. If the colour 
 *  reconnection is actually accepted, then the previous collections of "usual"
 *  clusters is first deleted and then the new one is created.
 *
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
  ColourReconnector() : _clreco( 0 ), _preco( 0.5 )
  {}
  //@}

  /**
   * Does the colour rearrangement, starting out from the list of particles in
   * the event record and the collection of "usual" clusters passed as
   * arguments. If the actual rearrangement is accepted, the initial collection of
   * clusters is overridden by the old ones.
   */
  void rearrange(EventHandler & ch,
                 ClusterVector & clusters);
    
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
   * Private and non-existent assignment operator.
   */
  ColourReconnector & operator=(const ColourReconnector &);

  /**
   * Do we do colour reconnections?
   */
  int _clreco;

  /**
   * Probability that a found reconnection possibility is actually accepted.
   */
  double _preco;


private:

  /**
   * Returns the Cluster (within the ClusterVector cv) where the sum of the
   * invariant Cluster masses becomes minimal in the case of a colour
   * reconnection with cl. If no reconnection partner can be found, a pointer to
   * the original Cluster cl is returned.
   */
  ClusterPtr _findRecoPartner(ClusterPtr cl, ClusterVector cv) const;

  /**
   * @return	true, if the two partons are splitting products of the same
   * 		gluon
   */
  bool _isColour8(tPPtr p1, tPPtr p2) const;

  /**
   * Reconnects the constituents of the given clusters to the (only) other
   * possible cluster combination.
   */
  pair <ClusterPtr,ClusterPtr> _reconnect(ClusterPtr c1, ClusterPtr c2) const;

};


}

#endif /* HERWIG_ColourReconnector_H */
