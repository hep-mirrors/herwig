// -*- C++ -*-
//
// MatchboxMECache.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxMECache_H
#define Herwig_MatchboxMECache_H
//
// This is the declaration of the MatchboxMECache class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/LastXCombInfo.h"

#include <boost/functional/hash.hpp>

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxMECache provides caching for matrix elements.
 *
 */
class MatchboxMECache: public HandlerBase, public LastXCombInfo<StandardXComb> {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxMECache();

  /**
   * The destructor.
   */
  virtual ~MatchboxMECache();
  //@}

public:

  /**
   * Set the XComb object.
   */
  virtual void setXComb(tStdXCombPtr xc) { 
    theLastXComb = xc;
  }

  /**
   * Flush the cache.
   */
  void flush() { theME2Cache.clear(); }

  /**
   * Caluclate a hash value for the phase space
   * point contained in meMomenta(). The hash value
   * is calculated using only the spatial components
   * of the outgoing partons, thus assuming momentum
   * conservation and definite and constant mass shells
   * for each outgoing particle.
   */
  size_t hashPhaseSpace() const;

  /**
   * Return true, if the matrix element needs to be 
   * recalculated for the given phase space point.
   * If not, return the cached value in the given reference.
   */
  bool calculateME2(double& xme2,
		    const pair<int,int>& corr = make_pair(0,0));

  /**
   * Cache a calculated matrix element
   * for the last phase space point.
   */
  void cacheME2(double xme2,
		const pair<int,int>& corr = make_pair(0,0));

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

private:

  /**
   * Helper for cache indexing.
   */
  struct MECacheKey {

    /**
     * The phasespace hash
     */
    size_t phasespaceHash;

    /**
     * The process.
     */
    const cPDVector& process;

    /**
     * Initialize the ref.
     */
    static const cPDVector& initProcess() {
      static cPDVector value; return value;
    }

    /**
     * Additional id to identify correlated ME's squared.
     */
    pair<int,int> corrIds;

    /**
     * Default constructor.
     */
    MECacheKey()
      : phasespaceHash(0), process(initProcess()),
	corrIds(0,0) {}

    /**
     * Standard constructor.
     */
    MECacheKey(size_t pshash, const cPDVector& proc,
	       const pair<int,int>& cids)
      : phasespaceHash(pshash), process(proc),
	corrIds(cids) {}

    /**
     * Compare for equality
     */
    inline bool operator==(const MECacheKey& x) const {
      return
	phasespaceHash == x.phasespaceHash &&
	process == x.process &&
	corrIds == x.corrIds;
    }

    /**
     * Compare for ordering
     */
    inline bool operator<(const MECacheKey& x) const {
      if ( phasespaceHash != x.phasespaceHash )
	return phasespaceHash < x.phasespaceHash;
      if ( process != x.process )
	return process < x.process;
      return corrIds < x.corrIds;
    }

  };

  /**
   * Map process and phase space hashes to already calculated matrix
   * elements. @todo This goes to XComb::meta ... but who does it?
   */
  map<MECacheKey,double> theME2Cache;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxMECache & operator=(const MatchboxMECache &);

};

}

#endif /* Herwig_MatchboxMECache_H */
