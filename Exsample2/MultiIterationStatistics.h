// -*- C++ -*-
//
// MultiIterationStatictis.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MultiIterationStatistics_H
#define Herwig_MultiIterationStatistics_H
//
// This is the declaration of the MultiIterationStatistics class.
//

#include "GeneralStatistics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Monte Carlo statistics for multiple iterations
 */
class MultiIterationStatistics: public Herwig::GeneralStatistics {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MultiIterationStatistics();

  /**
   * The destructor.
   */
  virtual ~MultiIterationStatistics();
  //@}

public:

  /**
   * Indicate the start of a new iteration
   */
  void nextIteration() {
    iterations().push_back(GeneralStatistics(*this));
    reset();
  }

  /**
   * Return the iterations done so far.
   */
  const vector<GeneralStatistics>& iterations() const {
    return theIterations;
  }

  /**
   * Access the iterations done so far.
   */
  vector<GeneralStatistics>& iterations() {
    return theIterations;
  }

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void put(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void get(PersistentIStream & is);
  //@}

private:

  /**
   * The currently accumulated iterations.
   */
  vector<GeneralStatistics> theIterations;

};

inline PersistentOStream& operator<<(PersistentOStream& os, const MultiIterationStatistics& s) {
  s.put(os); return os;
}

inline PersistentIStream& operator>>(PersistentIStream& is, MultiIterationStatistics& s) {
  s.get(is); return is;
}

}

#endif /* Herwig_MultiIterationStatistics_H */
