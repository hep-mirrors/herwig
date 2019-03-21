// -*- C++ -*-
//
// MatchboxXCombGroup.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxXCombGroup_H
#define Herwig_MatchboxXCombGroup_H
//
// This is the declaration of the MatchboxXCombGroup class.
//

#include "ThePEG/Handlers/StdXCombGroup.h"
#include "ThePEG/MatrixElement/MEGroup.h"
#include "Herwig/MatrixElement/Matchbox/Utility/MatchboxXCombData.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Matchbox extensions to StandardXComb
 */
class MatchboxXCombGroup: public StdXCombGroup, public MatchboxXCombData {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Standard constructor.
   */
  MatchboxXCombGroup(Energy newMaxEnergy, const cPDPair & inc,
		     tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
		     tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
		     const PBPair & newPartonBins, tCutsPtr newCuts, tMEGroupPtr newME,
		     const DiagramVector & newDiagrams, bool mir,
		     tStdXCombPtr newHead = tStdXCombPtr());

  /**
   * Default constructor.
   */
  MatchboxXCombGroup();

  /**
   * Destructor.
   */
  virtual ~MatchboxXCombGroup();

  //@}

public:

  /**
   * Reset all saved data about last generated phasespace point;
   */
  virtual void clean();

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxXCombGroup & operator=(const MatchboxXCombGroup &) = delete;

};

}

#endif /* Herwig_MatchboxXCombGroup_H */
