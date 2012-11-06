// -*- C++ -*-
//
// TreePhasespaceChannels.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_TreePhasespaceChannels_H
#define Herwig_TreePhasespaceChannels_H
//
// This is the declaration of the TreePhasespaceChannels class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/PhasespaceHelpers.h"

namespace Herwig {

using namespace ThePEG;

/**
 *
 * \ingroup Matchbox
 * \author Simon Platzer, Ken Arnold
 *
 * \brief Store channels for the tree phasespace.
 *
 * @see \ref TreePhasespaceChannelsInterfaces "The interfaces"
 * defined for TreePhasespaceChannels.
 */
class TreePhasespaceChannels: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  TreePhasespaceChannels();

  /**
   * The destructor.
   */
  virtual ~TreePhasespaceChannels();
  //@}

public:

  /**
   * Access the channel map
   */
  map<tStdXCombPtr,map<Ptr<Tree2toNDiagram>::ptr,pair <PhasespaceHelpers::PhasespaceTree, PhasespaceHelpers::PhasespaceTree> > >&
  channelMap() { return theChannelMap; }

  /**
   * Return the channel map
   */
  const map<tStdXCombPtr,map<Ptr<Tree2toNDiagram>::ptr,pair <PhasespaceHelpers::PhasespaceTree, PhasespaceHelpers::PhasespaceTree> > >&
  channelMap() const { return theChannelMap; }

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
   * Map xcomb's to channel vectors indexed by diagram id.
   */
  map<tStdXCombPtr,map<Ptr<Tree2toNDiagram>::ptr,pair<PhasespaceHelpers::PhasespaceTree,PhasespaceHelpers::PhasespaceTree> > > theChannelMap;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TreePhasespaceChannels & operator=(const TreePhasespaceChannels &);

};

}

#endif /* Herwig_TreePhasespaceChannels_H */
