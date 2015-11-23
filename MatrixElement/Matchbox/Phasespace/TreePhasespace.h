// -*- C++ -*-
//
// TreePhasespace.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_TreePhasespace_H
#define Herwig_TreePhasespace_H
//
// This is the declaration of the TreePhasespace class.
//

#include "Herwig/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/TreePhasespaceChannels.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer, Ken Arnold
 *
 * \brief TreePhasespace is a multichannel phasespace generator
 * adapting to singularity structures as determined from the matrix
 * elements diagrams.
 *
 * @see \ref TreePhasespaceInterfaces "The interfaces"
 * defined for TreePhasespace.
 */
class TreePhasespace: public MatchboxPhasespace {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  TreePhasespace();

  /**
   * The destructor.
   */
  virtual ~TreePhasespace();
  //@}

public:

  /**
   * Prepare a phase space generator for the given xcomb object.
   */
  virtual void setXComb(tStdXCombPtr);

  /**
   * Generate a phase space point and return its weight.
   */
  virtual double generateTwoToNKinematics(const double*,
					  vector<Lorentz5Momentum>& momenta);

  /**
   * Return the number of random numbers required to produce a given
   * multiplicity final state.
   */
  virtual int nDimPhasespace(int nFinal) const {
    if ( nFinal == 1 )
      return 1;
    return 3*(nFinal - 1); // one additional number needed for channel selection
  }

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

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}


private:

  /**
   * The object storing channel maps
   */
  Ptr<TreePhasespaceChannels>::ptr theChannelMap;

  /**
   * Map xcomb's to channel vectors indexed by diagram id.
   */
  map<tStdXCombPtr,
      map<Ptr<Tree2toNDiagram>::ptr,
	  pair <PhasespaceHelpers::PhasespaceTree, PhasespaceHelpers::PhasespaceTree> > >&
  channelMap() { return theChannelMap->channelMap(); }

  /**
   * The currently active channels.
   */
  map<tStdXCombPtr,
      map<Ptr<Tree2toNDiagram>::ptr,
	  pair <PhasespaceHelpers::PhasespaceTree, PhasespaceHelpers::PhasespaceTree> > >::iterator 
  lastChannelsIterator;

  /**
   * The phasespace info object to be used.
   */
  PhasespaceHelpers::PhasespaceInfo lastPhasespaceInfo;

  /**
   * Parameter steering from which on propagator virtualities are
   * sampled flat.
   */
  double x0;

  /**
   * Parameter steering at which virtuality singularities of
   * propagators are actually cut off.
   */
  double xc;

  /**
   * Parameter steering from which on propagator virtualities are
   * sampled flat.
   */
  Energy M0;

  /**
   * Parameter steering at which virtuality singularities of
   * propagators are actually cut off.
   */
  Energy Mc;

  /**
   * Choose whether to also use mirrored phasespace generation
   */
  bool theIncludeMirrored;
       
  /**
   * Return the currently active channels.
   */
  map<Ptr<Tree2toNDiagram>::ptr,
      pair <PhasespaceHelpers::PhasespaceTree, PhasespaceHelpers::PhasespaceTree> >& lastChannels() { 
    return lastChannelsIterator->second; 
  }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TreePhasespace & operator=(const TreePhasespace &);

};

}

#endif /* Herwig_TreePhasespace_H */
