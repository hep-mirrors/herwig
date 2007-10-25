// -*- C++ -*-
#ifndef HERWIG_ClusteringSelector_H
#define HERWIG_ClusteringSelector_H
//
// This is the declaration of the ClusteringSelector class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ClusteringSelector.fh"

#include "Clustering.h"
#include "ClusteringGuide.h"

namespace Herwig {

using namespace ThePEG;

  /**\ingroup CKKW
   *
   * The binary predicate used for comparing two
   * clusterings amongst their scale.
   *
   *@author Simon Plaetzer
   *
   */
  struct ClusteringScaleLess {

    /**
     * Return wether the first clustering scale is less than
     * the second clustering scale.
     */
    inline bool operator () (const pair<ClusteringPtr,tClusteringGuidePtr>& first,
			     const pair<ClusteringPtr,tClusteringGuidePtr>& second)
    { return first.first->scale() < second.first->scale() ; }

  };

/**\ingroup CKKW
 *
 * ClusteringSelector is used to order possible clusterings
 * in each step of a cascade history reconstruction. The default
 * implementation sorts the clusterings according to the scales
 * associated with each clustering. In case there are degenerate
 * clusterings from the same configuration with the same scale,
 * a clustering out of these is selested randomly according to
 * the weight set by the clusterer used.
 *
 *@author Simon Plaetzer
 *
 * @see Clusterer
 *
 * @see \ref ClusteringSelectorInterfaces "The interfaces"
 * defined for ClusteringSelector.
 */
class ClusteringSelector: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ClusteringSelector();

  /**
   * The destructor.
   */
  virtual ~ClusteringSelector();
  //@}

public:

  /**
   * Order the given clusterings.
   *
   * The default version orders
   * according to the scales associated with a clustering. In case
   * there are different clusterings with the same scale,
   * selectDegenerate is invoked.
   */
  virtual list<pair<ClusteringPtr, tClusteringGuidePtr> > sort
  (const list<pair<ClusteringPtr, tClusteringGuidePtr> >&) const;

  /**
   * Select a clustering out of a set of degenerate ones (in scale).
   * Default implementation selects randomly according to the weight.
   */
  virtual pair<ClusteringPtr,tClusteringGuidePtr> 
  selectDegenerate (const list<pair<ClusteringPtr,tClusteringGuidePtr> >&) const;

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ClusteringSelector> initClusteringSelector;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ClusteringSelector & operator=(const ClusteringSelector &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ClusteringSelector. */
template <>
struct BaseClassTrait<Herwig::ClusteringSelector,1> {
  /** Typedef of the first base class of ClusteringSelector. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ClusteringSelector class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ClusteringSelector>
  : public ClassTraitsBase<Herwig::ClusteringSelector> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ClusteringSelector"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ClusteringSelector is implemented. It may also include several, space-separated,
   * libraries if the class ClusteringSelector depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "ClusteringSelector.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ClusteringSelector.tcc"
#endif

#endif /* HERWIG_ClusteringSelector_H */
