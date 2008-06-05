// -*- C++ -*-
#ifndef HERWIG_NasonCKKWHandler_H
#define HERWIG_NasonCKKWHandler_H
//
// This is the declaration of the NasonCKKWHandler class.
//

#include "Herwig++/Shower/ShowerHandler.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Utilities/Interpolator.h"
#include "Herwig++/Utilities/GaussianIntegrator.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/Config/Pointers.h"
#include "NasonTree.h"

namespace Herwig {

  class NasonCKKWHandler;

}
//declaration of thepeg ptr
namespace ThePEG {

  ThePEG_DECLARE_POINTERS(Herwig::NasonCKKWHandler,NasonCKKWHandlerPtr);

}

namespace Herwig {

using namespace ThePEG;


/**
 * Here is the documentation of the NasonCKKWHandler class.
 *
 * @see \ref NasonCKKWHandlerInterfaces "The interfaces"
 * defined for NasonCKKWHandler.
 */
class NasonCKKWHandler: public ShowerHandler {

public:

  /**
   * The default constructor.
   */
  NasonCKKWHandler() : _npoint(200) {}

  /**
   *  Destructor
   */
  virtual ~NasonCKKWHandler();

  /**
   * Perform CKKW reweighting
   */
  virtual double reweightCKKW(int minMult, int maxMult);

  /**
   *  Main method for the cascade
   */
  virtual void cascade();

public:

  /**
   *  access to the nason tree object
   */
  inline NasonTreePtr getNasonTree(){
    return _theNasonTree;
  }

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

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
   * Clusters the partons and creates a branching history
   * by combining the 2 particles with smallest
   * jet measure out of all allowed pairings until we are left 
   * with \f$q\bar{q}\f$.
   */
  bool doClustering( ParticleVector theParts, PPtr vb );
  
  /**
   * Checks to see that the splitting is allowed.
   */
  bool splittingAllowed( ShowerParticlePtr part_i,
			 ShowerParticlePtr part_j,
			 int qq_pairs);
  
  
  /**
   * Checks to see that the splitting is allowed and finds the
   * Sudakov form factor for the splitting.
   */
  bool getSud( int * qq_pairs, long * emmitter_id,
	       SudakovPtr clusterSudakov,
	       ShowerParticlePtr part_i, 
	       ShowerParticlePtr part_j ) ;
  
  /**
   * Returns the durham jet measure, yij, for the two particles. 
   */
  double getJetMeasure(ShowerParticlePtr part_i, ShowerParticlePtr part_j);
  
  
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<NasonCKKWHandler> initNasonCKKWHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NasonCKKWHandler & operator=(const NasonCKKWHandler &);

private:

  /**
   *  The nason tree
   */
  NasonTreePtr _theNasonTree;

  /**
   *  Number of points for the interpolation tables
   */
  unsigned int _npoint;

  /**
   *  Centre of mass energy
   */
  Energy2 _s;

  /**
   *  Map containing the sudakovs for the final-state particles
   */
  multimap<long,pair<Interpolator<double,Energy>::Ptr,Interpolator<Energy,double>::Ptr> > _fbranchings;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NasonCKKWHandler. */
template <>
struct BaseClassTrait<Herwig::NasonCKKWHandler,1> {
  /** Typedef of the first base class of NasonCKKWHandler. */
  typedef Herwig::ShowerHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NasonCKKWHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NasonCKKWHandler>
  : public ClassTraitsBase<Herwig::NasonCKKWHandler> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::NasonCKKWHandler"; }
  /**
   * The name of a file containing the dynamic library where the class
   * NasonCKKWHandler is implemented. It may also include several, space-separated,
   * libraries if the class NasonCKKWHandler depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNasonShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_NasonCKKWHandler_H */


