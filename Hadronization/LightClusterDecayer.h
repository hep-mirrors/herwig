// -*- C++ -*-
#ifndef HERWIG_LightClusterDecayer_H
#define HERWIG_LightClusterDecayer_H

#include <ThePEG/Handlers/HandlerBase.h>
#include <ThePEG/EventRecord/Step.h>
#include "CluHadConfig.h"
#include "HadronSelector.h"


namespace Herwig {


using namespace ThePEG;

/** \ingroup Hadronization
 *  \class LightClusterDecayer
 *  \brief This class performs the decay of light clusters into a single hadron.
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *
 *  This is the class that performs the decay of light clusters into 
 *  only one hadron. The major difficulty is that a kinematical reshuffling
 *  is necessary, between the cluster under consideration and its
 *  "neighbouring" clusters, to conserve energy-momentum in one-body decay.
 *  Notice that, differently from what happens in Fortran Herwig,
 *  light (that is below the threshold for the production of the lightest 
 *  pair of hadrons with the proper flavours) fission products, produced 
 *  by the fission of heavy clusters in ClusterFissioner 
 *  have been already "decayed" into single hadron (the lightest one 
 *  with proper flavour) by the same latter class, without require 
 *  any reshuffling. Therefore the light clusters that are treated in 
 *  this LightClusterDecayer class are produced directly 
 *  (originally) by the ClusterFinder. 
 *	
 *  Notice:
 *  - The choice of the candidate cluster with whom to reshuffle momentum 
 *    is based on the minimal space-time distance from the light cluster.
 *  - An alternate choice of what is considered a "neighbour" could be
 *    implemented but was not considered for Herwig++.
 *
 *  @see HadronSelector
 */ 
class LightClusterDecayer: public ThePEG::HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline LightClusterDecayer();

  /**
   * Copy-constructor.
   */
  inline LightClusterDecayer(const LightClusterDecayer &);

  /**
   * Destructor.
   */
  virtual ~LightClusterDecayer();
  //@}

  /**
   * This method does the decay of light hadron in one hadron.
   *
   * This method requires a kinematical reshuffling for energy-momentum 
   * conservation. This is done explicitly by the (private) method 
   * reshuffling().
   */
  bool decay(const StepPtr &);

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<LightClusterDecayer> initLightClusterDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  LightClusterDecayer & operator=(const LightClusterDecayer &);

  /**
   * This (private) method, called by decay(), takes care of the kinematical
   * reshuffling necessary for energy-momentum conservation.
   */
  bool reshuffling( const long, tClusterPtr, tClusterPtr,
		    const StepPtr &, tClusterVector &) 
    throw (Veto, Stop, Exception); 
  
  /**
   * A pointer to a Herwig::HadronSelector object used for producing hadrons.
   */
  Ptr<HadronSelector>::pointer _hadronsSelector;

  /**
   * A parameter used for determining when b clusters are too light.
   *
   * This parameter is used for setting the lower threshold, \f$ t \f$ as
   * \f[ t' = t(1 + r B^1_{\rm lim}) \f]
   * where \f$ r \f$ is a random number [0,1].
   */
  double _B1Lim;

};

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {

template <>
/**
 * The following template specialization informs ThePEG about the
 * base class of LightClusterDecayer.
 */
struct BaseClassTrait<Herwig::LightClusterDecayer,1> {
  /** Typedef of the base class of LightClusterDecayer. */
  typedef ThePEG::HandlerBase NthBase;
};

template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
struct ClassTraits<Herwig::LightClusterDecayer>: 
    public ClassTraitsBase<Herwig::LightClusterDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig++::LightClusterDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwHadronization.so"; }
};

}

#endif // DOXYGEN STUFF

#include "LightClusterDecayer.icc"

#endif /* HERWIG_LightClusterDecayer_H */
