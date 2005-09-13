// -*- C++ -*-
#ifndef HERWIG_ClusterDecayer_H
#define HERWIG_ClusterDecayer_H

#include <ThePEG/Handlers/HandlerBase.h>
#include <ThePEG/EventRecord/Step.h>
#include "CluHadConfig.h"
#include "HadronSelector.h"
#include "Herwig++/Utilities/GlobalParameters.h"


namespace Herwig {
using namespace ThePEG;

class ThePEG::Particle;   // forward declaration

/** \ingroup Hadronization
 *  \class ClusterDecayer
 *  \brief This class decays the "normal" clusters
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *
 *  This class decays the "normal" clusters, e.g. ones that are not heavy
 *  enough for fission, and not too light to decay into one hadron. 
 *
 *  This class is directs the production of hadrons via 2-body cluster decays. 
 *  The selection of the hadron flavours is given by Herwig::HadronSelector.
 *
 *  @see HadronSelector
 */
class ClusterDecayer: public ThePEG::HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline ClusterDecayer();  
  
  /**
   * Copy constructor.
   */
  inline ClusterDecayer(const ClusterDecayer &);

  /**
   * Standard destructor.
   */
  virtual ~ClusterDecayer();
  //@}

  /** Decays all remaining clusters into hadrons. 
   * This routine decays the clusters that are left after
   * Herwig::ClusterFissioner::fission and
   * Herwig::LightClusterDecayer::decay have been called. These are all
   * the "normal" clusters which are not forced into hadrons by
   * the other functions.
   */
  void decay(const StepPtr&) 
    throw(Veto, Stop, Exception);

public:

  /**
   * Standard ThePEG function for writing a persistent stream.
   */
  void persistentOutput(PersistentOStream &) const;

  /**
   * Standard ThePEG function for reading from a persistent stream.
   */
  void persistentInput(PersistentIStream &, int);

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
  static ClassDescription<ClusterDecayer> initClusterDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  ClusterDecayer & operator=(const ClusterDecayer &);

public:

  /** Decays the cluster into two hadrons. 
   *
   *  This routine is used to take a given cluster and decay it into
   *  two hadrons which are returned. If one of the constituents is from
   *  the perturbative regime then the direction of the perturbative parton
   *  is remembered and the decay is preferentially in that direction. The 
   *  direction of the decay is given by
   *  \f[ \cos \theta = 1 + S \log r_1 \f]
   *  where \f$ S \f$ is a parameter of the model and \f$ r_1 \f$ is a random
   *  number [0,1].
   */
  pair<PPtr,PPtr> decayIntoTwoHadrons(tClusterPtr ptr) 
    throw(Veto, Stop, Exception);

private:

  /** Compute the positions of the new hadrons based on the clusters position.
   *
   *  This method calculates the positions of the children hadrons by a
   *  call to Herwig::Smearing::gaussianSmearing with width inversely 
   *  proportional to the cluster mass, around the parent cluster position.
   */
  void calculatePositions( const Lorentz5Momentum &, const LorentzPoint &, 
			   const Lorentz5Momentum &, const Lorentz5Momentum &,
			   LorentzPoint &, LorentzPoint &) const;

  /**
   * Pointer to a Herwig::HadronSelector for choosing decay types
   */
  Ptr<HadronSelector>::pointer _hadronsSelector;

  /**
   * Pointer to a Herwig::GlobalParameters for various global data
   */
  Ptr<GlobalParameters>::pointer _globalParameters;
  
  /**
   * Whether non-b clusters decay along the perturbative parton direction.
   */
  int _ClDir1;

  /**
   * Whether b clusters decay along the perturbative parton direction.
   */
  int _ClDir2;

  /**
   * The S parameter from decayIntoTwoHadrons for non-b clusters.
   */
  double _ClSmr1;

  /**
   * The S parameter from decayIntoTwoHadrons for b clusters.
   */
  double _ClSmr2;

  /**
   * Whether or not the hadrons produced should be on-shell
   * or generated used the MassGenerator
   */
  bool _onshell;

  /**
   * Number of tries to generate the masses of the decay products
   */
  unsigned int _masstry;

};


}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {

template <>
/** 
 * The following template specialization informs ThePEG about the
 * base class of ClusterDecayer.
 */
struct BaseClassTrait<Herwig::ClusterDecayer,1> {
  /** Typedef of the base class of ClusterDecayer. */
  typedef ThePEG::HandlerBase NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
struct ClassTraits<Herwig::ClusterDecayer>
  : public ClassTraitsBase<Herwig::ClusterDecayer> {
  /** Return the class name. */
  static string className() { return "Herwig++::ClusterDecayer"; }
  /** Return the name of the shared library to be loaded to get
   *  access to this class and every other class it uses
   *  (except the base class).
   */
  static string library() { return "libHwHadronization.so"; }
};

}

#endif // DOXYGEN

#include "ClusterDecayer.icc"

#endif /* HERWIG_ClusterDecayer_H */
