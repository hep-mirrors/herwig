// -*- C++ -*-
#ifndef HERWIG_ClusterDecayer_H
#define HERWIG_ClusterDecayer_H
/*! \class Herwig::ClusterDecayer ClusterDecayer.h "Herwig++\Hadronization\ClusterDecayer.h"
 *  \brief This class decays the "normal" clusters
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *  \ingroup Hadronization
 *
 * This class decays the "normal" clusters, e.g. ones that are not heavy
 * enough for fission, and not too light to decay into one hadron. 
 *
 * This class is directs the production of hadrons via 2-body cluster decays. 
 * The selection of the hadron flavours is given by Herwig::HadronSelector.
 *
 * See also:
 * HadronSelector.h.
 */

#include <ThePEG/Handlers/HandlerBase.h>
#include <ThePEG/EventRecord/Step.h>
#include "CluHadConfig.h"
#include "HadronSelector.h"
#include "Herwig++/Utilities/GlobalParameters.h"


namespace Herwig {


using namespace ThePEG;

class ThePEG::Particle;   // forward declaration

class ClusterDecayer: public ThePEG::HandlerBase {

public:

  inline ClusterDecayer();
  inline ClusterDecayer(const ClusterDecayer &);
  virtual ~ClusterDecayer();
  //!< Standard destructor.

  void decay(const StepPtr&) 
    throw(Veto, Stop, Exception);
  /*!< Decays all remaining clusters into hadrons. 
   *
   * This routine decays the clusters that are left after
   * Herwig::ClusterFissioner::fission and
   * Herwig::LightClusterDecayer::decay have been called. These are all
   * the "normal" clusters which are not forced into hadrons by
   * the other functions.
   */

public:

  void persistentOutput(PersistentOStream &) const;
  //!< Standard ThePEG function for writing a persistent stream.
  void persistentInput(PersistentIStream &, int);
  //!< Standard ThePEG function for reading from a persistent stream.

  static void Init();
  //!< Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  //!< Standard clone methods.

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  //!< Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  //!< Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  //!< Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<ClusterDecayer> initClusterDecayer;
  //!< Describe a concrete class with persistent data.

  ClusterDecayer & operator=(const ClusterDecayer &);
  //!<  Private and non-existent assignment operator.

public:
  pair<PPtr,PPtr> decayIntoTwoHadrons(tClusterPtr ptr) 
    throw(Veto, Stop, Exception);
  /*!< Decays the cluster into two hadrons. 
   *
   * This routine is used to take a given cluster and decay it into
   * two hadrons which are returned. If one of the constituents is from
   * the perturbative regime then the direction of the perturbative parton
   * is remembered and the decay is preferentially in that direction. The 
   * direction of the decay is given by
   * \f[ \cos \theta = 1 + S \log r_1 \f]
   * where \f$ S \f$ is a parameter of the model and \f$ r_1 \f$ is a random
   * number [0,1].
   */

private:
  void calculatePositions( const Lorentz5Momentum &, const LorentzPoint &, 
			   const Lorentz5Momentum &, const Lorentz5Momentum &,
			   LorentzPoint &, LorentzPoint &) const;
  /*!< Compute the positions of the new hadrons based on the clusters position.
   *
   * This method calculates the positions of the children hadrons by a
   * call to Herwig::Smearing::gaussianSmearing with width inversely 
   * proportional to the cluster mass, around the parent cluster position.
   */

  Ptr<HadronSelector>::pointer _hadronsSelector;
  //!< Pointer to a Herwig::HadronSelector for choosing decay types
  Ptr<GlobalParameters>::pointer _globalParameters;
  //!< Pointer to a Herwig::GlobalParameters for various global data
  
  int _ClDir1;
  //!< Whether non-b clusters decay along the perturbative parton direction.
  int _ClDir2;
  //!< Whether b clusters decay along the perturbative parton direction.
  double _ClSmr1;
  //!< The S parameter from decayIntoTwoHadrons for non-b clusters.
  double _ClSmr2;
  //!< The S parameter from decayIntoTwoHadrons for b clusters.

};


}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of ClusterDecayer.
template <>
struct BaseClassTrait<Herwig::ClusterDecayer,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ClusterDecayer>: public ClassTraitsBase<Herwig::ClusterDecayer> {
  static string className() { return "/Herwig++/ClusterDecayer"; }
  //!< Return the class name.
  static string library() { return "libHwHadronization.so"; }
  /*!< Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
};

}

#endif // DOXYGEN

#include "ClusterDecayer.icc"

#endif /* HERWIG_ClusterDecayer_H */
