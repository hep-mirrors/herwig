// -*- C++ -*-
#ifndef HERWIG_LightClusterDecayer_H
#define HERWIG_LightClusterDecayer_H
/*! \class Herwig::LightClusterDecayer LightClusterDecayer.h "Herwig++\Hadronization\LightClusterDecayer.h"
 * \brief This class performs the decay of light clusters into a single hadron.
 * \author Philip Stephens
 * \author Alberto Ribon
 * \ingroup Hadronization
 *
 * This is the class that performs the decay of light clusters into 
 * only one hadron. The major difficulty is that a kinematical reshuffling
 * is necessary, between the cluster under consideration and its
 * "neighbouring" clusters, to conserve energy-momentum in one-body decay.
 * Notice that, differently from what happens in Fortran Herwig,
 * light (that is below the threshold for the production of the lightest 
 * pair of hadrons with the proper flavours) fission products, produced 
 * by the fission of heavy clusters in ClusterFissioner 
 * have been already "decayed" into single hadron (the lightest one 
 * with proper flavour) by the same latter class, without require 
 * any reshuffling. Therefore the light clusters that are treated in 
 * this LightClusterDecayer class are produced directly 
 * (originally) by the ClusterFinder. 
 *	
 * Notice:
 * - The choice of the candidate cluster with whom to reshuffle momentum 
 *   is based on the minimal space-time distance from the light cluster.
 * - An alternate choice of what is considered a "neighbour" could be
 *   implemented but was not considered for Herwig++.
 *
 * See also:
 * HadronSelector.h.
 */ 

#include <ThePEG/Handlers/HandlerBase.h>
#include <ThePEG/EventRecord/Step.h>
#include "CluHadConfig.h"
#include "HadronSelector.h"


namespace Herwig {


using namespace ThePEG;

class LightClusterDecayer: public ThePEG::HandlerBase {

public:

  inline LightClusterDecayer();
  inline LightClusterDecayer(const LightClusterDecayer &);
  virtual ~LightClusterDecayer();

  bool decay(const StepPtr &);
  /*!< This method does the decay of light hadron in one hadron.
   *
   * This method requires a kinematical reshuffling for energy-momentum 
   *conservation. This is done explicitly by the (private) method 
   * reshuffling().
   */

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  //!< Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  //!< Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  //!< Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<LightClusterDecayer> initLightClusterDecayer;
  //!< Describe a concrete class with persistent data.

  LightClusterDecayer & operator=(const LightClusterDecayer &);
  //!< Private and non-existent assignment operator.

  bool reshuffling( const long, tClusterPtr, tClusterPtr,
		    const StepPtr &, tClusterVector &) 
    throw (Veto, Stop, Exception); 
  /*!< This (private) method, called by decay(), takes care of the kinematical
   * reshuffling necessary for energy-momentum conservation.
   */
  
  Ptr<HadronSelector>::pointer _hadronsSelector;
  //!< A pointer to a Herwig::HadronSelector object used for producing hadrons.
  double _B1Lim;
  /*!< A parameter used for determining when b clusters are too light.
   *
   * This parameter is used for setting the lower threshold, \f$ t \f$ as
   * \f[ t' = t(1 + r B^1_{\rm lim}) \f]
   * where \f$ r \f$ is a random number [0,1].
   */

};

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of LightClusterDecayer.
template <>
struct BaseClassTrait<Herwig::LightClusterDecayer,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::LightClusterDecayer>: public ClassTraitsBase<Herwig::LightClusterDecayer> {
  static string className() { return "/Herwig++/LightClusterDecayer"; }
  // Return the class name.
  static string library() { return "libHwHadronization.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#endif // DOXYGEN STUFF

#include "LightClusterDecayer.icc"

#endif /* HERWIG_LightClusterDecayer_H */
