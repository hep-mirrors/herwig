// -*- C++ -*-
#ifndef HERWIG_ColourReconnector_H
#define HERWIG_ColourReconnector_H

#include <ThePEG/Handlers/HandlerBase.h>
#include "CluHadConfig.h"


namespace Herwig {

using namespace ThePEG;

/** \ingroup Hadronization
 *  \class ColourReconnector
 *  \brief Class for changing colour reconnections of partons.
 *  \author Alberto Ribon
 * 
 *  This class does the nonperturbative colour rearrangement, after the 
 *  nonperturbative gluon splitting and the "normal" cluster formation. 
 *  It uses the list of particles in the event record, and the collections of
 *  "usual" clusters which is passed to the main method. If the colour 
 *  reconnection is actually accepted, then the previous collections of "usual"
 *  clusters is first deleted and then the new one is created.
 *
 *  Note: by default this class does nothing. It can be inherited and overridden
 *  in future hadronization models.
 */
class ThePEG::PartialCollisionHandler; // forward declaration


class ColourReconnector: public ThePEG::HandlerBase {

public:

  inline ColourReconnector();
  inline ColourReconnector(const ColourReconnector &);
  virtual ~ColourReconnector();

  /**
   * Does the colour rearrangment.
   *
   * Does the colour rearrangement, starting from the list of particles
   * in the event record, and the collection of "usual" clusters passed
   * in input. If the actual rearrangement is accepted, the new collection 
   * of clusters is overriden to the intial one.
   */
  void rearrange(PartialCollisionHandler & ch, const StepPtr & pstep,
                 ClusterVector & clusters) throw(Veto, Stop, Exception);
    
public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();

  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ColourReconnector> initColourReconnector;

  /**
   * Private and non-existent assignment operator.
   */
  ColourReconnector & operator=(const ColourReconnector &);

  /**
   * The numer of colours in the reconstruction.
   */
  int    _ClReco;

  /**
   * The probability of a reconstruction.
   */
  double _PReco;

};


}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ColourReconnector.
 */
template <>
struct BaseClassTrait<Herwig::ColourReconnector,1> {
  typedef ThePEG::HandlerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ColourReconnector>: public ClassTraitsBase<Herwig::ColourReconnector> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/ColourReconnector"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwHadronization.so"; }

};

}

#endif // DOXYGEN

#include "ColourReconnector.icc"

#endif /* HERWIG_ColourReconnector_H */
