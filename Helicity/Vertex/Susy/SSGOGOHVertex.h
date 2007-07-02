// -*- C++ -*-
#ifndef HERWIG_SSGOGOHVertex_H
#define HERWIG_SSGOGOHVertex_H
//
// This is the declaration of the SSGOGOHVertex class.
//

#include "Herwig++/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Models/Susy/MSSM.h"
#include "SSGOGOHVertex.fh"

namespace Herwig {
namespace Helicity {

/**
 * The is the coupling of higgs bosons in the MSSM to a pair
 * of SM fermions.
 *
 * @see \ref SSGOGOHVertexInterfaces "The interfaces"
 * defined for SSGOGOHVertex.
 */
class SSGOGOHVertex: public FFSVertex {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SSGOGOHVertex();

  /**
   * The destructor.
   */
  virtual ~SSGOGOHVertex();
  //@}

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
  
  /**
   * Calculate the coupling for the vertex
   * @param q2 The scale to at which evaluate the coupling.
   * @param particle1 The first particle in the vertex.
   * @param particle2 The second particle in the vertex.
   * @param particle3 The third particle in the vertex.
   * @param iint The incoming particle(only needed for vertices with
   * Majorana particles)
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			   tcPDPtr particle3, int iint);
  
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


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SSGOGOHVertex> initSSGOGOHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSGOGOHVertex & operator=(const SSGOGOHVertex &);

private:
  
  /**
   * A pointer to the MSSM object.
   */
  tMSSMPtr theMSSM;

  /**
   * The mass of the \f$W\f$.
   */
  Energy theMw;

  /**
   * The matrix \f$R_{ij}\f$ 
   */
  vector<vector<Complex> > theRij;

  /**
   * The matrix \f$Q_{ij}\f$ 
   */
  vector<vector<Complex> > theQij;

  /**
   * The matrix \f$Q_{ij}^{L'}\f$ 
   */
  vector<vector<Complex> > theQijLp;

  /**
   * The matrix \f$Q_{ij}^{R'}\f$ 
   */
  vector<vector<Complex> > theQijRp;

  /**
   * The matrix \f$R_{ij}^{''}\f$ 
   */
  vector<vector<Complex> > theRijdp;

  /**
   * The matrix \f$Q_{ij}^{''}\f$ 
   */
  vector<vector<Complex> > theQijdp; 

  /**
   * The value of \f$\sin\theta_W\f$
   */
  double theSw;

  /**
   * The value of \f$\sin\alpha\f$ 
   */
  double theSa;

  /**
   * The value of \f$\sin\beta\f$ 
   */
  double theSb;

  /**
   * The value of \f$\cos\alpha\f$ 
   */
  double theCa;

  /**
   * The value of \f$\cos\beta\f$ 
   */
  double theCb;  

  /**
   * The value of \f$\cos 2\beta\f$ 
   */
  double theC2b;  
  
  /**
   * The value of \f$\sin(\beta - \alpha)\f$ 
   */
  double theSba;

  /**
   * The value of \f$\cos(\beta - \alpha)\f$ 
   */
  double theCba;

  /**
   * The value of the coupling when it was last evaluated.
   */
  Complex theCoupLast;

  /**
   * The value of the left-coupling when it was last evaluated.
   */
  Complex theLLast;
  
  /**
   * The value of the right-coupling when it was last evaluated.
   */
  Complex theRLast;

  /**
   * The ID of the last higgs for which the vertex was evaluated
   */
  long theHLast;

  /**
   * The ID of the first gaugino when the coupling was las evaluated
   */
  long theID1Last;

  /**
   * The ID of the first gaugino when the coupling was las evaluated
   */
  long theID2Last;
};

}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSGOGOHVertex. */
template <>
struct BaseClassTrait<Herwig::Helicity::SSGOGOHVertex,1> {
  /** Typedef of the first base class of SSGOGOHVertex. */
  typedef Herwig::Helicity::FFSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSGOGOHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Helicity::SSGOGOHVertex>
  : public ClassTraitsBase<Herwig::Helicity::SSGOGOHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SSGOGOHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSGOGOHVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSGOGOHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so HwSusyVertex.so"; }
};

/** @endcond */

}

#include "SSGOGOHVertex.icc"

#endif /* HERWIG_SSGOGOHVertex_H */
