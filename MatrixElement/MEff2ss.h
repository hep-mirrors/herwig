// -*- C++ -*-
#ifndef HERWIG_MEff2ss_H
#define HERWIG_MEff2ss_H
//
// This is the declaration of the MEff2ss class.
//

#include "GeneralHardME.h"
#include "Herwig++/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/VSSVertex.h"
#include "Herwig++/Helicity/Vertex/Tensor/FFTVertex.h"
#include "Herwig++/Helicity/Vertex/Tensor/SSTVertex.h"
#include "MEff2ss.fh"

namespace Herwig {
using namespace ThePEG;
using Herwig::Helicity::FFSVertexPtr;
using Herwig::Helicity::FFVVertexPtr;
using Herwig::Helicity::VSSVertex;
using Herwig::Helicity::FFTVertexPtr;
using Herwig::Helicity::SSTVertexPtr;

/**
 * The MEff2ss class is designed to implement the matrix element for a
 * fermion-antifermion to scalar-scalar hard process. It inherits from 
 * GeneralHardME and implements the appropriate virtual functions for this 
 * specific spin combination.
 *
 * @see \ref MEff2ssInterfaces "The interfaces"
 * defined for MEff2ss.
 * @see GeneralHardME
 */
class MEff2ss: public GeneralHardME {
  
public: 
  
  /**
   * Convenient typedef for VSSVertex pointer
   */
  typedef Ptr<VSSVertex>::pointer VSSVertexPtr;

public:

  /**
   * The default constructor.
   */
  inline MEff2ss();

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
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


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  //@}

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEff2ss> initMEff2ss;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEff2ss & operator=(const MEff2ss &);

  private:

  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * fermion
   */
  vector<pair<FFSVertexPtr, FFSVertexPtr> > theFerm;

  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * vector
   */
  vector<pair<FFVVertexPtr, VSSVertexPtr> > theVec;
  
  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * tensor
   */
  vector<pair<FFTVertexPtr, SSTVertexPtr> > theTen;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/// \if TRAITSPECIALIZATIONS

/** This template specialization informs ThePEG about the
 *  base classes of MEff2ss. */
template <>
struct BaseClassTrait<Herwig::MEff2ss,1> {
  /** Typedef of the first base class of MEff2ss. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEff2ss class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEff2ss>
  : public ClassTraitsBase<Herwig::MEff2ss> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MEff2ss"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEff2ss is implemented. It may also include several, space-separated,
   * libraries if the class MEff2ss depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libHwGeneralME.so"; }
};

/// \endif

}

#include "MEff2ss.icc"

#endif /* HERWIG_MEff2ss_H */
