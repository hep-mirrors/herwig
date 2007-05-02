// -*- C++ -*-
#ifndef HERWIG_MEff2vv_H
#define HERWIG_MEff2vv_H
//
// This is the declaration of the MEff2vv class.
//

#include "GeneralHardME.h"
#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig++/Helicity/Vertex/Tensor/FFTVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/GeneralSVVVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/VVSVertex.h"
#include "Herwig++/Helicity/Vertex/Tensor/VVTVertex.h"
#include "Herwig++/Helicity/Vertex/Vector/VVVVertex.h"
#include "MEff2vv.fh"

namespace Herwig {
using namespace ThePEG;
using Herwig::Helicity::FFVVertexPtr;
using Herwig::Helicity::FFTVertexPtr;
using Herwig::Helicity::FFSVertexPtr;
using Herwig::Helicity::VVSVertexPtr;
using Herwig::Helicity::GeneralSVVVertexPtr;
using Herwig::Helicity::VVTVertexPtr;
using Herwig::Helicity::VVVVertexPtr;

/**
 * This class implements the matrix element calculation for a generic
 * \f$\Psi \Psi \rightarrow V^{\mu} V^{\nu}\f$ process. 
 *
 * @see \ref MEff2vvInterfaces "The interfaces"
 * defined for MEff2vv.
 */
class MEff2vv: public GeneralHardME {

public:

  /**
   * The default constructor.
   */
  inline MEff2vv();

  /** @name Virtual functions required by the GeneralHardME class. */
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
  inline virtual void doinit() throw(InitException);
  //@}
  
private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEff2vv> initMEff2vv;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEff2vv & operator=(const MEff2vv &);

private:

  /**
   * Storage for a dynamically cast vertices for a tchannel vector
   * intermediate
   */
  vector<pair<FFVVertexPtr, FFVVertexPtr> > theFerm;

  /**
   * Storage for a dynamically cast vertices for a schannel vector
   * intermediate
   */
  vector<pair<FFVVertexPtr, VVVVertexPtr> > theVec;

  /**
   * Storage for a dynamically cast vertices for a schannel scalar
   * intermediate
   */
  vector<pair<FFTVertexPtr, VVTVertexPtr> > theTen;

  /**
   * Storage for a dynamically cast vertices for a schannel scalar
   * intermediate for massless external vector bosons
   */
  vector<pair<FFSVertexPtr, GeneralSVVVertexPtr> > theSca1;

  /**
   * Storage for a dynamically cast vertices for a schannel scalar
   * intermediate for massive external vector bosons
   */
  vector<pair<FFSVertexPtr, VVSVertexPtr> > theSca2;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEff2vv. */
template <>
struct BaseClassTrait<Herwig::MEff2vv,1> {
  /** Typedef of the first base class of MEff2vv. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEff2vv class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEff2vv>
  : public ClassTraitsBase<Herwig::MEff2vv> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MEff2vv"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEff2vv is implemented. It may also include several, space-separated,
   * libraries if the class MEff2vv depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libHwGeneralME.so"; }
};

/** @endcond */

}

#include "MEff2vv.icc"

#endif /* HERWIG_MEff2vv_H */
