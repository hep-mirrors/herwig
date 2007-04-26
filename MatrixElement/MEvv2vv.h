// -*- C++ -*-
#ifndef HERWIG_MEvv2vv_H
#define HERWIG_MEvv2vv_H
//
// This is the declaration of the MEvv2vv class.
//

#include "Herwig++/MatrixElement/GeneralHardME.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/Vertex/Vector/VVVVertex.h"
#include "Herwig++/Helicity/Vertex/Tensor/VVTVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/VVSVertex.h"
#include "Herwig++/Helicity/Vertex/Vector/VVVVVertex.h"
#include "Herwig++/Helicity/Correlations/ProductionMatrixElement.h"
#include "MEvv2vv.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::VectorWaveFunction;
using Herwig::Helicity::VVSVertexPtr;
using Herwig::Helicity::VVVVertexPtr;
using Herwig::Helicity::VVTVertexPtr;
using Herwig::Helicity::VVVVVertexPtr;
using Herwig::Helicity::ProductionMatrixElement;

/**
 * This is the implementation of the matrix element for 
 * $2\ra 2$ massless vector-boson pair to vector-boson pair. It inherits from
 * GeneralHardME and implements the appropriate virtual member functions.
 *
 * @see \ref MEvv2vvInterfaces "The interfaces"
 * defined for MEvv2vv.
 */
class MEvv2vv: public GeneralHardME {

public:

  typedef vector<VectorWaveFunction> VBVector;

public:

  /**
   * The default constructor.
   */
  inline MEvv2vv();

public:

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

  /**
   * Compute the matrix element for \f$V\,V\ra V\,V\f$
   * @param vin VectorWaveFunctions for first incoming particle
   * @param vin2 VectorWaveFunctions for second incoming particle
   * @param vout VectorWaveFunctions for first outgoing particle
   * @param vout2  VectorWaveFunctions for outgoing particle
   * @param me2 colour averaged, spin summed ME
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement 
  vv2vvHeME(VBVector & vin1, VBVector & vin2, 
	    VBVector & vout1, VBVector & vout2,
	    double & me2) const;

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
  static ClassDescription<MEvv2vv> initMEvv2vv;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEvv2vv & operator=(const MEvv2vv &);

private:

  /**
   * Store the dynamically casted VVSVertex pointers
   */
  vector<pair<VVSVertexPtr, VVSVertexPtr> > theScaV;

  /**
   * Store the dynamically casted VVVVertex pointers
   */
  vector<pair<VVVVertexPtr, VVVVertexPtr> > theVecV;

  /**
   * Store the dynamically casted VVTVertex pointers
   */
  vector<pair<VVTVertexPtr, VVTVertexPtr> > theTenV;

  /**
   * Store the dynamically casted VVVVVertex pointer
   */
  VVVVVertexPtr theFPVertex;
  
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEvv2vv. */
template <>
struct BaseClassTrait<Herwig::MEvv2vv,1> {
  /** Typedef of the first base class of MEvv2vv. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEvv2vv class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEvv2vv>
  : public ClassTraitsBase<Herwig::MEvv2vv> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MEvv2vv"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEvv2vv is implemented. It may also include several, space-separated,
   * libraries if the class MEvv2vv depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwGeneralME.so"; }
};

/** @endcond */

}

#include "MEvv2vv.icc"

#endif /* HERWIG_MEvv2vv_H */
