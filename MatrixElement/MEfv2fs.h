// -*- C++ -*-
#ifndef HERWIG_MEfv2fs_H
#define HERWIG_MEfv2fs_H
//
// This is the declaration of the MEfv2fs class.
//

#include "GeneralHardME.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/Correlations/ProductionMatrixElement.h"
#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.fh"
#include "Herwig++/Helicity/Vertex/Scalar/FFSVertex.fh"
#include "Herwig++/Helicity/Vertex/Scalar/VSSVertex.fh"
#include "MEfv2fs.fh"

namespace Herwig {
using namespace ThePEG;
using Herwig::Helicity::FFVVertexPtr;
using Herwig::Helicity::FFSVertexPtr;
using Herwig::Helicity::VSSVertexPtr;
using Herwig::Helicity::SpinorWaveFunction;
using Herwig::Helicity::SpinorBarWaveFunction;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::ProductionMatrixElement;

/**
 * This class is designed to implement the matrix element for 
 * fermion-vector to fermion scalar. It inherits from GeneralHardME 
 * and implements the required virtual functions.
 *
 * @see @see \ref MEfv2fsInterfaces "The Interfaces"
 * defined for MEfv2fs.
 * @see GeneralHardME
 */
class MEfv2fs: public GeneralHardME {

  /** Vector of SpinorWaveFunctions. */
  typedef vector<SpinorWaveFunction> SpinorVector;

  /** Vector of SpinorBarWaveFunctions. */
  typedef vector<SpinorBarWaveFunction> SpinorBarVector;

  /** Vector of VectorWaveFunctions. */
  typedef vector<VectorWaveFunction> VecWFVector;

public:

  /**
   * The default constructor.
   */
  inline MEfv2fs();

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
   * Construct the vertex information for the spin correlations
   * @param subp Pointer to the relevent SubProcess
   */
  virtual void constructVertex(tSubProPtr subp);

private:

  /** @name Functions to calculate production matrix elements and me2. */
  //@{
  /**
   * Calculate me2 and the production matrix element for the normal mode.
   * @param spIn Vector of SpinorWaveFunction for the incoming fermion
   * @param vecIn Vector of VectorWaveFunction for incoming boson
   * @param spbOut Vector of SpinorBarWaveFunction for outgoing fermion
   * @param scaOut ScalarWaveFunction for outgoing scalar.
   * @param full_me The value of me2 calculation
   */
  ProductionMatrixElement fv2fbsHeME(const SpinorVector & spIn, 
				     const VecWFVector & vecIn,
				     const SpinorBarVector & spbOut,
				     const ScalarWaveFunction & scaOut,
				     double & full_me) const;
  
  /**
   * Calculate me2 and the production matrix element for the cc mode.
   * @param spbIn Vector of SpinorBarWaveFunction for the incoming fermion
   * @param vecIn Vector of VectorWaveFunction for incoming boson
   * @param spOut Vector of SpinorWaveFunction for outgoing fermion
   * @param scaOut ScalarWaveFunction for outgoing scalar.
   * @param full_me The value of me2 calculation
   */
  ProductionMatrixElement fbv2fsHeME(const SpinorBarVector & spbIn, 
				     const VecWFVector & vecIn,
				     const SpinorVector & spOut,
				     const ScalarWaveFunction & scaOut,
				     double & full_me) const;
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
  inline void doinit() throw(InitException);
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
  static ClassDescription<MEfv2fs> initMEfv2fs;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEfv2fs & operator=(const MEfv2fs &);

private:

  /**
   * Store a pair of  FFSVertex and VSSVertex pointers  
   */
  vector<pair<FFSVertexPtr, VSSVertexPtr> > theScaV;

  /**
   * Store a pair of  FFSVertex and FFVVertex pointers  
   */
  vector<pair<FFSVertexPtr, FFVVertexPtr> > theFermV;
  
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEfv2fs. */
template <>
struct BaseClassTrait<Herwig::MEfv2fs,1> {
  /** Typedef of the first base class of MEfv2fs. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEfv2fs class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEfv2fs>
  : public ClassTraitsBase<Herwig::MEfv2fs> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MEfv2fs"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEfv2fs is implemented. It may also include several, space-separated,
   * libraries if the class MEfv2fs depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libGeneralHardME.so"; }
};

/** @endcond */

}

#include "MEfv2fs.icc"

#endif /* HERWIG_MEfv2fs_H */
