// -*- C++ -*-
#ifndef HERWIG_NMSSM_H
#define HERWIG_NMSSM_H
//
// This is the declaration of the NMSSM class.
//

#include "Herwig++/Models/Susy/MSSM.h"
#include "NMSSM.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * Here is the documentation of the NMSSM class.
 *
 * @see \ref NMSSMInterfaces "The interfaces"
 * defined for NMSSM.
 */
class NMSSM: public MSSM {

public:

  /**
   * The default constructor.
   */
  inline NMSSM();

public:

  /**
   * Mixing matrix for the neutral CP-odd Higgs bosons
   */
  inline const MixingMatrixPtr & CPoddHiggsMix() const;

  /**
   *  The NMSSM couplings
   */
  //@{
  /**
   *  Superpotential \f$\lambda\f$ term
   */
  inline double lambda() const;

  /**
   *  Superpotential \f$\kappa\f$ coupling
   */
  inline double kappa() const;

  /**
   *  The V.E.V of the extra singlet field scaled
   * by \f$ lambda\f$, 
   */
  inline Energy lambdaVEV() const;
  
  /**
   * Soft trilinear \f$S\H_2 H_1\f$ coupling
   */
  inline Energy trilinearLambda() const;

  /**
   * Soft cubic \f$S\f$ coupling
   */
  inline Energy trilinearKappa() const;
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

public:

  /**
   * Pointer to the fermion-fermion-Higgs vertex
   */
  virtual inline tFFSVertexPtr vertexFFH() const;

  /**
   * Pointer to the two electroweak gauge boson Higgs vertex.
   */
  virtual inline tVVSVertexPtr vertexWWH() const;

  /**
   * Pointer to the electroweak gauge boson Higgs-Higgs vertex.
   */
  virtual inline tVSSVertexPtr vertexWHH() const;

  /**
   * Pointer to the higgs coupling to a pair of gauginos
   */
  virtual inline tFFSVertexPtr vertexGOGOH() const;

  /**
   * Pointer to the triple higgs vertex
   */
  virtual inline tSSSVertexPtr vertexHHH() const;

  /**
   * Pointer to higgs-sfermion-sfermion vertex 
   */
  virtual inline tSSSVertexPtr vertexHSS() const;
  
protected:

  /**
   *  Extract the parameters from the input blocks
   */
  virtual void extractParameters(bool checkModel=true);

  /**
   *  Create the mixing matrices for the model
   */
  virtual void createMixingMatrices();

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
  static ClassDescription<NMSSM> initNMSSM;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NMSSM & operator=(const NMSSM &);

private:

  /**
   *  Higgs mixing matrix
   */
  MixingMatrixPtr theHiggsAMix;

  /**
   *  The NMSSM couplings
   */
  //@{
  /**
   *  Superpotential \f$\lambda\f$ term
   */
  double _lambda;

  /**
   *  Superpotential \f$\kappa\f$ coupling
   */
  double _kappa;

  /**
   * Soft trilinear \f$S\H_2 H_1\f$ coupling
   */
  Energy _theAlambda;

  /**
   * Soft cubic \f$S\f$ coupling
   */
  Energy _theAkappa;

  /**
   *  The V.E.V of the extra singlet field scaled
   * by \f$ lambda\f$
   */
  Energy _lambdaVEV;
  //@}

  /** @name The NMSSM vertices.*/
  //@{
  /**
   * The fermion-fermion higgs vertex.
   */
  FFSVertexPtr _ffhvertex;

  /**
   * The vector-vector-higgs vertex.
   */
  VVSVertexPtr _wwhvertex;
  
  /**
   * The vector-higgs-higgs vertex
   */
  VSSVertexPtr _whhvertex;

  /**
   * The coupling of a pair of gauginos to the higgs 
   */
  FFSVertexPtr _gogohvertex;

  /**
   * The triple higgs coupling 
   */
  SSSVertexPtr _hhhvertex;

  /**
   * The higgs sfermion vertex 
   */
  SSSVertexPtr _hssvertex;

  /**
   * The neutralino-sfermion-fermion vertex 
   */
  FFSVertexPtr _nfsvertex;

  /**
   * The neutralino-neutralino-Z vertex 
   */
  FFVVertexPtr _nnzvertex;

  /**
   * The chargino-neutralino-W vertex 
   */
  FFVVertexPtr _cnwvertex;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NMSSM. */
template <>
struct BaseClassTrait<Herwig::NMSSM,1> {
  /** Typedef of the first base class of NMSSM. */
  typedef Herwig::MSSM NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NMSSM class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NMSSM>
  : public ClassTraitsBase<Herwig::NMSSM> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::NMSSM"; }
  /**
   * The name of a file containing the dynamic library where the class
   * NMSSM is implemented. It may also include several, space-separated,
   * libraries if the class NMSSM depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so HwNMSSM.so"; }
};

/** @endcond */

}

#include "NMSSM.icc"

#endif /* HERWIG_NMSSM_H */
