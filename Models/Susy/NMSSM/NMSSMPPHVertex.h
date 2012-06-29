// -*- C++ -*-
#ifndef HERWIG_NMSSMPPHVertex_H
#define HERWIG_NMSSMPPHVertex_H
//
// This is the declaration of the NMSSMPPHVertex class.
//

#include "Herwig++/Models/General/VVSLoopVertex.h"
#include "Herwig++/Models/StandardModel/StandardModel.fh"
#include "Herwig++/Models/Susy/MixingMatrix.h"

namespace Herwig {
using namespace ThePEG;

/**
 * This class implements the effective vertex for a higgs coupling
 * to a pair of photons in the NMSSM.
 *
 * @see \ref NMSSMPPHVertexInterfaces "The interfaces"
 * defined for NMSSMPPHVertex.
 */
class NMSSMPPHVertex: public VVSLoopVertex {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  NMSSMPPHVertex();
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
   * Calculate couplings
   *@param q2 Scale at which to evaluate coupling
   *@param p1 ParticleData pointer to first particle
   *@param p2 ParticleData pointer to second particle
   *@param p3 ParticleData pointer to third particle
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr p1, tcPDPtr p2,
			   tcPDPtr p3);

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<NMSSMPPHVertex> initNMSSMPPHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NMSSMPPHVertex & operator=(const NMSSMPPHVertex &);

private:

  /** @name Stored parameters. */
  //@{
  /**
   * The SM pointer 
   */
  tcHwSMPtr _theSM;
  
  /**
   * \f$ \sin\theta_W\f$ 
   */
  double _sw;

  /**
   * \f$ \cos\theta_W\f$ 
   */
  double _cw;

  /**
   * \f$ M_W\f$
   */
  Energy _mw;

  /**
   * \f$ M_Z \f$
   */
  Energy _mz;
  
   /**
   * The product \f$\lambda \langle S\rangle \f$.
   */
  
  Energy _lambdaVEV;
  
  /**
   *  The coefficient of the trilinear \f$SH_2 H_1\f$ term in the superpotential
   */
  
  double _lambda;
  
  /**
   * The value of the VEV of the higgs that couples to the up-type sector
   *  \f$ g*sqrt(2)M_W\cos\beta \f$
   */
  
  Energy _v1;

  /**
   * The value of the VEV of the higgs that couples to the down-type sector
   *  \f$ g*sqrt(2)M_W\cos\beta \f$
   */ 
  Energy _v2;
  
  /**
   * The top quark trilinear coupling
   */
  complex<Energy> _triTp;

  /**
   * The bottom quark trilinear coupling
   */
  complex<Energy> _triBt;

  /**
   * The tau quark trilinear coupling
   */
  complex<Energy> _triTa;
  
  /**
   * A pointer to the top quark
   */
  tcPDPtr _top;

  /**
   * A pointer to the bottom quark
   */
  tcPDPtr _bt;

  /**
   * A pointer to the tau lepton
   */
  tcPDPtr _tau;
  
  /**
   * CP-even Higgs mixing matrix 
   */
  MixingMatrixPtr _mixS;
  
  /**
   * CP-even Higgs mixing matrix 
   */
  MixingMatrixPtr _mixP;
  
  /**
   * \f$\tilde{t}\f$ mixing matrix  
   */
  MixingMatrixPtr _mixQt;
  
  /**
   * \f$\tilde{b}\f$ mixing matrix  
   */
  MixingMatrixPtr _mixQb;
  
  /**
   * \f$\tilde{\tau}\f$ mixing matrix  
   */
  MixingMatrixPtr _mixLt;

  /**
   * \f$ \sin\beta\f$ 
   */
  double _sb;

  /**
   * \f$ \cos\beta\f$ 
   */
  double _cb;

  /**
   *  The coefficient of the cubic singlet term in the superpotential
   */
  double _kappa;

  /**
   * The value of the VEV of the higgs that couples to the down-type sector
   *  \f$ g*sqrt(2)M_W\cos\beta \f$
   */
  Energy _vu;

  /**
   * The value of the VEV of the higgs that couples to the up-type sector
   * i.e. \f$ g*sqrt(2)M_W\sin\beta \f$
   */
  Energy _vd;
  
    /**
   * The value of the VEV of the singlet higgs
   */
  Energy _s;
  
  /**
   * The soft trilinear \f$SH_2 H_1\f$ coupling
   */
  Energy _theAl;

  /**
   * The U mixing matrix
   */
  tMixingMatrixPtr _mixU;

  /**
   * The V mixing matrix
   */
  tMixingMatrixPtr _mixV;
  
  /**
   * The top and bottom quark masses calculated at the last value
   * of \f$q^2\f$ 
   */
  pair<Energy, Energy> _masslast;

  /**
   * The scale at which the coupling was last evaluated 
   */
  Energy2 _q2last;

  /**
   * The value of the overall normalisation when the coupling was last 
   * evaluated.
   */
  double _couplast;

  /**
   * The value of the weak coupling was last 
   * evaluated.
   */
  double _coup;
  
  /**
   * The PDG code of the Higgs particle when the vertex was last evaluated 
   */
  long _hlast;

  /**
   * Whether the tensor coefficient need recalculating or not 
   */
  bool _recalc;
  //@}
};
  
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NMSSMPPHVertex. */
template <>
struct BaseClassTrait<Herwig::NMSSMPPHVertex,1> {
  /** Typedef of the first base class of NMSSMPPHVertex. */
  typedef Herwig::VVSLoopVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NMSSMPPHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NMSSMPPHVertex>
  : public ClassTraitsBase<Herwig::NMSSMPPHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::NMSSMPPHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * NMSSMPPHVertex is implemented. It may also include several, space-separated,
   * libraries if the class NMSSMPPHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so HwNMSSM.so"; }
};

/** @endcond */

}

#endif /* HERWIG_NMSSMPPHVertex_H */
