// -*- C++ -*-
#ifndef HERWIG_NMSSMHHHVertex_H
#define HERWIG_NMSSMHHHVertex_H
//
// This is the declaration of the NMSSMHHHVertex class.
//
#include "ThePEG/Helicity/Vertex/Scalar/SSSVertex.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Models/Susy/MixingMatrix.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/** \ingroup Helicity
 * The NMSSMHHHVertex defines the triple
 * Higgs coupling in the NMSSM.
 *
 * @see \ref NMSSMHHHVertexInterfaces "The interfaces"
 * defined for NMSSMHHHVertex.
 */
class NMSSMHHHVertex: public SSSVertex {

public:

  /**
   * The default constructor.
   */
  NMSSMHHHVertex();

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
   * Calculate the couplings. This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
			   tcPDPtr part3);

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
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<NMSSMHHHVertex> initNMSSMHHHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NMSSMHHHVertex & operator=(const NMSSMHHHVertex &) = delete;

private:

  /**
   * The mixing matrix combination \f$U^S_{ai}U^S_{bj}U^S_{ck}\f$
   * @param a The row element of the first CP-even mixing matrix   
   * @param b The column element of the first CP-even mixing matrix
   * @param c The row element of the second CP-even mixing matrix   
   * @param i The column element of the second CP-even mixing matrix
   * @param j The row element of the third CP-even mixing matrix   
   * @param k The column element of the third CP-even mixing matrix
   */
  Complex usMix(unsigned int a, unsigned int b, unsigned int c,
		unsigned int i, unsigned int j, unsigned int k) const {
    return (*_mixS)(a,i)*(*_mixS)(b,j)*(*_mixS)(c,k) +
           (*_mixS)(a,i)*(*_mixS)(c,j)*(*_mixS)(b,k) +
   	   (*_mixS)(b,i)*(*_mixS)(a,j)*(*_mixS)(c,k) +
   	   (*_mixS)(b,i)*(*_mixS)(c,j)*(*_mixS)(a,k) +
           (*_mixS)(c,i)*(*_mixS)(a,j)*(*_mixS)(b,k) +
           (*_mixS)(c,i)*(*_mixS)(b,j)*(*_mixS)(a,k);
  }
  
  /**
   * The mixing matrix combination \f$U^S_{ai}U^P_{bj}U^P_{ck}\f$
   * @param a The row element of the first CP-even mixing matrix   
   * @param b The column element of the first CP-even mixing matrix
   * @param c The row element of the second CP-even mixing matrix   
   * @param i The column element of the second CP-even mixing matrix
   * @param j The row element of the third CP-even mixing matrix   
   * @param k The column element of the third CP-even mixing matrix
   */
  Complex upMix(unsigned int a, unsigned int b, unsigned int c,
		unsigned int i, unsigned int j, unsigned int k) const {
    return (*_mixS)(a,i)*((*_mixP)(b,j)*(*_mixP)(c,k) + 
			  (*_mixP)(c,j)*(*_mixP)(b,k));
  }

private:

  /**
   * A pointer to the object containing the SM parameters 
   */
  tcHwSMPtr _theSM;

  /**
   * The \f$W\f$ mass 
   */
  Energy _mw;
  
  /**
   * The \f$Z\f$ mass 
   */
  Energy _mz;
     /**
   * The \f$b\f$ mass 
   */
  Energy _mb;
  
  /**
   * The \f$t\f$ mass 
   */
  Energy _mt;
  
  /**
   * \f$\sin^2\theta_W\f$
   */
  double _sw2;

  /**
   * \f$\cos\theta_W\f$
   */
  double _cw;

  /**
   * The CP-even Higgs mixing matrix 
   */
  MixingMatrixPtr _mixS;

  /**
   * The CP-odd Higgs mixing matrix 
   */
  MixingMatrixPtr _mixP;

  /**
   *  The coefficient of the trilinear \f$SH_2 H_1\f$ term in the superpotential
   */
  double _lambda;

  /**
   *  The coefficient of the cubic singlet term in the superpotential
   */
  double _kappa;
  
  /**
   * The product \f$\lambda \langle S\rangle \f$.
   */
  Energy _lambdaVEV;
  
  /**
   * The soft trilinear \f$SH_2 H_1\f$ coupling
   */
  Energy _theAl;
  
  /**
   * The soft cubic \f$S\f$ coupling
   */
  Energy _theAk;

  /**
   * \f$\sin\beta\f$
   */
  double _sb;

  /**
   * \f$\cos\beta\f$
   */
  double _cb;

  /**
   * \f$\sin2\beta\f$
   */
  double _s2b;

  /**
   * \f$\cos2\beta\f$
   */
  double _c2b;

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
   * The scale at which this vertex was last evaluated 
   */
  Energy2 _q2last;

  /**
   * The value of the EW coupling when it was last evaluated
   */
  double _glast;
        /**
   * left 3rd generation scalar quark mass
   */
  Energy _MQ3;

  /**
   *  right scalar top mass
   */
  Energy _MU2;

  /**
   *  Whether or onto to include the radiative terms
   */
  bool _includeRadiative;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NMSSMHHHVertex. */
template <>
struct BaseClassTrait<Herwig::NMSSMHHHVertex,1> {
  /** Typedef of the first base class of NMSSMHHHVertex. */
  typedef Helicity::SSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NMSSMHHHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NMSSMHHHVertex>
  : public ClassTraitsBase<Herwig::NMSSMHHHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::NMSSMHHHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * NMSSMHHHVertex is implemented. It may also include several, space-separated,
   * libraries if the class NMSSMHHHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so HwNMSSM.so"; }
};

/** @endcond */

}

#endif /* HERWIG_NMSSMHHHVertex_H */
