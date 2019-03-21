// -*- C++ -*-
#ifndef HERWIG_NMSSMHSFSFVertex_H
#define HERWIG_NMSSMHSFSFVertex_H
//
// This is the declaration of the NMSSMHSFSFVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/SSSVertex.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Models/Susy/MixingMatrix.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Helicity
 * This class defines the coupling of Higgs bosons to the sfermions
 * in the NMSSM.
 *
 * @see \ref NMSSMHSFSFVertexInterfaces "The interfaces"
 * defined for NMSSMHSFSFVertex.
 */
class NMSSMHSFSFVertex: public Helicity::SSSVertex {

public:

  /**
   * The default constructor.
   */
  NMSSMHSFSFVertex();

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NMSSMHSFSFVertex & operator=(const NMSSMHSFSFVertex &) = delete;

private:

  /**
   * Return the coupling of the charged higgs to the sfermions
   */
  Complex chargedHiggs(Energy2 q2, long id1, long id2);

private:
  
  /** @name Stored parameters for fast access.*/
  //@{
  /**
   * A pointer to the Standard Model object 
   */
  tcHwSMPtr _theSM;

  /**
   * The CP-even higgs mixing matrix
   */
  MixingMatrixPtr _mixS;

  /**
   * The CP-odd higgs mixing matrix
   */
  MixingMatrixPtr _mixP;

  /**
   * The \f$ \tilde{t}\f$ mixing matrix 
   */
  MixingMatrixPtr _mixTp;

  /**
   * The \f$ \tilde{b}\f$ mixing matrix 
   */
  MixingMatrixPtr _mixBt;

  /**
   * The \f$ \tilde{\tau}\f$ mixing matrix 
   */
  MixingMatrixPtr _mixTa;

  /**
   * The top quark trilinear coupling
   */
  complex<Energy> _triTp;

  /**
   * The bottom quark trilinear coupling
   */
  complex<Energy> _triBt;

  /**
   * The tau lepton trilinear coupling
   */
  complex<Energy> _triTa;

  /**
   * The value of \f$\lambda\f$.
   */
  double _lambda;
  
  /**
   * The value of \f$ \lambda <S> \f$, the 
   * V.E.V of the extra gauge singlet scaled
   * by \f$\lambda\f$
   */
  Energy _lambdaVEV;

  /**
   * The value of the V.E.V \f$ v_1 \f$ 
   */
  Energy _v1;

  /**
   * The value of the V.E.V \f$ v_2 \f$ 
   */
  Energy _v2;
  
  /**
   * The value of \f$ \sin\theta_W \f$
   */
  double _sw;

  /**
   * The value of \f$ \cos\theta_W \f$
   */
  double _cw;
  
  /**
   * The value of \f$ M_W \f$
   */
  Energy _mw;

  /**
   * The value of \f$ M_Z \f$
   */
  Energy _mz;

  /**
   * The value of \f$ \sin\beta \f$
   */
  double _sb;

  /**
   * The value of \f$ \cos\beta \f$
   */
  double _cb;

  /**
   * The value of \f$ \tan\beta \f$
   */
  double _tb;

  //@}
  
  /** @name Store previously calculated values for speed. */
  //@{
  /**
   * The scale at which the last calculation took place. 
   */
  Energy2 _q2last;
  
  /**
   * The value of the dimensionless coupling \f$g_W\f$ when
   * last calculated.
   */
  double _couplast;

  /**
   * The value of mass of the counterpart SM fermion when
   * last calculated.
   */
  pair<Energy, Energy> _masslast;

  /**
   * The PDG codes of the particles in the vertex when it was last evaluated
   */
  pair<long, long> _idlast;
  //@}
};
}


#endif /* HERWIG_NMSSMHSFSFVertex_H */
