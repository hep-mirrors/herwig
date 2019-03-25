// -*- C++ -*-
#ifndef Herwig_RPVFFSVertex_H
#define Herwig_RPVFFSVertex_H
//
// This is the declaration of the RPVFFSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "RPV.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the RPVFFSVertex class.
 *
 * @see \ref RPVFFSVertexInterfaces "The interfaces"
 * defined for RPVFFSVertex.
 */
class RPVFFSVertex: public Helicity::FFSVertex {

public:

  /**
   * The default constructor.
   */
  RPVFFSVertex();
  
  /**
   * Calculate the coupling for the vertex
   * @param q2 The scale to at which evaluate the coupling.
   * @param particle1 The first particle in the vertex.
   * @param particle2 The second particle in the vertex.
   * @param particle3 The third particle in the vertex.
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			   tcPDPtr particle3);

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
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

protected:

  /**
   *    Neutralino-sfermion-fermion
   */
  void neutralinoSfermionCoupling(Energy2 q2, tcPDPtr fermion, tcPDPtr gaugino,
				  tcPDPtr sfermion);

  /**
   *    Chargino-sfermion-fermion
   */
  void charginoSfermionCoupling(Energy2 q2, tcPDPtr fermion, tcPDPtr gaugino,
				tcPDPtr sfermion);

  /**
   *  Higgs to SM fermions
   */
  void higgsFermionCoupling(Energy2 q2, tcPDPtr f1, tcPDPtr f2, tcPDPtr higgs);

  /**
   *  Higgs to gauginos (general slepton case with mixing)
   */
  void higgsGauginoCoupling(Energy2 q2, tcPDPtr f1, tcPDPtr f2, tcPDPtr higgs);

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RPVFFSVertex & operator=(const RPVFFSVertex &) = delete;

private:

  /**
   *  Which interactions to include
   */
  unsigned int interactions_;

  /**
   *  Mixing Matrices
   */
  //@{  
  /**
   * Pointer to stop mixing matrix
   */
  tMixingMatrixPtr _stop;

  /**
   * Pointer to sbottom mixing matrix
   */
  tMixingMatrixPtr _sbot;

  /**
   * Pointer to stau mixing matrix 
   */
  tMixingMatrixPtr _stau;

  /**
   * Pointer to U chargino mixing matrix 
   */
  tMixingMatrixPtr _umix;

  /**
   * Pointer to V chargino mixing matrix
   */
  tMixingMatrixPtr _vmix;
  
  /**
   * Pointer to the neutralino mixing matrix
   */
  tMixingMatrixPtr _nmix;

  /**
   *  Neutral scalar Higgs mixing matrix
   */
  MixingMatrixPtr _mixH;

  /**
   *  Neutral pseudoscalar Higgs mixing matrix
   */
  MixingMatrixPtr _mixP;

  /**
   *  Charged Higgs mixing matrix
   */
  MixingMatrixPtr _mixC;
  //@}

  /**
   * The mass of the \f$W\f$.
   */
  Energy mw_;

  /**
   * The energy scale at which the coupling 
   * was last evaluated 
   */
  Energy2 _q2last;

  /**
   * The value of the coupling at the scale last evaluated
   */
  Complex _couplast;
  
  /**
   * Store the value of the left coupling when it was last evaluated
   */
  Complex _leftlast;
  
  /**
   * Store the value of the right coupling when it was last evaluated
   */
  Complex _rightlast;

  /**
   * Store the id of the last gaugino to be evaluated
   */
  long _id1last;
  
  /**
   * Store the id of the last SM fermion to be evaluated
   */
  long _id2last;

  /**
   * Store the id of the last scalar to be evaluated
   */
  long _id3last;

  /**
   *  Include Yukawa's ? in neutralino/chargino interactions
   */
  bool yukawa_;
 
  /**
   * Pointer to the Susy Model object
   */
  tRPVPtr model_;

  /**
   * \f$\sin(\theta_w)\f$
   */
  double _sw;

  /**
   * \f$\cos(\theta_w)\f$
   */
  double _cw;

  /**
   * \f$\sin(\beta)\f$
   */
  double _sb;

  /**
   * \f$\cos(\beta)\f$
   */
  double _cb;

  /**
   *  VEV for down type Higgs
   */
  Energy vd_;

  /**
   *  VEV for up   type Higgs
   */
  Energy vu_;

  /**
   *  Values of the masses
   */
  pair<Energy,Energy> _massLast;

  /**
   *  OCCHL
   */
  vector<vector<vector<Complex> > > OCCHL_;

  /**
   *  ONNHL
   */
  vector<vector<vector<Complex> > > ONNHL_;

  /**
   *  OCCAL
   */
  vector<vector<vector<Complex> > > OCCAL_;

  /**
   *  ONNAL
   */
  vector<vector<vector<Complex> > > ONNAL_;

  /**
   *  OCNSL
   */
  vector<vector<vector<Complex> > > OCNSL_;

  /**
   *   OCNSR
   */
  vector<vector<vector<Complex> > > OCNSR_;
};

}

#endif /* Herwig_RPVFFSVertex_H */
