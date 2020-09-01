// -*- C++ -*-
#ifndef Herwig_MEee2eeX_H
#define Herwig_MEee2eeX_H
//
// This is the declaration of the MEee2eeX class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "GammaGammaAmplitude.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEee2eeX class implements the processes \f$e^+e^-\toe^+e^- \gamma\gamma\f$ followed by
 * \f$\gamma\gamma\to X\f$ using the GammaGammaAmplitude class
 *
 * @see \ref MEee2eeXInterfaces "The interfaces"
 * defined for MEee2eeX.
 */
class MEee2eeX: public HwMEBase {

public:
  
  /**
   * The default constructor.
   */
  MEee2eeX() : Q2_1min_(ZERO), Q2_1max_(Constants::MaxEnergy2),
	       Q2_2min_(ZERO), Q2_2max_(Constants::MaxEnergy2)
  {}
  
public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const {
    return amp_->orderInAlphaS();
  }

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const {
    return amp_->orderInAlphaS()+2;
  }

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics();

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given nDim() uniform
   * random numbers in the interval \f$ ]0,1[ \f$. To help the phase space
   * generator, the dSigHatDR should be a smooth function of these
   * numbers, although this is not strictly necessary.
   * @param r a pointer to the first of nDim() consecutive random numbers.
   * @return true if the generation succeeded, otherwise false.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

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

protected:

  /**
   * Calculation of the leptonic currents
   */
  //@{
  /**
   *  Calculation of the electron current
   */
  vector<VectorWaveFunction> electronCurrent() const;

  /**
   *  Calculation of the positron current
   */
  vector<VectorWaveFunction> positronCurrent() const;
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
  virtual void doinit();
  //@}

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2eeX & operator=(const MEee2eeX &) = delete;

private:

  /**
   *   fermion-fermion-photon vertex
   */
  AbstractFFVVertexPtr FFPVertex_;

  /**
   *  Pointer to the particle data object for the photon
   */
  PDPtr gamma_;

  /**
   *  Pointer to the amplitude for the \f$\gamma\gamma$ process
   */
  GammaGammaAmpPtr amp_;

private:

  /**
   *   Cuts
   */
  //@{
  /**
   *  Minimum value of \f$Q_1^2\f$
   */
  Energy2 Q2_1min_;
  
  /**
   *  Maximum value of \f$Q_1^2\f$
   */
  Energy2 Q2_1max_;
  
  /**
   *  Minimum value of \f$Q_2^2\f$
   */
  Energy2 Q2_2min_;
  
  /**
   *  Maximum value of \f$Q_2^2\f$
   */
  Energy2 Q2_2max_;
  //@}

private:

  /**
   *  Kinematic quantities stored for accuracy
   */
  //@{
  /**
   *  Outgoing electron
   */
  double cHalf1_,sHalf1_,phi1_;
  Energy2 t1_;
  
  /**
   *  Outgoing positron
   */
  double cHalf2_,sHalf2_,phi2_;
  Energy2 t2_;
  //@}

};

}

#endif /* Herwig_MEee2eeX_H */



// // -*- C++ -*-
// #ifndef Herwig_MEee2eePseudoScalar_H
// #define Herwig_MEee2eePseudoScalar_H
// //
// // This is the declaration of the MEee2eePseudoScalar class.
// //

// #include "Herwig/MatrixElement/HwMEBase.h"
// #include "Herwig/Models/StandardModel/StandardModel.h"
// #include "Herwig/MatrixElement/ProductionMatrixElement.h"
// #include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
// #include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

// namespace Herwig {

// using namespace ThePEG;

// /**
//  * Here is the documentation of the MEee2eePseudoScalar class.
//  *
//  * @see \ref MEee2eePseudoScalarInterfaces "The interfaces"
//  * defined for MEee2eePseudoScalar.
//  */
// class MEee2eePseudoScalar: public HwMEBase {

// public:

//   /**
//    * The default constructor.
//    */
//   MEee2eePseudoScalar() : Q2_1min_(ZERO), Q2_1max_(Constants::MaxEnergy2),
// 			  Q2_2min_(ZERO), Q2_2max_(Constants::MaxEnergy2),
// 			  F00_({0.274/GeV,0.274/GeV,0.344/GeV}),
// 			  LambdaP2_({0.6*GeV2,0.6*GeV2,0.6*GeV2})
//   {}

// public:

//   /** @name Virtual functions required by the MEBase class. */
//   //@{
//   /**
//    * Return the order in \f$\alpha_S\f$ in which this matrix
//    * element is given.
//    */
//   virtual unsigned int orderInAlphaS() const {
//     return 0;
//   }

//   /**
//    * Return the order in \f$\alpha_{EW}\f$ in which this matrix
//    * element is given.
//    */
//   virtual unsigned int orderInAlphaEW() const {
//     return 4;
//   }

//   /**
//    * The matrix element for the kinematical configuration
//    * previously provided by the last call to setKinematics(), suitably
//    * scaled by sHat() to give a dimension-less number.
//    * @return the matrix element scaled with sHat() to give a
//    * dimensionless number.
//    */
//   virtual double me2() const;

//   /**
//    * Return the scale associated with the last set phase space point.
//    */
//   virtual Energy2 scale() const;

//   /**
//    * Set the typed and momenta of the incoming and outgoing partons to
//    * be used in subsequent calls to me() and colourGeometries()
//    * according to the associated XComb object. If the function is
//    * overridden in a sub class the new function must call the base
//    * class one first.
//    */
//   virtual void setKinematics();

//   /**
//    * The number of internal degrees of freedom used in the matrix
//    * element.
//    */
//   virtual int nDim() const;

//   /**
//    * Generate internal degrees of freedom given nDim() uniform
//    * random numbers in the interval \f$ ]0,1[ \f$. To help the phase space
//    * generator, the dSigHatDR should be a smooth function of these
//    * numbers, although this is not strictly necessary.
//    * @param r a pointer to the first of nDim() consecutive random numbers.
//    * @return true if the generation succeeded, otherwise false.
//    */
//   virtual bool generateKinematics(const double * r);

//   /**
//    * Return the matrix element squared differential in the variables
//    * given by the last call to generateKinematics().
//    */
//   virtual CrossSection dSigHatDR() const;

//   /**
//    * Add all possible diagrams with the add() function.
//    */
//   virtual void getDiagrams() const;

//   /**
//    * Get diagram selector. With the information previously supplied with the
//    * setKinematics method, a derived class may optionally
//    * override this method to weight the given diagrams with their
//    * (although certainly not physical) relative probabilities.
//    * @param dv the diagrams to be weighted.
//    * @return a Selector relating the given diagrams to their weights.
//    */
//   virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

//   /**
//    * Return a Selector with possible colour geometries for the selected
//    * diagram weighted by their relative probabilities.
//    * @param diag the diagram chosen.
//    * @return the possible colour geometries weighted by their
//    * relative probabilities.
//    */
//   virtual Selector<const ColourLines *>
//   colourGeometries(tcDiagPtr diag) const;

//   /**
//    *  Construct the vertex of spin correlations.
//    */
//   virtual void constructVertex(tSubProPtr);
//   //@}
// public:

//   /** @name Functions used by the persistent I/O system. */
//   //@{
//   /**
//    * Function used to write out object persistently.
//    * @param os the persistent output stream written to.
//    */
//   void persistentOutput(PersistentOStream & os) const;

//   /**
//    * Function used to read in object persistently.
//    * @param is the persistent input stream read from.
//    * @param version the version number of the object when written.
//    */
//   void persistentInput(PersistentIStream & is, int version);
//   //@}

//   /**
//    * The standard Init function used to initialize the interfaces.
//    * Called exactly once for each class by the class description system
//    * before the main function starts or
//    * when this class is dynamically loaded.
//    */
//   static void Init();

// protected:

//   /** @name Clone Methods. */
//   //@{
//   /**
//    * Make a simple clone of this object.
//    * @return a pointer to the new object.
//    */
//   virtual IBPtr clone() const;

//   /** Make a clone of this object, possibly modifying the cloned object
//    * to make it sane.
//    * @return a pointer to the new object.
//    */
//   virtual IBPtr fullclone() const;
//   //@}

// protected:

//   /** @name Standard Interfaced functions. */
//   //@{
//   /**
//    * Initialize this object after the setup phase before saving an
//    * EventGenerator to disk.
//    * @throws InitException if object could not be initialized properly.
//    */
//   virtual void doinit();
//   //@}

// private:

//   /**
//    * The assignment operator is private and must never be called.
//    * In fact, it should not even be implemented.
//    */
//   MEee2eePseudoScalar & operator=(const MEee2eePseudoScalar &) = delete;

// private:

//   /**
//    *   fermion-fermion-photon vertex
//    */
//   AbstractFFVVertexPtr FFPVertex_;

//   /**
//    *  Pointer to the particle data object for the photon
//    */
//   PDPtr gamma_;

// private:

//   /**
//    *   Cuts
//    */
//   //@{
//   /**
//    *  Minimum value of \f$Q_1^2\f$
//    */
//   Energy2 Q2_1min_;
  
//   /**
//    *  Maximum value of \f$Q_1^2\f$
//    */
//   Energy2 Q2_1max_;
  
//   /**
//    *  Minimum value of \f$Q_2^2\f$
//    */
//   Energy2 Q2_2min_;
  
//   /**
//    *  Maximum value of \f$Q_2^2\f$
//    */
//   Energy2 Q2_2max_;
//   //@}

// private:

//   /**
//    *   Parameters for the form-factors
//    */
//   //@{
//   /**
//    *    Form factors at $Q_1^2=Q_2^2=0$
//    */
//   vector<InvEnergy> F00_;

//   /**
//    * Pole mass squared parameter for the form factors
//    */
//   vector<Energy2> LambdaP2_;
//   //@}

// private:

//   /**
//    *  Kinematic quantities stored for accuracy
//    */
//   //@{
//   /**
//    *  Outgoing electron
//    */
//   double cHalf1_,sHalf1_,phi1_;
//   Energy2 t1_;
  
//   /**
//    *  Outgoing positron
//    */
//   double cHalf2_,sHalf2_,phi2_;
//   Energy2 t2_;
//   //@}
  
// };

// }

// #endif /* Herwig_MEee2eePseudoScalar_H */
