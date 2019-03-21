// -*- C++ -*-
#ifndef HERWIG_MEDiffraction_H
#define HERWIG_MEDiffraction_H
//
// This is the declaration of the MEDiffraction class.
//

#include "Herwig/MatrixElement/HwMEBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEDiffraction class provides a simple colour singlet exchange matrix element
 * to be used in the soft component of the multiple scattering model of the 
 * underlying event
 *
 * @see \ref MEDiffractionInterfaces "The interfaces"
 * defined for MEDiffraction.
 */
class MEDiffraction: public HwMEBase {

public:

  MEDiffraction();

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

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

  /**
   * Expect the incoming partons in the laboratory frame
   */
  /* virtual bool wantCMS() const { return false; } */

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

  /**
   * Initialize this object. Called in the run phase just before a run begins.
   */
  virtual void doinitrun();
  //@}

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

  /* The matrix element squared */
  double theme2;
  
  /* Use only delta as excited state */
  bool deltaOnly; 

  
  /* Direction of the excited proton */
  unsigned int diffDirection; 
  
  /* Number of clusters the dissociated proton decays into */
  unsigned int dissociationDecay; 
  
  /* The mass of the consitutent quark */
  Energy mq() const {return Energy(0.325*GeV);}
  
  /* The mass of the constituent diquark */
  Energy mqq() const {return Energy(0.650*GeV);}
  
  /* The proton-pomeron slope */
  double theprotonPomeronSlope;
  
  /* The soft pomeron intercept */
  double thesoftPomeronIntercept;
  
  /* The soft pomeron slope */
  double thesoftPomeronSlope;
 
  
  /**
   * Sample the diffractive mass squared M2 and the momentum transfer t
   */
  pair<pair<Energy2,Energy2>,Energy2> diffractiveMassAndMomentumTransfer() const;
  

  /**
   * Random value for the diffractive mass squared M2 according to (M2/s0)^(-intercept)
   */
  Energy2 randomM2() const;

  /**
   * Random value for t according to exp(diffSlope*t)
   */
  Energy2 randomt(Energy2 M2) const;
  
  /**
   * Random value for t according to exp(diffSlope*t) for double diffraction
   */
  
  Energy2 doublediffrandomt(Energy2 M12, Energy2 M22) const;
  
  
  /**
   * Returns the momenta of the two-body decay of momentum pp
   */
  pair<Lorentz5Momentum,Lorentz5Momentum> twoBodyDecayMomenta(Lorentz5Momentum pp) const;

  /**
   * Returns the proton-pomeron slope
   */
  InvEnergy2 protonPomeronSlope() const; 
  
  /**
   * Returns the soft pomeron intercept
   */  
  double softPomeronIntercept() const; 
 
  //M12 and M22 are masses squared of 
  //outgoing particles
  
  /**
   * Returns the minimal possible value of momentum transfer t given the center
   * of mass energy and diffractive masses
   */
  Energy2 tminfun(Energy2 s, Energy2 M12, Energy2 M22) const;
 
  /**
   * Returns the maximal possible value of momentum transfer t given the center
   * of mass energy and diffractive masses
   */
  Energy2 tmaxfun(Energy2 s , Energy2 M12, Energy2 M22) const; 

  /**
   * Returns the minimal possible value of diffractive mass
   */
  //lowest possible mass given the constituent masses of quark and diquark
  Energy2 M2min() const{return sqr(getParticleData(2212)->mass()+mq()+mqq());}
  
  /**
   * Returns the maximal possible value of diffractive mass
   */
  Energy2 M2max() const{
    return sqr(generator()->maximumCMEnergy()-getParticleData(2212)->mass());
  }//TODO:modify to get proper parameters

  InvEnergy2 softPomeronSlope() const;
  
   

  /* Kallen function */
  template<int L, int E, int Q, int DL, int DE, int DQ>
  Qty<2*L,2*E,2*Q,DL,DE,DQ> kallen(Qty<L,E,Q,DL,DE,DQ> a,
                                   Qty<L,E,Q,DL,DE,DQ> b,
                                   Qty<L,E,Q,DL,DE,DQ> c) const {
    return a*a + b*b + c*c - 2.0*(a*b + b*c + c*a);
  }
  
  

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEDiffraction> initMEDiffraction;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEDiffraction & operator=(const MEDiffraction &) = delete;

  bool isInRunPhase; 
  
 
  /* The proton mass */
  Energy theProtonMass;
  
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEDiffraction. */
template <>
struct BaseClassTrait<Herwig::MEDiffraction,1> {
  /** Typedef of the first base class of MEDiffraction. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEDiffraction class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEDiffraction>
  : public ClassTraitsBase<Herwig::MEDiffraction> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEDiffraction"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEDiffraction is implemented. It may also include several, space-separated,
   * libraries if the class MEDiffraction depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() {return "HwMEHadron.so";}
};

/** @endcond */

}

#endif /* HERWIG_MEDiffraction_H */
