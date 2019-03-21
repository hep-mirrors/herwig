// -*- C++ -*-
//
// MEHiggsPair.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2009-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEHiggsPairOL_H
#define HERWIG_MEHiggsPairOL_H
//
// This is the declaration of the MEHiggsPairOL class.
//
// The implementation of this process is based upon hep-ph/0112161 by G.F. Giudice, R. Rattazzi, J.D. Wells.

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "HiggsPair.h"


namespace Herwig {
using namespace ThePEG;


/* 
 * Interface to external FORTRAN hpair functions (T.Plehn, M.Spira & P.Zerwas, arXiv:hep-ph/9603205)
 * original Program found at http://people.web.psi.ch/spira/hpair/ 
 */
extern "C" { 
  complex<double> eta_(complex<double>* C1, complex<double>* C2);
  complex<double> d04_(double* P1, double* P2, double* P3, double* P4, double* P12,double* P23, double* M1, double* M2, double* M3, double* M4);
  complex<double> c03_(double* P1,double* P2,double* P3,double* M1,double* M2,double* M3); 
  complex<double> cspen_(double* Z);
  complex<double> sqe_(double* A, double* B, double* C);
  complex<double> etas_(complex<double>*Y,complex<double>* R,complex<double>*RS);
  void formfac_(double* AMQ, double* S, double* T, double* U, double* M1, double* M2, complex<double>* C0AB,complex<double>* C0AC,complex<double>* C0AD,complex<double>* C0BC,complex<double>* C0BD,complex<double>* C0CD,complex<double>* D0ABC,complex<double>* D0BAC,complex<double>* D0ACB);
}

extern "C" {
  void ol_setparameter_int(const char* param, int val);
  void ol_setparameter_double(const char* param, double val);
  void ol_setparameter_string(const char* param, char* val);
  int ol_register_process(const char* process, int amptype);
  int ol_n_external(int id);
  void ol_phase_space_point(int id, double sqrt_s, double* pp);
  void ol_start();
  void ol_finish();
  void ol_evaluate_tree(int id, double* pp, double* m2_tree);
  void ol_evaluate_loop(int id, double* pp, double* m2_tree, double* m2_loop, double* acc);
  void ol_evaluate_loop2(int id, double* pp, double* m2_loop2, double* acc);
}

extern "C" { 
  extern struct{
      complex<double> A1;
      complex<double> A2;
      complex<double> H1;
      complex<double> H2;
      complex<double> Z1;
      complex<double> Z2; 
      complex<double> AA1;
      complex<double> AA2;
      complex<double> HH1;
      complex<double> HH2; 
      complex<double> AH1;
      complex<double> AH2;
  } form_;
}   

/**
 * The MEHiggsPairOL class implements the matrix elements for
 * HiggsPairian \f$2\to2\f$ scattering process
 */
class MEHiggsPairOL: public HwMEBase {


public:

  /**
   * The default constructor.
   */
  MEHiggsPairOL();

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const { return 0; }

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const { return 0; }

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

 /**
   * Generate internal degrees of freedom 
   */

  virtual bool generateKinematics(const double * r);


  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;


  

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

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;
  //@}


  virtual void setKinematics();

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;


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
   * Initialize this object before the run phase
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinitrun();


  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

  /**
   *  Access to the higgs data
   */ 
  PDPtr higgs() const { return _higgs; }

  /**
   *  Set the higgs data
   */ 
  void higgs(PDPtr in) {_higgs =in;}


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

  /** @name Helper functions for me2. */
  //@{
  /**
   */
 
private:
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEHiggsPairOL> initMEHiggsPairOL;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEHiggsPairOL & operator=(const MEHiggsPairOL &) = delete;


  /**
   * process identifiers
   */ 
  int _id, _idtri, _idbox, _idint;

  /**
   *  The higgs boson
   */
  PDPtr _higgs;

  /* 
   * The top mass
   */
  Energy _topmass;

  /* 
   * The bottom mass
   */
 Energy _bottommass;

  /*
   * The Z boson mass
   */ 
  Energy _zmass;
  

  /*					       
   * Higgs boson mass(es) 
   */
  Energy _m1, _m2;

  /*					       
   * Heavy H mass
   */
  Energy _heavyHmass;


  /*					       
   * Heavy H width
   */
  Energy _heavyHwidth;




private:
  /**
   * Pointer to the model.
   */
  tcSMPtr _theModel;

  /**
   *  The mass generator for the Higgs
   */
  GenericMassGeneratorPtr _hmass;

  /**
   *  multiplier for the SM triple-coupling
   */
  double _selfcoupling;


  /**
   *  Processes to include
   */
  unsigned int _process;

 /**
   *  On-shell mass for the higgs
   */
  Energy _mh;


  /**
   *  On-shell width for the higgs
   */
  Energy _wh;

  /*
   * Fix alphaS
   */ 
  unsigned int _fixedalphaS;

 /*
   * Value of AlphaS if fixed using option 2 above
   */ 
  double _alphasfixedvalue;

  /*
   * Fix scale of whole process
   */ 
  unsigned int _fixedscale;
  
  /*
   * Scale to use for alpha_S if fixed
   */
  Energy _alphascale;

  /*
   * Choose between the OpenLoops and HPAIR implementations
   */ 
  unsigned int _implementation;

  /*
   * include the widths
   */
  unsigned int _includeWidths;

  /*
   * include the b quark loops
   */
  unsigned int _includeBquark;


  /* 
   * Base scale to use if chosen to be fixed
   */ 
  Energy _basescale;

  /* 
   * scale multiplier 
   */ 
  double _scalemultiplier;


  /* 
   * Masses of fermions 
   */
  double Mass_E;
  double Mass_M;
  double Mass_L;
  double Mass_T;
  double Mass_U;
  double Mass_C;
  double Mass_D;
  double Mass_S;
  double Mass_B;
  double Mass_Z;
  double Mass_W;
  double Mass_H;
  double Width_C;
  double Width_B;
  double Width_T;
  double Width_W;
  double Width_Z;
  double Width_H;


  double getMass_E() const { return Mass_E; }
  double getMass_M() const { return Mass_M; }
  double getMass_L() const { return Mass_L; }

  double getMass_T() const { return Mass_T; }
  double getMass_U() const { return Mass_U; }  
  double getMass_C() const { return Mass_C; }
  double getMass_D() const { return Mass_D; }
  double getMass_S() const { return Mass_S; }
  double getMass_B() const { return Mass_B; }

  /*
   * Masses of bosons
   */
  double getMass_Z() const { return Mass_Z; }
  double getMass_W() const { return Mass_W; }
  double getMass_H() const { return Mass_H; }

  /*
   * Widths of particles
   */ 
  double getWidth_C() const { return Width_C; }
  double getWidth_B() const { return Width_B; }
  double getWidth_T() const { return Width_T; }
  double getWidth_W() const { return Width_W; }
  double getWidth_Z() const { return Width_Z; }
  double getWidth_H() const { return Width_H; }


  /*
   * QCD and QED couplings 
   */
  double Coupl_Alpha_QED; 
  double Coupl_Alpha_QCD; 

  
  /* Scalar integral initialization from hpair.f (T.Plehn, M.Spira & P.Zerwas, arXiv:hep-ph/9603205)
     Program found at http://people.web.psi.ch/spira/hpair/ */
  
  virtual vector<Complex> iniscal(double AMQ, double S, double T,double U, double M1, double M2) const;
  
  /* Matrix element calculation from hpair */

  virtual double MATRIX(double S, double T,double U, double M1, double M2) const;


 /* OpenLoops initialization */
  virtual void InitOpenLoops();

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEHiggsPairOL. */
template <>
struct BaseClassTrait<Herwig::MEHiggsPairOL,1> {
  /** Typedef of the first base class of MEHiggsPairOL. */
  typedef Herwig::HwMEBase NthBase;

};

/** This template specialization informs ThePEG about the name of
 *  the MEHiggsPairOL class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEHiggsPairOL>
  : public ClassTraitsBase<Herwig::MEHiggsPairOL> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEHiggsPairOL"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEHiggsPairOL is implemented. It may also include several, space-separated,
   * libraries if the class MEQCD2to2Fast depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HiggsPair.so"; }
};

/** @endcond */

}


#endif /* HERWIG_MEHiggsPairOL_H */
