// -*- C++ -*-
//
// MRST.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MRST_H
#define HERWIG_MRST_H

#include <ThePEG/PDF/PDFBase.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace Herwig {

using namespace ThePEG;

/** \ingroup PDF
 *
 *  Implementation of the MRST PDFs
 *
 * @see \ref MRSTInterfaces "The interfaces"
 * defined for MRST.
 */
class MRST : public PDFBase {
  /**
   *  Enumeration to storage the types of partons
   */
  enum PDFFlavour { upValence = 1, dnValence, glu, upSea, chm, str, bot, dnSea };

  /**
   *  Enum for type of pdf to return
   */
  enum PDFType {Sea,Valence,Total};

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  MRST();
  //@}
 
public:

  /** @name Virtual functions from PDFBase */
  //@{
  /**
   * Return true if this PDF can handle the extraction of parton from the
   * given particle ie. if the particle is a proton or neutron.
   * @param particle The particle
   */
  virtual bool canHandleParticle(tcPDPtr particle) const;

  /**
   * Return the parton types which are described by these parton
   * densities.
   * @param p The particle
   */
  virtual cPDVector partons(tcPDPtr p) const;

  /**
   * Return x times the pdf for the given parameters
   * @param particle The beam particle
   * @param parton The parton for which to return the PDF.
   * @param partonScale The scale at which to evaluate the PDF.
   * @param x The momentum fraction
   * @param eps ??? an unknown parameter from ThePEG.
   * @param particleScale The scale for the particle
   */
  virtual double xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                     double x, double eps = 0.0,
                     Energy2 particleScale = ZERO) const;

  /**
   * Return x times the valence pdf for the given parameters
   * @param particle The beam particle
   * @param parton The parton for which to return the PDF.
   * @param partonScale The scale at which to evaluate the PDF.
   * @param x The momentum fraction
   * @param eps ??? an unknown parameter from ThePEG.
   * @param particleScale The scale for the particle
   */
  virtual double xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                      double x, double eps = 0.0,
                      Energy2 particleScale = ZERO) const;

  /**
   * The sea density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and momentum fraction \a x. The \a particle is
   * assumed to have a virtuality \a particleScale. If not overidden
   * by a sub class this implementation will assume that the
   * difference between a quark and anti-quark distribution is due do
   * valense quarks.
   */
  virtual double xfsx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		      double x, double eps = 0.0,
		      Energy2 particleScale = ZERO) const;
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
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MRST> initMRST;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MRST & operator=(const MRST &);

private:

  /**
   * This function interpolates the PDF value for the given 
   * bin numbers and fractional positions inside the bin. It is a helper function for pdfValue() 
   * 
   * @param i The PDF flavour
   * @param n The x bin index
   * @param m The q^2 bin index
   * @param u The fractional position along the current x bin
   * @param t The fractional position along the current q^2 bin
   */
  inline double lookup(PDFFlavour i, int n, int m, double u, double t) const;


  /**
   * This function calculates the PDF value for the given particles and a given x and q
   * @param x The energy fraction
   * @param q2 The scale
   * @param particle The beam particle
   * @param parton The parton for which to return the PDF.
   * @param type Type of PDF, sea, valence or total.
   */
  double pdfValue(double x, Energy2 q2, 
		  tcPDPtr particle, tcPDPtr parton,PDFType type) const;

  /**
   * Returns an integer j such that x lies inbetween xx[j] and xx[j+1].
   * @param xx The x values
   * @param n  The number of values
   * @param x  The x value
   */
  inline int locate(double xx[],int n,double x) const;

  /**
   *  Return  the estimate of the derivative at \f$x_2\f$ obtained by a polynomial 
   * interpolatio using the three points \f$(x_i,y_i)\f$
   * @param x1 The \f$x\f$ value at the first point
   * @param x2 The \f$x\f$ value at the second point
   * @param x3 The \f$x\f$ value at the third point
   * @param y1 The \f$y\f$ value at the first point
   * @param y2 The \f$y\f$ value at the second point
   * @param y3 The \f$y\f$ value at the third point
   */
  inline double polderivative(double x1, double x2, double x3,
			      double y1, double y2, double y3) const;

  /**
   *  Read the data from the file
   */
  virtual void readSetup(istream &);

  /**
   *  Initialize the data
   * @param reread Whether or not to reread the data
   */
  void initialize(bool reread = true);

private:

  /**
   *  Parameters for the MRST data tables
   */
  //@{
  /**
   *  Number of distributions to interpolate
   */
  static const int np=8;
  /**
   *  Number of points in \f$x\f$ for the interpolation
   */
  static const int nx=49;
  
  /**
   *  Number of points in \f$q^2\f$ for the interpolation
   */
  static const int nq=37;

  /**
   *  \f$q^2\f$ bin where charm introduced
   */
  static const int nqc0=2; 
  
  /**
   *  \f$q^2\f$ bin where bottom introduced
   */
  static const int nqb0=11;

  /**
   *  Parameter for the FORTRAN interpolation
   */
  static const int ntenth=23;
  
  /**
   *  Minimum value of \f$x\f$
   */
  static const double xmin;

  /**
   *  Maximum value of \f$x\f$
   */
  static const double xmax;
  
  /**
   *  Minimum value of \f$q^2\f$.
   */
  static const Energy2 qsqmin;
  
  /**
   *  Maximum value of \f$q^2\f$.
   */
  static const Energy2 qsqmax;
  
  /**
   *  Mass squared of the charm quark
   */
  static const Energy2 mc2;
  
  /**
   *  Mass squared of the bottom quark
   */
  static const Energy2 mb2;
  //@}

  /**
   *  Use FORTRAN or C++ MRST interpolation
   */
  unsigned _inter;

  /**
   *  X value to switch from cubic to linear
   */
  double _xswitch;

  /**
   *  The name of the file
   */
  string _file;

  /**
   *  Array containing the data to be interpolated
   */
  //  double data[np+1][nx+1][nq+1];
  vector<vector<vector<double> > > data;

  /**
   *  Array containing the data to be interpolated
   */
  //  double data[np+1][nx+1][nq+1];
  vector<vector<vector<double> > > fdata;

  /**
   *  The \f$x\f$ values for interpolation
   */
  static double xx[nx+1];

  /**
   *  The \f$x\f$ values for interpolation
   */
  static double lxx[nx+1];

  /**
   *  The \f$x\f$ values for interpolation
   */
  static double lxxb[nx+1];

  /**
   *  The \f$q^2\f$ values for interpolation
   */
  static double qq[nq+1];

  /**
   *  The \f$q^2\f$ values for interpolation
   */
  static double lqq[nq+1];

  /**
   *  Coefficients used for interpolation
   */
  double c[np+1][nx][nq][5][5];

  /**
   *  The powers n0 for the FORTRAN interpolation
   */
  static double n0[np+1];

  /**
   *  where or not initialized
   */
  static bool initialized;
};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MRST. */
template <>
struct BaseClassTrait<Herwig::MRST,1> {
  /** Typedef of the first base class of MRST. */
  typedef PDFBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MRST class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MRST>: public ClassTraitsBase<Herwig::MRST> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MRST"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MRST is implemented. It may also include several, space-separated,
   * libraries if the class MRST depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMRST.so"; }
};

/** @endcond */

}

#include "MRST.icc"

#endif
