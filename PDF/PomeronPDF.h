// -*- C++ -*-
//
// PomeronPDF.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PomeronPDF_H
#define HERWIG_PomeronPDF_H

#include <ThePEG/PDF/PDFBase.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace Herwig {

using namespace ThePEG;

/** \ingroup PDF
 *
 *  Implementation of the PomeronPDF PDFs
 *
 * @see \ref PomeronPDFInterfaces "The interfaces"
 * defined for PomeronPDF.
 */
class PomeronPDF : public PDFBase {
public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  PomeronPDF();
 

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
   *  Enumeration to storage the types of partons
   */
  enum PDFFlavour {charm, gluon, singlet};

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PomeronPDF & operator=(const PomeronPDF &) = delete;

private:

  /**
   * This function calculates the PDF value for the given particles and a given x and q
   * @param flPDF Flavour of the PDF function
   * @param x     The x of the pomeron
   * @param qq    The scale
   *
   */
  double getPDFValue(PDFFlavour flPDF, double x, Energy2 qq) const;

  /**
   * Load the PDF tables from file to the vectors.
   */
  void loadTables() const;

private:

  /**
   *  Number of PDF flavours
   */
  static const int nPDFFlavour_;

  
 /**
  * Vector of matrixes of PDF functions for different flavours:
  * pdfTable [Flavours][Q^2][log x]
  */
  mutable vector<vector<vector<double> > >     pdfTable_;

 /**
  * Vector of grid with log x values: lxGrid [Flavours] [log x] 
  */
  mutable vector<vector<double> >              lxGrid_;

 /**
  * Vector of grid with log qq values: lxGrid [Flavours] [log qq] 
  */
  mutable vector<vector<double> >              lqqGrid_;
 
 /**
  * Vector of names of PDF data files for each flavour 
  */
  vector<string> fileName_;

  /**
   *  Base of the filename
   */
  string rootName_;

   /**
   *  Number of ln x points
   */
  int nxPoints_;
  
  /**
   *  Number of qq flavours
   */
  int nqPoints_;

  /**
   *  Switch between different PDF fits. Possible values: 0 fit 2007, 1 fit A 2006, 2 fit B 2006 
   */
  int PDFFit_;

  /**
   *  Switch between different aproaches when the values are out of PDF range: 0 Extrapolate, 1 Freeze
   */
  int boundary_;


};

}

#endif
