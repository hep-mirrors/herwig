// -*- C++ -*-
//
// MatchboxNLOME.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxNLOME_H
#define HERWIG_MatchboxNLOME_H
//
// This is the declaration of the MatchboxNLOME class.
//

#include "ThePEG/MatrixElement/MEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/InsertionOperators/MatchboxInsertionOperator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxNLOME assembles a MatchboxMEBase
 * and several MatchboxInsertionOperator objects to calculate
 * the Born + finite virtual cross section in a NLO
 * QCD calculation.
 *
 * @see \ref MatchboxNLOMEInterfaces "The interfaces"
 * defined for MatchboxNLOME.
 */
class MatchboxNLOME: public MEBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxNLOME();

  /**
   * The destructor.
   */
  virtual ~MatchboxNLOME();
  //@}

public:

  /** @name Subprocess and diagram information. */
  //@{

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const {
    theMatrixElement->diagrams();
    useDiagrams(theMatrixElement);
  }

  /**
   * Return true, if this matrix element does not want to
   * make use of mirroring processes; in this case all
   * possible partonic subprocesses with a fixed assignment
   * of incoming particles need to be provided through the diagrams
   * added with the add(...) method.
   */
  virtual bool noMirror () const { return true; }

  /**
   * With the information previously supplied with the
   * setKinematics(...) method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const  { 
    return matrixElement()->diagrams(dv); 
  }

  /**
   * Select a diagram. Default version uses diagrams(const
   * DiagramVector &) to select a diagram according to the
   * weights. This is the only method used that should be outside of
   * MEBase.
   */
  virtual DiagramIndex diagram(const DiagramVector & dv) const {
    return matrixElement()->diagram(dv); 
  }

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const { return matrixElement()->colourGeometries(diag); }

  /**
   * Select a ColpurLines geometry. The default version returns a
   * colour geometry selected among the ones returned from
   * colourGeometries(tcDiagPtr).
   */
  virtual const ColourLines &
  selectColourGeometry(tcDiagPtr diag) const { return matrixElement()->selectColourGeometry(diag); }

  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const {
    return 
      virtuals().empty() ? theMatrixElement->orderInAlphaS() :
      theMatrixElement->orderInAlphaS() + 1;
  }

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 0.
   */
  virtual unsigned int orderInAlphaEW() const { return theMatrixElement->orderInAlphaEW(); }

  //@}

  /** @name Phasespace generation */
  //@{

  /**
   * Set the XComb object to be used in the next call to
   * generateKinematics() and dSigHatDR().
   */
  virtual void setXComb(tStdXCombPtr);

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics() { matrixElement()->setKinematics(); }

  /**
   * Generate internal degrees of freedom given nDim() uniform random
   * numbers in the interval ]0,1[. To help the phase space generator,
   * the 'dSigHatDR' should be a smooth function of these numbers,
   * although this is not strictly necessary. The return value should
   * be true of the generation succeeded. If so the generated momenta
   * should be stored in the meMomenta() vector.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * The number of internal degreed of freedom used in the matrix
   * element. This default version returns 0;
   */
  virtual int nDim() const { if ( theNDim < 0 ) getNDim(); return theNDim; }

  /**
   * Get the ranom numbers needed
   */
  void getNDim() const;

  /**
   * Return true, if this matrix element will generate momenta for the
   * incoming partons itself.  The matrix element is required to store
   * the incoming parton momenta in meMomenta()[0,1]. No mapping in
   * tau and y is performed by the PartonExtractor object, if a
   * derived class returns true here. The phase space jacobian is to
   * include a factor 1/(x1 x2).
   */
  virtual bool haveX1X2() const { return theMatrixElement->haveX1X2(); }

  /**
   * Return true, if this matrix element expects
   * the incoming partons in their center-of-mass system
   */
  virtual bool wantCMS () const { return theMatrixElement->wantCMS(); }

  /**
   * Return true, if the XComb steering this matrix element
   * should keep track of the random numbers used to generate
   * the last phase space point
   */
  virtual bool keepRandomNumbers() const { return true; }

  /**
   * Clear the information previously provided by a call to
   * setKinematics(...).
   */
  virtual void clearKinematics() { matrixElement()->clearKinematics(); }

  //@}

  /** @name Scale choices, couplings and PDFs */
  //@{

  /**
   * Return the scale associated with the phase space point provided
   * by the last call to setKinematics().
   */
  virtual Energy2 scale() const { return matrixElement()->scale(); }

  /**
   * Return the value of \f$\alpha_S\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaS(scale()).
   */
  virtual double alphaS() const { return matrixElement()->alphaS(); }

  /**
   * Return the value of \f$\alpha_EM\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaEM(scale()).
   */
  virtual double alphaEM() const { return matrixElement()->alphaEM(); }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the first incoming parton itself.
   */
  virtual bool havePDFWeight1() const { return theMatrixElement->havePDFWeight1(); }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the second incoming parton itself.
   */
  virtual bool havePDFWeight2() const { return theMatrixElement->havePDFWeight2(); }

  //@}

  /** @name Matrix elements and evaluation */
  //@{

  /**
   * Return the Born matrix element
   */
  Ptr<MatchboxMEBase>::tcptr matrixElement() const { return theMatrixElement; }

  /**
   * Return the Born matrix element
   */
  Ptr<MatchboxMEBase>::tptr matrixElement() { return theMatrixElement; }

  /**
   * Set the Born matrix element
   */
  void matrixElement(Ptr<MatchboxMEBase>::ptr b) { theMatrixElement = b; }

  /**
   * Return the virtual corrections
   */
  const vector<Ptr<MatchboxInsertionOperator>::ptr>& virtuals() const {
    return theVirtuals;
  }

  /**
   * Return the virtual corrections
   */
  vector<Ptr<MatchboxInsertionOperator>::ptr>& virtuals() {
    return theVirtuals;
  }

  /**
   * Return true, if cancellationn of epsilon poles should be checked.
   */
  bool checkPoles() const { return theCheckPoles; }

  /**
   * Switch on checking of epsilon pole cancellation.
   */
  void doCheckPoles() { theCheckPoles = true; }

  /**
   * Perform the check of epsilon pole cancellation. Results will be
   * written to the log file, one per phasespace point.
   */
  void logPoles() const;

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   */
  virtual double me2() const;

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;

  //@}

  /** @name Caching and diagnostic information */
  //@{

  /**
   * Inform this matrix element that a new phase space
   * point is about to be generated, so all caches should
   * be flushed.
   */
  virtual void flushCaches();

  /**
   * Dump the setup to an ostream
   */
  void print(ostream&) const;

  /**
   * Print information on last event generated.
   */
  void printLastEvent(ostream& os) const;

  //@}

  /** @name Methods used to setup NLO ME objects */
  //@{

  /**
   * Clone the dependencies, using a given prefix.
   */
  void cloneDependencies(const std::string& prefix = "");

  //@}

  /** @name Methods required to setup the event record */
  //@{

  /**
   * construct the spin information for the interaction
   */
  virtual void constructVertex(tSubProPtr sub) { matrixElement()->constructVertex(sub); }

  /**
   * Comlete a SubProcess object using the internal degrees of freedom
   * generated in the last generateKinematics() (and possible other
   * degrees of freedom which was intergated over in dSigHatDR(). This
   * default version does nothing. Will be made purely virtual in the
   * future.
   */
  virtual void generateSubCollision(SubProcess & sub) { matrixElement()->generateSubCollision(sub); }

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

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The Born matrix element
   */
  Ptr<MatchboxMEBase>::ptr theMatrixElement;

  /**
   * The virtual corrections.
   */
  vector<Ptr<MatchboxInsertionOperator>::ptr> theVirtuals;

  /**
   * The number of random variables needed
   */
  mutable int theNDim;

  /**
   * True, if cancellationn of epsilon poles should be checked.
   */
  bool theCheckPoles;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxNLOME & operator=(const MatchboxNLOME &);

};

}

#endif /* HERWIG_MatchboxNLOME_H */
