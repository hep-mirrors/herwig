// -*- C++ -*-
//
// TraceBasis.h is a part of ColorFull
// Copyright (C) 2010-2011 Simon Platzer & Malin Sjodahl
//
// ColorFull is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef COLORFULL_TraceBasis_H
#define COLORFULL_TraceBasis_H
//
// This is the declaration of the TraceBasis class.
//

#include "Herwig/MatrixElement/Matchbox/Utility/ColourBasis.h"
#include "Trace_basis.h"

namespace ColorFull {

using namespace ThePEG;
using namespace Herwig;

/**
 * TraceBasis implements the trace colour basis.
 *
 * @see \ref TraceBasisInterfaces "The interfaces"
 * defined for TraceBasis.
 */
class TraceBasis: public Herwig::ColourBasis {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  TraceBasis();

  /**
   * The destructor.
   */
  virtual ~TraceBasis();
  //@}

public:

  /**
   * Clear this colour basis
   */
  virtual void clear();

  /**
   * Return a map of basis tensor indices to vectors identifying a
   * certain ordering corresponding to the given colour structure. May
   * not be supported by all colour basis implementations.
   */
  virtual map<size_t,vector<vector<size_t> > > basisList(const vector<PDT::Colour>&) const;

  /**
   * Prepare the basis for the normal ordered legs and return the
   * dimensionality of the basis.
   */
  virtual size_t prepareBasis(const vector<PDT::Colour>&);

  /**
   * Gather any implementation dependend details when reading a basis
   */
  virtual void readBasisDetails(const vector<PDT::Colour>&);

  /**
   * Return the scalar product of basis tensors labelled a and b in
   * the basis used for the given normal ordered legs.
   */
  virtual double scalarProduct(size_t a, size_t b,
			       const vector<PDT::Colour>& abBasis) const;

  /**
   * Return the matrix element of a colour charge
   * <c_{n+1,a}|T_i|c_{n,b}> between basis tensors a and b, with
   * respect to aBasis and bBasis
   */
  virtual double tMatrixElement(size_t i, size_t a, size_t b,
				const vector<PDT::Colour>& aBasis,
				const vector<PDT::Colour>& bBasis) const;

  /**
   * Return true, if the colour basis is capable of assigning colour
   * flows.
   */
  virtual bool haveColourFlows() const { return true; }

  /**
   * Return true, if a large-N colour connection exists for the
   * given external legs and basis tensor.
   */
  virtual bool colourConnected(const cPDVector&,
			       const vector<PDT::Colour>&,
			       const pair<int,bool>&, 
			       const pair<int,bool>&, 
			       size_t) const;

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The color functions object to be used
   */
  mutable Col_functions colorFunctions;

  /**
   * Map legs to known basis vectors.
   */
  mutable map<vector<PDT::Colour>,Trace_basis> theBasisMap;

  /**
   * Memorize scalar product intermediate results.
   */
  mutable map<string,Polynomial> theScalarProducts;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TraceBasis & operator=(const TraceBasis &);

};

}

#endif /* COLORFULL_TraceBasis_H */
