// -*- C++ -*-
//
// SimpleColourBasis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_SimpleColourBasis_H
#define Herwig_SimpleColourBasis_H
//
// This is the declaration of the SimpleColourBasis class.
//

#include "Herwig/MatrixElement/Matchbox/Utility/ColourBasis.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief SimpleColourBasis implements the colour algebra needed for
 * electroweak boson and electroweak boson + jet production at NLO. It
 * mainly serves as an example for the general ColourBasis interface.
 *
 */
class SimpleColourBasis: public ColourBasis {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SimpleColourBasis();

  /**
   * The destructor.
   */
  virtual ~SimpleColourBasis();
  //@}

public:

  /**
   * Prepare the basis for the normal ordered legs and return the
   * dimensionality of the basis.
   */
  virtual size_t prepareBasis(const vector<PDT::Colour>&);

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
   * Return true, if a large-N colour connection exists for the
   * given external legs and basis tensor.
   */
  virtual bool colourConnected(const cPDVector&,
			       const vector<PDT::Colour>&,
			       const pair<int,bool>&, 
			       const pair<int,bool>&, 
			       size_t) const;

  /**
   * Return true, if the colour basis is capable of assigning colour
   * flows.
   */
  virtual bool haveColourFlows() const { return true; }

  /**
   * Create ids for bases
   */
  void makeIds() const;

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
   * id for 88
   */
  mutable vector<PDT::Colour> id88;

  /**
   * id for 33bar
   */
  mutable vector<PDT::Colour> id33bar;

  /**
   * id for 888
   */
  mutable vector<PDT::Colour> id888;

  /**
   * id for 33bar8
   */
  mutable vector<PDT::Colour> id33bar8;

  /**
   * id for 8888
   */
  mutable vector<PDT::Colour> id8888;

  /**
   * id for 33bar88
   */
  mutable vector<PDT::Colour> id33bar88;

  /**
   * id for 33bar33bar
   */
  mutable vector<PDT::Colour> id33bar33bar;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SimpleColourBasis & operator=(const SimpleColourBasis &);

};

}

#endif /* Herwig_SimpleColourBasis_H */
