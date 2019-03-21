// -*- C++ -*-
//
// ColourFlowBasis.h is a part of CVolver
// Copyright (C) 2013-2017 Simon Platzer, The Herwig Collaboration
//
// CVolver is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef CVOLVER_ColourFlowBasis_H
#define CVOLVER_ColourFlowBasis_H
//
// This is the declaration of the ColourFlowBasis class.
//

#include "Herwig/MatrixElement/Matchbox/Utility/ColourBasis.h"
#include "ColourFlows.h"

namespace CVolver {

using namespace ThePEG;
using namespace Herwig;

/**
 * Specify particle data traits for ThePEG
 */
template<>
struct ParticleDataTraits<ThePEG::PDT::Colour> {

  /**
   * Return true, if singlet
   */
  static bool isSinglet(const ThePEG::PDT::Colour& pd) { 
    return pd == PDT::Colour0;
  }

  /**
   * Return true, if anti-fundamental
   */
  static bool isAntiFundamental(const ThePEG::PDT::Colour& pd) {
    return pd == PDT::Colour3bar;
  }

  /**
   * Return true, if fundamental
   */
  static bool isFundamental(const ThePEG::PDT::Colour& pd) {
    return pd == PDT::Colour3;
  }

  /**
   * Return true, if adjoint
   */
  static bool isAdjoint(const ThePEG::PDT::Colour& pd) {
    return pd == PDT::Colour8;
  }

};

/**
 * ColourFlowBasis implements the colour flow basis.
 *
 * @see \ref ColourFlowBasisInterfaces "The interfaces"
 * defined for ColourFlowBasis.
 */
class ColourFlowBasis: public Herwig::ColourBasis {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ColourFlowBasis();

  /**
   * The destructor.
   */
  virtual ~ColourFlowBasis();
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
   * Map colour signatures to colour crossings
   */
  map<vector<PDT::Colour>,ColourFlowCrossing> theCrossings;

  /**
   * Colour flows indexed by basis size; note this is a unique
   * assignment
   */
  map<size_t,vector<ColourFlow> > theFlows;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ColourFlowBasis & operator=(const ColourFlowBasis &) = delete;

};

}

#endif /* CVOLVER_ColourFlowBasis_H */
