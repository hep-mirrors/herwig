// -*- C++ -*-
#ifndef PYTHIA7_MEee2gZ2bb_H
#define PYTHIA7_MEee2gZ2bb_H
//
// This is the declaration of the <!id>MEee2gZ2bb<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// The <!id>MEee2gZ2qq<!!id> class implements the
// <i>e<sup>+</sup>e<sup>-</sup> -&gt; gamma/Z<sup>0</sup> -&gt;
// b+bbar</i> matrix element. Both the continuum and
// <i>Z<sup>0</sup></i> pole term as well as the interference term is
// included. Although not a strict QCD matrix element the cass
// inherits from <!class>ME2to2Base<!!class>, mainly to inherit the
// parameter for the number of active quark flavours.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:ME2to2QCD.html">ME2to2QCD.h</a>.
// 

#include "Pythia7/MatrixElement/ME2to2QCD.h"

namespace Pythia7 {

class MEee2gZ2bb: public ME2to2QCD {

public:

  MEee2gZ2bb();
  MEee2gZ2bb(const MEee2gZ2bb &);
  virtual ~MEee2gZ2bb();
  // Standard ctors and dtor

public:

  virtual unsigned int orderInAlphaS() const;
  virtual unsigned int orderInAlphaEW() const;
  // Return the order in respective couplings in which this matrix
  // element is given. Returns 0 and 2 respectively.

  virtual double me2() const;
  // Return the matrix element for the kinematical configuation
  // previously provided by the last call to setKinematics().

  virtual void getDiagrams() const;
  // Add all possible diagrams with the add() function.

  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
  // Return a Selector with possible colour geometries for the selected
  // diagram weighted by their relative probabilities.

  virtual Selector<DiagramIndex> diagrams(const DiagramVector &) const;
  // Weight the given diagrams with their relative probabilities.

  virtual Energy2 scale() const;
  // Return the scale associated with the last set phase space point.

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interface.

protected:

  virtual IBPtr clone() const;
  virtual IBPtr fullclone() const;
  // Standard clone methods

  virtual void doinit() throw(InitException);
  // Standard Interfaced virtual functions.

protected:

  vector<double> coefs;
  // Constants for the different terms set from the StandardModel in
  // the init() function.

  Energy2 mZ2;
  Energy2 GZ2;
  // The mass squared and width squared of the Z0.

  mutable double lastCont;
  mutable double lastBW;
  // The last continuum and Breit-Wigner terms to be used to select
  // primary diagram.

private:

  static ClassDescription<MEee2gZ2bb> initMEee2gZ2bb;

  MEee2gZ2bb & operator=(const MEee2gZ2bb &);
  //  Private and non-existent assignment operator.

};

template <>
struct BaseClassTrait<MEee2gZ2bb,1> {
  typedef ME2to2QCD NthBase;
};

template <>
struct ClassTraits<MEee2gZ2bb>: public ClassTraitsBase<MEee2gZ2bb> {
  static string className() { return "/Herwig++/MEee2gZ2bb"; }
  static string library() { return "HwME.so"; }
};

}

#endif /* PYTHIA7_MEee2gZ2qq_H */
