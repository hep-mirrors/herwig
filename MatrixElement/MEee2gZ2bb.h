// -*- C++ -*-
#ifndef HERWIG_MEee2gZ2bb_H
#define HERWIG_MEee2gZ2bb_H
//
// This is the declaration of the MEee2gZ2bb class.


#include <ThePEG/MatrixElement/ME2to2QCD.h>

namespace ThePEG {

/** \ingroup MatrixElement
 * 
 *  The MEee2gZ2qq class implements the e+ e- -> gamma/Z0 b bbar 
 *  matrix element. Both the continuum and Z0 pole term as well as the 
 *  interference term is included. Although not a strict QCD matrix 
 *  element the cass inherits from ME2to2Base, mainly to inherit the
 *  parameter for the number of active quark flavours.
 *
 *  @see ME2to2QCD
 */
class MEee2gZ2bb: public ME2to2QCD {

public:

  /**
   * Standard ctors and dtor.
   */
  MEee2gZ2bb();
  MEee2gZ2bb(const MEee2gZ2bb &);
  virtual ~MEee2gZ2bb();

public:

  /**
   * Return the order in respective couplings in which this matrix
   * element is given. Returns 0 and 2 respectively.
   */
  virtual unsigned int orderInAlphaS() const;
  virtual unsigned int orderInAlphaEW() const;

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics().
   */
  virtual double me2() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   * Weight the given diagrams with their relative probabilities.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector &) const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

public:

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

  /**
   * Standard clone methods.
   */
  virtual IBPtr clone() const;
  virtual IBPtr fullclone() const;

  /**
   * Standard Interfaced virtual functions.
   */
  virtual void doinit() throw(InitException);

protected:
 
  /**
   * Constants for the different terms set from the StandardModel in
   * the init() function.
   */
  vector<double> coefs;

  /**
   * The mass squared and width squared of the Z0.
   */
  Energy2 mZ2;
  Energy2 GZ2;

  /**
   * The last continuum and Breit-Wigner terms to be used to select
   * primary diagram.
   */
  mutable double lastCont;
  mutable double lastBW;

private:

  static ClassDescription<MEee2gZ2bb> initMEee2gZ2bb;

  /**
   * Private and non-existent assignment operator.
   */
  MEee2gZ2bb & operator=(const MEee2gZ2bb &);

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
