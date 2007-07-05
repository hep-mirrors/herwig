// -*- C++ -*-
#ifndef HERWIG_TwoBodyDecayConstructor_H
#define HERWIG_TwoBodyDecayConstructor_H
//
// This is the declaration of the TwoBodyDecayConstructor class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/Decay/General/GeneralTwoBodyDecayer.fh"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "TwoBodyDecayConstructor.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::VertexType;
using Helicity::VertexBasePtr;

/**
 * The TwoBodyDecayConstructor class inherits from the dummy base class
 * NBodyDecayConstructorBase and implements the necessary functions in
 * order to create the 2 body decaymodes for a given set of vertices
 * stored in a Model class.
 *
 * @see NBodyDecayConstructor
 **/
class TwoBodyDecayConstructor: public NBodyDecayConstructorBase {

public:

  /**
   * The default constructor.
   */
  inline TwoBodyDecayConstructor();

  /**
   * Function used to determine allowed decaymodes
   *@param part vector of ParticleData pointers containing particles in model
   */
  virtual void DecayList(const PDVector & part);
  
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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<TwoBodyDecayConstructor> initTwoBodyDecayConstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TwoBodyDecayConstructor & operator=(const TwoBodyDecayConstructor &);

private:
  
  /** @name Functions to create decayers and decaymodes. */
  //@{
  /**
   * Function to create decays
   * @param inpart Incoming particle 
   * @param vert The vertex to create decays for
   * @param ilist Which list to search
   * @param iv Row number in _theExistingDecayers member
   * @return vector of ParticleData ptrs
   */
  PDVector createModes(tPDPtr inpart, VertexBasePtr vert,
		       unsigned int ilist,
		       unsigned int iv);

  /**
   * Function to create decayer for specific vertex
   * @param vert Pointer to vertex 
   * @param icol Integer referring to the colmun in _theExistingDecayers
   * @param ivert Integer referring to the row in _theExistingDecayers
   * member variable
   */
  void createDecayer(VertexBasePtr vert, unsigned int icol,
		     unsigned int ivert);

  /**
   * Create decay mode(s) from given part and decay modes
   * @param inpart pointer to incoming particle
   * @param decays list of allowed interactions
   * @param decayer The decayer responsible for this decay
   */
  void createDecayMode(tPDPtr inpart,
		       const PDVector & decays,
		       GeneralTwoBodyDecayerPtr decayer);

  /**
   * Set the branching ratio of this mode. This requires 
   * calculating a new width for the decaying particle and reweighting
   * the current branching fractions.
   * @param dm The decaymode for which to set the branching ratio
   * @param pwidth The calculated width of the mode
   */
    void setBranchingRatio(tDMPtr dm, Energy pwidth);

  /**
   * Set the interfaces on the decayers to initialise them
   * @param name Fullname of the decayer in the EventGenerator
   * including the path
   */
  void initializeDecayers(string name) const;
  //@}

private:
  
  /**
   * Existing decayers
   */
   vector<vector<GeneralTwoBodyDecayerPtr> > _theExistingDecayers;

  /**
   * Model Pointer
   */
  Ptr<Herwig::StandardModel>::pointer _theModel;

  /**
   * Whether to initialize decayers or not
   */
  bool _init;
  
  /**
   * Number of iterations if initializing (default 1)
   */
  int _iteration;

  /**
   * Number of points to do in initialization
   */
  int _points;
};
  
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TwoBodyDecayConstructor. */
template <>
struct BaseClassTrait<Herwig::TwoBodyDecayConstructor,1> {
  /** Typedef of the first base class of TwoBodyDecayConstructor. */
  typedef Herwig::NBodyDecayConstructorBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TwoBodyDecayConstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TwoBodyDecayConstructor>
  : public ClassTraitsBase<Herwig::TwoBodyDecayConstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::TwoBodyDecayConstructor"; }
  /** Return the name of the shared library be loaded to get
   *  access to the TwoBodyDecayConstructor class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwModelGenerator.so"; }
};

/** @endcond */

}

#include "TwoBodyDecayConstructor.icc"

#endif /* HERWIG_TwoBodyDecayConstructor_H */
