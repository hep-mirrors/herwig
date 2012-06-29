// -*- C++ -*-
#ifndef HERWIG_HardProcessConstructor_H
#define HERWIG_HardProcessConstructor_H
//
// This is the declaration of the HardProcessConstructor class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "HPDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "Herwig++/MatrixElement/General/GeneralHardME.h"
#include "HardProcessConstructor.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HardProcessConstructor class.
 *
 * @see \ref HardProcessConstructorInterfaces "The interfaces"
 * defined for HardProcessConstructor.
 */
class HardProcessConstructor: public Interfaced {

public:

  /** Vector of HPDiagrams. */
  typedef vector<HPDiagram> HPDVector;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  HardProcessConstructor() : debug_(false) {}
  //@}

  /**
   * The main function to create diagrams etc for the processes
   */
  virtual void constructDiagrams() = 0;

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

  /** Functions to set up colour flows and matrix elements. */
  //@{
  /**
   * Determine whether the ordering of the outgoing states is the same
   * as the ordering in the matrix elements
   * @param diag The diagram to question
   */
  void fixFSOrder(HPDiagram & diag);

  /**
   * Assign a diagram to the appropriate colour flow(s).
   * @param diag The diagram to assign
   */
  void assignToCF(HPDiagram & diag);

  /**
   * Assign a $s$-channel diagram to the appropriate colour flow(s).
   * @param diag The diagram to assign
   */
  void sChannelCF(HPDiagram & diag);

  /**
   * Assign a $t$-channel diagram to the appropriate colour flow(s).
   * @param diag The diagram to assign
   */
  void tChannelCF(HPDiagram & diag);

  /**
   * Assign a $u$-channel diagram to the appropriate colour flow(s).
   * @param diag The diagram to assign
   */
  void uChannelCF(HPDiagram & diag);

  /**
   * Assign a $u$-channel diagram to the appropriate colour flow(s).
   * @param diag The diagram to assign
   */
  void fourPointCF(HPDiagram & diag);
  //@}

  /**
   * Pointer to the model being used
   */
  tHwSMPtr model() const {return model_;}
  
  /**
   * Pointer to the sub process handler
   */
  tSubHdlPtr subProcess() const {return subProcess_;}

  /**
   * Whether to print the debug information with the matrix 
   * element. This is here solely so it can be passed to 
   * a matrix element that is created here.
   */
  bool debug() const {return debug_;}

  /**
   * Get the correct colour factor matrix.
   * @param extpart Vector of external ParticleData pointers
   */
  GeneralHardME::ColourStructure colourFlow(const tcPDVector & extpart) const;

  /**
   * Search for a diagram that has already been created
   * @param diagram The diagram to search for
   * @param group The group of diagrams to search through 
   */
  bool duplicate(const HPDiagram & diagram, 
		 const HPDVector & group) const;

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<HardProcessConstructor> initHardProcessConstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HardProcessConstructor & operator=(const HardProcessConstructor &);

private:

  /**
   * Pointer to the model being used
   */
  tHwSMPtr model_;
  
  /**
   * Pointer to the sub process handler
   */
   tSubHdlPtr subProcess_;

  /**
   * Whether to print the debug information with the matrix 
   * element. This is here solely so it can be passed to 
   * a matrix element that is created here.
   */
  bool debug_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HardProcessConstructor. */
template <>
struct BaseClassTrait<Herwig::HardProcessConstructor,1> {
  /** Typedef of the first base class of HardProcessConstructor. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HardProcessConstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HardProcessConstructor>
  : public ClassTraitsBase<Herwig::HardProcessConstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HardProcessConstructor"; }
};

/** @endcond */

}

#endif /* HERWIG_HardProcessConstructor_H */
