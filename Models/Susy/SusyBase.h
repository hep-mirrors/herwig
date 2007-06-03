// -*- C++ -*-
#ifndef HERWIG_SusyBase_H
#define HERWIG_SusyBase_H
//
// This is the declaration of the SusyBase class.
//

#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Models/Susy/MixingMatrix.h"
#include "Herwig++/Helicity/Vertex/Scalar/VSSVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/SSSVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/VVSSVertex.h"
#include "SusyBase.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Models
   * This class is designed to be a base class for SUSY models. It implements
   * the readSetup function from InterfacedBase to open a LesHouches file and
   * contains other member functions to read in the data.
   *
   * @see StandardModel
   */  

class SusyBase: public StandardModel {

public:
  
  /*@name Some convenient typedefs. */
  //@{
  /** Map to hold key, parameter pairs. */
  typedef map<long, double> ParamMap;
  //@}
  
public:

  /**
   * The default constructor.
   */
  inline SusyBase();

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
  
  /**
   * Function to read information from setup file
   * @param is istream object to be read from 
   */
  virtual void readSetup(istream & is) throw(SetupException);

public:

  /** @name Access functions. */
  //@{
  /**
   * Value of \f$\tan\beta\f$.
   */
  inline double tanBeta() const;

  /**
   * Value of \f$\mu\f$ parameter.
   */
  inline double muParameter() const;

  /**
   * Value of Higgs mixing angle \f$\alpha\f$.
   */
  inline double higgsMixingAngle() const;

  /**
   * Value of up-type trilinear couplings
   */
  inline const vector<Complex> & upTypeTrilinear() const;

  /**
   * Value of down-type trilinear couplings
   */
  inline const vector<Complex> & downTypeTrilinear() const;

  /**
   * Value of lepton trilinear couplings
   */
  inline const vector<Complex> & leptonTypeTrilinear() const;

  /**
   * The neutralino mixing matrix
   */
  inline const MixingMatrixPtr & neutralinoMix() const;

  /**
   * The U-type chargino mixing matrix
   */
  inline const MixingMatrixPtr & charginoUMix() const;

  /**
   * The V-type chargino mixing matrix
   */
  inline const MixingMatrixPtr & charginoVMix() const;

  /**
   * The stop mixing matrix
   */
  inline const MixingMatrixPtr & stopMix() const;

  /**
   * The sbottom chargino mixing matrix
   */
  inline const MixingMatrixPtr & sbottomMix() const;

  /**
   * The stau mixing matrix
   */
  inline const MixingMatrixPtr & stauMix() const;
  //@}

private:
  
  /**@name Functions to help file read-in. */
  //@{
  /**
   * Read block from LHA file
   * @param ifs input stream containg data
   */
  ParamMap readBlock(ifstream & ifs) const throw(SetupException);

  /**
   * Function to read mixing matrix from LHA file
   * @param ifs input stream containg data
   * @param size The size of the matrix read-in
   */
  const vector<Complex> readMatrix(ifstream & ifs, unsigned int & size)
    throw(SetupException);

  /**
   * Read decaymodes from LHA file
   * @param ifs input stream containg data
   * @param decay string containing name of parent and value of total width
   */
  void readDecay(ifstream & ifs, string decay) const throw(SetupException);

  /**
   * Create a DecayMode object in the repository
   * @param tag string containing first portion of tag including '->'
   * @param products Decay products
   * @param brat Branching ratio of this mode 
   */
  void createDecayMode(string tag, vector<long> products,
		       double brat) const;

  /**
   * Create a object MixingMatrix in the repository
   * @param name Name of the mixing matrix, i.e. nmix, umix...
   * @param values Value of each entry in the matrix
   * @param size The size of the matrix
   */
  void createMixingMatrix(string name, vector<Complex> & values,
			  unsigned int size);

  /**
   * Reset masses in the repository to values read from LHA file.
   */
  void resetRepositoryMasses();

  /**
   * Adjust row of Mixing Matrix if a negative mass occurs in LHA file
   * @param id The PDG code of the particle with a negative mass
   */
  void adjustMixingMatrix(long id);
  //@}
  
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
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SusyBase> initSusyBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SusyBase & operator=(const SusyBase &);

private:

  /*
   * Storage of the parameters.
   */
  //@{
  /**
   *  The minpar parameter block
   */
   ParamMap theMinPar;

  /**
   * The hmix parameter block
   */
   ParamMap theHMix;
  //@}

  /**
   * Storage of the masses (gets cleared once they have been set in the 
   * repository) .
   */
   ParamMap theMasses;

  /**
   * Store pointers to the gaugino mixing matrices
   */
  //@{
  /**
   * The neutralino mixing matrix
   */
  MixingMatrixPtr theNMix; 

  /**
   * The \f$U\f$ mixing matrix for the charginos
   */
  MixingMatrixPtr theUMix; 

  /**
   * The \f$V\f$ mixing matrix for the charginos
   */
  MixingMatrixPtr theVMix; 
  //@}

  /**
   * Store pointers to the squark/slepton mixing matrices
   */
  //@{
  /**
   *  The \f$\tilde{t}\f$ mixing matrix
   */
  MixingMatrixPtr theStopMix;

  /**
   *  The \f$\tilde{b}\f$ mixing matrix
   */
  MixingMatrixPtr theSbotMix;

  /**
   *  The \f$\tilde{\tau}\f$ mixing matrix
   */
  MixingMatrixPtr theStauMix;
  //@}
  /**
   * Value of higgs mixing angle.
   */
  double theAlpha;

  /**
   * Trilinear couplings stored as vector of complex numbers to make use
   * of routine already available to read complex matrices
   */
  //@{
  /**
   *  For the up type squarks
   */
  vector<Complex> theAu;

  /**
   *  For the down type squarks
   */
  vector<Complex> theAd;

  /**
   *  For the charged sleptons
   */
  vector<Complex> theAe;
  //@}

  /**@name Vertex pointers. */
  //@{
  /**
   * Pointer to the Z-sfermion-sfermion vertex
   **/
  VSSVertexPtr theZSFSFVertex;
  
  /**
   * Pointer to the photon-sfermion-sfermion vertex
   */
  VSSVertexPtr thePSFSFVertex;
  
  /**
   * Pointer to the W-sfermion-sfermion vertex
   */
  VSSVertexPtr theWSFSFVertex;
  
  /**
   * Pointer to the neutralino-fermion-sfermion vertex
   */
  FFSVertexPtr theNFSFVertex;
  
  /**
   * Pointer to the gluino-fermion-sfermion coupling
   */
  FFSVertexPtr theGFSFVertex;

  /**
   * Pointer to the \f$h^0\f$-sfermion-sfermion vertex
   */
  SSSVertexPtr theH1SFSFVertex;
  
  /**
   * Pointer to the \f$H^0\f$-sfermion-sfermion vertex
   */
  SSSVertexPtr theH2SFSFVertex;

  /**
   * Pointer to the \f$H^0\f$-sfermion-sfermion vertex
   */
  SSSVertexPtr theH3SFSFVertex;

  /**
   * Pointer to the \f$\tilde{\chi}^+\f$-fermion-sfermion vertex
   */
  FFSVertexPtr theCFSFVertex;

  /**
   * Pointer to the gluon-sfermion-sfermion vertex
   */
  VSSVertexPtr theGSFSFVertex;

  /**
   * Pointer to the gluon-gluon-squark-squark vertex;
   */
  VVSSVertexPtr theGGSQSQVertex;

  /**
   * Pointer to the gluon-gluino-gluino vertex
   */
  FFVVertexPtr theGSGSGVertex; 

  /**
   * Pointer to the neutralino-neutralino-Z vertex
   */
  FFVVertexPtr theNNZVertex;

  /**
   * Pointer to the  vertex chargino-chargino-photon vertex
   */
  FFVVertexPtr theCCPVertex;

  /**
   * Pointer to the  vertex chargino-chargino-Z vertex
   */
  FFVVertexPtr theCCZVertex;
  
  /**
   * Pointer to the  vertex chargino-neutralino-Z vertex
   */
  FFVVertexPtr theCNWVertex;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SusyBase. */
template <>
struct BaseClassTrait<Herwig::SusyBase,1> {
  /** Typedef of the first base class of SusyBase. */
  typedef Herwig::StandardModel NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SusyBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SusyBase>
  : public ClassTraitsBase<Herwig::SusyBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SusyBase"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SusyBase is implemented. It may also include several, space-separated,
   * libraries if the class SusyBase depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#include "SusyBase.icc"

#endif /* HERWIG_SusyBase_H */
