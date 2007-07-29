// -*- C++ -*-
#ifndef HERWIG_SusyBase_H
#define HERWIG_SusyBase_H
//
// This is the declaration of the SusyBase class.
//

#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Models/Susy/MixingMatrix.h"
#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/SSSVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/VVSSVertex.h"
#include "SusyBase.fh"

namespace Herwig {
using namespace ThePEG;

  /*@name Some convenient typedefs. */
  //@{
  /** 
   * Map to hold key, parameter pairs. 
   */
  typedef map<long, double> ParamMap;
  //@}

/** \ingroup Models
   * This class is designed to be a base class for SUSY models. There is
   * an interface to set the name of the spectrum file to read in 
   * the necessary parameters for a SUSY model.
   *
   * @see \ref SusyBaseInterfaces "The interfaces"
   * defined for SusyBase.
   * @see StandardModel
   */  

class SusyBase: public StandardModel {
  
public:

  /**
   * The default constructor.
   */
  inline SusyBase();

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
  inline Energy muParameter() const;

  /**
   * Value of soft breaking mass for the bino
   */
  inline Energy softMOne() const;

  /**
   * Value of soft breaking mass for the wino
   */
  inline Energy softMTwo() const;

  /**
   * Value of soft breaking mass for the gluino
   */
  inline Energy softMThree() const;

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

  /**
   * Tell EventGenerator::init() that this is an object that wants to
   * create new interfaced objects in the read step
   */
  virtual bool preInitialize() const;

  /**@name Functions to access specific vertices.*/
  //@{
  /**
   * Pointer to the MSSM fermion-antifermion-higgs vertex 
   */
  virtual inline tFFSVertexPtr vertexFFH() const;
  
  /**
   * Pointer to the MSSM double gauge boson-higgs vertex 
   */
  virtual inline tVVSVertexPtr vertexWWH() const;
  
  /**
   * Pointer to the MSSM effective higgs-gluon-gluon vertex
   */
  virtual inline tGeneralSVVVertexPtr vertexHGG() const;
  //@}


private:
  
  /**
   * Read the SLHA file with the name set by an interface.
   */
  void readSLHA() throw(InitException);
  
private:
  
  /**@name Functions to help file read-in. */
  //@{
  /**
   * Read block from LHA file
   * @param ifs input stream containg data
   * @param name The name of the block
   */
  void readBlock(ifstream & ifs,string name) throw(InitException);

  /**
   * Function to read mixing matrix from LHA file
   * @param ifs input stream containg data
   * @param row Number of rows
   * @param col Number of columns
   */
  const MixingVector readMatrix(ifstream & ifs, unsigned int & row,
				unsigned int & col) throw(InitException);

  /**
   * Read decaymodes from LHA file
   * @param ifs input stream containg data
   * @param decay string containing name of parent and value of total width
   */
  void readDecay(ifstream & ifs, string decay) const throw(InitException);

  /**
   * Create a DecayMode object in the repository
   * @param tag string containing first portion of tag including '->'
   * @param products Decay products
   * @param brat Branching ratio of this mode 
   */
  void createDecayMode(string tag, vector<long> products,
		       double brat) const;

protected:

  /**
   *  Create the mixing matrices for the model
   */
  virtual void createMixingMatrices();

  /**
   *  Extract the parameters from the input blocks
   */
  virtual void extractParameters(bool checkmodel=true);

  /**
   * Create a object MixingMatrix in the repository
   * @param Pointer to the mixing matrix
   * @param name Name of the mixing matrix, i.e. nmix, umix...
   * @param values Value of each entry in the matrix
   * @param size The size of the matrix
   */
  void createMixingMatrix(MixingMatrixPtr & matrix, string name, 
			  const MixingVector & values,
			  MatrixSize size);

  /**
   * Reset masses in the repository to values read from LHA file.
   */
  void resetRepositoryMasses();

  /**
   * Adjust row of Mixing Matrix if a negative mass occurs in LHA file
   * @param id The PDG code of the particle with a negative mass
   */
  virtual void adjustMixingMatrix(long id);
  //@}

  /**
   *  Access to the mixings and parameters for the inheriting classes
   */
  //@{
  /**
   *  Parameter blocks
   */
  inline const map<string,ParamMap> & parameters() const;

  /**
   *  Mixing blocks
   */
  inline const map<string,pair<MatrixSize,MixingVector> > & mixings() const;
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
   *  Parameter blocks
   */
  map<string,ParamMap> _parameters;

  /**
   *  Mixing blocks
   */
  map<string,pair<MatrixSize, MixingVector> > _mixings;

  /**
   *  \f$\tan\beta\f$
   */
  double _tanbeta;

  /**
   *  \f$\mu\f$
   */
  Energy _mu;

  /**
   * The bilinear breaking mass term for the bino
   */
  Energy theMone;
  
  /**
   * The bilinear breaking mass term for the wino
   */
  Energy theMtwo;

  /**
   * The bilinear breaking mass term for the gluinos
   */
  Energy theMthree;  
  //@}

  /**
   * The name of the SLHA file
   */
  string theSLHAName;
  
  /**
   *  Neutralino and Chargino mixing matrices
   */
  //@{
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
   * Pointer to the Higgs-sfermion-sfermion vertex
   */
  SSSVertexPtr theHSFSFVertex;

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

  /**
   * Pointer to the vertex fermion-antifermion-higgs vertex
   */
  FFSVertexPtr theSSFFHVertex;

  /**
   * Pointer to the vertex gaugino-gaugino-higgs vertex
   */
  FFSVertexPtr theGOGOHVertex;
  
  /**
   * Pointer to the vertex for a pair of gauge bosons and higgs
   */
  VVSVertexPtr theSSWWHVertex;
  
  /**
   * Pointer to the vertex for a gauge boson and higgs
   */
  VSSVertexPtr theWHHVertex;

  /**
   * Pointer to triple higgs vertex
   */
  SSSVertexPtr theHHHVertex;
  
  /**
   * The effective coupling of the higgs to a pai of gluons in the MSSM
   */
  GeneralSVVVertexPtr theSSHGGVertex;
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
  static string className() { return "Herwig::SusyBase"; }
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
