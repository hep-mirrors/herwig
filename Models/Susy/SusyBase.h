// -*- C++ -*-
//
// SusyBase.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SusyBase_H
#define HERWIG_SusyBase_H
//
// This is the declaration of the SusyBase class.
//

#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Models/Susy/MixingMatrix.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractSSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSSVertex.h"
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
  SusyBase();

public:

  /** @name Access functions. */
  //@{
  /**
   * Value of \f$\tan\beta\f$.
   */
  inline double tanBeta() const { return _tanbeta; }

  /**
   * Value of \f$\mu\f$ parameter.
   */
  inline Energy muParameter() const { return _mu; }

  /**
   * Value of soft breaking mass for the bino
   */
  inline Energy softMOne() const { return theMone; }

  /**
   * Value of soft breaking mass for the wino
   */
  inline Energy softMTwo() const { return theMtwo; }

  /**
   * Value of soft breaking mass for the gluino
   */
  inline Energy softMThree() const { return theMthree; }

  /**
   * The neutralino mixing matrix
   */
  inline const MixingMatrixPtr & neutralinoMix() const { 
    return theNMix;
  }

  /**
   * The U-type chargino mixing matrix
   */
  inline const MixingMatrixPtr & charginoUMix() const {
    return theUMix;
  }

  /**
   * The V-type chargino mixing matrix
   */
  inline const MixingMatrixPtr & charginoVMix() const {
    return theVMix;
  }
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

  /**@name Functions to access specific vertices.*/
  //@{
  /**
   * Pointer to the electroweak gauge boson Higgs-Higgs vertex.
   */
  virtual inline tAbstractVSSVertexPtr vertexWHH() const {
    return theWHHVertex;
  }

  /**
   * Pointer to the higgs coupling to a pair of gauginos
   */
  virtual inline tAbstractFFSVertexPtr vertexGOGOH() const {
    return theGOGOHVertex;
  }

  /**
   * Pointer to the triple higgs vertex
   */
  virtual inline tAbstractSSSVertexPtr vertexHHH() const {
    return theHHHVertex;
  }

  /**
   * Pointer to higgs-sfermion-sfermion vertex 
   */
  virtual inline tAbstractSSSVertexPtr vertexHSS() const {
    return theHSFSFVertex;
  }
  //@}

protected:

  /**
   * Function to read information from a setup file.
   * @param is istream object to read file.
   */
  virtual void readSetup(istream & is);

private:
  
  /**@name Functions to help file read-in. */
  //@{
  /**
   * Read block from LHA file
   * @param ifs input stream containg data
   * @param name The name of the block
   */
  void readBlock(ifstream & ifs,string name);

  /**
   * Function to read mixing matrix from LHA file
   * @param ifs input stream containg data
   * @param row Number of rows
   * @param col Number of columns
   */
  const MixingVector readMatrix(ifstream & ifs, unsigned int & row,
				unsigned int & col);

  /**
   * Read decaymodes from LHA file
   * @param ifs input stream containg data
   * @param decay string containing name of parent and value of total width
   */
  void readDecay(ifstream & ifs, string decay) const;

  /**
   * Create a DecayMode object in the repository
   * @param tag The tag identifying the decay mode including the prefix
   * 'decaymode'
   * @param brat Branching ratio of this mode 
   */
  void createDecayMode(string tag, double brat) const;

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
   * @param matrix Pointer to the mixing matrix
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
  inline const map<string,ParamMap> & parameters() const {
    return _parameters;
  }

  /**
   *  Mixing blocks
   */
  inline const map<string,pair<MatrixSize,MixingVector> > & mixings() const {
    return _mixings;
  }
  //@}

  /**
   * Reset neutralino mixing matrix
   */
  inline void neutralinoMix(MixingMatrixPtr nm) { theNMix = nm; }

  /**
   * Reset the U-type chargino mixing matrix
   */
  inline void charginoUMix(MixingMatrixPtr um) { theUMix = um; }

  /**
   *  Reset the V-type chargino mixing matrix
   */
  inline void charginoVMix(MixingMatrixPtr vm) { theVMix = vm; }
  
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

  /**
   *  Whether or not the SLHA fiel has been read
   */
  bool _readFile;

  /**
   * Whether or not to replace the top decay modes with those from
   * the SLHA files
   */
  bool _topModesFromFile;

  /**
   *  Tolerance for branching ratios
   */
  double _tolerance;

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
   * Pointer to the gauge boson sfermion-sfermion vertex
   */
  AbstractVSSVertexPtr theWSFSFVertex;
  
  /**
   * Pointer to the neutralino-fermion-sfermion vertex
   */
  AbstractFFSVertexPtr theNFSFVertex;
  
  /**
   * Pointer to the gluino-fermion-sfermion coupling
   */
  AbstractFFSVertexPtr theGFSFVertex;

  /**
   * Pointer to the Higgs-sfermion-sfermion vertex
   */
  AbstractSSSVertexPtr theHSFSFVertex;

  /**
   * Pointer to the \f$\tilde{\chi}^+\f$-fermion-sfermion vertex
   */
  AbstractFFSVertexPtr theCFSFVertex;

  /**
   * Pointer to the gluon-sfermion-sfermion vertex
   */
  AbstractVSSVertexPtr theGSFSFVertex;

  /**
   * Pointer to the gluon-gluon-squark-squark vertex;
   */
  AbstractVVSSVertexPtr theGGSQSQVertex;

  /**
   * Pointer to the gluon-gluino-gluino vertex
   */
  AbstractFFVVertexPtr theGSGSGVertex; 

  /**
   * Pointer to the gluino-neutralino-gluon vertex
   */
  AbstractFFVVertexPtr theGNGVertex;

  /**
   * Pointer to the neutralino-neutralino-Z vertex
   */
  AbstractFFVVertexPtr theNNZVertex;

  /**
   * Pointer to the neutralino-neutralino-photon vertex
   */
  AbstractFFVVertexPtr theNNPVertex;

  /**
   * Pointer to the  vertex chargino-chargino-Z vertex
   */
  AbstractFFVVertexPtr theCCZVertex;
  
  /**
   * Pointer to the  vertex chargino-neutralino-Z vertex
   */
  AbstractFFVVertexPtr theCNWVertex;

  /**
   * Pointer to the vertex gaugino-gaugino-higgs vertex
   */
  AbstractFFSVertexPtr theGOGOHVertex;
  
  /**
   * Pointer to the vertex for a gauge boson and higgs
   */
  AbstractVSSVertexPtr theWHHVertex;

  /**
   * Pointer to triple higgs vertex
   */
  AbstractSSSVertexPtr theHHHVertex;
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


#endif /* HERWIG_SusyBase_H */
