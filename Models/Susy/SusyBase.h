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
#include "MixingMatrix.h"
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
  double tanBeta() const { return tanBeta_; }

  /**
   * Value of \f$\mu\f$ parameter.
   */
  Energy muParameter() const { return mu_; }

  /**
   * The neutralino mixing matrix
   */
  const MixingMatrixPtr & neutralinoMix() const { 
    return NMix_;
  }

  /**
   * The U-type chargino mixing matrix
   */
  const MixingMatrixPtr & charginoUMix() const {
    return UMix_;
  }

  /**
   * The V-type chargino mixing matrix
   */
  const MixingMatrixPtr & charginoVMix() const {
    return VMix_;
  }

  /**
   *  The phase for gluino vertices
   */
  const Complex & gluinoPhase() const {return gluinoPhase_;}
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

public:

  /**
   *  Soft breaking parameters
   */
  //@{
  /**
   * The bilinear breaking mass term for the bino
   */
  const Energy & M1() const {return M1_;}
  
  /**
   * The bilinear breaking mass term for the wino
   */
  const Energy & M2() const {return M2_;}

  /**
   * The bilinear breaking mass term for the gluinos
   */
  const Energy & M3() const {return M3_;}

  /**
   *  The soft breaking mass squared for \f$H_1\f$
   */
  const Energy2 & Mh12() const {return mH12_;}

  /**
   *  The soft breaking mass squared for \f$H_2\f$
   */
  const Energy2 & Mh22() const {return mH22_;}

  /**
   *  Soft breaking mass for the first generation lepton doublet
   */
  const Energy & MeL() const {return meL_;}

  /**
   *  Soft breaking mass for the second generation lepton doublet
   */
  const Energy & MmuL() const {return mmuL_;}

  /**
   *  Soft breaking mass for the third generation lepton doublet
   */
  const Energy & MtauL() const {return mtauL_;} 

  /**
   *  Soft breaking mass for the first generation lepton singlet
   */
  const Energy & MeR() const {return meR_;}

  /**
   *  Soft breaking mass for the second generation lepton singlet
   */
  const Energy & MmuR() const {return mmuR_;}

  /**
   *  Soft breaking mass for the third generation lepton singlet
   */
  const Energy & MtauR() const {return mtauR_;} 

  /**
   *  Soft breaking mass for the first generation quark doublet
   */
  const Energy & Mq1L() const {return mq1L_;}

  /**
   *  Soft breaking mass for the second generation quark doublet
   */
  const Energy & Mq2L() const {return mq2L_;}

  /**
   *  Soft breaking mass for the third generation quark doublet
   */
  const Energy & Mq3L() const {return mq3L_;}

  /**
   *  Soft breaking mass for the down singlet
   */ 
  const Energy & MdR() const {return mdR_;}

  /**
   *  Soft breaking mass for the up singlet
   */ 
  const Energy & MuR() const {return muR_;}

  /**
   *  Soft breaking mass for the strange singlet
   */ 
  const Energy & MsR() const {return msR_;}

  /**
   *  Soft breaking mass for the charm singlet
   */ 
  const Energy & McR() const {return mcR_;}

  /**
   *  Soft breaking mass for the bottom singlet
   */ 
  const Energy & MbR() const {return mbR_;}

  /**
   *  Soft breaking mass for the top singlet
   */ 
  const Energy & MtR() const {return mtR_;}
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
  void readBlock(ifstream & ifs,string name,string line);

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
  const map<string,ParamMap> & parameters() const {
    return parameters_;
  }

  /**
   *  Mixing blocks
   */
  const map<string,pair<MatrixSize,MixingVector> > & mixings() const {
    return mixings_;
  }
  //@}

  /**
   * Reset neutralino mixing matrix
   */
  void neutralinoMix(MixingMatrixPtr nm) { NMix_ = nm; }

  /**
   * Reset the U-type chargino mixing matrix
   */
  void charginoUMix(MixingMatrixPtr um) { UMix_ = um; }

  /**
   *  Reset the V-type chargino mixing matrix
   */
  void charginoVMix(MixingMatrixPtr vm) { VMix_ = vm; }
  
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
  bool readFile_;

  /**
   * Whether or not to replace the top decay modes with those from
   * the SLHA files
   */
  bool topModesFromFile_;

  /**
   *  Tolerance for branching ratios
   */
  double tolerance_;

  /*
   * Storage of the parameters.
   */
  //@{
  /**
   *  Parameter blocks
   */
  map<string,ParamMap> parameters_;

  /**
   *  Mixing blocks
   */
  map<string,pair<MatrixSize, MixingVector> > mixings_;

  /**
   *  \f$\tan\beta\f$
   */
  double tanBeta_;

  /**
   *  \f$\mu\f$
   */
  Energy mu_;
  //@}

  /**
   *  Soft breaking parameters
   */
  //@{
  /**
   * The bilinear breaking mass term for the bino
   */
  Energy M1_;
  
  /**
   * The bilinear breaking mass term for the wino
   */
  Energy M2_;

  /**
   * The bilinear breaking mass term for the gluinos
   */
  Energy M3_;

  /**
   *  The soft breaking mass squared for \f$H_1\f$
   */
  Energy2 mH12_;

  /**
   *  The soft breaking mass squared for \f$H_2\f$
   */
  Energy2 mH22_;

  /**
   *  Soft breaking mass for the first generation lepton doublet
   */
  Energy meL_;

  /**
   *  Soft breaking mass for the second generation lepton doublet
   */
  Energy mmuL_;

  /**
   *  Soft breaking mass for the third generation lepton doublet
   */
  Energy mtauL_; 

  /**
   *  Soft breaking mass for the first generation lepton singlet
   */
  Energy meR_;

  /**
   *  Soft breaking mass for the second generation lepton singlet
   */
  Energy mmuR_;

  /**
   *  Soft breaking mass for the third generation lepton singlet
   */
  Energy mtauR_; 

  /**
   *  Soft breaking mass for the first generation quark doublet
   */
  Energy mq1L_;

  /**
   *  Soft breaking mass for the second generation quark doublet
   */
  Energy mq2L_;

  /**
   *  Soft breaking mass for the third generation quark doublet
   */
  Energy mq3L_;

  /**
   *  Soft breaking mass for the down singlet
   */ 
  Energy mdR_;

  /**
   *  Soft breaking mass for the up singlet
   */ 
  Energy muR_;

  /**
   *  Soft breaking mass for the strange singlet
   */ 
  Energy msR_;

  /**
   *  Soft breaking mass for the charm singlet
   */ 
  Energy mcR_;

  /**
   *  Soft breaking mass for the bottom singlet
   */ 
  Energy mbR_;

  /**
   *  Soft breaking mass for the top singlet
   */ 
  Energy mtR_;
  //@}

  /**
   *  Phase for the gluino
   */
  Complex gluinoPhase_;

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
  MixingMatrixPtr NMix_; 

  /**
   * The \f$U\f$ mixing matrix for the charginos
   */
  MixingMatrixPtr UMix_; 

  /**
   * The \f$V\f$ mixing matrix for the charginos
   */
  MixingMatrixPtr VMix_; 
  //@}

  /**@name Vertex pointers. */
  //@{
  /**
   * Pointer to the gauge boson sfermion-sfermion vertex
   */
  AbstractVSSVertexPtr WSFSFVertex_;
  
  /**
   * Pointer to the neutralino-fermion-sfermion vertex
   */
  AbstractFFSVertexPtr NFSFVertex_;
  
  /**
   * Pointer to the gluino-fermion-sfermion coupling
   */
  AbstractFFSVertexPtr GFSFVertex_;

  /**
   * Pointer to the Higgs-sfermion-sfermion vertex
   */
  AbstractSSSVertexPtr HSFSFVertex_;

  /**
   * Pointer to the \f$\tilde{\chi}^+\f$-fermion-sfermion vertex
   */
  AbstractFFSVertexPtr CFSFVertex_;

  /**
   * Pointer to the gluon-sfermion-sfermion vertex
   */
  AbstractVSSVertexPtr GSFSFVertex_;

  /**
   * Pointer to the gluon-gluon-squark-squark vertex;
   */
  AbstractVVSSVertexPtr GGSQSQVertex_;

  /**
   * Pointer to the gluon-gluino-gluino vertex
   */
  AbstractFFVVertexPtr GSGSGVertex_; 

  /**
   * Pointer to the gluino-neutralino-gluon vertex
   */
  AbstractFFVVertexPtr GNGVertex_;

  /**
   * Pointer to the neutralino-neutralino-Z vertex
   */
  AbstractFFVVertexPtr NNZVertex_;

  /**
   * Pointer to the neutralino-neutralino-photon vertex
   */
  AbstractFFVVertexPtr NNPVertex_;

  /**
   * Pointer to the  vertex chargino-chargino-Z vertex
   */
  AbstractFFVVertexPtr CCZVertex_;
  
  /**
   * Pointer to the  vertex chargino-neutralino-Z vertex
   */
  AbstractFFVVertexPtr CNWVertex_;

  /**
   * Pointer to the vertex gaugino-gaugino-higgs vertex
   */
  AbstractFFSVertexPtr GOGOHVertex_;
  
  /**
   * Pointer to the vertex for a gauge boson and higgs
   */
  AbstractVSSVertexPtr WHHVertex_;

  /**
   *  Pointer to the vertex for flavour changing stop decay
   */
  AbstractFFSVertexPtr NCTVertex_;
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
