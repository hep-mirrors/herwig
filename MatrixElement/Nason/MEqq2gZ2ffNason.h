// -*- C++ -*-
#ifndef HERWIG_MEqq2gZ2ffNason_H
#define HERWIG_MEqq2gZ2ffNason_H
//
// This is the declaration of the MEqq2gZ2ffNason class.
//

#include "ThePEG/MatrixElement/ME2to2Base.h"
#include "Herwig++/MatrixElement/General/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Utilities/Statistic.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "MEqq2gZ2ffNason.fh"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEqq2gZ2ffNason class implements the products of Standard Model
 * fermion antifermion pairs via the \f$Z^0\f$ resonance including
 * photon interference terms.
 *
 * @see \ref MEqq2gZ2ffNasonInterfaces "The interfaces"
 * defined for MEqq2gZ2ffNason.
 */
class MEqq2gZ2ffNason: public ME2to2Base {

public:

  /**
   * The default constructor.
   */
  inline MEqq2gZ2ffNason();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  inline virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   *  Construct the vertex of spin correlations.
   */
  virtual void constructVertex(tSubProPtr);
  /**
   * The number of internal degreed of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given 'nDim()' uniform
   * random numbers in the interval ]0,1[. To help the phase space
   * generator, the 'dSigHatDR()' should be a smooth function of these
   * numbers, although this is not strictly necessary. Return
   * false if the chosen points failed the kinematical cuts.
   */
  virtual bool generateKinematics(const double * r);
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

protected:

  /**
   * Matrix element for \f$q\bar{q}\to \gamma/Z \to f\bar{f}\f$.
   * @param fin  Spinors for incoming quark
   * @param ain  Spinors for incoming antiquark
   * @param fout Spinors for incoming quark
   * @param aout Spinors for incoming antiquark
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double qqbarME(vector<SpinorWaveFunction>    & fin ,
		 vector<SpinorBarWaveFunction> & ain ,
		 vector<SpinorBarWaveFunction> & fout,
		 vector<SpinorWaveFunction>    & aout,
		 bool me) const;

  /**
   *  Calculate the correction weight
   */
  double NLOweight() const;
  
  mutable double _max_wgt;
 
  mutable double  _xb_a,_xb_b,_alphaS,_TF,_CF;
  mutable Energy2 _mll2,_mu2;
  mutable tcPDPtr _parton_a,_parton_b,_gluon;
  mutable Ptr<BeamParticleData>::transient_const_pointer _hadron_A,_hadron_B;
   
    double x(double xt, double v) const;
    double x_a(double x, double v) const;
    double x_b(double x, double v) const;
    double xbar(double v) const;
    double Ltilde_qg(double x, double v) const;
    double Ltilde_gq(double x, double v) const;
    double Ltilde_qq(double x, double v) const;
    double Vtilde_qq() const;
    double Ccalbar_qg(double x) const;
    double Fcal_qg(double x, double v) const;
    double Fcal_gq(double x, double v) const;
    double Fcal_qq(double x, double v) const;
    double Ftilde_qg(double xt, double v) const;
    double Ftilde_gq(double xt, double v) const;
    double Ftilde_qq(double xt, double v) const;
    double Ctilde_qg(double x, double v) const;
    double Ctilde_gq(double x, double v) const;
    double Ctilde_qq(double x, double v) const;

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

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEqq2gZ2ffNason> initMEqq2gZ2ffNason;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEqq2gZ2ffNason & operator=(const MEqq2gZ2ffNason &);

private:

  /**
   *  Pointer to the vertices for the helicity calculations
   */
  //@{
  /**
   *  Pointer to the Z vertex
   */
  AbstractFFVVertexPtr _theFFZVertex;

  /**
   *  Pointer to the photon vertex
   */
  AbstractFFVVertexPtr _theFFPVertex;
  //@}

  /**
   *  Pointers to the intermediate resonances
   */
  //@{
  /**
   *  Pointer to the Z ParticleData object
   */
  tcPDPtr _z0;

  /**
   *  Pointer to the photon ParticleData object
   */
  tcPDPtr _gamma;
  //@}

  /**
   *  Switches to control the particles in the hard process
   */
  //@{
  /**
   *  Allowed flavours for the incoming quarks
   */
  unsigned int _maxflavour;

  /**
   *  Whether to include both \f$Z^0\f$ and \f$\gamma\f$ or only one
   */
  unsigned int _gammaZ;
  
  /**
   *  Which processes to include
   */
  unsigned int _process;
  //@}

  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement _me;

  /**
   *  Parameters for the NLO weight
   */
  //@{
  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int _contrib;

  /**
   *  The magnitude of the correction term to reduce the negative contribution
   */
  double _a;

  /**
   *  The power of the correction term to reduce the negative contribution
   */
  double _p;

  /**
   * Histograms implemented as vectors of statistics to take into 
   * account the acdc sampling
   */
  mutable vector<Statistic> x_h, v_h, x_pos_h, v_pos_h, x_neg_h, v_neg_h;

  /**
   *  The cut-off on the pdfs
   */
  double _eps;

  //@}
  /**
   *  Radiation variables
   */
  //@{
  /**
   *   The \f$\tilde{x}\f$ variable
   */
  double _xt;

  /**
   *  The \f$v\f$ angular variable
   */
  double _v;

  //@}

  /**
   * Statistics & Histograms for testing.
   */
  mutable int no_wgts;
  mutable int no_negwgts;
  mutable double maxy, miny;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEqq2gZ2ffNason. */
template <>
struct BaseClassTrait<Herwig::MEqq2gZ2ffNason,1> {
  /** Typedef of the first base class of MEqq2gZ2ffNason. */
  typedef ME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEqq2gZ2ffNason class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEqq2gZ2ffNason>
  : public ClassTraitsBase<Herwig::MEqq2gZ2ffNason> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEqq2gZ2ffNason"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MEqq2gZ2ffNason class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwNasonME.so"; }
};

/** @endcond */

}

#include "MEqq2gZ2ffNason.icc"

#endif /* HERWIG_MEqq2gZ2ffNason_H */
