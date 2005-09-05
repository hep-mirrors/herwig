#ifndef HERWIG_MRST_H
#define HERWIG_MRST_H

#include <ThePEG/PDF/PDFBase.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "MRSTData.h"

namespace Herwig {

using namespace ThePEG;

static const int np=8;
static const int nx=49;
static const int nq=37;
static const int nqc0=2; // Charm not introduced before 2nd bin in q^2
static const int nqb0=11; // Bottom not introduced before 11th bin in q^2
static const double xmin=1E-5;
static const double xmax=1.0;
static const double qsqmin=1.25;
static const double qsqmax=1E7;
static const double mc2=2.045;
static const double mb2=18.5;

typedef Ptr<MRSTData>::pointer MRSTDatPtr;

/** \ingroup PDF
 *
 *  Some comment should be provided!
 */
class MRST : public PDFBase {
 
public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  MRST();

  /**
   * Copy-constructor.
   */
  MRST(const MRST &);

  /**
   * Destructor.
   */
  virtual ~MRST();
  //@}


  /**
   * Return true if this PDF can handle the extraction of parton from the
   * given particle ie. if the particle is a proton or neutron.
   */
  virtual bool canHandleParticle(tcPDPtr particle) const;

  /**
   * Return the parton types which are described by these parton
   * densities.
   */
  virtual cPDVector partons(tcPDPtr p) const;

  /**
   * Return the true pdf for the given parameters, with the momentum
   * fraction given as l=log(1/x) and simply x respectively.
   */
  virtual double xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                     double l, Energy2 particleScale = 0.0*GeV2) const;

  /**
   * 
   */
  virtual double xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                     double x, double eps = 0.0,
                     Energy2 particleScale = 0.0*GeV2) const;

  virtual double xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                     double l, Energy2 particleScale = 0.0*GeV2) const;
  virtual double xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                      double x, double eps = 0.0,
                      Energy2 particleScale = 0.0*GeV2) const;

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();
  
  virtual IBPtr clone() const;
  virtual IBPtr fullclone() const;

 protected:

  /**
   * Standard Interfaced virtual functions.
   */
  virtual void doupdate() throw(UpdateException);
  virtual void doinit() throw(InitException);
  //virtual void doinitrun();
  virtual void dofinish();

  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  virtual void rebind(const TranslationMap & trans) throw(RebindException);

  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  virtual IVector getReferences();

private:

  string _file;
  //MRSTDatPtr dataPtr;

  void initialize(bool reread = true);
  virtual void readSetup(istream &) throw(SetupException);

  static ClassDescription<MRST> initMRST;

  /**
   * Private and non-existent assignment operator.
   */
  MRST & operator=(const MRST &);

  double data[np+1][nx+1][nq+1];

 private:
  static double xx[nx+1];
  static double qq[nq+1];
  double c[np+1][nx][nq][5][5]; //coefficients used for interpolation
  double table[np+1]; // The values for each parton at a given scale

  enum { upValence = 1, dnValence, glu, upSea, chm, str, bot, dnSea };
	 
  int locate(double xx[],int n,double x);
  double polderivative(double x1, double x2, double x3,
		       double y1, double y2, double y3);

  /**
   * This function calculates the values for the given x and q
   */
  void update(double x, double q);
};

}

namespace ThePEG {

template <>
struct BaseClassTrait<Herwig::MRST,1> {
  typedef PDFBase NthBase;
};

template <>
struct ClassTraits<Herwig::MRST>: public ClassTraitsBase<Herwig::MRST> {
  static string className() { return "/Herwig++/PDF/MRST"; }
  static string library() { return "MRST.so"; }
};

}

#endif
