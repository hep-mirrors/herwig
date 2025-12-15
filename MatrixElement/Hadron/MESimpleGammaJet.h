// -*- C++ -*-
#ifndef Herwig_MESimpleGammaJet_H
#define Herwig_MESimpleGammaJet_H
//
// This is the declaration of the MESimpleGammaJet class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MESimpleGammaJet class.
 *
 * @see \ref MESimpleGammaJetInterfaces "The interfaces"
 * defined for MESimpleGammaJet.
 */
class MESimpleGammaJet: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MESimpleGammaJet();

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  virtual unsigned int orderInAlphaS() const;
  virtual unsigned int orderInAlphaEW() const;
  virtual double me2() const;
  virtual Energy2 scale() const;
  virtual void getDiagrams() const;
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
  virtual void constructVertex(tSubProPtr);
  //@}

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  void persistentOutput(PersistentOStream & os) const;
  void persistentInput (PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /** @name Standard Interfaced functions. */
  //@{
  virtual void doinit();
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  virtual IBPtr clone() const;
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Matrix elements for the different subprocesses */
  //@{
  double qqbarME(std::vector<SpinorWaveFunction> & fin,
                 std::vector<SpinorBarWaveFunction> & ain,
                 std::vector<VectorWaveFunction> & gout,
                 std::vector<VectorWaveFunction> & Aout,
                 bool me=false) const;

  double qgME(std::vector<SpinorWaveFunction> & fin,
              std::vector<VectorWaveFunction> & gin,
              std::vector<SpinorBarWaveFunction> & fout,
              std::vector<VectorWaveFunction> & Aout,
              bool me=false) const;

  double qbargME(std::vector<SpinorBarWaveFunction> & fin,
                 std::vector<VectorWaveFunction> & gin,
                 std::vector<SpinorWaveFunction> & fout,
                 std::vector<VectorWaveFunction> & Aout,
                 bool me=false) const;
  //@}

private:

  MESimpleGammaJet & operator=(const MESimpleGammaJet &);

private:

  /** @name Vertices/pointers used in the helicity ME */
  //@{
  AbstractFFVVertexPtr _theFFGammaVertex;
  AbstractFFVVertexPtr _theQQGVertex;
  tcPDPtr _gamma;
  //@}

  /** @name Switches/parameters controlling allowed subprocesses and flavours */
  //@{
  unsigned int _process;
  unsigned int _maxflavour;
  //@}

  /** Matrix element for spin correlations */
  mutable ProductionMatrixElement _me;

};

}

#endif /* Herwig_MESimpleGammaJet_H */
