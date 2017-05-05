// -*- C++ -*-
//
// ColourBasis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ColourBasis_H
#define HERWIG_ColourBasis_H
//
// This is the declaration of the ColourBasis class.
//

#include "ThePEG/Handlers/HandlerBase.h"

#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/MatrixElement/MEBase.h"

#include "Herwig/MatrixElement/Matchbox/Utility/MatchboxXComb.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.fh"

#include <iterator>

namespace Herwig {

using std::iterator_traits;
using std::distance;

using namespace ThePEG;

using boost::numeric::ublas::matrix;
using boost::numeric::ublas::symmetric_matrix;
using boost::numeric::ublas::compressed_matrix;
using boost::numeric::ublas::upper;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief ColourBasis is an interface to a colour basis
 * implementation.
 *
 */
class ColourBasis: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ColourBasis();

  /**
   * The destructor.
   */
  virtual ~ColourBasis();
  //@}

public:

  /**
   * Return the factory which produced this matrix element
   */
  Ptr<MatchboxFactory>::tptr factory() const;

  /**
   * Set the factory which produced this matrix element
   */
  void factory(Ptr<MatchboxFactory>::tptr f);

  /**
   * Clone this colour basis.
   */
  Ptr<ColourBasis>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<ColourBasis>::ptr>(clone());
  }

  /**
   * Clear this colour basis
   */
  virtual void clear();

  /**
   * Prepare for the given sub process and return the basis
   * dimensionality.
   */
  size_t prepare(const cPDVector&, bool);

  /**
   * Prepare for the given diagrams.
   */
  size_t prepare(const MEBase::DiagramVector&, bool);

  /**
   * Return the index map.
   */
  const map<cPDVector,map<size_t,size_t> >& indexMap() const { return theIndexMap; }

  /**
   * Return a map of basis tensor indices to vectors identifying a
   * certain ordering corresponding to the given colour structure. May
   * not be supported by all colour basis implementations.
   */
  virtual map<size_t,vector<vector<size_t> > > basisList(const vector<PDT::Colour>&) const {
    return map<size_t,vector<vector<size_t> > >();
  }

  /**
   * Given a physical subprocess, a colour to amplitude label map and
   * a basis tensor index, return an identifier of the ordering
   * coresponding to the given colour structure. This will only return
   * sensible results for colour bases which implement the basisList
   * query.
   */
  const string& orderingString(const cPDVector& sub, 
			       const map<size_t,size_t>& colourToAmplitude,
			       size_t tensorId);

  /**
   * Given a physical subprocess, a colour to amplitude label map and
   * a basis tensor index, return an identifier of the ordering
   * coresponding to the given colour structure. This will only return
   * sensible results for colour bases which implement the basisList
   * query.
   */
  const set<vector<size_t> >& ordering(const cPDVector& sub, 
				       const map<size_t,size_t>& colourToAmplitude,
				       size_t tensorId, size_t shift = 0);

  /**
   * For the given subprocess and amplitude vectors
   * calculate the amplitude squared.
   */
  double me2(const cPDVector&, const map<vector<int>,CVector>&) const;
 
  /**
   * For the given subprocess and amplitude vectors
   * calculate the interference.
   */
  double interference(const cPDVector&, 
		      const map<vector<int>,CVector>&,
		      const map<vector<int>,CVector>&) const;

  /**
   * For the given subprocess and amplitude vector
   * calculate the colour correlated amplitude.
   */
  double colourCorrelatedME2(const pair<size_t,size_t>&,
			     const cPDVector&, 
			     const map<vector<int>,CVector>&) const;

  /**
   * For the given subprocess and amplitude vector
   * calculate the amplitude squared.
   */
  Complex interference(const cPDVector&, 
		       const CVector&, const CVector&) const;

  /**
   * For the given subprocess and amplitude vector
   * calculate the colour correlated amplitude.
   */
  Complex colourCorrelatedInterference(const pair<size_t,size_t>&,
				       const cPDVector&, 
				       const CVector&, const CVector&) const;

  /**
   * For the given subprocess and amplitude given as amp amp^\dagger
   * calculate the amplitude squared.
   */
  double me2(const cPDVector&, const matrix<Complex>&) const;

  /**
   * For the given subprocess and amplitude given as amp amp^\dagger
   * calculate the colour correlated amplitude.
   */
  double colourCorrelatedME2(const pair<size_t,size_t>&,
			     const cPDVector&, 
			     const matrix<Complex>&) const;

  /**
   * Return the scalar product matrix for the given process.
   */
  const symmetric_matrix<double,upper>& scalarProducts(const cPDVector&) const;

  /**
   * Return the matrix representation of a colour charge.
   */
  const compressed_matrix<double>& charge(const cPDVector&, size_t) const;

  /**
   * Return the non-vanishing elements of a colour charge.
   */
  const vector<pair<size_t,size_t> >& chargeNonZero(const cPDVector&, size_t) const;

  /**
   * Return the correlator matrix for the given process.
   */
  const symmetric_matrix<double,upper>& correlator(const cPDVector&,
						   const pair<size_t,size_t>&) const;

  /**
   * Return true, if the colour basis is capable of assigning colour
   * flows.
   */
  virtual bool haveColourFlows() const { return false; }

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  Selector<const ColourLines *> colourGeometries(tcDiagPtr diag,
						 const map<vector<int>,CVector>& amps);

  /**
   * Return the colour tensor used for the selected colour flow
   */
  size_t tensorIdFromFlow(tcDiagPtr diag, const ColourLines * cl);

  /**
   * Match colour representation.
   */
  struct matchRep {
    PDT::Colour m;
    matchRep(PDT::Colour n)
      : m(n) {}
    bool operator()(PDT::Colour c) const {
      return c == m;
    }
  };

  /**
   * Return true, if this basis is running in large-N mode
   */
  virtual bool largeN() const { return theLargeN; }

  /**
   * Switch to large n
   */
  void doLargeN(bool yes = true) { theLargeN = yes; }

  /**
   * Convert particle data to colour information
   */
  vector<PDT::Colour> projectColour(const cPDVector&) const;

  /**
   * Perform a normal ordering of the external legs. This default
   * implementation assumes normal ordered legs as 3 3bar ... 3 3bar 8 ... 8
   * while removing all non-coloured particles.
   */
  virtual vector<PDT::Colour> normalOrder(const vector<PDT::Colour>&) const;

  /**
   * Determine the mapping of process to colour indices and return the
   * normal ordered vector of colour indices
   */
  vector<PDT::Colour> normalOrderMap(const cPDVector& sub);

  /**
   * Get the normal ordered legs
   */
  const vector<PDT::Colour>& normalOrderedLegs(const cPDVector& sub) const;

  /**
   * Convert the legs to a string.
   */
  string file(const vector<PDT::Colour>&) const;

  /**
   * Calculate T_i^\dagger X T_j
   */
  void chargeProduct(const compressed_matrix<double>& ti,
		     const vector<pair<size_t,size_t> >& tiNonZero,
		     const symmetric_matrix<double,upper>& X,
		     const compressed_matrix<double>& tj,
		     const vector<pair<size_t,size_t> >& tjNonZero,
		     symmetric_matrix<double,upper>& result) const;

  /**
   * Calculate T_i X T_j^\dagger
   */
  void chargeProductAdd(const compressed_matrix<double>& ti,
			const vector<pair<size_t,size_t> >& tiNonZero,
			const matrix<Complex>& X,
			const compressed_matrix<double>& tj,
			const vector<pair<size_t,size_t> >& tjNonZero,
			matrix<Complex>& result,
			double factor = 1.) const;

public:

  /**
   * Find a coloured path from a to b within the given diagram.
   */
  static list<pair<int,bool> > colouredPath(pair<int,bool> a, pair<int,bool> b,
					    Ptr<Tree2toNDiagram>::tcptr);

  /**
   * Get all colour flows for the given diagram.
   */
  static list<list<list<pair<int,bool> > > > colourFlows(Ptr<Tree2toNDiagram>::tcptr);

  /**
   * Convert a flow to a string representation appropriate for
   * ColourLines
   */
  static string cfstring(const list<list<pair<int,bool> > >&);

protected:

  /**
   * Prepare the basis for the normal ordered legs and return the
   * dimensionality of the basis.
   */
  virtual size_t prepareBasis(const vector<PDT::Colour>&) = 0;

  /**
   * Return the scalar product of basis tensors labelled a and b in
   * the basis used for the given normal ordered legs.
   */
  virtual double scalarProduct(size_t a, size_t b,
			       const vector<PDT::Colour>& abBasis) const = 0;

  /**
   * Return the matrix element of a colour charge
   * <c_{n+1,a}|T_i|c_{n,b}> between basis tensors a and b, with
   * respect to aBasis and bBasis
   */
  virtual double tMatrixElement(size_t i, size_t a, size_t b,
				const vector<PDT::Colour>& aBasis,
				const vector<PDT::Colour>& bBasis) const = 0;

  /**
   * Return true, if a large-N colour connection exists for the
   * given external legs and basis tensor.
   */
  virtual bool colourConnected(const cPDVector&,
			       const vector<PDT::Colour>&,
			       const pair<int,bool>&, 
			       const pair<int,bool>&, 
			       size_t) const;

  /**
   * Return true, if a large-N colour connection exists for the
   * given external legs and basis tensor.
   */
  virtual bool colourConnected(const vector<PDT::Colour>&,
			       int, int, size_t) const {
    return false;
  }

  /**
   * Match up colour flows for given diagram to basis tensors.
   */
  vector<string> makeFlows(Ptr<Tree2toNDiagram>::tcptr, size_t) const;

  /**
   * Return the colour line map.
   */
  map<Ptr<Tree2toNDiagram>::tcptr,vector<ColourLines*> >&
  colourLineMap();

  /**
   * Update the colour line map for a given diagram.
   */
  void updateColourLines(Ptr<Tree2toNDiagram>::tcptr);

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The factory which produced this matrix element
   */
  Ptr<MatchboxFactory>::tptr theFactory;

  typedef map<vector<PDT::Colour>,symmetric_matrix<double,upper> >
  ScalarProductMap;

  typedef map<vector<PDT::Colour>,map<size_t,compressed_matrix<double> > > ChargeMap;
  typedef map<vector<PDT::Colour>,map<size_t,vector<pair<size_t,size_t > > > > ChargeNonZeroMap;

  typedef map<vector<PDT::Colour>,map<pair<size_t,size_t>,symmetric_matrix<double,upper> > > CorrelatorMap;

  /**
   * True, if this basis is running in large-N mode
   */
  bool theLargeN;

  /**
   * Map external legs to normal ordered versions
   */
  map<cPDVector,vector<PDT::Colour> > theNormalOrderedLegs;

  /**
   * Index mappings to normal order from given leg assignments,
   * indexed by the original leg assignment.
   */
  map<cPDVector,map<size_t,size_t> > theIndexMap;

  /**
   * The scalar product matrix S_n = <c_{n,a}|c_{n,b}> , indexed
   * by normal ordered leg assignments.
   */
  ScalarProductMap theScalarProducts;

  /**
   * The colour charge matrices <c_{n+1,a}|T_i|c_{n,b}> indexed by
   * the `n' normal ordered legs and the index i.
   */
  ChargeMap theCharges;

  /**
   * The nonzero elements of the charge matrices.
   */
  ChargeNonZeroMap theChargeNonZeros;

  /**
   * The correlator matrices T_i\cdot T_j -> T_i^\dagger S_{n+1} T_j
   * with T_i = <c_{n+1,a}|T_i|c_{n,b}> indexed by the `n' basis
   * normal ordered legs and indices i,j
   */
  CorrelatorMap theCorrelators;

  /**
   * Map diagrams to colour flows indexed by basis tensor.
   */
  map<Ptr<Tree2toNDiagram>::tcptr,vector<string> > theFlowMap;

  /**
   * Map diagrams to colour line objects.
   */
  map<Ptr<Tree2toNDiagram>::tcptr,vector<ColourLines*> > theColourLineMap;

  /**
   * Store ordering identifiers
   */
  map<cPDVector,map<size_t,string> > theOrderingStringIdentifiers;

  /**
   * Store ordering identifiers
   */
  map<cPDVector,map<size_t,set<vector<size_t> > > > theOrderingIdentifiers;

  /**
   * Write out yet unknown basis computations.
   */
  void writeBasis(const string& prefix = "") const;

  /**
   * Read in the basis computation which are supposed to be known.
   */
  void readBasis();

  /**
   * Read in the basis computation which are supposed to be known.
   */
  bool readBasis(const vector<PDT::Colour>&);

  /**
   * Gather any implementation dependend details when reading a basis
   */
  virtual void readBasisDetails(const vector<PDT::Colour>&) {}

  /**
   * Write out symmetric matrices.
   */
  void write(const symmetric_matrix<double,upper>&, ostream&) const;

  /**
   * Read in symmetric matrices.
   */
  void read(symmetric_matrix<double,upper>&, istream&);

  /**
   * Write out compressed matrices.
   */
  void write(const compressed_matrix<double>&, ostream&,
	     const vector<pair<size_t,size_t> >&) const;

  /**
   * Read in compressed matrices.
   */
  void read(compressed_matrix<double>&, istream&,
	    vector<pair<size_t,size_t> >&);

  /**
   * True, if an attempt to read in basis information has been
   * completed.
   */
  bool didRead;

  /**
   * True, if an attempt to write out basis information has been
   * completed.
   */
  mutable bool didWrite;

  /**
   * Temporary storage.
   */
  matrix<double> tmp;

  /**
   * The search path
   */
  string theSearchPath;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ColourBasis & operator=(const ColourBasis &);

};

}

#endif /* HERWIG_ColourBasis_H */
