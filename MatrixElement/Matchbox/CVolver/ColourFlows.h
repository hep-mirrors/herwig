// -*- C++ -*-

//
// ColourFlowBasis.hpp is part of CVolver, (C) 2013 Simon Plätzer -- simon.plaetzer@desy.de
// CVolver is licenced under version 2 of the GPL, see COPYING for details.
//

#ifndef CVOLVER_ColourFlowBasis_hpp_included
#define CVOLVER_ColourFlowBasis_hpp_included

#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <algorithm>

namespace CVolver {

  /**
   * Generate the identical permutation.
   */
  template<class ForwardIterator, class T>
  void iota(ForwardIterator begin, ForwardIterator end, T value) {
    while ( begin != end ) {
      *begin = value; 
      value += 1;
      ++begin;
    }
  }

  /**
   * Return the identical permutation
   */
  inline std::vector<std::size_t> identicalPermutation(const std::size_t n) {
    std::vector<std::size_t> res(n);
    CVolver::iota(res.begin(),res.end(),0);
    return res;
  }

  /**
   * Return a random permutation
   */
  template<class Rnd>
  std::vector<std::size_t> randomPermutation(const std::size_t n, Rnd& rnd) {
    std::vector<std::size_t> pick = identicalPermutation(n);
    std::vector<std::size_t> res(n);
    std::size_t count = 0;
    while ( !pick.empty() ) {
      std::size_t i = rnd.template uniform_int_distribution<std::size_t>(0,pick.size()-1);
      res[count] = pick[i];
      pick.erase(pick.begin()+i);
      ++count;
    }
    return res;
  }

  /**
   * A colour flow (basis tensor).
   */
  class ColourFlow {

  public:

    /**
     * Default constructor
     */
    ColourFlow() {}

    /**
     * Construct a colour flow given a permutation
     */
    explicit ColourFlow(const std::vector<size_t>& perm)
      : thePermutation(perm) {}

    /**
     * Generate a random colour flow
     */
    template<class Rnd>
    static ColourFlow randomFlow(const std::size_t& n, Rnd& rnd) {
      return ColourFlow(randomPermutation(n,rnd));
    }

    /**
     * Generate all colour flows
     */
    static std::set<ColourFlow> allFlows(const std::size_t& n);

    /**
     * Compare for equality
     */
    bool operator==(const ColourFlow& other) const {
      return thePermutation == other.thePermutation;
    }

    /**
     * Compare for inequality
     */
    bool operator!=(const ColourFlow& other) const {
      return thePermutation != other.thePermutation;
    }

    /**
     * Compare for ordering
     */
    bool operator<(const ColourFlow& other) const {
      return thePermutation < other.thePermutation;
    }

    /**
     * Conjugate this basis tensor
     */
    ColourFlow& conjugate() {
      std::vector<std::size_t> tmp = thePermutation;
      for ( std::size_t k = 0; k < tmp.size(); ++k ) {
	thePermutation[tmp[k]] = k;
      }
      return *this;
    }

    /**
     * Return the scalar product with another basis tensor as the
     * resulting power of N
     */
    std::size_t scalarProduct(const ColourFlow& other) const;

    /**
     * Return the permutation
     */
    const std::vector<std::size_t>& permutation() const {
      return thePermutation;
    }

    /**
     * Return the anti-colour index connected to the given colour index
     */
    const std::size_t& antiColour(const std::size_t& i) const {
      assert(i < thePermutation.size());
      return thePermutation[i];
    }

    /**
     * Return the colour index connected to the given anti-colour index
     */
    std::size_t colour(const std::size_t& i) const {
      std::vector<std::size_t>::const_iterator k =
	std::find(thePermutation.begin(),thePermutation.end(),i);
      assert(k != thePermutation.end());
      return std::distance(thePermutation.begin(),k);
    }

    /**
     * Return the swapped indices, if this corresponds to a
     * transposition of the given one, or (0,0)
     */
    std::pair<std::size_t,std::size_t> getTranspositionOf(const ColourFlow& other) const;

    /**
     * Act a transposition on this index
     */
    ColourFlow& swap(const std::size_t& i, const std::size_t& j) {
      assert(i < thePermutation.size() && j < thePermutation.size());
      std::swap(thePermutation[i],thePermutation[j]);
      return *this;
    }

    /**
     * Add another flow to this basis tensor; this correspons to
     * emitting a `singlet' gluon.
     */
    ColourFlow& emitSinglet() {
      thePermutation.push_back(thePermutation.size());
      return *this;
    }

    /**
     * Emit a gluon on a given colour line
     */
    ColourFlow& emitFromColour(const std::size_t& i) {
      assert(i < thePermutation.size());
      thePermutation.push_back(thePermutation[i]);
      thePermutation[i] = thePermutation.size()-1;
      return *this;
    }

    /**
     * Emit a gluon on a given anti-colour line
     */
    ColourFlow& emitFromAntiColour(const std::size_t& i) {
      return emitFromColour(colour(i));
    }

    /**
     * Return true, if this colour flow is non-vanishing for the given
     * colours and anticolours
     */
    bool isNonZero(const std::vector<std::size_t>& colours,
		   const std::vector<std::size_t>& antiColours) const;

    /**
     * Return the number of coloured legs
     */
    size_t nLegs() const { return thePermutation.size(); }

  private:

    /**
     * The vector representing the permutation of anti-fundamental
     * w.r.t. fundamental indices.
     */
    std::vector<std::size_t> thePermutation;

  };

  /**
   * ParticleData traits
   */
  template<class ParticleData>
  struct ParticleDataTraits {

    /**
     * Return true, if singlet
     */
    static bool isSinglet(const ParticleData&) { return true; }

    /**
     * Return true, if anti-fundamental
     */
    static bool isAntiFundamental(const ParticleData&) { return false; }

    /**
     * Return true, if fundamental
     */
    static bool isFundamental(const ParticleData&) { return false; }

    /**
     * Return true, if adjoint
     */
    static bool isAdjoint(const ParticleData&) { return false; }

  };

  /**
   * The crossing of a physical process to a colour flow basis
   */
  class ColourFlowCrossing {

  private:

    /**
     * Add colour leg mapping
     */
    void addColourCrossing(const std::size_t& leg,
			   std::size_t& count,
			   const double& sign) {
      theColourMap[count] = leg;
      theColourCrossingSigns[count] = sign;
      theReverseColourMap[leg] = count;
      ++count;
    }

    /**
     * Add anti-colour leg mapping
     */
    void addAntiColourCrossing(const std::size_t& leg,
			       std::size_t& count,
			       const double& sign) {
      theAntiColourMap[count] = leg;
      theAntiColourCrossingSigns[count] = sign;
      theReverseAntiColourMap[leg] = count;
      ++count;
    }

  public:

    /**
     * Default constructor
     */
    ColourFlowCrossing()
      : theNFlows(0) {}

    /**
     * Construct for the given process
     */
    template<class ParticleData>
    explicit ColourFlowCrossing(const std::vector<ParticleData>& proc,
				bool signs = true) 
      : theNFlows(0) {
      typedef ParticleDataTraits<ParticleData> Traits;
      std::size_t colourCounter = 0;
      std::size_t antiColourCounter = 0;
      for ( std::size_t k = 0; k < proc.size(); ++k ) {
	if ( Traits::isSinglet(proc[k]) )
	  continue;
	double sign = k > 1 ? 1. : -1.;
	if ( Traits::isAntiFundamental(proc[k]) ) {
	  if ( k > 1 )
	    addAntiColourCrossing(k,antiColourCounter,signs ? sign : 1.0);
	  else
	    addColourCrossing(k,colourCounter,signs ? sign : 1.0);
	}
	if ( Traits::isFundamental(proc[k]) ) {
	  if ( k > 1 )
	    addColourCrossing(k,colourCounter,signs ? sign : 1.0);
	  else
	    addAntiColourCrossing(k,antiColourCounter,signs ? sign : 1.0);
	}
	if ( Traits::isAdjoint(proc[k]) ) {
	  addColourCrossing(k,colourCounter,1.0);
	  addAntiColourCrossing(k,antiColourCounter,1.0);
	}
      }
      theNFlows = colourCounter;
    }

    /**
     * Return the number of colour flows
     */
    const std::size_t& nFlows() const { return theNFlows; }

    /**
     * Return the external leg for the given colour
     */
    std::size_t colourLeg(const std::size_t& i) const {
      std::map<std::size_t,std::size_t>::const_iterator l =
	theColourMap.find(i);
      assert(l != theColourMap.end());
      return l->second;
    }

    /**
     * Return the external leg for the given anti-colour
     */
    std::size_t antiColourLeg(const std::size_t& i) const {
      std::map<std::size_t,std::size_t>::const_iterator l =
	theAntiColourMap.find(i);
      assert(l != theAntiColourMap.end());
      return l->second;
    }

    /**
     * Return the crossing sign for the given colour
     */
    std::size_t colourCrossingSign(const std::size_t& i) const {
      std::map<std::size_t,double>::const_iterator l =
	theColourCrossingSigns.find(i);
      assert(l != theColourCrossingSigns.end());
      return l->second;
    }

    /**
     * Return the crossing sign for the given anti-colour
     */
    double antiColourCrossingSign(const std::size_t& i) const {
      std::map<std::size_t,double>::const_iterator l =
	theAntiColourCrossingSigns.find(i);
      assert(l != theAntiColourCrossingSigns.end());
      return l->second;
    }

    /**
     * Return true, if the external line carries colour
     */
    bool coloured(const std::size_t& i) const {
      return theReverseColourMap.find(i) != theReverseColourMap.end();
    }

    /**
     * Return the colour line for the given external leg
     */
    std::size_t colourLine(const std::size_t& i) const {
      std::map<std::size_t,std::size_t>::const_iterator l =
	theReverseColourMap.find(i);
      assert(l != theReverseColourMap.end());
      return l->second;
    }

    /**
     * Return true, if the external line carries anti-colour
     */
    bool antiColoured(const std::size_t& i) const {
      return theReverseAntiColourMap.find(i) != theReverseAntiColourMap.end();
    }

    /**
     * Return the anti-colour line for the given external leg
     */
    std::size_t antiColourLine(const std::size_t& i) const {
      std::map<std::size_t,std::size_t>::const_iterator l =
	theReverseAntiColourMap.find(i);
      assert(l != theReverseAntiColourMap.end());
      return l->second;
    }

  private:

    /**
     * The number of colour flows
     */
    std::size_t theNFlows;

    /**
     * Map colour legs to external legs
     */
    std::map<std::size_t,std::size_t> theColourMap;

    /**
     * Map anti-colour legs to external legs
     */
    std::map<std::size_t,std::size_t> theAntiColourMap;

    /**
     * Map external legs to colour legs
     */
    std::map<std::size_t,std::size_t> theReverseColourMap;

    /**
     * Map external legs to anti-colour legs
     */
    std::map<std::size_t,std::size_t> theReverseAntiColourMap;

    /**
     * Map colour legs to crossing signs
     */
    std::map<std::size_t,double> theColourCrossingSigns;

    /**
     * Map anti-colour legs to crossing signs
     */
    std::map<std::size_t,double> theAntiColourCrossingSigns;

  };

}

#endif // CVOLVER_ColourFlowBasis_hpp_included
