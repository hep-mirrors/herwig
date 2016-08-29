// -*- C++ -*-
//
// Tree2toNGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_Tree2toNGenerator_H
#define Herwig_Tree2toNGenerator_H
//
// This is the declaration of the Tree2toNGenerator class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Generate Tree2toNDiagrams for a given process.
 *
 * @see \ref Tree2toNGeneratorInterfaces "The interfaces"
 * defined for Tree2toNGenerator.
 */
class Tree2toNGenerator: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  Tree2toNGenerator();

  /**
   * The destructor.
   */
  virtual ~Tree2toNGenerator();
  //@}

public:

  /**
   * Generate all diagrams for the given process.
   */
  vector<Ptr<Tree2toNDiagram>::ptr> generate(const PDVector&,
					     unsigned int orderInGs,
					     unsigned int orderInGem);

  typedef vector<Ptr<Helicity::VertexBase>::ptr> VertexVector;

  /**
   * Access the vertices
   */
  VertexVector& vertices() { return theVertices; }

  /**
   * Return the vertices
   */
  const VertexVector& vertices() const { return theVertices; }

  /**
   * Access the particles to be excluded from internal lines
   */
  PDVector& excludeInternal() { return theExcludeInternal; }

  /**
   * Return the particles to be excluded from internal lines
   */
  const PDVector& excludeInternal() const { return theExcludeInternal; }

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
   * A node in internally used trees.
   */
  struct Vertex {

    /**
     * The outgoing particles. If this is a spacelike node, the first
     * child is considered the next spacelike (or second incoming)
     * line. If children are empty, this is an external line.
     */
    vector<Vertex> children;

    /**
     * The incoming line at this node.
     */
    PDPtr parent;

    /**
     * True, if this is spacelike node.
     */
    bool spacelike;

    /**
     * The external leg id.
     */
    int externalId;

    /**
     * The parent diagram id.
     */
    int parentId;

    /**
     * The default constructor.
     */
    Vertex()
      : spacelike(false), externalId(-1), parentId(-1) {}

    /**
     * Debug printout.
     */
    void print(ostream& os, const string& prefix = "") const {
      os << prefix << parent->PDGName()
	 << "[" << (spacelike ? "s" : "t") << "] (";
      if ( externalId < 0 )
	os << "x)\n";
      else
	os << externalId << ")\n";
      if ( !children.empty() ) {
	os << prefix << "|__\n";
	children[0].print(os,prefix + "|  ");
	os << prefix << "|__\n";
	children[1].print(os,prefix + "|  ");
      }
    }

    /**
     * Count the number of spacelike lines
     */
    int nspace() const {
      if ( children.empty() )
	return 1;
      int ret = 1;
      ret += children[0].nspace();
      return ret;
    }

    /**
     * Update diagram returning a map of external ids to diagram id
     * parents.
     */
    void update(Tree2toNDiagram& diag,
		map<int,pair<int,PDPtr> >& outgoing,
		int& lastUsed) {
      if ( externalId == 0 ) {
	assert(lastUsed==0);
	++lastUsed;
	diag.operator,(parent);
	children[0].parentId = lastUsed;
	children[1].parentId = lastUsed;
	children[0].update(diag,outgoing,lastUsed);
	children[1].update(diag,outgoing,lastUsed);
	for ( map<int,pair<int,PDPtr> >::iterator out =
		outgoing.begin(); out != outgoing.end(); ++out ) {
	  diag.operator,(out->second.first);
	  diag.operator,(out->second.second);
	}
	return;
      }
      if ( spacelike ) {
	++lastUsed;
	diag.operator,(parent);
	if ( externalId == 1 )
	  return;
	children[0].parentId = lastUsed;
	children[1].parentId = lastUsed;
	children[0].update(diag,outgoing,lastUsed);
	children[1].update(diag,outgoing,lastUsed);
	return;
      }
      if ( children.empty() ) {
	outgoing[externalId] =
	  make_pair(parentId,parent);
	return;
      }
      diag.operator,(parentId);
      diag.operator,(parent);
      ++lastUsed;
      children[0].parentId = lastUsed;
      children[1].parentId = lastUsed;
      children[0].update(diag,outgoing,lastUsed);
      children[1].update(diag,outgoing,lastUsed);
    }

    /**
     * Generate a diagram of given id.
     */
    Tree2toNDiagram generate(int id) {
      int nsp = nspace();
      Tree2toNDiagram res(nsp);
      int diagid = 0;
      map<int,pair<int,PDPtr> > out;
      update(res,out,diagid);
      res.operator,(-id);
      return res;
    }

  };

  /**
   * For the given set of trees determine all allowed clusterings.
   */
  list<vector<Vertex> > cluster(const vector<Vertex>& children,
				unsigned int orderInGs,
				unsigned int orderInGem) const;
  
  /**
   * For the given set of outgoing lines cluster recursively.
   */
  list<vector<Vertex> > clusterAll(const list<vector<Vertex> >& current,
				   unsigned int orderInGs,
				   unsigned int orderInGem) const;

  /**
   * For the given set of outgoing lines cluster recursively.
   */
  list<vector<Vertex> > clusterAll(const PDVector& external,
				   unsigned int orderInGs,
				   unsigned int orderInGem);

  /**
   * Helper for topology restrictions
   */
  struct LineMatcher {

    /**
     * The group of lines to be considered
     */
    set<tcPDPtr> particles;

    /**
     * The range allowed
     */
    pair<int,int> range;

    /**
     * The current count
     */
    int count;

    /**
     * Default constructor
     */
    LineMatcher()
      : range(0,0), count(0) {}

    /**
     * Construct given particles and a range
     */
    LineMatcher(const PDVector& p,
		const pair<int,int>& r)
      : range(r), count(0) {
      copy(p.begin(),p.end(),inserter(particles,particles.begin()));
    }

    /**
     * Rebind the particle data pointers
     */
    void rebind(Tree2toNGenerator* g) {
      set<tcPDPtr> oldp = particles;
      particles.clear();
      for ( set<tcPDPtr>::const_iterator p = oldp.begin();
	    p != oldp.end(); ++p )
	particles.insert(g->getParticleData((**p).id()));
    }

    /**
     * Reset this matcher
     */
    void reset() {
      count = 0;
    }

    /**
     * Count the given multiplicity
     */
    void add(tcPDPtr p, int n) {
      if ( particles.find(p) == particles.end() )
	return;
      count += n;
    }

    /**
     * Ceck if restrictions are met
     */
    bool check() const {
      return
	count >= range.first && count <= range.second;
    }

  };
  
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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The vertices to be used.
   */
  VertexVector theVertices;

  /**
   * The particles to be excluded from internal lines
   */
  PDVector theExcludeInternal;

  /**
   * Maximum order in gs to consider.
   */
  unsigned int maxOrderGs;

  /**
   * Maximum order in gem to consider.
   */
  unsigned int maxOrderGem;

  /**
   * Wether or not the generator has been prepared
   */
  bool prepared;

  /**
   * The vertices to be excluded.
   */
  VertexVector theExcludeVertices;

  /**
   * Minimal and maximal occurences of spacelike internal lines
   */
  vector<LineMatcher> spaceLikeAllowed;

  /**
   * Minimal and maximal occurences of timelike internal lines
   */
  vector<LineMatcher> timeLikeAllowed;

  /**
   * The next particle for which internal lines need to be restricted
   */
  PDVector theRestrictLines;

  /**
   * Command to set an allowed range of spacelike internal lines
   */
  string doSpaceLikeRange(string);

  /**
   * Command to set an allowed range of timelike internal lines
   */
  string doTimeLikeRange(string);

  /**
   * Command to clear the restrict lines container
   */
  string doClearRestrictLines(string);

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Tree2toNGenerator & operator=(const Tree2toNGenerator &);

};

inline PersistentOStream& operator<<(PersistentOStream& os, const Tree2toNGenerator::LineMatcher& m) {
  os << m.particles << m.range << m.count;
  return os;
}

inline PersistentIStream& operator>>(PersistentIStream& is, Tree2toNGenerator::LineMatcher& m) {
  is >> m.particles >> m.range >> m.count;
  return is;
}

}

#endif /* Herwig_Tree2toNGenerator_H */
