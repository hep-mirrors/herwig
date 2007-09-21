// -*- C++ -*-

#ifndef HERWIG_Partitioner_H
#define HERWIG_Partitioner_H

#include <vector>
#include <map>
#include <algorithm>

namespace Herwig {

  using namespace std;

  template<class T, class S>
  void appendvector (vector<T>& to, const vector<S>& from) {
    for (unsigned int f = 0; f < from.size(); ++f)
      to.push_back(from[f]);
  }

  /**\ingroup CKKW
   *
   * Partitioner returns all possible ordered
   * partitions from a given vector. It stores information
   * which can be reused in building partitions.
   *
   *@author Simon Plaetzer
   *
   */
  class Partitioner {

  public:

    /**
     * Return all partitions of size degree out of the
     * input container in
     */
    template<class T>
    vector<vector<T> > partitions (const vector<T>& in, unsigned int degree) {
      unsigned int size = in.size();
      map<pair<unsigned int, unsigned int>, vector<vector<unsigned int> > >::iterator
	iset = _knownIndexSets.find(make_pair(size,degree));
      if (iset == _knownIndexSets.end()) {
	vector<vector<unsigned int> > newSet;
	index_partitions(0,size,degree,newSet);
	iset = _knownIndexSets.insert(make_pair(make_pair(size,degree),newSet)).first;
      }
      vector<vector<T> > result;
      for (vector<vector<unsigned int> >::iterator is = iset->second.begin();
	   is != iset->second.end(); ++is)
	result.push_back(project(in,*is));
      return result;
    }

    /**
     * For an index range [begin,end), get all ordered
     * partitions of size degree, giving the possible
     * index sets for partitions to be done.
     */
    void index_partitions (const unsigned int& begin, const unsigned int& end,
			   const unsigned int& degree,
			   vector<vector<unsigned int> >& partitions,
			   const unsigned int& current = 1,
			   const vector<unsigned int>& indices
			   = vector<unsigned int>()) {
      unsigned int dummy = 0;
      vector<unsigned int> temp = indices;
      temp.push_back(dummy);
      unsigned int next_current = current; next_current++;

      if (current < degree) {
	for (unsigned int i=begin; i<end; i++) {
	  temp.back() = i;
	  if (find(temp.begin(),--temp.end(),i) == --temp.end()) {
	    index_partitions (i, end, degree, partitions, next_current, temp);
	  }
	}
      }
      else {
	for(unsigned int i=begin; i<end; i++) {
	  temp.back() = i;
	  if (find(temp.begin(),--temp.end(),i)==--temp.end()) {
	    partitions.push_back(temp);
	  }
	}
      }
    }

    /**
     * Return all members of given random access container
     * at positions listed in indexSet
     */
    template<class T>
    vector<T> project (const vector<T>& in, const vector<unsigned int>& indexSet) const {
      vector<T> temp;
      for (unsigned int i = 0; i< in.size(); ++i)
	if (find(indexSet.begin(),indexSet.end(),i)!=indexSet.end())
	  temp.push_back(in[i]);
      return temp;
    }

    /**
     * Return all members of given random access container
     * at positions not listed in indexSet
     */
    template<class T>
    vector<T> projectComplement (const vector<T>& in,
				 const vector<unsigned int>& indexSet) const {
      vector<T> temp;
      for (unsigned int i = 0; i< in.size(); ++i)
	if (find(indexSet.begin(),indexSet.end(),i)==indexSet.end())
	  temp.push_back(in[i]);
      return temp;  
    }

  private:

    /**
     * Cache already computed index sets by size of input
     * and size of partitions.
     */
    map<pair<unsigned int, unsigned int>, vector<vector<unsigned int> > >
    _knownIndexSets;


  };

}

#endif // HERWIG_Partitioner_H
