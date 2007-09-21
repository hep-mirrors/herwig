// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusteringSelector class.
//

#include "ClusteringSelector.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Utilities/Selector.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/UseRandom.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ClusteringSelector.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ClusteringSelector::~ClusteringSelector() {}

void ClusteringSelector::persistentOutput(PersistentOStream & ) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void ClusteringSelector::persistentInput(PersistentIStream & , int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<ClusteringSelector> ClusteringSelector::initClusteringSelector;
// Definition of the static class description member.


list<pair<ClusteringPtr,tClusteringGuidePtr> >
ClusteringSelector::sort (const list<pair<ClusteringPtr,tClusteringGuidePtr> >& clusterings) const {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "== ClusteringSelector::sort" << endl
		     << "got " << clusterings.size() << " clusterings" << endl;
#endif

  // first sort according to scales
  list<pair<ClusteringPtr,tClusteringGuidePtr> > sorted = clusterings;
  ClusteringScaleLess less;
  sorted.sort<ClusteringScaleLess>(less);
  // now look, if we got clusterings with equal scales

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "looking for degenerate clusterings" << endl;
#endif

  list<list<pair<ClusteringPtr,tClusteringGuidePtr> > > possiblyDegenerate;

  while (!sorted.empty()) {

    Energy2 currentScale = sorted.front().first->scale();
    list<pair<ClusteringPtr,tClusteringGuidePtr> > current;
    current.push_back(sorted.front());
    sorted.pop_front();

    while (sorted.front().first->scale() == currentScale && !sorted.empty()) {
      current.push_back(sorted.front());
      sorted.pop_front();
    }

    possiblyDegenerate.push_back(current);
    current.clear();

  }

  for (list<list<pair<ClusteringPtr,tClusteringGuidePtr> > >::iterator d= possiblyDegenerate.begin();
       d != possiblyDegenerate.end(); ++d) {
    if (d->size() == 1)
      sorted.push_back(d->front());
    else {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
      generator()->log() << "got degenerate set (clustering -> guide)" << endl;
      for(list<pair<ClusteringPtr,tClusteringGuidePtr> >::iterator dc = d->begin();
	  dc != d->end(); ++dc)
	generator()->log() << dc->first << " -> " << dc->second << "\n";
      generator()->log() << endl;
#endif
      pair<ClusteringPtr,tClusteringGuidePtr> select = selectDegenerate(*d);
      
#ifdef HERWIG_DEBUG_CKKW_EXTREME

      generator()->log() << "selected " << select.first << " -> "  << select.second << endl;

#endif

      sorted.push_back(select);
    }
  }
  return sorted;
}

pair<ClusteringPtr,tClusteringGuidePtr>
ClusteringSelector::selectDegenerate
(const list<pair<ClusteringPtr,tClusteringGuidePtr> >& deg) const {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "== ClusteringSelector::selectDegenerate" << endl;
#endif

  // normalize the weights

  double sum = 0.;

  for (list<pair<ClusteringPtr,tClusteringGuidePtr> >::const_iterator d= deg.begin(); d!=deg.end(); ++d) {
    if ((*d).first->weight() < 0.) 
      throw Exception () << "CKKW : ClusteringSelector::selectDegenerate : negative clustering weight."
			 << Exception::eventerror;
    sum += (*d).first->weight();
  }

  if (sum <= 0.) throw Exception () << "CKKW : ClusteringSelector::selectDegenerate : zero or negative sum"
				    << " of weights." << Exception::eventerror;

  for (list<pair<ClusteringPtr,tClusteringGuidePtr> >::const_iterator d= deg.begin(); d!=deg.end(); ++d)
    (*d).first->weight((*d).first->weight()/sum);

  // select a clustering

  Selector<pair<ClusteringPtr,tClusteringGuidePtr> > selector;
  for (list<pair<ClusteringPtr,tClusteringGuidePtr> >::const_iterator d= deg.begin(); d!=deg.end(); ++d)
    selector.insert((*d).first->weight(),*d);
  return selector.select(UseRandom::rnd());
}



void ClusteringSelector::Init() {

  static ClassDocumentation<ClusteringSelector> documentation
    ("ClusteringSelector is used to order clusterings occuring in a "
     "parton shower history reconstruction. The default version orders "
     "according to the scales associated with a clustering. In case "
     "there are different clusterings with the same scale, one is "
     "selected out of them accordin to their relative weights.");

}

