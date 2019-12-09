// -*- C++ -*-
//
// Remapper.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <cstdlib>

#include "Remapper.h"

using namespace Herwig;
using namespace std;

Remapper::Remapper()
  : minSelection(0.0), smooth(false) {}

Remapper::Remapper(unsigned int nBins,
		   double nMinSelection,
		   bool nSmooth) 
  : minSelection(nMinSelection), smooth(nSmooth) {
  assert(minSelection > 0.0);
  double step = 1./nBins;
  for ( unsigned int k = 1; k <= nBins; ++k )
    weights[k*step] = 0.0;
}

void Remapper::fill(double x, double w) {
  map<double,double>::iterator k =
    weights.upper_bound(x);
  assert(k != weights.end());
  k->second += abs(w);
}

void Remapper::finalize() {
  map<double,double> nweights = weights;
  double step = nweights.begin()->first;
  double norm = 0.0;
  for ( map<double,double>::const_iterator k =
	  nweights.begin(); k != nweights.end(); ++k )
    norm += k->second;
  double sum = 0.0;
  for ( map<double,double>::iterator k =
	  nweights.begin(); k != nweights.end(); ++k ) {
    k->second /= norm;
    k->second = max(k->second,minSelection);
    sum += k->second;
  }
  if ( smooth ) {
    assert(nweights.size() >= 2);
    map<double,double> nnweights = nweights;
    nnweights.begin()->second =
      0.5*(nweights.begin()->second + (++nweights.begin())->second);
    (--nnweights.end())->second =
      0.5*((--nweights.end())->second + (--(--nweights.end()))->second);
    sum = nnweights.begin()->second + (--nnweights.end())->second;
    map<double,double>::iterator nb = ++nnweights.begin();
    map<double,double>::const_iterator b = ++nweights.begin();
    for ( ; b != --nweights.end(); ++b, ++nb ) {
      map<double,double>::const_iterator bm = b; --bm;
      map<double,double>::const_iterator bp = b; ++bp;
      nb->second = (bm->second + b->second + bp->second)/3.;
      sum += nb->second;
    }
    nweights = nnweights;
  }
  norm = 0.0;
  for ( map<double,double>::const_iterator k =
	  nweights.begin(); k != nweights.end(); ++k ) {
    norm += k->second;
    SelectorEntry s;
    s.lower = k->first - step;
    s.upper = k->first;
    s.value = k->second / sum / step;
    selector[norm/sum] = s;
  }
}

pair<double,double> Remapper::generate(double r) const {
  map<double,SelectorEntry>::const_iterator bin
    = selector.upper_bound(r);
  // happens, if there is truely non-zero cross section
  // then let the integrator pick up this result downstream
  if ( bin == selector.end() )
    return make_pair(r,1.);
  const SelectorEntry& binInfo = bin->second;
  double intUp = bin->first;
  double intLow = 
    bin != selector.begin() ? (--bin)->first : 0.0;
  double x = 
    ((binInfo.upper-binInfo.lower)/(intUp-intLow))*
    (r-(intLow*binInfo.upper-intUp*binInfo.lower)/(binInfo.upper-binInfo.lower));
  return pair<double,double>(x,binInfo.value);
}

void Remapper::fromXML(const XML::Element& elem) {

  elem.getFromAttribute("MinSelection",minSelection);
  elem.getFromAttribute("Smooth",smooth);
  size_t nbins = 0;
  elem.getFromAttribute("NBins",nbins);

  list<XML::Element>::const_iterator cit;

  cit = elem.findFirst(XML::ElementTypes::Element,"BinData");
  if ( cit == elem.children().end() )
    throw runtime_error("[ExSample::Remapper] Expected a BinData element.");

  const XML::Element& bindata = *cit;
  cit = bindata.findFirst(XML::ElementTypes::ParsedCharacterData,"");
  if ( cit == bindata.children().end() )
    throw runtime_error("[ExSample::Remapper] Expected bin data.");
  istringstream bdata(cit->content());
  for ( size_t k = 0; k < nbins; ++k ) {
    double x,w;
    bdata >> w >> x;
    weights[w] = x;
  }

  cit = elem.findFirst(XML::ElementTypes::Element,"SelectorData");
  if ( cit == elem.children().end() )
    throw runtime_error("[ExSample::Remapper] Expected a SelectorData element.");

  const XML::Element& selectordata = *cit;
  cit = selectordata.findFirst(XML::ElementTypes::ParsedCharacterData,"");
  if ( cit == selectordata.children().end() )
    throw runtime_error("[ExSample::Remapper] Expected selector data.");
  istringstream sdata(cit->content());
  for ( size_t k = 0; k < nbins; ++k ) {
    double w; SelectorEntry s;
    sdata >> w >> s.lower >> s.upper >> s.value;
    selector[w] = s;
  }

}

XML::Element Remapper::toXML() const {

  XML::Element res(XML::ElementTypes::Element,"Remapper");
  res.appendAttribute("MinSelection",minSelection);
  res.appendAttribute("Smooth",smooth);
  res.appendAttribute("NBins",weights.size());

  XML::Element bindata(XML::ElementTypes::Element,"BinData");

  ostringstream bdata;
  bdata << setprecision(17);
  for ( map<double,double>::const_iterator b = weights.begin();
	b != weights.end(); ++b )
    bdata << b->first << " " << b->second << " ";

  XML::Element belem(XML::ElementTypes::ParsedCharacterData,bdata.str());
  bindata.append(belem);
  
  res.append(bindata);

  XML::Element selectordata(XML::ElementTypes::Element,"SelectorData");

  ostringstream sdata;
  sdata << setprecision(17);
  for ( map<double,SelectorEntry>::const_iterator b = selector.begin();
	b != selector.end(); ++b )
    sdata << b->first << " " << b->second.lower << " " 
	  << b->second.upper << " " << b->second.value << " ";

  XML::Element selem(XML::ElementTypes::ParsedCharacterData,sdata.str());
  selectordata.append(selem);
  
  res.append(selectordata);

  return res;

}

inline double sqr(double x) {
  return x*x;
}

void Remapper::test(size_t n, std::ostream& os) {

  double sumwFlat = 0.0;
  double sumw2Flat = 0.0;

  for ( size_t k = 0; k < n; ++k ) {
    double x = drand48();
    double fx = x < 0.7 ? 5.*pow(x,0.4)*pow(0.7-x,2.4) : ( x > 0.8 ? x*x : 0.0 );
    sumwFlat += fx;
    sumw2Flat += fx*fx;
    fill(x,fx);
  }

  finalize();

  Remapper check(weights.size(),0.001,false);
  double sumw = 0.0;
  double sumw2 = 0.0;

  for ( size_t k = 0; k < n; ++k ) {
    double r = drand48();
    pair<double,double> rw = generate(r);
    double x = rw.first;
    double fx = x < 0.7 ? 5.*pow(x,0.4)*pow(0.7-x,2.4) : ( x > 0.8 ? x*x : 0.0 );
    fx /= rw.second;
    sumw += fx;
    sumw2 += fx*fx;
    check.fill(x,1.);
  }

  cerr << setprecision(6) 
       << "int flat   = "
       << (sumwFlat/n) << " +/- " << sqrt(abs(sqr(sumwFlat/n)-sumw2Flat/n)/(n-1.)) << "\n"
       << "int mapped = "
       << (sumw/n) << " +/- " << sqrt(abs(sqr(sumw/n)-sumw2/n)/(n-1.)) << "\n"
       << "int exact  = 0.353848\n" << flush;

  double sum = 0.0;
  double sumCheck = 0.0;
  map<double,double>::const_iterator w = weights.begin();
  map<double,double>::const_iterator wx = check.weights.begin();
  for ( ; w != weights.end(); ++w, ++wx ) {
    sum += w->second;
    sumCheck += wx->second;
  }

  double step = weights.begin()->first;

  map<double,SelectorEntry>::const_iterator s = selector.begin();
  w = weights.begin();
  wx = check.weights.begin();

  for ( ; s != selector.end(); ++s, ++w, ++wx )
    os << s->second.lower << " " << s->second.upper << " "
       << w->second/sum/step << " " << s->second.value << " "
       << s->first << " " << wx->second/sumCheck/step << "\n"
       << flush;

}



