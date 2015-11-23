// -*- C++ -*-
//
// Interpolator.tcc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Interpolator class.
//

#include "Interpolator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

template <typename ValT, typename ArgT>
void Interpolator<ValT,ArgT>::persistentOutput(PersistentOStream & os) const {
  os << _xval << _fun << _order 
     << ounit(_funit,TypeTraits<ValT>::baseunit) 
     << ounit(_xunit,TypeTraits<ArgT>::baseunit);
}

template <typename ValT, typename ArgT>
void Interpolator<ValT,ArgT>::persistentInput(PersistentIStream & is, int) {
  is >> _xval >> _fun >> _order 
     >> iunit(_funit,TypeTraits<ValT>::baseunit) 
     >> iunit(_xunit,TypeTraits<ArgT>::baseunit);
}

/**
 * Macro for Interpolator user classes to use. Only they know what the
 * template arguments are going to be.
 */
#define HERWIG_INTERPOLATOR_CLASSDESC(Name,ValT,ArgT) \
/**                                                   \
* Register the Interpolator with ThePEG               \
*/                                                    \
DescribeClass<Interpolator<ValT,ArgT>,Interfaced>     \
describeHerwigInterpolatorFor##Name("Herwig::Interpolator<"#ValT","#ArgT">","");\



template <typename ValT, typename ArgT>
void Interpolator<ValT,ArgT>::Init() {

  static ClassDocumentation<Interpolator<ValT,ArgT> > documentation
    ("The Interpolator class is design to interpolate a table of values");

  static Parameter<Interpolator<ValT,ArgT>,unsigned int> interfaceOrder
    ("Order",
     "Order of the interpolation",
     &Interpolator::_order, 3, 1, 10,
     false, false, Interface::limited);

  static ParVector<Interpolator<ValT,ArgT>,double> interfaceXValues
    ("XValues",
     "The x values for the interpolation",
     &Interpolator::_xval, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<Interpolator<ValT,ArgT>,double> interfaceFunctionValues
    ("FunctionValues",
     "The function values for the interpolation",
     &Interpolator::_fun, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static Parameter<Interpolator<ValT,ArgT>,ValT> interfaceValueType
    ("ValueType",
     "The unit of the function values",
     &Interpolator<ValT,ArgT>::_funit, 
     TypeTraits<ValT>::baseunit, 
     1.0*TypeTraits<ValT>::baseunit, 
     0*TypeTraits<ValT>::baseunit, 
     0*TypeTraits<ValT>::baseunit,
     false, true, Interface::nolimits);

  static Parameter<Interpolator<ValT,ArgT>,ArgT> interfaceArgType
    ("ArgType",
     "The unit of the function arguments",
     &Interpolator<ValT,ArgT>::_xunit, 
     TypeTraits<ArgT>::baseunit, 
     1.0*TypeTraits<ArgT>::baseunit, 
     0*TypeTraits<ArgT>::baseunit, 
     0*TypeTraits<ArgT>::baseunit,
     false, true, Interface::nolimits);

}

template <typename ValT, typename ArgT>
ValT Interpolator<ValT,ArgT>::operator ()(ArgT xpt) const {
  const double xpoint = xpt / _xunit;
  // size of the vectors
  unsigned int isize(_xval.size());
  // workout the numer of points we need
  unsigned int m(std::min(_order,isize)),mp(m+1),ix,iy;
  // search for the point if the function increases
  int mid,iupp=isize,ilow=0;
  if(_xval[0]>_xval[_xval.size()-1]) {
    do {
      mid=(iupp+ilow)/2;
      if(xpoint>_xval[mid]){iupp=mid;}
      else{ilow=mid;}
    }
    while(iupp-ilow>1);
  }
  // search for the point if the function decreases 
  else {
    do {
      mid=(iupp+ilow)/2;
      if(xpoint<_xval[mid]){iupp=mid;}
      else{ilow=mid;}
    }
    while(iupp-ilow>1);
  }
  // ilow is now the midpoint
  mid=ilow;
  // copy the re-ordered interpolation points 
  // number of points
  unsigned int npoints(_order+2-_order%2),icopy,j(0);
  int iloc(0),i(0);
  do {
    icopy=mid+iloc;
    if(icopy>isize-1) npoints=mp;
    else {
      _copyx[j]   = _xval[icopy];
      _copyfun[j] = _fun [icopy];
      ++j;
    }
    iloc=-iloc;
    if(iloc>=0){++iloc;}
  }
  while(j<npoints);
  // do this interpolation
  bool extra(npoints!=mp);
  for(ix=0;ix<m;++ix) {
    if(extra) {
      icopy=m-ix-1;
      _copyfun[m+1]=(_copyfun[m+1]-_copyfun[m-1])/(_copyx[m+1]-_copyx[icopy]);
    }
    i=m;
    for(iy=ix;iy<m;++iy) {
      icopy=i-ix-1;
      _copyfun[i]=(_copyfun[i]-_copyfun[i-1])/(_copyx[i]-_copyx[icopy]);
      --i;
    }
  }
  double sum(_copyfun[m]);
  if(extra) sum=0.5*(sum+_copyfun[m+1]);
  i=m-1;
  for(ix=0;ix<m;++ix) {
    sum=_copyfun[i]+(xpoint-_copyx[i])*sum;
    --i;
  }
  return sum * _funit;
}
