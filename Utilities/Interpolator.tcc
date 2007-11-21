// -*- C++ -*-
//
// Interpolator.tcc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include <cassert>

using namespace Herwig;

template <typename ValT, typename ArgT>
Interpolator<ValT,ArgT>::Interpolator(vector<ValT> f, 
				      vector<ArgT> x, 
				      unsigned int order) 
  : _fun(f.size(),0.0),_xval(x.size(),0.0),_order(order),
    _funit(TypeTraits<ValT>::baseunit), 
    _xunit(TypeTraits<ArgT>::baseunit)
{
  // check the size of the vectors is the same
  if(_fun.size()!=_xval.size())
    throw Exception() << "Interpolator: The size of the vectors containing " 
		      << "the x and function values are different" 
		      << Exception::runerror;
  if(_order<1)
    throw Exception() << "Interpolator: The order of interpolation is too low" 
		      << " using linear interpolation" 
		      << Exception::runerror;
  assert(x.size() == f.size());
  for (size_t i = 0; i < f.size(); ++i) {
    _fun [i] = f[i] / _funit;
    _xval[i] = x[i] / _xunit;
  }
}

template <typename ValT, typename ArgT>
Interpolator<ValT,ArgT>::Interpolator(size_t size, 
				      double f[], ValT funit,
				      double x[], ArgT xunit,
				      unsigned int order)
  : _fun(size,0.0),_xval(size,0.0),_order(order),
    _funit(funit),_xunit(xunit)
{
  // check the size of the vectors is the same
  if(_order<1)
    throw Exception() << "Interpolator: The order of interpolation is too low" 
		      << " using linear interpolation" 
		      << Exception::runerror;
  for (size_t i = 0; i < size; ++i) {
    _fun [i] = f[i];
    _xval[i] = x[i];
  }
}

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

template <typename ValT, typename ArgT>
ClassDescription<Interpolator<ValT,ArgT> > 
Interpolator<ValT,ArgT>::initInterpolator;
// Definition of the static class description member.

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
     &Interpolator<ValT,ArgT>::xfunit, 
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
  vector<double> copyx,copyfun;
  // number of points
  unsigned int npoints(_order+2-_order%2),icopy;
  int iloc(0),i;
  do {
    icopy=mid+iloc;
    if(icopy>isize-1){npoints=mp;}
    else {
      copyx.push_back(_xval[icopy]);
      copyfun.push_back(_fun[icopy]);
    }
    iloc=-iloc;
    if(iloc>=0){++iloc;}
  }
  while(copyx.size()<npoints);
  // do this interpolation
  bool extra(npoints!=mp);
  for(ix=0;ix<m;++ix) {
    if(extra) {
      icopy=m-ix-1;
      copyfun[m+1]=(copyfun[m+1]-copyfun[m-1])/(copyx[m+1]-copyx[icopy]);
    }
    i=m;
    for(iy=ix;iy<m;++iy) {
      icopy=i-ix-1;
      copyfun[i]=(copyfun[i]-copyfun[i-1])/(copyx[i]-copyx[icopy]);
      --i;
    }
  }
  double sum(copyfun[m]);
  if(extra) sum=0.5*(sum+copyfun[m+1]);
  i=m-1;
  for(ix=0;ix<m;++ix) {
    sum=copyfun[i]+(xpoint-copyx[i])*sum;
    --i;
  }
  return sum * _funit;
}
