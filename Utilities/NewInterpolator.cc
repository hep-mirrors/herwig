// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NewInterpolator class.
//

#include "NewInterpolator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

NewInterpolator::NewInterpolator(vector<double> f, vector<double> x, int order) 
  : _xval(x),_fun(f),_order(order)
{
  // check the size of the vectors is the same
  if(_fun.size()!=_xval.size())
    throw Exception() << "NewInterpolator: The size of the vectors containing " 
		      << "the x and function values are different" 
		      << Exception::runerror;
  if(_order<1)
    throw Exception() << "NewInterpolator: The order of interpolation is too low" 
		      << " using linear interpolation" 
		      << Exception::runerror;
}

void NewInterpolator::persistentOutput(PersistentOStream & os) const {
  os << _xval << _fun << _order;
}

void NewInterpolator::persistentInput(PersistentIStream & is, int) {
  is >> _xval >> _fun >> _order;
}

ClassDescription<NewInterpolator> NewInterpolator::initNewInterpolator;
// Definition of the static class description member.

void NewInterpolator::Init() {

  static ClassDocumentation<NewInterpolator> documentation
    ("The NewInterpolator class is design to interpolate a table of values");

  static Parameter<NewInterpolator,unsigned int> interfaceOrder
    ("Order",
     "Order of the interpolation",
     &NewInterpolator::_order, 3, 1, 10,
     false, false, Interface::limited);

  static ParVector<NewInterpolator,double> interfaceXValues
    ("XValues",
     "The x values for the interpolation",
     &NewInterpolator::_xval, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<NewInterpolator,double> interfaceFunctionValues
    ("FunctionValues",
     "The function values for the interpolation",
     &NewInterpolator::_fun, -1, 0., 0, 0,
     false, false, Interface::nolimits);
}

double NewInterpolator::operator ()(double xpoint) const {
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
      icopy=m-ix;
      copyfun[m+1]=(copyfun[m+1]-copyfun[m-1])/(copyx[m+1]-copyx[m-1]);
    }
    i=m;
    for(iy=ix;iy<m;++iy) {
      icopy=i-1;
      copyfun[i]=(copyfun[i]-copyfun[i-1])/(copyx[i]-copyx[icopy]);
      --i;
    }
  }
  double sum(copyfun[m]);
  if(extra) sum=0.5*(sum+copyfun[m+1]);
  i=m-1;
  for(ix=0;ix<m;++ix) sum=copyfun[i]+(xpoint-copyx[i])*sum;--i;
  return sum;
}
