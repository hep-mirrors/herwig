// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Interpolator class.
//
//  Author: Peter Richardson
//

#include "Interpolator.h"
#include "ThePEG/Interface/ClassDocumentation.h"


namespace Herwig{
  
using namespace Genfun;
using std::endl;

FUNCTION_OBJECT_IMP(Interpolator)
  
Interpolator::Interpolator(const Interpolator & right) 
{  }

Interpolator::Interpolator(vector<double> f, vector<double> x, int order)
{
  //
  _xval=x;
  _fun=f;
  _order=order;
  // check the size of the vectors is the same
  if(_fun.size()!=_xval.size())
    {
      std::cerr << "Interpolator: The size of the vectors containing " 
		<< "the x and function values are different" << std::endl;
      // make them the same (use the smallest)
      if(_fun.size()<_xval.size())
	{_xval.resize(_fun.size());}
      else if(_xval.size()<_fun.size())
	{_fun.resize(_xval.size());}
    }
  if(_order<1)
    {
      std::cerr << "Interpolator: The order of interpolation is too low" 
		<< " using linear interpolation" << std::endl;
      _order=1;
    }
}

// destructor
Interpolator::~Interpolator() {}
  
double Interpolator::operator ()(double xpoint) const {
  // size of the vectors
  unsigned int isize(_xval.size());
  // workout the numer of points we need
  unsigned int m(std::min(_order,isize)),mp(m+1),ix,iy;
  // search for the point if the function increases
  int mid,iupp=isize,ilow=0;
  if(_xval[0]>_xval[_xval.size()-1])
    {
      do
	{
	  mid=(iupp+ilow)/2;
	  if(xpoint>_xval[mid]){iupp=mid;}
	  else{ilow=mid;}
	}
      while(iupp-ilow>1);
    }
  // search for the point if the function decreases 
  else
    {
      do
	{
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
  do
    {
      icopy=mid+iloc;
      if(icopy>isize-1){npoints=mp;}
      else
	{
	  copyx.push_back(_xval[icopy]);
	  copyfun.push_back(_fun[icopy]);
	}
      iloc=-iloc;
      if(iloc>=0){++iloc;}
    }
  while(copyx.size()<npoints);
  // do this interpolation
  bool extra(npoints!=mp);
  for(ix=0;ix<m;++ix)
    {
      if(extra)
	{
	  icopy=m-ix-1;
	  copyfun[m+1]=(copyfun[m+1]-copyfun[m-1])/(copyx[m+1]-copyx[icopy]);
	}
      i=m;
      for(iy=ix;iy<m;++iy)
	{
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
  return sum;
}

}
