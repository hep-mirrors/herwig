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

#include "Interpolator2d.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace std;

void Interpolator2d::persistentOutput(PersistentOStream & os) const {
}

void Interpolator2d::persistentInput(PersistentIStream & is, int) {
}

//Interpolator2d::initInterpolator2d;
// Definition of the static class description member.

void Interpolator2d::Init() {
  static ClassDocumentation< Interpolator2d > documentation
    ("The Interpolator2d class is designed to interpolate a 2d table of values");
}

double Interpolator2d::operator ()(double x1, double x2) const {
  //find which grid (i,j) we are in and get the coefficients of that ij
  int i = int( ( x1 - _x1_points[0] ) / _dx1 );
  int j = int( ( x2 - _x2_points[0] ) / _dx2 );

  //grid square is labelled by i,j and include i+1,j i+1,j+1, i,j+1
  //so there are i-1 * j-1 squares

  if( i > _the_weights.size() - 2 || i <  0 
      || j > _the_weights[0].size() - 2 || j < 0 ){
    cerr<< "out of range at ("<<x1<<","<<x2<<") at ("<<i<<","<<j<<")\n";
    cerr<<"weights x size = "<< _the_weights.size()<<"\n"; 
    cerr<<"weights y size = "<< _the_weights[0].size()<<"\n"; 
    cerr<<"x1 size = "<<_x1_points.size()<<"\n";
    cerr<<"x2 size = "<<_x2_points.size()<<"\n";
    cerr<<"dx1 = "<< _dx1 <<"\n";
    cerr<<"dx2 = "<< _dx2 <<"\n";
    cerr<<"x1 min = "<<  _x1_points[0]<<", x1 max = "<<_x1_points[ _x1_points.size() - 1]<<"\n";
    cerr<<"x2 min = "<<  _x2_points[0]<<", x2 max = "<<_x2_points[ _x2_points.size() - 1]<<"\n";

    return 0.;
  }
  
  double t;
  double u;
  //find co-ordinate parameters t and u
  if( i == _the_weights.size() - 1 )
    t = ( x1 - _x1_points[i] )/ ( 0. -  _x1_points[i] );
  else
    t = ( x1 - _x1_points[i] )/ ( _x1_points[i+1] -  _x1_points[i] );

  if( i == _the_weights[0].size() - 1 )
    u = ( x2 - _x2_points[j] )/ ( 0. -  _x2_points[j] );
  else
    u = ( x2 - _x2_points[j] )/ ( _x2_points[j+1] -  _x2_points[j] );
  
  //return c_kl t^k u^l
  double result = 0.;
  //
  for(int k = 0; k < 4; k++ ){
    for( int l = 0; l < 4; l++ )
      result += _coefficients[i][j][k][l] * pow( t, k ) * pow( u, l );
  }
  return result;
}


void Interpolator2d::gridHists(){
  ofstream gridOut( "grid.out" );

   //make a test histogram
  for( unsigned int ix = 0; ix < _x1_points.size(); ++ix ) {
    gridOut <<"NEW FRAME \nSET WINDOW X 1.6 8 Y 3.5 9\n"
	    <<"SET FONT DUPLEX\n"
	    <<"SET ORDER X Y DX \n";
    for( unsigned int jx = 0; jx <  _x2_points.size(); ++jx ){
      gridOut<< _x2_points[jx] <<"\t"<< _the_weights[ix][jx] <<"\t"<<_dx2/2.<<"\n";
    }
    gridOut<<"HIST BLACK \n";
  }
  gridOut.close();

}


Interpolator2d::Interpolator2d( const vector< vector< double > > & the_weights, 
				const vector< double > & x1_points,
				const vector< double > & x2_points ){
  _the_weights = the_weights;
  _x1_points = x1_points;
  _x2_points = x2_points; 
  _dx1 = x1_points[1] - x1_points[0];
  _dx2 = x2_points[1] - x2_points[0];

  //create the magic matrix
  int wt_d[16][16] =
    { { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
      {-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0 },
      { 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 },
      { 0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1 },
      { 0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1 },
      {-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0 },
      { 9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2 },
      {-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2 },
      { 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0 },
      {-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1 },
      { 4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1 } };
  //convert this to a vector - easier and better
  //initialise
  vector< int > dummy( 16, 0 );
  vector< vector< int > > wt( 16, dummy );
  _wt = wt;  
  //fill the 2d vector with elements from the array
  for( int i = 0; i < 16; i++ ){
    for( int j = 0; j < 16; j++ ){
      _wt[i][j] = wt_d[i][j];
    }
  } 

  //fill the derivative matrices - note edges are not filled
  //the edge points have all derivatives set to zero.
  _theDerivX1 = findDerivX1();
  _theDerivX2 = findDerivX2();
  _theDerivX1X2 = findDerivX1X2();

  //loop over all grid squares and find the Cij coefficients of each
  //the grid square index is determined by the bottom left point
  //so the grid square indices go from i (j) = 0 -> m-1 (n-1) 
  //initialise coefficients
  vector< double > dummy1( 4, 0. );
  vector< vector< double > > dummyMatrix( 4, dummy1 );
  vector< vector< vector< double > > > dummy2( _the_weights[0].size(), dummyMatrix );
  vector< vector< vector< vector< double > > > > 
    coefficients( _the_weights.size(), dummy2 );
  
  _coefficients = coefficients;

  //fill coefficients
  for( int i = 0; i < _the_weights.size() - 1; i++ ){
    for( int j = 0; j < _the_weights[0].size() - 1; j++ ){
      _coefficients[i][j] = getCoefficient( i, j );
    }
  }
}


vector< vector< double > > Interpolator2d::getCoefficient( unsigned int i, unsigned int j ){
  //initialise the coefficient 4*4 matrix that is to be returned
  vector< double > dummy1( 4, 0 );
  vector< vector< double > > theCoeff( 4, dummy1 );
  
  //fill the y vals and derivs for grid square
  vector< double > y;
  y.push_back( _the_weights[i][j] );
  y.push_back( _the_weights[i+1][j] );
  y.push_back( _the_weights[i+1][j+1] );
  y.push_back( _the_weights[i][j+1] );
  //fill derivs of grid square
  vector< double > yx1;
  yx1.push_back( _theDerivX1[i][j] );
  yx1.push_back( _theDerivX1[i+1][j] );
  yx1.push_back( _theDerivX1[i+1][j+1] );
  yx1.push_back( _theDerivX1[i][j+1] );
  vector< double > yx2;
  yx2.push_back( _theDerivX2[i][j] );
  yx2.push_back( _theDerivX2[i+1][j] );
  yx2.push_back( _theDerivX2[i+1][j+1] );
  yx2.push_back( _theDerivX2[i][j+1] );
  vector< double > yx12;
  yx12.push_back( _theDerivX1X2[i][j] );
  yx12.push_back( _theDerivX1X2[i+1][j] );
  yx12.push_back( _theDerivX1X2[i+1][j+1] );
  yx12.push_back( _theDerivX1X2[i][j+1] );
  
  //make temporary vectors used in calc
  vector< double > x(16), cl(16);
  for( int i = 0; i < 4; i++ ){
    x[i] = y[i];
    x[i+4] = yx1[i] * _dx1;
    x[i+8] = yx2[i] * _dx2;
    x[i+12] = yx12[i] * _dx1 * _dx2;
  }
  //matrix multiply by wt[i][j] store result in cl
  for(int i = 0; i < 16; i++ ){
    double element = 0.;
    for( int j = 0; j < 16; j++ )
      element += _wt[i][j]*x[j];
    cl[i] = element;
  }
  int l = 0;
  for( int i = 0; i < 4; i++ ){
    for( int j = 0; j < 4; j++ ) {
      theCoeff[i][j] = cl[l++];
    }
  }
  return theCoeff;
}

vector< vector< double > > Interpolator2d::findDerivX1(){
  vector< double > dummy( _the_weights[0].size(), 0. );
  vector< vector< double > > deriv( _the_weights.size(), dummy );
  for( unsigned int i = 0; i < _the_weights.size(); i++ ){
    for( unsigned int j = 0; j < _the_weights[i].size(); j++ ){
      //special cases - use less accurate approxs if on the edge of limits
      if( i == 0 ){
	deriv[i][j] = ( _the_weights[i+1][j] - _the_weights[i][j] ) 
	  / _dx1;
      }
      else if( i == _the_weights.size() - 1 ){
	deriv[i][j] = ( _the_weights[i][j] - _the_weights[i-1][j] ) 
	  / _dx1;
      }
      else{
	deriv[i][j] = ( _the_weights[i+1][j] - _the_weights[i-1][j] ) 
	  / 2. / _dx1;		  
      }
    }
  }
  return deriv;
}

vector< vector< double > > Interpolator2d::findDerivX2(){
  vector< double > dummy( _the_weights[0].size(), 0. );
  vector< vector< double > > deriv( _the_weights.size(), dummy );
  //loop over all of grid but not edges (cannot find derivs here
  for( unsigned int i = 0; i < _the_weights.size(); i++ ){
    for( unsigned int j = 0; j < _the_weights[i].size(); j++ ){
      //special cases - use less accurate approxs if on the edge of limits
      if( j == 0 ){
	deriv[i][j] = ( _the_weights[i][j+1] - _the_weights[i][j] ) 
	  / _dx2;
      }
      else if( j == _the_weights[i].size() - 1 ){
	deriv[i][j] = ( _the_weights[i][j] - _the_weights[i][j-1] ) 
	  / _dx2;
      }
      else{
	deriv[i][j] = ( _the_weights[i][j+1] - _the_weights[i][j-1] ) 
	  / 2. / _dx2;      }
    }
  }
  return deriv;
}

vector< vector< double > > Interpolator2d::findDerivX1X2(){
  vector< double > dummy( _the_weights[0].size(), 0. );
  vector< vector< double > > deriv( _the_weights.size(), dummy );
  //loop over all of grid but not edges (cannot find derivs here
  for( unsigned int i = 0; i < _the_weights.size(); i++ ){
    for( unsigned int j = 0; j < _the_weights[i].size(); j++ ){
       //special cases - use less accurate approxs if on the edge of limits
      if( i == 0 ){
	if( j != _the_weights.size()-1 ){
	  deriv[i][j] = - ( _the_weights[i+1][j-1] - _the_weights[i][j] ) 
	    / _dx1 / _dx2;
	}
	else{
	deriv[i][j] = ( _the_weights[i+1][j+1] - _the_weights[i][j] ) 
	  / _dx1 / _dx2;
	}
      } 
	else if( i == _the_weights.size()-1 ){
	  if( j != 0 ){
	  deriv[i][j] = ( _the_weights[i][j] - _the_weights[i-1][j-1] ) 
	    / _dx1 / _dx2;
	}
	else{
	  //is the sign correct here??
	  deriv[i][j] = - ( _the_weights[i][j] - _the_weights[i-1][j+1] ) 
	    / _dx1 / _dx2;
	}
      }
      else{
	deriv[i][j] = ( _the_weights[i+1][j+1] - _the_weights[i+1][j-1] -
			_the_weights[i-1][j+1] + _the_weights[i-1][j-1] ) 
	  / 2. / _dx1 / 2. / _dx2;
      }
    }
  }
  return deriv;
}
