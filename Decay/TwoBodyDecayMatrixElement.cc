// -*- C++ -*-
//
// TwoBodyMatrixElement.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoBodyMatrixElement class.
//
// Author: Peter Richardson
//

#include "TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG;
    
// calculate the decay matrix for this decay
RhoDMatrix TwoBodyDecayMatrixElement::
calculateDMatrix(const vector<RhoDMatrix> & rhoout) const {
  // rhomatrix to be returned
  RhoDMatrix output(inspin(), false);
  // loop over all helicity components of the matrix element
  for(int ihel1=0;ihel1<outspin()[0];++ihel1) {
    for(int jhel1=0;jhel1<outspin()[0];++jhel1) {
      if(rhoout[0](ihel1,jhel1)==0.)continue;
      for(int ihel2=0;ihel2<outspin()[1];++ihel2) {
        for(int jhel2=0;jhel2<outspin()[1];++jhel2) {
          if(rhoout[1](ihel2,jhel2)==0.)continue;
  	  for(int ihel0=0;ihel0<inspin();++ihel0) {
            for(int jhel0=0;jhel0<inspin();++jhel0) {
      
	      output(ihel0,jhel0) += 
		     matrixElement_[ihel0][ihel1][ihel2]*
		conj(matrixElement_[jhel0][jhel1][jhel2])*
		rhoout[0](ihel1,jhel1)*rhoout[1](ihel2,jhel2);
	    }
	  }
	}
      }
    }
  }
  // ensure unit trace for the matrix
  output.normalize();
  // return the answer
  return output;
}

// calculate the rho matrix for a given outgoing particle
RhoDMatrix TwoBodyDecayMatrixElement::
calculateRhoMatrix(int id,const RhoDMatrix & rhoin,
		   const vector<RhoDMatrix> & rhoout) const {
  // rhomatrix to be returned
  RhoDMatrix output(outspin()[id], false);
  // loop over all helicity components of the matrix element
  if(id==0) {
    for(int ihel0=0;ihel0<inspin();++ihel0) {
      for(int jhel0=0;jhel0<inspin();++jhel0) {
        if (rhoin(ihel0,jhel0)==0.)continue;
        for(int ihel2=0;ihel2<outspin()[1];++ihel2) {
          for(int jhel2=0;jhel2<outspin()[1];++jhel2) {
            if(rhoout[0](ihel2,jhel2)==0.)continue;
            for(int ihel1=0;ihel1<outspin()[0];++ihel1) {
              if(matrixElement_[ihel0][ihel1][ihel2]==0.)continue;
              for(int jhel1=0;jhel1<outspin()[0];++jhel1) {
                output(ihel1,jhel1) +=
                matrixElement_[ihel0][ihel1][ihel2]*
                conj(matrixElement_[jhel0][jhel1][jhel2])*
                rhoin(ihel0,jhel0)*
                rhoout[0](ihel2,jhel2);
              }
            }
          }
        }
      }
    }
  }
  else {
    for(int ihel0=0;ihel0<inspin();++ihel0) {
      for(int jhel0=0;jhel0<inspin();++jhel0) {
        if (rhoin(ihel0,jhel0)==0.)continue;
	for(int ihel1=0;ihel1<outspin()[0];++ihel1) {
	  for(int jhel1=0;jhel1<outspin()[0];++jhel1) {
            if (rhoout[0](ihel1,jhel1)==0.)continue;
	    for(int ihel2=0;ihel2<outspin()[1];++ihel2) {
	      for(int jhel2=0;jhel2<outspin()[1];++jhel2) {
		output(ihel2,jhel2) += 
		  matrixElement_[ihel0][ihel1][ihel2]*
		  conj(matrixElement_[jhel0][jhel1][jhel2])*
		  rhoin(ihel0,jhel0)*
		  rhoout[0](ihel1,jhel1);
	      }
	    }
	  }
	}
      }
    }
  }
  // return the answer
  output.normalize();
  return output;
}

// contract the matrix element with the rho matrix of the incoming particle
Complex TwoBodyDecayMatrixElement::contract(const RhoDMatrix & in) const {
  Complex me=0.;
  for(int ihel0=0;ihel0<inspin();++ihel0) {
    for(int jhel0=0;jhel0<inspin();++jhel0) {
      for(int ihel1=0;ihel1<outspin()[0];++ihel1) {
	  for(int ihel2=0;ihel2<outspin()[1];++ihel2) {
	      me +=  matrixElement_[ihel0][ihel1][ihel2]*
		conj(matrixElement_[jhel0][ihel1][ihel2])*in(ihel0,jhel0);
	  }
      }
    }
  }
  return me;
}

// contract the matrix element with the rho matrix of the incoming particle
Complex TwoBodyDecayMatrixElement::contract(const TwoBodyDecayMatrixElement & con, 
					    const RhoDMatrix & in) {
  Complex me=0.;
  for(int ihel0=0;ihel0<inspin();++ihel0) {
    for(int jhel0=0;jhel0<inspin();++jhel0) {
      for(int ihel1=0;ihel1<outspin()[0];++ihel1) {
	for(int ihel2=0;ihel2<outspin()[1];++ihel2) {
	  // compute the term
	  me += matrixElement_[ihel0][ihel1][ihel2]*
	    conj(con.matrixElement_[jhel0][ihel1][ihel2])*in(ihel0,jhel0);
	}
      }
    }
  }
  return me;
}
