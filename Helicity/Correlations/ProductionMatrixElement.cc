// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ProductionMatrixElement class.
//
// Author: Peter Richardson
//

#include "ProductionMatrixElement.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ProductionMatrixElement.tcc"
#endif

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
ProductionMatrixElement::~ProductionMatrixElement() {}
    
NoPIOClassDescription<ProductionMatrixElement> ProductionMatrixElement::initProductionMatrixElement;
// Definition of the static class description member.
    
void ProductionMatrixElement::Init() {
  
  static ClassDocumentation<ProductionMatrixElement> documentation
    ("The \\classname{ProductionMatrixElement} class is designed to store the "
     "matrix element for a hard interaction.");
  
}
  
// calculate a decay matrix for one of the incoming particles
RhoDMatrix ProductionMatrixElement::calculateDMatrix(int id, RhoDMatrix rhoin,
						     vector<RhoDMatrix> rhoout)
{
  // vectors for the helicities
  vector<unsigned int> ihel1(_outspin.size()+2),ihel2(_outspin.size()+2);
  // rhomatrix to be returned
  RhoDMatrix output(_inspin[id]); output.zero();
  // loop over all helicity components of the matrix element
  // outer loop
  Complex temp;
  unsigned int ix,iy;
  int ixa,iya;
  for(ix=0;ix<_matrixelement.size();++ix)
    {
      // map the vector index to the helicities
      for(ixa=_outspin.size()+1;ixa>=0;--ixa)
	{ihel1[ixa]=(ix%_constants[ixa])/_constants[ixa+1];}
      // inner loop
      for(iy=0;iy<_matrixelement.size();++iy)
	{
	  // map the vector index to the helicities	   
	  for(iya=_outspin.size()+1;iya>=0;--iya)
	    {ihel2[iya]=(iy%_constants[iya])/_constants[iya+1];}
	  // matrix element piece
	  temp=_matrixelement[ix]*conj(_matrixelement[iy]);
	  // spin density matrices for the outgoing particles
	  for(unsigned int iz=0;iz<_outspin.size();++iz)
	    {temp*=rhoout[iz](ihel1[iz+2],ihel2[iz+2]);}
	  // construct the spin density matrix
	  if(id==0)
	    {
	      temp*=rhoin(ihel1[1],ihel2[1]);
	      output(ihel1[0],ihel2[0])+=temp;
	    }
	  else
	    {
	      temp*=rhoin(ihel1[0],ihel2[0]);
	      output(ihel1[1],ihel2[1])+=temp;
	    }
	}
    }
  // return the answer
  return output;
}

// calculate the rho matrix for a given outgoing particle
RhoDMatrix ProductionMatrixElement::calculateRhoMatrix(int id,
						       RhoDMatrix rhoin0,
						       RhoDMatrix rhoin1,
						       vector<RhoDMatrix>rhoout)
{
  unsigned int ix,iy;
  int ixa,iya;
  // vectors for the helicities
  vector<unsigned int> ihel1(_outspin.size()+2),ihel2(_outspin.size()+2);
  // rhomatrix to be returned
  RhoDMatrix output(_outspin[id]); output.zero();
  // loop over all helicity components of the matrix element
  // outer loop
  Complex temp;
  for(ix=0;ix<_matrixelement.size();++ix)
    {
      // map the vector index to the helicities
      for(ixa=_outspin.size()+1;ixa>=0;--ixa)
	{ihel1[ixa]=(ix%_constants[ixa])/_constants[ixa+1];}
      // inner loop
      for(iy=0;iy<_matrixelement.size();++iy)
	{
	  // map the vector index to the helicities	   
	  for(iya=_outspin.size()+1;iya>=0;--iya)
	    {ihel2[iya]=(iy%_constants[iya])/_constants[iya+1];}
	  // matrix element piece
	  temp=_matrixelement[ix]*conj(_matrixelement[iy]);
	  // spin denisty matrix for the incoming particles
	  temp*=rhoin0(ihel1[0],ihel2[0]);
	  temp*=rhoin1(ihel1[1],ihel2[1]);
	  // spin density matrix for the outgoing particles
	  for(unsigned int iz=0;iz<_outspin.size()-1;++iz)
	    {
	      if(int(iz)<id){temp*=rhoout[iz](ihel1[iz+2],ihel2[iz+2]);}
	      else{temp*=rhoout[iz](ihel1[iz+3],ihel2[iz+3]);}
	    }
	  output(ihel1[id+2],ihel2[id+2])+=temp;
	}
    }      
  // normalise the matrix so it has unit trace
  output.normalize();
  // return the answer
  return output;
}

}
}

