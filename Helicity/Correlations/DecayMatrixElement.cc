// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayMatrixElement class.
//
// Author: Peter Richardson
//

#include "DecayMatrixElement.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DecayMatrixElement.tcc"
#endif

namespace Herwig {
namespace Helicity {

using namespace ThePEG;
    
DecayMatrixElement::~DecayMatrixElement() {}
    
NoPIOClassDescription<DecayMatrixElement> DecayMatrixElement::initDecayMatrixElement; 
// Definition of the static class description member.
    
void DecayMatrixElement::Init() {
  
  static ClassDocumentation<DecayMatrixElement> documentation
    ("The\\classname{DecayMatrixElement} class is designed to store the helicity "
     "amplitude expression for the matrix element of a decay.");
  
}
    
// calculate the decay matrix for this decay
RhoDMatrix DecayMatrixElement::calculateDMatrix(vector<RhoDMatrix> rhoout)
{
  // vectors for the helicities
  vector<int> ihel1(_outspin.size()+1),ihel2(_outspin.size()+1);
  // rhomatrix to be returned
  RhoDMatrix output(_inspin);
  // make sure that this is zeroed
  for(int ix=0;ix<_inspin;++ix)
    {
      int ixa=ix-_inspin/2;if(_inspin%2==0&&ixa>=0){++ixa;}
      for(int iy=0;iy<_inspin;++iy)
	{
	  int iya=iy-_inspin/2;if(_inspin%2==0&&iya>=0){++iya;}
	  output(ixa,iya)=0.;
	}
    }
  // loop over all helicity components of the matrix element
  // outer loop
  Complex temp;
  for(unsigned int ix=0;ix<_matrixelement.size();++ix)
    {
      // map the vector index to the helicities
      for(int ixa=_outspin.size();ixa>0;--ixa)
	{
	  ihel1[ixa]=(ix%_constants[ixa])/_constants[ixa+1]-int(_outspin[ixa-1]/2);
	  if(_outspin[ixa-1]%2==0&&ihel1[ixa]>=0){++ihel1[ixa];}
	}
      ihel1[0]=(ix%_constants[0])/_constants[1]-int(_inspin/2);
      if(_inspin%2==0&&ihel1[0]>=0){++ihel1[0];}
      // inner loop
      for(unsigned int iy=0;iy<_matrixelement.size();++iy)
	{
	  // map the vector index to the helicities	   
	  for(unsigned int iya=_outspin.size();iya>0;--iya)
	    {
	      ihel2[iya]=(iy%_constants[iya])/_constants[iya+1]
		-int(_outspin[iya-1]/2);
	      if(_outspin[iya-1]%2==0&&ihel2[iya]>=0){++ihel2[iya];}
	    }
	  ihel2[0]=(iy%_constants[0])/_constants[1]-int(_inspin/2);
	  if(_inspin%2==0&&ihel2[0]>=0){++ihel2[0];}
	  // matrix element piece
	  temp=_matrixelement[ix]*conj(_matrixelement[iy]);
	  // spin density matrices for the outgoing particles
	  for(unsigned int iz=0;iz<_outspin.size();++iz)
	    {temp*=rhoout[iz](ihel1[iz+1],ihel2[iz+1]);}
	  output(ihel1[0],ihel2[0])+=temp;
	}
    }
  // ensure unit trace for the matrix
  output.normalize();
  // return the answer
  return output;
}

// calculate the rho matrix for a given outgoing particle
RhoDMatrix DecayMatrixElement::calculateRhoMatrix(int id,RhoDMatrix rhoin,
						  vector<RhoDMatrix>rhoout)
{
  // vectors for the helicities
  vector<int> ihel1(_outspin.size()+1),ihel2(_outspin.size()+1);
  // rhomatrix to be returned
  RhoDMatrix output(_outspin[id]); output.zero();
  // loop over all helicity components of the matrix element
  // outer loop
  Complex temp;
  for(unsigned int ix=0;ix<_matrixelement.size();++ix)
    {
      // map the vector index to the helicities
      for(unsigned int ixa=_outspin.size();ixa>0;--ixa)
	{
	  ihel1[ixa]=(ix%_constants[ixa])/_constants[ixa+1]-int(_outspin[ixa-1]/2);
	  if(_outspin[ixa-1]%2==0&&ihel1[ixa]>=0){++ihel1[ixa];}
	}
      ihel1[0]=(ix%_constants[0])/_constants[1]-int(_inspin/2);
      if(_inspin%2==0&&ihel1[0]>=0){++ihel1[0];}
      // inner loop
      for(unsigned int iy=0;iy<_matrixelement.size();++iy)
	{
	  // map the vector index to the helicities	   
	  for(unsigned int iya=_outspin.size();iya>0;--iya)
	    {
	      ihel2[iya]=(iy%_constants[iya])/_constants[iya+1]
		-int(_outspin[iya-1]/2);
	      if(_outspin[iya-1]%2==0&&ihel2[iya]>=0){++ihel2[iya];}
	    }
	  ihel2[0]=(iy%_constants[0])/_constants[1]-int(_inspin/2);
	  if(_inspin%2==0&&ihel2[0]>=0){++ihel2[0];}
	  // matrix element piece
	  temp=_matrixelement[ix]*conj(_matrixelement[iy]);
	  // spin denisty matrix for the incoming particle
	  temp*=rhoin(ihel1[0],ihel2[0]);
	  // spin density matrix for the outgoing particles
	  for(unsigned int iz=0;iz<_outspin.size()-1;++iz)
	    {
	      if(int(iz)<id){temp*=rhoout[iz](ihel1[iz+1],ihel2[iz+1]);}
	      else{temp*=rhoout[iz](ihel1[iz+2],ihel2[iz+2]);}
	    }
	  // add to the rho matrix
	  output(ihel1[id+1],ihel2[id+1])+=temp;
	}
    }
  // return the answer
  output.normalize();
  return output;
}

}
}

