// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2vv class.
//

#include "MEff2vv.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

using namespace Herwig;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::SpinorWaveFunction;
using Herwig::Helicity::SpinorBarWaveFunction;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::TensorWaveFunction;
using Herwig::Helicity::FFVVertexPtr;
using Herwig::Helicity::FFTVertexPtr;
using Herwig::Helicity::FFSVertexPtr;
using Herwig::Helicity::VVSVertexPtr;
using Herwig::Helicity::GeneralSVVVertexPtr;
using Herwig::Helicity::VVTVertex;
using Herwig::Helicity::VVVVertexPtr;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

typedef Ptr<VVTVertex>::pointer VVTVertexPtr;

double MEff2vv::me2() const {
  //Define wavefunctions
  SpinorWaveFunction spIn(meMomenta()[0], mePartonData()[0], incoming);
  SpinorBarWaveFunction spbIn(meMomenta()[1], mePartonData()[1], incoming);
  VectorWaveFunction vec1(meMomenta()[2], mePartonData()[2], outgoing);
  VectorWaveFunction vec2(meMomenta()[3], mePartonData()[3], outgoing);
  //Massless vectors?
  bool masslC  = (mePartonData()[2]->mass() == 0.);
  bool masslD  = (mePartonData()[3]->mass() == 0.);
  //Define vectors to store diagrams and square elements
  const HPCount ndiags(numberOfDiags());
  const size_t ncf(numberOfFlows());
  vector<Complex> diag(ndiags, Complex(0.));
  vector<double> me(ndiags, 0.);
  double full_me(0.);
  const Energy2 m2 = scale();
  const vector<vector<double> > cfactors = getColourFactors();
  //Intermediate wavefunctions
  ScalarWaveFunction interS; SpinorWaveFunction interF;
  VectorWaveFunction interV; TensorWaveFunction interT;
  //Loop over helicities and calculate diagrams
  for(unsigned int fhel1 = 0; fhel1 < 2; ++fhel1) {
    spIn.reset(fhel1);
    for(unsigned int fhel2 = 0; fhel2 < 2; ++fhel2) {
      spbIn.reset(fhel2);
      for(unsigned int vhel1 = 0; vhel1 < 3; ++vhel1) {
	if(vhel1 == 1 && masslC) ++vhel1;
	vec1.reset(vhel1);
	for(unsigned int vhel2 = 0; vhel2 < 3; ++vhel2) {
	  if(vhel2 == 1 && masslD) ++vhel2;
	  vec2.reset(vhel2);
	  //loop and calculate diagrams
	  vector<Complex> flows = vector<Complex>(2, Complex(0.));
	  for(HPCount ix = 0; ix < ndiags; ++ix) {
	    HPDiagram current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel) {
	      if(current.intermediate->iSpin() == PDT::Spin1Half) {
		if(current.ordered.second) {
		  interF = theFerm[ix].first->evaluate(m2, 3, offshell, spIn, 
						       vec1);
		  diag[ix] = theFerm[ix].second->evaluate(m2, interF, spbIn,
							  vec2);
		}
		else {
		  interF = theFerm[ix].first->evaluate(m2, 3, offshell, spIn, 
						       vec2);
		  diag[ix] = theFerm[ix].second->evaluate(m2, interF, spbIn,  
							  vec1);
		}
	      }	      
	    }
	    else {
	      if(current.intermediate->iSpin() == PDT::Spin0) {
		if(masslC && masslD) {
		  interS = theSca1[ix].first->evaluate(m2, 1, offshell, 
						       spIn, spbIn);
		  diag[ix] = theSca1[ix].second->evaluate(m2, interS, vec1, vec2);
 		}
		else {
		  interS = theSca2[ix].first->evaluate(m2, 1, offshell, 
						       spIn, spbIn);
		  diag[ix] = theSca2[ix].second->evaluate(m2, vec1, vec2, interS);
		}
	      }
	      else if(current.intermediate->iSpin() == PDT::Spin1) {
		interV = theVec[ix].first->evaluate(m2, 1, offshell, spIn, spbIn);
		diag[ix] = theVec[ix].second->evaluate(m2, vec1, vec2, interV);
	      }
	      else { 
		interT = theTen[ix].first->evaluate(m2, 1, offshell, spIn, spbIn);
		diag[ix] = theTen[ix].second->evaluate(m2, vec1, vec2, interT);
	      }
	    }
	    //compute individual diagram squared and save for diagram sel
	    me[ix] += norm(diag[ix]);
	    //find out which colourflow this diag belongs to and add to it
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy)
	      flows[current.colourFlow[iy].first - 1] +=
		current.colourFlow[iy].second*diag[ix];

	  }//end-of-diagram loop
	  //Add results, with appropriate colour factors to full_me
	  for(size_t ii = 0; ii < ncf; ++ii)
	    for(size_t ij = 0; ij < ncf; ++ij)
	      full_me += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();

	}
      }
    }
  }
  const double cAvg = mePartonData()[0]->coloured() ? 1./9. : 1.;
  const double ifact = mePartonData()[2]->id() == mePartonData()[3]->id() ?
    0.5 : 1;
  DVector save(ndiags);
  for(DVector::size_type idx = 0; idx < ndiags; ++idx) 
    save[idx] = 0.25*cAvg*ifact*me[idx];
  meInfo(save);
  full_me *= 0.25*cAvg*ifact;
  return full_me;
}

Selector<const ColourLines *>
MEff2vv::colourGeometries(tcDiagPtr diag) const {
  static  vector<ColourLines> cf(9);
  //33b->11
  cf[0] = ColourLines("1 2 -3");
  cf[1] = ColourLines("1 -2");
  //33b->88
  cf[2] = ColourLines("1 4, -4 2 5, -5 -3");
  cf[3] = ColourLines("1 5, -5 2 4, -4 -3");
  cf[4] = ColourLines("1 3 4, -5 -3 -2, -4 5");
  cf[5] = ColourLines("1 3 5, -4 -3 -2, -5 4");
  cf[6] = ColourLines("1 -2, 4 -5, 5 -4");
  //33b->81
  cf[7] = ColourLines("1 4, -4 2 -3");
  cf[8] = ColourLines("1 2 5, -5 -3");

HPDiagram current = getProcessInfo()[abs(diag->id()) - 1];
  vector<ColourLines>::size_type cl(0);
  if(mePartonData()[2]->iColour() == PDT::Colour8 &&
     mePartonData()[3]->iColour() == PDT::Colour8 ) {
    if(current.channelType == HPDiagram::tChannel)
      cl = current.ordered.second ? 2 : 3;
    else {
      if(current.intermediate->iColour() == PDT::Colour0)
	cl = 6;
      else
	cl = 4 + rnd2(0.5, 0.5);
    }
  }
  else if(mePartonData()[2]->iColour() == PDT::Colour0 &&
	  mePartonData()[3]->iColour() == PDT::Colour0 )
    cl = (current.channelType == HPDiagram::tChannel) ? 0 : 1;
  else if(mePartonData()[2]->iColour() == PDT::Colour8 &&
	  mePartonData()[3]->iColour() == PDT::Colour0 )
    cl = (current.channelType == HPDiagram::tChannel) ? 7 : 8;
  else if(mePartonData()[2]->iColour() == PDT::Colour0 &&
	  mePartonData()[3]->iColour() == PDT::Colour8 )
    cl = (current.channelType == HPDiagram::tChannel) ? 8 : 7;
  else {
    return Selector<const ColourLines *>();
}

  Selector<const ColourLines *> sel;
  sel.insert(1., &cf[cl]);
  return sel;
}


void MEff2vv::persistentOutput(PersistentOStream & os) const {
  os << theFerm << theVec << theTen << theSca1 << theSca2;
}

void MEff2vv::persistentInput(PersistentIStream & is, int) {
  is >> theFerm >> theVec >> theTen >> theSca1 >> theSca2;
}

ClassDescription<MEff2vv> MEff2vv::initMEff2vv;
// Definition of the static class description member.

void MEff2vv::Init() {

  static ClassDocumentation<MEff2vv> documentation
    ("This class implements the matrix element calculation of the 2->2 "
     "process, fermion-antifermion -> vector vector");

}

