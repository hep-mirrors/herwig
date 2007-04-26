// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2ss class.
//

#include "MEff2ss.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

using namespace Herwig;
using Herwig::Helicity::SpinorWaveFunction;
using Herwig::Helicity::SpinorBarWaveFunction;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::TensorWaveFunction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

double MEff2ss::me2() const {
  //first setup  wavefunctions for external particles
  SpinorWaveFunction sp(meMomenta()[0], mePartonData()[0], incoming);
  SpinorBarWaveFunction spbar(meMomenta()[1], mePartonData()[1], incoming);
  ScalarWaveFunction sca1(meMomenta()[2], mePartonData()[2],
			  Complex(1.), outgoing);
  ScalarWaveFunction sca2(meMomenta()[3], mePartonData()[3],
			  Complex(1.), outgoing);
  //Define factors
  const Energy2 m2(scale());
  const vector<vector<double> > cfactors = getColourFactors();
  const HPCount ndiags = numberOfDiags();
  const size_t ncf = numberOfFlows();
  vector<double> me(ndiags, 0.);
  vector<Complex> diag(ndiags, Complex(0.)), flows(ncf, Complex(0.));
  double full_me(0.);
  ScalarWaveFunction interS; VectorWaveFunction interV; 
  SpinorBarWaveFunction interFB; TensorWaveFunction interT;
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) {
    sp.reset(ihel1);
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      spbar.reset(ihel2);
      flows = vector<Complex>(ncf, Complex(0.));
      for(HPCount ix = 0; ix < ndiags; ++ix) {
	HPDiagram current = getProcessInfo()[ix];
	tcPDPtr internal(current.intermediate);	
	if(current.channelType == HPDiagram::tChannel &&
	   internal->iSpin() == PDT::Spin1Half) {
	  if(current.ordered.second) {
	    interFB = theFerm[ix].second->evaluate(m2, 3, internal, spbar, sca2);
	    diag[ix] = theFerm[ix].first->evaluate(m2, sp, interFB, sca1);
	  }
	  else {
	    interFB = theFerm[ix].second->evaluate(m2, 3, internal, spbar, sca1);
	    diag[ix] = theFerm[ix].first->evaluate(m2, sp, interFB, sca2);
	  }
	}
	else if(current.channelType == HPDiagram::sChannel) {
	  if(internal->iSpin() == PDT::Spin1) {
	    interV = theVec[ix].first->evaluate(m2, 1, internal, sp, spbar);
	    diag[ix] = theVec[ix].second->evaluate(m2, interV, sca2, sca1);
	  }
	  else if(internal->iSpin() == PDT::Spin2) {
	    interT = theTen[ix].first->evaluate(m2, 1, internal, sp, spbar);
	    diag[ix] = theTen[ix].second ->evaluate(m2, sca2, sca1, interT);
	  }
	  else 
	    diag[ix] = 0.;
	}
	me[ix] += norm(diag[ix]);
	//colourflows
	for(unsigned int iy = 0; iy < current.colourFlow.size(); ++iy)
	  flows[current.colourFlow[iy].first - 1] += 
	    current.colourFlow[iy].second * diag[ix];

      }//end of diag loop
      for(unsigned int ii = 0; ii < ncf; ++ii) 
	for(unsigned int ij = 0; ij < ncf; ++ij)
	  full_me += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();
    }
  }
  const double identFact = mePartonData()[2]->id() == mePartonData()[3]->id() 
    ? 0.5 : 1;
  const double colourAvg = mePartonData()[0]->iColour() == PDT::Colour3 
    ? 1./9. : 1.;
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = 0.25*identFact*colourAvg*me[ix];
  meInfo(save);
  full_me *= 0.25*identFact*colourAvg;
  return full_me;
}

Selector<const ColourLines *>
MEff2ss::colourGeometries(tcDiagPtr diag) const {
  static vector<ColourLines> cf(16);
  //33b->33b
  cf[0] = ColourLines("1 2 -3, 4 -2 -5");
  cf[1] = ColourLines("1 3 4, -2 -3 -5");
  cf[2] = ColourLines("1 4, -3 -5");	 
  cf[3] = ColourLines("1 -2, 4 -5");
  //33->33
  cf[4] = ColourLines("1 2 5, 3 -2 4");  
  cf[5] = ColourLines("1 2 4, 3 -2 5");
  cf[6] = ColourLines("1 4, 3 5");	
  cf[7] = ColourLines("1 5, 3 4");
  //3b3b->3b3b
  cf[8] = ColourLines("-1 -2 -5, -3 2 -4");
  cf[9] = ColourLines("-1 -2 -4, -3 2 -5");
  cf[10] = ColourLines("-1 -4, -3 -5");
  cf[11] = ColourLines("-1 -5, -3 -4");  
  //33b->11
  cf[12] = ColourLines("1 2 -3");
  cf[13] = ColourLines("1 -2");
  //11->11
  cf[14] = ColourLines("");
  //11->33b
  cf[15] = ColourLines("4 -5");
  HPDiagram current = getProcessInfo()[abs(diag->id()) - 1];
  vector<ColourLines>::size_type cl(0);
  PDT::Colour inac(mePartonData()[0]->iColour());
  PDT::Colour inbc(mePartonData()[1]->iColour());
  PDT::Colour outac(mePartonData()[2]->iColour());
  if(inac == PDT::Colour0 && inbc == PDT::Colour0) {
    cl = outac == PDT::Colour0 ?  14 : 15;
  }
  else if(inac == PDT::Colour3 && inbc == PDT::Colour3) {
    if(current.intermediate->iColour() == PDT::Colour8)
      cl = current.ordered.second ? 4 : 5;
    else
      cl = current.ordered.second ? 6 : 7;
  }
  else if(inac == PDT::Colour3bar && inbc == PDT::Colour3bar) {
    if(current.intermediate->iColour() == PDT::Colour8)
      cl = current.ordered.second ? 8 : 9;
    else
      cl = current.ordered.second ? 10 : 11 ;
  }
  else {
    if(outac == PDT::Colour3) {
      if(current.intermediate->iColour() == PDT::Colour8)
	cl = current.channelType == HPDiagram::tChannel ? 0 : 1;
      else
	cl = current.channelType == HPDiagram::tChannel ? 2 : 3;
    }
    else
      cl = current.channelType == HPDiagram::tChannel ? 12 : 13;
  }
  Selector<const ColourLines *> sel;
  sel.insert(1., &cf[cl]);
  return sel;
}

void MEff2ss::persistentOutput(PersistentOStream & os) const {
  os << theFerm << theVec << theTen;
}

void MEff2ss::persistentInput(PersistentIStream & is, int) {
  is >> theFerm >> theVec >> theTen;
}

ClassDescription<MEff2ss> MEff2ss::initMEff2ss;
// Definition of the static class description member.

void MEff2ss::Init() {

  static ClassDocumentation<MEff2ss> documentation
    ("MEff2ss implements the ME calculation of the fermion-antifermion "
     "to scalar-scalar hard process.");

}

