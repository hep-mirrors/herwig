// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEfv2vf class.
//

#include "MEfv2vf.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include <numeric>

using namespace Herwig;
using Herwig::Helicity::FFVVertexPtr;
using Herwig::Helicity::VVVVertexPtr;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;
using Herwig::Helicity::SpinorWaveFunction;
using Herwig::Helicity::SpinorBarWaveFunction;
using Herwig::Helicity::VectorWaveFunction;

double MEfv2vf::me2() const {
  //wavefunctions
  SpinorVector sp(2);
  VBVector vecIn(2), vecOut(3);
  SpinorBarVector spb(2);
  double fullme(0.);
  if(mePartonData()[0]->id() > 0) {
    for(unsigned int i = 0; i < 2; ++i) {
      sp[i] = SpinorWaveFunction(meMomenta()[0], mePartonData()[0], i, incoming);
      vecIn[i] = VectorWaveFunction(meMomenta()[1], mePartonData()[1], 2*i, 
				    incoming);
      vecOut[2*i] = VectorWaveFunction(meMomenta()[2], mePartonData()[2], 2*i, 
				     outgoing);
      spb[i] = SpinorBarWaveFunction(meMomenta()[3], mePartonData()[3], i, 
				     outgoing);
    }
    if(mePartonData()[2]->mass() > 0.0*MeV)
      vecOut[1] = VectorWaveFunction(meMomenta()[2], mePartonData()[2], 1, 
				     outgoing);
    fv2vfHeME(sp, vecIn, vecOut, spb, fullme);
  }
  else {
    for(unsigned int i = 0; i < 2; ++i) {
      spb[i] = SpinorBarWaveFunction(meMomenta()[0], mePartonData()[0], i, 
				     incoming);
      vecIn[i] = VectorWaveFunction(meMomenta()[1], mePartonData()[1], 2*i, 
				    incoming);
      vecOut[2*i] = VectorWaveFunction(meMomenta()[2], mePartonData()[2], 2*i, 
				     outgoing);
      sp[i] = SpinorWaveFunction(meMomenta()[3], mePartonData()[3], i, outgoing);
    }
    if(mePartonData()[2]->mass() > 0.0*MeV)
      vecOut[1] = VectorWaveFunction(meMomenta()[2], mePartonData()[2], 1, 
				     outgoing);
    fbv2vfbHeME(spb, vecIn, vecOut, sp, fullme);
  }
  return fullme;
}

ProductionMatrixElement
MEfv2vf::fv2vfHeME(SpinorVector & spIn,  VBVector & vecIn, 
		   VBVector & vecOut, SpinorBarVector & spbOut, 
		   double & mesq) const {
  const HPCount ndiags = numberOfDiags();
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors(getColourFactors());
  const Energy2 q2(scale());
  const bool masslessC = mePartonData()[2]->mass() == 0.0*MeV;
  ProductionMatrixElement prodME(PDT::Spin1Half, PDT::Spin1, PDT::Spin1,
				 PDT::Spin1Half);
  vector<Complex> diag(ndiags, Complex(0., 0.));
  vector<double> me(ndiags, 0.);
 
  //intermediate wave functions
  SpinorBarWaveFunction interFB; VectorWaveFunction interV;
  SpinorWaveFunction interF;
  //loop over helicities
  for(unsigned int ifh = 0; ifh < 2; ++ifh) {
    for(unsigned int ivh = 0; ivh < 2; ++ivh) {
      for(unsigned int ovh = 0; ovh < 3; ++ovh) {
	if(masslessC && ovh == 1) ++ovh;
	for(unsigned int ofh = 0; ofh < 2; ++ofh) {
	  //store the flows
	  vector<Complex> cflows(ncf, Complex(0. ,0.));
	  for(HPCount ix = 0; ix < ndiags; ++ix) {
	    HPDiagram diagram = getProcessInfo()[ix];
	    tcPDPtr offshell = diagram.intermediate;
	    if(diagram.channelType == HPDiagram::tChannel) {
	      //t-chan spin-1/2
	      if(offshell->iSpin() == PDT::Spin1Half) {
		interFB = theFerm[ix].second->evaluate(q2, 3, offshell, 
						       spbOut[ofh], vecIn[ivh]);
		diag[ix] = theFerm[ix].first->evaluate(q2, spIn[ifh], interFB, 
						       vecOut[ovh]);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		interV = theVec[ix].second->evaluate(q2, 3, offshell, 
						     vecIn[ivh], vecOut[ovh]);
		diag[ix] = theVec[ix].first->evaluate(q2, spIn[ifh], 
						      spbOut[ofh], interV);
	      }
	      else
		diag[ix] = 0.0;
	    }
	    else if(diagram.channelType == HPDiagram::sChannel) {
	      interFB = theFerm[ix].second->evaluate(q2, 1, offshell, 
						     spbOut[ofh], vecOut[ovh]);
	      diag[ix] = theFerm[ix].first->evaluate(q2, spIn[ifh], interFB, 
						     vecIn[ivh]);
	    }
	    me[ix] += norm(diag[ix]);
	    //Compute flows
 	    for(size_t iy = 0; iy < diagram.colourFlow.size(); ++iy)
	      cflows[diagram.colourFlow[iy].first - 1] += 
		diagram.colourFlow[iy].second * diag[ix];
	    
	  }//end of diag loop
	  //Now add flows to me2 with appropriate colour factors
	  for(size_t ii = 0; ii < ncf; ++ii)
	    for(size_t ij = 0; ij < ncf; ++ij)
	      mesq += cfactors[ii][ij]*(cflows[ii]*conj(cflows[ij])).real();

	  prodME(ifh, 2*ivh, ovh, ofh) = 
	    std::accumulate(cflows.begin(), cflows.end(), Complex(0., 0.));
	}
      }
    }
  }
  const double colfact = mePartonData()[0]->iColour() == PDT::Colour3 ? 
    1./24. : 1;
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = 0.25*colfact*me[ix];
  meInfo(save);
  mesq *= 0.25*colfact;
  return prodME;
}

ProductionMatrixElement
MEfv2vf::fbv2vfbHeME(SpinorBarVector & spbIn,  VBVector & vecIn, 
		     VBVector & vecOut, SpinorVector & spOut, 
		     double & mesq) const {
  const HPCount ndiags = numberOfDiags();
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors(getColourFactors());
  const Energy2 q2(scale());
  const bool masslessC = mePartonData()[2]->mass() == 0.0*MeV;
  ProductionMatrixElement prodME(PDT::Spin1Half, PDT::Spin1, PDT::Spin1,
				 PDT::Spin1Half);
  vector<Complex> diag(ndiags, Complex(0., 0.));
  vector<double> me(ndiags, 0.);

  //intermediate wave functions
  SpinorBarWaveFunction interFB; VectorWaveFunction interV;
  SpinorWaveFunction interF;
  //loop over helicities
  for(unsigned int ifh = 0; ifh < 2; ++ifh) {
    for(unsigned int ivh = 0; ivh < 2; ++ivh) {
      for(unsigned int ovh = 0; ovh < 3; ++ovh) {
	if(masslessC && ovh == 1) ++ovh;
	for(unsigned int ofh = 0; ofh < 2; ++ofh) {
	  //store the flows
	  vector<Complex> cflows(ncf, Complex(0. ,0.));
	  for(HPCount ix = 0; ix < ndiags; ++ix) {
	    HPDiagram diagram = getProcessInfo()[ix];
	    tcPDPtr offshell = diagram.intermediate;
	    if(diagram.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin1Half) {
		interFB = theFerm[ix].first->evaluate(q2, 3, offshell, 
						       spbIn[ifh], vecOut[ovh]);
		diag[ix] = theFerm[ix].second->evaluate(q2, spOut[ofh], interFB, 
						       vecIn[ivh]);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		interV = theVec[ix].first->evaluate(q2, 3, offshell, 
						     spOut[ofh], spbIn[ifh]);
		diag[ix] = theVec[ix].second->evaluate(q2, vecIn[ivh], interV, 
						      vecOut[ovh]);
	      }
	      else diag[ix] = 0.0;
	    }
	    else if(diagram.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin1Half) {
		interFB = theFerm[ix].first->evaluate(q2, 1, offshell, spbIn[ifh], 
						      vecIn[ivh]);
		diag[ix] = theFerm[ix].second->evaluate(q2, spOut[ofh], interFB, 
						       vecOut[ovh]);
	      }
	    }
	    me[ix] += norm(diag[ix]);
	    //Compute flows
 	    for(size_t iy = 0; iy < diagram.colourFlow.size(); ++iy)
	      cflows[diagram.colourFlow[iy].first - 1] += 
		diagram.colourFlow[iy].second * diag[ix];

	  }//end of diag loop
	  //Now add flows to me2 with appropriate colour factors
	  for(size_t ii = 0; ii < ncf; ++ii)
	    for(size_t ij = 0; ij < ncf; ++ij)
	      mesq += cfactors[ii][ij]*(cflows[ii]*conj(cflows[ij])).real();
	  
	  prodME(ifh, ivh, ovh, ofh) = 
	    std::accumulate(cflows.begin(), cflows.end(), Complex(0., 0.));
	}
      }
    }
  }
  const double colfact = mePartonData()[0]->iColour() == PDT::Colour3bar ? 
    1./24. : 1;
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = 0.25*colfact*me[ix];
  meInfo(save);
  mesq *= 0.25*colfact;
  return prodME;
}

Selector<const ColourLines *>
MEfv2vf::colourGeometries(tcDiagPtr diag) const {
  static vector<ColourLines> cf(10);
  //3 8->8 3
  cf[0] = ColourLines("1 4, -4 2 -3, 3 5");
  cf[1] = ColourLines("1 2 -3, -4 -2 5, 3 4");
  cf[2] = ColourLines("1 -2, 2 3 4, -4 5");
  //3b 8 -> 8 3b
  cf[3] = ColourLines("-4 -1, 3 2 4, -5 -3");
  cf[4] = ColourLines("3 2 -1, -5 -2 4, -4 -3");
  cf[5] = ColourLines("2 -1, -4 -3 -2, -5 4");
  //3 8 -> 0 3
  cf[6] = ColourLines("1 2 -3, 3 5");
  cf[7] = ColourLines("1 -2, 2 3 5");
  //3b 8 -> 0 3b
  cf[8] = ColourLines("3 2 -1, -3 -5");
  cf[9] = ColourLines("2 -1, -5 -3 -2");
  
  HPDiagram current = getProcessInfo()[abs(diag->id()) - 1];
  PDT::Colour offcolour = current.intermediate->iColour();
  bool octetC = (mePartonData()[2]->iColour() == PDT::Colour8);
  bool cc(mePartonData()[0]->id() < 0);
  
  Selector<const ColourLines *> select;
  int icf(-1);
  if(current.channelType == HPDiagram::sChannel)
    if(octetC) icf = cc ? 5 : 2;
    else icf = cc ? 9 : 7;
  else if(current.channelType == HPDiagram::tChannel) {
    if(offcolour == PDT::Colour3 || offcolour == PDT::Colour3bar) {
      if(octetC) icf = cc ? 3 : 0; 
      else  icf = cc ? 8 : 6;
    }
    else if(offcolour == PDT::Colour8) 
      icf = cc ? 4 : 1;
    else
      throw MEException() << "MEfv2vf::colourGeometries - There is an incorrect "
			  << "coloured object in a t-channel diagram. "
			  << "No colour lines set."
			  << Exception::warning;
  }
  else
    throw MEException() << "MEfv2vf::colourGeometries - Incorrect diagram type "
			<< "encountered. No colour lines set." 
			<< Exception::warning;
  if(icf >= 0) select.insert(1., &cf[icf]);
  return select;
}


void MEfv2vf::persistentOutput(PersistentOStream & os) const {
  os << theFerm << theVec;
}

void MEfv2vf::persistentInput(PersistentIStream & is, int) {
  is >> theFerm >> theVec;
}

ClassDescription<MEfv2vf> MEfv2vf::initMEfv2vf;
// Definition of the static class description member.

void MEfv2vf::Init() {

  static ClassDocumentation<MEfv2vf> documentation
    ("This is the implementation of the matrix element for a fermion-vector boson"
     "to a vector-fermion.");

}

