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
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "HardVertex.h"
#include<numeric>

using namespace Herwig;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::FFVVertexPtr;
using ThePEG::Helicity::FFTVertexPtr;
using ThePEG::Helicity::FFSVertexPtr;
using ThePEG::Helicity::VVSVertexPtr;
using ThePEG::Helicity::GeneralSVVVertexPtr;
using ThePEG::Helicity::VVTVertex;
using ThePEG::Helicity::VVVVertexPtr;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

typedef Ptr<VVTVertex>::pointer VVTVertexPtr;

double MEff2vv::me2() const {
  SpinorVector sp(2);
  SpinorBarVector sbar(2);
  // vector 
  bool mc  = !(mePartonData()[2]->mass() > 0.*MeV);
  bool md  = !(mePartonData()[3]->mass() > 0.*MeV);
  VBVector v1(3), v2(3);  
  for( unsigned int i = 0; i < 2; ++i ) {
    sp[i] = SpinorWaveFunction(meMomenta()[0], mePartonData()[0], i, incoming);
    sbar[i] = SpinorBarWaveFunction(meMomenta()[1], mePartonData()[1], i, 
				    incoming);
    v1[2*i] = VectorWaveFunction(meMomenta()[2], mePartonData()[2],2*i , 
				 outgoing);
    v2[2*i] = VectorWaveFunction(meMomenta()[3], mePartonData()[3], 2*i, 
				 outgoing);
  }
  if( !mc ) v1[1] = VectorWaveFunction(meMomenta()[2], mePartonData()[2], 1, 
				       outgoing);
  if( !md ) v2[1] = VectorWaveFunction(meMomenta()[3], mePartonData()[3], 1, 
				       outgoing);
  double full_me(0.);
  ff2vvME(sp, sbar, v1, mc, v2, md, full_me);

#ifndef NDEBUG
  if(  debugME() ) debug(full_me);
#endif
  return full_me;
}

ProductionMatrixElement 
MEff2vv::ff2vvME(const SpinorVector & sp, const SpinorBarVector sbar, 
		 const VBVector & v1, bool m1, const VBVector & v2, bool m2,
		 double & me2) const {
    //Define vectors to store diagrams and square elements
  const HPCount ndiags(numberOfDiags());
  const size_t ncf(numberOfFlows());
  vector<Complex> diag(ndiags, Complex(0.));
  vector<double> me(ndiags, 0.);
  const Energy2 q2 = scale();
  const vector<vector<double> > cfactors = getColourFactors();
  //Intermediate wavefunctions
  ScalarWaveFunction interS; SpinorWaveFunction interF;
  VectorWaveFunction interV; TensorWaveFunction interT;
  ProductionMatrixElement pme(PDT::Spin1Half, PDT::Spin1Half,
			      PDT::Spin1, PDT::Spin1);
  //Loop over helicities and calculate diagrams
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      for(unsigned int vh1 = 0; vh1 < 3; ++vh1) {
 	if( vh1 == 1 && m1 ) ++vh1;
	for(unsigned int vh2 = 0; vh2 < 3; ++vh2) {
	  if( vh2 == 1 && m2 ) ++vh2;
	  //loop and calculate diagrams
	  vector<Complex> flows = vector<Complex>(ncf, Complex(0.));
	  for(HPCount ix = 0; ix < ndiags; ++ix) {
	    HPDiagram current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel) {
	      if(current.intermediate->iSpin() == PDT::Spin1Half) {
		if(current.ordered.second) {
		  interF = theFerm[ix].first->evaluate(q2, 3, offshell, sp[if1], 
						       v1[vh1]);
		  diag[ix] = theFerm[ix].second->evaluate(q2, interF, sbar[if2],
							  v2[vh2]);
		}
		else {
		  interF = theFerm[ix].first->evaluate(q2, 3 , offshell, sp[if1], 
						       v2[vh2]);
		  diag[ix] = theFerm[ix].second->evaluate(q2, interF, sbar[if2],  
							  v1[vh1]);
		}	      
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(current.intermediate->iSpin() == PDT::Spin0) {
		if( m1 && m2 ) {
		  interS = theSca1[ix].first->evaluate(q2, 1, offshell, 
						       sp[if1], sbar[if2]);
		  diag[ix] = theSca1[ix].second->evaluate(q2, interS, 
							  v1[vh1], v2[vh2]);
 		}
		else {
		  interS = theSca2[ix].first->evaluate(q2, 1, offshell, 
						       sp[if1], sbar[if2]);
		  diag[ix] = theSca2[ix].second->evaluate(q2, v1[if1], v2[if2], 
							  interS);
		}
	      }
	      else if(current.intermediate->iSpin() == PDT::Spin1) {
		interV = theVec[ix].first->evaluate(q2, 5, offshell, sp[if1], 
						    sbar[if2]);
 		diag[ix] = theVec[ix].second->evaluate(q2, v1[vh1], 
						       v2[vh2], interV);
	      }
	      else if(current.intermediate->iSpin() == PDT::Spin2) {
		interT = theTen[ix].first->evaluate(q2, 1, offshell, sp[if1],
						    sbar[if2]);
		diag[ix] = theTen[ix].second->evaluate(q2, v1[vh1], v2[vh2], 
						       interT);
	      }
	    }
	    //compute individual diagram squared and save for diagram sel
	    me[ix] += norm(diag[ix]);
	    //find out which colourflow this diag belongs to and add to it
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy)
	      flows[current.colourFlow[iy].first - 1] +=
		current.colourFlow[iy].second*diag[ix];

	  }//end-of-diagram loop

	  //set appropriate element in the ProductionMatrixElement
	  pme(if1, if2, vh1, vh2) = 
	    std::accumulate(flows.begin(), flows.end(), Complex(0.0, 0.0));
	  
	  //Add results, with appropriate colour factors to me2
	  for(size_t ii = 0; ii < ncf; ++ii)
	    for(size_t ij = 0; ij < ncf; ++ij)
	      me2 += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();
	  
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
  me2 *= 0.25*cAvg*ifact;
  return pme;
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

void MEff2vv::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  if( (mePartonData()[0]->id() != 1 && mePartonData()[0]->id() != 2) ||
      (mePartonData()[1]->id() != -1 && mePartonData()[1]->id() != -2) ||
      mePartonData()[2]->id() != 5100021 || 
      mePartonData()[3]->id() != 5100021 ) return;
  
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  Energy2 s(sHat());
  Energy2 mf2 = meMomenta()[2].m2();
  Energy2 t3(tHat() - mf2), u4(uHat() - mf2);
  double analytic = gs4*( mf2*( (57.*s/t3/u4)  - (4.*s*s*s/t3/t3/u4/u4) 
				- (108./s) )  
			  + (20.*s*s/t3/u4) - 93. + (108.*t3*u4/s/s) )/27.;
  
  double diff = abs( analytic - me2 );
  if( diff  > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << "   ratio: " 
      << analytic/me2 << '\n';
  }


}
