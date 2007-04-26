// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2ff class.
//

#include "MEff2ff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/Helicity/Correlations/HardVertex.h"
#include <numeric>

using namespace Herwig;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::TensorWaveFunction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;
using Herwig::Helicity::HardVertex;
using Herwig::Helicity::HardVertexPtr;
using ThePEG::Helicity::SpinfoPtr;

double MEff2ff::me2() const {
  // Three different cases to handle, Psi,Psibar->Psi,Psibar, Psi,Psi->PsiPsi &
  // Psi,Psibar->lambda,lambda
  tcPDPtr outa(mePartonData()[2]), outb(mePartonData()[3]);
  bool majorana(false);
  if( (!outa->CC() && !outb->CC() ) || 
      ((abs(outa->id()) > 1000000 && abs(outa->id()) < 2000000) &&
       (abs(outb->id()) > 1000000 && abs(outb->id()) < 2000000)) )
    majorana = true;
  
  double full_me(0.);
  vector<SpinorWaveFunction> spA(2), spB(2);
  vector<SpinorBarWaveFunction> spbA(2), spbB(2);
  if(mePartonData()[1]->id() < 0) {
    for(unsigned int ih = 0; ih < 2; ++ih) {
      spA[ih] = SpinorWaveFunction(meMomenta()[0], mePartonData()[0], ih, 
				   incoming);
      spbA[ih] = SpinorBarWaveFunction(meMomenta()[1], mePartonData()[1], ih, 
				       incoming);
      spB[ih] = SpinorWaveFunction(meMomenta()[3], outb, ih, outgoing);
      spbB[ih] = SpinorBarWaveFunction(meMomenta()[2], outa, ih, outgoing);
    }
    if(majorana) {
      vector<SpinorWaveFunction> spC(2);
      vector<SpinorBarWaveFunction> spbC(2);
      for(unsigned int ih = 0; ih < 2; ++ih) {
	spC[ih] = SpinorWaveFunction(meMomenta()[2], outa, ih, outgoing);
	spbC[ih] = SpinorBarWaveFunction(meMomenta()[3], outb, ih, outgoing);
      }
      ffb2mfmfHeME(spA, spbA, spbB, spB, spC, spbC, full_me);
      SpinorWaveFunction spOut2(meMomenta()[2], outa, outgoing);
      SpinorBarWaveFunction spbarOut2(meMomenta()[3], outb, outgoing);
     }
    else 
      ffb2ffbHeME(spA, spbA, spbB, spB, full_me); 
  }
  else {
    SpinorVector spA(2), spB(2);
    SpinorBarVector spbA(2), spbB(2);
    for(unsigned int ih = 0; ih < 2; ++ih) {
      spA[ih] = SpinorWaveFunction(meMomenta()[0], mePartonData()[0], ih,
				   incoming);
      spB[ih] = SpinorWaveFunction(meMomenta()[1], mePartonData()[1], ih,
				   incoming);
      spbA[ih] = SpinorBarWaveFunction(meMomenta()[2], outa, ih,
				       outgoing);
      spbB[ih] = SpinorBarWaveFunction(meMomenta()[3], outb, ih,
				       outgoing);
    }
    ff2ffHeME(spA, spB, spbA, spbB, full_me);
  }
  return full_me;
}

ProductionMatrixElement 
MEff2ff::ffb2ffbHeME(SpinorVector & fin, SpinorBarVector & fbin,
		     SpinorBarVector & fbout, SpinorVector & fout,
		     double & me2) const {
  const HPCount ndiags = numberOfDiags();
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors(getColourFactors());
  const Energy2 m2(scale());
  vector<Complex> diag(ndiags, Complex(0.)), flows(ncf, Complex(0.));
  vector<double> me(ndiags, 0.);
  ScalarWaveFunction interS; VectorWaveFunction interV;
  TensorWaveFunction interT;
  ProductionMatrixElement prodME(PDT::Spin1Half, PDT::Spin1Half, PDT::Spin1Half,
				 PDT::Spin1Half);
  for(unsigned int ifhel1 = 0; ifhel1 < 2; ++ifhel1) {
    for(unsigned int ifhel2 = 0; ifhel2 < 2; ++ifhel2) {
      for(unsigned int ofhel1 = 0; ofhel1 < 2; ++ofhel1) {
	for(unsigned int ofhel2 = 0; ofhel2 < 2; ++ofhel2) {
	  flows = vector<Complex>(ncf, Complex(0.));
	  for(HPCount ix = 0; ix < ndiags; ++ix) {
	    HPDiagram current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		interS = theScaV[ix].second->evaluate(m2, 3, offshell, 
						      fout[ofhel2], fbin[ifhel2]);
		diag[ix] = theScaV[ix].first->evaluate(m2, fin[ifhel1],
						       fbout[ofhel1], interS);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		interV = theVecV[ix].second->evaluate(m2, 3, offshell, 
						      fout[ofhel2], fbin[ifhel2]);
		diag[ix] = -theVecV[ix].first->evaluate(m2, fin[ifhel1],
						       fbout[ofhel1], interV);
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		interT = theTenV[ix].second->evaluate(m2, 3, offshell, 
						      fout[ofhel2], fbin[ifhel2]);
		diag[ix] = theTenV[ix].first->evaluate(m2, fin[ifhel1],
						       fbout[ofhel1], interT);
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		interS = theScaV[ix].second->evaluate(m2, 1, offshell,
						      fout[ofhel2],fbout[ofhel1]);
		diag[ix] = theScaV[ix].first->evaluate(m2, fin[ifhel1], 
						       fbin[ifhel2], interS);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		interV = theVecV[ix].second->evaluate(m2, 1, offshell,
						      fout[ofhel2],fbout[ofhel1]);
		diag[ix] = theVecV[ix].first->evaluate(m2, fin[ifhel1], 
						       fbin[ifhel2], interV);
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		interT = theTenV[ix].second->evaluate(m2, 1, offshell,
						      fout[ofhel2],fbout[ofhel1]);
		diag[ix] = theTenV[ix].first->evaluate(m2, fin[ifhel1], 
						       fbin[ifhel2], interT);
	      }
	    }
	    else {
	      throw MEException() << "Incorrect diagram type in matrix element "
				  << fullName()
				  << Exception::warning;
	      diag[ix] = 0.;
	    }
	    me[ix] += norm(diag[ix]);
	    //Compute flows
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy)
	      flows[current.colourFlow[iy].first - 1] += 
		current.colourFlow[iy].second * diag[ix];
	    
	  }//diagram loop

	  //Now add flows to me2 with appropriate colour factors
	  for(size_t ii = 0; ii < ncf; ++ii)
	    for(size_t ij = 0; ij < ncf; ++ij)
	      me2 += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();

	  prodME(ifhel1, ifhel2, ofhel1, ofhel2) = 
	    std::accumulate(flows.begin(), flows.end(), Complex(0., 0.));
	}//end of first helicity loop
      }
    }
  }
  const double identfact = mePartonData()[2]->id() == mePartonData()[3]->id() 
    ? 0.5 : 1.;
  const double colfact = mePartonData()[0]->iColour() == PDT::Colour3 ? 
    1./9. : 1;
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = 0.25*identfact*colfact*me[ix];
  meInfo(save);
  me2 *= 0.25*identfact*colfact;
  
  return prodME;
}

ProductionMatrixElement 
MEff2ff:: ff2ffHeME(SpinorVector & fin, SpinorVector & fin2,
		    SpinorBarVector & fbout, SpinorBarVector & fbout2,
		    double & me2) const {
  const HPCount ndiags = getProcessInfo().size();
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors(getColourFactors());
  const Energy2 q2(scale());
  vector<Complex> diag(ndiags, Complex(0.)), flows(ncf, Complex(0.));
  vector<double> me(ndiags, 0.);
  ScalarWaveFunction interS; VectorWaveFunction interV;
  TensorWaveFunction interT;
  ProductionMatrixElement prodME(PDT::Spin1Half, PDT::Spin1Half, PDT::Spin1Half,
				 PDT::Spin1Half);
  for(unsigned int ifhel1 = 0; ifhel1 < 2; ++ifhel1) {
    for(unsigned int ifhel2 = 0; ifhel2 < 2; ++ifhel2) {
      for(unsigned int ofhel1 = 0; ofhel1 < 2; ++ofhel1) {
	for(unsigned int ofhel2 = 0; ofhel2 < 2; ++ofhel2) {
	  flows = vector<Complex>(ncf, Complex(0.));
	  for(HPCount ix = 0; ix < ndiags; ++ix) {
	    HPDiagram current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		if(current.ordered.second) {
		  interS = theScaV[ix].second->evaluate(q2, 3, offshell, 
							fin2[ifhel2], 
							fbout2[ofhel2]);
		  diag[ix] = theScaV[ix].first->evaluate(q2, fin[ifhel1], 
							 fbout[ofhel1], interS);
		}
		else {
		  interS = theScaV[ix].second->evaluate(q2, 3, offshell, 
							fin2[ifhel2], 
							fbout[ofhel1]);
		  diag[ix] = theScaV[ix].first->evaluate(q2, fin[ifhel1], 
							 fbout2[ofhel2], interS);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  interV = theVecV[ix].second->evaluate(q2, 3, offshell, 
							fin2[ifhel2],
							fbout2[ofhel2]);
		  diag[ix] = theVecV[ix].first->evaluate(q2, fin[ifhel1], 
							 fbout[ofhel1], interV);
		}
		else {
		  interV = theVecV[ix].second->evaluate(q2, 3, offshell, 
							fin2[ifhel2], 
							fbout[ofhel1]);
		  diag[ix] = -theVecV[ix].first->evaluate(q2, fin[ifhel1], 
							  fbout2[ofhel2], interV);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		if(current.ordered.second) {
		  interT = theTenV[ix].second->evaluate(q2, 3, offshell, 
							fin2[ifhel2], 
							fbout2[ofhel2]);
		  diag[ix] = theTenV[ix].first->evaluate(q2, fin[ifhel1], 
							 fbout[ofhel1], interT);
		}
		else {
		  interT = theTenV[ix].second->evaluate(q2, 3, offshell, 
							fin2[ifhel2], 
							fbout[ofhel1]);
		  diag[ix] = theTenV[ix].first->evaluate(q2, fin[ifhel1], 
							 fbout2[ofhel2],
							 interT);
		}
	      }
	    }
	    else {
	      throw MEException() << "Incorrect diagram type in matrix element "
				  << fullName()
				  << Exception::warning;
	      diag[ix] = 0.;
	    }
	    me[ix] += norm(diag[ix]);
	    //Compute flows
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy)
	      flows[current.colourFlow[iy].first - 1] += 
		current.colourFlow[iy].second * diag[ix];
 	    
	  }//diagram loop
	  //Now add flows to me2 with appropriate colour factors
	  for(size_t ii = 0; ii < ncf; ++ii)
	    for(size_t ij = 0; ij < ncf; ++ij)
	      me2 += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();

	  prodME(ifhel1, ifhel2, ofhel1, ofhel2) = 
	    std::accumulate(flows.begin(), flows.end(), Complex(0., 0.));

	}
      }
    }
  }
  const double identfact = mePartonData()[2]->id() == mePartonData()[3]->id() 
    ? 0.5 : 1.;
  const double colfact = mePartonData()[0]->iColour() == PDT::Colour3 ? 1./9. : 1;
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = 0.25*identfact*colfact*me[ix];
  meInfo(save);
  me2 = 0.25*identfact*colfact*me2;
  return prodME;
}

ProductionMatrixElement 
MEff2ff::ffb2mfmfHeME(SpinorVector & fin, SpinorBarVector & fbin, 
		      SpinorBarVector & fbout, SpinorVector & fout,
		      SpinorVector & fout2, SpinorBarVector & fbout2,
		      double & me2) const {
  //useful constants
  const HPCount ndiags = numberOfDiags();
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors(getColourFactors());
  const Energy2 m2(scale());
  //store results
  vector<Complex> diag(ndiags, Complex(0.)), flows(ncf, Complex(0.));
  vector<double> me(ndiags, 0.);
  //intermediate wavefunctions
  ScalarWaveFunction interS; VectorWaveFunction interV;
  TensorWaveFunction interT;
  //ProductionMatrixElement object
  ProductionMatrixElement prodME(PDT::Spin1Half, PDT::Spin1Half,
				 PDT::Spin1Half, PDT::Spin1Half);
  for(unsigned int ifhel1 = 0; ifhel1 < 2; ++ifhel1) {
    for(unsigned int ifhel2 = 0; ifhel2 < 2; ++ifhel2) {
      for(unsigned int ofhel1 = 0; ofhel1 < 2; ++ofhel1) {
	for(unsigned int ofhel2 = 0; ofhel2 < 2; ++ofhel2) {
	  flows = vector<Complex>(ncf, Complex(0.));
	  for(size_t ix = 0; ix < ndiags; ++ix) {
	    HPDiagram current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		if(current.ordered.second) {
		  interS = theScaV[ix].second->evaluate(m2, 3, offshell, 
							fout[ofhel2], 
							fbin[ifhel2]);
		  diag[ix] = theScaV[ix].first->evaluate(m2, fin[ifhel1], 
							 fbout[ofhel1], 
							 interS);
		}
		else {
		  interS = theScaV[ix].second->evaluate(m2, 3, offshell, 
							fout2[ofhel1], 
							fbin[ifhel2]);
		  diag[ix] = -theScaV[ix].first->evaluate(m2, fin[ifhel1],
							  fbout2[ofhel2],
							  interS);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  interV = theVecV[ix].second->evaluate(m2, 3, offshell, 
							fout[ofhel2], 
							fbin[ifhel2]);
		  diag[ix] = theVecV[ix].first->evaluate(m2, fin[ifhel1], 
							 fbout[ofhel1], 
							 interV);
		}
		else {
		  interV = theVecV[ix].second->evaluate(m2, 3, offshell, 
							fout2[ofhel1], 
							fbin[ifhel2]);
		  diag[ix] = theVecV[ix].first->evaluate(m2, fin[ifhel1],
							 fbout2[ofhel2],
							 interV);
 		}
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		if(current.ordered.second) {
		  interT = theTenV[ix].second->evaluate(m2, 3, offshell, 
							fout[ofhel2], 
							fbin[ifhel2]);
		  diag[ix] = theTenV[ix].first->evaluate(m2, fin[ifhel1], 
							 fbout[ofhel1], 
							 interT);
		}
		else {
		  interT = theTenV[ix].second->evaluate(m2, 3, offshell, 
							fout2[ofhel1], 
							fbin[ifhel2]);
		  diag[ix] = theTenV[ix].first->evaluate(m2, fin[ifhel1],
							 fbout2[ofhel2],
							 interT);
 		}
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		interS = theScaV[ix].second->evaluate(m2, 1, offshell,
						      fout[ofhel2],fbout[ofhel1]);
		diag[ix] = theScaV[ix].first->evaluate(m2, fin[ifhel1], 
						       fbin[ifhel2], interS);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		interV = theVecV[ix].second->evaluate(m2, 1, offshell,
						      fout[ofhel2],fbout[ofhel1]);
		diag[ix] = theVecV[ix].first->evaluate(m2, fin[ifhel1], 
						       fbin[ifhel2], interV);
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		interT = theTenV[ix].second->evaluate(m2, 1, offshell,
						      fout[ofhel2],fbout[ofhel1]);
		diag[ix] = theTenV[ix].first->evaluate(m2, fin[ifhel1], 
						       fbin[ifhel2], interT);
	      }
	    }
	    else {
	      throw MEException() << "Incorrect diagram type in matrix element "
				  << fullName()
				  << Exception::warning;
	      diag[ix] = 0.;
	    }
	    me[ix] += norm(diag[ix]);
	    // Compute flows
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy)
	      flows[current.colourFlow[iy].first - 1] += 
		current.colourFlow[iy].second * diag[ix];
	    
	  }//diagram loop	  
	  //Now add flows to me2 with appropriate colour factors
	  for(size_t ii = 0; ii < ncf; ++ii)
	    for(size_t ij = 0; ij < ncf; ++ij)
	      me2 += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();
	  
	  prodME(ifhel1, ifhel2, ofhel1, ofhel2) = 
	    std::accumulate(flows.begin(), flows.end(), Complex(0., 0.));
	    
	}// ofhel2 loop
      }
    }
  }	  
  return prodME;
}

Selector<const ColourLines *>
MEff2ff::colourGeometries(tcDiagPtr diag) const {
  static vector<ColourLines> cf(24);
  //33b->11
  cf[0] = ColourLines("1 2 -3");
  cf[1] = ColourLines("1 -2");
  //33b->18
  cf[2] = ColourLines("1 2 5, -3 -5");
  cf[3] = ColourLines("1 4, -4 2 -3");
  //33b->33bar
  cf[4] = ColourLines("1 2 -3, 4 -2 -5");
  cf[5] = ColourLines("1 3 4, -2 -3 -5");
  cf[6] = ColourLines("1 4, -3 -5");
  cf[7] = ColourLines("1 -2, 4 -5");     
  //33b->88
  cf[8] = ColourLines("1 4, -4 2 5, -5 -3");
  cf[9] = ColourLines("1 5, -5 2 4, -4 -3");
  cf[10] = ColourLines("1 3 4, -5 -3 -2, -4 5");
  cf[11] = ColourLines("1 3 5, -4 -3 -2, -5 4");
  //33->33
  cf[12] = ColourLines("1 2 5, 3 -2 4");
  cf[13] = ColourLines("1 2 4, 3 -2 5");
  cf[14] = ColourLines("1 4, 3 5");
  cf[15] = ColourLines("1 5, 3 4");
  //3b->3b
  cf[16] = ColourLines("-1 -2 -5, -3 2 -4");
  cf[17] = ColourLines("-1 -2 -4, -3 2 -5");
  cf[18] = ColourLines("-1 -4, -3 -5");
  cf[19] = ColourLines("-1 -5, -3 -4");  
  //33b->81
  cf[20]= ColourLines("1 4, -4 2 -3");
  cf[21]= ColourLines("-3 -4, 1 2 4");
  //11->33bar
  cf[22] = ColourLines("4 -5");
  //11->11
  cf[23] = ColourLines("");
  
  HPDiagram current = getProcessInfo()[abs(diag->id()) - 1];
  //select appropriate set of diagrams
  PDT::Colour inac(mePartonData()[0]->iColour());
  PDT::Colour inbc(mePartonData()[1]->iColour());
  PDT::Colour outac(mePartonData()[2]->iColour());
  PDT::Colour outbc(mePartonData()[3]->iColour());
  vector<ColourLines>::size_type cl(0);
  if(inac == PDT::Colour3 && inbc == PDT::Colour3) {
    if(current.intermediate->iColour() == PDT::Colour8)
      cl = current.ordered.second ? 12 : 13;
    else
      cl = current.ordered.second ? 14 : 15;
  }
  else if(inac == PDT::Colour3bar && inbc == PDT::Colour3bar) {
    if(current.intermediate->iColour() == PDT::Colour8)
      cl = current.ordered.second ? 16 : 17;
    else
      cl = current.ordered.second ? 18 : 19;
  }
  else if(inac == PDT::Colour0 && inbc == PDT::Colour0)
    cl = (outac == PDT::Colour0) ? 23 : 22;
  else {
    if(outac == PDT::Colour0 || outbc == PDT::Colour0 ) {
      if(current.channelType == HPDiagram::tChannel) {
	if(outac == outbc)
	  cl = current.ordered.second ? 0 : 1;
	else if(outac == PDT::Colour0 && outbc == PDT::Colour8)
	  cl = current.ordered.second ? 2 : 3;
	else
	  cl = current.ordered.second ? 20 : 21;
      }
    }
    else if(outbc == PDT::Colour3bar) {
      if(current.channelType == HPDiagram::tChannel)
	cl = (current.intermediate->iColour() == PDT::Colour8) ? 4 : 6;
      else
	cl = (current.intermediate->iColour() == PDT::Colour8) ? 5 : 7;
    }
    else {
      if(current.channelType == HPDiagram::tChannel)
	cl = current.ordered.second ? 8 : 9;
      else
	cl = 10 + rnd2(0.5, 0.5);
    }
  }
  Selector<const ColourLines *> sel;
  sel.insert(1., &cf[cl]);
  return sel;
}

void MEff2ff::constructVertex(tSubProPtr subp) {
  // Hard proces external particles
  ParticleVector hardpro(4);
  hardpro[0] = subp->incoming().first; 
  hardpro[1] = subp->incoming().second;
  hardpro[2] = subp->outgoing()[0]; 
  hardpro[3] = subp->outgoing()[1];
  //particle ordering
  if(hardpro[0] < 0 && (hardpro[0]->id() != hardpro[1]->id())) 
    swap(hardpro[0], hardpro[1]);
  if(hardpro[2]->id() < 0 && (hardpro[2]->id() != hardpro[3]->id())) 
    swap(hardpro[2], hardpro[3]);
  //pick which process we are doing
  //qqbar initial state
  if( hardpro[1]->id() < 0 && (hardpro[0]->id() != hardpro[1]->id()) ) {
    SpinorVector spA, spB;
    SpinorBarVector spbA, spbB;
    SpinorWaveFunction(spA, hardpro[0], incoming, false, true);
    SpinorBarWaveFunction(spbA, hardpro[1], incoming, false, true);
    SpinorBarWaveFunction(spbB, hardpro[2], outgoing,true, true);
    SpinorWaveFunction(spB, hardpro[3], outgoing, true, true);
    //majorana
    if(!hardpro[2]->dataPtr()->CC() || hardpro[2]->id() == 1000024 || 
       hardpro[2]->id() == 1000037) {
      SpinorVector spC;
      SpinorBarVector spbC;
      for(unsigned int ix=0;ix<2;++ix) {
	spC.push_back(SpinorWaveFunction(-spbB[ix].getMomentum(),spbB[ix].getParticle(),
					 spbB[ix].wave().bar().conjugate(),
					 spbB[ix].direction()));
	spbC.push_back(SpinorBarWaveFunction(-spB[ix].getMomentum(),spB[ix].getParticle(),
					     spB[ix].wave().bar().conjugate(),
					     spB[ix].direction()));
      }
      double me;
      ProductionMatrixElement prodME = ffb2mfmfHeME(spA, spbA, spbB, spB, spC, 
						    spbC, me);
      HardVertexPtr hardvertex = new_ptr(HardVertex());
      hardvertex->ME(prodME);
      for(ParticleVector::size_type i = 0; i < 4; ++i) {
	dynamic_ptr_cast<SpinfoPtr>(hardpro[i]->spinInfo())->
	  setProductionVertex(hardvertex);
      }
    }
    //ffbar->ffbar
    else {
      double dummy;
      ProductionMatrixElement prodME = ffb2ffbHeME(spA, spbA, spbB, spB,
						   dummy);
      HardVertexPtr hardvertex = new_ptr(HardVertex());
      hardvertex->ME(prodME);
      for(ParticleVector::size_type i = 0; i < 4; ++i)
	dynamic_ptr_cast<SpinfoPtr>(hardpro[i]->spinInfo())->
	  setProductionVertex(hardvertex);
    }
  }   
  else {
    SpinorVector spA, spB;
    SpinorBarVector spbA, spbB;
    SpinorWaveFunction(spA,hardpro[0],incoming, false, true);
    SpinorWaveFunction(spB,hardpro[1],incoming, false, true);
    SpinorBarWaveFunction(spbA, hardpro[2], outgoing, true, true);
    SpinorBarWaveFunction(spbB, hardpro[3], outgoing, true, true);
    double dummy;
    ProductionMatrixElement prodME = ff2ffHeME(spA, spB, spbA, spbB,
					       dummy);
    HardVertexPtr hardvertex = new_ptr(HardVertex());
    hardvertex->ME(prodME);
    for(ParticleVector::size_type i = 0; i < 4; ++i)
      dynamic_ptr_cast<SpinfoPtr>(hardpro[i]->spinInfo())->
	setProductionVertex(hardvertex);
  }

}

void MEff2ff::persistentOutput(PersistentOStream & os) const {
  os << theScaV << theVecV << theTenV;
}

void MEff2ff::persistentInput(PersistentIStream & is, int) {
  is >> theScaV >> theVecV >> theTenV;
}

ClassDescription<MEff2ff> MEff2ff::initMEff2ff;
// Definition of the static class description member.

void MEff2ff::Init() {

  static ClassDocumentation<MEff2ff> documentation
    ("This is the implementation of the matrix element for fermion-"
     "antifermion -> fermion-antifermion.");

}
