// -*- C++ -*-
//
// MEff2ff.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2ff class.
//

#include "MEff2ff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;
using ThePEG::Helicity::SpinfoPtr;

void MEff2ff::doinit() throw(InitException) {
  GeneralHardME::doinit();
  size_t ndiags = numberOfDiags();
  theScaV.resize(ndiags);
  theVecV.resize(ndiags);
  theTenV.resize(ndiags);
  for(size_t ix = 0;ix < ndiags; ++ix) {
    HPDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      theScaV[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      theVecV[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractFFTVertexPtr vert1 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.first);
      AbstractFFTVertexPtr vert2 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.second);
      theTenV[ix] = make_pair(vert1, vert2);
    }
  }
}

double MEff2ff::me2() const {
  tcPDPtr ina(mePartonData()[0]), inb(mePartonData()[1]),
    outa(mePartonData()[2]), outb(mePartonData()[3]);
  bool majorana(false);
  if( (!outa->CC() && !outb->CC() ) || 
      ((abs(outa->id()) > 1000000 && abs(outa->id()) < 2000000) &&
       (abs(outb->id()) > 1000000 && abs(outb->id()) < 2000000)) )
    majorana = true;
  
  double full_me(0.);
  vector<SpinorWaveFunction> spA(2), spB(2);
  vector<SpinorBarWaveFunction> spbA(2), spbB(2);
  if( ina->id() > 0 && inb->id() < 0) {
    for(unsigned int ih = 0; ih < 2; ++ih) {
      spA[ih] = SpinorWaveFunction(rescaledMomenta()[0], ina, ih, 
				   incoming);
      spbA[ih] = SpinorBarWaveFunction(rescaledMomenta()[1], inb, ih, 
				       incoming);
      spB[ih] = SpinorWaveFunction(rescaledMomenta()[3], outb, ih, outgoing);
      spbB[ih] = SpinorBarWaveFunction(rescaledMomenta()[2], outa, ih, outgoing);
    }
    if(majorana) {
      vector<SpinorWaveFunction> spC(2);
      vector<SpinorBarWaveFunction> spbC(2);
      for(unsigned int ih = 0; ih < 2; ++ih) {
	spC[ih] = SpinorWaveFunction(rescaledMomenta()[2], outa, ih, outgoing);
	spbC[ih] = SpinorBarWaveFunction(rescaledMomenta()[3], outb, ih, outgoing);
      }
      ffb2mfmfHeME(spA, spbA, spbB, spB, spC, spbC, full_me);
      SpinorWaveFunction spOut2(rescaledMomenta()[2], outa, outgoing);
      SpinorBarWaveFunction spbarOut2(rescaledMomenta()[3], outb, outgoing);
     }
    else {
      ffb2ffbHeME(spA, spbA, spbB, spB, full_me); 
    }
  }
  else if( ina->id() > 0 && inb->id() > 0 ) {
    SpinorVector spA(2), spB(2);
    SpinorBarVector spbA(2), spbB(2);
    for(unsigned int ih = 0; ih < 2; ++ih) {
      spA[ih] = SpinorWaveFunction(rescaledMomenta()[0], ina, ih,
				   incoming);
      spB[ih] = SpinorWaveFunction(rescaledMomenta()[1], inb, ih,
				   incoming);
      spbA[ih] = SpinorBarWaveFunction(rescaledMomenta()[2], outa, ih,
				       outgoing);
      spbB[ih] = SpinorBarWaveFunction(rescaledMomenta()[3], outb, ih,
				       outgoing);
    }
    ff2ffHeME(spA, spB, spbA, spbB, full_me);
  }
  else if( ina->id() < 0 && inb->id() < 0 ) {
    SpinorVector spA(2), spB(2);
    SpinorBarVector spbA(2), spbB(2);
    for(unsigned int ih = 0; ih < 2; ++ih) {
      spbA[ih] = SpinorBarWaveFunction(rescaledMomenta()[0], ina, ih,
				       incoming);
      spbB[ih] = SpinorBarWaveFunction(rescaledMomenta()[1], inb, ih,
				   incoming);
      spA[ih] = SpinorWaveFunction(rescaledMomenta()[2], outa, ih,
				   outgoing);
      spB[ih] = SpinorWaveFunction(rescaledMomenta()[3], outb, ih,
				   outgoing);
    }
    fbfb2fbfbHeME(spbA, spbB, spA, spB, full_me);
  }
  else 
    throw MEException() 
      << "MEff2ff::me2() - Cannot find correct function to deal with process " 
      << ina->PDGName() << "," << inb->PDGName() << "->" << outa->PDGName() 
      << "," << outb->PDGName() << "\n";

#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif

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
MEff2ff::fbfb2fbfbHeME(SpinorBarVector & fbin, SpinorBarVector & fbin2,
		       SpinorVector & fout, SpinorVector & fout2,
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
							fout2[ofhel2],
							fbin2[ifhel2]);
		  diag[ix] = theScaV[ix].first->evaluate(q2, fout[ofhel1],
							 fbin[ifhel1], interS);
		}
		else {
		  interS = theScaV[ix].second->evaluate(q2, 3, offshell, 
							fout[ofhel1],
							fbin2[ifhel2]);
		  diag[ix] = -theScaV[ix].first->evaluate(q2, fout2[ofhel2],
							 fbin[ifhel1], interS);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  interV = theVecV[ix].second->evaluate(q2, 3, offshell, 
							fout2[ofhel2],
							fbin2[ifhel2]);
		  diag[ix] = theVecV[ix].first->evaluate(q2, fout[ofhel1],
							 fbin[ifhel1], interV);
		}
		else {
		  interV = theVecV[ix].second->evaluate(q2, 3, offshell, 
							fout[ofhel1],
							fbin2[ifhel2]);
		  diag[ix] = -theVecV[ix].first->evaluate(q2, fout2[ofhel2],
							 fbin[ifhel1], interV);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		if(current.ordered.second) {
		  interT = theTenV[ix].second->evaluate(q2, 3, offshell, 
							fout2[ofhel2],
							fbin2[ifhel2]);
		  diag[ix] = theTenV[ix].first->evaluate(q2, fout[ofhel1],
							 fbin[ifhel1], interT);
		}
		else {
		  interT = theTenV[ix].second->evaluate(q2, 3, offshell, 
							fout[ofhel1],
							fbin2[ifhel2]);
		  diag[ix] = -theTenV[ix].first->evaluate(q2, fout2[ofhel2],
							 fbin[ifhel1], interT);
		}
	      }
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
  const double colfact = (mePartonData()[0]->iColour() == PDT::Colour3bar) 
    ? 1./9. : 1;
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = 0.25*identfact*colfact*me[ix];
  meInfo(save);
  me2 *= 0.25*identfact*colfact;
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
  const double identfact = mePartonData()[2]->id() == mePartonData()[3]->id()  
    ? 0.5 : 1.; 
  const double colfact = mePartonData()[0]->coloured() ? 1./9. : 1; 
  DVector save(ndiags); 
  for(DVector::size_type ix = 0; ix < ndiags; ++ix) 
    save[ix] = 0.25*identfact*colfact*me[ix]; 
  meInfo(save); 
  me2 *= 0.25*identfact*colfact;
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
  cf[3] = ColourLines("1 5, -5 2 -3");
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
	  cl = 0;
	else if( outbc == PDT::Colour8 )
	  cl = current.ordered.second ? 2 : 3;
	else
	  cl = current.ordered.second ? 20 : 21;
      }
      else
	cl = 1;
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

  //ensure particle ordering is the same as it was when
  //the diagrams were created
  if( hardpro[0]->id() != getIncoming().first )
    swap(hardpro[0], hardpro[1]);
  if( hardpro[2]->id() != getOutgoing().first )
    swap(hardpro[2], hardpro[3]);

  //Need to use rescale momenta to calculate matrix element
  cPDVector data(4);
  vector<Lorentz5Momentum> momenta(4);
  for( size_t i = 0; i < 4; ++i ) {
    data[i] = hardpro[i]->dataPtr();
    momenta[i] = hardpro[i]->momentum();
  }
  rescaleMomenta(momenta, data);

  double dummy(0.);
  //pick which process we are doing
  if( hardpro[0]->id() > 0) {
    //common spinors
    SpinorVector spA;
    SpinorBarVector spbB;
    SpinorWaveFunction(spA, hardpro[0], incoming, false, true);
    SpinorBarWaveFunction(spbB, hardpro[2], outgoing,true, true);

    //ME spinors
    SpinorWaveFunction sp1r(rescaledMomenta()[0], data[0], incoming);
    SpinorBarWaveFunction spb2r(rescaledMomenta()[2], data[2], outgoing);

    //majorana
    if(!hardpro[2]->dataPtr()->CC() || hardpro[2]->id() == 1000024 || 
       hardpro[2]->id() == 1000037) {
      SpinorVector spB, spC;
      SpinorBarVector spbA, spbC;
      SpinorBarWaveFunction(spbA, hardpro[1], incoming, false, true);
      SpinorWaveFunction(spB, hardpro[3], outgoing, true, true);
      //ME spinors
      SpinorBarWaveFunction spb1r(rescaledMomenta()[1], data[1], incoming);
      SpinorWaveFunction sp2r(rescaledMomenta()[3], data[3], outgoing);

      for( unsigned int ihel = 0; ihel < 2; ++ihel ) {
	sp1r.reset(ihel);
	spA[ihel] = sp1r;
	spb1r.reset(ihel);
	spbA[ihel] = spb1r;
	
	spb2r.reset(ihel);
	spbB[ihel] = spb2r;
	sp2r.reset(ihel);
	spB[ihel] = sp2r;

	//extra spinors
	spC.push_back(SpinorWaveFunction(-spbB[ihel].getMomentum(),
					 spbB[ihel].getParticle(),
					 spbB[ihel].wave().bar().conjugate(),
					 spbB[ihel].direction()));
	spbC.push_back(SpinorBarWaveFunction(-spB[ihel].getMomentum(),
					     spB[ihel].getParticle(),
					     spB[ihel].wave().bar().conjugate(),
					     spB[ihel].direction()));
      }
      ProductionMatrixElement prodME = ffb2mfmfHeME(spA, spbA, spbB, spB, spC, 
						    spbC, dummy);
      HardVertexPtr hardvertex = new_ptr(HardVertex());
      hardvertex->ME(prodME);
      for(ParticleVector::size_type i = 0; i < 4; ++i) {
	dynamic_ptr_cast<SpinfoPtr>(hardpro[i]->spinInfo())->
	  setProductionVertex(hardvertex);
      }
    }
    //ffbar->ffbar
    else if( hardpro[1]->id() < 0 ) {
      SpinorVector spB;
      SpinorBarVector spbA;
      SpinorBarWaveFunction(spbA, hardpro[1], incoming, false, true);
      SpinorWaveFunction(spB, hardpro[3], outgoing, true, true);

      //ME spinors
      SpinorBarWaveFunction spb1r(rescaledMomenta()[1], data[1], incoming);
      SpinorWaveFunction sp2r(rescaledMomenta()[3], data[3], outgoing);

      for( unsigned int ihel = 0; ihel < 2; ++ihel ) {
	sp1r.reset(ihel);
	spA[ihel] = sp1r;
	spb1r.reset(ihel);
	spbA[ihel] = spb1r;
	
	spb2r.reset(ihel);
	spbB[ihel] = spb2r;
	sp2r.reset(ihel);
	spB[ihel] = sp2r;
      }
      
      ProductionMatrixElement prodME = ffb2ffbHeME(spA, spbA, spbB, spB,
						   dummy);
      HardVertexPtr hardvertex = new_ptr(HardVertex());
      hardvertex->ME(prodME);
      for(ParticleVector::size_type i = 0; i < 4; ++i)
	dynamic_ptr_cast<SpinfoPtr>(hardpro[i]->spinInfo())->
	  setProductionVertex(hardvertex);
    }
    //ff2ff
    else {
      SpinorVector spB;
      SpinorBarVector spbA;
      SpinorWaveFunction(spB,hardpro[1],incoming, false, true);
      SpinorBarWaveFunction(spbA, hardpro[3], outgoing, true, true);

      SpinorWaveFunction sp2r(rescaledMomenta()[1], data[1], incoming);
      SpinorBarWaveFunction spb1r(rescaledMomenta()[3], data[3], outgoing);

      for( unsigned int ihel = 0; ihel < 2; ++ihel ) {
	sp1r.reset(ihel);
	spA[ihel] = sp1r;
	sp2r.reset(ihel);
	spB[ihel] = sp2r;

	spb2r.reset(ihel);
	spbB[ihel] = spb2r;
	spb1r.reset(ihel);
	spbA[ihel] = spb1r;
	
      }

      ProductionMatrixElement prodME = ff2ffHeME(spA, spB, spbB, spbA,
						 dummy);
      HardVertexPtr hardvertex = new_ptr(HardVertex());
      hardvertex->ME(prodME);
      for(ParticleVector::size_type i = 0; i < 4; ++i)
	dynamic_ptr_cast<SpinfoPtr>(hardpro[i]->spinInfo())->
	  setProductionVertex(hardvertex);
    }
  } 
  //fbarfbar->fbarfbar
  else {
    SpinorVector spA, spB;
    SpinorBarVector spbA, spbB;
    SpinorBarWaveFunction(spbA,hardpro[0],incoming, false, true);
    SpinorBarWaveFunction(spbB,hardpro[1],incoming, false, true);
    SpinorWaveFunction(spA, hardpro[2], outgoing, true, true);
    SpinorWaveFunction(spB, hardpro[3], outgoing, true, true);

    //ME spinors
    SpinorBarWaveFunction spb1r(rescaledMomenta()[0], data[0], incoming);
    SpinorBarWaveFunction spb2r(rescaledMomenta()[1], data[1], incoming);
    SpinorWaveFunction sp1r(rescaledMomenta()[2], data[2], outgoing);
    SpinorWaveFunction sp2r(rescaledMomenta()[3], data[3], outgoing);
    
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {
      spb1r.reset(ihel);
      spbA[ihel] = spb1r;
      spb2r.reset(ihel);
      spbB[ihel] = spb2r;
      
      sp1r.reset(ihel);
      spA[ihel] = sp1r;
      sp2r.reset(ihel);
      spB[ihel] = sp2r;
    }
    
    ProductionMatrixElement prodME = fbfb2fbfbHeME(spbA, spbB, spA, spB,
						   dummy);
    HardVertexPtr hardvertex = new_ptr(HardVertex());
    hardvertex->ME(prodME);
    for(ParticleVector::size_type i = 0; i < 4; ++i)
      dynamic_ptr_cast<SpinfoPtr>(hardpro[i]->spinInfo())->
	setProductionVertex(hardvertex);
  }
  
#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif

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

void MEff2ff::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id1 = mePartonData()[0]->id();
  long id2 = mePartonData()[1]->id();
  long id3 = mePartonData()[2]->id();
  long id4 = mePartonData()[3]->id();
  long aid1 = abs(mePartonData()[0]->id());
  long aid2 = abs(mePartonData()[1]->id());
  long aid3 = abs(mePartonData()[2]->id());
  long aid4 = abs(mePartonData()[3]->id());
  if( (aid1 != 1 && aid1 != 2) || (aid2 != 1 && aid2 != 2) ) return;
  double analytic(0.);
  if( id3 == id4 && id3 == 1000021 ) {
    tcSMPtr sm = generator()->standardModel();
    double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
    int Nc = sm->Nc();
    double Cf = (sqr(Nc) - 1.)/2./Nc;
    Energy2 mgo2 = meMomenta()[3].m2();
    long squark = (aid1 == 1) ? 1000001 : 1000002;
    Energy2 muL2 = sqr(getParticleData(squark)->mass());
    Energy2 deltaL = muL2 - mgo2;
    Energy2 muR2 = sqr(getParticleData(squark + 1000000)->mass());
    Energy2 deltaR = muR2 - mgo2;
    Energy2 s(sHat());
    Energy2 m3s = meMomenta()[2].m2();
    Energy2 m4s = meMomenta()[3].m2();
    Energy4 spt2 = uHat()*tHat() - m3s*m4s;
    Energy2 t3(tHat() - m3s), u4(uHat() - m4s);
    
    double Cl = 2.*spt2*( (u4*u4 - deltaL*deltaL) + (t3*t3 - deltaL*deltaL)
			  - (s*s/Nc/Nc) )/s/s/(u4 - deltaL)/(t3 - deltaL);
    Cl += deltaL*deltaL*( (1./sqr(t3 - deltaL)) + (1./sqr(u4 - deltaL))
			  - ( sqr( (1./(t3 - deltaL)) - 
				   (1./(u4 - deltaL)) )/Nc/Nc ) );
    
    double Cr = 2.*spt2*( (u4*u4 - deltaR*deltaR) + (t3*t3 - deltaR*deltaR)
			  - (s*s/Nc/Nc) )/s/s/(u4 - deltaR)/(t3 - deltaR);
    Cr += deltaR*deltaR*( (1./sqr(t3 - deltaR)) + (1./sqr(u4 - deltaR))
			  - ( sqr( (1./(t3 - deltaR)) 
				   - (1./(u4 - deltaR)) )/Nc/Nc ) );
    analytic = gs4*Cf*(Cl + Cr)/4.;
  }
  else if( (aid3 == 5100001 || aid3 == 5100002 ||
	    aid3 == 6100001 || aid3 == 6100002) &&
	   (aid4 == 5100001 || aid4 == 5100002 ||
	    aid4 == 6100001 || aid4 == 6100002) ) {
    tcSMPtr sm = generator()->standardModel();
    double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
    Energy2 s(sHat());
    Energy2 mf2 = meMomenta()[2].m2();
    Energy2 t3(tHat() - mf2), u4(uHat() - mf2);
    Energy4 s2(sqr(s)), t3s(sqr(t3)), u4s(sqr(u4));
    
    bool iflav = (aid2 - aid1 == 0);
    int alpha(aid3/1000000), beta(aid4/1000000);
    bool oflav = ((aid3 - aid1) % 10  == 0);
    if( alpha != beta ) {
      if( ( id1 > 0 && id2 > 0) ||
	  ( id1 < 0 && id2 < 0) ) {
	if( iflav )
	  analytic = gs4*( mf2*(2.*s2*s/t3s/u4s - 4.*s/t3/u4) 
			   + 2.*sqr(s2)/t3s/u4s - 8.*s2/t3/u4 + 5. )/9.;
	else
	  analytic = gs4*( -2.*mf2*(1./t3 + u4/t3s) + 0.5 + 2.*u4s/t3s)/9.;
      }
      else
	analytic = gs4*( 2.*mf2*(1./t3 + u4/t3s) + 5./2. + 4.*u4/t3 
			 + 2.*u4s/t3s)/9.;
    }
    else {
      if( ( id1 > 0 && id2 > 0) ||
	  ( id1 < 0 && id2 < 0) ) {
	if( iflav ) {
	  analytic = gs4*( mf2*(6.*t3/u4s + 6.*u4/t3s - s/t3/u4) 
			   + 2.*(3.*t3s/u4s + 3.*u4s/t3s 
				 + 4.*s2/t3/u4 - 5.) )/27.;
	}
	else
	  analytic = 2.*gs4*( -mf2*s/t3s + 0.25 + s2/t3s )/9.;
      }
      else {
	if( iflav ) {
	  if( oflav )	  
	    analytic = gs4*( 2.*mf2*(4./s + s/t3s - 1./t3) + 23./6.+ 2.*s2/t3s
			     + 8.*s/3./t3 + 6.*t3/s + 8.*t3s/s2 )/9.;
	  else
	    analytic = 4.*gs4*( 2.*mf2/s + (t3s + u4s)/s2)/9.;
	}
	else 
	  analytic = gs4*(4.*mf2*s/t3s + 5. + 4.*s2/t3s + 8.*s/t3 )/18.;
      }
    }
    if( id3 == id4 ) analytic /= 2.;

  }
  else return;
  double diff = abs(analytic - me2);
  if( diff  > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff  << "  ratio: " << analytic/me2 
      << '\n';
  }
}
