// -*- C++ -*-
//
// MEfv2vf.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEfv2vf class.
//

#include "MEfv2vf.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::VectorWaveFunction;

void MEfv2vf::doinit() throw(InitException) {
  GeneralHardME::doinit();
  HPCount ndiags = numberOfDiags();
  theFerm.resize(ndiags);
  theVec.resize(ndiags);
  for(HPCount ix = 0; ix < ndiags; ++ix) {
    HPDiagram diagram = getProcessInfo()[ix];
    PDT::Spin offspin = diagram.intermediate->iSpin();
    if(diagram.channelType == HPDiagram::sChannel ||
       ( diagram.channelType == HPDiagram::tChannel 
	 && offspin == PDT::Spin1Half)) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
		(diagram.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(diagram.vertices.second);
      theFerm[ix] = make_pair(vert1, vert2);
    }
    else {
      if(offspin == PDT::Spin1) {
	AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	  (diagram.vertices.first);
	AbstractVVVVertexPtr vert2 = dynamic_ptr_cast<AbstractVVVVertexPtr>
	  (diagram.vertices.second);
	theVec[ix] = make_pair(vert1, vert2);
      }
    }
  }
}

double MEfv2vf::me2() const {
  //wavefunctions
  SpinorVector sp(2);
  VBVector vecIn(2), vecOut(3);
  SpinorBarVector spb(2);
  double fullme(0.);
  bool mc = !(mePartonData()[2]->mass() > 0.*MeV);
  if(mePartonData()[0]->id() > 0) {
    for(unsigned int i = 0; i < 2; ++i) {
      sp[i] = SpinorWaveFunction(rescaledMomenta()[0], mePartonData()[0], i, 
				 incoming);
      vecIn[i] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*i, 
				    incoming);
      vecOut[2*i] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 2*i, 
				     outgoing);
      spb[i] = SpinorBarWaveFunction(rescaledMomenta()[3], mePartonData()[3], i, 
				     outgoing);
    }
    if( !mc )
      vecOut[1] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 1, 
				     outgoing);
    fv2vfHeME(sp, vecIn, vecOut, mc, spb, fullme);
  }
  else {
    for(unsigned int i = 0; i < 2; ++i) {
      spb[i] = SpinorBarWaveFunction(rescaledMomenta()[0], mePartonData()[0], i, 
				     incoming);
      vecIn[i] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*i, 
				    incoming);
      vecOut[2*i] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 2*i, 
				     outgoing);
      sp[i] = SpinorWaveFunction(rescaledMomenta()[3], mePartonData()[3], i, 
				 outgoing);
    }
    if( !mc )
      vecOut[1] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 1, 
				     outgoing);
    fbv2vfbHeME(spb, vecIn, vecOut, mc, sp, fullme);
  }

#ifndef NDEBUG
  if( debugME() ) debug(fullme);
#endif

  return fullme;
}

ProductionMatrixElement
MEfv2vf::fv2vfHeME(const SpinorVector & spIn,  const VBVector & vecIn, 
		   const VBVector & vecOut, bool mc, 
		   const SpinorBarVector & spbOut, double & mesq) const {
  const HPCount ndiags = numberOfDiags();
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors(getColourFactors());
  const Energy2 q2(scale());
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
	if(mc && ovh == 1) ++ovh;
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
MEfv2vf::fbv2vfbHeME(const SpinorBarVector & spbIn, const VBVector & vecIn, 
		     const VBVector & vecOut, bool mc, 
		     const SpinorVector & spOut, double & mesq) const {
  const HPCount ndiags = numberOfDiags();
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors(getColourFactors());
  const Energy2 q2(scale());
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
	if(mc && ovh == 1) ++ovh;
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

void MEfv2vf::constructVertex(tSubProPtr sub) {
  ParticleVector ext(4);
  ext[0] = sub->incoming().first;
  ext[1] = sub->incoming().second;
  ext[2] = sub->outgoing()[0];
  ext[3] = sub->outgoing()[1];

  if( ext[0]->data().iSpin() > ext[1]->data().iSpin() ) swap(ext[0], ext[1]);
  if( ext[2]->data().iSpin() < ext[3]->data().iSpin() ) swap(ext[2], ext[3]);

  VBVector v1, v3;
  VectorWaveFunction(v1, ext[1], incoming, false, true, true);
  //function to calculate me2 expects massless incoming vectors
  // and this constructor sets the '1' polarisation at element [2] 
  //in the vector
  v1[1] = v1[2];
  bool mc = !(ext[2]->data().mass() > 0.*MeV);
  VectorWaveFunction(v3, ext[2], outgoing, true, mc, true);
  SpinorVector sp;  SpinorBarVector sbar;
  double dummy(0.);
  HardVertexPtr hv = new_ptr(HardVertex());
  //Need to use rescale momenta to calculate matrix element
  cPDVector data(4);
  vector<Lorentz5Momentum> momenta(4);
  for( size_t i = 0; i < 4; ++i ) {
    data[i] = ext[i]->dataPtr();
    momenta[i] = ext[i]->momentum();
  }
  rescaleMomenta(momenta, data);
  
  VectorWaveFunction vir(rescaledMomenta()[1], data[1], incoming);
  VectorWaveFunction vor(rescaledMomenta()[2], data[2], outgoing);

  if( ext[0]->id() > 0 ) {
    SpinorWaveFunction(sp, ext[0], incoming, false, true);
    SpinorBarWaveFunction(sbar, ext[3], outgoing, true, true);

    SpinorWaveFunction spr(rescaledMomenta()[0], data[0], incoming);
    SpinorBarWaveFunction sbr(rescaledMomenta()[3], data[3], outgoing);
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      spr.reset(ihel);
      sp[ihel] = spr;
      vir.reset(2*ihel);
      v1[ihel] = vir;
      vor.reset(2*ihel);
      v3[2*ihel] = vor;
      sbr.reset(ihel);
      sbar[ihel] = sbr;
    }
    if( !mc ) {
      vor.reset(1);
      v3[1] = vor;
    }
    ProductionMatrixElement pme = fv2vfHeME(sp, v1, v3, mc, sbar, dummy);
    hv->ME(pme);
    for(unsigned int i = 0; i < 4; ++i) 
      dynamic_ptr_cast<SpinfoPtr>(ext[i]->spinInfo())->setProductionVertex(hv);
  }
  else {
    SpinorBarWaveFunction(sbar, ext[0], incoming, false, true);
    SpinorWaveFunction(sp, ext[3], outgoing, true, true);

    SpinorBarWaveFunction sbr(rescaledMomenta()[0], data[0], incoming);
    SpinorWaveFunction spr(rescaledMomenta()[3], data[3], outgoing);
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      sbr.reset(ihel);
      sbar[ihel] = sbr;
      vir.reset(2*ihel);
      v1[ihel] = vir;
      vor.reset(2*ihel);
      v3[2*ihel] = vor;
      spr.reset(ihel);
      sp[ihel] = spr;
    }
    if( !mc ) {
      vor.reset(1);
      v3[1] = vor;
    }
    ProductionMatrixElement pme = fbv2vfbHeME(sbar, v1, v3, mc, sp, dummy);
    hv->ME(pme);
    for(unsigned int i = 0; i < 4; ++i) 
      dynamic_ptr_cast<SpinfoPtr>(ext[i]->spinInfo())->setProductionVertex(hv);
  }
  
#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif
}

void MEfv2vf::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id1 = abs(mePartonData()[0]->id());
  long id4 = abs(mePartonData()[3]->id());
  if( (id1 != 1 && id1 != 2) || mePartonData()[1]->id() != 21 ||
      mePartonData()[2]->id() != 5100021 || 
      (id4 != 5100001 && id4 != 5100002 &&
       id4 != 6100001 && id4 != 6100002) ) return;
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  Energy2 s(sHat());
  Energy2 mf2 = meMomenta()[2].m2();
  //  Energy4 spt2 = uHat()*tHat() - sqr(mf2);
  //swap t and u as formula defines process vf->vf
  Energy2 t3(uHat() - mf2), u4(tHat() - mf2);
  Energy4 s2(sqr(s)), t3s(sqr(t3)), u4s(sqr(u4));

  double analytic = -gs4*( 5.*s2/12./t3s + s2*s/t3s/u4 + 11.*s*u4/6./t3s
			   + 5.*u4s/12./t3s + u4s*u4/s/t3s)/3.;
  double diff = abs(analytic - me2);
  if( diff > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << "  ratio: " << analytic/me2 << '\n';
  }


}
