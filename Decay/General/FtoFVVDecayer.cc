// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FtoFVVDecayer class.
//

#include "FtoFVVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"
#include <numeric>

using namespace Herwig;

IBPtr FtoFVVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FtoFVVDecayer::fullclone() const {
  return new_ptr(*this);
}

void FtoFVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _sca << _fer << _vec << _ten;
}

void FtoFVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _sca >> _fer >> _vec >> _ten;
}

ClassDescription<FtoFVVDecayer> FtoFVVDecayer::initFtoFVVDecayer;
// Definition of the static class description member.

void FtoFVVDecayer::Init() {

  static ClassDocumentation<FtoFVVDecayer> documentation
    ("The FtoFVVDecayer class implements the general decay of a fermion to "
     "a fermion and a pair of vectors.");

}

void FtoFVVDecayer::doinit() throw(InitException) {
  GeneralThreeBodyDecayer::doinit();
  unsigned int ndiags = getProcessInfo().size();
  _sca.resize(ndiags);
  _fer.resize(ndiags);
  _vec.resize(ndiags);
  _ten.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractVVSVertexPtr vert2 = dynamic_ptr_cast<AbstractVVSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in FtoFVVDecayer::doinit()"
	<< Exception::runerror;
      _sca[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1Half) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in FtoFVVDecayer::doinit()"
	<< Exception::runerror;
      _fer[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.first);
      AbstractVVVVertexPtr vert2 = dynamic_ptr_cast<AbstractVVVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in FtoFVVDecayer::doinit()"
	<< Exception::runerror;
      _vec[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractFFTVertexPtr vert1 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.first);
      AbstractVVTVertexPtr vert2 = dynamic_ptr_cast<AbstractVVTVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a tensor diagram in FtoFVVDecayer::doinit()"
	<< Exception::runerror;
      _ten[ix] = make_pair(vert1, vert2);
    }
  }
}

double  FtoFVVDecayer::me2(bool vertex, const int ichan, const Particle & inpart,
			   const ParticleVector & decay) const {
  // spin density matrix
  RhoDMatrix rhoin(PDT::Spin1Half);
  
  //Set up wave-functions
  bool ferm( inpart.id() > 0 );
  //incoming
  vector<SpinorWaveFunction> fwave;
  vector<SpinorBarWaveFunction> fbwave;
  if( ferm ) {
    SpinorWaveFunction(fwave, rhoin, const_ptr_cast<tPPtr>(&inpart),
		       Helicity::incoming, true, vertex);
    if( fwave[0].wave().Type() != u_spinortype )
      fwave[0].conjugate();
    if( fwave[1].wave().Type() != u_spinortype )
      fwave[1].conjugate();
  }
  else {
    SpinorBarWaveFunction(fbwave, rhoin, const_ptr_cast<tPPtr>(&inpart),
			  Helicity::incoming, true, vertex);
    if( fbwave[0].wave().Type() != v_spinortype )
      fbwave[0].conjugate();
    if( fbwave[1].wave().Type() != v_spinortype )
      fbwave[1].conjugate();
  }

  // outgoing, keep track of fermion and first occurrence of vector positions
  int isp(-1), ivec(-1);
  pair<vector<VectorWaveFunction>, vector<VectorWaveFunction> > vwave;
  // outgoing particles
  for(int ix = 0; ix < 3; ++ix) {
    tPPtr p = decay[ix];
    if( p->dataPtr()->iSpin() == PDT::Spin1Half ) {
      isp = ix;
      if( ferm ) {
	SpinorBarWaveFunction(fbwave, p, Helicity::outgoing, true, vertex);
	if( fbwave[0].wave().Type() != u_spinortype )
	  fbwave[0].conjugate();
	if( fbwave[1].wave().Type() != u_spinortype )
	  fbwave[1].conjugate();
      }
      else {
	SpinorWaveFunction(fwave, p, Helicity::outgoing, true, vertex);
	if( fwave[0].wave().Type() != v_spinortype )
	  fwave[0].conjugate();
	if( fwave[1].wave().Type() != v_spinortype )
	  fwave[1].conjugate();
      }
    }
    else if( p->dataPtr()->iSpin() == PDT::Spin1 ) {
      if( ivec < 0 ) {
	ivec = ix;
	VectorWaveFunction(vwave.first, p, Helicity::outgoing, true, false, 
			   vertex);
      }
      else {
	VectorWaveFunction(vwave.second, p, Helicity::outgoing, true, false, 
			   vertex);
      }
    }
  }
  assert(isp >= 0 && ivec >= 0);
  Energy2 scale(sqr(inpart.mass()));
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  const size_t ncf(numberOfFlows());
  vector<Complex> flows(ncf, Complex(0.)), largeflows(ncf, Complex(0.)); 
  vector<DecayMatrixElement> 
    mes(ncf,DecayMatrixElement(PDT::Spin1Half, 
			       (isp == 0) ? PDT::Spin1Half : PDT::Spin1,
			       (isp == 1) ? PDT::Spin1Half : PDT::Spin1,
			       (isp == 2) ? PDT::Spin1Half : PDT::Spin1));
  vector<DecayMatrixElement> 
    mel(ncf,DecayMatrixElement(PDT::Spin1Half, 
			       (isp == 0) ? PDT::Spin1Half : PDT::Spin1,
			       (isp == 1) ? PDT::Spin1Half : PDT::Spin1,
			       (isp == 2) ? PDT::Spin1Half : PDT::Spin1));
  //Helicity calculation
  for( unsigned int if1 = 0; if1 < 2; ++if1 ) {
    for( unsigned int if2 = 0; if2 < 2; ++if2 ) {
      for( unsigned int iv1 = 0; iv1 < 3; ++iv1 ) {
	for( unsigned int iv2 = 0; iv2 < 3; ++iv2 ) {
	flows = vector<Complex>(ncf, Complex(0.));
	largeflows = vector<Complex>(ncf, Complex(0.));
	unsigned int idiag(0);
	for(vector<TBDiagram>::const_iterator dit = getProcessInfo().begin();
	    dit != getProcessInfo().end(); ++dit) {
	  // If we are selecting a particular channel
	  if( ichan >= 0 && diagramMap()[ichan] != idiag ) {
	    ++idiag;
	    continue;
	  }
	  tcPDPtr offshell = (*dit).intermediate;
	  if( offshell->CC() ) offshell = offshell->CC();
	  Complex diag;
	  if( offshell->iSpin() == PDT::Spin1Half ) {
	    //Need make sure we connect the correct particles 
	    VectorWaveFunction vw1, vw2;
	    if( (*dit).channelType == TBDiagram::channel23 ) {
	      vw1 = vwave.first[iv1];
	      vw2 = vwave.second[iv2];
	    }
	    else if( (*dit).channelType == TBDiagram::channel12 ) {
	      vw1 = vwave.second[iv2];
	      vw2 = vwave.first[iv1];
	    }
	    else {
	      if( ivec < isp ) {
		vw1 = vwave.second[iv2];
		vw2 = vwave.first[iv1];
	      }
	      else {
		vw1 = vwave.first[iv1];
		vw2 = vwave.second[iv2];
	      }
	    }
	    if( ferm ) {
	      SpinorWaveFunction inters = 
		_fer[idiag].first->evaluate(scale, widthOption(), offshell,
					    fwave[if1], vw1);
	      diag = _fer[idiag].second->evaluate(scale, inters, fbwave[if2],
						  vw2);
	    }
	    else {
	      SpinorBarWaveFunction inters = 
		_fer[idiag].first->evaluate(scale, widthOption(), offshell,
					    fbwave[if2], vw1);
	      diag = _fer[idiag].second->evaluate(scale, fwave[if1], inters, 
						  vw2);
	    }
	  }
	  else if( offshell->iSpin() == PDT::Spin0 ) {
	    ScalarWaveFunction inters = 
	      _sca[idiag].first->evaluate(scale, widthOption(), offshell, 
					  fwave[if1], fbwave[if2]);
	    diag = _sca[idiag].second->evaluate(scale, vwave.first[iv1],
						vwave.second[iv2], inters);
	  }
	  else if( offshell->iSpin() == PDT::Spin1 ) {
	    VectorWaveFunction interv = 
	      _vec[idiag].first->evaluate(scale, widthOption(), offshell, 
					  fwave[if1], fbwave[if2]);
	    diag = _vec[idiag].second->evaluate(scale, vwave.first[iv1],
						vwave.second[iv2], interv);
	  } 
	  else if( offshell->iSpin() == PDT::Spin2 ) {
	    TensorWaveFunction intert = 
	      _ten[idiag].first->evaluate(scale, widthOption(), offshell, 
					  fwave[if1], fbwave[if2]);
	    diag = _ten[idiag].second->evaluate(scale, vwave.first[iv1],
						vwave.second[iv2], intert);
	  }
	  else 
	    throw Exception()
	      << "Unknown intermediate in FtoFVVDecayer::me2()" 
	      << Exception::runerror;
	  //NO sign
	  if( !ferm ) diag *= -1;

	  if(ichan<0) {
	    for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
	      flows[dit->colourFlow[iy].first - 1] += 
		dit->colourFlow[iy].second * diag;
	    }
	    for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
	      largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		dit->largeNcColourFlow[iy].second * diag;
	    }
	  }
	  else {
	    for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
	      if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
	      flows[dit->colourFlow[iy].first - 1] += 
		dit->colourFlow[iy].second * diag;
	    }
	    for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
	      if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
	      largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		dit->largeNcColourFlow[iy].second * diag;
	    }
	  }
	  ++idiag;
	}// end diagram loop

	// now add the flows to the me2 
	unsigned int h1(if1), h2(if2);
	if( !ferm ) swap(h1,h2);
	  for(unsigned int ix = 0; ix < ncf; ++ix) {
	    if(isp == 0) {
	      mes[ix](h1, h2, iv1, iv2) = flows[ix];
	      mel[ix](h1, h2, iv1, iv2) = largeflows[ix];
	    }
	    else if(isp == 1) { 
	      mes[ix](h1, iv1, h2, iv2) = flows[ix];
	      mel[ix](h1, iv1, h2, iv2) = largeflows[ix];
	    }
	    else if(isp == 2) { 
	      mes[ix](h1, iv1, iv2, h2) = flows[ix];
	      mel[ix](h1, iv1, h2, iv2) = largeflows[ix];
	    }
	  }

	}//v2
      }//v1
    }//f2
  }//f1
  
  //Finally, work out me2. This depends on whether we are selecting channels
  //or not
  double me2(0.);
  if(ichan<0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix].contract(mes[iy],rhoin)).real();
	me2 += con;
	if(ix==iy) {
	  con = nfactors[ix][iy]*(mel[ix].contract(mel[iy],rhoin)).real();
	  pflows[ix] += con;
	}
      }
    }
    if(vertex) {
      double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
      ptotal *= UseRandom::rnd();
      for(unsigned int ix=0;ix<pflows.size();++ix) {
	if(ptotal<=pflows[ix]) {
	  colourFlow(ix);
	  ME(mes[ix]);
	  break;
	}
	ptotal-=pflows[ix];
      }
    }
  }
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow].contract(mel[iflow],rhoin)).real();
  }
  // return the matrix element squared

  return me2;
}

WidthCalculatorBasePtr FtoFVVDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<FtoFVVDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),
		  outgoing()[2]->mass()));
}
