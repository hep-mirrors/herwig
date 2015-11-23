// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FtoFVVDecayer class.
//

#include "FtoFVVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
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

void FtoFVVDecayer::doinit() {
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

double  FtoFVVDecayer::me2(const int ichan, const Particle & inpart,
			   const ParticleVector & decay,
			   MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != inpart.id();
  // special handling or first/last call
  //Set up wave-functions
  bool ferm( inpart.id() > 0 );
  if(meopt==Initialize) {
    if( ferm ) {
      SpinorWaveFunction::
	calculateWaveFunctions(_fwave,_rho,const_ptr_cast<tPPtr>(&inpart),
			       Helicity::incoming);
      if( _fwave[0].wave().Type() != u_spinortype )
	_fwave[0].conjugate();
      if( _fwave[1].wave().Type() != u_spinortype )
	_fwave[1].conjugate();
    }
    else {
      SpinorBarWaveFunction::
	calculateWaveFunctions(_fbwave, _rho, const_ptr_cast<tPPtr>(&inpart),
			       Helicity::incoming);
      if( _fbwave[0].wave().Type() != v_spinortype )
	_fbwave[0].conjugate();
      if( _fbwave[1].wave().Type() != v_spinortype )
	_fbwave[1].conjugate();
    }
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(ferm)
      SpinorWaveFunction::constructSpinInfo(_fwave,
					    const_ptr_cast<tPPtr>(&inpart),
					    Helicity::incoming,true);
    else
      SpinorBarWaveFunction::constructSpinInfo(_fbwave,
					       const_ptr_cast<tPPtr>(&inpart),
					       Helicity::incoming,true);
    int ivec(-1);
    // outgoing particles
    for(int ix = 0; ix < 3; ++ix) {
      tPPtr p = decay[ix];
      if( p->dataPtr()->iSpin() == PDT::Spin1Half ) {
	if( ferm ) {
	  SpinorBarWaveFunction::
	    constructSpinInfo(_fbwave, p, Helicity::outgoing,true);
	}
	else {
	  SpinorWaveFunction::
	    constructSpinInfo(_fwave , p, Helicity::outgoing,true);
	}
      }
      else if( p->dataPtr()->iSpin() == PDT::Spin1 ) {
	if( ivec < 0 ) {
	  ivec = ix;
	  VectorWaveFunction::
	    constructSpinInfo(_vwave.first , p, Helicity::outgoing, true, false);
	}
	else {
	  VectorWaveFunction::
	    constructSpinInfo(_vwave.second, p, Helicity::outgoing, true, false);
	}
      }
    }
    return 0.;
  }
  // outgoing, keep track of fermion and first occurrence of vector positions
  int isp(-1), ivec(-1);
  // outgoing particles
  pair<bool,bool> mass = make_pair(false,false);
  for(int ix = 0; ix < 3; ++ix) {
    tPPtr p = decay[ix];
    if( p->dataPtr()->iSpin() == PDT::Spin1Half ) {
      isp = ix;
      if( ferm ) {
	SpinorBarWaveFunction::
	  calculateWaveFunctions(_fbwave, p, Helicity::outgoing);
	if( _fbwave[0].wave().Type() != u_spinortype )
	  _fbwave[0].conjugate();
	if( _fbwave[1].wave().Type() != u_spinortype )
	  _fbwave[1].conjugate();
      }
      else {
	SpinorWaveFunction::
	  calculateWaveFunctions(_fwave, p, Helicity::outgoing);
	if( _fwave[0].wave().Type() != v_spinortype )
	  _fwave[0].conjugate();
	if( _fwave[1].wave().Type() != v_spinortype )
	  _fwave[1].conjugate();
      }
    }
    else if( p->dataPtr()->iSpin() == PDT::Spin1 ) {
      bool massless = p->id() == ParticleID::gamma || p->id() == ParticleID::g;
      if( ivec < 0 ) {
	ivec = ix;
	VectorWaveFunction::
	  calculateWaveFunctions(_vwave.first , p, Helicity::outgoing, massless);
	mass.first = massless;
      }
      else {
	VectorWaveFunction::
	  calculateWaveFunctions(_vwave.second, p, Helicity::outgoing, massless);
	mass.second = massless;
      }
    }
  }
  assert(isp >= 0 && ivec >= 0);
  Energy2 scale(sqr(inpart.mass()));
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  const size_t ncf(numberOfFlows());
  vector<Complex> flows(ncf, Complex(0.)), largeflows(ncf, Complex(0.)); 
  vector<GeneralDecayMEPtr> 
    mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half, 
					      (isp == 0) ? PDT::Spin1Half : PDT::Spin1,
					      (isp == 1) ? PDT::Spin1Half : PDT::Spin1,
					      (isp == 2) ? PDT::Spin1Half : PDT::Spin1)));
  vector<GeneralDecayMEPtr> 
    mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half, 
					      (isp == 0) ? PDT::Spin1Half : PDT::Spin1,
					      (isp == 1) ? PDT::Spin1Half : PDT::Spin1,
					      (isp == 2) ? PDT::Spin1Half : PDT::Spin1)));
  //Helicity calculation
  for( unsigned int if1 = 0; if1 < 2; ++if1 ) {
    for( unsigned int if2 = 0; if2 < 2; ++if2 ) {
      for( unsigned int iv1 = 0; iv1 < 3; ++iv1 ) {
	if ( mass.first && iv1 == 1 ) continue;
	for( unsigned int iv2 = 0; iv2 < 3; ++iv2 ) {
	  if ( mass.second && iv2 == 1 ) continue;
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
	  if(cc&&offshell->CC()) offshell=offshell->CC();
	  Complex diag;
	  if( offshell->iSpin() == PDT::Spin1Half ) {
	    // Make sure we connect the correct particles 
	    VectorWaveFunction vw1, vw2;
	    if( (*dit).channelType == TBDiagram::channel23 ) {
	      vw1 = _vwave.first[iv1];
	      vw2 = _vwave.second[iv2];
	    }
	    else if( (*dit).channelType == TBDiagram::channel12 ) {
	      vw1 = _vwave.second[iv2];
	      vw2 = _vwave.first[iv1];
	    }
	    else {
	      if( ivec < isp ) {
		vw1 = _vwave.second[iv2];
		vw2 = _vwave.first[iv1];
	      }
	      else {
		vw1 = _vwave.first[iv1];
		vw2 = _vwave.second[iv2];
	      }
	    }
	    if( ferm ) {
	      SpinorWaveFunction inters = 
		_fer[idiag].first->evaluate(scale, widthOption(), offshell,
					    _fwave[if1], vw1);
	      diag = _fer[idiag].second->evaluate(scale, inters, _fbwave[if2],
						  vw2);
	    }
	    else {
	      SpinorBarWaveFunction inters = 
		_fer[idiag].first->evaluate(scale, widthOption(), offshell,
					    _fbwave[if2], vw1);
	      diag = _fer[idiag].second->evaluate(scale, _fwave[if1], inters, 
						  vw2);
	    }
	  }
	  else if( offshell->iSpin() == PDT::Spin0 ) {
	    ScalarWaveFunction inters = 
	      _sca[idiag].first->evaluate(scale, widthOption(), offshell, 
					  _fwave[if1], _fbwave[if2]);
	    diag = _sca[idiag].second->evaluate(scale, _vwave.first[iv1],
						_vwave.second[iv2], inters);
	  }
	  else if( offshell->iSpin() == PDT::Spin1 ) {
	    VectorWaveFunction interv = 
	      _vec[idiag].first->evaluate(scale, widthOption(), offshell, 
					  _fwave[if1], _fbwave[if2]);
	    diag = _vec[idiag].second->evaluate(scale, _vwave.first[iv1],
						_vwave.second[iv2], interv);
	  } 
	  else if( offshell->iSpin() == PDT::Spin2 ) {
	    TensorWaveFunction intert = 
	      _ten[idiag].first->evaluate(scale, widthOption(), offshell, 
					  _fwave[if1], _fbwave[if2]);
	    diag = _ten[idiag].second->evaluate(scale, _vwave.first[iv1],
						_vwave.second[iv2], intert);
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
	      (*mes[ix])(h1, h2, iv1, iv2) = flows[ix];
	      (*mel[ix])(h1, h2, iv1, iv2) = largeflows[ix];
	    }
	    else if(isp == 1) { 
	      (*mes[ix])(h1, iv1, h2, iv2) = flows[ix];
	      (*mel[ix])(h1, iv1, h2, iv2) = largeflows[ix];
	    }
	    else if(isp == 2) { 
	      (*mes[ix])(h1, iv1, iv2, h2) = flows[ix];
	      (*mel[ix])(h1, iv1, h2, iv2) = largeflows[ix];
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
	double con = cfactors[ix][iy]*(mes[ix]->contract(*mes[iy],_rho)).real();
	me2 += con;
	if(ix==iy) {
	  con = nfactors[ix][iy]*(mel[ix]->contract(*mel[iy],_rho)).real();
	  pflows[ix] += con;
	}
      }
    }
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
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow]->contract(*mel[iflow],_rho)).real();
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
		  outgoing()[2]->mass(),relativeError()));
}
