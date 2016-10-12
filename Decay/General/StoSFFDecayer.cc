// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StoSFFDecayer class.
//

#include "StoSFFDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

IBPtr StoSFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr StoSFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void StoSFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _sca << _fer << _vec << _ten;
}

void StoSFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _sca >> _fer >> _vec >> _ten;
}

ClassDescription<StoSFFDecayer> StoSFFDecayer::initStoSFFDecayer;
// Definition of the static class description member.

void StoSFFDecayer::Init() {

  static ClassDocumentation<StoSFFDecayer> documentation
    ("The StoSFFDecayer class implements the general decay of a scalar to "
     "a scalar and two fermions.");

}

WidthCalculatorBasePtr StoSFFDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<StoSFFDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),outgoing()[2]->mass(),
		  relativeError()));
}

void StoSFFDecayer::doinit() {
  GeneralThreeBodyDecayer::doinit(); 
  unsigned int ndiags = getProcessInfo().size();
  _sca.resize(ndiags);
  _fer.resize(ndiags);
  _vec.resize(ndiags);
  _ten.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if( offshell->CC() ) offshell = offshell->CC();
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractSSSVertexPtr vert1 = dynamic_ptr_cast<AbstractSSSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in StoSFFDecayer::doinit()"
	<< Exception::runerror;
      _sca[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1Half) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a fermion diagram in StoSFFDecayer::doinit()"
	<< Exception::runerror;
      _fer[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractVSSVertexPtr vert1 = dynamic_ptr_cast<AbstractVSSVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in StoSFFDecayer::doinit()"
	<< Exception::runerror;
      _vec[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractSSTVertexPtr vert1 = dynamic_ptr_cast<AbstractSSTVertexPtr>
	(current.vertices.first);
      AbstractFFTVertexPtr vert2 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a tensor diagram in StoSFFDecayer::doinit()"
	<< Exception::runerror;
      _ten[ix] = make_pair(vert1, vert2);
    }
  }
}

double StoSFFDecayer::me2(const int ichan, const Particle & inpart,
			  const ParticleVector & decay,
			  MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != inpart.id();
  // special handling or first/last call
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),
			     Helicity::incoming);
    _swave = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),
				Helicity::incoming);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
			Helicity::incoming,true);
    for(unsigned int ix=0;ix<decay.size();++ix) {
      if(decay[ix]->dataPtr()->iSpin()==PDT::Spin0) {
	ScalarWaveFunction::constructSpinInfo(decay[ix],Helicity::outgoing,true);
      }
      else {
	SpinorWaveFunction::
	  constructSpinInfo(_outspin[ix].first,decay[ix],Helicity::outgoing,true);
      }
    }
    return 0.;
  }
  // get the wavefunctions for all the particles
  ScalarWaveFunction outScalar;
  unsigned int isca(0);
  for(unsigned int ix=0;ix<decay.size();++ix) {
    if(decay[ix]->dataPtr()->iSpin()==PDT::Spin0) {
      isca = ix;
      outScalar = ScalarWaveFunction(decay[ix]->momentum(),
				     decay[ix]->dataPtr(),Helicity::outgoing);
    }
    else {
      SpinorWaveFunction::
	calculateWaveFunctions(_outspin[ix].first,decay[ix],Helicity::outgoing);
      _outspin[ix].second.resize(2);
      if(_outspin[ix].first[0].wave().Type() == SpinorType::u) {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  _outspin[ix].second[iy] = _outspin[ix].first[iy].bar();
	  _outspin[ix].first[iy].conjugate();
	}
      }
      else {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  _outspin[ix].second[iy] = _outspin[ix].first[iy].bar();
	  _outspin[ix].second[iy].conjugate();
	}
      }
    }
  }
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  Energy2 scale(sqr(inpart.mass()));  
  const size_t ncf(numberOfFlows());
  vector<Complex> flows(ncf, Complex(0.)), largeflows(ncf, Complex(0.)); 
  vector<GeneralDecayMEPtr> 
    mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
					      isca==0 ? PDT::Spin0 : PDT::Spin1Half,
					      isca==1 ? PDT::Spin0 : PDT::Spin1Half,
					      isca==2 ? PDT::Spin0 : PDT::Spin1Half)));
  vector<GeneralDecayMEPtr> 
    mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
					      isca == 0 ? PDT::Spin0 : PDT::Spin1Half,
					      isca == 1 ? PDT::Spin0 : PDT::Spin1Half,
					      isca == 2 ? PDT::Spin0 : PDT::Spin1Half)));
  static const unsigned int out2[3]={1,0,0},out3[3]={2,2,1};
  for(unsigned int s1 = 0;s1 < 2; ++s1) {
    for(unsigned int s2 = 0;s2 < 2; ++s2) {
      flows = vector<Complex>(ncf, Complex(0.));
      largeflows = vector<Complex>(ncf, Complex(0.));
      unsigned int idiag(0);
      for(vector<TBDiagram>::const_iterator dit = getProcessInfo().begin();
	  dit != getProcessInfo().end(); ++dit) {
	// channels if selecting
	if( ichan >= 0 && diagramMap()[ichan] != idiag ) {
	  ++idiag;
	  continue;
	}
	tcPDPtr offshell = dit->intermediate;
	if(cc&&offshell->CC()) offshell=offshell->CC();
	Complex diag;
	double sign = out3[dit->channelType] < out2[dit->channelType] ? 1. : -1.;
	// intermediate scalar
	if     (offshell->iSpin() == PDT::Spin0) {
	  ScalarWaveFunction inters = _sca[idiag].first->
	    evaluate(scale, widthOption(), offshell, _swave, outScalar);
	  unsigned int h1(s1),h2(s2);
	  if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	  if(decay[out2[dit->channelType]]->id()<0&&
	     decay[out3[dit->channelType]]->id()>0) {
	    diag =-sign*_sca[idiag].second->
	      evaluate(scale,
		       _outspin[out2[dit->channelType]].first [h1],
		       _outspin[out3[dit->channelType]].second[h2],inters);
	  }
	  else {
	    diag = sign*_sca[idiag].second->
	      evaluate(scale,
		       _outspin[out3[dit->channelType]].first [h2],
		       _outspin[out2[dit->channelType]].second[h1],inters);
	  }
	}
	// intermediate fermion
	else if(offshell->iSpin() == PDT::Spin1Half) {
	  int iferm = 
	    decay[out2[dit->channelType]]->dataPtr()->iSpin()==PDT::Spin1Half
	    ? out2[dit->channelType] : out3[dit->channelType];
	  unsigned int h1(s1),h2(s2);
	  if(dit->channelType>iferm) swap(h1,h2);
	  sign = iferm<dit->channelType ? 1. : -1.;


	  if((decay[dit->channelType]->id() < 0 &&decay[iferm]->id() > 0 ) ||
	     (decay[dit->channelType]->id()*offshell->id()>0)) {
	    SpinorWaveFunction    inters = _fer[idiag].first->
	      evaluate(scale,widthOption(),offshell,
		       _outspin[dit->channelType].first [h1],_swave);
	    diag = -sign*_fer[idiag].second->
	      evaluate(scale,inters,_outspin[iferm].second[h2],outScalar);
	  }
	  else {
	    SpinorBarWaveFunction inters = _fer[idiag].first->
	      evaluate(scale,widthOption(),offshell,
		       _outspin[dit->channelType].second[h1],_swave);
	    diag =  sign*_fer[idiag].second->
	      evaluate(scale,_outspin[iferm].first [h2],inters,outScalar);
	  }
	}
	// intermediate vector
	else if(offshell->iSpin() == PDT::Spin1) {
	  VectorWaveFunction interv = _vec[idiag].first->
	    evaluate(scale, widthOption(), offshell, _swave, outScalar);
	  unsigned int h1(s1),h2(s2);
	  if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	  if(decay[out2[dit->channelType]]->id()<0&&
	     decay[out3[dit->channelType]]->id()>0) {
	    diag =-sign*_vec[idiag].second->
	      evaluate(scale,
		       _outspin[out2[dit->channelType]].first [h1],
		       _outspin[out3[dit->channelType]].second[h2],interv);
	  }
	  else {
	    diag = sign*_vec[idiag].second->
	      evaluate(scale,
		       _outspin[out3[dit->channelType]].first [h2],
		       _outspin[out2[dit->channelType]].second[h1],interv);
	  }
	}
	// intermediate tensor
	else if(offshell->iSpin() == PDT::Spin2) {
 	  TensorWaveFunction intert = _ten[idiag].first->
	    evaluate(scale, widthOption(), offshell, _swave, outScalar);
	  unsigned int h1(s1),h2(s2);
	  if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	  if(decay[out2[dit->channelType]]->id()<0&&
	     decay[out3[dit->channelType]]->id()>0) {
	    diag =-sign*_ten[idiag].second->
	      evaluate(scale,
		       _outspin[out2[dit->channelType]].first [h1],
		       _outspin[out3[dit->channelType]].second[h2],intert);
	  }
	  else {
	    diag = sign*_ten[idiag].second->
	      evaluate(scale,
		       _outspin[out3[dit->channelType]].first [h2],
		       _outspin[out2[dit->channelType]].second[h1],intert);
	  }
	}
	// unknown
	else throw Exception()
	  << "Unknown intermediate in StoSFFDecayer::me2()" 
	  << Exception::runerror;
	// matrix element for the different colour flows
	if(ichan < 0) {
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
	    if(dit->colourFlow[iy].first - 1 != colourFlow()) continue;
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
      }
      for(unsigned int ix = 0; ix < ncf; ++ix) {
	if(isca == 0) {
	  (*mes[ix])(0, 0, s1, s2) = flows[ix];
	  (*mel[ix])(0, 0, s1, s2) = largeflows[ix];
	}
	else if(isca == 1 ) { 
	  (*mes[ix])(0, s1, 0, s2) = flows[ix];
	  (*mel[ix])(0, s1, 0, s2) = largeflows[ix];
	}
	else if(isca == 2) { 
	  (*mes[ix])(0, s1,s2, 0) = flows[ix];
	  (*mel[ix])(0, s1,s2, 0) = largeflows[ix] ;
	}
      }
    }
  }
  double me2(0.);
  if(ichan < 0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix]->contract(*mes[iy],_rho)).real();
	me2 += con;
	if(ix == iy) {
	  con = nfactors[ix][iy]*(mel[ix]->contract(*mel[iy],_rho)).real();
	  pflows[ix] += con;
	}
      }
    }
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *= UseRandom::rnd();
    for(unsigned int ix = 0;ix < pflows.size(); ++ix) {
      if(ptotal <= pflows[ix]) {
	colourFlow(ix);
	ME(mes[ix]);
	break;
      }
      ptotal -= pflows[ix];
    }
  }
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow]->contract(*mel[iflow],_rho)).real();
  }
  // return the matrix element squared
  return me2;
}
