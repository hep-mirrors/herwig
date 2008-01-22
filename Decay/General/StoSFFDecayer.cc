// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StoSFFDecayer class.
//

#include "StoSFFDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

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
		  outgoing()[0]->mass(),outgoing()[1]->mass(),outgoing()[2]->mass()));
}

void StoSFFDecayer::doinit() throw(InitException) {
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

double StoSFFDecayer::me2(bool vertex, const int ichan, const Particle & inpart,
			  const ParticleVector & decay) const {
  const vector<vector<double> > cfactors(getColourFactors());
  Energy2 scale(sqr(inpart.mass()));
  // spin density matrix
  RhoDMatrix rhoin(PDT::Spin0);
  rhoin.average();
  // get the wavefunctions for all the particles
  ScalarWaveFunction inScalar(const_ptr_cast<tPPtr>(&inpart),rhoin,
			      Helicity::incoming,true,vertex);
  ScalarWaveFunction outScalar;
  pair<vector<SpinorWaveFunction>,vector<SpinorBarWaveFunction> > outspin[3];
  unsigned int isca(0);
  for(unsigned int ix=0;ix<decay.size();++ix) {
    if(decay[ix]->dataPtr()->iSpin()==PDT::Spin0) {
      isca = ix;
      outScalar = ScalarWaveFunction(decay[ix],Helicity::outgoing,true,vertex);
    }
    else {
      SpinorWaveFunction(outspin[ix].first,decay[ix],Helicity::outgoing,true,vertex);
      if(outspin[ix].first[0].wave().Type() == u_spinortype) {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  outspin[ix].second.push_back(outspin[ix].first[iy].bar());
	  outspin[ix].first[iy].conjugate();
	}
      }
      else {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  outspin[ix].second.push_back(outspin[ix].first[iy].bar());
	  outspin[ix].second[iy].conjugate();
	}
      }
    }
  }
  const size_t ncf(numberOfFlows());
  vector<Complex> flows(ncf, Complex(0.));
  vector<DecayMatrixElement> 
    mes(ncf,DecayMatrixElement(PDT::Spin0,
			       isca==0 ? PDT::Spin0 : PDT::Spin1Half,
			       isca==1 ? PDT::Spin0 : PDT::Spin1Half,
			       isca==2 ? PDT::Spin0 : PDT::Spin1Half));
  static const unsigned int out2[3]={1,0,0},out3[3]={2,2,1};
  for(unsigned int s1=0;s1<2;++s1) {
    for(unsigned int s2=0;s2<2;++s2) {
      flows = vector<Complex>(ncf, Complex(0.));
      int nchan=-1;
      unsigned int idiag=0;
      for(vector<TBDiagram>::const_iterator dit=getProcessInfo().begin();
	  dit!=getProcessInfo().end();++dit) {
	// channels if selecting
	++nchan;
	if(ichan>0&&ichan!=nchan) continue;
	tcPDPtr offshell = dit->intermediate;
	Complex diag;
	double sign = out3[dit->channelType]<out2[dit->channelType] ? 1. : -1.;
	// intermediate scalar
	if     (offshell->iSpin() == PDT::Spin0) {
	  ScalarWaveFunction inters = _sca[idiag].first->
	    evaluate(scale, widthOption(), offshell, inScalar, outScalar);
	  unsigned int h1(s1),h2(s2);
	  if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	  if(decay[out2[dit->channelType]]->id()<0&&
	     decay[out3[dit->channelType]]->id()>0) {
	    diag =-sign*_sca[idiag].second->
	      evaluate(scale,
		       outspin[out2[dit->channelType]].first [h1],
		       outspin[out3[dit->channelType]].second[h2],inters);
	  }
	  else {
	    diag = sign*_sca[idiag].second->
	      evaluate(scale,
		       outspin[out3[dit->channelType]].first [h2],
		       outspin[out2[dit->channelType]].second[h1],inters);
	  }
	}
	// intermediate fermion
	else if(offshell->iSpin() == PDT::Spin1Half) {
	  int iferm = decay[out2[dit->channelType]]->dataPtr()->iSpin()==PDT::Spin1Half
	    ? out2[dit->channelType] : out3[dit->channelType];
	  unsigned int h1(s1),h2(s2);
	  if(dit->channelType>iferm) swap(h1,h2);
	  sign = iferm<dit->channelType ? 1. : -1.;
	  if(decay[dit->channelType]->id()<0&&decay[iferm]->id()>0) {
	    SpinorWaveFunction    inters = _fer[idiag].first->
	      evaluate(scale,widthOption(),offshell,
		       outspin[dit->channelType].first [h1],inScalar);
	    diag = -sign*_fer[idiag].second->
	      evaluate(scale,inters,outspin[iferm].second[h2],outScalar);
	  }
	  else {
	    SpinorBarWaveFunction inters = _fer[idiag].first->
	      evaluate(scale,widthOption(),offshell,
		       outspin[dit->channelType].second[h1],inScalar);
	    diag =  sign*_fer[idiag].second->
	      evaluate(scale,outspin[iferm].first [h2],inters,outScalar);
	  }
	}
	// intermediate vector
	else if(offshell->iSpin() == PDT::Spin1) {
	  VectorWaveFunction interv = _vec[idiag].first->
	    evaluate(scale, widthOption(), offshell, inScalar, outScalar);
	  unsigned int h1(s1),h2(s2);
	  if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	  if(decay[out2[dit->channelType]]->id()<0&&
	     decay[out3[dit->channelType]]->id()>0) {
	    diag =-sign*_vec[idiag].second->
	      evaluate(scale,
		       outspin[out2[dit->channelType]].first [h1],
		       outspin[out3[dit->channelType]].second[h2],interv);
	  }
	  else {
	    diag = sign*_vec[idiag].second->
	      evaluate(scale,
		       outspin[out3[dit->channelType]].first [h2],
		       outspin[out2[dit->channelType]].second[h1],interv);
	  }
	}
	// intermediate tensor
	else if(offshell->iSpin() == PDT::Spin2) {
 	  TensorWaveFunction intert = _ten[idiag].first->
	    evaluate(scale, widthOption(), offshell, inScalar, outScalar);
	  unsigned int h1(s1),h2(s2);
	  if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	  if(decay[out2[dit->channelType]]->id()<0&&
	     decay[out3[dit->channelType]]->id()>0) {
	    diag =-sign*_ten[idiag].second->
	      evaluate(scale,
		       outspin[out2[dit->channelType]].first [h1],
		       outspin[out3[dit->channelType]].second[h2],intert);
	  }
	  else {
	    diag = sign*_ten[idiag].second->
	      evaluate(scale,
		       outspin[out3[dit->channelType]].first [h2],
		       outspin[out2[dit->channelType]].second[h1],intert);
	  }
	}
	// unknown
	else throw Exception()
	  << "Unknown intermediate in FtoFFFDecayer::me2()" 
	  << Exception::runerror;
	// matrix element for the different colour flows
	for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
	  flows[dit->colourFlow[iy].first - 1] += 
	    dit->colourFlow[iy].second * diag;
	}
	++idiag;
      }
      // now add the flows to the me2 with appropriate colour factors
      for(unsigned int ix = 0; ix < ncf; ++ix) {
	if     (isca==0) mes[ix](0, 0,s1,s2) = flows[ix];
	else if(isca==1) mes[ix](0,s1, 0,s2) = flows[ix];
	else if(isca==2) mes[ix](0,s1,s2, 0) = flows[ix];
      }
    }
  }
  vector<double> pflows(ncf,0.);
  double me2(0.);
  for(unsigned int ix = 0; ix < ncf; ++ix) {
    for(unsigned int iy = 0; iy < ncf; ++ iy) {
      double con = cfactors[ix][iy]*(mes[ix].contract(mes[iy],rhoin)).real();
      if(ix==iy) pflows[ix] += con;
      me2 += con;
    }
  }
  // select the matrix element according to the colour flow
  if(vertex) {
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *=UseRandom::rnd();
    for(unsigned int ix=0;ix<pflows.size();++ix) {
      if(ptotal<=pflows[ix]) {
	ME(mes[ix]);
	break;
      }
      ptotal-=pflows[ix];
    }
  }
  // make the colour connections
  colourConnections(inpart, decay);
  // return the matrix element squared
  return me2;
}
