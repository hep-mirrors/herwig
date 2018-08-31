// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FtoFFFDecayer class.
//

#include "FtoFFFDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;

IBPtr FtoFFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FtoFFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void FtoFFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << sca_ << vec_ << ten_;
}

void FtoFFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> sca_ >> vec_ >> ten_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FtoFFFDecayer,GeneralThreeBodyDecayer>
describeHerwigFtoFFFDecayer("Herwig::FtoFFFDecayer", "Herwig.so");

void FtoFFFDecayer::Init() {

  static ClassDocumentation<FtoFFFDecayer> documentation
    ("The FtoFFFDecayer class implements the general decay of a fermion to "
     "three fermions.");

}

void FtoFFFDecayer::setupDiagrams(bool kinCheck) {
  GeneralThreeBodyDecayer::setupDiagrams(kinCheck);
  if(outgoing().empty()) return;
  unsigned int ndiags = getProcessInfo().size();
  sca_.resize(ndiags);
  vec_.resize(ndiags);
  ten_.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if( offshell->CC() ) offshell = offshell->CC();
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in FtoFFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      sca_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in FtoFFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      vec_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractFFTVertexPtr vert1 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.first);
      AbstractFFTVertexPtr vert2 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a tensor diagram in FtoFFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      ten_[ix] = make_pair(vert1, vert2);
    }
  }
}

double  FtoFFFDecayer::me2(const int ichan, const Particle & inpart,
			   const ParticleVector & decay,
			   MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != inpart.id();
  // special handling or first/last call
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  const size_t ncf(numberOfFlows());
  Energy2 scale(sqr(inpart.mass()));
  if(meopt==Initialize) {
    SpinorWaveFunction::
      calculateWaveFunctions(inwave_.first,rho_,const_ptr_cast<tPPtr>(&inpart),
			    Helicity::incoming);
    inwave_.second.resize(2);
    if(inwave_.first[0].wave().Type() == SpinorType::u) {
      for(unsigned int ix = 0; ix < 2; ++ix) {
	inwave_.second[ix] = inwave_.first[ix].bar();
	inwave_.second[ix].conjugate();
      }
    }
    else {
      for(unsigned int ix = 0; ix < 2; ++ix) {
	inwave_.second[ix] = inwave_.first[ix].bar();
	inwave_.first[ix].conjugate();
      }
    }
    // fix rho if no correlations
    fixRho(rho_);
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()<0) 
      SpinorWaveFunction::constructSpinInfo(inwave_.first,
					    const_ptr_cast<tPPtr>(&inpart),
					    Helicity::incoming,true);
    else
      SpinorBarWaveFunction::constructSpinInfo(inwave_.second,
					       const_ptr_cast<tPPtr>(&inpart),
					       Helicity::incoming,true);
    // outgoing particles
    for(unsigned int ix = 0; ix < 3; ++ix) {
      SpinorWaveFunction::
      constructSpinInfo(outwave_[ix].first,decay[ix],Helicity::outgoing,true);
    }
  }
  // outgoing particles
  for(unsigned int ix = 0; ix < 3; ++ix) {
    SpinorWaveFunction::
      calculateWaveFunctions(outwave_[ix].first,decay[ix],Helicity::outgoing);
    outwave_[ix].second.resize(2);
    if(outwave_[ix].first[0].wave().Type() == SpinorType::u) {
      for(unsigned int iy = 0; iy < 2; ++iy) {
	outwave_[ix].second[iy] = outwave_[ix].first[iy].bar();
	outwave_[ix].first[iy].conjugate();
      }
    }
    else {
      for(unsigned int iy = 0; iy < 2; ++iy) {
	outwave_[ix].second[iy] = outwave_[ix].first[iy].bar();
	outwave_[ix].second[iy].conjugate();
      }
    }
  }
  bool ferm = inpart.id()>0; 
  vector<Complex> flows(ncf, Complex(0.)),largeflows(ncf, Complex(0.)); 
  static const unsigned int out2[3]={1,0,0},out3[3]={2,2,1};
  vector<GeneralDecayMEPtr> mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
								      PDT::Spin1Half,PDT::Spin1Half)));
  vector<GeneralDecayMEPtr> mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
								      PDT::Spin1Half,PDT::Spin1Half)));
  unsigned int ihel[4];
  for(ihel[0] = 0; ihel[0] < 2; ++ihel[0]) {
    for(ihel[1] = 0; ihel[1] < 2; ++ihel[1]) {
      for(ihel[2] = 0; ihel[2] < 2; ++ihel[2]) {
	for(ihel[3] = 0; ihel[3] < 2; ++ihel[3]) {
	  flows = vector<Complex>(ncf, Complex(0.));
	  largeflows = vector<Complex>(ncf, Complex(0.));
	  unsigned int idiag=0;
	  for(vector<TBDiagram>::const_iterator dit=getProcessInfo().begin();
	      dit!=getProcessInfo().end();++dit) {
	    if(ichan>=0&&diagramMap()[ichan]!=idiag) {
	      ++idiag;
	      continue;
	    }
	    // the sign from normal ordering
	    double sign = ferm ? 1. : -1;
	    // outgoing wavefunction and NO sign
	    if     (dit->channelType==TBDiagram::channel23) sign *= -1.;
	    else if(dit->channelType==TBDiagram::channel13) sign *=  1.;
	    else if(dit->channelType==TBDiagram::channel12) sign *= -1.;
	    else throw Exception()
	      << "Unknown diagram type in FtoFFFDecayer::me2()" << Exception::runerror;
	    // wavefunctions
	    SpinorWaveFunction    w0,w3;
	    SpinorBarWaveFunction w1,w2;
	    // incoming wavefunction
	    if(ferm) {
	      w0 = inwave_.first [ihel[0]];
	      w1 = outwave_[dit->channelType].second[ihel[dit->channelType+1]];
	    }
	    else {
	      w0 = outwave_[dit->channelType].first [ihel[dit->channelType+1]];
	      w1 = inwave_.second[ihel[0]];
	    }
	    if(decay[out2[dit->channelType]]->id()<0&&
	       decay[out3[dit->channelType]]->id()>0) {
	      w2 = outwave_[out3[dit->channelType]].second[ihel[out3[dit->channelType]+1]];
	      w3 = outwave_[out2[dit->channelType]].first [ihel[out2[dit->channelType]+1]];
	      sign *= -1.;
	    }
	    else {
	      w2 = outwave_[out2[dit->channelType]].second[ihel[out2[dit->channelType]+1]];
	      w3 = outwave_[out3[dit->channelType]].first [ihel[out3[dit->channelType]+1]];
	    }
	    tcPDPtr offshell = dit->intermediate;
	    if(cc&&offshell->CC()) offshell=offshell->CC();
	    Complex diag(0.);
	    // intermediate scalar
	    if     (offshell->iSpin() == PDT::Spin0) { 
	      ScalarWaveFunction inters = sca_[idiag].first->
		evaluate(scale, widthOption(), offshell, w0, w1);
	      diag = sca_[idiag].second->evaluate(scale,w3,w2,inters);
	    }
	    // intermediate vector
	    else if(offshell->iSpin() == PDT::Spin1) {
	      VectorWaveFunction interv = vec_[idiag].first->
		evaluate(scale, widthOption(), offshell, w0, w1);
	      diag = vec_[idiag].second->evaluate(scale,w3,w2,interv);
	    }
	    // intermediate tensor
	    else if(offshell->iSpin() == PDT::Spin2) {
	      TensorWaveFunction intert = ten_[idiag].first->
		evaluate(scale, widthOption(), offshell, w0, w1);
	      diag = ten_[idiag].second->evaluate(scale,w3,w2,intert);
	    }
	    // unknown
	    else throw Exception()
	      << "Unknown intermediate in FtoFFFDecayer::me2()" 
	      << Exception::runerror;
	    // apply NO sign
	    diag *= sign;
	    // matrix element for the different colour flows
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
	  }
	  // now add the flows to the me2 with appropriate colour factors
	  for(unsigned int ix = 0; ix < ncf; ++ix) {
	    (*mes[ix])(ihel[0],ihel[1],ihel[2],ihel[3]) =      flows[ix];
	    (*mel[ix])(ihel[0],ihel[1],ihel[2],ihel[3]) = largeflows[ix];
	  }
	}
      }
    }
  }
  double me2(0.);
  if(ichan<0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix]->contract(*mes[iy],rho_)).real();
	me2 += con;
	if(ix==iy) {
	  con = nfactors[ix][iy]*(mel[ix]->contract(*mel[iy],rho_)).real();
	  pflows[ix] += con;
	}
      }
    }
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *=UseRandom::rnd();
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
    me2 = nfactors[iflow][iflow]*(mel[iflow]->contract(*mel[iflow],rho_)).real();
  }
  // return the matrix element squared
  return me2;
}

WidthCalculatorBasePtr FtoFFFDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<FtoFFFDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),outgoing()[2]->mass(),
		 relativeError()));
}
