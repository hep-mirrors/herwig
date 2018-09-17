// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StoFFFFDecayer class.
//

#include "StoFFFFDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG::Helicity;


namespace {

  inline bool isParticle(tPPtr part) {
    return part->id()>0 && part->dataPtr()->CC();
  }
  inline bool isAntiParticle(tPPtr part) {
    return part->id()<0 && part->dataPtr()->CC();
  }
  inline bool isMajorana(tPPtr part) {
    return !part->dataPtr()->CC();
  }

}

IBPtr StoFFFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr StoFFFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void StoFFFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << firstVVS_ << firstVSS_ << firstSSS_ << firstFFS_
     << secondFFV_ << secondFFS_
     << thirdFFV_ << thirdFFS_ << sign_;
}

void StoFFFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> firstVVS_ >> firstVSS_ >> firstSSS_ >> firstFFS_
     >> secondFFV_ >> secondFFS_
     >> thirdFFV_ >> thirdFFS_ >> sign_;
}

DescribeClass<StoFFFFDecayer,GeneralFourBodyDecayer>
describeStoFFFFDecayer("Herwig::StoFFFFDecayer", "Herwig.so");

void StoFFFFDecayer::Init() {

  static ClassDocumentation<StoFFFFDecayer> documentation
    ("The StoFFFFDecayer class performs the 4-fermion "
     "decays of scalar particles in BSM models");

}

void StoFFFFDecayer::doinit() {
  GeneralFourBodyDecayer::doinit(); 
  unsigned int ndiags = getProcessInfo().size();
  firstVVS_ .resize(ndiags);
  firstVSS_ .resize(ndiags);
  firstSSS_ .resize(ndiags);
  firstFFS_ .resize(ndiags);
  secondFFV_.resize(ndiags);
  secondFFS_.resize(ndiags);  
  thirdFFV_ .resize(ndiags);
  thirdFFS_ .resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    const NBDiagram & current = getProcessInfo()[ix];
    // first vertex
    firstVVS_[ix] = dynamic_ptr_cast<AbstractVVSVertexPtr>(current.vertex);
    firstVSS_[ix] = dynamic_ptr_cast<AbstractVSSVertexPtr>(current.vertex);
    firstSSS_[ix] = dynamic_ptr_cast<AbstractSSSVertexPtr>(current.vertex);
    firstFFS_[ix] = dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertex);
    // get the other vertices
    const NBVertex & second = current.vertices.begin()->second.incoming ?
      current.vertices.begin()->second : (++current.vertices.begin())->second;
    // get the other vertices
    const NBVertex & third = current.vertices.begin()->second.incoming ?
      (++current.vertices.begin())->second : (++second.vertices.begin())->second;
    // second vertex
    secondFFV_[ix] = dynamic_ptr_cast<AbstractFFVVertexPtr>(second.vertex);
    secondFFS_[ix] = dynamic_ptr_cast<AbstractFFSVertexPtr>(second.vertex);  
    // third vertex
    thirdFFV_ [ix] = dynamic_ptr_cast<AbstractFFVVertexPtr>(third .vertex);
    thirdFFS_ [ix] = dynamic_ptr_cast<AbstractFFSVertexPtr>(third .vertex); 
    assert( ( firstVVS_[ix] ||  firstVSS_[ix] || 
	      firstSSS_[ix] ||  firstFFS_[ix] ) && 
	    (secondFFV_[ix] || secondFFS_[ix] ) &&
	    ( thirdFFV_[ix] ||  thirdFFS_[ix] ));
    // NO sign
    int order = 
      current.channelType[0]*1000+current.channelType[1]*100+
      current.channelType[2]*10  +current.channelType[3];
    switch(order) {
    case 1234: case 1342: case 1423: 
    case 2143: case 2314: case 2431:
    case 3124: case 3241: case 3412:
    case 4132: case 4213: case 4321:
      sign_.push_back( 1.);
      break;
    case 1243: case 1324: case 1432:
    case 2134: case 2341: case 2413:
    case 3142: case 3214: case 3421:
    case 4123: case 4231: case 4312:
      sign_.push_back(-1.);
      break;
    default:
      assert(false);
    }
  }
}

double StoFFFFDecayer::me2(const int ichan, const Particle & inpart,
			   const ParticleVector & decay, MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming->id() != inpart.id();
  // special handling or first/last call
  const vector<vector<double> > & cfactors(getColourFactors());
  const vector<vector<double> > & nfactors(getLargeNcColourFactors());
  const size_t ncf(numberOfFlows());
  Energy2 scale(sqr(inpart.mass()));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&inpart),
			     Helicity::incoming);
    swave_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),
				Helicity::incoming);
    // fix rho if no correlations
    fixRho(rho_);
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
			Helicity::incoming,true);
    // outgoing particles
    for(unsigned int ix = 0; ix < 4; ++ix) {
      SpinorWaveFunction::
      constructSpinInfo(outwave_[ix].first,decay[ix],Helicity::outgoing,true);
    }
  }
  // outgoing particles
  for(unsigned int ix = 0; ix < 4; ++ix) {
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
  // matrix element for the colour flows
  vector<Complex> flows(ncf, Complex(0.)),largeflows(ncf, Complex(0.));
  vector<GeneralDecayMEPtr> mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
								      PDT::Spin1Half,PDT::Spin1Half,
								      PDT::Spin1Half,PDT::Spin1Half)));
  vector<GeneralDecayMEPtr> mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
								      PDT::Spin1Half,PDT::Spin1Half,
								      PDT::Spin1Half,PDT::Spin1Half)));
  // calculate the matrix element
  unsigned int ihel[4];
  for(ihel[0] = 0; ihel[0] < 2; ++ihel[0]) {
    for(ihel[1] = 0; ihel[1] < 2; ++ihel[1]) {
      for(ihel[2] = 0; ihel[2] < 2; ++ihel[2]) {
	for(ihel[3] = 0; ihel[3] < 2; ++ihel[3]) {
	  flows = vector<Complex>(ncf, Complex(0.));
	  largeflows = vector<Complex>(ncf, Complex(0.));
	  unsigned int idiag=0;
	  for(vector<NBDiagram>::const_iterator dit=getProcessInfo().begin();
	      dit!=getProcessInfo().end();++dit) {
	    if(ichan>=0&&diagramMap()[ichan]!=idiag) {
	      ++idiag;
	      continue;
	    }
	    // location of the particles
	    int iloc[4];
	    for(unsigned int ix=0;ix<4;++ix) iloc[ix] = dit->channelType[ix]-1;
	    // NO sign
	    double sign = sign_[idiag];
	    // work out the type of topology
	    bool topo = dit->vertices.begin()->second.incoming;
	    const NBVertex & second = topo ? 
   	          dit->vertices.begin()  ->second :
	      (++(dit->vertices.begin()))->second;
	    const NBVertex & third = topo ?
	      (++(dit->  vertices.begin()))->second :
	      (++(second.vertices.begin()))->second;
	    // extract the intermediates
	    tPDPair inter = make_pair(second.incoming,
				      third .incoming);
	    if(    inter.second->CC()) inter.second = inter.second->CC();
	    if(cc&&inter.first ->CC()) inter.first  = inter.first ->CC();
	    if(cc&&inter.second->CC()) inter.second = inter.second->CC();
	    // value for the diagram
	    Complex diag(0.);
	    // first compute the last part of the diagram
	    VectorWaveFunction offVector2;
	    ScalarWaveFunction offScalar2;
	    // intermediate scalar
	    if(inter.second->iSpin()==PDT::Spin0) {
	      if( (isAntiParticle(decay[iloc[2]]) || isMajorana(decay[iloc[2]])) &&
		   (isParticle   (decay[iloc[3]]) || isMajorana(decay[iloc[3]])) ) {
		sign *= -1.;
		offScalar2 = thirdFFS_[idiag]->
		  evaluate(scale, widthOption(),inter.second,
			   outwave_[iloc[2]].first [ihel[iloc[2]]],
			   outwave_[iloc[3]].second[ihel[iloc[3]]]);
	      }
	      else {
		offScalar2 = thirdFFS_[idiag]->
		  evaluate(scale, widthOption(),inter.second,
			   outwave_[iloc[3]].first [ihel[iloc[3]]],
			   outwave_[iloc[2]].second[ihel[iloc[2]]]);
	      }
	    }
	    // intermediate vector
	    else if(inter.second->iSpin()==PDT::Spin1) {
	      if( (isAntiParticle(decay[iloc[2]]) || isMajorana(decay[iloc[2]])) &&
		   (isParticle(decay[iloc[3]])||isMajorana(decay[iloc[3]]))) {
		sign *= -1.;
		offVector2 = thirdFFV_[idiag]->
		  evaluate(scale, widthOption(),inter.second,
			   outwave_[iloc[2]].first [ihel[iloc[2]]],
			   outwave_[iloc[3]].second[ihel[iloc[3]]]);
	      }
	      else {
		offVector2 = thirdFFV_[idiag]->
		  evaluate(scale, widthOption(),inter.second,
			   outwave_[iloc[3]].first [ihel[iloc[3]]],
			   outwave_[iloc[2]].second[ihel[iloc[2]]]);
	      }
	    }
	    // unknown
	    else
	      assert(false);
	    // first topology
	    if(topo) {
	      // first intermediate
	      if(inter.first->CC()) inter.first = inter.first->CC();
	      VectorWaveFunction offVector1;
	      ScalarWaveFunction offScalar1;
	      // intermeidate scalar
	      if(inter.first->iSpin()==PDT::Spin0) {
		if(decay[iloc[0]]->id()<0&&
		   decay[iloc[1]]->id()>0) {
		  sign *= -1.;
		  offScalar1 = secondFFS_[idiag]->
		    evaluate(scale, widthOption(),inter.first,
			     outwave_[iloc[0]].first [ihel[iloc[0]]],
			     outwave_[iloc[1]].second[ihel[iloc[1]]]);
		}
		else {
		  offScalar1 = secondFFS_[idiag]->
		    evaluate(scale, widthOption(),inter.first,
			     outwave_[iloc[1]].first [ihel[iloc[1]]],
			     outwave_[iloc[0]].second[ihel[iloc[0]]]);
		}
	      }
	      // intermediate vector
	      else if(inter.first->iSpin()==PDT::Spin1) {
		if(decay[iloc[0]]->id()<0&&
		   decay[iloc[1]]->id()>0) {
		  sign *= -1.;
		  offVector1 = secondFFV_[idiag]->
		    evaluate(scale, widthOption(),inter.first,
			     outwave_[iloc[0]].first [ihel[iloc[0]]],
			     outwave_[iloc[1]].second[ihel[iloc[1]]]);
		}
		else {
		  offVector1 = secondFFV_[idiag]->
		    evaluate(scale, widthOption(),inter.first,
			     outwave_[iloc[1]].first [ihel[iloc[1]]],
			     outwave_[iloc[0]].second[ihel[iloc[0]]]);
		}
	      }
	      // unknown
	      else
		assert(false);
	      // put it all together	      
	      if(inter.first->iSpin()==PDT::Spin0) {	      
		if(inter.second->iSpin()==PDT::Spin0) {
		  diag = firstSSS_[idiag]->
		    evaluate(scale,swave_,offScalar1,offScalar2);
		}
		else if(inter.second->iSpin()==PDT::Spin1) {
		  diag = firstVSS_[idiag]->
		    evaluate(scale,offVector2,offScalar1,swave_);
		}
	      }
	      else if(inter.first->iSpin()==PDT::Spin1) {	      
		if(inter.second->iSpin()==PDT::Spin0) {
		  diag = firstVSS_[idiag]->
		    evaluate(scale,offVector1,offScalar2,swave_);
		}
		else if(inter.second->iSpin()==PDT::Spin1) {
		  diag = firstVVS_[idiag]->
		    evaluate(scale,offVector1,offVector2,swave_);
		}
	      }
	    }
	    // second topology
	    else {
	      if(((isAntiParticle(decay[iloc[0]]) || isMajorana(decay[iloc[0]]))&&
		  (isParticle    (decay[iloc[1]]) || isMajorana(decay[iloc[1]])))) {
		sign *= -1.;
		SpinorWaveFunction    inters = firstFFS_[idiag]->
		  evaluate(scale,widthOption(),inter.first,
			   outwave_[iloc[0]].first [ihel[iloc[0]]],swave_);
		if(inter.second->iSpin()==PDT::Spin0) {
		  diag = secondFFS_[idiag]->
		    evaluate(scale,inters,outwave_[iloc[1]].second[ihel[iloc[1]]],
			     offScalar2);
		}
		else {
		  diag = secondFFV_[idiag]->
		    evaluate(scale,inters,outwave_[iloc[1]].second[ihel[iloc[1]]],
			     offVector2);
		}
	      }
	      else {
		SpinorBarWaveFunction inters = firstFFS_[idiag]->
		  evaluate(scale,widthOption(),inter.first,
			   outwave_[iloc[0]].second[ihel[iloc[0]]],swave_);
		if(inter.second->iSpin()==PDT::Spin0) {
		  diag = secondFFS_[idiag]->
		    evaluate(scale,outwave_[iloc[1]].first [ihel[iloc[1]]],inters,
			     offScalar2);
		}
		else {
		  diag = secondFFV_[idiag]->
		    evaluate(scale,outwave_[iloc[1]].first [ihel[iloc[1]]],inters,
			     offVector2);
		}
	      }
	    }
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
	    (*mes[ix])(0,ihel[0],ihel[1],ihel[2],ihel[3]) =      flows[ix];
	    (*mel[ix])(0,ihel[0],ihel[1],ihel[2],ihel[3]) = largeflows[ix];
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
  return me2*scale*UnitRemoval::InvE2;
}

// OLD TESTING CODE

// extracted from StandardModel.h
// public:

//   virtual void StopCouplings(Energy2 &,tcPDPtr, tcPDPtr,
// 			     double &, double &, double &,
// 			     Complex &, Complex &,
// 			     vector<Complex> &, vector<Complex> &,
// 			     vector<Complex> &, 
// 			     vector<Complex> &, vector<Complex> &, 
// 			     vector<Complex> &, vector<Complex> &, 
// 			     vector<Complex> &, vector<Complex> &, 
// 			     vector<vector<Complex> > &, vector<vector<Complex> > &, 
// 			     vector<vector<Complex> > &, vector<vector<Complex> > &, 
// 			     vector<Complex> &, vector<Complex> &, 
// 			     vector<Complex> &, vector<Complex> &,
// 			     double &, double &) {
//     assert(false);
//   }

// extracted from RunningMass.cc
  // if(id==5) return 4.8787783899999999*GeV;
  // if(id==15) return 1.7770999999999999*GeV;

// extracted from MSSM.h

// public:

//   virtual void StopCouplings(Energy2 &, tcPDPtr, tcPDPtr,
// 			     double & g, double & sw, double & cw,
// 			     Complex & aL, Complex & aR,
// 			     vector<Complex> & cL, vector<Complex> & cR,
// 			     vector<Complex> & d,
// 			     vector<Complex> & bL, vector<Complex> & bR,
// 			     vector<Complex> & kL, vector<Complex> & kR,
// 			     vector<Complex> & fL, vector<Complex> & fR,
// 			     vector<vector<Complex> > & eL, vector<vector<Complex> > & eR,
// 			     vector<vector<Complex> > & gL, vector<vector<Complex> > & gR,
// 			     vector<Complex> & hL, vector<Complex> & hR,
// 			     vector<Complex> & lL, vector<Complex> & lR,
// 			     double & ytop, double & ytau);

// extracted from MSSM.cc

// void MSSM::StopCouplings(Energy2 & scale, tcPDPtr ferm, tcPDPtr anti, double & g, double & sw, double & cw,
// 			 Complex & aL, Complex & aR,
// 			 vector<Complex> & cL, vector<Complex> & cR,
// 			 vector<Complex> & d,
// 			 vector<Complex> & bL, vector<Complex> & bR,
// 			     vector<Complex> & kL, vector<Complex> & kR,
// 			 vector<Complex> & fL, vector<Complex> & fR,
// 			 vector<vector<Complex> > & eL, vector<vector<Complex> > & eR,
// 			 vector<vector<Complex> > & gL, vector<vector<Complex> > & gR,
// 			 vector<Complex> & hL, vector<Complex> & hR,
// 			 vector<Complex> & lL, vector<Complex> & lR,
// 			 double & ytop, double & ytau) {
//   MixingMatrixPtr stop = stopMix();
//   MixingMatrixPtr sbot = sbottomMix();
//   MixingMatrixPtr stau = stauMix();
//   MixingMatrixPtr neut = neutralinoMix();
//   MixingMatrixPtr vmix = charginoVMix();
//   MixingMatrixPtr umix = charginoUMix();
//   sw = sqrt(   sin2ThetaW());
//   cw = sqrt(1.-sin2ThetaW());
//   g  =  sqrt(4.0*Constants::pi*alphaEMMZ()/sin2ThetaW());
//   Energy mb = mass(scale,getParticleData(ParticleID::b));
//   Energy mt = mass(scale,getParticleData(ParticleID::t));
//   Energy mw = getParticleData(ParticleID::Wplus)->mass();
//   double tb = tanBeta();
//   double sb = tb/sqrt(1 + sqr(tb));
//   double cb = sqrt(1 - sqr(sb));
//   Complex n1prime = (*neut)(0,0)*cw + (*neut)(0,1)*sw;
//   Complex n2prime = (*neut)(0,1)*cw - (*neut)(0,0)*sw;
//   double yt =  double(mt/mw)/sb/sqrt(2.);
//   Complex bracketl = 2./ 3.*sw*( conj(n1prime) - sw*conj(n2prime)/cw );
//   Complex bracketr = 2./3.*sw*n1prime - n2prime*(-0.5 + 2./3.*sqr(sw))/cw;
//   ytop = mt/tb/mw;
//   aL = -yt*conj((*neut)(0,3))*(*stop)(0,0) + sqrt(2.)*(*stop)(0,1)*bracketl;
//   aR = -yt*     (*neut)(0,3) *(*stop)(0,1) - sqrt(2.)*(*stop)(0,0)*bracketr;
//   cL.resize(2); cR.resize(2); d.resize(2.);
//   kL.resize(2); kR.resize(2);
//   bL.resize(2); bR.resize(2);
//   double yb =  double(mb/mw)/cb/sqrt(2.);
//   bracketl = -1./3.*sw*( conj(n1prime) - sw*conj(n2prime)/cw );
//   bracketr = -1./3.*sw*n1prime - n2prime*(0.5  -1./3.*sqr(sw))/cw;
//   for(unsigned int k=0;k<2;++k) {
//     cL[k] =-yb*conj((*neut)(0,2))*(*sbot)(k,0) + sqrt(2.)*(*sbot)(k,1)*bracketl;
//     cR[k] =-yb*     (*neut)(0,2) *(*sbot)(k,1) - sqrt(2.)*(*sbot)(k,0)*bracketr;
//     d[k] = (*stop)(0,0)*(*sbot)(k,0);
//     bL[k] = mb*conj((*umix)(k,1))/sqrt(2.)/mw/cb*(*stop)(0,0);
//     bR[k] = -(*vmix)(k,0)*(*stop)(0,0)+mt*(*vmix)(k,1)/sqrt(2.)/mw/sb*(*stop)(0,1);
//     kR[k] = -     (*neut)(0,3)*conj((*vmix)(k,1))/sqrt(2.)
//       +     (*neut)(0,1) *conj((*vmix)(k,0));
//     kL[k] =  conj((*neut)(0,2))*    (*umix)(k,1) /sqrt(2.)
//       +conj((*neut)(0,1))*     (*umix)(k,0) ;
//   }
//   fL.resize(2); fR.resize(2);
//   double qf = ferm->charge()/eplus;
//   Energy mf = (abs(ferm->id())<5||(abs(ferm->id())>=11&&abs(ferm->id())<=14)) ? ZERO : mass(scale, ferm);
//   Energy ma = (abs(anti->id())<5||(abs(anti->id())>=11&&abs(anti->id())<=14)) ? ZERO : mass(scale, anti);
//   fL[0] = 0.;
//   fR[0] = - sqrt(2.)*(qf*sw*n1prime - n2prime*(-0.5 + qf*sqr(sw))/cw);
//   fL[1] = + sqrt(2.)*qf*sw*( conj(n1prime) - sw*conj(n2prime)/cw );
//   fR[1] = 0.;
//   eL.resize(2,vector<Complex>(2,0.));
//   eR.resize(2,vector<Complex>(2,0.));
//   for(unsigned int i=0;i<2;++i) {
//     eR[i][0] = ma*conj((*umix)(i,1))/sqrt(2.)/mw/cb;
//     eL[i][0] = -(*vmix)(i,0); 
//     eL[i][1] = 0.;
//     eR[i][1] = 0.;
//   }
//   hL.resize(2); hR.resize(2);
//   double ya =  double(ma/mw)/cb/sqrt(2.);
//   qf =-anti->charge()/eplus;
//   bracketl = qf*sw*( conj(n1prime) - sw*conj(n2prime)/cw );
//   bracketr = qf*sw*n1prime - n2prime*(0.5 +qf*sqr(sw))/cw;
//   if(abs(anti->id())==ParticleID::tauminus) {
//     for(unsigned int k=0;k<2;++k) {
//       hR[k] =-ya*conj((*neut)(0,2))*(*stau)(k,0) + sqrt(2.)*(*stau)(k,1)*bracketl;
//       hL[k] =-ya*     (*neut)(0,2) *(*stau)(k,1) - sqrt(2.)*(*stau)(k,0)*bracketr;
//     }
//   }
//   else {
//     hR[0] = 0.;
//     hL[0] = - sqrt(2.)*bracketr;
//     hR[1] = + sqrt(2.)*bracketl;
//     hL[1] = 0.;
//   }
//   gL.resize(2,vector<Complex>(2,0.));
//   gR.resize(2,vector<Complex>(2,0.)); 
//   double y = ma/mw/sqrt(2.)/cb;
//   ytau = ma/mw*tb;
//   for(unsigned int i=0;i<2;++i) {
//     if(abs(anti->id())==ParticleID::tauminus) {
//       for(unsigned int j=0;j<2;++j) {
// 	gL[i][j] = 0.;
// 	gR[i][j] = -(*umix)(i,0)*(*stau)(j,0) + ya*(*stau)(j,1)*(*umix)(i,1);
//       } 
//     }
//     else {
//       gL[i][0] = 0.;
//       gR[i][0] = -(*umix)(i,0); 
//       gL[i][1] = 0.;
//       gR[i][1] = 0.;
//     }
//   }
//   double tw = sw/cw;
//   lL.resize(2);
//   lR.resize(2);
//   for(unsigned int j = 0; j < 2; ++j) {
//     lL[j] = -(conj((*neut)(0, 3)*(*vmix)(j,0) + ((*neut)(0,1) + (*neut)(0,0)*tw)*(*vmix)(j,1)/sqrt(2)))*cb;
//     lR[j] = -(     (*neut)(0, 2)*(*umix)(j,0) - ((*neut)(0,1) + (*neut)(0,0)*tw)*(*umix)(j,1)/sqrt(2) )*sb;
//   }
// }

// extracted from SSFFHVertex.cc
// if(particle1->id()!=ParticleID::b) theMassLast.first  = theMSSM->mass(q2,particle1);
// if(particle2->id()!=ParticleID::b) theMassLast.second = theMSSM->mass(q2,particle2);

// extracted from FourBodyDecayConstructor.cc
// from createDecayMode
  // unsigned int nferm=0;
  // tcPDPtr bottom,ferm,anti,chi;
  // for(OrderedParticles::const_iterator it=diagrams[0].outgoing.begin();
  //     it!=diagrams[0].outgoing.end();++it) {
  //   if((**it).iSpin()==PDT::Spin1Half) ++nferm;
  //   if(abs((**it).id())==ParticleID::b)
  //     bottom = *it;
  //   else if(abs((**it).id())>1000000)
  //     chi = *it;
  //   else if((**it).id()<0)
  //     anti = *it;
  //   else
  //     ferm = *it;
  // }
  // if(!bottom||!chi||!ferm||!anti) return;
  // if(ferm->id()-abs(anti->id())!=1) return;
  //if(anti->id()!=ParticleID::tauplus) return;
// from decayList
  // set<PDPtr> new_particles;
  // for(set<PDPtr>::iterator it=particles.begin();it!=particles.end();++it) {
  //   if((**it).id()==ParticleID::SUSY_t_1) new_particles.insert(*it);
  // }
  // NBodyDecayConstructorBase::DecayList(new_particles);

// extracted from StoFFFFDecayer.h

// private :

//   InvEnergy2 stopMatrixElement(const Particle & inpart,
// 			       const ParticleVector & decay,
// 			       InvEnergy2 me2) const;

// #include "Herwig/Models/StandardModel/StandardModel.h"
// InvEnergy2 StoFFFFDecayer::stopMatrixElement(const Particle & inpart,
// 					     const ParticleVector & decay,
// 					     InvEnergy2 me2) const {
//   // extract the momenta and check the process
//   bool found[4]={false,false,false,false};
//   Lorentz5Momentum pb,pf,pfp,pChi;
//   double col = 1.;
//   tcPDPtr ferm,anti;
//   for(unsigned int ix=0;ix<decay.size();++ix) {
//     long id = decay[ix]->id();
//     if(id==ParticleID::b) {
//       found[0] = true;
//       pb    = decay[ix]->momentum();
//     }
//     else if(id==ParticleID::SUSY_chi_10) {
//       found[1] = true;
//       pChi  = decay[ix]->momentum();
//     }
//     else if(abs(id)%2==0) {
//       if(decay[ix]->dataPtr()->coloured()) col = 3.;
//       found[2] = true;
//       pf   = decay[ix]->momentum();
//       ferm = decay[ix]->dataPtr();
//     }
//     else {
//       found[3] = true;
//       pfp  = decay[ix]->momentum();
//       anti = decay[ix]->dataPtr();
//     }
//   }
//   // check the process
//   if(!found[0]||!found[1]||!found[2]||!found[3]) return ZERO;
//   // extract the couplings we need
//   HwSMPtr model = dynamic_ptr_cast<HwSMPtr>(generator()->standardModel());
//   double sw(0.),cw(0.),g(0.);
//   Energy mb = getParticleData(ParticleID::b)->mass();
//   Energy mt = getParticleData(ParticleID::t)->mass();
//   Energy mChi = getParticleData(ParticleID::SUSY_chi_10)->mass();
//   Energy mw = getParticleData(ParticleID::Wplus)->mass();
//   Energy mbt[2] = {getParticleData(ParticleID::SUSY_b_1)->mass(),
// 		   getParticleData(ParticleID::SUSY_b_2)->mass()};
//   Energy mP[2]  = {getParticleData(ParticleID::SUSY_chi_1plus)->mass(),
// 		   getParticleData(ParticleID::SUSY_chi_2plus)->mass()}; 
//   Energy msf[2]={ZERO,ZERO};
//   tcPDPtr sf = getParticleData(1000000+abs(ferm->id()));
//   msf[0] = sf->mass();
//   sf = getParticleData(2000000+abs(ferm->id()));
//   msf[1] = sf ? sf->mass() : 1e30*GeV;
//   Energy msfp[2]={getParticleData(1000000+abs(anti->id()))->mass(),
// 		  getParticleData(2000000+abs(anti->id()))->mass()};
//   Complex aL(0.),aR(0.);
//   vector<Complex> cL,cR,d,bL,bR,kL,kR,fL,fR,hL,hR,lL,lR;
//   vector<vector<Complex> > eL,eR,gL,gR;
//   double ytop,ytau;
//   Energy2 scale = sqr(inpart.mass());
//   model->StopCouplings(scale,ferm,anti,g,sw,cw,aL,aR,cL,cR,d,bL,bR,kL,kR,fL,fR,eL,eR,gL,gR,hL,hR,
// 		       lL,lR,ytop,ytau);
//   // compute the matrix element
//   Lorentz5Momentum pw   = pf+pfp; pw.rescaleMass();
//   Lorentz5Momentum ptop = pw+pb; ptop.rescaleMass();
//   Lorentz5Momentum ptt  = inpart.momentum(); 
//   Lorentz5Momentum pbt  = inpart.momentum()-pw; pbt.rescaleMass();
//   Lorentz5Momentum pChiP= inpart.momentum()-pb; pb.rescaleMass();
//   Lorentz5Momentum psf  = pChi+pf;psf.rescaleMass();
//   Lorentz5Momentum psfp = pChi+pfp;psfp.rescaleMass();
  
//   Energy2 ptpt = ptop*ptop;
//   Energy2 pChipfp  = pChi*pfp;
//   Energy2 ptopptop = ptop.m2();
//   Energy2 pbpf     = pb*pf; 
//   Energy2 pbpfp    = pb*pfp; 
//   Energy2 ptoppfp  = ptop*pfp; 
//   Energy2 pChipt   = pChi*ptop;
//   Energy2 pChiptt   = pChi*ptt;
//   Energy2 pfpfp    = pf*pfp;
//   Energy2 pChipb   = pChi*pb;
//   Energy2 pChipf   = pChi*pf;
//   Energy2 pttptt   = ptt.m2();
//   Energy2 pbptt  = pb*ptt;
//   Energy2 pfpptt  = pfp*ptt;
//   Energy2 pfppt  = pfp*ptop;
//   Energy2 pfptt   = pf*ptt; 
//   Energy2 ptptt   = ptop*ptt; 
//   Energy2 pfpt   = pf*ptop; 
//   Energy2 pbpt   = pb*ptop; 
//   Energy2 pChiPpChiP=pChiP*pChiP;
//   Energy2 pChipChiP=pChi*pChiP;
//   Energy2 pbpChiP = pb*pChiP;
//   Energy2 pfppChiP = pfp*pChiP;
//   Energy2 ptpChiP = ptop*pChiP;
//   Energy2 pttpChiP = ptt*pChiP;
//   Energy2 pfpChiP = pf*pChiP;

//   Energy mf = pf.mass();
//   Energy mfp = pfp.mass();
//   assert(model);
//   InvEnergy2 pTop = 1./(ptop.m2()-sqr(mt));
//   InvEnergy2 pW   = 1./(pw  .m2()-sqr(mw));
//   InvEnergy2 pBT[2]  = {1./(pbt.m2()-sqr(mbt[0])),1./(pbt.m2()-sqr(mbt[1]))};
//   InvEnergy2 pP[2]   = {1./(pChiP.m2()-sqr(mP[0])),1./(pChiP.m2()-sqr(mP[1]))};
//   InvEnergy2 pSF [2]  = {1./(psf .m2()-sqr(msf [0])),1./(psf .m2()-sqr(msf [1]))};
//   if(abs(ferm->id())==ParticleID::nu_e||abs(ferm->id())==ParticleID::nu_mu||abs(ferm->id())==ParticleID::nu_tau)
//     pSF[1] = ZERO;
//   InvEnergy2 pSFP[2]  = {1./(psfp.m2()-sqr(msfp[0])),1./(psfp.m2()-sqr(msfp[1]))};
//   // top squared
//   InvEnergy2 mett = pow(g,6)*sqr(pTop*pW)*
//     ( norm(aR) * ( -4.*pChipfp*pbpf*ptopptop + 8.*pChipt*pbpf*ptoppfp ) +
//       norm(aL) * (  4.*pChipfp*pbpf*sqr(mt))
//       + real(aL*conj(aR)) * (  - 8*pbpf*ptoppfp*mChi*mt )
//       );
//   // colour factors
//   mett *=col;
//   // sbottom squared
//   complex<InvEnergy2> mebb(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     for(unsigned int j=0;j<2;++j) {
//       mebb += pow(g,6)*d[i]*d[j]*pBT[i]*pBT[j]*sqr(pW)*
// 	(pChipb*(cR[i]*cR[j]+cL[i]*cL[j])- mChi*mb*(cL[j]*cR[i]+cL[i]*cR[j]))*
// 	( - 4*pfpfp*pttptt + 8*pfptt*pfpptt + pfpfp*sqr(mfp) + pfpfp*sqr(mf)
// 	  - 4*pfptt*sqr(mfp) - 4*pfpptt*sqr(mf) + 2*sqr(mf)*sqr(mfp) );
//     }
//   }
//   // colour factors
//   mebb *=col;
//   // stop sbottom
//   complex<InvEnergy2> metb(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     metb += pow(g,6)*d[i]*pTop*pBT[i]*sqr(pW)*cR[i]*
//       (
//        + aR*(   - 2*pChipb*pfpfp*ptptt + 2*pChipb*pfpt*
//          pfpptt + 2*pChipb*pfptt*pfppt + 2*pChipf*pbpfp*ptptt - 
//          2*pChipf*pbpt*pfpptt - 2*pChipf*pbptt*pfppt - 2*pChipfp
//          *pbpf*ptptt - 2*pChipfp*pbpt*pfptt + 2*pChipfp*pbptt*
//          pfpt + 2*pChipt*pbpf*pfpptt + 2*pChipt*pbpfp*pfptt - 2*
//          pChipt*pbptt*pfpfp + 2*pChiptt*pbpf*pfppt - 2*pChiptt*
//          pbpfp*pfpt + 2*pChiptt*pbpt*pfpfp)
//        + aL*mChi*mt*(  - 2*pbpf*pfpptt - 2*pbpfp*pfptt  + 2*pbptt*pfpfp ));
//   }
//   // colour factors
//   metb *=col;
//   complex<InvEnergy2> mecc(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     for(unsigned int j=0;j<2;++j) {
//       mecc += pow(g,6)*pP[i]*pP[j]*sqr(pW)*
// 	(+ bR[i]*bR[j]*kR[i]*kR[j] * ( 16*pChipb*pChipf*pChipfp + 16*pChipb*pChipf*pfpfp + 16*pChipb*pChipf*sqr(mfp) + 16*pChipf*
//          pChipfp*pbpf + 16*pChipf*pbpf*pfpfp + 16*pChipf*pbpf*sqr(mfp) + 8*pChipf*pbpfp*sqr(mfp) - 8*pChipf*pbpfp*sqr(mf) - 16*sqr(pChipf)*pbpfp )

//        + bR[i]*bR[j]*kL[i]*kL[j] * ( 8*pChipfp*pbpf*mP[i]*mP[j] )

//        + bL[i]*bL[j]*kR[i]*kR[j] * ( 8*pChipf*pbpfp*mP[i]*mP[j] )

//        + bL[i]*bL[j]*kL[i]*kL[j] * ( 16*pChipb*pChipf*pChipfp + 16*pChipb*
//          pChipfp*pfpfp + 16*pChipb*pChipfp*sqr(mf) + 16*pChipf*
//          pChipfp*pbpfp - 8*pChipfp*pbpf*sqr(mfp) + 8*pChipfp*pbpf*
//          sqr(mf) + 16*pChipfp*pbpfp*pfpfp + 16*pChipfp*pbpfp*sqr(mf) - 
// 				     16*sqr(pChipfp)*pbpf )

//        + mb*bL[j]*bR[i]*kR[i]*kR[j] * (  - 8*pChipf*pChipfp*mP[j] - 8*pChipf*pfpfp*mP[j] - 8*pChipf*sqr(mfp)*mP[j] )

//        + mb*bL[j]*bR[i]*kL[i]*kL[j] * (  - 8*pChipf*pChipfp*mP[i] - 8*pChipfp*pfpfp*mP[i] - 8*pChipfp*sqr(mf)*mP[i] )

//        + mb*bL[i]*bR[j]*kR[i]*kR[j] * (  - 8*pChipf*pChipfp*mP[i] - 8*pChipf*pfpfp*mP[i] - 8*pChipf*sqr(mfp)*mP[i] )

//        + mb*bL[i]*bR[j]*kL[i]*kL[j] * (  - 8*pChipf*pChipfp*mP[j] - 8*pChipfp*pfpfp*mP[j] - 8*pChipfp*sqr(mf)*mP[j] )

//        + mChi*bR[i]*bR[j]*kR[i]*kL[j] * (  - 4*pChipb*pfpfp*mP[j] + 4*pChipf*pbpfp*mP[j] - 4*pChipfp*pbpf*mP[j] - 8*pbpf*pfpfp*mP[j] - 4*pbpf*sqr(mfp)*mP[j] + 4*pbpfp*sqr(mf)*mP[j] )

//        + mChi*bR[i]*bR[j]*kL[i]*kR[j] * (  - 4*pChipb*pfpfp*mP[i] + 4*pChipf*pbpfp*mP[i] - 4*pChipfp*pbpf*mP[i] - 8*pbpf*pfpfp*mP[i] - 4*
//          pbpf*sqr(mfp)*mP[i] + 4*pbpfp*sqr(mf)*mP[i] )

//        + mChi*bL[i]*bL[j]*kR[i]*kL[j] * (  - 4*pChipb*pfpfp*mP[i] - 4*pChipf*pbpfp*mP[i] + 4*pChipfp*pbpf*mP[i] + 4*pbpf*sqr(mfp)*mP[i] - 8*pbpfp*pfpfp*mP[i] - 4*pbpfp*sqr(mf)*mP[i] )

//        + mChi*bL[i]*bL[j]*kL[i]*kR[j] * (  - 4*pChipb*pfpfp*mP[j] - 4*pChipf*pbpfp*mP[j] + 4*pChipfp*pbpf*mP[j] + 4*pbpf*sqr(mfp)*mP[j] - 8*
//          pbpfp*pfpfp*mP[j] - 4*pbpfp*sqr(mf)*mP[j] )

//        + mChi*mb*bL[j]*bR[i]*kR[i]*kL[j] * ( 8*pChipf*pfpfp + 8*pChipfp*
// 					     pfpfp + 4*pfpfp*sqr(mfp) + 4*pfpfp*sqr(mf) + 8*sqr(pfpfp) )

//        + mChi*mb*bL[j]*bR[i]*kL[i]*kR[j] * ( 4*pfpfp*mP[i]*mP[j] )

//        + mChi*mb*bL[i]*bR[j]*kR[i]*kL[j] * ( 4*pfpfp*mP[i]*mP[j] )

//        + mChi*mb*bL[i]*bR[j]*kL[i]*kR[j] * ( 8*pChipf*pfpfp + 8*pChipfp*
// 					 pfpfp + 4*pfpfp*sqr(mfp) + 4*pfpfp*sqr(mf) + 8*sqr(pfpfp) )

//        + sqr(mChi)*bR[i]*bR[j]*kR[i]*kR[j] * (  - 8*pChipf*pbpfp )

//        + sqr(mChi)*bL[i]*bL[j]*kL[i]*kL[j] * (  - 8*pChipfp*pbpf )

//        + mChi*sqr(mChi)*mb*bL[j]*bR[i]*kR[i]*kL[j] * ( 4*pfpfp )

//        + mChi*sqr(mChi)*mb*bL[i]*bR[j]*kL[i]*kR[j] * ( 4*pfpfp ));

//     }
//   }
//   mecc *=col;
//   complex<InvEnergy2> metc(ZERO);
//   for(unsigned int j=0;j<2;++j) {
//     metc +=  pow(g,6)*pP[j]*pTop*sqr(pW)/sqrt(2.)*(
//       + aR*bR[j]*kR[j] * (  - 8*pChipb*pChipf*pbpfp + 8*pChipb*pChipf*pfpfp - 8*pChipb*
// 			pChipfp*pbpf - 8*pChipb*pChipfp*sqr(mf) + 8*pChipb*pbpf*pfpfp - 8*pChipb*pbpfp*
//          sqr(mf) - 4*pChipb*pfpfp*sqr(mfp) - 4*
// 			    pChipb*pfpfp*sqr(mf) - 8*pChipb*sqr(mf)*sqr(mfp) + 8*sqr(pChipb)*pfpfp + 8*pChipf*pChipfp*
//          pbpf + 8*pChipf*pbpf*pbpfp + 16*
//          pChipf*pbpf*pfpfp + 16*pChipf*pbpf*sqr(mfp) + 4*pChipf*pbpfp*sqr(mfp) - 4*pChipf*pbpfp*
// 			    sqr(mf) - 8*sqr(pChipf)*pbpfp + 4*pChipfp*
//          pbpf*sqr(mfp) - 4*pChipfp*pbpf*sqr(mf) - 8*
// 			    pChipfp*sqr(pbpf) )

//        + aR*mb*bL[j]*kR[j] * (  - 4*pChipb*pfpfp*mP[j] - 4*pChipf*pbpfp*mP[j] - 8*pChipf*pfpfp*mP[j]
//           - 4*pChipf*sqr(mfp)*mP[j] + 4*pChipfp*
//          pbpf*mP[j] + 4*pChipfp*sqr(mf)*mP[j] )

//       + aR*sqr(mb)*bR[j]*kR[j] * ( 8*pChipf*pChipfp + 4*pChipf*sqr(mfp) + 4*pChipfp*sqr(mf) )

//        + aR*mChi*bR[j]*kL[j] * (  - 8*pbpf*pbpfp*mP[j] - 8*pbpf*pfpfp*mP[j] - 8*pbpf*sqr(mfp)*mP[j] )

//        + aR*mChi*mb*bL[j]*kL[j] * ( 8*sqr(mf)*sqr(mfp) + 8*
//          pChipf*pbpfp + 8*pChipf*pfpfp + 8*
//          pChipf*sqr(mfp) + 8*pbpfp*pfpfp + 8*
//          pbpfp*sqr(mf) + 8*pfpfp*sqr(mfp) + 8*pfpfp*
// 				    sqr(mf) + 8*sqr(pfpfp) )

//        + aR*sqr(mChi)*bR[j]*kR[j] * ( 8*pbpf*pbpfp + 4
//          *pbpf*sqr(mfp) + 4*pbpfp*sqr(mf) )

//       + aR*sqr(mChi)*sqr(mb)*bR[j]*kR[j] * (  - 4*pfpfp )

//        + aL*mt*bR[j]*kL[j] * ( 8*pChipfp*pbpf*mP[j] )

//        + aL*mb*mt*bL[j]*kL[j] * (  - 8*pChipf*pChipfp - 8*pChipfp*pfpfp - 8*pChipfp*sqr(mf) )

//        + aL*mChi*mt*bR[j]*kR[j] * (  - 4*pChipb*pfpfp + 4*pChipf*pbpfp - 4*pChipfp*pbpf - 8*pbpf*pfpfp - 4*pbpf*sqr(mfp) + 4*pbpfp*sqr(mf) )

//       + aL*mChi*mb*mt*bL[j]*kR[j] * ( 4*pfpfp*mP[j] ));
//   }
//   metc *= col;
  
//   complex<InvEnergy2> mebc(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     for(unsigned int j=0;j<2;++j) {
//       mebc += pow(g,6)*pP[j]*sqr(pW)*d[i]*pBT[i]/sqrt(2.)*
// 	(+ cR[i]*bR[j]*kR[j] * (  - 8*pChipb*pChipf*pfpptt + 4
//          *pChipb*pChipf*sqr(mfp) - 8*pChipb*
//          pChipfp*pfptt + 4*pChipb*pChipfp*sqr(mf) + 8*pChipb*pChiptt*pfpfp + 2*pChipb*
//          pfpfp*sqr(mfp) + 2*pChipb*pfpfp*sqr(mf) - 4
//          *pChipb*pfptt*sqr(mfp) - 4*pChipb*pfpptt*sqr(mf) + 4
//          *pChipb*sqr(mf)*sqr(mfp) + 8*pChipf*pbpfp*
//          pfptt + 2*pChipf*pbpfp*sqr(mfp) - 2*
//          pChipf*pbpfp*sqr(mf) - 8*pChipf*pbptt*pfpfp - 4*pChipf*pbptt*sqr(mfp) - 8*pChipfp*pbpf*
//          pfptt - 2*pChipfp*pbpf*sqr(mfp) + 2*
//          pChipfp*pbpf*sqr(mf) + 4*pChipfp*pbptt*sqr(mf) + 8*pChiptt*pbpf*pfpfp + 4*pChiptt*pbpf*
//          sqr(mfp) - 4*pChiptt*pbpfp*sqr(mf))

//        + cL[i]*bL[j]*kL[j] * (  - 8*pChipb*pChipf*pfpptt + 4
//          *pChipb*pChipf*sqr(mfp) - 8*pChipb*
//          pChipfp*pfptt + 4*pChipb*pChipfp*sqr(mf) + 8*pChipb*pChiptt*pfpfp + 2*pChipb*
//          pfpfp*sqr(mfp) + 2*pChipb*pfpfp*sqr(mf) - 4
//          *pChipb*pfptt*sqr(mfp) - 4*pChipb*pfpptt*sqr(mf) + 4
//          *pChipb*sqr(mf)*sqr(mfp) - 8*pChipf*pbpfp*
//          pfpptt + 2*pChipf*pbpfp*sqr(mfp) - 2*
//          pChipf*pbpfp*sqr(mf) + 4*pChipf*pbptt*sqr(mfp) + 8*pChipfp*pbpf*pfpptt - 2*pChipfp*pbpf
//          *sqr(mfp) + 2*pChipfp*pbpf*sqr(mf) - 8*
//          pChipfp*pbptt*pfpfp - 4*pChipfp*pbptt*sqr(mf) - 4*pChiptt*pbpf*sqr(mfp) + 8*pChiptt*
//          pbpfp*pfpfp + 4*pChiptt*pbpfp*sqr(mf))

//        + mb*cR[i]*bL[j]*kR[j] * ( 4*pChipf*pfpptt*mP[j] - 2*pChipf*sqr(mfp)*mP[j] + 4*pChipfp*pfptt*mP[j]
//           - 2*pChipfp*sqr(mf)*mP[j] - 4*pChiptt*pfpfp*mP[j] )

//        + mb*cL[i]*bR[j]*kL[j] * ( 4*pChipf*pfpptt*mP[j] - 2*pChipf*sqr(mfp)*mP[j] + 4*pChipfp*pfptt*mP[j]
//           - 2*pChipfp*sqr(mf)*mP[j] - 4*pChiptt*pfpfp*mP[j] )

//        + mChi*cR[i]*bR[j]*kL[j] * (  - 4*pbpf*pfpptt*mP[j] + 2*pbpf*sqr(mfp)*mP[j] - 4*pbpfp*pfptt*mP[j]
//           + 2*pbpfp*sqr(mf)*mP[j] + 4*pbptt*pfpfp*mP[j] )

//        + mChi*cL[i]*bL[j]*kR[j] * (  - 4*pbpf*pfpptt*mP[j] + 2*pbpf*sqr(mfp)*mP[j] - 4*pbpfp*pfptt*mP[j] + 2*pbpfp*sqr(mf)*mP[j] + 4*pbptt*pfpfp*mP[j] )

//        + mChi*mb*cR[i]*bL[j]*kL[j] * (  - 4*sqr(mf)*sqr(mfp) + 4*pChipf*pfpptt - 2*pChipf*sqr(mfp) + 4*pChipfp*pfptt - 2*pChipfp*sqr(mf) - 4*pChiptt*pfpfp - 2*pfpfp*sqr(mfp) - 2*pfpfp*sqr(mf) + 4*pfptt*sqr(mfp) + 4*pfpptt*sqr(mf) )

//        + mChi*mb*cL[i]*bR[j]*kR[j] * (  - 4*sqr(mf)*sqr(mfp) + 4*pChipf*pfpptt - 2*pChipf*sqr(mfp) + 4*pChipfp*pfptt - 2*pChipfp*sqr(mf) - 4*pChiptt*pfpfp - 2*pfpfp*sqr(mfp) - 2*pfpfp*sqr(mf) + 4*pfptt*sqr(mfp) + 4*pfpptt*sqr(mf) )

//        + sqr(mChi)*cR[i]*bR[j]*kR[j] * ( 4*pbpf*pfpptt - 2*pbpf*sqr(mfp) + 4*pbpfp*pfptt - 2*pbpfp*sqr(mf) - 4*pbptt*pfpfp )

// 	 + sqr(mChi)*cL[i]*bL[j]*kL[j] * ( 4*pbpf*pfpptt - 2*pbpf*sqr(mfp) + 4*pbpfp*pfptt - 2*pbpfp*sqr(mf) - 4*pbptt*pfpfp ));
//     }
//   }
//   mebc *= col;
//   complex<InvEnergy2> meff(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     for(unsigned int j=0;j<2;++j) {
//       for(unsigned int k=0;k<2;++k) {
// 	for(unsigned int l=0;l<2;++l) {
// 	  meff += pow(g,6)*pP[i]*pP[j]*pSF[k]*pSF[l]*pChipf*(fR[k]*fR[l]+fL[k]*fL[l])*
// 	    (
// 	     + ( eR[i][k]*eR[j][l]*bR[i]*bR[j] + eL[i][k]*eL[j][l]*bL[i]*bL[j] )* ( 4.*pbpfp*mP[i]*mP[j] )
// 	     + ( eR[i][k]*eR[j][l]*bL[i]*bL[j] + eL[i][k]*eL[j][l]*bR[i]*bR[j] )* (  - 4*pbpfp*pChiPpChiP + 
// 	  								    8*pbpChiP*pfppChiP )
// 	     +mb*mP[i]*( eR[i][k]*eR[j][l]*bL[j]*bR[i] + eL[i][k]*eL[j][l]*bL[i]*bR[j]  ) * (  - 4*pfppChiP )
// 	     +mb*mP[j]*( eR[i][k]*eR[j][l]*bL[i]*bR[j] + eL[i][k]*eL[j][l]*bL[j]*bR[i]  ) * (  - 4*pfppChiP ));
// 	}
//       }
//     }
//   }
//   meff *= col;
//   complex<InvEnergy2> mepp(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     for(unsigned int j=0;j<2;++j) {
//       for(unsigned int k=0;k<2;++k) {
// 	for(unsigned int l=0;l<2;++l) {
// 	  mepp += pow(g,6)*pP[i]*pP[j]*pSFP[k]*pSFP[l]*pChipfp*(hR[k]*hR[l]+hL[k]*hL[l])*
// 	    ((gR[i][k]*gR[j][l]*bR[i]*bR[j] + gL[i][k]*gL[j][l]*bL[i]*bL[j]) * ( 4*pbpf*mP[i]*mP[j] ) +
// 	     (gR[i][k]*gR[j][l]*bL[i]*bL[j] + gL[i][k]*gL[j][l]*bR[i]*bR[j]) * 
// 	     (  - 4*pbpf*pChiPpChiP + 8*pbpChiP*pfpChiP ));
// 	}
//       }
//     }
//   }
//   mepp *= col;
//   complex<InvEnergy2> metf(ZERO);
//   for(unsigned int j=0;j<2;++j) {
//     for(unsigned int l=0;l<2;++l) {
//       metf +=  pow(g,6)*pTop*pW*pP[j]*pSF[l]*
// 	(+ aR*eL[j][l]*fR[l]*bR[j] * ( 2*pChipb*pfpfp*ptpChiP - 2*pChipb*pfpt*pfppChiP 
// 				       - 2*pChipb*pfpChiP*pfppt - 2*pChipf*pbpfp*ptpChiP 
// 				       + 2*pChipf*pbpt*pfppChiP + 2*pChipf*pbpChiP*pfppt 
// 				       - 2*pChipfp*pbpf*ptpChiP + 2*pChipfp*pbpt*pfpChiP
// 				       - 2*pChipfp*pbpChiP*pfpt + 2*pChipt*pbpf*pfppChiP 
// 				       - 2*pChipt*pbpfp*pfpChiP + 2*pChipt*pbpChiP*pfpfp 
// 				       + 2*pChipChiP*pbpf*pfppt + 2*pChipChiP*pbpfp*pfpt 
// 				       - 2*pChipChiP*pbpt*pfpfp )
// 	 + aL*mChi*mt*eL[j][l]*fR[l]*bR[j] * (  - 2*pbpf*pfppChiP + 2*pbpfp*pfpChiP - 2*pbpChiP*pfpfp )
	 
// 	 );
//     }
//   }
//   metf *= col;
//   // sbottom sfermion interference
//   complex<InvEnergy2> mebf(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     for(unsigned int j=0;j<2;++j) {
//       for(unsigned int l=0;l<2;++l) {
// 	mebf += 2.*pow(g,6)*d[i]*pBT[i]*pW*pP[j]*pSF[l]*
// 	  (cR[i]*eL[j][l]*fR[l]*bR[j] * ( pChipb*pfpfp*pttpChiP - pChipb*pfptt*pfppChiP - pChipb*pfpChiP*pfpptt +
// 				      pChipf*pbpfp*pttpChiP - pChipf*pbptt*pfppChiP - pChipf*pbpChiP*pfpptt -
// 				      pChipfp*pbpf*pttpChiP + pChipfp*pbptt*pfpChiP - pChipfp*pbpChiP*pfptt +
// 				      pChiptt*pbpf*pfppChiP - pChiptt*pbpfp*pfpChiP + pChiptt*pbpChiP*pfpfp +
// 				      pChipChiP*pbpf*pfpptt + pChipChiP*pbpfp*pfptt - pChipChiP*pbptt*pfpfp)
// 	   + mChi*mP[j]*cL[i]*eL[j][l]*fR[l]*bL[j] * ( - pbpf*pfpptt - pbpfp*pfptt + pbptt*pfpfp ));
//       }
//     }
//   }
//   mebf *= col;
//   // chi W sfermion interference
//   complex<InvEnergy2> mecf(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     for(unsigned int j=0;j<2;++j) {
//       for(unsigned int l=0;l<2;++l) {
// 	mecf -= pow(g,6)*pP[i]*pW*pP[j]*pSF[l]*eL[j][l]*fR[l]/sqrt(2.)*
// 	  (+ bR[i]*bR[j]*kR[i] * ( 8*pChipf*pbpfp*pChiPpChiP - 16*pChipf*pbpChiP*pfppChiP )
// 	   + bR[i]*bR[j]*kL[i]*mChi*mP[i] *  ( 4*pbpf*pfppChiP - 4*pbpfp*pfpChiP + 4*pbpChiP*pfpfp)
// 	   + bL[i]*bL[j]*kR[i]*mP[i]*mP[j] * (  - 8*pChipf*pbpfp )
// 	   + bL[i]*bL[j]*kL[i]*mChi*mP[j]  * (-4*pbpf*pfppChiP + 4*pbpfp*pfpChiP + 4*pbpChiP*pfpfp)
// 	   );
//       }
//     }
//   }
//   mecf *= col;

//   complex<InvEnergy2> mebp(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     for(unsigned int j=0;j<2;++j) {
//       for(unsigned int l=0;l<2;++l) {
// 	mebp += pow(g,6)*d[i]*pBT[i]*pW*pP[j]*pSFP[l]*
// 	  (
// 	   + cL[i]*gR[j][l]*hL[l]*bL[j] * 
// 	   (  - 2*pChipb*pfpfp*pttpChiP + 2*pChipb*pfptt*pfppChiP 
// 	      + 2*pChipb*pfpChiP*pfpptt + 2*pChipf*pbpfp*pttpChiP 
// 	      - 2*pChipf*pbptt*pfppChiP + 2*pChipf*pbpChiP*pfpptt
// 	      - 2*pChipfp*pbpf*pttpChiP + 2*pChipfp*pbptt*pfpChiP
// 	      + 2*pChipfp*pbpChiP*pfptt + 2*pChiptt*pbpf*pfppChiP
// 	      - 2*pChiptt*pbpfp*pfpChiP - 2*pChiptt*pbpChiP*pfpfp
// 	      - 2*pChipChiP*pbpf*pfpptt - 2*pChipChiP*pbpfp*pfptt
// 	      + 2*pChipChiP*pbptt*pfpfp)
				  
// 	   + mChi*cR[i]*gR[j][l]*hL[l]*bR[j]*mP[j] * 
// 	   ( 2*pbpf*pfpptt + 2*pbpfp*pfptt - 2*pbptt*pfpfp ));
//       }
//     }
//   }
//   mebp *= col;
//   // chi W anti sfermion interference
//   complex<InvEnergy2> mecp(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     for(unsigned int j=0;j<2;++j) {
//       for(unsigned int l=0;l<1;++l) {
// 	mecp +=pow(g,6)*pP[i]*pW*pP[j]*pSFP[l]/sqrt(2.)*gR[j][l]*hL[l]*
// 	  (
// 	   + bR[i]*bR[j]*kL[i] * (  - 8*pChipfp*pbpf*mP[i]*mP[j] )
// 	   + bL[i]*bL[j]*kL[i] * ( 8*pChipfp*pbpf*pfppChiP + 8*pChipfp*pbpf*pChiPpChiP - 8*pChipfp*pbpfp*pfpChiP - 8*pChipfp*pbpChiP*pfpfp - 16*pChipfp*pbpChiP*pfpChiP)
// 	   + mChi*bR[i]*bR[j]*kR[i]*mP[j]*( 4*pbpf*pfppChiP - 4*pbpfp*pfpChiP + 4*pbpChiP*pfpfp)
// 	   + mChi*bL[i]*bL[j]*kR[i]*mP[i]*(  - 4*pbpf*pfppChiP + 8*pbpfp*pfpfp + 4*pbpfp*pfpChiP + 4*pbpChiP*pfpfp));
//       }
//     }
//   }
//   mecp *= col;
//   // sfermion antisfermion interferences
//   complex<InvEnergy2> mefp(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     for(unsigned int j=0;j<2;++j) {
//       for(unsigned int k=0;k<1;++k) {
// 	for(unsigned int l=0;l<2;++l) {
// 	  mefp +=pow(g,6)*pP[i]*pP[j]*pSF[k]*pSFP[l]*fR[k]*gR[j][l]*
// 	    (
// 	     + eR[i][k]*hR[l]*bR[i]*bR[j]*mP[i]*mP[j] * 
// 	     ( 2*pChipb*pfpfp - 2*pChipf*pbpfp - 2*pChipfp*pbpf )
	     
// 	     + eR[i][k]*hR[l]*bL[i]*bL[j]* 
// 	     (  - 2*pChipb*pfpfp*pChiPpChiP + 2*pChipf*pbpfp*pChiPpChiP 
// 		- 4*pChipf*pbpChiP*pfppChiP + 2*pChipfp*pbpf*pChiPpChiP 
// 		- 4*pChipfp*pbpChiP*pfpChiP + 4*pChipChiP*pbpChiP*pfpfp)
	     
// 	     + eL[i][k]*hL[l]*bL[i]*bL[j]*mChi*mP[i] * 
// 	     (-2*pbpf*pfppChiP + 2*pbpfp*pfpChiP + 2*pbpChiP*pfpfp )
	     
// 	     + eL[i][k]*hL[l]*bR[i]*bR[j]*mChi*mP[j] *
// 	     ( 2*pbpf*pfppChiP - 2*pbpfp*pfpChiP + 2*pbpChiP*pfpfp ));
// 	}
//       }
//     }
//   }
//   mefp *= col;
//   // top higgs
//   Energy mHiggs = getParticleData(ParticleID::Hplus)->mass();
//   InvEnergy2 pH = 1./(pw.m2()-sqr(mHiggs));
//   InvEnergy2 meht = 0.25*sqr(ytop*ytau)*pow(g,6)*sqr(pTop*pH)*pfpfp*
//     ( norm(aR) *sqr(mt)*4*pChipb +
//       norm(aL) * (  - 4*pChipb*ptpt + 8*pChipt*pbpt )
//       + real(aL*conj(aR)) * mChi*mt * (  - 8*pbpt ));
//   // colour factors
//   meht *=col;
//   // chargino higgs
//   complex<InvEnergy2> mehc(ZERO);
//   for(unsigned int i=0;i<2;++i) {
//     for(unsigned int j=0;j<2;++j) {
//       mehc += pow(g,6)*pP[i]*pP[j]*sqr(pH)*sqr(ytau)*pfpfp*
//       (    + (bR[i]*bR[j]*lR[i]*lR[j]+bL[i]*bL[j]*lL[i]*lL[j])*mP[i]*mP[j]*2*pChipb
//            + (bR[i]*bR[j]*lL[i]*lL[j]+bL[i]*bL[j]*lR[i]*lR[j])*(-2*pChipb*pChiPpChiP+4*pChipChiP*pbpChiP)
//            + (bR[i]*bR[j]*lR[i]*lL[j]+bL[i]*bL[j]*lL[i]*lR[j])*mChi*mP[i]*2*pbpChiP
//            + (bR[i]*bR[j]*lL[i]*lR[j]+bL[i]*bL[j]*lR[i]*lL[j])*mChi*mP[j]*2*pbpChiP);
//     }
//   }


//    // hehbc =
//    //     + cR1*d1*bR2*kR2 * ( 2*pChi.pb*pf.pfp*mP2 )

//    //     + cL1*d1*bL2*kL2 * ( 2*pChi.pb*pf.pfp*mP2 )

//    //     + mChi*cR1*d1*bR2*kL2 * ( 2*pb.pChiP*pf.pfp )

//    //     + mChi*cL1*d1*bL2*kR2 * ( 2*pb.pChiP*pf.pfp );

//   // InvEnergy2 meTotal = meff.real()+mecc.real()+2.*mecf.real();
//   InvEnergy2 meTotal = abs(mehc.real());
//   // if(abs(anti->id())==ParticleID::tauminus) {
//   //   cerr << "testing inter " << (me2-mepp.real()-meff.real())*GeV2 << " " << 2.*mefp.real()*GeV2
//   // 	 << " " << 0.5*(me2-mepp.real()-meff.real())/mefp.real() << "\n";
//   //   cerr << "testing the matrix element " << meTotal*GeV2 << " "
//   // 	 << me2*GeV2 << " " << meTotal/me2 << "\n";
//   // }
//   return meTotal;
// }

// extracted from main diagram loop


	    // //\todo remove testing
	    // top
	    // if(!(abs(inter.first ->id())==ParticleID::t&&
	    // 	 abs(inter.second->id())==ParticleID::Wplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // sbottom
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_b_1||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_b_2)&&
	    // 	 abs(inter.second->id())==ParticleID::Wplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // chargino W
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    // 	 abs(inter.second->id())==ParticleID::Wplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // sneutrino
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    // 	 abs(inter.second->id())!=ParticleID::Wplus&&
	    // 	 inter.second->id()<0)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // sneutrino
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    //  	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    // 	 abs(inter.second->id())!=ParticleID::Wplus&&
	    // 	 inter.second->id()<0)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // charged slepton
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    // 	 abs(inter.second->id())!=ParticleID::Wplus&&
	    // 	 inter.second->id()>0)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // slepton sneutrino
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    // 	 abs(inter.second->id())!=ParticleID::Wplus&&
	    // 	 abs(inter.second->id())!=ParticleID::Hplus&&
	    // 	 inter.second->id()>0) &&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    //  	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    //  	 abs(inter.second->id())!=ParticleID::Wplus&&
	    //  	 abs(inter.second->id())!=ParticleID::Hplus&&
	    //  	 inter.second->id()<0)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // chi W and charged slepton
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    // 	 abs(inter.second->id())==ParticleID::Wplus) &&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    // 	 abs(inter.second->id())!=ParticleID::Wplus&&
	    // 	 inter.second->id()>0) ) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // chi W and sneutrino
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    // 	 abs(inter.second->id())==ParticleID::Wplus) &&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    // 	 abs(inter.second->id())!=ParticleID::Wplus&&
	    // 	 inter.second->id()<0) ) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // top/sbottom interference
	    // if(!(abs(inter.first ->id())==ParticleID::t&&
	    // 	 abs(inter.second->id())==ParticleID::Wplus)&&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_b_1||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_b_2)&&
	    // 	 abs(inter.second->id())==ParticleID::Wplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // top chiW interference
	    // if(!(abs(inter.first ->id())==ParticleID::t&&
	    // 	 abs(inter.second->id())==ParticleID::Wplus)&&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    //  	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    //  	 abs(inter.second->id())==ParticleID::Wplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // bottom chiW interference
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_b_1||
	    //  	  abs(inter.first ->id())==ParticleID::SUSY_b_2)&&
	    // 	 abs(inter.second->id())==ParticleID::Wplus)&&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    //   	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    //   	 abs(inter.second->id())==ParticleID::Wplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // top charged slepton
	    // if(!(abs(inter.first ->id())==ParticleID::t&&
	    // 	 abs(inter.second->id())==ParticleID::Wplus) &&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    //  	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    //  	 abs(inter.second->id())!=ParticleID::Wplus&&
	    //  	 inter.second->id()>0)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // top sneutrino
	    // if(!(abs(inter.first ->id())==ParticleID::t&&
	    // 	 abs(inter.second->id())==ParticleID::Wplus) &&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    //  	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    //  	 abs(inter.second->id())!=ParticleID::Wplus&&
	    //  	 inter.second->id()<0)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // ~b sneutrino
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_b_1||
	    //   	  abs(inter.first ->id())==ParticleID::SUSY_b_2)&&
	    //  	 abs(inter.second->id())==ParticleID::Wplus)&&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    //  	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    //  	 abs(inter.second->id())!=ParticleID::Wplus&&
	    //  	 inter.second->id()<0)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // ~b slepton
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_b_1||
	    //   	  abs(inter.first ->id())==ParticleID::SUSY_b_2)&&
	    //  	 abs(inter.second->id())==ParticleID::Wplus)&&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    //  	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    //  	 abs(inter.second->id())!=ParticleID::Wplus&&
	    //  	 inter.second->id()>0)) {
	    //   ++idiag;
	    //   continue;
	    // }


	    // top H
	    // if(!(abs(inter.first ->id())==ParticleID::t&&
	    // 	 abs(inter.second->id())==ParticleID::Hplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // sbottom H
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_b_1||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_b_2)&&
	    // 	 abs(inter.second->id())==ParticleID::Hplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // chargino H
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    // 	 abs(inter.second->id())==ParticleID::Hplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // top/sbottom interference
	    // if(!(abs(inter.first ->id())==ParticleID::t&&
	    // 	 abs(inter.second->id())==ParticleID::Hplus)&&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_b_1||
	    // 	  abs(inter.first ->id())==ParticleID::SUSY_b_2)&&
	    // 	 abs(inter.second->id())==ParticleID::Hplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // top chiH interference
	    // if(!(abs(inter.first ->id())==ParticleID::t&&
	    // 	 abs(inter.second->id())==ParticleID::Hplus)&&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    //  	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    //  	 abs(inter.second->id())==ParticleID::Hplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // bottom chiH interference
	    // if(!((abs(inter.first ->id())==ParticleID::SUSY_b_1||
	    //  	  abs(inter.first ->id())==ParticleID::SUSY_b_2)&&
	    // 	 abs(inter.second->id())==ParticleID::Hplus)&&
	    //    !((abs(inter.first ->id())==ParticleID::SUSY_chi_1plus||
	    //   	  abs(inter.first ->id())==ParticleID::SUSY_chi_2plus)&&
	    //   	 abs(inter.second->id())==ParticleID::Hplus)) {
	    //   ++idiag;
	    //   continue;
	    // }
	    // all Higgs
	    // if(abs(inter.second->id())==ParticleID::Hplus) {
	    //   ++idiag;
	    //   continue;
	    // }
//reomved from end of me2()
//InvEnergy2 output = stopMatrixElement(inpart,decay,me2*UnitRemoval::InvE2);
  // return output*scale;
