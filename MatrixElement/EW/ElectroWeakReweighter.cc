// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ElectroWeakReweighter class.
//

#include "ElectroWeakReweighter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/operation.hpp"
#include "EWProcess.h"
#include "HighEnergyMatching.h"
#include "ElectroWeakMatching.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;

tEWCouplingsPtr ElectroWeakReweighter::staticEWCouplings_ = tEWCouplingsPtr();


ElectroWeakReweighter::ElectroWeakReweighter() {}

ElectroWeakReweighter::~ElectroWeakReweighter() {}

IBPtr ElectroWeakReweighter::clone() const {
  return new_ptr(*this);
}

IBPtr ElectroWeakReweighter::fullclone() const {
  return new_ptr(*this);
}

void ElectroWeakReweighter::persistentOutput(PersistentOStream & os) const {
  os << EWCouplings_ << collinearSudakov_ << softSudakov_;
}

void ElectroWeakReweighter::persistentInput(PersistentIStream & is, int) {
  is >> EWCouplings_ >> collinearSudakov_ >> softSudakov_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ElectroWeakReweighter,ReweightBase>
  describeHerwigElectroWeakReweighter("Herwig::ElectroWeakReweighter", "HwMEEW.so");

void ElectroWeakReweighter::Init() {

  static ClassDocumentation<ElectroWeakReweighter> documentation
    ("There is no documentation for the ElectroWeakReweighter class");

  static Reference<ElectroWeakReweighter,EWCouplings> interfaceEWCouplings
    ("EWCouplings",
     "The object to calculate the electroweak couplings",
     &ElectroWeakReweighter::EWCouplings_, false, false, true, false, false);

  static Reference<ElectroWeakReweighter,CollinearSudakov> interfaceCollinearSudakov
    ("CollinearSudakov",
     "The collinear Sudakov",
     &ElectroWeakReweighter::collinearSudakov_, false, false, true, false, false);

  static Reference<ElectroWeakReweighter,SoftSudakov> interfaceSoftSudakov
    ("SoftSudakov",
     "The soft Sudakov",
     &ElectroWeakReweighter::softSudakov_, false, false, true, false, false);

}

namespace {

void axpy_prod_local(const boost::numeric::ublas::matrix<Complex>  & A,
		     const boost::numeric::ublas::matrix<complex<InvEnergy2> > & B,
		     boost::numeric::ublas::matrix<complex<InvEnergy2> > & C) {
  assert(A.size2()==B.size1());
  C.resize(A.size1(),B.size2());
  for(unsigned int ix=0;ix<A.size1();++ix) {
    for(unsigned int iy=0;iy<B.size2();++iy) {
      C(ix,iy) = ZERO;
      for(unsigned int iz=0;iz<A.size2();++iz) {
	C(ix,iy) += A(ix,iz)*B(iz,iy);
      }
    }
  }
}

void axpy_prod_local(const boost::numeric::ublas::matrix<complex<InvEnergy2> >  & A,
		     const boost::numeric::ublas::vector<complex<Energy2> >      & B,
		     boost::numeric::ublas::vector<Complex > & C) {
  assert(A.size2()==B.size());
  C.resize(A.size1());
  for(unsigned int ix=0;ix<A.size1();++ix) {
    C(ix) = ZERO;
      for(unsigned int iz=0;iz<A.size2();++iz) {
	C(ix) += A(ix,iz)*B(iz);
      }
  }
}

void axpy_prod_local(const boost::numeric::ublas::matrix<complex<InvEnergy2> > & A,
		     const boost::numeric::ublas::matrix<Complex> & B,
		     boost::numeric::ublas::matrix<complex<InvEnergy2> > & C) {
  assert(A.size2()==B.size1());
  C.resize(A.size1(),B.size2());
  for(unsigned int ix=0;ix<A.size1();++ix) {
    for(unsigned int iy=0;iy<B.size2();++iy) {
      C(ix,iy) = ZERO;
      for(unsigned int iz=0;iz<A.size2();++iz) {
	C(ix,iy) += A(ix,iz)*B(iz,iy);
      }
    }
  }
}

}

double ElectroWeakReweighter::weight() const {
  EWCouplings_->initialize();
  staticEWCouplings_ = EWCouplings_;
  // cerr << "aEM\n";
  // for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.1) {
  //   cerr << scale/GeV << " " 
  // 	 << EWCouplings_->aEM(scale) << "\n";
  // }
  // cerr << "aS\n";
  // for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.4) {
  //   cerr << scale/GeV << " " 
  // 	 << EWCouplings_->aS(scale) << "\n";
  // }
  // cerr << "y_t\n";
  // for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.4) {
  //   cerr << scale/GeV << " " 
  // 	 << EWCouplings_->y_t(scale) << "\n";
  // }
  // cerr << "lambda\n";
  // for(Energy scale=91.2*GeV; scale<10*TeV; scale *= 1.4) {
  //   cerr << scale/GeV << " " 
  // 	 << EWCouplings_->lambda(scale) << "\n";
  // }
  // cerr << "vev\n";
  // for(Energy scale=91.2*GeV; scale<10*TeV; scale *= 1.4) {
  //   cerr << scale/GeV << " " 
  // 	 << EWCouplings_->vev(scale)/GeV << "\n";
  // }
  // collinearSudakov_->makePlots();
  // Energy2 s = sqr(5000.*GeV);
  // Energy2 t = -0.25*s;
  // Energy2 u = -0.75*s;
  // testEvolution(s,t,u);




  cerr <<  subProcess() << "\n";
  cerr << *subProcess() << "\n";
  cerr << subProcess()->outgoing()[0] << *subProcess()->outgoing()[0] << "\n";
  cerr << subProcess()->outgoing()[0]->spinInfo() << "\n";
  cerr << subProcess()->outgoing()[0]->spinInfo()->productionVertex() << "\n";
  if(subProcess()->outgoing().size()!=2)
    return 1.;
  if(subProcess()->incoming().first->id()==ParticleID::g &&
     subProcess()->incoming().second->id()==ParticleID::g) {
    if(subProcess()->outgoing()[0]->id()==ParticleID::g &&
       subProcess()->outgoing()[1]->id()==ParticleID::g)
      return 1.;
    else if(abs(subProcess()->outgoing()[0]->id())<=5 &&
	    subProcess()->outgoing()[0]->id()==-subProcess()->outgoing()[1]->id()) {
      return reweightggqqbar();
    }
    else
      assert(false);
  }
  else if(abs(subProcess()->incoming().first->id())<=5 &&
	  subProcess()->incoming().first->id()==-subProcess()->incoming().second->id()) {
    if(subProcess()->outgoing()[0]->id()==ParticleID::g &&
       subProcess()->outgoing()[1]->id()==ParticleID::g)
      return reweightqqbargg();
    else
      assert(false);
  }
  else
    assert(false);
  assert(false);
  staticEWCouplings_ = tEWCouplingsPtr();
}

void ElectroWeakReweighter::testEvolution(Energy2 s,Energy2 t, Energy2 u) const {
  Energy highScale = sqrt(s);
  Energy ewScale = coupling()->mZ();
  Energy lowScale = 50.0*GeV;
  for (unsigned int i=0; i<45;++i) {
    EWProcess::Process process = (EWProcess::Process)i;
    cerr << "process " << process << "\n";
    // result for all EW and QCD SCET contributions:
    boost::numeric::ublas::matrix<complex<InvEnergy2> > highMatch_val 
      = HighEnergyMatching::highEnergyMatching(highScale,s,t,u,process,true,true);
    boost::numeric::ublas::matrix<Complex> highRunning_val
      = softSudakov_->highEnergyRunning(highScale,ewScale,s,t,u,process);
    boost::numeric::ublas::matrix<Complex> ewMatch_val = 
      ElectroWeakMatching::electroWeakMatching(ewScale,s,t,u,process,true);
    boost::numeric::ublas::matrix<Complex> lowRunning_val = 
      softSudakov_->lowEnergyRunning(ewScale,lowScale,s,t,u,process);
    boost::numeric::ublas::matrix<Complex> collinearHighRunning_val =
      collinearSudakov_->highEnergyRunning(highScale,ewScale,s,process,false);
    boost::numeric::ublas::matrix<Complex> collinearEWMatch_val =
      collinearSudakov_->electroWeakMatching(ewScale,s,process,true);
    boost::numeric::ublas::matrix<Complex> collinearLowRunning_val =
      collinearSudakov_->lowEnergyRunning(ewScale,lowScale,s,process);
    boost::numeric::ublas::matrix<Complex> lowMatchTemp_val = 
      boost::numeric::ublas::zero_matrix<Complex>(ewMatch_val.size1(),ewMatch_val.size2());
    for (unsigned int ii=0; ii<ewMatch_val.size1(); ++ii) {
      for (unsigned int jj=0; jj<ewMatch_val.size2(); ++jj) {
	lowMatchTemp_val(ii,jj) = collinearEWMatch_val(ii,jj)*ewMatch_val(ii,jj);
      }
    }
    boost::numeric::ublas::matrix<Complex> temp(highRunning_val.size1(),collinearHighRunning_val.size2());
    boost::numeric::ublas::axpy_prod(highRunning_val,collinearHighRunning_val,temp);
    boost::numeric::ublas::matrix<Complex> temp2(collinearLowRunning_val.size1(),lowRunning_val.size2());
    boost::numeric::ublas::axpy_prod(collinearLowRunning_val,lowRunning_val,temp2);
    boost::numeric::ublas::matrix<Complex> temp3(temp2.size1(),lowMatchTemp_val.size2());
    boost::numeric::ublas::axpy_prod(temp2,lowMatchTemp_val,temp3);
    temp2.resize(temp3.size1(),temp.size2());
    boost::numeric::ublas::axpy_prod(temp3,temp,temp2);
    boost::numeric::ublas::matrix<complex<InvEnergy2> >  result(temp2.size1(),highMatch_val.size2());
    axpy_prod_local(temp2,highMatch_val,result);
    for(unsigned int ix=0;ix<result.size1();++ix) {
      for(unsigned int iy=0;iy<result.size2();++iy) {
	cerr << s*result(ix,iy) << " ";
      }
      cerr << "\n";
    }
  }
}

namespace {

void SackGluonPolarizations(Lorentz5Momentum &p1,
			    Lorentz5Momentum &p2,
			    Lorentz5Momentum &p3,
			    Lorentz5Momentum &p4,
			    Energy2 s, Energy2 t, Energy2 u, Energy2 m2,
			    vector<LorentzVector<Complex> > & eps3,
			    vector<LorentzVector<Complex> > & eps4,
			    unsigned int iopt) {
  static const Complex I(0.,1.);
  // p1 is p-, p2 is p+
  // p3 is k-, p4 is k+
  // both final-state
  if(iopt==0) {
    // swap t and u due Aneesh's defn
    Energy3 den1 = sqrt((u*t-sqr(m2))*(s-4.*m2));
    Energy3 den2 = sqrt(s*(u*t-sqr(m2)));
    LorentzVector<Complex> eps3Para = (m2+t)/den1*p1 -(m2+u)/den1*p2 +(u-t)/den1*p3;
    LorentzVector<Complex> eps3Perp = 2./den2*epsilon(p1,p2,p3);
    LorentzVector<Complex> eps4Para = (m2+t)/den1*p2 -(m2+u)/den1*p1 +(u-t)/den1*p4;
    LorentzVector<Complex> eps4Perp = 2./den2*epsilon(p1,p2,p4);
    eps3.push_back(sqrt(0.5)*(eps3Para+I*eps3Perp));
    eps3.push_back(sqrt(0.5)*(eps3Para-I*eps3Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para+I*eps4Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para-I*eps4Perp));
    if(m2!=ZERO) assert(false);
  }
  // both initial-state
  else if(iopt==1) {
    if(m2!=ZERO) assert(false);
    LorentzVector<Complex> eps3Para( 1., 0.,0.,0.);
    LorentzVector<Complex> eps3Perp( 0.,-1.,0.,0.);
    LorentzVector<Complex> eps4Para(-1.,0.,0., 0.);
    LorentzVector<Complex> eps4Perp( 0., 1.,0.,0.);
    eps3.push_back(sqrt(0.5)*(eps3Para+I*eps3Perp));
    eps3.push_back(sqrt(0.5)*(eps3Para-I*eps3Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para+I*eps4Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para-I*eps4Perp));
  }
}

}

double ElectroWeakReweighter::reweightqqbargg() const {
  // momenta and invariants
  Lorentz5Momentum p1   = subProcess()->incoming().first ->momentum();
  tcPDPtr          q    = subProcess()->incoming().first ->dataPtr();
  Lorentz5Momentum p2   = subProcess()->incoming().second->momentum();
  tcPDPtr          qbar = subProcess()->incoming().second->dataPtr();
  if(subProcess()->incoming().first->id()<0) {
    swap(p1,p2  );
    swap(q ,qbar);
  }
  Lorentz5Momentum p3 = subProcess()->outgoing()[0]->momentum();
  Lorentz5Momentum p4 = subProcess()->outgoing()[1]->momentum();
  tcPDPtr           g = subProcess()->outgoing()[1]->dataPtr();
  Energy2 s = (p1+p2).m2();
  Energy2 t = (p1-p4).m2();
  Energy2 u = (p1-p3).m2();
  // // boost to partonci rest frame
  // Lorentz5Momentum psum=p1+p2;
  // LorentzRotation boost(-psum.boostVector());
  // p1 *= boost;
  // p2 *= boost;
  // p3 *= boost;
  // p4 *= boost;
  // cerr << "testing momenta in reweight A " << p1/GeV << "\n";
  // cerr << "testing momenta in reweight B " << p2/GeV << "\n";
  // cerr << "testing momenta in reweight C " << p3/GeV << "\n";
  // cerr << "testing momenta in reweight D " << p4/GeV << "\n";
  // LO matrix element coefficents
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    bornQQGGweights,bornRRGGweights;
  // quark left doublet
  if(q->id()!=5) {
    bornQQGGweights = evaluateRunning(EWProcess::QQGG,s,t,u,true);
  }
  else {
    bornQQGGweights = evaluateRunning(EWProcess::QtQtGG,s,t,u,true);
  }
  // quark right singlet
  if(abs(subProcess()->incoming().first->id())%2==0)
    bornRRGGweights = evaluateRunning(EWProcess::UUGG,s,t,u,true);
  else
    bornRRGGweights = evaluateRunning(EWProcess::DDGG,s,t,u,true);
  // EW corrected matrix element coefficients
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    EWQQGGweights,EWRRGGweights;
  // quark left doublet
  if(q->id()!=5) {
    EWQQGGweights = evaluateRunning(EWProcess::QQGG,s,t,u,false);
  }
  else {
    EWQQGGweights = evaluateRunning(EWProcess::QtQtGG,s,t,u,false);
  }
  // quakr right singlet
  if(abs(subProcess()->incoming().first->id())%2==0)
    EWRRGGweights = evaluateRunning(EWProcess::UUGG,s,t,u,false);
  else
    EWRRGGweights = evaluateRunning(EWProcess::DDGG,s,t,u,false);
  // cerr << "testing matrices\n";
  // for(unsigned int ix=0;ix<bornRRGGweights.size1();++ix) {
  //   for(unsigned int iy=0;iy<bornRRGGweights.size2();++iy) {
  //     cerr << bornRRGGweights(ix,iy)*GeV2 << " ";
  //   }
  //   cerr << "\n";
  // }
  // for(unsigned int ix=0;ix<bornQQGGweights.size1();++ix) {
  //   for(unsigned int iy=0;iy<bornQQGGweights.size2();++iy) {
  //     cerr << bornQQGGweights(ix,iy)*GeV2 << " ";
  //   }
  //   cerr << "\n";
  // }

  SpinorWaveFunction       qw(p1,q   ,incoming);
  SpinorBarWaveFunction qbarw(p2,qbar,incoming);
  // VectorWaveFunction      g1w(p3,getParticleData(ParticleID::g),outgoing);
  // VectorWaveFunction      g2w(p4,getParticleData(ParticleID::g),outgoing);
  vector<LorentzVector<Complex> > eps3,eps4;
  SackGluonPolarizations(p1,p2,p3,p4,s,t,u,ZERO,eps3,eps4,0);
  boost::numeric::ublas::matrix<Complex>
    bornME = boost::numeric::ublas::zero_matrix<Complex>(3,3),
    EWME   = boost::numeric::ublas::zero_matrix<Complex>(3,3);
    // testME = boost::numeric::ublas::zero_matrix<Complex>(3,3);
  for(unsigned int iq=0;iq<2;++iq) {
    if(iq==0) {
      qw.reset   (0);
      qbarw.reset(1);
    }
    else {
      qw.reset   (1);
      qbarw.reset(0);
    }
    LorentzVector<complex<Energy> > current  = iq==0 ?
      qw.dimensionedWave(). leftCurrent(qbarw.dimensionedWave()) :
      qw.dimensionedWave().rightCurrent(qbarw.dimensionedWave());
    for(unsigned int i1=0;i1<2;++i1) {
      complex<Energy> d31 = eps3[i1].dot(p1);
      for(unsigned int i2=0;i2<2;++i2) {
	// g1w.reset(2*i1);
	// g2w.reset(2*i2);
	boost::numeric::ublas::vector<complex<Energy2> > M(5);
	Complex         d34 = eps3[i1].dot(eps4[i2]);
	complex<Energy> d42 = eps4[i2].dot(p2);
	// M0 in paper
	M(0) = qw.dimensionedWave().slash(eps3[i1])
	  .slash(p4-p2).vectorCurrent(qbarw.dimensionedWave()).dot(eps4[i2]);

	// // really t channel
	// Complex MEU = -Complex(0.,1.)*qw.dimensionedWave().slash(g1w.wave())
	//   .slash(p4-p2).vectorCurrent(qbarw.dimensionedWave()).dot(g2w.wave())/u;
	// // really u channel
	// Complex MET = -Complex(0.,1.)*qw.dimensionedWave().slash(g2w.wave())
	//   .slash(p3-p2).vectorCurrent(qbarw.dimensionedWave()).dot(g1w.wave())/t;
	// // s channel
	// Complex MES = -Complex(0.,1.)*(current.dot(p3-p4)*g1w.wave().dot(g2w.wave())
	// 			       -2.*current.dot(g1w.wave())*(g2w.wave().dot(p3))
	// 			       -2.*current.dot(g2w.wave())*(g1w.wave().dot(p4)))/s;

	// // really t channel
	// Complex MEU = -Complex(0.,1.)*qw.dimensionedWave().slash(eps3[i1])
	//   .slash(p4-p2).vectorCurrent(qbarw.dimensionedWave()).dot(eps4[i2])/u;
	// // really u channel
	// Complex MET = -Complex(0.,1.)*qw.dimensionedWave().slash(eps4[i2])
	//   .slash(p3-p2).vectorCurrent(qbarw.dimensionedWave()).dot(eps3[i1])/t;
	// // s channel
	// Complex MES = -Complex(0.,1.)*(current.dot(p3-p4)*eps3[i1].dot(eps4[i2])
	// 			       -2.*current.dot(eps3[i1])*(eps4[i2].dot(p3))
	// 			       -2.*current.dot(eps4[i2])*(eps3[i1].dot(p4)))/s;
	// cerr << "NEW flows " << MEU+MES << " " << MET-MES << "\n";




	// cerr << "testing new U " << MET << "\n";
	// M4 in paper
	M(2) =  current.dot(eps4[i2])*d31;
	// M5 in paper
	M(3) = -current.dot(eps3[i1])*d42;
	// M1 in paper (missing factor)
	M(1) =  current.dot(p4);
	// M6 in paper
	M(4) = M(1)*d31*d42/GeV2;
	// M1 final factor
	M(1) *= d34;
	// coefficient of different contributions
	boost::numeric::ublas::vector<Complex> Cborn(3),CEW(3),Ctest(3);


	// Ctest(0) = 1./6.*( MEU+MET);
	// Ctest(1) = 0.5*( MEU+MET);
	// Ctest(2) = 0.5*(MEU+MES-MET+MES);
	if(iq==0) {
	  axpy_prod_local(bornQQGGweights,M,Cborn);
	  axpy_prod_local(EWQQGGweights  ,M,CEW  );
	}
	else {
	  axpy_prod_local(bornRRGGweights,M,Cborn);
	  axpy_prod_local(EWRRGGweights  ,M,CEW  );
	}
	unsigned int ioff = (Cborn.size()==6 && q->id()%2!=0) ? 3 : 0;
	for(unsigned int ix=0;ix<3;++ix) {
	  for(unsigned int iy=0;iy<3;++iy) {
	    bornME(ix,iy) += Cborn(ix+ioff)*conj(Cborn(iy+ioff));
	    EWME  (ix,iy) += CEW  (ix+ioff)*conj(CEW  (iy+ioff));
	    // testME  (ix,iy) += Ctest  (ix)*conj(Ctest  (iy));
	  }
	}
      }
    }
  }
  double born = 24.*real(bornME(0,0))+20./3.*real(bornME(1,1))+12.*real(bornME(2,2));
  double EW   = 24.*real(EWME(0,0))+20./3.*real(EWME(1,1))+12.*real(EWME(2,2));
  // double test = 24.*real(testME(0,0))+20./3.*real(testME(1,1))+12.*real(testME(2,2));
  
  // double gs2 = 4.*Constants::pi*ElectroWeakReweighter::coupling()->a3(sqrt(s));

  // cerr << "testing born A " << 0.125*born/sqr(gs2)/9. << "\n";
  // cerr << "testing born B " << 0.125*test/9. << "\n";

  return EW/born;
}

boost::numeric::ublas::matrix<complex<InvEnergy2> > 
ElectroWeakReweighter::evaluateRunning(EWProcess::Process process, Energy2 s,
				       Energy2 t, Energy2 u, bool born) const {
  using namespace boost::numeric::ublas;
  bool SU3save = coupling()->SU3();
  bool EWsave  = coupling()-> EW();
  Energy highScale = sqrt(s);
  Energy   ewScale = coupling()->mZ();
  Energy  lowScale = ewScale;
  // result for all EW and QCD SCET contributions:
  // high energy matching
  matrix<complex<InvEnergy2> > highMatch_val 
    = HighEnergyMatching::highEnergyMatching(highScale,s,t,u,process,!born,false);
  // low energy matching
  matrix<Complex> 
    ewMatch_val = ElectroWeakMatching::electroWeakMatching(ewScale,s,t,u,process,!born);
  matrix<Complex> collinearEWMatch_val =
    collinearSudakov_->electroWeakMatching(ewScale,s,process,!born);
  matrix<Complex> highRunning_val,lowRunning_val,
    collinearHighRunning_val,collinearLowRunning_val;
  if(born) {
    highRunning_val = identity_matrix<Complex>(softSudakov_->numberGauge(process));
    lowRunning_val  = identity_matrix<Complex>(softSudakov_->numberBrokenGauge(process));
    collinearHighRunning_val = identity_matrix<Complex>(softSudakov_->numberGauge(process));
    collinearLowRunning_val  = identity_matrix<Complex>(softSudakov_->numberBrokenGauge(process));
  }
  else {
    coupling()->SU3(false);
    coupling()-> EW( true);
    highRunning_val = softSudakov_->highEnergyRunning(highScale,ewScale,s,t,u,process);
    lowRunning_val = softSudakov_->lowEnergyRunning(ewScale,lowScale,s,t,u,process);
    collinearHighRunning_val = collinearSudakov_->highEnergyRunning(highScale,ewScale,s,process,false);
    collinearLowRunning_val = collinearSudakov_->lowEnergyRunning(ewScale,lowScale,s,process);
  };
  matrix<Complex> lowMatchTemp_val = 
    zero_matrix<Complex>(ewMatch_val.size1(),ewMatch_val.size2());
  for (unsigned int ii=0; ii<ewMatch_val.size1(); ++ii) {
    for (unsigned int jj=0; jj<ewMatch_val.size2(); ++jj) {
      lowMatchTemp_val(ii,jj) = collinearEWMatch_val(ii,jj)*ewMatch_val(ii,jj);
    }
  }
  // perform all the multiplications
  matrix<Complex> temp(highRunning_val.size1(),collinearHighRunning_val.size2());
  axpy_prod(highRunning_val,collinearHighRunning_val,temp);
  matrix<Complex> temp2(collinearLowRunning_val.size1(),lowRunning_val.size2());
  axpy_prod(collinearLowRunning_val,lowRunning_val,temp2);
  matrix<Complex> temp3(temp2.size1(),lowMatchTemp_val.size2());
  axpy_prod(temp2,lowMatchTemp_val,temp3);
  temp2.resize(temp3.size1(),temp.size2());
  axpy_prod(temp3,temp,temp2);
  matrix<complex<InvEnergy2> >  result(temp2.size1(),highMatch_val.size2());
  axpy_prod_local(temp2,highMatch_val,result);
  // for(unsigned int ix=0;ix<result.size1();++ix) {
  //   for(unsigned int iy=0;iy<result.size2();++iy) {
  // 	cerr << s*result(ix,iy) << " ";
  //   }
  //   cerr << "\n";
  // }
  // reset the couplings
  coupling()->SU3(SU3save);
  coupling()-> EW( EWsave);
  // return the answer
  return result;
}


double ElectroWeakReweighter::reweightggqqbar() const {
  // momenta and invariants
  Lorentz5Momentum p1   = subProcess()->incoming().first ->momentum();
  Lorentz5Momentum p2   = subProcess()->incoming().second->momentum();
  Lorentz5Momentum p3   = subProcess()->outgoing()[0]->momentum();
  Lorentz5Momentum p4   = subProcess()->outgoing()[1]->momentum();
  tcPDPtr          qbar = subProcess()->outgoing()[0]->dataPtr();
  tcPDPtr          q    = subProcess()->outgoing()[1]->dataPtr();
  if(q->id()<0) {
    swap(p3,p4  );
    swap(q ,qbar);
  }
  Energy2 s = (p1+p2).m2();
  Energy2 t = (p1-p4).m2();
  Energy2 u = (p1-p3).m2();
  // boost to partonic rest frame and rescale momenta of outgoing
  // so zero mass
  Lorentz5Momentum psum=p1+p2;
  LorentzRotation boost(-psum.boostVector());
  p1 *= boost;
  p2 *= boost;
  p3 *= boost;
  p4 *= boost;
  p3.setMass(ZERO);
  p3.rescaleRho();
  p4.setMass(ZERO);
  p4.rescaleRho();
  // LO matrix element coefficents
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    bornQQGGweights,bornRRGGweights;
  // quark left doublet
  if(q->id()<5) {
    bornQQGGweights = evaluateRunning(EWProcess::QQGG,s,t,u,true);
  }
  else {
    bornQQGGweights = evaluateRunning(EWProcess::QtQtGG,s,t,u,true);
  }
  // quark right singlet
  if(q->id()==0) {
    if(q->id()==6)
      bornRRGGweights = evaluateRunning(EWProcess::tRtRGG,s,t,u,true);
    else
      bornRRGGweights = evaluateRunning(EWProcess::UUGG,s,t,u,true);
  }
  else
    bornRRGGweights = evaluateRunning(EWProcess::DDGG,s,t,u,true);
  // EW corrected matrix element coefficients
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    EWQQGGweights,EWRRGGweights;
  // quark left doublet
  if(q->id()<5) {
    EWQQGGweights = evaluateRunning(EWProcess::QQGG,s,t,u,false);
  }
  else {
    EWQQGGweights = evaluateRunning(EWProcess::QtQtGG,s,t,u,false);
  }
  // quark right singlet
  if(q->id()%2==0) {
    if(q->id()==6)
      EWRRGGweights = evaluateRunning(EWProcess::tRtRGG,s,t,u,false);
    else
      EWRRGGweights = evaluateRunning(EWProcess::UUGG,s,t,u,false);
  }
  else
    EWRRGGweights = evaluateRunning(EWProcess::DDGG,s,t,u,false);
  SpinorWaveFunction       qw(p4,qbar,incoming);
  SpinorBarWaveFunction qbarw(p3,q   ,incoming);
  vector<LorentzVector<Complex> > eps1,eps2;
  SackGluonPolarizations(p1,p2,p3,p4,s,t,u,ZERO,eps1,eps2,1);
  boost::numeric::ublas::matrix<Complex>
    bornME = boost::numeric::ublas::zero_matrix<Complex>(3,3),
    EWME   = boost::numeric::ublas::zero_matrix<Complex>(3,3);
  // helicities of outgoing quarks
  for(unsigned int iq=0;iq<2;++iq) {
    if(iq==0) {
      qw.reset   (0);
      qbarw.reset(1);
    }
    else {
      qw.reset   (1);
      qbarw.reset(0);
    }
    LorentzVector<complex<Energy> > current  = iq==0 ?
      qw.dimensionedWave(). leftCurrent(qbarw.dimensionedWave()) :
      qw.dimensionedWave().rightCurrent(qbarw.dimensionedWave());
    for(unsigned int i1=0;i1<2;++i1) {
      complex<Energy> d31 = eps1[i1].dot(p3);
      for(unsigned int i2=0;i2<2;++i2) {
  	boost::numeric::ublas::vector<complex<Energy2> > M(5);
  	Complex         d34 = eps1[i1].dot(eps2[i2]);
  	complex<Energy> d42 = eps2[i2].dot(p4);
  	// M0 in paper
  	M(0) = qw.dimensionedWave().slash(eps1[i1])
  	  .slash(p2-p4).vectorCurrent(qbarw.dimensionedWave()).dot(eps2[i2]);
  	// M4 in paper
  	M(2) =  current.dot(eps2[i2])*d31;
  	// M5 in paper
  	M(3) = -current.dot(eps1[i1])*d42;
  	// M1 in paper (missing factor)
  	M(1) =  current.dot(p2);
  	// M6 in paper
  	M(4) = M(1)*d31*d42/GeV2;
  	// M1 final factor
  	M(1) *= d34;
  	// coefficient of different contributions
 	boost::numeric::ublas::vector<Complex> Cborn(3),CEW(3);
  	if(iq==0) {
  	  axpy_prod_local(bornQQGGweights,M,Cborn);
  	  axpy_prod_local(EWQQGGweights  ,M,CEW  );
  	}
  	else {
  	  axpy_prod_local(bornRRGGweights,M,Cborn);
  	  axpy_prod_local(EWRRGGweights  ,M,CEW  );
  	}
	unsigned int ioff = (Cborn.size()==6 && q->id()%2!=0) ? 3 : 0;
  	for(unsigned int ix=0;ix<3;++ix) {
  	  for(unsigned int iy=0;iy<3;++iy) {
  	    bornME(ix,iy) += Cborn(ix+ioff)*conj(Cborn(iy+ioff));
  	    EWME  (ix,iy) += CEW  (ix+ioff)*conj(CEW  (iy+ioff));
  	  }
  	}
      }
    }
  }
  double born = 24.*real(bornME(0,0))+20./3.*real(bornME(1,1))+12.*real(bornME(2,2));
  double EW   = 24.*real(EWME(0,0))+20./3.*real(EWME(1,1))+12.*real(EWME(2,2));
  return EW/born;
}
