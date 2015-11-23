// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2gZ2qqPowheg class.
//

#include "MEee2gZ2llPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig/Utilities/Maths.h"

using namespace Herwig;

IBPtr MEee2gZ2llPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEee2gZ2llPowheg::fullclone() const {
  return new_ptr(*this);
}

void MEee2gZ2llPowheg::persistentOutput(PersistentOStream & os) const {
  os << contrib_ << yPow_ << zPow_;
}

void MEee2gZ2llPowheg::persistentInput(PersistentIStream & is, int) {
  is >> contrib_ >> yPow_ >> zPow_;
}

ClassDescription<MEee2gZ2llPowheg> MEee2gZ2llPowheg::initMEee2gZ2llPowheg;
// Definition of the static class description member.

void MEee2gZ2llPowheg::Init() {

  static ClassDocumentation<MEee2gZ2llPowheg> documentation
    ("The MEee2gZ2llPowheg class implements the next-to-leading order "
     "matrix element for e+e- > q qbar in the POWHEG scheme");

  static Switch<MEee2gZ2llPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEee2gZ2llPowheg::contrib_, 1, false, false);
  static SwitchOption interfaceContributionLeadingOrder
    (interfaceContribution,
     "LeadingOrder",
     "Just generate the leading order cross section",
     0);
  static SwitchOption interfaceContributionPositiveNLO
    (interfaceContribution,
     "PositiveNLO",
     "Generate the positive contribution to the full NLO cross section",
     1);
  static SwitchOption interfaceContributionNegativeNLO
    (interfaceContribution,
     "NegativeNLO",
     "Generate the negative contribution to the full NLO cross section",
     2);

  static Parameter<MEee2gZ2llPowheg,double> interfacezPower
    ("zPower",
     "The sampling power for z",
     &MEee2gZ2llPowheg::zPow_, 0.5, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<MEee2gZ2llPowheg,double> interfaceyPower
    ("yPower",
     "The sampling power for y",
     &MEee2gZ2llPowheg::yPow_, 0.9, 0.0, 1.0,
     false, false, Interface::limited);

}

int MEee2gZ2llPowheg::nDim() const {
  return MEee2gZ2ll::nDim() + ( contrib_>0 ? 3 : 0 );
}

bool MEee2gZ2llPowheg::generateKinematics(const double * r) {
  // Generate the radiative integration variables:
  if(contrib_>0) {
    z_   = r[nDim()-1];
    y_   = r[nDim()-2];
    phi_ = r[nDim()-3]*Constants::twopi;
  }
  // Continue with lo matrix element code:
  return MEee2gZ2ll::generateKinematics(r);
}

double MEee2gZ2llPowheg::me2() const {
  // if leading order just return the LO matrix element
  if(contrib_==0) return MEee2gZ2ll::me2();
  // cast the vertices
  tcFFVVertexPtr Zvertex = dynamic_ptr_cast<tcFFVVertexPtr>(FFZVertex());
  tcFFVVertexPtr Pvertex = dynamic_ptr_cast<tcFFVVertexPtr>(FFPVertex());
  // compute the spinors
  vector<SpinorWaveFunction> fin,aout;
  vector<SpinorBarWaveFunction>  ain,fout;
  SpinorWaveFunction    ein  (rescaledMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction pin  (rescaledMomenta()[1],mePartonData()[1],incoming);
  SpinorBarWaveFunction qkout(rescaledMomenta()[2],mePartonData()[2],outgoing);
  SpinorWaveFunction    qbout(rescaledMomenta()[3],mePartonData()[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    ein.reset(ix)  ;
    fin.push_back( ein  );
    pin.reset(ix)  ;
    ain.push_back( pin  );
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  // me to be returned
  ProductionMatrixElement   output(PDT::Spin1Half,PDT::Spin1Half,
				   PDT::Spin1Half,PDT::Spin1Half);
  ProductionMatrixElement  gammaME(PDT::Spin1Half,PDT::Spin1Half,
				   PDT::Spin1Half,PDT::Spin1Half);
  ProductionMatrixElement ZbosonME(PDT::Spin1Half,PDT::Spin1Half,
				   PDT::Spin1Half,PDT::Spin1Half);
  ProductionMatrixElement    nloME(PDT::Spin1Half,PDT::Spin1Half,
				   PDT::Spin1Half,PDT::Spin1Half);
  // wavefunctions for the intermediate particles
  VectorWaveFunction interZ,interG;
  // temporary storage of the different diagrams
  // sum over helicities to get the matrix element
  unsigned int inhel1,inhel2,outhel1,outhel2;
  double total[4]={0.,0.,0.,0.};
  LorentzPolarizationVector momDiff = 
    (rescaledMomenta()[2]-rescaledMomenta()[3])/2./
    (rescaledMomenta()[2].mass()+rescaledMomenta()[3].mass());
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      // intermediate Z
      interZ = FFZVertex()->evaluate(scale(),1,Z0(),fin[inhel1],ain[inhel2]);
      // intermediate photon
      interG = FFPVertex()->evaluate(scale(),1,gamma(),fin[inhel1],ain[inhel2]);
      // scalars
      Complex scalar1 = interZ.wave().dot(momDiff);
      Complex scalar2 = interG.wave().dot(momDiff);
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {		
	  // first the Z exchange diagram
	  Complex diag1 = FFZVertex()->evaluate(scale(),aout[outhel2],fout[outhel1],
						interZ);
	  // then the photon exchange diagram
	  Complex diag2 = FFPVertex()->evaluate(scale(),aout[outhel2],fout[outhel1],
						interG);
	  // extra stuff for NLO
	  LorentzPolarizationVector left  = 
	    aout[outhel2].wave(). leftCurrent(fout[outhel1].wave());
	  LorentzPolarizationVector right = 
	    aout[outhel2].wave().rightCurrent(fout[outhel1].wave());
	  Complex scalar = 
	    aout[outhel2].wave().scalar(fout[outhel1].wave());
	  // nlo specific pieces
	  Complex diag3 =
	           Complex(0.,1.)*Zvertex->norm()*
	    (Zvertex->right()*( left.dot(interZ.wave())) +
	     Zvertex-> left()*(right.dot(interZ.wave())) -
	     ( Zvertex-> left()+Zvertex->right())*scalar1*scalar);
	  diag3 += Complex(0.,1.)*Pvertex->norm()*
	    (Pvertex->right()*( left.dot(interG.wave())) +
	     Pvertex-> left()*(right.dot(interG.wave())) -
	     ( Pvertex-> left()+Pvertex->right())*scalar2*scalar);
	  // add up squares of individual terms
	  total[1] += norm(diag1);
	  ZbosonME(inhel1,inhel2,outhel1,outhel2) = diag1;
	  total[2] += norm(diag2);
	  gammaME (inhel1,inhel2,outhel1,outhel2) = diag2;
	  // the full thing including interference
	  diag1 += diag2;
	  total[0] += norm(diag1);
	  output(inhel1,inhel2,outhel1,outhel2)=diag1;
	  // nlo piece
	  total[3] += real(diag1*conj(diag3) + diag3*conj(diag1));
	  nloME(inhel1,inhel2,outhel1,outhel2)=diag3;
	}
      }
    }
  }
  // spin average
  for(int ix=0;ix<4;++ix) total[ix] *= 0.25;
  // special for polarization beams if needed
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = 
      {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
       beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    total[0] = output  .average(rho[0],rho[1]);
    total[1] = ZbosonME.average(rho[0],rho[1]);
    total[2] = gammaME .average(rho[0],rho[1]);
    total[3] = real(output.average(nloME ,rho[0],rho[1]) +
		    nloME .average(output,rho[0],rho[1]));
  }
  // save the stuff for diagram selection
  DVector save;
  save.push_back(total[2]);
  save.push_back(total[1]);
  meInfo(save);
  // now for the NLO bit
  double mu2 = 0.25*sqr(rescaledMomenta()[2].mass()+rescaledMomenta()[3].mass())/sHat();
  double mu = sqrt(mu2);
  double mu4 = sqr(mu2);
  double lmu = log(mu);
  double v = sqrt(1.-4.*mu2),v2(sqr(v));
  double omv = 4.*mu2/(1.+v);
  double f1,f2,fNS,VNS;
  double r = omv/(1.+v),lr(log(r));
  // normal form
  if(mu>1e-4) {
    f1 = 
      ( +1. + 3.*log(0.5*(1.+v)) - 1.5*log(0.5*(1.+v2)) + sqr(Constants::pi)/6.
	- 0.5*sqr(lr) - (1.+v2)/v*(lr*log(1.+v2) + sqr(Constants::pi)/12. 
				       -0.5*log(4.*mu2)*lr + 0.25*sqr(lr)));
    fNS =  -0.5*(1.+2.*v2)*lr/v + 1.5*lr - 2./3.*sqr(Constants::pi) + 0.5*sqr(lr)
      + (1.+v2)/v*(Herwig::Math::ReLi2(r) + sqr(Constants::pi)/3. - 0.25*sqr(lr) + lr*log((2.*v/ (1.+v))));
    VNS = 1.5*log(0.5*(1.+v2)) 
      + 0.5*(1.+v2)/v*( 2.*lr*log(2.*(1.+v2)/sqr(1.+v))  + 2.*Herwig::Math::ReLi2(sqr(r)) 
			    - 2.*Herwig::Math::ReLi2(2.*v/(1.+v)) - sqr(Constants::pi)/6.)
      + log(1.-mu) - 2.*log(1.-2.*mu) - 4.*mu2/(1.+v2)*log(mu/(1.-mu)) - mu/(1.-mu)
      + 4.*(2.*mu2-mu)/(1.+v2) + 0.5*sqr(Constants::pi); 
    f2 = mu2*lr/v;
  }
  // small mass limit
  else {
    f1 = -1./6.*
      ( - 6. - 24.*lmu*mu2 - 15.*mu4 - 12.*mu4*lmu - 24.*mu4*sqr(lmu) 
	+ 2.*mu4*sqr(Constants::pi) - 12.*mu2*mu4 - 96.*mu2*mu4*sqr(lmu) 
	+ 8.*mu2*mu4*sqr(Constants::pi) - 80.*mu2*mu4*lmu);
    fNS = - mu2/18.*( + 36.*lmu - 36. - 45.*mu2 + 216.*lmu*mu2 - 24.*mu2*sqr(Constants::pi) 
		      + 72.*mu2*sqr(lmu) - 22.*mu4 + 1032.*mu4 * lmu
		      - 96.*mu4*sqr(Constants::pi) + 288.*mu4*sqr(lmu));
    VNS = - mu2/1260.*(-6930. + 7560.*lmu + 2520.*mu - 16695.*mu2 + 1260.*mu2*sqr(Constants::pi) 
		       + 12600.*lmu*mu2 + 1344.*mu*mu2 - 52780.*mu4 + 36960.*mu4*lmu 
		       + 5040.*mu4*sqr(Constants::pi) - 12216.*mu*mu4);
    f2 = mu2*( 2.*lmu + 4.*mu2*lmu + 2.*mu2 + 12.*mu4*lmu + 7.*mu4);
  }
  // add up bits for f1
  f1 += fNS+VNS;
  // now for the real correction
  double jac = 1.;
  // generate y
  double yminus = 0.; 
  double yplus  = 1.-2.*mu*(1.-mu)/(1.-2*mu2);
  double rhoymax = pow(yplus-yminus,1.-yPow_);
  double rhoy = y_*rhoymax;
  double y = yminus+pow(rhoy,1./(1.-yPow_));
  jac *= pow(y-yminus,yPow_)*rhoymax/(1.-yPow_);
  // generate z 
  double vt = sqrt(max(sqr(2.*mu2+(1.-2.*mu2)*(1.-y))-4.*mu2,0.))/(1.-2.*mu2)/(1.-y);
  double zplus  = (1.+vt)*(1.-2.*mu2)*y/2./(mu2 +(1.-2.*mu2)*y);
  double zminus = (1.-vt)*(1.-2.*mu2)*y/2./(mu2 +(1.-2.*mu2)*y);
  double rhozmax = pow(zplus-zminus,1.-zPow_);
  double rhoz = z_*rhozmax;
  double z = zminus+pow(rhoz,1./(1.-zPow_));
  jac *= pow(z-zminus,zPow_)*rhozmax/(1.-zPow_);
  // calculate x1,x2,x3 and xT 
  double x2 = 1. - y*(1.-2.*mu2);
  double x1 = 1. - z*(x2-2.*mu2);
  double x3 = 2.-x1-x2;
  double xT = sqrt(max(0.,sqr(x3) -0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1))/(sqr(x2)-4.*mu2)));
  // calculate the momenta
  Energy M = sqrt(sHat());
  Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*mu2,0.)),0.5*M*x2,M*mu); 
  Lorentz5Momentum pemit (-0.5*M*xT*cos(phi_),-0.5*M*xT*sin(phi_),
			  0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*mu2,0.)),0.5*M*x1,M*mu);
  Lorentz5Momentum pgluon( 0.5*M*xT*cos(phi_), 0.5*M*xT*sin(phi_),
			   0.5*M*sqrt(max(sqr(x3)-sqr(xT),0.)),0.5*M*x3,ZERO);
  if(abs(pspect.z()+pemit.z()-pgluon.z())/M<1e-6) 
    pgluon.setZ(-pgluon.z());
  else if(abs(pspect.z()-pemit.z()+pgluon.z())/M<1e-6) 
    pemit .setZ(- pemit.z());
  // loop over the possible emitting partons
  vector<cPDPtr> partons(mePartonData());
  partons.push_back(cPDPtr());
  double realwgt=0.;
  for(unsigned int iemit=0;iemit<2;++iemit) {
    // boost and rotate momenta
    LorentzRotation eventFrame( ( rescaledMomenta()[2] +
				  rescaledMomenta()[3] ).findBoostToCM() );
    Lorentz5Momentum spectator = eventFrame*rescaledMomenta()[2+iemit];
    eventFrame.rotateZ( -spectator.phi() );
    eventFrame.rotateY( -spectator.theta()  );
    eventFrame.invert();
    vector<Lorentz5Momentum> momenta(rescaledMomenta());
    if(iemit==0) {
      momenta[3] = eventFrame*pspect;
      momenta[2] = eventFrame*pemit ;
    }
    else {
      momenta[2] = eventFrame*pspect;
      momenta[3] = eventFrame*pemit ;
    }
    momenta.push_back(eventFrame*pgluon);
    // calculate the weight
    if(1.-x1>1e-5 && 1.-x2>1e-5) {
      partons[4] = gamma();
      realwgt += meRatio(partons,momenta,iemit,true);
    }
  }
  // total real emission contribution
  double realFact = 0.25*(1.-y)*jac*sqr(1.-2.*mu2)/sqrt(1.-4.*mu2);
  realwgt *= realFact;
  // coupling prefactors
  double charge = sqr(double(mePartonData()[2]->iCharge())/3.);
  double aEM =    SM().alphaEM(scale())/Constants::pi*charge;
  // coupling factors
  // correction for real emission
  realwgt *= aEM;
  // and virtual
  f1 *= aEM;
  f2 *= aEM;
  // the born + virtual + real
  total[0] = total[0]*(1. + f1 +realwgt) + f2*total[3];
  // return the answer
  if(contrib_==2) total[0] *=-1.;
  return max(total[0],0.);
}
