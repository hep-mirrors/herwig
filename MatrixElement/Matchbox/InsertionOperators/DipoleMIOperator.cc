// -*- C++ -*-
//
// DipoleMIOperator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleMIOperator class.
//

#include "DipoleMIOperator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "ThePEG/PDT/DecayMode.h"

#include "Herwig++/MatrixElement/Matchbox/MatchboxFactory.h"

#include <gsl/gsl_sf_dilog.h>

using namespace Herwig;
using Constants::pi;

DipoleMIOperator::DipoleMIOperator() 
  : MatchboxInsertionOperator(),
    CA(-1.0), CF(-1.0), 
    gammaQuark(-1.0), gammaGluon(-1.0),
    betaZero(-1.),
    KQuark(-1.0), KGluon(-1.0) {}

DipoleMIOperator::~DipoleMIOperator() {}

IBPtr DipoleMIOperator::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleMIOperator::fullclone() const {
  return new_ptr(*this);
}

void DipoleMIOperator::setXComb(tStdXCombPtr xc) {
  MatchboxInsertionOperator::setXComb(xc);
  if ( CA < 0. ) {
    CA = SM().Nc();
    CF = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc());
    gammaQuark = (3./2.)*CF;
    // gammaGluon = (11./6.)*CA - (1./3.)*lastBorn()->nLight();
    gammaGluon = (11./6.)*CA - (1./3.)*lastBorn()->nLightJetVec().size();
    // betaZero = gammaGluon;
    betaZero = (11./6.)*CA - (1./3.)*(lastBorn()->nLightJetVec().size()+lastBorn()->nHeavyJetVec().size());
    KQuark = (7./2.-sqr(pi)/6.)*CF;
    // KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*lastBorn()->nLight();
    KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*lastBorn()->nLightJetVec().size();
  }
}

bool DipoleMIOperator::apply(const cPDVector& pd) const {

  // DipoleMIOperator should only apply, if massive
  // partons exist in the given process (aka in the
  // final state):
  // This applies also to gluons in the final state of
  // the LO process and further splittings g->Q\bar{Q}
  // inside the jet.
  bool mFSet = false;
  if ( NHeavyJetVec().size() != 0 ) {
    mFSet = true;
  }

  // Partons in the initial state are assumed massless in
  // the CS approach:
  // The following loop checks for at least one existing
  // combination (note that the single apply function is
  // not checking for massless condition) 'n in addition
  // for at least one massive parton in the final state,
  // 'n for only massless partons in the initial state.
  bool first = false;
  bool second = false;
  bool finalmass = false;
  bool initialmass = false;
  int idp = 0;
  for ( cPDVector::const_iterator p = pd.begin();
	p != pd.end(); ++p, ++idp ) {
    if ( (*p)->coloured() && (*p)->mass()!=ZERO && idp > 1 ) {
      finalmass = true;
    }
    if ( (*p)->coloured() && (*p)->mass()!=ZERO && idp < 2 ) {
      initialmass = true;
    }
    if ( !first ) {
      if ( apply(*p) ) {
	first = true;
      }
    } else {
      if ( apply(*p) ) {
	second = true;
      }
    }
  }

  if ( first && second && (finalmass || mFSet) && !initialmass ) {
    cout << "DipoleMIOperator::apply (master apply): Return true!" << endl;
    cout << endl;
    cout << "     !!!!! Attention !!!!!" << endl;
    cout << "     Number of massless flavours in jet particle group (aka n_f) = " << NLightJetVec().size() << endl;
    cout << "     Number of massive flavours in jet particle group (aka n_F or n_{f,h}) = " << NHeavyJetVec().size() << endl;
    cout << "     Ensure consistent usage!" << endl;
    cout << endl;
  }
  if ( !( first && second && (finalmass || mFSet) && !initialmass ) )
    cout << "DipoleMIOperator::apply (master apply): Return false!" << endl;

  return first && second && (finalmass || mFSet) && !initialmass;

}

bool DipoleMIOperator::apply(tcPDPtr pd) const {
  return
    // pd->mass() == ZERO &&
    (abs(pd->id()) < 7 || pd->id() == ParticleID::g);
}

vector<int> DipoleMIOperator::NLightJetVec() const {

  const map<string,PDVector>& theParticleGroups = MatchboxFactory::currentFactory()->particleGroups();
  map<string,PDVector>::const_iterator theIt = theParticleGroups.find("j");
  if ( theIt == theParticleGroups.end() )
    throw Exception() << "DipoleMIOperator::NLightJetVec(): Could not find a jet particle group named 'j'" << Exception::abortnow;

  const PDVector& theJetConstitutents = theIt->second;
  vector<int> theNLightJetVec;

  for ( PDVector::const_iterator theP = theJetConstitutents.begin();
        theP != theJetConstitutents.end(); ++theP ) {
    if ( (**theP).id() > 0 && (**theP).id() < 7 && (**theP).mass() == ZERO )
      theNLightJetVec.push_back( (**theP).id() );
  }

  return theNLightJetVec;

}

vector<int> DipoleMIOperator::NHeavyJetVec() const {

  const map<string,PDVector>& theParticleGroups = MatchboxFactory::currentFactory()->particleGroups();
  map<string,PDVector>::const_iterator theIt = theParticleGroups.find("j");
  if ( theIt == theParticleGroups.end() )
    throw Exception() << "DipoleMIOperator::NHeavyJetVec(): Could not find a jet particle group named 'j'" << Exception::abortnow;

  const PDVector& theJetConstitutents = theIt->second;
  vector<int> theNHeavyJetVec;

  for ( PDVector::const_iterator theP = theJetConstitutents.begin();
        theP != theJetConstitutents.end(); ++theP ) {
    if ( (**theP).id() > 0 && (**theP).id() < 7 && (**theP).mass() != ZERO )
      theNHeavyJetVec.push_back( (**theP).id() );
  }

  return theNHeavyJetVec;

}

vector<int> DipoleMIOperator::NLightBornVec() const {

  // For the moment just count all quark and antiquark
  // constituents in the Born process.

  vector<int> theNLightBornVec;

  for ( cPDVector::const_iterator j = mePartonData().begin();
	j != mePartonData().end(); ++j ) {
    // if ( (**j).id() > 0 && (**j).id() < 7 && (**j).mass() == ZERO )
    if ( abs((**j).id()) < 7 && (**j).mass() == ZERO )
      theNLightBornVec.push_back( (**j).id() );
  }

  return theNLightBornVec;

}

vector<int> DipoleMIOperator::NHeavyBornVec() const {

  // For the moment just count all quark and antiquark
  // constituents in the Born process.

  vector<int> theNHeavyBornVec;

  for ( cPDVector::const_iterator j = mePartonData().begin();
	j != mePartonData().end(); ++j ) {
    // if ( (**j).id() > 0 && (**j).id() < 7 && (**j).mass() != ZERO )
    if ( abs((**j).id()) < 7 && (**j).mass() != ZERO )
      theNHeavyBornVec.push_back( (**j).id() );
  }

  return theNHeavyBornVec;

}

double DipoleMIOperator::me2() const {

  if ( !isExpanded() )
    throw InitException() << "DipoleMIOperator only implemented in the expanded convention.";
  if ( isDR() )
    throw InitException() << "DipoleMIOperator not implemented for dimensional reduction.";

  Energy2 mu2 = lastBorn()->mu2();

  double kappa=0.;

  double res = 0.;

  int idj = 0; int idk = 0;

  // j is emitter, k is spectator

  for ( cPDVector::const_iterator j = mePartonData().begin();
	j != mePartonData().end(); ++j, ++idj ) {

    if ( !apply(*j) ) {
      continue;
    }

    if ( apply(*j) && idj < 2 && (**j).mass() != ZERO )
      throw InitException() << "DipoleMIOperator: Initial state partons must not be massive!";

    idk = 0;

    for ( cPDVector::const_iterator k = mePartonData().begin();
	  k != mePartonData().end(); ++k, ++idk ) {

      if ( !apply(*k) ) {
        continue;
      }

      if ( j == k || lastBorn()->noDipole(idj,idk) ) {
        continue;
      }

      if ( apply(*k) && idk < 2 && (**k).mass() != ZERO )
        throw InitException() << "DipoleMIOperator: Initial state partons must not be massive!";

      double delta = 0.0;

      Energy2 sjk = 2.*meMomenta()[idj]*meMomenta()[idk];
      
      if ( idj > 1 ) { // Involves idk > 1 as well as idk < 2
        delta +=
	  ( ((**j).id() == ParticleID::g ? CA : CF) *
	    ( Vj(**j,**k,sjk,kappa) - sqr(pi)/3. ) +
	    ((**j).id() == ParticleID::g ? GammaGluon() : GammaQuark(**j)) +
	    ((**j).id() == ParticleID::g ? gammaGluon : gammaQuark) * (1 + log(mu2/sjk)) +
	    ((**j).id() == ParticleID::g ? KGluon : KQuark) );
      }

      if ( idj < 2 && idk > 1 ) {
        delta +=
	  ( ((**j).id() == ParticleID::g ? CA : CF) *
	    ( Vj(**j,**k,sjk,2./3.,true) - sqr(pi)/3. ) +
	    // ((**j).id() == ParticleID::g ? GammaGluon() : GammaQuark(**j)) + // Actually the finite parts of GammaGluon and GammaQuark are zero for m_j=0
	    ((**j).id() == ParticleID::g ? gammaGluon : gammaQuark) * (1 + log(mu2/sjk)) +
	    ((**j).id() == ParticleID::g ? KGluon : KQuark) );
      }

      // If j and k are incoming, same contribution as in DipoleIOperator.
      // Master apply prevents the DipoleIOperator though from applying in
      // case of at least one massive parton. So, add expanded finite term
      // for additional initial-initial correlation here.
      if ( idj < 2 && idk < 2 ) {
        delta += ( ((**j).id() == ParticleID::g ? CA : CF) * ( -1.*sqr(pi)/3. + 1./2.*log(mu2/sjk)*log(mu2/sjk) ) +
                 ((**j).id() == ParticleID::g ? gammaGluon : gammaQuark) * ( 1. + log(mu2/sjk) ) +
                 ((**j).id() == ParticleID::g ? KGluon : KQuark) );
      }

      delta *= lastBorn()->colourCorrelatedME2(make_pair(idj,idk));

      res += delta;

    }
  }

  // NOTE: The accountance for the alpha_s running has yet to be checked for the
  // massive case, or rather to be made consistent within the whole framework in
  // general. So far proceed similarly to the massless case here (no DR though),
  // however betaZero contains not only the light flavours (see setXComb()) now,
  // and also consider arXiv::hep-ph/0011222v3 eq.(27).

//   Energy2 muR2 = 
//     lastBorn()->renormalizationScale()*
//     sqr(lastBorn()->renormalizationScaleFactor());
//   if ( muR2 != mu2 ) {
//     // throw InitException() << "DipoleMIOperator::me2(): muR2 != mu2. Scale dependence restoration for alpha_s not yet understood in the massive case. Please use MatchboxSHatScale as renormalization scale for now.";
//     res -=
//       betaZero *
//       lastBorn()->orderInAlphaS() * 
//       log(muR2/mu2) * lastBorn()->me2();
//       for ( cPDVector::const_iterator j = mePartonData().begin();
//             j != mePartonData().end(); ++j ) {
//         if ( !(**j).coloured() ) continue;
//         res += 
//           ( (**j).id()==ParticleID::g ? gammaGluon-1./3.*lastBorn()->nHeavyJetVec().size() : gammaQuark ) * 
//           log(muR2/mu2) * lastBorn()->me2();
//       }
//   }

  Energy2 muR2 = 
    lastBorn()->renormalizationScale()*
    sqr(lastBorn()->renormalizationScaleFactor());
  if ( muR2 != mu2 ) {
    res -=
      betaZero *
      lastBorn()->orderInAlphaS() * log(muR2/mu2) *
      lastBorn()->me2();
  }

//   // include the finite renormalization for DR here; ATTENTION this
//   // has to be mentioned in the manual!  see hep-ph/9305239 for
//   // details; this guarantees an expansion in alpha_s^\bar{MS} when
//   // using dimensional reduction
//   if ( isDR() && isDRbar() )
//     res -= (CA/6.)*lastBorn()->orderInAlphaS()*lastBorn()->me2();

  res *= ( - lastBorn()->lastAlphaS() / (2.*pi) );
  
  return res;

}

double DipoleMIOperator::oneLoopDoublePole() const {

  if ( !isExpanded() )
    throw InitException() << "DipoleMIOperator only implemented in the expanded convention.";
  if ( isDR() )
    throw InitException() << "DipoleMIOperator not implemented for dimensional reduction.";

  double res = 0.;

  int idj = 0; int idk = 0;

  // j is emitter, k is spectator

  for ( cPDVector::const_iterator j = mePartonData().begin();
	j != mePartonData().end(); ++j, ++idj ) {

    if ( !apply(*j) ) {
      continue;
    }

    if ( apply(*j) && idj < 2 && (**j).mass() != ZERO )
      throw InitException() << "DipoleMIOperator: Initial state partons must not be massive!";

    idk = 0;

    for ( cPDVector::const_iterator k = mePartonData().begin();
	  k != mePartonData().end(); ++k, ++idk ) {

      if ( !apply(*k) ) {
        continue;
      }

      if ( j == k || lastBorn()->noDipole(idj,idk) ) {
        continue;
      }

      if ( apply(*k) && idk < 2 && (**k).mass() != ZERO )
        throw InitException() << "DipoleMIOperator: Initial state partons must not be massive!";

      double delta = 0.0;

      if (idj>1) { // Involves idk > 1 as well as idk < 2
        delta += ( (**j).id() == ParticleID::g ? CA : CF ) * VsDoublePole(**j,**k);
      }
      else if (idj<2 && idk>1) {
        delta += ( (**j).id() == ParticleID::g ? CA : CF ) * VsDoublePole(**j,**k);
      }
      else if (idj<2 && idk<2) {
        delta += ( (**j).id() == ParticleID::g ? CA : CF );
      }

      delta *= lastBorn()->colourCorrelatedME2(make_pair(idj,idk));

      res += delta;

    }
  }

  res *= ( - lastBorn()->lastAlphaS() / (2.*pi) );
  
  return res;

}

double DipoleMIOperator::oneLoopSinglePole() const {

  if ( !isExpanded() )
    throw InitException() << "DipoleMIOperator only implemented in the expanded convention.";
  if ( isDR() )
    throw InitException() << "DipoleMIOperator not implemented for dimensional reduction.";

  Energy2 mu2 = lastBorn()->mu2();

  double res = 0.;

  int idj = 0; int idk = 0;

  // j is emitter, k is spectator

  for ( cPDVector::const_iterator j = mePartonData().begin();
	j != mePartonData().end(); ++j, ++idj ) {

    if ( !apply(*j) ) {
      continue;
    }

    if ( apply(*j) && idj < 2 && (**j).mass() != ZERO )
      throw InitException() << "DipoleMIOperator: Initial state partons must not be massive!";

    idk = 0;

    for ( cPDVector::const_iterator k = mePartonData().begin();
	  k != mePartonData().end(); ++k, ++idk ) {

      if ( !apply(*k) ) {
        continue;
      }

      if ( j == k || lastBorn()->noDipole(idj,idk) ) {
        continue;
      }

      if ( apply(*k) && idk < 2 && (**k).mass() != ZERO )
        throw InitException() << "DipoleMIOperator: Initial state partons must not be massive!";

      double delta = 0.0;

      Energy2 sjk = 2.*meMomenta()[idj]*meMomenta()[idk];

      if (idj>1) { // Involves idk > 1 as well as idk < 2
        delta += ( (**j).id() == ParticleID::g ? CA : CF ) * VsSinglePole(**j,**k,sjk);
        delta += ( (**j).id() == ParticleID::g ? GammaGluonSinglePole() : GammaQuarkSinglePole(**j) );
      }
      else if (idj<2 && idk>1) {
        delta += ( (**j).id() == ParticleID::g ? CA : CF ) * VsSinglePole(**j,**k,sjk);
        delta += ( (**j).id() == ParticleID::g ? gammaGluon : gammaQuark );
      }
      else if (idj<2 && idk<2) {
        delta += ( (**j).id() == ParticleID::g ? CA : CF ) * log(mu2/sjk);
        delta += ( (**j).id() == ParticleID::g ? gammaGluon : gammaQuark );
      }

      delta *= lastBorn()->colourCorrelatedME2(make_pair(idj,idk));

      res += delta;

    }
  }

  res *= ( - lastBorn()->lastAlphaS() / (2.*pi) );

  return res;

}

double DipoleMIOperator::Vj(const ParticleData& j, const ParticleData& k,
			    Energy2 sjk, double kappa, bool mFSetEmpty) const {
  
  Energy2 mu2 = lastBorn()->mu2();

  double res = 0.;
  
  // sjk is being handed over as input parameter to DipoleMIOperator::Vj()
  // kappa is being handed over as input parameter to DipoleMIOperator::Vj()

  Energy2 mj2 = sqr(j.mass()), mk2 = sqr(k.mass());
  Energy mj = j.mass(), mk = k.mass();
  Energy2 Qjk2 = sjk + mj2 + mk2;
  Energy Qjk = sqrt(Qjk2);
  
  double vjk = rootOfKallen(Qjk2,mj2,mk2) / sjk;
  double rho = sqrt( abs(1.-vjk)/(1.+vjk) ); // abs() because for small mass 1.-vjk can get O(-1.e-16)
  double rhoj = sqrt( ( 1. - vjk + 2.*mj2/Qjk2 / (1.-mj2/Qjk2-mk2/Qjk2) ) /
    ( 1. + vjk + 2.*mj2/Qjk2 / (1.-mj2/Qjk2-mk2/Qjk2) ) );
  double rhok = sqrt( ( 1. - vjk + 2.*mk2/Qjk2 / (1.-mj2/Qjk2-mk2/Qjk2) ) /
    ( 1. + vjk + 2.*mk2/Qjk2 / (1.-mj2/Qjk2-mk2/Qjk2) ) );
  
  //////////////////////////////////////
  // Finite terms of S part (6.20)    //
  // Expanded convention              //
  //////////////////////////////////////
  
  ParticleData l = ( mj2 == ZERO ? k : j );
  
  // both masses zero
  if( mj2 == ZERO && mk2 == ZERO ) {
    res += 0.0;
    // Expanded
    res += 1./2.*log(mu2/sjk)*log(mu2/sjk);
  }

  // one mass zero
  else if( mj2 == ZERO || mk2 == ZERO ) {
    Energy2 m2 = sqr(l.mass());
    res += -1./4.*sqr(log(m2/sjk)) - sqr(pi)/12. -
      1./2.*log(m2/sjk)*log(sjk/Qjk2) - 1./2.*log(m2/Qjk2)*log(sjk/Qjk2);
    // Expanded
    res += 1./2.*log(mu2/sjk)*log(m2/sjk) +
        1./4.*log(mu2/sjk)*log(mu2/sjk);
  }

  // no mass zero
  else if( mj2 != ZERO && mk2 != ZERO ) {
    res += 1./vjk * ( -1./4.*sqr(log(rhoj*rhoj)) - 1./4.*sqr(log(rhok*rhok)) -
      sqr(pi)/6. + ( rho==0. ? 0. : log(rho)*log(Qjk2/sjk) ) );
    // Expanded
    res += 1./vjk * ( rho==0. ? 0. : log(rho)*log(mu2/sjk) );
  }

  //////////////////////////////////////
  // NS part (6.21)-(6.26)            //
  //////////////////////////////////////
  
  // V_q (j is quark)

  // j is massive quark
  if( mj2 != ZERO ) {
    assert( abs(j.id()) < 7);
    res += gammaQuark/CF * log(sjk/Qjk2); // part common iff j is massive quark and k either massive quark or massless parton (6.21),(6.22)

    // k is massive quark (6.21)
    if( mk2 != ZERO ) {
      assert( abs(k.id()) < 7);
      res += 1./vjk * ( ( rho==0. ? 0. : log(rho*rho)*log(1.+rho*rho) ) + 2.*gsl_sf_dilog(rho*rho) -
	gsl_sf_dilog(1.-rhoj*rhoj) - gsl_sf_dilog(1.-rhok*rhok) - sqr(pi)/6. ) +
	log((Qjk-mk)/Qjk) - 2.*log((sqr(Qjk-mk)-mj2)/Qjk2) - 2.*mj2/sjk*log(mj/(Qjk-mk)) -
	mk/(Qjk-mk) + 2.*mk*(2.*mk-Qjk)/sjk + sqr(pi)/2.;
    }

    // k is massless parton (6.22)
    else {
      res += sqr(pi)/6. - gsl_sf_dilog(sjk/Qjk2) -
	2.*log(sjk/Qjk2) - mj2/sjk*log(mj2/Qjk2);
    }

  }

  // V_j (j either massless quark (6.23) or gluon (6.24),(6.26))
  else {

    // k is massless parton
    if( mk == ZERO ) {
      // only contributes if j is gluon (6.26)
      if( j.id() == ParticleID::g ) {
	// sum over all quark flavours
	if( !mFSetEmpty )
	  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) { // only heavy quarks in jet (aka g->QQbar at NLO)
	    Energy2 mF2 = sqr( getParticleData(lastBorn()->nHeavyJetVec()[f])->mass() );
	    // sum only over quarks which meet special condition
	    if( sjk <= 4.*sqrt(mF2)*(sqrt(mF2)+mk) )
	      continue;
	    double rho1 = sqrt( 1. - 4.*mF2 / sqr(Qjk-mk) );
// 	    res += 2./3./CA * ( log((1.+rho1)/2.) - rho1/3.*(3.+sqr(rho1)) - 0.5*log(mF2/Qjk2) );
	    res += 2./3./CA * ( log((1.+rho1)/2.) - rho1/3.*(3.+sqr(rho1)) - 0.5*log(mF2/sjk) );
	  }
	  // The last term with Q_{aux} in (6.26) cancels against a similar term in GammaGluon().
      }
    }

    // k is massive quark
    else {
      assert( abs(k.id()) < 7);
      // part common to j massless quark or gluon (6.23),(6.24)
      res += sqr(pi)/6. - gsl_sf_dilog(sjk/Qjk2);
      // j is massless quark (6.23)
      if( abs(j.id()) < 7)
	res += gammaQuark/CF * ( log(sjk/Qjk2) - 2.*log((Qjk-mk)/Qjk) - 2.*mk/(Qjk+mk) );
      // j is gluon (6.24)
      else if( j.id() == ParticleID::g ) {
	// part independent of other heavy quark flavours
	res += gammaGluon/CA * ( log(sjk/Qjk2) - 2.*log((Qjk-mk)/Qjk) - 2.*mk/(Qjk+mk) ) +
	  (kappa-2./3.) * mk2/sjk * (1./CA*lastBorn()->nLight()-1.) * log(2.*mk/(Qjk+mk));
	// part containing other heavy quark flavours
	if( !mFSetEmpty )
	  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) { // only heavy quarks in jet (aka g->QQbar at NLO)
	    Energy2 mF2 = sqr( getParticleData(lastBorn()->nHeavyJetVec()[f])->mass() );
	    // sum only over quarks which meet special condition
	    if( sjk <= 4.*sqrt(mF2)*(sqrt(mF2)+mk) )
	      continue;
	    double rho1 = sqrt( 1. - 4.*mF2 / sqr(Qjk-mk) );
	    double rho2 = sqrt( 1. - 4.*mF2 / (Qjk2-mk2) );
	    res += 2./3./CA * ( log((Qjk-mk)/Qjk) + mk*rho1*rho1*rho1/(Qjk+mk) + log((1.+rho1)/2.) -
	      rho1/3.*(3.+sqr(rho1)) - 1./2.*log(mF2/Qjk2) ) +
	      1./CA * ( rho2*rho2*rho2*log((rho2-rho1)/(rho2+rho1)) - log((1.-rho1)/(1.+rho1)) -
	      8.*rho1*mF2/sjk ) * (kappa-2./3.) * mk2/sjk;
	  }
	  // The term with Q_{aux} in (6.24) cancels against a similar term in GammaGluon().
      }
    }
  }

  return res;

}

double DipoleMIOperator::VsDoublePole(const ParticleData& j, const ParticleData& k) const {

  double res = 0.;
  
  ////////////////////////////////////////////////
  // Double pole coefficient of S part (6.20)   //
  // Expanded convention                        //
  ////////////////////////////////////////////////
  
  Energy2 mj2 = sqr(j.mass()), mk2 = sqr(k.mass());
  
  // both masses zero
  if( mj2 == ZERO && mk2 == ZERO ) {
    res += 1.0;
  }

  // one mass zero
  else if( mj2 == ZERO || mk2 == ZERO ) {
    res += 1./2.;
  }

  // no mass zero
  else if( mj2 != ZERO && mk2 != ZERO ) {
    res += 0.0;
  }

  return res;

}

double DipoleMIOperator::VsSinglePole(const ParticleData& j, const ParticleData& k, 
                                      Energy2 sjk) const {

  Energy2 mu2 = lastBorn()->mu2();

  double res = 0.;
  
  // sjk is being handed over as input parameter to DipoleMIOperator::VsDoublePole()

  Energy2 mj2 = sqr(j.mass()), mk2 = sqr(k.mass());
  Energy2 Qjk2 = sjk + mj2 + mk2;
  
  double vjk = rootOfKallen(Qjk2,mj2,mk2) / sjk;
  double rho = sqrt( abs(1.-vjk)/(1.+vjk) ); // abs() because for small mass 1.-vjk can get O(-1.e-16)
  
  ////////////////////////////////////////////////
  // Single pole coefficient of S part (6.20)   //
  // Expanded convention                        //
  ////////////////////////////////////////////////
  
  ParticleData l = ( mj2 == ZERO ? k : j );
  
  // both masses zero
  if( mj2 == ZERO && mk2 == ZERO ) {
    res += log(mu2/sjk);
  }

  // one mass zero
  else if( mj2 == ZERO || mk2 == ZERO ) {
    Energy2 m2 = sqr(l.mass());
//     res += 1./2.*(1.0 + log(m2/sjk));
    res += 1./2.*(log(mu2/sjk) + log(m2/sjk));
  }

  // no mass zero
  else if( mj2 != ZERO && mk2 != ZERO ) {
    res += 1./vjk * ( rho==0. ? 0. : log(rho) );
  }

  return res;

}

double DipoleMIOperator::GammaQuark(const ParticleData& j) const {
  if ( j.mass() == ZERO )
    return 0.;
  Energy2 mu2 = lastBorn()->mu2();
  return CF * ( 0.5*log(sqr(j.mass())/mu2) - 2. );
}

double DipoleMIOperator::GammaGluon() const {
  // Finite contribution cancels with similar contribution in VjNS.
  return 0.;
}

double DipoleMIOperator::GammaQuarkSinglePole(const ParticleData& j) const {
  if ( j.mass() == ZERO ) {
    return gammaQuark;
  }
  return CF;
}

double DipoleMIOperator::GammaGluonSinglePole() const {
  return gammaGluon;
}

void DipoleMIOperator::persistentOutput(PersistentOStream & os) const {
  os << CA << CF << gammaQuark << gammaGluon << betaZero 
     << KQuark << KGluon;
}

void DipoleMIOperator::persistentInput(PersistentIStream & is, int) {
  is >> CA >> CF >> gammaQuark >> gammaGluon >> betaZero 
     >> KQuark >> KGluon;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<DipoleMIOperator,MatchboxInsertionOperator>
describeHerwigDipoleMIOperator("Herwig::DipoleMIOperator", "Herwig.so");

void DipoleMIOperator::Init() {

  static ClassDocumentation<DipoleMIOperator> documentation
    ("DipoleMIOperator");

  DipoleRepository::registerInsertionOperator<0,DipoleMIOperator>("MassiveIOperator");

}

