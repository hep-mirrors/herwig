// -*- C++ -*-
//
// PhasespaceHelpers.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "PhasespaceHelpers.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

using namespace Herwig;
using namespace Herwig::PhasespaceHelpers;
using namespace Herwig::RandomHelpers;

Energy PhasespaceInfo::generateMass(tcPDPtr data, 
                                    const pair<Energy,Energy>& range) {

  double xlow = sqr(range.first)/sHat;
  if ( range.first < ZERO )
    xlow = -xlow;

  double xup = sqr(range.second)/sHat;
  if ( range.second < ZERO )
    xup = -xup;

  double mu = sqr(data->hardProcessMass())/sHat;
  double gamma = sqr(data->hardProcessWidth())/sHat;

  if ( gamma < 1e-14 )
    gamma = 0.0;

  if ( M0 != ZERO )
    x0 = M0/sqrtSHat;
  if ( Mc != ZERO )
    xc = Mc/sqrtSHat;

  double r = rnd();

  pair<double,double> event;
  
  if ( gamma == 0. ) {

    if ( mu < xlow || mu > xup ) {
      if ( abs(xlow-mu) < xc )
        xlow = mu;
      if ( abs(xup-mu) < xc )
        xup = mu;
    }

    if ( mu < xlow || mu > xup ) {
      event =
      generate(inverse(mu,xlow,xup),r);
    } else {
      pair<double,double> pLeft(xlow,xlow < mu-x0 ? mu-x0 : xlow);
      pair<double,double> pRight(xup > mu+x0 ? mu+x0 : xup,xup);
      pair<double,double> fLeft(pLeft.second < mu-x0 ? mu-x0 : pLeft.second,mu-xc);
      if ( fLeft.first >= fLeft.second ) fLeft.first = fLeft.second;
      pair<double,double> fRight(mu+xc,pRight.first > mu+x0 ? mu+x0 : pRight.first);
      if ( fRight.first >= fRight.second ) fRight.second = fRight.first;

      
      const bool pL= abs( pLeft.first - pLeft.second ) < 1e-14;
      const bool fL= abs( fLeft.first - fLeft.second ) < 1e-14;
      const bool fR= abs( fRight.first - fRight.second ) < 1e-14;
      const bool pR= abs( pRight.first - pRight.second ) < 1e-14;


      if ( !pL && !fL && !fR && !pR ) {
        event =
        generate((piecewise(),
                  inverse(mu,pLeft.first,pLeft.second),
                  match(flat(fLeft.first,fLeft.second))) +
                 match((piecewise(),
                        flat(fRight.first,fRight.second),
                        match(inverse(mu,pRight.first,pRight.second)))),
                 r);
      } else if ( pL && !fL && !fR && !pR ) {
        event =
        generate(flat(fLeft.first,fLeft.second) +
                 match((piecewise(),
                        flat(fRight.first,fRight.second),
                        match(inverse(mu,pRight.first,pRight.second)))), r);
      } else if ( !pL && !fL && !fR && pR ) {
        event =
        generate((piecewise(),
                  inverse(mu,pLeft.first,pLeft.second),
                  match(flat(fLeft.first,fLeft.second))) +
                 match(flat(fRight.first,fRight.second)),
                 r);
      } else if ( pL && fL && !fR && !pR ) {
        event =
        generate((piecewise(),flat(fRight.first,fRight.second),
                  match(inverse(mu,pRight.first,pRight.second))), r);
      } else if ( !pL && !fL && fR && pR ) {
        event =
        generate((piecewise(),
                  inverse(mu,pLeft.first,pLeft.second),
                  match(flat(fLeft.first,fLeft.second))),r);
      } else if ( pL && !fL && !fR && pR ) {
        event = generate(flat(fLeft.first,fLeft.second) +
                 match(flat(fRight.first,fRight.second)),r);
      } else if ( pL && fL && !fR && pR ) {
        event = generate(flat(fRight.first,fRight.second),r);
      } else if ( pL && !fL && fR && pR ) {
        event = generate(flat(fLeft.first,fLeft.second),r);
      } else if ( pL && fL && fR && pR ) {
        throw Veto();
      } else assert(false);
    }
  } else {
    event = generate(breitWigner(mu,gamma,xlow,xup),r);
  }

  if ( abs(event.first) < xc )
    throw Veto();

  weight *= event.second;

  Energy res = sqrt(abs(event.first)*sHat);
  if ( event.first < 0. )
    res = -res;

  return res;
}

Lorentz5Momentum PhasespaceInfo::generateKt(const Lorentz5Momentum& p1,
                                            const Lorentz5Momentum& p2,
                                            Energy pt) {

  double phi = 2.*Constants::pi*rnd();
  weight *= 2.*Constants::pi;

  Lorentz5Momentum P = p1 + p2;

  Energy2 Q2 = abs(P.m2());

  Lorentz5Momentum Q =  Lorentz5Momentum(ZERO,ZERO,ZERO,sqrt(Q2),sqrt(Q2));

  bool boost =
    abs((P-Q).vect().mag2()/GeV2) > 1e-8 ||
    abs((P-Q).t()/GeV) > 1e-4;
  boost &= (P*Q-Q.mass2())/GeV2 > 1e-8;

  Lorentz5Momentum inFrame1;
  if ( boost )
    inFrame1 = p1 + ((P*p1-Q*p1)/(P*Q-Q.mass2()))*(P-Q);
  else
    inFrame1 = p1;

  Energy ptx = inFrame1.x();
  Energy pty = inFrame1.y();
  Energy q = 2.*inFrame1.z();

  Energy Qp = sqrt(4.*(sqr(ptx)+sqr(pty))+sqr(q));
  Energy Qy = sqrt(4.*sqr(pty)+sqr(q));

  double cPhi = cos(phi);
  double sPhi = sqrt(1.-sqr(cPhi));
  if ( phi > Constants::pi )
    sPhi = -sPhi;

  Lorentz5Momentum kt;

  kt.setT(ZERO);
  kt.setX(pt*Qy*cPhi/Qp);
  kt.setY(-pt*(4*ptx*pty*cPhi/Qp+q*sPhi)/Qy);
  kt.setZ(2.*pt*(-ptx*q*cPhi/Qp + pty*sPhi)/Qy);

  if ( boost )
    kt = kt + ((P*kt-Q*kt)/(P*Q-Q.mass2()))*(P-Q);
  kt.setMass(-pt);
  kt.rescaleRho();

  return kt;

}

void PhasespaceTree::setup(const Tree2toNDiagram& diag, 
                           int pos) {
  doMirror = false;

  pair<int,int> dchildren =  diag.children(pos);

  data = diag.allPartons()[pos];

  spacelike = pos < diag.nSpace();

  if ( pos == 0 )
    externalId = 0;

  if ( dchildren.first == -1 ) {
    externalId = diag.externalId(pos);
    leafs.insert(externalId);
    return;
  }

  children.push_back(PhasespaceTree());
  children.back().setup(diag,dchildren.first);
  children.push_back(PhasespaceTree());
  children.back().setup(diag,dchildren.second);

  if ( !children[0].children.empty() &&
      children[1].children.empty() &&
      !spacelike )
    swap(children[0],children[1]);
  if ( spacelike &&
      !children[0].spacelike )
    swap(children[0],children[1]);

  copy(children[0].leafs.begin(),children[0].leafs.end(),
       inserter(leafs,leafs.begin()));
  copy(children[1].leafs.begin(),children[1].leafs.end(),
       inserter(leafs,leafs.begin()));

}

void PhasespaceTree::setupMirrored(const Tree2toNDiagram& diag, 
                                   int pos) {

  doMirror = true;

  spacelike = pos < diag.nSpace();

  pair<int,int> dchildren;
  if (pos != 0 && spacelike)
    dchildren = {diag.parent(pos), diag.children(diag.parent(pos)).second};
  else if ( !spacelike ) dchildren = diag.children(pos);
  else                   dchildren = {-1,-1};

  data = diag.allPartons()[pos];

  if ( pos == diag.nSpace() - 1 )
    externalId = 1;

  if ( dchildren.first == -1 ) {
    externalId = diag.externalId(pos);
    leafs.insert(externalId);
    return;
  }
  


  children.push_back(PhasespaceTree());
  children.back().setupMirrored(diag,dchildren.first);
  children.push_back(PhasespaceTree());
  children.back().setupMirrored(diag,dchildren.second);

  if ( !children[0].children.empty() &&
      children[1].children.empty() &&
      !spacelike )
    swap(children[0],children[1]);
  if ( spacelike &&
      !children[0].spacelike ) {
    assert (false);
  }

  copy(children[0].leafs.begin(),children[0].leafs.end(),
       inserter(leafs,leafs.begin()));
  copy(children[1].leafs.begin(),children[1].leafs.end(),
       inserter(leafs,leafs.begin()));

}

void PhasespaceTree::print(int in) {
  for (int i = 0; i != in; i++)
    cerr << "   ";
  cerr << " |- "  << data->PDGName() << " " << externalId << "\n" << flush;
  if ( !children.empty() ) {
    children[1].print(in+1);
    children[0].print(in+int(!spacelike));
  }
  else {
    cerr << "\n";
  }
}

void PhasespaceTree::init(const vector<Lorentz5Momentum>& meMomenta) {

  if ( children.empty() ) {
    massRange.first = meMomenta[externalId].mass();
    massRange.second = massRange.first;
    if ( !doMirror && externalId == 1 )
      momentum = meMomenta[1];
    if ( doMirror && externalId == 0 )
      momentum = meMomenta[0];
    momentum.setMass(meMomenta[externalId].mass());
    return;
  }

  children[0].init(meMomenta);
  children[1].init(meMomenta);

  if ( !children[0].spacelike &&
      !children[1].spacelike ) {
    massRange.first =
    children[0].massRange.first +
    children[1].massRange.first;
  }

}

void PhasespaceTree::generateKinematics(PhasespaceInfo& info,
                                        vector<Lorentz5Momentum>& meMomenta) {
  
  if ( !doMirror && externalId == 0 ) {
    init(meMomenta);
    Energy2 s = (meMomenta[0]+meMomenta[1]).m2();
    double sign = meMomenta[0].z() >= ZERO ? 1. : -1;
    momentum = Lorentz5Momentum(ZERO,ZERO,sign*sqrt(s)/2.,sqrt(s)/2.,ZERO);
    backwardMomentum = Lorentz5Momentum(ZERO,ZERO,-sign*sqrt(s)/2.,sqrt(s)/2.,ZERO);
  }
  else if ( doMirror && externalId == 1) {
    init(meMomenta);
    Energy2 s = (meMomenta[0]+meMomenta[1]).m2();
    double sign = meMomenta[0].z() >= ZERO ? 1. : -1;
    momentum = Lorentz5Momentum(ZERO,ZERO,-sign*sqrt(s)/2.,sqrt(s)/2.,ZERO);
    backwardMomentum = Lorentz5Momentum(ZERO,ZERO,sign*sqrt(s)/2.,sqrt(s)/2.,ZERO);
  }

  if ( children.empty() ) {
    if ( ( !doMirror && externalId != 1 )
        || ( doMirror && externalId !=0 ) )
      meMomenta[externalId] = momentum;
    return;
  }

    // s-channel
  if ( ( !doMirror && externalId == 0 &&
          children[0].externalId == 1 )
     || ( doMirror && externalId == 1 &&
          children[0].externalId == 0 ) ) {
        children[1].momentum = meMomenta[0] + meMomenta[1];
        children[1].momentum.setMass(info.sqrtSHat);
        children[1].momentum.rescaleEnergy();
        children[1].generateKinematics(info,meMomenta);
        return;
      }

  if ( !spacelike ) {

    Energy mij = momentum.mass();
    Energy mi,mj;
    
      // work out the mass for the first child
    if ( !children[0].children.empty() ) {
      Energy sumOthers = ZERO;
      for ( size_t k = 2; k < meMomenta.size(); ++k )
        if ( children[1].leafs.find(k) != children[1].leafs.end() )
          sumOthers += meMomenta[k].mass();
      children[0].massRange.second = momentum.mass() - sumOthers;
      if ( children[0].massRange.second < children[0].massRange.first )
        throw Veto();
      if ( children[0].massRange.second > momentum.mass() )
        throw Veto();
      mi = info.generateMass(children[0].data,children[0].massRange);
      children[0].momentum.setMass(mi);
    } else {
      mi = children[0].momentum.mass();
    }

      // work out the mass for the second child
    if ( !children[1].children.empty() ) {
      children[1].massRange.second = momentum.mass()-children[0].momentum.mass();
      if ( children[1].massRange.second < children[1].massRange.first )
        throw Veto();
      mj = info.generateMass(children[1].data,children[1].massRange);
      children[1].momentum.setMass(mj);
    } else {
      mj = children[1].momentum.mass();
    }

    Energy2 mij2 = sqr(mij);
    Energy2 mi2 = sqr(mi);
    Energy2 mj2 = sqr(mj);

      // perform the decay
    Energy4 lambda2 = sqr(mij2-mi2-mj2)-4.*mi2*mj2;
    if ( lambda2 <= ZERO )
      throw Veto();
    Energy2 lambda = sqrt(lambda2);
    double phi = 2.*Constants::pi*info.rnd();
    double cosPhi = cos(phi);
    double sinPhi = sqrt(1.-sqr(cosPhi));
    if ( phi > Constants::pi )
      sinPhi = -sinPhi;
    info.weight *= Constants::pi*lambda/(2.*mij2);
    double cosTheta = 2.*info.rnd() - 1.;
    double sinTheta = sqrt(1.-sqr(cosTheta));
    Energy p = lambda/(2.*mij);
    children[0].momentum.setX(p*cosPhi*sinTheta);
    children[0].momentum.setY(p*sinPhi*sinTheta);
    children[0].momentum.setZ(p*cosTheta);
    children[0].momentum.rescaleEnergy();
    if ( momentum.m2() <= ZERO )
      throw Veto();
    Boost out = momentum.boostVector();
    if ( out.mag2() > Constants::epsilon ) {
      children[0].momentum.boost(out);
    }
    children[1].momentum = momentum - children[0].momentum;
    children[1].momentum.setMass(mj);
    children[1].momentum.rescaleEnergy();

      // go on with next branchings
    children[0].generateKinematics(info,meMomenta);
    children[1].generateKinematics(info,meMomenta);

    return;

  }

    // get the minimum mass of the `W' system
  Energy Wmin = ZERO;
  PhasespaceTree* current = &children[0];
  while ( !(current->children.empty()) ) {
    Wmin += current->children[1].massRange.first;
    current = &(current->children[0]);
  }

    // get the CM energy avaialble
  Energy2 s = (momentum+backwardMomentum).m2();
  if ( s <= ZERO )
    throw Veto();

    // generate a mass for the timelike child
  Energy mi;
  if ( !children[1].children.empty() ) {
    children[1].massRange.second = sqrt(s)-Wmin;
    if ( children[1].massRange.second < children[1].massRange.first )
      throw Veto();
    mi = info.generateMass(children[1].data,children[1].massRange);
    children[1].momentum.setMass(mi);
  } else {
    mi = children[1].momentum.mass();
  }
  Energy2 mi2 = sqr(mi);
 
    // wether or not this is the last 2->2 scatter
  bool lastScatter = children[0].children[0].children.empty();

    // `W' mass relevant for the other boundaries
  Energy MW = Wmin;

    // generate a mass for second outgoing leg, if needed
  if ( lastScatter )
    if ( !children[0].children[1].children.empty() ) {
        // get the maximum `W' mass
      Energy Wmax = sqrt(s)-children[1].momentum.mass();
      children[0].children[1].massRange.second = Wmax;
      if ( children[0].children[1].massRange.second <
          children[0].children[1].massRange.first )
        throw Veto();
      MW = info.generateMass(children[0].children[1].data,
                             children[0].children[1].massRange);
      children[0].children[1].momentum.setMass(MW);
    }
  Energy2 MW2 = sqr(MW);

  Energy ma = momentum.mass();
  Energy2 ma2 = sqr(ma);
  if ( ma < ZERO )
    ma2 = -ma2;
  Energy mb = backwardMomentum.mass();
  Energy2 mb2 = sqr(mb);
  if ( mb < ZERO )
    mb2 = -mb2;

    // pick the ys variable
  Energy2 ys = ZERO;
  if ( !lastScatter ) {
    ys = info.rnd()*(sqr(sqrt(s)-mi)-MW2);
    info.weight *= (sqr(sqrt(s)-mi)-MW2)/info.sHat;
  }

  Energy4 lambda2 = sqr(s-ma2-mb2)-4.*ma2*mb2;
  if ( lambda2 <= ZERO ) {
    throw Veto();
  }
  Energy2 lambda = sqrt(lambda2);
  info.weight *= info.sHat/(4.*lambda);

    // get the boundaries on the momentum transfer
  Energy4 rho2 = sqr(s-ys-MW2-mi2)-4.*mi2*(ys+MW2);
  if ( rho2 < ZERO )
    throw Veto();
  Energy2 rho = sqrt(rho2);
  Energy4 tau2 =
  ys*(ma2-mb2+s)
  - sqr(s)+s*(ma2+mb2+mi2+MW2)-(mi2-MW2)*(ma2-mb2);
  pair<Energy2,Energy2> tBounds
  ((tau2-rho*lambda)/(2.*s),(tau2+rho*lambda)/(2.*s));
  children[0].massRange.first = sqrt(abs(tBounds.first));
  if ( tBounds.first < ZERO )
    children[0].massRange.first = -children[0].massRange.first;
  children[0].massRange.second = sqrt(abs(tBounds.second));
  if ( tBounds.second < ZERO )
    children[0].massRange.second = -children[0].massRange.second;

    // generate a momentum transfer
  Energy mai = info.generateMass(children[0].data,children[0].massRange);
  children[0].momentum.setMass(mai);
  Energy2 t = sqr(mai);
  if ( mai < ZERO )
    t = -t;

  Energy2 u = -s -t + ys + ma2 + mb2 + mi2 + MW2;
  Energy2 st = s - ma2 - mb2;
  Energy2 tt = t - mi2 - ma2;
  Energy2 ut = u - mi2 - mb2;

    // get the timelike momentum
  double xa = (-st*ut+2.*mb2*tt)/lambda2;
  double xb = (-st*tt+2.*ma2*ut)/lambda2;
  Energy2 pt2 = (st*tt*ut-ma2*sqr(ut)-mb2*sqr(tt)-mi2*sqr(st)+4.*ma2*mb2*mi2)/lambda2;
  if ( pt2 < ZERO )
    throw Veto();
  Energy pt = sqrt(pt2);

  children[1].momentum =
  xa*momentum + xb*backwardMomentum 
  + info.generateKt(momentum,backwardMomentum,pt);
  children[1].momentum.setMass(mi);
  children[1].momentum.rescaleEnergy();

  children[0].momentum =
  momentum - children[1].momentum;
  children[0].momentum.setMass(mai);
  bool changeSign = false;
  if ( children[0].momentum.t() < ZERO && mai > ZERO ) changeSign = true;
  if ( mai < ZERO )
    children[0].momentum.rescaleRho();
  else
    children[0].momentum.rescaleEnergy();
  if ( changeSign ) children[0].momentum.setT(-children[0].momentum.t());
 
  children[1].generateKinematics(info,meMomenta);

  if ( !lastScatter ) {
    children[0].backwardMomentum = backwardMomentum;
    children[0].generateKinematics(info,meMomenta);
  } else {
    children[0].children[1].momentum =
    backwardMomentum + children[0].momentum;
    children[0].children[1].momentum.setMass(MW);
    children[0].children[1].momentum.rescaleEnergy();
    children[0].children[1].generateKinematics(info,meMomenta);
  }

}


void PhasespaceTree::put(PersistentOStream& os) const {
  os << children.size();
  if ( !children.empty() ) {
    children[0].put(os);
    children[1].put(os);
  }
  os << data << externalId << leafs << spacelike << doMirror;
}

void PhasespaceTree::get(PersistentIStream& is) {
  size_t nc; is >> nc;
  assert(nc == 0 || nc == 2);
  if ( nc == 2 ) {
    children.resize(2,PhasespaceTree());
    children[0].get(is);
    children[1].get(is);
  }
  is >> data >> externalId >> leafs >> spacelike >> doMirror;
}
