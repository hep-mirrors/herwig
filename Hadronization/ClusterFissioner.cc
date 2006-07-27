// -*- C++ -*-
//
// Thisk is the implementation of the non-inlined, non-templated member
// functions of the ClusterFissioner class.
//

#include "ClusterFissioner.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/EventRecord/Collision.h>
#include "Herwig++/Utilities/Kinematics.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Herwig++/Utilities/CheckId.h"
#include "Cluster.h"
#include <iomanip>

using namespace Herwig;

void ClusterFissioner::persistentOutput(PersistentOStream & os) const {
  os << _hadronsSelector << _globalParameters
     << _clMax << _clPow << _pSplit1 << _pSplit2 << _btClM << _iopRem;
}

void ClusterFissioner::persistentInput(PersistentIStream & is, int) {
  is >> _hadronsSelector >> _globalParameters
     >> _clMax >> _clPow >> _pSplit1 >> _pSplit2 >> _btClM >> _iopRem;
}

ClassDescription<ClusterFissioner> ClusterFissioner::initClusterFissioner;
// Definition of the static class description member.


void ClusterFissioner::Init() {

  static ClassDocumentation<ClusterFissioner> documentation
    ("Class responsibles for chopping up the clusters");

  static Reference<ClusterFissioner,HadronSelector> 
    interfaceHadronSelector("HadronSelector", 
                             "A reference to the HadronSelector object", 
                             &Herwig::ClusterFissioner::_hadronsSelector,
			     false, false, true, false);

  static Reference<ClusterFissioner,GlobalParameters> 
    interfaceGlobalParameters("GlobalParameters", 
			      "A reference to the GlobalParameters object", 
			      &Herwig::ClusterFissioner::_globalParameters,
			      false, false, true, false);
  
  static Parameter<ClusterFissioner,Energy>
    interfaceClMax ("ClMax","cluster max mass  (unit [GeV])",
                    &ClusterFissioner::_clMax, GeV, 3.35*GeV, 0.0*GeV, 10.0*GeV,
		    false,false,false);

  static Parameter<ClusterFissioner,double>
    interfaceClPow ("ClPow","cluster mass exponent",
                    &ClusterFissioner::_clPow, 0, 2.0, 0.0, 10.0,false,false,false);

  static Parameter<ClusterFissioner,double>
    interfacePSplt1 ("PSplt1","cluster mass splitting param for u,d,s,c",
                    &ClusterFissioner::_pSplit1, 0, 1.0, 0.0, 10.0,false,false,false);

  static Parameter<ClusterFissioner,double>
    interfacePSplt2 ("PSplt2","cluster mass splitting param for b",
                    &ClusterFissioner::_pSplit2, 0, 1.0, 0.0, 10.0,false,false,false);

  static Switch<ClusterFissioner,int> interfaceRemnantOption
    ("RemnantOption",
     "Option for the treatment of remnant clusters",
     &ClusterFissioner::_iopRem, 1, false, false);
  static SwitchOption interfaceRemnantOptionSoft
    (interfaceRemnantOption,
     "Soft",
     "Both clusters produced in the fission of the beam cluster"
     " are treated as soft clusters.",
     0);
  static SwitchOption interfaceRemnantOptionHard
    (interfaceRemnantOption,
     "Hard",
     "Only the cluster containing the remnant is treated as a soft cluster.",
     1);

  static Parameter<ClusterFissioner,Energy> interfaceBTCLM
    ("BTCLM",
     "Parameter for the mass spectrum of remnant clusters",
     &ClusterFissioner::_btClM, GeV, 1.*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

void ClusterFissioner::fission(const StepPtr &pstep) {
  /*****************
   * Loop over the (input) collection of cluster pointers, and store in 
   * the vector  splitClusters  all the clusters that need to be split
   * (these are beam clusters, if soft underlying event is off, and 
   *  heavy non-beam clusters).
   ********************/
  vector<tClusterPtr> splitClusters; 
  
  ClusterVector clusters; 
  for (ParticleSet::iterator it = pstep->particles().begin();
       it!= pstep->particles().end(); it++) { 
    if((*it)->id() == ExtraParticleID::Cluster) 
      clusters.push_back(dynamic_ptr_cast<ClusterPtr>(*it));
  }
  
  for(ClusterVector::iterator it = clusters.begin() ; 
      it != clusters.end() ; ++it) {
    
    /**************
     * Skip 3-component clusters that have been redefined (as 2-component 
     * clusters) or not available clusters. The latter check is indeed 
     * redundant now, but it is used for possible future extensions in which, 
     * for some reasons, some of the clusters found by ClusterFinder are tagged
     * straight away as not available.
     **************/
    if((*it)->isRedefined() || !(*it)->isAvailable()) continue;
    // if the cluster is a beam cluster and Underlying event is off
    // add it to the vector of clusters to be split
    if((*it)->isBeamCluster()) {
      if (_globalParameters->isSoftUnderlyingEventON()) {
	(*it)->isAvailable(false);
      } else {
	splitClusters.push_back(*it);
      }
      continue;
    }

    // If the cluster is heavy add it to the vector of clusters to be split.
    if(pow((*it)->mass() , _clPow) > 
       pow(_clMax, _clPow) + pow((*it)->sumConstituentMasses(), _clPow))
      splitClusters.push_back(*it);
  }
  // split the clusters
  vector<tClusterPtr>::const_iterator iter;
  for(iter = splitClusters.begin(); 
      iter != splitClusters.end() ; 
      ++iter) {
    cut(*iter, pstep, clusters);
  }
}

void ClusterFissioner::cut(tClusterPtr cluster, const StepPtr &pstep, 
   			   ClusterVector &clusters) {
  /**************************************************
   * This method does the splitting of the cluster pointed by  cluPtr
   * and "recursively" by all of its cluster children, if heavy. All of these
   * new children clusters are added (indeed the pointers to them) to the
   * collection of cluster pointers  collecCluPtr. The method works as follows.
   * Initially the vector vecCluPtr contains just the input pointer to the
   * cluster to be split. Then it will be filled "recursively" by all
   * of the cluster's children that are heavy enough to require, in their turn,
   * to be split. In each loop, the last element of the vector vecCluPtr is 
   * considered (only once because it is then removed from the vector).
   * This approach is conceptually recursive, but avoid the overhead of
   * a concrete recursive function. Furthermore it requires minimal changes
   * in the case that the fission of an heavy cluster could produce more
   * than two cluster children as assumed now. 
   *
   * Draw the masses: for normal, non-beam clusters a power-like mass dist
   * is used, whereas for beam clusters a fast-decreasing exponential mass 
   * dist is used instead (to avoid many iterative splitting which could 
   * produce an unphysical large transverse energy from a supposed soft beam
   * remnant process).
   ****************************************/
  vector<tClusterPtr> clusterStack; 
  clusterStack.push_back(cluster);
  // Here we recursively loop over clusters in the stack and cut them
  while (!clusterStack.empty()) {
    // take the last element of the vector
    tClusterPtr iCluster = clusterStack.back(); 
    clusterStack.pop_back();   
    // split it                 
    cutType ct = cut(iCluster);
    // There are cases when we don't want to split, even if it fails mass test
    if(!ct.first.first || !ct.second.first)
      {
	// if an unsplit beam cluster leave if for the underlying event
// 	if(iCluster->isBeamCluster()
// 	   &&_globalParameters->isSoftUnderlyingEventON()) {
// 	  iCluster->isAvailable(false);
// 	}	
	continue;
      }
    // check if clusters
    ClusterPtr one = dynamic_ptr_cast<ClusterPtr>(ct.first.first);
    ClusterPtr two = dynamic_ptr_cast<ClusterPtr>(ct.second.first);
    // is a beam cluster must be split into two hadrons
//     if(iCluster->isBeamCluster()&&(!one||!two)
//        &&_globalParameters->isSoftUnderlyingEventON())
//       {
// 	iCluster->isAvailable(false);
// 	continue;
//       }

    // There should always be a intermediate quark(s) from the splitting, but
    // in case there isn't
    if(ct.first.second) {
      pstep->addIntermediate(ct.first.second);
      ct.first.second->addChild(ct.first.first);
    }
    if(ct.second.second) {
      pstep->addIntermediate(ct.second.second);
      ct.second.second->addChild(ct.second.first);
    }
    pstep->addDecayProduct(iCluster, ct.first.first);
    pstep->addDecayProduct(iCluster, ct.second.first);

    // Sometimes the clusters decay C -> H + C' rather then C -> C' + C''
    if(one) {
      clusters.push_back(one);
//       if(one->isBeamCluster()&&
// 	 _globalParameters->isSoftUnderlyingEventON()) one->isAvailable(false);
      if(pow(one->mass(), _clPow) > 
	 pow(_clMax, _clPow) + pow(one->sumConstituentMasses(), _clPow)
	 &&one->isAvailable()) {
	clusterStack.push_back(one);
      } 
    }
    if(two) {
      clusters.push_back(two);
//       if(two->isBeamCluster()&&
// 	 _globalParameters->isSoftUnderlyingEventON()) two->isAvailable(false);
      if(pow(two->mass(), _clPow) > 
	 pow(_clMax, _clPow) + pow(two->sumConstituentMasses(), _clPow)
	 && two->isAvailable()) {
	clusterStack.push_back(two);
      } 
    }
  }
}

ClusterFissioner::cutType ClusterFissioner::cut(tClusterPtr &cluster) {
  // Get the actual particles making up the cluster
  long idQ1 = 0, idQ2 = 0;
  tPPtr ptrQ1 = cluster->particle(0), ptrQ2 = cluster->particle(1);
  if(ptrQ1) idQ1 = ptrQ1->id();
  if(ptrQ2) idQ2 = ptrQ2->id(); 
  
  // And check if those particles are from a beam remnant
  bool rem1 = cluster->isBeamRemnant(0);
  bool rem2 = cluster->isBeamRemnant(1);
  // workout which distribution to use
  bool soft1=rem1||(_iopRem==0&&rem2);
  bool soft2=rem2||(_iopRem==0&&rem1);
  // Initialization for the exponential ("soft") mass distribution.
  static const InvEnergy b = 2.0 / _btClM;
  static const int max_loop = 1000;
  int counter = 0;
  bool succeeded=false;
  Energy Mc1 = Energy(), Mc2 = Energy(),m1=Energy(),m2=Energy(),m=Energy();
  bool toHadron1,toHadron2;
  long idNew;
  do
    {
      succeeded=false;
      ++counter;
      // Draw new flavour
      idNew = drawNewFlavour();      // draw the new flavour (idNew > 0)
      if(!CheckId::canBeMeson(idQ1,-idNew) && !CheckId::canBeBaryon(idQ1,-idNew))
	idNew = -idNew;
      // Check that new clusters can produce particles and there is enough
      // phase space to choose the drawn flavour
      Energy Mc = cluster->mass(); 
      m1 = ptrQ1->data().constituentMass();
      m2 = ptrQ2->data().constituentMass();
      m  = getParticleData(abs(idNew))->constituentMass();
      // Do not split in the case there is no phase space available
      // (it happens sometimes for clusters with b-flavour)
      if(Mc <  m1+m + m2+m) continue;
      // power for splitting
      double exp1=_pSplit1, exp2=_pSplit1;
      if(CheckId::hasBeauty(idQ1)) exp1 = _pSplit2;
      if(CheckId::hasBeauty(idQ2)) exp2 = _pSplit2;
      // If, during the drawing of candidate masses, too many attempts fail 
      // (because the phase space available is tiny) 
      //_hadronsSelector->lightestHadron(idQ2, idNew)then give up (the cluster 
      // is not split).
      Mc1 = Energy();
      Mc2 = Energy();
      drawChildMass(Mc,m1,m2,m,Mc1,exp1,b,soft1);
      drawChildMass(Mc,m2,m1,m,Mc2,exp2,b,soft2);
      if(Mc1<m1+m || Mc2<m+m2 || Mc1+Mc2>Mc) continue;
      /**************************
       * New (not present in Fortran Herwig):
       * check whether the fragment masses  Mc1  and  Mc2  are above the 
       * threshold for the production of the lightest pair of hadrons with the 
       * right flavours. If not, then set by hand the mass to the lightest 
       * single hadron with the right flavours, in order to solve correctly
       * the kinematics, and (later in this method) create directly such hadron
       * and add it to the children hadrons of the cluster that undergoes the
       * fission (i.e. the one pointed by iCluPtr). Notice that in this special
       * case, the heavy cluster that undergoes the fission has one single 
       * cluster child and one single hadron child. We prefer this approach,
       * rather than to create a light cluster, with the mass set equal to
       * the lightest hadron, and let then the class LightClusterDecayer to do 
       * the job to decay it to that single hadron, for two reasons: 
       * First, because the sum of the masses of the two constituents can be, 
       * in this case, greater than the mass of that hadron, hence it would
       * be impossible to solve the kinematics for such two components, and
       * therefore we would have a cluster whose components are undefined.
       * Second, the algorithm is faster, because it avoids the reshuffling
       * procedure that would be necessary if we used LightClusterDecayer
       * to decay the light cluster to the lightest hadron.   
       ****************************/
      toHadron1 = false;
      if(Mc1 < _hadronsSelector->massLightestHadronPair(idQ1,-idNew)) { 
	Mc1 = _hadronsSelector->massLightestHadron(idQ1,-idNew);          
	toHadron1 = true;
      }
      
      toHadron2 = false;
      if(Mc2 < _hadronsSelector->massLightestHadronPair(idQ2,idNew)) { 
	Mc2 = _hadronsSelector->massLightestHadron(idQ2,idNew);           
	toHadron2 = true;
      }

      // Check if the decay kinematics is still possible: if not then 
      // force the one-hadron decay for the other cluster as well.
      if(Mc1 + Mc2  >  Mc) {
	if(!toHadron1) {
	  Mc1 = _hadronsSelector->massLightestHadron(idQ1, -idNew); 
	  toHadron1 = true;	
	} else if(!toHadron2) {
	  Mc2 = _hadronsSelector->massLightestHadron(idQ2, idNew);
	  toHadron2 = true;
	}
      }
      succeeded = (Mc>=Mc1+Mc2); 
    }
  while (!succeeded && counter < max_loop);

  if(counter >= max_loop)
    return cutType(PPair(PPtr(),PPtr()),PPair(PPtr(),PPtr()));



  // Determined the (5-components) momenta (all in the LAB frame)
  Lorentz5Momentum pClu = cluster->momentum(); // known
  Lorentz5Momentum p0Q1 = ptrQ1->momentum(); // known (mom Q1 before fission)
  Lorentz5Momentum pClu1, pClu2, pQ1, pQone, pQtwo, pQ2; //unknown
  pClu1.setMass(Mc1);
  pClu2.setMass(Mc2);
  pQ1.setMass(m1);
  pQ2.setMass(m2);
  pQone.setMass(m); 
  pQtwo.setMass(m);
    
  calculateKinematics(pClu,p0Q1,toHadron1,toHadron2,
		      pClu1,pClu2,pQ1,pQone,pQtwo,pQ2);                // out

  /******************
   * The previous methods have determined the kinematics and positions
   * of C -> C1 + C2. 
   * In the case that one of the two product is light, that means either
   * decayOneHadronClu1 or decayOneHadronClu2 is true, then the momenta
   * of the components of that light product have not been determined,
   * and a (light) cluster will not be created: the heavy father cluster
   * decays, in this case, into a single (not-light) cluster and a
   * single hadron. In the other, "normal", cases the father cluster
   * decays into two clusters, each of which has well defined components.
   * Notice that, in the case of components which point to particles, the
   * momenta of the components is properly set to the new values, whereas
   * we do not change the momenta of the pointed particles, because we 
   * want to keep all of the information (that is the new momentum of a
   * component after the splitting, which is contained in the _momentum
   * member of the Component class, and the (old) momentum of that component
   * before the splitting, which is contained in the momentum of the
   * pointed particle). Please not make confusion of this only apparent
   * inconsistency!
   ********************/
  LorentzPoint posC,pos1,pos2;
  posC = cluster->labVertex();
  calculatePositions(pClu,posC,pClu1,pClu2,pos1,pos2);
  cutType rval;
  if(toHadron1) rval.first = produceHadron(idQ1,-idNew,pClu1,pos1);
  else rval.first = produceCluster(ptrQ1,-idNew, pClu1, pos1,
				   pQ1, pQone, rem1);

  if(toHadron2) rval.second = produceHadron(idQ2,idNew,pClu2,pos2);
  else rval.second = produceCluster(ptrQ2, idNew, pClu2, pos2,
				    pQ2, pQtwo, rem2);
  return rval;
}

ClusterFissioner::PPair ClusterFissioner::produceHadron(long id1, long id2, 
		                                        Lorentz5Momentum &a,
				                        LorentzPoint &b) const {
  PPair rval;
  rval.first = getParticle(_hadronsSelector->lightestHadron(id1, id2));
  rval.second = getParticle(id2);
  rval.first->set5Momentum(a);
  rval.first->setLabVertex(b);
  return rval;
}

ClusterFissioner::PPair ClusterFissioner::produceCluster(tPPtr &p1, long id, 
		                                         Lorentz5Momentum &a,
				                         LorentzPoint &b, 
							 Lorentz5Momentum &c,
				                         Lorentz5Momentum &d, 
							 bool isRem) const {
  PPair rval;
  rval.second = getParticle(id);
  ClusterPtr cluster = new_ptr(Cluster(p1,rval.second));
  rval.first = cluster;
  rval.first->set5Momentum(a);
  rval.first->setLabVertex(b);
  if(cluster->particle(0)->id() == p1->id()) {
    cluster->particle(0)->set5Momentum(c);
    cluster->particle(1)->set5Momentum(d);
    cluster->setBeamRemnant(0,isRem);
  } else {
    cluster->particle(0)->set5Momentum(d);
    cluster->particle(1)->set5Momentum(c);
    cluster->setBeamRemnant(1,isRem);
  }
  return rval;
}

long ClusterFissioner::drawNewFlavour() const {

  // Flavour is assumed to be only  u, d, s,  with weights
  // (which are not normalized probabilities) given
  // by the same weights as used in HadronsSelector for
  // the decay of clusters into two hadrons. 
  double prob_d = _hadronsSelector->pwtDquark();
  double prob_u = _hadronsSelector->pwtUquark();
  double prob_s = _hadronsSelector->pwtSquark();
  double sum = prob_d + prob_u + prob_s;
  prob_d = prob_d / sum;
  prob_u = prob_u / sum;
  prob_s = prob_s / sum;
  int choice = rnd3(prob_u, prob_d, prob_s);
  long idNew = 0;
  switch (choice) {
  case 0: idNew = ThePEG::ParticleID::u; break;  
  case 1: idNew = ThePEG::ParticleID::d; break;
  case 2: idNew = ThePEG::ParticleID::s; break;
  }
  return idNew;
}


void ClusterFissioner::drawChildMass(const Energy M, const Energy m1,
				     const Energy m2, const Energy m,
                                     Energy & Mclu, const double expt, 
                                     const InvEnergy b, const bool soft) const {

  /***************************
   * This method, given in input the cluster mass Mclu of an heavy cluster C,
   * made of consituents of masses m1 and m2, draws the masses Mclu1 and Mclu2
   * of, respectively, the children cluster C1, made of constituent masses m1 
   * and m, and cluster C2, of mass Mclu2 and made of constituent masses m2 
   * and m. The mass is extracted from one of the two following mass 
   * distributions:
   *   --- power-like ("normal" distribution)
   *                        d(Prob) / d(M^exponent) = const
   *       where the exponent can be different from the two children C1 (exp1)
   *       and C2 (exponent2).
   *   --- exponential ("soft" distribution) 
   *                        d(Prob) / d(M^2) = exp(-b*M) 
   *       where b = 2.0 / average.
   * Such distributions are limited below by the masses of
   * the constituents quarks, and above from the mass of decaying cluster C. 
   * The choice of which of the two mass distributions to use for each of the
   * two cluster children is dictated by  iRemnant  (see below).
   * If the number of attempts to extract a pair of mass values that are 
   * kinematically acceptable is above some fixed number (max_loop, see below)
   * the method gives up and returns false; otherwise, when it succeeds, it
   * returns true. 
   *
   * These distributions have been modified from HERWIG:   
   * Before these were:
   *      Mclu1 = m1 + (Mclu - m1 - m2)*pow( rnd(), 1.0/exponent1 );
   * The new one coded here is a more efficient version, same density 
   * but taking into account 'in phase space from' beforehand      
   ***************************/
  Energy max = M-m1-m2-2.0*m;
  double rmin = exp(-b*max);
  if(rmin > 50.0) rmin = 0.0;

  // hard cluster
  if(!soft)
    {Mclu = pow(rnd(pow(M-m1-m2-m, expt), pow(m, expt)), 1./expt) + m1;}
  // Otherwise it uses a soft mass distribution
  else 
    { 
      double r1 = rnd(rmin, 1.0-rmin) * rnd(rmin, 1.0-rmin);
      if(r1 > rmin) Mclu = m1 + m - log(r1)/b;
      else Mclu = Energy();
    }
}


void ClusterFissioner::calculateKinematics(const Lorentz5Momentum & pClu,
					   const Lorentz5Momentum & p0Q1,
					   const bool toHadron1,
					   const bool toHadron2,
					   Lorentz5Momentum & pClu1, 
					   Lorentz5Momentum & pClu2, 
					   Lorentz5Momentum & pQ1,
					   Lorentz5Momentum & pQbar,
					   Lorentz5Momentum & pQ,     
					   Lorentz5Momentum & pQ2bar) const {

  /******************
   * This method solves the kinematics of the two body cluster decay:
   *    C (Q1 Q2bar)  --->  C1 (Q1 Qbar)  +  C2 (Q Q2bar)
   * In input we receive the momentum of C, pClu, and the momentum 
   * of the quark Q1 (constituent of C), p0Q1, both in the LAB frame.
   * Furthermore, two boolean variables inform whether the two fission
   * products (C1, C2) decay immediately into a single hadron (in which
   * case the cluster itself is identify with that hadron) and we do 
   * not have to solve the kinematics of the components (Q1,Qbar) for
   * C1 and (Q,Q2bar) for C2.
   * The output is given by the following momenta (all 5-components, 
   * and all in the LAB frame):
   *   pClu1 , pClu2   respectively of   C1 , C2  
   *   pQ1 , pQbar     respectively of   Q1 , Qbar  in  C1
   *   pQ  , pQ2bar    respectively of   Q  , Q2    in  C2
   * The assumption, suggested from the string model, is that, in C frame,
   * C1 and its constituents Q1 and Qbar are collinear, and collinear to
   * the direction of Q1 in C (that is before cluster decay); similarly,
   * (always in the C frame) C2 and its constituents Q and Q2bar are
   * collinear (and therefore anti-collinear with C1,Q1,Qbar).
   * The solution is then obtained by using Lorentz boosts, as follows.
   * The kinematics of C1 and C2 is solved in their parent C frame,
   * and then boosted back in the LAB. The kinematics of Q1 and Qbar
   * is solved in their parent C1 frame and then boosted back in the LAB; 
   * similarly, the kinematics of Q and Q2bar is solved in their parent 
   * C2 frame and then boosted back in the LAB. In each of the three
   * "two-body decay"-like cases, we use the fact that the direction 
   * of the motion of the decay products is known in the rest frame of 
   * their parent. This is obvious for the first case in which the 
   * parent rest frame is C; but it is also true in the other two cases
   * where the rest frames are C1 and C2. This is because C1 and C2 
   * are boosted w.r.t. C in the same direction where their components, 
   * respectively (Q1,Qbar) and (Q,Q2bar) move in C1 and C2 rest frame
   * respectively.      
   * Of course, although the notation used assumed that C = (Q1 Q2bar)
   * where Q1 is a quark and Q2bar an antiquark, indeed everything remain
   * unchanged also in all following cases:
   *  Q1 quark, Q2bar antiquark;           --> Q quark;
   *  Q1 antiquark , Q2bar quark;          --> Q antiquark;  
   *  Q1 quark, Q2bar diquark;             --> Q quark
   *  Q1 antiquark, Q2bar anti-diquark;    --> Q antiquark
   *  Q1 diquark, Q2bar quark              --> Q antiquark
   *  Q1 anti-diquark, Q2bar antiquark;    --> Q quark
   **************************/

  // Calculate the unit three-vector, in the C frame, along which
  // all of the constituents and children clusters move.
  Lorentz5Momentum u(p0Q1);
  u.boost( -pClu.boostVector() );        // boost from LAB to C
  // the unit three-vector is then  u.vect().unit()
						       
  // Calculate the momenta of C1 and C2 in the (parent) C frame first, 
  // where the direction of C1 is u.vect().unit(), and then boost back in the 
  // LAB frame.

  if (pClu.m() < pClu1.mass() + pClu2.mass() ) {
    throw Exception() << "Impossible Kinematics in ClusterFissioner::calculateKinematics() (A)" 
		      << Exception::eventerror;
  }
  Kinematics::twoBodyDecay(pClu, pClu1.mass(), pClu2.mass(), 
			   u.vect().unit(), pClu1, pClu2);

  // In the case that cluster1 does not decay immediately into a single hadron,
  // calculate the momenta of Q1 (as constituent of C1) and Qbar in the
  // (parent) C1 frame first, where the direction of Q1 is u.vect().unit(), 
  // and then boost back in the LAB frame. 
  if(!toHadron1) {
    if (pClu1.m() < pQ1.mass() + pQbar.mass() ) {
      throw Exception() << "Impossible Kinematics in ClusterFissioner::calculateKinematics() (B)" 
			<< Exception::eventerror;
    }
    Kinematics::twoBodyDecay(pClu1, pQ1.mass(), pQbar.mass(), 
			     u.vect().unit(), pQ1, pQbar);
  }

  // In the case that cluster2 does not decay immediately into a single hadron,
  // Calculate the momenta of Q and Q2bar (as constituent of C2) in the
  // (parent) C2 frame first, where the direction of Q is u.vect().unit(), 
  // and then boost back in the LAB frame. 
  if(!toHadron2) {
    if (pClu2.m() < pQ.mass() + pQ2bar.mass() ) {
      throw Exception() << "Impossible Kinematics in ClusterFissioner::calculateKinematics() (C)" 
			<< Exception::eventerror;
    }
    Kinematics::twoBodyDecay(pClu2, pQ.mass(), pQ2bar.mass(), 
			     u.vect().unit(), pQ, pQ2bar);
  }
}


void ClusterFissioner::calculatePositions(const Lorentz5Momentum & pClu,
					  const LorentzPoint & positionClu,
					  const Lorentz5Momentum & pClu1,
					  const Lorentz5Momentum & pClu2,
					  LorentzPoint & positionClu1,
					  LorentzPoint & positionClu2) const 
{
  // Determine positions of cluster children.
  // See Marc Smith's thesis, page 127, formulas (4.122) and (4.123).
  Energy Mclu  = pClu.m();
  Energy Mclu1 = pClu1.m();
  Energy Mclu2 = pClu2.m();

  // Calculate the unit three-vector, in the C frame, along which 
  // children clusters move.
  Lorentz5Momentum u(pClu1);
  u.boost( -pClu.boostVector() );        // boost from LAB to C frame
  // the unit three-vector is then  u.vect().unit()

  Energy pstarChild = Kinematics::pstarTwoBodyDecay(Mclu,Mclu1,Mclu2);

  // First, determine the relative positions of the children clusters
  // in the parent cluster reference frame.
  double GeV2mm = _globalParameters->conversionFactorGeVtoMillimeter();
  Length x1 = GeV2mm * (Mclu*0.25 + 0.5
                     *(pstarChild + (sqr(Mclu2) - sqr(Mclu1))/(2.0*Mclu)))/GeV;
  Length t1 = ((Mclu/GeV) * GeV2mm - x1); 
  LorentzDistance distanceClu1( x1 * u.vect().unit(), t1 );
  Length x2 = GeV2mm * (-Mclu*0.25 + 0.5
                     *(-pstarChild + (sqr(Mclu2) - sqr(Mclu))/(2.0*Mclu)))/GeV;
  Length t2 = ((Mclu/GeV) * GeV2mm + x2);
  LorentzDistance distanceClu2(x2 * u.vect().unit(), t2);

  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Hadronization ) {
    generator()->log() << "ClusterFissioner::calculatePositions : *** extreme debugging ***" << endl
                       << "\t distanceClu1 = " << distanceClu1 
                       << "\t invariant length = " << distanceClu1.mag() << "  [mm] " << endl
                       << "\t distanceClu2 = " << distanceClu2 
                       << "\t invariant length = " << distanceClu2.mag() << "  [mm] " << endl;
  }

  // Then, transform such relative positions from the parent cluster
  // reference frame to the Lab frame.
  distanceClu1.boost( pClu.boostVector() );
  distanceClu2.boost( pClu.boostVector() );

  // Finally, determine the absolute positions in the Lab frame.
  positionClu1 = positionClu + distanceClu1;
  positionClu2 = positionClu + distanceClu2;

}

