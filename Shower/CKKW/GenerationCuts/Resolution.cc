#include "Resolution.h"

#ifdef useResolution_KtTilde

#include "Herwig++/Shower/CKKW/Clustering/Partitioner.h"
#include <utility>
#include <cmath>

Herwig::Partitioner partitioner;

bool resolvable (const FiveMomentum& firstEmission,
		 const FiveMomentum& secondEmission,
		 const FiveMomentum& spectator) {

  // currently cannot handle initial state, so we consider
  // these resolvable
  if (firstEmission.initial ||
      secondEmission.initial ||
      spectator.initial)
    return true;
  
  FiveMomentum p = add(firstEmission,secondEmission);
  double q2 = dotProduct(p,p);

  FiveMomentum n = spectator;
  n.E = std::sqrt(sqr(n.px)+sqr(n.py)+sqr(n.pz));
  n.invMass2 = 0.;

  double z = dotProduct(n,firstEmission)/dotProduct(n,p);
  double qt2 = q2/(z*(1.-z));

  return (min(sqr(z),sqr(1-z))*qt2 > resolution);

}

void checkCut_ (int * passFlag) {

  *passFlag = 1;

  std::vector<std::vector<FiveMomentum> > parts
    = partitioner.partitions(jets,3);

  for (std::vector<std::vector<FiveMomentum> >::iterator c = parts.begin();
       c != parts.end(); ++c)
    if (!resolvable((*c)[0],(*c)[1],(*c)[2]) ||
	!resolvable((*c)[0],(*c)[2],(*c)[1]) ||
	!resolvable((*c)[1],(*c)[2],(*c)[0])) {
      *passFlag = 0;
      return;
    }

}

#else

void checkCut_ (int * passFlag) {
  *passFlag = 1;
}

#endif

