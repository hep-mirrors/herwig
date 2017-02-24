#ifndef HERWIG_RPV_HELPER_H
#define HERWIG_RPV_HELPER_H

namespace Herwig {

namespace RPV_helper {

  unsigned int neutralinoIndex(long id) {
    if(id> 1000000) {
      return id<1000025 ? id-1000022 : (id-1000005)/10;
    }
    else if(abs(id)<=16) {
      return (abs(id)-4)/2;
    }
    else {
      return id-13;
    }
  }
  
  unsigned int charginoIndex(long id) {
    return abs(id)>1000000 ? (abs(id)-1000024)/13 : (abs(id)-7)/2;
  }

}

}

#endif
