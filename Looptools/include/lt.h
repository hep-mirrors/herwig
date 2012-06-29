* lt.h
* internal common blocks for the LoopTools routines
* this file is part of LoopTools
* last modified 21 Dec 10 th


#include "ff.h"

* the cache-pointer structure is (see cache.c):
* 1. int valid
* 2. Node *last
* 3. Node *first
* 4. (not used)

	integer ncaches
	parameter (ncaches = 8)

	integer*8 cacheptr(4,KIND,ncaches)
	integer*8 savedptr(2,ncaches)
	double precision maxdev
	integer warndigits, errdigits
	integer serial, versionkey
	integer debugkey, debugfrom, debugto

	common /ltvars/
     &    cacheptr, savedptr,
     &    maxdev,
     &    warndigits, errdigits,
     &    serial, versionkey,
     &    debugkey, debugfrom, debugto

	integer cmpbits

	common /ltcache/ cmpbits

	double complex cache(2,ncaches)
	equivalence (cacheptr, cache)

#ifndef sig
#define sig(c) int(sign(1D0,DBLE(r))
#define DEBUGLEVEL ibits(debugkey,8,2)
#endif

