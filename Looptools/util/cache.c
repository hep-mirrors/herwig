/*
	cache.c
		caching of tensor coefficients in
		dynamically allocated memory
		this file is part of LoopTools
		last modified 9 Dec 10 th
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define cachelookup_ ljcachelookup_

#define cachelookup cachelookup_
#define ltcache ltcache_

#ifndef BIGENDIAN
#define BIGENDIAN 0
#endif

#ifndef KIND
#define KIND 1
#endif

#if KIND == 2
#define MSB (1-BIGENDIAN)
#else
#define MSB 0
#endif


typedef long long dblint;

typedef unsigned long long udblint;

typedef struct { dblint part[KIND]; } Real;

typedef struct { Real re, im; } Complex;

typedef long Integer;


extern struct {
  long cmpbits;
} ltcache;


/* (a < 0) ? -1 : 0 */
#define NegQ(a) ((a) >> (sizeof(a)*8 - 1))

/* (a < 0) ? 0 : a */
#define IDim(a) ((a) & NegQ(-(a)))


static int SignBit(const dblint i)
{
  return (udblint)i >> (8*sizeof(i) - 1);
}


static Integer PtrDiff(const void *a, const void *b)
{
  return (char *)a - (char *)b;
}


static dblint CmpPara(const Real *para1, const Real *para2, int n,
  const dblint mask)
{
  while( n-- ) {
    const dblint c = (mask & para1->part[MSB]) -
                     (mask & para2->part[MSB]);
    if( c ) return c;
    ++para1;
    ++para2;
  }
  return 0;
}


#if KIND == 2

static dblint CmpParaLo(const Real *para1, const Real *para2, int n,
  const dblint mask)
{
  while( n-- ) {
    dblint c = para1->part[MSB] - para2->part[MSB];
    if( c ) return c;
    c = (mask & para1->part[1-MSB]) - (mask & para2->part[1-MSB]);
    if( c ) return c;
    ++para1;
    ++para2;
  }
  return 0;
}

#endif


Integer cachelookup(const Real *para, double *base,
  void (*calc)(const Real *, Real *, const long *),
  const int *pnpara, const int *pnval)
{

  const long one = 1;
  const Integer C_size = sizeof(Complex);
  const int npara = *pnpara, nval = *pnval;

  typedef struct node {
    struct node *next[2], *succ;
    int serial;
    Real para[2];
  } Node;

#define base_valid (int *)&base[0]
#define base_last (Node ***)&base[1]
#define base_first (Node **)&base[2]

  const int valid = *base_valid;
  Node **last = *base_last;
  Node **next = base_first;
  Node *node;

  if( last == NULL ) last = next;

  if( ltcache.cmpbits > 0 ) {
    dblint mask = -(1ULL << IDim(64 - ltcache.cmpbits));
#if KIND == 2
    dblint (*cmp)(const Real *, const Real *, int, const dblint) = CmpPara;
    if( ltcache.cmpbits >= 64 ) {
      mask = -(1ULL << IDim(128 - ltcache.cmpbits));
      cmp = CmpParaLo;
    }
#else
#define cmp CmpPara
#endif

    while( (node = *next) && node->serial < valid ) {
      const dblint i = cmp(para, node->para, npara, mask);
      if( i == 0 ) {
        goto found;
      }
      next = &node->next[SignBit(i)];
    }
  }

  node = *last;

  if( node == NULL ) {
	/* The "Real para[2]" bit in Node is effectively an extra
	   Complex for alignment so that node can be reached with
	   an integer index into base */
    node = malloc(sizeof(Node) + npara*sizeof(Real) + nval*sizeof(Complex));
    if( node == NULL ) {
      fputs("Out of memory for LoopTools cache.\n", stderr);
      exit(1);
    }
    node = (Node *)((char *)node +
      (PtrDiff(base, &node->para[npara]) & (sizeof(Complex) - 1)));
    node->succ = NULL;
    node->serial = valid;
    *last = node;
  }

  *next = node;
  *base_last = &node->succ;
  *base_valid = valid + 1;

  node->next[0] = NULL;
  node->next[1] = NULL;

  memcpy(node->para, para, npara*sizeof(Real));
  calc(node->para, &node->para[npara], &one);

found:
  return PtrDiff(&node->para[npara], base)/C_size;
}

