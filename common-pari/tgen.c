#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "c-target.h"

/*
  tgen: congruence conditions modulo p for a tame ramification pattern, for
  any prime p dividing none of the ramification indices (degree 8, which
  has no hand-written functions; the method works for any degree).

  If the pattern has the index e exactly k_e times, the characteristic
  polynomials mod p of the integral elements of the algebras with this
  pattern are exactly the products  prod_e G_e(x)^e  with G_e monic of
  degree k_e (any factorization of G_e is realized by some splitting into
  residue degrees).  Write f = R*H with R = prod_{e>=2} G_e^e of degree s
  and H = G_1 arbitrary of degree DEG - s.

  In the search's loop order (A_1, A_2, A_DEG, A_3, ..., A_{DEG-1}) the
  first DEG - s coefficients are free.  At the next coefficient the
  function runs through every R; for each, f = 0 mod R is a linear system
  for the remaining s coefficients mod p, solved by elimination (all
  solutions are enumerated when the system is degenerate).  This is the
  same scheme as the hand-written functions of the lower degrees: for
  pattern 2,1,...,1 it enumerates the double root a.

  The work per value of the free coefficients is p^(number of ramified
  parts).  When that exceeds TGEN_MAXENUM (default 50; environment
  variable of the same name), enumerating the residue classes costs more
  than it saves, since the coefficient ranges are then smaller than p.
  In that case no congruence is imposed on the loops, and instead the
  modular filter (c-qfilter.c) requires p | disc(f), which holds for every
  ramified pattern.
*/

#define pp plist[pnum]
#define MAXD 10

static int Lord[MAXD];          /* loop order of the coefficients */

typedef struct {
  int    s, k1, nram;           /* deg R, number of free coefs, ramified parts */
  int    ek[MAXD+1];            /* ek[e] = k_e for e >= 2 */
  long   Rcount, Ridx;          /* the R's: Ridx in [0, p^nram) */
  long   freev[MAXD+1];         /* residues of the free coefficients */
  int    det[MAXD];             /* the determined coefficients */
  long   part[MAXD];            /* particular solution */
  long   nullb[MAXD][MAXD];     /* basis of the solution space */
  int    ndim;
  long   combo[MAXD];           /* current combination of the basis */
} tgstate;

static tgstate tg[10];
int tgen_filter[10];            /* 1: filter mode (p | disc(f)) for this prime */

static long md(long a, long q) { a %= q; return a < 0 ? a + q : a; }
static long pinv(long a, long p) {
  long r = 1, b = md(a, p), e = p - 2;
  while(e) { if(e & 1) r = r*b % p; b = b*b % p; e >>= 1; }
  return r;
}

static void setup_order(void) {
  int i, k = 0;
  Lord[k++] = 1;
  if(DEG >= 2) Lord[k++] = 2;
  if(DEG >= 3) Lord[k++] = DEG;
  for(i=3; i<DEG; i++) Lord[k++] = i;
}

/* R for index idx, coefficients low -> high, degree s, monic.  */
static void build_R(long pnum, long idx, long *R) {
  tgstate *T = &tg[pnum];
  long p = pp, cur[MAXD+1], tmp[2*MAXD+2];
  int e, k, i, j, t, deg = 0;
  cur[0] = 1;
  for(e=2; e<=DEG; e++) {
    long G[MAXD+1];
    int kk = T->ek[e];
    if(kk == 0) continue;
    for(i=0; i<kk; i++) { G[i] = idx % p; idx /= p; }
    G[kk] = 1;
    for(t=0; t<e; t++) {               /* cur *= G */
      for(i=0; i<=deg+kk; i++) tmp[i] = 0;
      for(i=0; i<=deg; i++) for(j=0; j<=kk; j++)
        tmp[i+j] = (tmp[i+j] + cur[i]*G[j]) % p;
      deg += kk;
      for(i=0; i<=deg; i++) cur[i] = tmp[i];
    }
  }
  for(k=0; k<=deg; k++) R[k] = cur[k];
}

/* Solve f = 0 mod R for the determined coefficients.  Returns 0 if there is
   no solution. */
static int solve(long pnum) {
  tgstate *T = &tg[pnum];
  long p = pp, R[MAXD+1], rm[MAXD+1][MAXD], M[MAXD][MAXD+1];
  int s = T->s, i, j, m, row, col, piv[MAXD], npiv = 0, isdet[MAXD+1];

  build_R(pnum, T->Ridx, R);
  /* rm[m] = x^m mod R, as a vector of length s */
  for(m=0; m<=DEG; m++) {
    for(i=0; i<s; i++) rm[m][i] = 0;
    if(m < s) { rm[m][m] = 1; continue; }
    /* x * rm[m-1] */
    {
      long top = rm[m-1][s-1];
      for(i=s-1; i>=1; i--) rm[m][i] = rm[m-1][i-1];
      rm[m][0] = 0;
      for(i=0; i<s; i++) rm[m][i] = md(rm[m][i] - top*R[i], p);
    }
  }
  for(i=0; i<=DEG; i++) isdet[i] = 0;
  for(j=0; j<s; j++) isdet[T->det[j]] = 1;
  /* f = x^DEG + sum A_i x^(DEG-i); known part */
  for(i=0; i<s; i++) {
    long v = rm[DEG][i];
    int a;
    for(a=1; a<=DEG; a++) if(!isdet[a]) v += T->freev[a] * rm[DEG-a][i];
    M[i][s] = md(-v, p);
    for(j=0; j<s; j++) M[i][j] = rm[DEG - T->det[j]][i];
  }
  /* row reduce [M | rhs] */
  row = 0;
  for(col=0; col<s && row<s; col++) {
    int r;
    for(r=row; r<s && M[r][col]==0; r++) ;
    if(r == s) continue;
    for(j=0; j<=s; j++) { long t = M[row][j]; M[row][j] = M[r][j]; M[r][j] = t; }
    {
      long iv = pinv(M[row][col], p);
      for(j=0; j<=s; j++) M[row][j] = M[row][j]*iv % p;
    }
    for(r=0; r<s; r++) if(r != row && M[r][col]) {
      long f = M[r][col];
      for(j=0; j<=s; j++) M[r][j] = md(M[r][j] - f*M[row][j], p);
    }
    piv[npiv++] = col;
    row++;
  }
  for(i=row; i<s; i++) if(M[i][s]) return 0;          /* inconsistent */
  for(j=0; j<s; j++) T->part[j] = 0;
  for(i=0; i<npiv; i++) T->part[piv[i]] = M[i][s];
  /* basis of the null space: one vector per free column */
  T->ndim = 0;
  for(col=0; col<s; col++) {
    int isp = 0;
    for(i=0; i<npiv; i++) if(piv[i] == col) isp = 1;
    if(isp) continue;
    for(j=0; j<s; j++) T->nullb[T->ndim][j] = 0;
    T->nullb[T->ndim][col] = 1;
    for(i=0; i<npiv; i++) T->nullb[T->ndim][piv[i]] = md(-M[i][col], p);
    T->ndim++;
  }
  for(i=0; i<T->ndim; i++) T->combo[i] = 0;
  return 1;
}

static void set_residues(long pnum) {
  tgstate *T = &tg[pnum];
  long p = pp;
  int j, i;
  for(j=0; j<T->s; j++) {
    long v = T->part[j];
    for(i=0; i<T->ndim; i++) v += T->combo[i]*T->nullb[i][j];
    pdata[pnum][T->det[j]][0] = p;
    pdata[pnum][T->det[j]][1] = md(v, p);
  }
}

static int first_R(long pnum) {
  tgstate *T = &tg[pnum];
  for(T->Ridx = 0; T->Ridx < T->Rcount; T->Ridx++)
    if(solve(pnum)) { set_residues(pnum); return 1; }
  return 0;
}

static long advance(long pnum) {
  tgstate *T = &tg[pnum];
  int i;
  for(i=0; i<T->ndim; i++) {
    if(++T->combo[i] < pp) { set_residues(pnum); return 1; }
    T->combo[i] = 0;
  }
  for(T->Ridx++; T->Ridx < T->Rcount; T->Ridx++)
    if(solve(pnum)) { set_residues(pnum); return 1; }
  return 0;
}

int tgen_init(long pnum, long which, long what) {
  tgstate *T = &tg[pnum];
  int t, j;

  if(which == 0) {
    const char *pat = pattern_string(cases[pnum][what].start);
    const char *q;
    setup_order();
    if(pat == NULL) { fprintf(stderr, "tgen: unknown pattern\n"); exit(1); }
    memset(T, 0, sizeof(*T));
    for(q = pat; *q; q++) {
      int e = *q - '0';
      if(e % pp == 0) { fprintf(stderr, "tgen: pattern is wild at %ld\n", pp); exit(1); }
      if(e >= 2) { T->ek[e]++; T->s += e; T->nram++; }
    }
    T->k1 = DEG - T->s;
    T->Rcount = 1;
    for(j=0; j<T->nram; j++) T->Rcount *= pp;
    for(j=0; j<T->s; j++) T->det[j] = Lord[T->k1 + j];
    for(j=1; j<=DEG; j++) {
      pflags[pnum][j] = CALL_NONE;
      pdata[pnum][j][0] = 1;
      pdata[pnum][j][1] = 0;
    }
    {
      const char *ev = getenv("TGEN_MAXENUM");
      long lim = ev ? atol(ev) : 50;
      tgen_filter[pnum] = (T->Rcount > lim);
    }
    pdata[pnum][0][0] = DEG - (int) strlen(pat);   /* discriminant exponent */
    if(tgen_filter[pnum]) return 1;   /* no congruences on the loops */
    for(j=0; j<T->s; j++) pdata[pnum][T->det[j]][0] = pp;
    for(t=1; t<T->k1; t++) pflags[pnum][Lord[t]] = CALL_INIT;
    pflags[pnum][Lord[T->k1]] = CALL_BOTH;
    return 1;
  }
  if(tgen_filter[pnum]) return 1;
  /* capture the value of the previous coefficient in the loop order */
  for(t=1; t<=T->k1; t++)
    if(Lord[t] == which) T->freev[Lord[t-1]] = md(what, pp);
  if(which == Lord[T->k1]) return first_R(pnum);
  return 1;
}

long tgen_next(long pnum, long which) {
  if(tgen_filter[pnum] || which != Lord[tg[pnum].k1]) return 0;
  return advance(pnum);
}
