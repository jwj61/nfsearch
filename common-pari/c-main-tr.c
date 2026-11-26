/*
  Totally real fields of degree DEG with discriminant at most B, by
  Voight's method (J. Voight, "Enumeration of totally real number fields of bounded
  root discriminant", ANTS VIII, 2008).

  Shared by c4-tr (4/c4-main-tr.c) and c5-tr (5-pari/c5-main-tr.c), which
  are links to this file.

  usage:  cN-tr [-jN] [-s] B

    -jN  (or -j N) use N processes
    -s   keep only fields whose discriminant is a (positive) perfect square

  Output goes to TRSearch_B (TRSquare_B with -s), one polredabs polynomial
  per line; a field can appear more than once.

  The polynomials are f = x^n + A1 x^(n-1) + ... + An, with the
  coefficients chosen from the top down.

    A1   0 <= A1 <= n/2 (translation x -> x+k and x -> -x), as in cN-stand.

    A2   For f totally real, T2 = sum of the squares of the roots
         = A1^2 - 2 A2 exactly, so Hunter's bound T2 <= A1^2/n + gamma
         (the same bound, with the same slack, as cN-stand) and n <= T2 (AM-GM, |N| >= 1),
         T2 > A1^2/n (the roots are distinct) give the range of A2.

    Aj   (j >= 2)  g = f^(n-j), the (n-j)-th derivative, has degree j and
         is known except for its constant term C = (n-j)! Aj.  Its
         derivative g' = f^(n-j+1) is known, with real roots
         b_1 < ... < b_(j-1) (found at the previous level).  By Rolle the
         roots of g interlace them, so (-1)^(j-i) g(b_i) > 0 for each i:
         a lower or upper bound on C at every b_i.  Also, by the
         Laguerre-Samuelson inequality all roots of g lie in [L,U], where
         L,U = mu -+ sqrt(j-1) sigma and mu, sigma^2 are the mean and the
         variance of the roots of g (known from its top three
         coefficients when j >= 3), so g(U) >= 0 and (-1)^j g(L) >= 0.

    An   in addition 1 <= |An| <= (T2/n)^(n/2) (AM-GM).

  The roots of g, needed at the next level, are found by bisection in the
  intervals [L,b_1], [b_1,b_2], ..., [b_(j-1),U], each of which contains
  exactly one of them.  Since g'(b_i) = 0, an error in b_i changes g(b_i)
  only to second order; the bounds on C are widened by a small margin to
  allow for rounding, and every polynomial that is written out is checked
  to be totally real with Sturm's theorem.

  With A1 = 0, x -> -x negates A3, so A3 >= 0 
  (c5-stand takes A5 >= 0 instead; either choice is valid, not both).

  In degree 4 Hunter's bound gives every primitive field but can miss
  imprimitive ones (C4, V4, D4), as in c4-stand.  In degree 5 every field
  is primitive, so the search is complete.

  With -s, before the exact discriminant D of f is computed, D is tested
  for being a square modulo each of the moduli in sq_mods (see sq_setup).
  The discriminant of the generic polynomial x^n + a1 x^(n-1) + ... + an
  is computed once, as a polynomial in a1, ..., an; for each value of
  A_(n-1) its reduction mod m as a polynomial in An is formed (for each
  modulus only when first needed), and each An costs a Horner evaluation
  and a table lookup, usually for only one or two moduli.

  With -jN the work is shared out by the values of A_(n-1): the workers
  walk the same loops, and claim each A_(n-1) through par_claim.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include "c-target.h"
#include "c-doout.h"

typedef long double ldbl;

/* added to Hunter's bound on T2, as in c4-stand and c5-stand */
#if DEG==4
#define HUNTER_SLACK (9./7.)
#elif DEG==5
#define HUNTER_SLACK (5./4.)
#else
#define HUNTER_SLACK (5./4.)
#endif

static ldbl  cf[DEG+1];          /* cf[i] = coefficient of x^i in f */
static ldbl  rts[DEG+1][DEG+1];  /* rts[j]: the j roots of f^(n-j), increasing */
static ldbl  fact[DEG+1];
static long long Acur[DEG+1];    /* Acur[j] = A_j */
static double t2max;             /* Hunter bound on T2 for the current A1 */
static GEN    bigbgen;
static int    squareonly = 0;
static FILE  *outfile;
static long   ltop;
static long   nrejsturm = 0;     /* passed the bounds, not totally real */

#define MARGIN(v) (1e-9L * (1.0L + fabsl(v)))

/* g = f^(k), the k-th derivative; coefficient of x^(i-k) is
   cf[i] i!/(i-k)!.  Value at x without the constant term (with=0) or
   with it (with=1). */
static ldbl deriv_eval(int k, ldbl x, int with)
{
  int i;
  ldbl v = 0;
  for(i=DEG; i>k; i--) v = v*x + cf[i]*fact[i]/fact[i-k];
  v *= x;
  if(with) v += cf[k]*fact[k];
  return v;
}

/* Laguerre-Samuelson interval [*L,*U] for the roots of f^(k) (degree
   m = n-k >= 2); needs cf[k..n].  Returns 0 if the variance is not
   positive (then f^(k), hence f, cannot have distinct real roots). */
static int lag_interval(int k, ldbl *L, ldbl *U)
{
  int m = DEG - k;
  ldbl g2 = cf[DEG]*fact[DEG]/fact[DEG-k];        /* x^m */
  ldbl g1 = cf[DEG-1]*fact[DEG-1]/fact[DEG-1-k];  /* x^(m-1) */
  ldbl g0 = cf[DEG-2]*fact[DEG-2]/fact[DEG-2-k];  /* x^(m-2) */
  ldbl s1 = -g1/g2, s2 = s1*s1 - 2*g0/g2;
  ldbl mu = s1/m, var = s2/m - mu*mu, w;
  if(var <= 0) return 0;
  w = sqrtl((m-1)*var);
  *L = mu - w - MARGIN(mu) - MARGIN(w);
  *U = mu + w + MARGIN(mu) + MARGIN(w);
  return 1;
}

/* The one root of f^(k) in [lo,hi], by bisection. */
static ldbl bisect_root(int k, ldbl lo, ldbl hi)
{
  int it;
  ldbl slo = deriv_eval(k, lo, 1);
  for(it=0; it<200 && hi-lo > 0; it++) {
    ldbl mid = (lo+hi)/2, sm;
    if(mid <= lo || mid >= hi) break;
    sm = deriv_eval(k, mid, 1);
    if((sm > 0) == (slo > 0)) { lo = mid; slo = sm; }
    else hi = mid;
  }
  return (lo+hi)/2;
}

/* After A_j is set: the j roots of g = f^(n-j) into rts[j].  Returns 0
   if g cannot have j distinct real roots. */
static int find_roots(int j)
{
  int k = DEG - j, i;
  ldbl L, U;
  if(j == 1) { rts[1][0] = -cf[DEG-1]*fact[DEG-1]/fact[DEG]; return 1; }
  if(!lag_interval(k, &L, &U)) return 0;
  for(i=0; i<j; i++) {
    ldbl lo = (i == 0)   ? L : rts[j-1][i-1];
    ldbl hi = (i == j-1) ? U : rts[j-1][i];
    rts[j][i] = bisect_root(k, lo, hi);
  }
  return 1;
}

/* Bounds on C, the constant term of g = f^(n-j), from interlacing with
   the roots of g' and from the Laguerre-Samuelson interval.  Returns 0 if
   the bounds are empty. */
static int rolle_bounds(int j, ldbl *clo, ldbl *chi)
{
  int k = DEG - j, i;
  ldbl lo = -1e30L, hi = 1e30L, h, L, U;
  for(i=0; i<j-1; i++) {        /* (-1)^(j-1-i) g(b_i) > 0 (b_i = rts[j-1][i], 0-based) */
    h = deriv_eval(k, rts[j-1][i], 0);
    if(((j-1-i) & 1) == 0) { if(-h > lo) lo = -h; }    /* g(b) > 0 */
    else                   { if(-h < hi) hi = -h; }    /* g(b) < 0 */
  }
  if(j >= 3) {
    if(!lag_interval(k, &L, &U)) return 0;
    h = deriv_eval(k, U, 0);    /* g(U) >= 0 */
    if(-h > lo) lo = -h;
    h = deriv_eval(k, L, 0);    /* (-1)^j g(L) >= 0 */
    if(j & 1) { if(-h < hi) hi = -h; }
    else      { if(-h > lo) lo = -h; }
  }
  *clo = lo - MARGIN(lo);
  *chi = hi + MARGIN(hi);
  return *clo <= *chi;
}

static void set_coef(int j, long long v)
{
  cf[DEG-j] = (ldbl) v;
  Acur[j] = v;
  if(j < DEG) myaffsi(v, coefs[j]);   /* An: in test_poly, after the cheap tests */
  switch(j) {
  case 1: A1 = v; break;
  case 2: A2 = v; break;
  case 3: A3 = v; break;
  case 4: A4 = v; break;
  case 5: A5 = v; break;
  case 6: A6 = v; break;
  case 7: A7 = v; break;
  case 8: A8 = v; break;
  }
}

/* ---- -s: modular test that the discriminant can be a square ---- */

/* the order matters a little: most polynomials are rejected by the first
   one or two moduli */
static const long sq_mods[] = {64, 63, 65, 11, 17, 19, 23, 29, 31, 37, 41, 43, 47};
#define NSQMODS ((int)(sizeof(sq_mods)/sizeof(sq_mods[0])))
#define SQMODMAX 65

static unsigned char sq_issq[NSQMODS][SQMODMAX];   /* squares mod m */
static long  sq_nmono = 0, sq_maxe = 0, sq_maxpow = 0;
static long (*sq_exp)[DEG+1] = NULL;  /* exponents of a1..an in each monomial of disc */
static long *sq_coef[NSQMODS];        /* its coefficient mod m */
static long  sq_dc[NSQMODS][2*DEG];   /* disc mod m as a polynomial in An, this A_(n-1) */
static int   sq_valid[NSQMODS];       /* sq_dc[mi] computed for this A_(n-1) */
static long  sq_rej = 0;              /* polynomials rejected by the test */
static long  sq_vars[DEG+1];          /* PARI variable numbers of a1..an */

static void sq_walk(GEN p, long *e)
{
  long k, i;
  if(typ(p) == t_INT) {
    if(!signe(p)) return;
    for(i=1; i<=DEG; i++) {
      sq_exp[sq_nmono][i] = e[i];
      if(i < DEG && e[i] > sq_maxpow) sq_maxpow = e[i];
    }
    if(e[DEG] > sq_maxe) sq_maxe = e[DEG];
    for(k=0; k<NSQMODS; k++) sq_coef[k][sq_nmono] = smodis(p, sq_mods[k]);
    sq_nmono++;
    return;
  }
  if(typ(p) != t_POL) pari_err_TYPE("sq_walk", p);
  for(i=1; i<=DEG && sq_vars[i] != varn(p); i++) ;
  if(i > DEG) pari_err_TYPE("sq_walk", p);
  for(k=0; k<=degpol(p); k++) { e[i] = k; sq_walk(gel(p, k+2), e); }
  e[i] = 0;
}

static long sq_count(GEN p)
{
  long k, n = 0;
  if(typ(p) == t_INT) return signe(p) ? 1 : 0;
  for(k=0; k<=degpol(p); k++) n += sq_count(gel(p, k+2));
  return n;
}

/* the discriminant of x^n + a1 x^(n-1) + ... + an, as a list of monomials */
static void sq_setup(void)
{
  char str[400], *t;
  long e[DEG+1], k, i;
  pari_sp av = avma;
  GEN D;
  for(k=0; k<NSQMODS; k++) {
    memset(sq_issq[k], 0, SQMODMAX);
    for(i=0; i<sq_mods[k]; i++) sq_issq[k][(i*i) % sq_mods[k]] = 1;
  }
  t = str + sprintf(str, "poldisc(x^%d", DEG);
  for(i=1; i<=DEG; i++) t += sprintf(t, "+trsq_a%ld*x^%ld", i, DEG-i);
  sprintf(t, ",x)");
  for(i=1; i<=DEG; i++) {
    char nm[20];
    sprintf(nm, "trsq_a%ld", i);
    sq_vars[i] = varn(gp_read_str(nm));
  }
  D = gp_read_str(str);
  k = sq_count(D);
  sq_exp = malloc(k * sizeof(*sq_exp));
  for(i=0; i<NSQMODS; i++) sq_coef[i] = malloc(k * sizeof(long));
  if(sq_exp == NULL) { perror("malloc"); exit(1); }
  for(i=0; i<=DEG; i++) e[i] = 0;
  sq_walk(D, e);
  if(sq_maxpow >= 2*DEG*DEG || sq_maxe >= 2*DEG) { printf("sq_setup: degree too large\n"); exit(1); }
  avma = av;
}

static void sq_newnode(void)
{
  int k;
  for(k=0; k<NSQMODS; k++) sq_valid[k] = 0;
}

/* disc mod m as a polynomial in An, from the current A1..A_(n-1) */
static void sq_compute(int k)
{
  long m = sq_mods[k], i, t, x;
  long pw[DEG][2*DEG*DEG];     /* pw[i][t] = A_i^t mod m */
  for(i=1; i<DEG; i++) {
    x = (long)(Acur[i] % m); if(x < 0) x += m;
    pw[i][0] = 1 % m;
    for(t=1; t<=sq_maxpow; t++) pw[i][t] = pw[i][t-1]*x % m;
  }
  for(t=0; t<=sq_maxe; t++) sq_dc[k][t] = 0;
  for(t=0; t<sq_nmono; t++) {
    long v = sq_coef[k][t];
    for(i=1; i<DEG; i++) v = v*pw[i][sq_exp[t][i]] % m;
    sq_dc[k][sq_exp[t][DEG]] += v;
  }
  for(t=0; t<=sq_maxe; t++) {
    sq_dc[k][t] %= m;
    if(sq_dc[k][t] < 0) sq_dc[k][t] += m;
  }
  sq_valid[k] = 1;
}

/* 0 if the discriminant of f (with An = Acur[DEG]) is not a square mod
   one of the moduli */
static int sq_maybe_square(void)
{
  int k;
  long m, x, v, t;
  for(k=0; k<NSQMODS; k++) {
    m = sq_mods[k];
    if(!sq_valid[k]) sq_compute(k);
    x = (long)(Acur[DEG] % m); if(x < 0) x += m;
    v = sq_dc[k][sq_maxe];
    for(t=sq_maxe-1; t>=0; t--) v = (v*x + sq_dc[k][t]) % m;
    if(!sq_issq[k][v]) return 0;
  }
  return 1;
}

/* f is set: test it and write out the field if it qualifies. */
static void test_poly(void)
{
  polcount++;
  if(squareonly && !sq_maybe_square()) { sq_rej++; return; }
  myaffsi(Acur[DEG], coefs[DEG]);
  td = discsr(a);
  if(signe(td) <= 0) { avma = ltop; return; }     /* totally real: D > 0 */
  if(squareonly && !Z_issquare(td)) { avma = ltop; return; }
  /* necessary condition for d_K <= B: D = d_K i^2 has a square divisor >= D/B */
  if(!disc_may_be_small(td, bigbgen)) { avma = ltop; return; }
  postdisc++;
  if(polisirreducible(a)) {
    if(ZX_sturm(a) != DEG) { nrejsturm++; avma = ltop; return; }
    td = discf(a);
    if(gcmp(td, bigbgen) <= 0) {
      discfcount++;
      pari_fprintf(outfile, "%Ps\n", polredabs0(a, 0));
    }
  }
  avma = ltop;
}

/* Choose A_j, then go on to A_(j+1). */
static void level(int j)
{
  long long lo, hi, v;
  ldbl clo, chi, fk = fact[DEG-j];

  if(j == 2) {
    double t2lo = (A1*A1)/(double) DEG;
    if(t2lo < DEG) t2lo = DEG;
    /* T2 = A1^2 - 2 A2 in [t2lo, t2max] */
    lo = (long long) ceil((A1*A1 - t2max)/2.0 - 1e-9);
    hi = (long long) floor((A1*A1 - t2lo)/2.0 + 1e-9);
  } else {
    if(!rolle_bounds(j, &clo, &chi)) return;
    lo = (long long) ceill(clo/fk);
    hi = (long long) floorl(chi/fk);
  }

  if(j == 3 && A1 == 0 && lo < 0) lo = 0;   /* x -> -x */

  if(j == DEG) {
    /* 1 <= |An| <= (T2/n)^(n/2) */
    double t2 = (double)(A1*A1 - 2*A2);
    long long nb = (long long) floor(pow(t2/DEG, DEG/2.0) + 1e-6);
    if(lo < -nb) lo = -nb;
    if(hi >  nb) hi =  nb;
  }

  for(v=lo; v<=hi; v++) {
    if(j == DEG && v == 0) continue;
    set_coef(j, v);
    if(j == DEG) { test_poly(); continue; }
    if(j == DEG-1) {
      if(!par_claim()) continue;               /* another worker has it */
      if(squareonly) sq_newnode();
    }
    if(!find_roots(j)) continue;
    if(j == 2 && verbosity > 2) {
      printf("A1 = %lld, A2 = %lld\n", A1, A2);
      fflush(stdout);
    }
    level(j+1);
  }
}

int main(int argc, char *argv[])
{
  struct timeval time1, time2;
  struct timezone tzn;
  FILE  *afile;
  long   timetot;
  long long bigB;
  char   s[512];
  char  *pargs[20], *rest[20];
  int    pjobs, nrest, k;

  pari_init(2000000, 500000);
  argc = par_jobs_option(argc, argv, pargs, &pjobs);   /* -jN */
  nrest = 0;
  for(k=0; k<argc; k++) {
    if(k > 0 && strcmp(pargs[k], "-s") == 0) squareonly = 1;
    else rest[nrest++] = pargs[k];
  }
  if(nrest != 2 || atoll(rest[1]) < 1) {
    printf("usage: %s [-jN] [-s] B\n"
           "  totally real fields of degree %d with discriminant <= B\n"
           "  -jN  use N processes\n"
           "  -s   keep only fields whose discriminant is a perfect square\n",
           pargs[0], DEG);
    exit(1);
  }

  welcome();
  myinitprimes();
  initmygens();
  if(squareonly) sq_setup();

  bigB = atoll(rest[1]);
  bigbgen = stoi((long) bigB);
  verbosity = 1;

  fact[0] = 1;
  for(k=1; k<=DEG; k++) fact[k] = fact[k-1]*k;
  cf[DEG] = 1;
  for(k=0; k<10; k++) { arc[k] = 1; shift[k] = 0; }

  sprintf(s, "%s_%lld", squareonly ? "TRSquare" : "TRSearch", bigB);
  snprintf(outfilename, sizeof(outfilename), "%.200s", s);

  afile = fopen("Current", "w");
  fprintf(afile, "Working on %s.\n", outfilename);
  fclose(afile);

  outfile = fopen(s, "a");
  if(outfile == NULL) { perror(s); exit(1); }

  polcount = postdisc = t2count = discfcount = 0;
  gettimeofday(&time1, &tzn);

  /* Hunter's bound, with the slack of cN-stand */
  Dmaxroot = 1.0 * bigB/(1.0*DEG);
  Dmaxroot *= mygamma[DEG-1];
  Dmaxroot = fpow(Dmaxroot, 1.0/(DEG-1));
  Dmaxroot += HUNTER_SLACK;
  printf("Totally real, degree %d, disc <= %lld%s, %d process%s\n",
         DEG, bigB, squareonly ? ", square discriminant" : "",
         pjobs, pjobs > 1 ? "es" : "");
  printf("Max T2 (A1 = 0) is %5.3f\n", Dmaxroot);

  ltop = avma;

  if(par_start(pjobs, outfile)) {   /* serial run or parallel worker */
    for(A1=0; 2*A1<=DEG; A1++) {
      set_coef(1, A1);
      t2max = Dmaxroot1;
      find_roots(1);
      level(2);
    }
    t2count = nrejsturm;            /* passed on to the parent in t2count */
    qf_rejected = sq_rej;           /* and this in qf_rejected */
    par_finish(outfile);            /* a worker exits here */
  }
  nrejsturm = t2count;
  sq_rej = qf_rejected;

  gettimeofday(&time2, &tzn);
  timetot = time2.tv_sec - time1.tv_sec;
  fclose(outfile);

  if(squareonly)
    printf("Modular square test: %ld of %ld polynomials rejected (%ld moduli, %ld monomials in the discriminant).\n",
           sq_rej, polcount, (long) NSQMODS, sq_nmono);
  printf("Searched %ld polynomials, %ld past the discriminant filter, "
         "%ld fields written to %s, %ld seconds.\n",
         polcount, postdisc, discfcount, s, timetot);
  if(nrejsturm)
    printf("(%ld irreducible polynomials inside the bounds were not totally real)\n",
           nrejsturm);

  afile = fopen("Current", "a");
  fprintf(afile, "Done.\n");
  fclose(afile);
  return 0;
}
