\\ count.gp -- count the solvable degree-p fields unramified outside S by ray class
\\ groups, independently of the Kummer construction in solvable.gp.
\\
\\ Usage:  read("solvable_prime_fields/claude/count.gp"); solvcount(p, S)   (p odd)
\\ returns [N_1, ..., N_{p-1}], N_d = number of fields with Galois group C_p : C_d
\\ (0 unless d | p-1).  solvcount(p, S, 1) certifies, without GRH, the p-parts that
\\ are used (bnfcertify(,1) + p-saturation, sc_psat); cheaper than a full bnfcertify.
\\
\\ A field F with group C_p : C_d <-> its Galois closure L, a cyclic degree p
\\ extension of the cyclic field K = L^{C_p} of degree d, Galois over Q, with
\\ Gal(K/Q) acting on Gal(L/K) by a faithful character chi.  K is unramified outside S,
\\ L/K is tamely ramified above q in S, q != p, and its conductor exponent above
\\ P | p is at most floor(p e_P/(p-1)) + 1 (the upper ramification break of a cyclic
\\ degree p extension of a p-adic field with absolute ramification index e is at
\\ most p e/(p-1)).  So L is cut out by a quotient Cl_m(K) -> Z/p on which a generator
\\ sigma of Gal(K/Q) acts by k = chi(sigma), m = prod_{q in S, q != p} q *
\\ prod_{P | p} P^{floor(p e_P/(p-1))+1} (no real places: p is odd), and the number of
\\ such L for (K, chi) is (p^r - 1)/(p - 1), r = dim of the k-eigenspace of sigma on
\\ Cl_m(K)/p (semisimple: d is prime to p).  K runs over the kernels of the characters
\\ (Z/f)^* -> F_p^* with f = prod_{q in S} q^{a_q} (a_p = 1, a_q = 1 + v_q(p-1),
\\ a_2 = 2 + v_2(p-1)), the characters with a given kernel being the faithful chi.

\\ the subfield of Q(zeta_f) fixed by H (HNF on G.gen); for large f by class field theory
\\ over Q (the cost of galoissubcyclo grows with phi(f))
sc_field(G, H, f) =
{
  if(f < 10^5, return(galoissubcyclo(G, H)));
  my(bnr = bnrinit(bnfinit(y, 1), [f, [1]]), M);
  M = Mat(vector(#G.gen, i, bnrisprincipal(bnr, G.gen[i], 0)));
  bnrclassfield(bnr, mathnf(matconcat([M * H, matdiagonal(bnr.cyc)])), 2);
}

\\ residue at the degree-1 prime modP of x (nf element or famat); 0 if x is not a P-unit
sc_psatred(nf, x, P, modP) =
{
  if(type(x) != "t_MAT", x = Mat([x, 1]));
  for(i = 1, #x~, if(nfeltval(nf, x[i,1], P), return(0)));
  prod(i = 1, #x~, nfmodpr(nf, x[i,1], modP)^x[i,2]);
}
\\ unconditional proof that the fundamental units, zeta (if p | w) and alpha_j with
\\ (alpha_j) = g_j^{d_j} (p | d_j) are independent in K^*/K^*p: p-th power residue
\\ characters at degree-1 primes q = 1 mod p.  Returns 0 if it gives up.
sc_psat(bnf, p, maxprimes = 2000) =
{
  my(nf = bnf.nf, G = List(), U = bnfunits(bnf)[1], cyc = bnf.cyc, M = [;], q = 1, cnt = 0, n);
  for(i = 1, #U - 1, listput(G, U[i]));
  if(bnf.tu[1] % p == 0, listput(G, bnf.tu[2]));
  for(j = 1, #cyc, if(cyc[j] % p == 0,
    listput(G, bnfisprincipal(bnf, idealpow(nf, bnf.gen[j], cyc[j]), 4)[2])));
  n = #G;
  if(!n, return(1));
  while(matrank(M) < n,
    if(cnt++ > maxprimes, return(0));
    q = nextprime(q + 1);
    if(q % p != 1 || nf.disc % q == 0, next);
    foreach(idealprimedec(nf, q), P,
      if(P.f != 1, next);
      my(modP = nfmodprinit(nf, P), zp, row = vector(n), ok = 1);
      zp = ffprimroot(nfmodpr(nf, 1, modP))^((q - 1) / p);
      for(i = 1, n, my(r = sc_psatred(nf, G[i], P, modP));
        if(r == 0, ok = 0; break);
        row[i] = fflog(r^((q - 1) / p), zp, p));
      if(ok, M = matconcat([M; Mod(row, p)]))));
  1;
}

solvcount(p, S, flag = 0) =
{
  my(t0 = getabstime(), f, G, cyc, g0, steps, chis, groups, GM, cnt = vector(p - 1));
  S = Set(S);
  if(#S == 0, return(cnt));
  f = prod(i = 1, #S, my(q = S[i]); q^if(q == p, 1, if(q == 2, 2 + valuation(p - 1, 2), 1 + valuation(p - 1, q))));
  G = znstar(f, 1); cyc = G.cyc; g0 = znprimroot(p);
  steps = vector(#cyc, i, (p - 1) / gcd(cyc[i], p - 1));
  groups = Map();
  forvec(a = vector(#cyc, i, [0, (p - 1) / steps[i] - 1]), my(e = vector(#cyc, i, a[i] * steps[i]), K0, H);
    K0 = matkerint(Mat(concat(e, [p - 1])));
    H = mathnf(matconcat([K0[1..#cyc, ], matdiagonal(cyc)]));
    if(mapisdefined(groups, H), my(l = mapget(groups, H)); listput(l, e); mapput(groups, H, l),
                                 mapput(groups, H, List([e]))));
  GM = Mat(groups);
  for(iK = 1, #GM[,1], my(H = GM[iK,1], cl = Vec(GM[iK,2]), d = matdet(H), polK, bnf, xa, sigma, m, bnr, ip, A, Am);
    polK = if(d == 1, y, subst(polredbest(sc_field(G, H, f)), x, y));
    xa = 0; forstep(j = 1, oo, 1, if(gcd(j, f) == 1 && znorder(g0^((cl[1] * znlog(j, G)) % (p - 1))) == d, xa = j; break));
    if(d == 1, sigma = y,
      my(D = poldisc(polK), cj = nfgaloisconj(polK), dc = lcm(apply(s -> denominator(content(s)), Vec(cj))), ell = xa, fr);
      while(!isprime(ell) || D % ell == 0 || dc % ell == 0, ell += f);
      fr = lift(Mod(Mod(1, ell) * y, Mod(1, ell) * polK)^ell);
      sigma = [s | s <- cj, Mod(1, ell) * s == fr]; if(#sigma != 1, error("Frobenius")); sigma = sigma[1]);
    bnf = bnfinit(polK, 1);
    \\ GRH-free: the factor base generates Cl(K) (bnfcertify(,1), Zimmert bound) and the p-Selmer
    \\ generators are independent mod p-th powers (sc_psat).  Then Cl(K) -> Cl_true and the unit
    \\ index are of order prime to p, so Cl_m(K)/p (all that is used below) is correct.
    if(flag && d > 1 && !bnfcertify(bnf, 1), error("bnfcertify(,1) failed for ", polK));
    if(flag && d > 1 && !sc_psat(bnf, p), error("p-saturation failed for ", polK));
    \\ the modulus m, as a factorization matrix
    my(Pr = List(), Ex = List());
    for(i = 1, #S, my(q = S[i], dec = idealprimedec(bnf, q));
      for(j = 1, #dec, listput(Pr, dec[j]); listput(Ex, if(q == p, (p * dec[j].e) \ (p - 1) + 1, 1))));
    m = matrix(#Pr, 2, i, j, if(j == 1, Pr[i], Ex[i]));
    bnr = bnrinit(bnf, m, 1);
    ip = [i | i <- [1..#bnr.cyc], bnr.cyc[i] % p == 0];       \\ Cl_m / p = F_p^#ip
    if(#ip == 0, next);
    A = Mat(vector(#bnr.gen, j, bnrisprincipal(bnr, nfgaloisapply(bnf, sigma, bnr.gen[j]), 0)));
    Am = Mod(matrix(#ip, #ip, i, j, A[ip[i], ip[j]]), p);
    for(ic = 1, #cl, my(k = g0^((cl[ic] * znlog(xa, G)) % (p - 1)), r = #matker(Am - k));
      cnt[d] += (p^r - 1) / (p - 1)));
  printf("count p = %d, S = %s: %s (total %d), %.2f s\n", p, S, cnt, vecsum(cnt), (getabstime() - t0) / 1000.);
  cnt;
}
