read("count4.gp");
\\ solvable.gp -- all number fields of prime degree p with solvable Galois group
\\ that are unramified outside a finite set S of primes, by Kummer descent.
\\
\\ Usage:   read("solvable_prime_fields/claude/solvable.gp");
\\          R = solvfields(7, [2,3,7]);   \\ R = [[pol, d], ...], sorted by |disc|
\\ pol is a polredabs defining polynomial; Gal(pol) = C_p : C_d (order p*d, d | p-1).
\\ solvfields(p, S, 1) also verifies every polynomial (nfdisc supported on S,
\\ polgalois order p*d, pairwise distinct polredabs).  solvfields(p, S, 2): certify
\\ Sel(M) for the auxiliary fields M (bnfcertify(,1) + p-saturation of the T'-units); otherwise GRH is assumed
\\ there.  Flags add.  Timings of the stages are printed and kept in the global SVTIME.
\\
\\ Theory (p odd).  A degree-p field F has solvable Galois group iff that group is
\\ C_p : C_d <= AGL(1,p), d | p-1.  Let L be the Galois closure, K = L^{C_p} (cyclic of
\\ degree d), chi : Gal(K/Q) -> F_p^* the (faithful) action on Gal(L/K), omega the
\\ mod p cyclotomic character, psi = omega/chi, all viewed as characters of
\\ (Z/f)^*, where f = prod_{q in S u {p}} q^{a_q} (a_p = 1, a_q = 1+v_q(p-1) for odd
\\ q != p, a_2 = 2+v_2(p-1)); f is divisible by the conductor of every cyclic field of
\\ degree | p-1 unramified outside S u {p}.  Kummer theory over K(zeta_p) and descent
\\ along the group of order prime to p (Hilbert 90 and H^1 = 0 on mu_p) give a bijection
\\   {F up to isomorphism} <-> {(psi, line in Sel(M_psi)[psi] satisfying (*))}
\\ where psi runs over Hom((Z/f)^*, F_p^*) (with psi = omega on (Z/p)^* if p is not in S,
\\ which says K is unramified at p), M = M_psi is the cyclic field cut out by ker psi,
\\ and Sel(M) = {a in M^*/M^*p : p | v_P(a) for all P not above S u {p}}, the
\\ psi-eigenspace taken for a generator sigma of Gal(M/Q).  (*) is empty if p in S;
\\ otherwise M has one prime P above p, with e = p-1, and (*) says that a is a unit
\\ times a p-th power at P whose unit part is a p-th power modulo P^c, c = p if
\\ zeta_p is in M_P and c = p+1 otherwise (if zeta_p is not in M_P an unramified Kummer
\\ class over M_P(zeta_p) coming from M_P is trivial; if it is, unramified classes
\\ are U^(pe/(p-1)) mod p-th powers).  From a line <alpha>, sigma(alpha) = alpha^k gamma^p,
\\ k = psi(sigma); E = M(t), t^p = alpha, has the automorphism t -> gamma t^k extending
\\ sigma (unique when zeta_p is not in M), its fixed field F' is a conjugate of F and
\\ F' = Q(theta), theta = sum_j sigma~^j(t) (a sum of distinct powers of t, so not in M);
\\ the relative characteristic polynomial of theta over M has rational coefficients
\\ and defines F.
\\ p = 2 is done directly (quadratic fields).

\\ (set default(parisizemax, ...) in the calling script: changing it inside read() aborts the read)

\\ ---------- small helpers ----------
sv_order(k, p) = znorder(Mod(k, p));
\\ order of the S-class group Cl(M)/<classes of the primes in T>
sv_hS(bnf, T) = if(#bnf.cyc == 0, 1, matdet(mathnf(matconcat([matdiagonal(bnf.cyc), Mat(vector(#T, i, bnfisprincipal(bnf, T[i], 0)))]))));
\\ nonzero vectors of F_p^e up to scalars (first nonzero coordinate 1)
sv_lines(e, p) =
{
  my(L = List());
  forvec(v = vector(e, i, [0, p-1]), my(i = 1); while(i <= e && v[i] == 0, i++);
    if(i <= e && v[i] == 1, listput(L, v)));
  Vec(L);
}
\\ apply an automorphism to a compact (famat) or plain element
sv_galfamat(nf, tau, z) =
{
  if(type(z) != "t_MAT", nfgaloisapply(nf, tau, z),
    matconcat([apply(e -> nfgaloisapply(nf, tau, nfalgtobasis(nf, e)), z[,1]), z[,2]]));
}
\\ real vector of embeddings (r1 real values, then Re, Im of the r2 complex ones)
sv_toreal(r1, r2, e) = my(v = vector(r1, i, real(e[i]))); for(i = 1, r2, v = concat(v, [real(e[r1+i]), imag(e[r1+i])])); v;
\\ log vector: log|z|_v at archimedean places (complex doubled), v_P(z) log N(P) for P in T
sv_logvec(nf, z, T) =
{
  if(type(z) == "t_MAT", return(sum(i = 1, #z[,1], z[i,2] * sv_logvec(nf, z[i,1], T))));
  my(em = nfeltembed(nf, z), r1 = nf.r1);
  concat(vector(#em, i, if(i <= r1, 1, 2) * log(abs(em[i]))),
         vector(#T, i, nfeltval(nf, z, T[i]) * log(idealnorm(nf, T[i]))));
}
\\ prod_i gens[i]^e[i] as one compact representation (equal factors merged)
sv_famatprod(nf, gens, e) =
{
  my(m = Map(), F, E);
  for(i = 1, #gens, if(e[i], my(g = gens[i]);
    if(type(g) != "t_MAT", g = matrix(1, 2, a, b, if(b == 1, g, 1)));
    for(l = 1, #g[,1], my(z = nfalgtobasis(nf, g[l,1]), v = 0);
      mapisdefined(m, z, &v); mapput(m, z, v + e[i] * g[l,2]))));
  F = Mat(m); if(#F == 0, return(matrix(1, 2, a, b, 1)));
  my(idx = [i | i <- [1..#F[,1]], F[i,2] != 0]);
  if(#idx == 0, matrix(1, 2, a, b, 1), matrix(#idx, 2, a, b, F[idx[a], b]));
}
sv_logs(v, digits) = localprec(digits); log(1. * v);
sv_nfprec(nf, digits) = localprec(digits); nfnewprec(nf);
sv_logvecp(nf, z, T, digits) = localprec(digits); sv_logvec(nf, z, T);
\\ w -> w + p m (free coordinates) making prod gens^w small (Babai rounding in p*Log)
sv_reduce(LV, w, p) =
{
  my(r = #LV, w2 = Vec(w));
  if(r == 0, return(w2));
  my(B = p * LV, U = qflll(B), Br = B * U, t = LV * vector(r, j, w[j])~, m);
  m = round(matsolve(Br~ * Br, -(Br~ * t)));
  for(j = 1, r, w2[j] += p * (U * m)[j]);
  w2;
}
\\ complex logarithms of the embeddings of a compact or plain element, at the
\\ precision of nfh (1. * : nfeltembed of a rational factor is exact)
sv_logemb(nfh, z) =
{
  if(type(z) != "t_MAT", return(log(1. * nfeltembed(nfh, z))));
  sum(i = 1, #z[,1], z[i,2] * log(1. * nfeltembed(nfh, z[i,1])));
}
\\ the log-embedding matrix of the generators at precision prec + (digits of the
\\ largest exponent sum in their compact forms): exp(LE * w) is then prod gens^w with
\\ relative precision about prec digits (the exponents of compact S-units can be huge)
sv_LE(bnf, gens, prec) =
{
  my(ex = vecmax(concat([0], vector(#gens, j, my(g = gens[j]); if(type(g) == "t_MAT", vecsum(abs(g[,2])), 1)))));
  localprec(prec + ceil(log(1 + 1000 * ex) / log(10)) + 10);
  my(nfh = nfnewprec(bnf.nf), cache = Map(), LE = matrix(nfh.r1 + nfh.r2, #gens));
  \\ the compact representations of different generators share most factors: log each
  \\ distinct factor once
  for(j = 1, #gens, my(g = gens[j]);
    if(type(g) != "t_MAT", g = matrix(1, 2, a, b, if(b == 1, g, 1)));
    for(i = 1, #g[,1], my(zb = nfalgtobasis(bnf, g[i,1]), l);
      if(!mapisdefined(cache, zb, &l), l = log(1. * nfeltembed(nfh, zb)); mapput(cache, zb, l));
      for(t = 1, #l, LE[t, j] += g[i,2] * l[t])));
  [LE, nfh];
}
\\ alpha = N^p * prod gens^w * b^p as an algebraic integer (currently always b = 1; h is
\\ the size of prod gens^w in digits).
\\ alpha is recovered from its embeddings exp(LE w + p log b) (LEh = sv_LE(...)) and
\\ checked: rounding error < 2^-40, and |norm(alpha)| equal to the norm predicted by the
\\ exact valuations of the factored form.  Returns [alpha, rational primes dividing N b],
\\ or [0, 0, 0, 0] if the precision was too low.
sv_makealpha(bnf, gens, w, T, p, h, LEh) =
{
  my(r = #gens, A = List(), Pl, N = 1, v, c, err, al, nrm, ps, b = 1, nfh = LEh[2]);
  for(j = 1, r, if(w[j], my(g = gens[j]);
    if(type(g) != "t_MAT", listput(A, [g, w[j]]), for(i = 1, #g[,1], listput(A, [g[i,1], g[i,2] * w[j]])))));
  if(#A == 0, return([1, [], 1, 1]));
  A = Vec(A); A = matrix(#A, 2, i, j, A[i][j]);
  Pl = T;
  \\ (no reduction by p-th powers of non-S-units: idealredmodpower does not accept elements
  \\ in factored form (p = 17, S = {2,17}); a large alpha just costs precision, b = 1)
  ps = Set(apply(P -> P.p, Pl));
  my(vals = vector(#Pl, i, nfeltval(bnf, A, Pl[i])));
  for(i = 1, #ps, my(q = ps[i], a = 0);
    for(j = 1, #Pl, if(Pl[j].p == q, a = max(a, ceil(-vals[j] / (p * Pl[j].e)))));
    N *= q^a);
  localprec(precision(nfh.roots[1]));
  v = LEh[1] * vector(r, j, w[j])~;
  if(b != 1, v += p * sv_logemb(nfh, b)~);
  v = exp(v~) * N^p;
  my(Emb = Mat(vector(#nfh.zk, j, sv_toreal(nfh.r1, nfh.r2, nfeltembed(nfh, nfh.zk[j]))~)));
  c = round(matsolve(Emb, sv_toreal(nfh.r1, nfh.r2, v)~), &err);
  if(err > -40, return([0, 0, 0, 0]));
  al = nfbasistoalg(bnf, c);
  if(al == 0, return([0, 0, 0, 0]));
  nrm = N^(p * poldegree(bnf.pol)) * prod(i = 1, #Pl, idealnorm(bnf, Pl[i])^vals[i]);
  if(abs(nfeltnorm(bnf, al)) != nrm, return([0, 0, 0, 0]));
  [lift(al), ps, N, b];
}
\\ gamma with sigma(alpha) = alpha^k gamma^p, mu_p not in M: alpha = N^p prod gens^w b^p,
\\ sigma(gens) = gens^Tint * (torsion), so gamma = N^(1-k) prod gens^((Tint-k)w/p)
\\ sigma(b)/b^k * (a root of unity of M, of order prime to p).  Evaluated from the
\\ log-embeddings, rounded (denominator from that of beta = gamma^p), and checked
\\ exactly: gamma^p = beta.  Returns 0 on failure.
sv_gamma(bnf, al, beta, sigma, k, w, Tint, LEh, Nn, b, p) =
{
  my(nfh = LEh[2], m, lg, d, D = 1, Emb, mu = bnf.tu, fa);
  m = (Tint - k) * vector(#Tint, j, w[j])~ / p;
  if(denominator(m) != 1, error("w is not in the eigenspace"));
  localprec(precision(nfh.roots[1]));
  lg = LEh[1] * m + (1 - k) * log(Nn) * vector(nfh.r1 + nfh.r2, i, 1)~;
  if(b != 1, lg += sv_logemb(nfh, nfgaloisapply(bnf, sigma, b))~ - k * sv_logemb(nfh, b)~);
  d = denominator(nfalgtobasis(bnf, beta)); fa = factor(d);
  for(i = 1, #fa[,1], D *= fa[i,1]^ceil(fa[i,2] / p));
  Emb = Mat(vector(#nfh.zk, j, sv_toreal(nfh.r1, nfh.r2, nfeltembed(nfh, nfh.zk[j]))~));
  my(zl = 1. * nfeltembed(nfh, mu[2]), v0 = exp(lg~) * D);
  for(j = 0, mu[1] - 1, my(c, err, g);
    c = round(matsolve(Emb, sv_toreal(nfh.r1, nfh.r2, vector(#v0, t, v0[t] * zl[t]^j))~), &err);
    if(err > -20, next);
    g = nfbasistoalg(bnf, c) / D;
    if(nfeltpow(bnf, g, p) == nfalgtobasis(bnf, beta), return(lift(g))));
  0;
}
\\ unconditional proof that gens (compact T'-units, + zeta if p | w) are independent in
\\ M^*/M^*p, via p-th power residue characters at degree-1 primes outside T'.  Together with
\\ bnfcertify(bnf, 1) this makes Sel(M) GRH-free.  Returns 0 if it gives up.
sv_psatred(nf, x, P, modP) =
{
  if(type(x) != "t_MAT", x = Mat([x, 1]));
  for(i = 1, #x~, if(nfeltval(nf, x[i,1], P), return(0)));
  prod(i = 1, #x~, nfmodpr(nf, x[i,1], modP)^x[i,2]);
}
sv_psat(bnf, gens, Tp, p, maxprimes = 2000) =
{
  my(nf = bnf.nf, n = #gens, M = [;], q = 1, cnt = 0, Tq = Set(apply(P -> P.p, Tp)));
  while(matrank(M) < n,
    if(cnt++ > maxprimes, return(0));
    q = nextprime(q + 1);
    if(q % p != 1 || nf.disc % q == 0 || setsearch(Tq, q), next);
    foreach(idealprimedec(nf, q), P,
      if(P.f != 1, next);
      my(modP = nfmodprinit(nf, P), zp, row = vector(n), ok = 1);
      zp = ffprimroot(nfmodpr(nf, 1, modP))^((q - 1) / p);
      for(i = 1, n, my(r = sv_psatred(nf, gens[i], P, modP));
        if(r == 0, ok = 0; break);
        row[i] = fflog(r^((q - 1) / p), zp, p));
      if(ok, M = matconcat([M; Mod(row, p)]))));
  1;
}

\\ the S-unit exponent coordinates (mod p) of z; torsion coordinate appended if tor
sv_coord(bnf, U, rf, z, tor, p) =
{
  my(v = bnfisunit(bnf, z, U)); if(#v == 0, error("not an S-unit"));
  vector(rf + tor, j, if(j <= rf, v[j], lift(v[rf + 1])));
}

\\ the cyclic field fixed by the subgroup H (HNF on G.gen) of G = znstar(f,1), [G:H] = n.
\\ Small f: galoissubcyclo.  Large f (its cost grows with phi(f)): bnrclassfield over Q (Kummer theory over Q(zeta_l), l | n), the ray class
\\ group mod f*oo being (Z/f)^* via (a) -> a for a > 0, so the Artin symbol of (l) is
\\ Frob_l, the same convention as galoissubcyclo.
sv_fieldM(G, H, n, f) =
{
  if(f < 10^5, return(galoissubcyclo(G, H)));
  my(bnr = bnrinit(bnfinit(y, 1), [f, [1]]), M);
  \\ the isomorphism G -> Cl_{f oo}(Q), a -> class of (a) (a > 0), on the generators of G
  M = Mat(vector(#G.gen, i, bnrisprincipal(bnr, G.gen[i], 0)));
  bnrclassfield(bnr, mathnf(matconcat([M * H, matdiagonal(bnr.cyc)])), 2);
}
\\ ---------- p = 2 ----------
sv_quadratic(S) =
{
  my(g = concat([-1], S), L = List());
  forvec(e = vector(#g, i, [0, 1]), my(d = core(prod(i = 1, #g, g[i]^e[i])));
    if(d == 1 || (d % 4 != 1 && !setsearch(Set(S), 2)), next);
    listput(L, [polredabs(x^2 - d), 1]));
  Vec(L);
}

\\ ---------- main ----------
solvfields(p, S, flag = 0) =
{
  SVNFROOTS = 0; SVFALLBACK = 0;
  my(cnts,cntsfound,cnts2);
  my(t0 = getabstime(), tm = vector(6), Sp, f, G, cyc, g0, om, psis, groups, res = List(), nM = 0, stats = List());
  S = Set(S);
  if(p == 2, res = vecsort(sv_quadratic(S), r -> [abs(nfdisc(r[1])), Vec(r[1])]); SVTIME = [getabstime() - t0];
    if(bitand(flag, 1), sv_verify(res, p, S)); return(res));
  Sp = setunion(S, [p]);
  f = prod(i = 1, #Sp, my(q = Sp[i]); q^if(q == p, 1, if(q == 2, 2 + valuation(p - 1, 2), 1 + valuation(p - 1, q))));
  G = znstar(f, 1); cyc = G.cyc;
  g0 = znprimroot(p);
  \\ omega on the generators, as exponents of g0
  om = vector(#cyc, i, znlog(Mod(G.gen[i], p), g0));
  \\ all psi : (Z/f)^* -> F_p^*, psi(gen_i) = g0^e_i
  my(steps = vector(#cyc, i, (p - 1) / gcd(cyc[i], p - 1)), xp = 0, phi = 0);
  if(!setsearch(S, p),
    xp = lift(chinese(Mod(lift(g0), p), Mod(1, f / p)));      \\ generator of the inertia (Z/p)^*
    phi = lift(chinese(Mod(1, p), Mod(p, f / p))));            \\ Frobenius at p on the prime-to-p part
  psis = List();
  forvec(a = vector(#cyc, i, [0, (p - 1) / steps[i] - 1]), my(e = vector(#cyc, i, a[i] * steps[i]));
    if(xp && (e * znlog(xp, G)) % (p - 1) != 1, next);    \\ psi(xp) = g0: psi = omega on (Z/p)^*
    listput(psis, e));
  \\ group psi by kernel
  groups = Map();
  for(i = 1, #psis, my(e = psis[i], K0, H);
    K0 = matkerint(Mat(concat(e, [p - 1])));
    H = mathnf(matconcat([K0[1..#cyc, ], matdiagonal(cyc)]));
    if(mapisdefined(groups, H), my(l = mapget(groups, H)); listput(l, e); mapput(groups, H, l),
                                 mapput(groups, H, List([e]))));
  my(GM = Mat(groups));
  for(iM = 1, #GM[,1], my(H = GM[iM,1], pl = Vec(GM[iM,2]), n = matdet(H), polM, bnf, ta, xa, ell, tau,
       T, Tp, extra, q, U, rf, tor, gens, LV, Tm, rows, P, c, sigma, tt);
    tt = getabstime();
    nM++;
    \\ the field M
    polM = if(n == 1, y, subst(polredbest(sv_fieldM(G, H, n, f)), x, y));
    \\ a generator x_a of (Z/f)^*/H: psi(x_a) of order n
    xa = 0; forstep(j = 1, oo, 1, if(gcd(j, f) == 1 && sv_order(g0^((pl[1] * znlog(j, G)) % (p - 1)), p) == n, xa = j; break));
    \\ sigma = Frob_ell, ell = x_a mod f
    if(n == 1, sigma = y,
      my(D = poldisc(polM), cj = nfgaloisconj(polM), dc);
      dc = lcm(apply(s -> denominator(content(s)), Vec(cj)));   \\ (content of a t_POL vector is a t_POL)
      ell = xa; while(!isprime(ell) || D % ell == 0 || dc % ell == 0, ell += f);
      my(fr = lift(Mod(Mod(1, ell) * y, Mod(1, ell) * polM)^ell));
      sigma = [s | s <- cj, Mod(1, ell) * s == fr];
      if(#sigma != 1, error("Frobenius not identified for ", polM));
      sigma = sigma[1]);
    bnf = bnfinit(polM, 1);
    \\ flag 2 (GRH-free): the factor base generates Cl(M) (Zimmert bound), so the true
    \\ class group is a quotient of bnf.cyc; with the p-saturation test below this proves Sel(M)
    if(bitand(flag, 2) && n > 1 && !bnfcertify(bnf, 1), error("bnfcertify(,1) failed for ", polM));
    tm[1] += getabstime() - tt; tt = getabstime();
    \\ T = primes above S u {p}; enlarge until the T'-class number is prime to p
    T = concat([idealprimedec(bnf, q) | q <- Sp]);
    Tp = T; extra = []; q = 1;
    if(bnf.no % p == 0, while(sv_hS(bnf, Tp) % p == 0, q = nextprime(q + 1); if(setsearch(Sp, q), next);
      my(Pq = idealprimedec(bnf, q)); Tp = concat(Tp, Pq); extra = concat(extra, Pq)));
    U = bnfunits(bnf, Tp); rf = #U[1] - 1;
    tor = (bnf.tu[1] % p == 0);
    gens = vector(rf, j, U[1][j]);
    if(tor, gens = concat(gens, [nfalgtobasis(bnf, bnf.tu[2])]));
    if(bitand(flag, 2) && n > 1 && #gens && !sv_psat(bnf, gens, Tp, p), error("T'-units not p-saturated for ", polM));
    my(r = rf + tor, LEh, LEprec = 100, Tint, nfL, VT, lN, ra, perA, perT);
    \\ log-embeddings of the generators (one high-precision computation per M, reused for
    \\ every alpha); the log vectors LV: log|g|_v (complex places doubled), v_P(g) log N(P)
    LEh = sv_LE(bnf, gens, LEprec); nfL = LEh[2];
    VT = matrix(#Tp, rf, i, j, nfeltval(bnf, gens[j], Tp[i]));
    lN = sv_logs(vector(#Tp, i, idealnorm(bnf, Tp[i])), precision(nfL.roots[1]));
    ra = nfL.r1 + nfL.r2;
    LV = matrix(ra + #Tp, rf, i, j, if(i <= ra, if(i <= nfL.r1, 1, 2) * real(LEh[1][i, j]), VT[i - ra, j] * lN[i - ra]));
    \\ LLL-reduce the S-unit basis in the log lattice.  bnfunits may return generators of
    \\ enormous height (log|g| ~ 10^35 when the class group is large, e.g. p = 17, S = {2}),
    \\ which makes the log lattice numerically degenerate.  New generators g'_j =
    \\ prod_i g_i^Ul[i,j] (Ul unimodular), in compact form; coordinates c on the old
    \\ generators become Ul^-1 c on the new ones.
    my(Ul = matid(rf), Uli = matid(rf));
    if(rf, Ul = qflll(LV);
      if(matsize(Ul) != [rf, rf] || abs(matdet(Ul)) != 1, error("LLL of the S-unit lattice failed"));
      Uli = Ul^-1;
      gens = concat(vector(rf, j, sv_famatprod(bnf, gens[1..rf], Ul[, j])), gens[rf + 1..#gens]);
      LV = LV * Ul; VT = VT * Ul;
      LEh = [LEh[1] * matrix(#gens, #gens, i, j, if(i <= rf && j <= rf, Ul[i, j], i == j)), LEh[2]]);
    \\ the matrix of sigma on O_T'^* / p.  If mu_p is not in M: sigma permutes the archimedean
    \\ places and the primes of T', so Log(sigma g) is a permutation of Log(g), and the
    \\ exponents of sigma(g_j) on the free generators are read off in the log lattice (exact
    \\ up to torsion, which is a p-th power).  If mu_p is in M, bnfisunit (torsion coordinate).
    if(tor, Tm = Mod(1, p) * Mat(vector(r, j, my(c = sv_coord(bnf, U, rf, sv_galfamat(bnf, sigma, gens[j]), tor, p)~);
                                                 concat(Uli * c[1..rf], c[rf + 1..r]))),
      \\ |g(sigma(rho_t))| = |g(rho_perA[t])|, and sigma(Tp[i]) = Tp[perT[i]]
      perA = vector(ra, t, my(z = subst(lift(sigma), y, nfL.roots[t]), dd = vector(ra, u, min(abs(z - nfL.roots[u]), abs(z - conj(nfL.roots[u])))));
        vecsort(dd, , 1)[1]);
      perT = vector(#Tp, i, my(Q = idealhnf(bnf, nfgaloisapply(bnf, sigma, Tp[i])), jj = [j | j <- [1..#Tp], idealhnf(bnf, Tp[j]) == Q]);
        if(#jj != 1, error("T' is not Galois stable")); jj[1]);
      if(#Set(perA) != ra || #Set(perT) != #Tp, error("not a permutation"));
      my(G2 = LV~ * LV);
      Tint = Mat(vector(r, j, my(Ls = vector(ra + #Tp), cr, cc);
        for(t = 1, ra, Ls[t] = LV[perA[t], j]);
        for(i = 1, #Tp, Ls[ra + perT[i]] = LV[ra + i, j]);
        Ls = Ls~;
        cr = matsolve(G2, LV~ * Ls); cc = round(cr);
        if(normlp(cr - cc) > 10^-6 || normlp(LV * cc - Ls) > 10^-6,     \\ ill-conditioned: exact fallback
          SVFALLBACK++; cc = Uli * sv_coord(bnf, U, rf, sv_galfamat(bnf, sigma, gens[j]), 0, p)~);
        cc)); Tm = Mod(1, p) * Tint);
    \\ linear conditions: valuations at the extra primes, and (*) at p if p is not in S
    rows = matrix(0, r);
    if(#extra, rows = matconcat([rows; Mod(1, p) * matrix(#extra, r, i, j, if(j <= rf, nfeltval(bnf, gens[j], extra[i]), 0))]));
    if(!setsearch(S, p),
      P = idealprimedec(bnf, p); if(#P != 1 || P[1].e != p - 1, error("unexpected splitting of p in ", polM));
      P = P[1];
      my(pi = P.gen[2], bid, ucond, vrow);
      if(nfeltval(bnf, pi, P) != 1, error("uniformizer"));
      c = if((pl[1] * znlog(phi, G)) % (p - 1) == 0, p, p + 1);     \\ psi(Frob) = 1 <=> zeta_p in M_P
      bid = idealstar(bnf, idealpow(bnf, P, c), 1);
      vrow = vector(r, j, if(j <= rf, nfeltval(bnf, gens[j], P), 0));
      ucond = [i | i <- [1..#bid.cyc], bid.cyc[i] % p == 0];
      my(lg = vector(r, j, my(z = gens[j], m);
        if(type(z) != "t_MAT", z = matrix(1, 2, a, b, if(b == 1, z, 1)));
        if(vrow[j], m = #z[,1]; z = matrix(m + 1, 2, a, b, if(a <= m, z[a, b], if(b == 1, pi, -vrow[j]))));
        ideallog(bnf, z, bid)));
      rows = matconcat([rows; Mod(1, p) * Mat(vrow); Mod(1, p) * matrix(#ucond, r, i, j, lg[j][ucond[i]])]));
    tm[2] += getabstime() - tt;
    \\ each psi in this kernel class
    for(ip = 1, #pl, my(e = pl[ip], k = lift(g0^((e * znlog(xa, G)) % (p - 1))), Ek, dim, chiexp, d);
      tt = getabstime();
      Ek = matker(Tm - k);
      if(#Ek && matsize(rows)[1], Ek = Ek * matker(rows * Ek));
      dim = #Ek;
      chiexp = (om - e) % (p - 1); d = (p - 1) / gcd(concat(chiexp, [p - 1]));
      listput(stats, [polM, k, d, dim]);
      tm[3] += getabstime() - tt;
      if(dim == 0, next);
      my(LN = sv_lines(dim, p));
      for(il = 1, #LN, my(w, al, N, Nn, bb, prec, beta, gam, mon, th, fpol, ok, listP);
        tt = getabstime();
        w = lift(Ek * LN[il]~);
        w = sv_reduce(LV, Vec(w), p);
        my(h = if(#LV, vecmax(abs(LV * vector(rf, j, w[j])~)) / log(10), 0));   \\ size of prod gens^w (digits)
        prec = max(LEprec, 100); al = 0;
        while(al == 0, if(LEprec < prec, LEh = sv_LE(bnf, gens, prec); LEprec = prec);
          [al, N, Nn, bb] = sv_makealpha(bnf, gens, w, Tp, p, h, LEh);
          if(al == 0, prec *= 2; if(prec > 10^5, error("precision"))));
        tm[4] += getabstime() - tt; tt = getabstime();
        \\ gamma: sigma(alpha) = alpha^k gamma^p
        if(n == 1, gam = 1,
          beta = nfeltdiv(bnf, nfgaloisapply(bnf, sigma, al), nfeltpow(bnf, al, k));
          gam = if(tor, 0, sv_gamma(bnf, al, beta, sigma, k, w, Tint, LEh, Nn, bb, p));
          if(gam == 0, SVNFROOTS++;              \\ fallback (and the case mu_p in M)
            gam = nfroots(bnf, x^p - nfbasistoalg(bnf, beta));
            if(#gam == 0, error("no p-th root: alpha not in the eigenspace?"));
            gam = lift(gam[1])));
        \\ theta = sum_{j<n} sigma~^j(t), with monomials [c, e] = c t^e, t^p = alpha
        ok = 0;
        for(iz = 0, if(tor, p - 1, 0), my(g1 = if(iz, lift(Mod(gam * bnf.tu[2]^(iz * bnf.tu[1] / p), polM)), gam));   \\ gamma zeta_p^iz
          mon = [[1, 1]]; th = 0;
          for(j = 1, n, my(m = mon[#mon], c2, e2);
            c2 = Mod(subst(lift(m[1]), y, Mod(sigma, polM)) * g1^m[2], polM); e2 = k * m[2];
            c2 *= Mod(al, polM)^(e2 \ p); e2 %= p;
            mon = concat(mon, [[lift(c2), e2]]));
          if(mon[n + 1] == [1, 1], ok = 1; break));
        if(!ok, error("sigma~ is not an automorphism of order n"));
        th = sum(j = 1, n, Mod(mon[j][1], polM) * x^mon[j][2]);
        fpol = rnfcharpoly(bnf, x^p - Mod(al, polM), th, x);
        fpol = apply(cf -> my(z = lift(cf)); if(poldegree(z, y) > 0, error("charpoly not rational: M=", polM, " k=", k, " alpha=", al, " gamma=", gam, " theta=", th, " f=", fpol)); z, fpol);
        tm[5] += getabstime() - tt; tt = getabstime();
        listP = setunion(Sp, N);                \\ N = primes dividing alpha (possible index primes)
        fpol = polredbest([fpol, listP]);
        fpol = polredabs([fpol, Sp]);
        tm[6] += getabstime() - tt;
        listput(res, [fpol, d]))));
  res = Vec(res);
  res = vecsort(res, r -> [abs(nfdisc([r[1], Sp])), Vec(r[1])]);
  SVTIME = concat([getabstime() - t0], tm);
  SVSTATS = Vec(stats);
  printf("p = %d, S = %s: %d fields; %d fields M, %d characters psi; time %.2f s\n", p, S, #res, nM, #psis, (getabstime() - t0) / 1000.);
  if(SVNFROOTS || SVFALLBACK, printf("  (nfroots used for gamma %d times; bnfisunit fallback for sigma %d times)\n", SVNFROOTS, SVFALLBACK));
  printf("  [bnfinit M %.2f, S-units %.2f, eigenspaces %.2f, alpha %.2f, charpoly %.2f, polred %.2f]\n", tm[1]/1000., tm[2]/1000., tm[3]/1000., tm[4]/1000., tm[5]/1000., tm[6]/1000.);
  if(bitand(flag, 1), sv_verify(res, p, S));
  sv_verify(res,p,S);
  cnts=solvcount(p,S,flag==0);
  cnts2=[];
  for(j=1,#cnts,if(cnts[j]>0,cnts2=concat(cnts2,[[j,cnts[j]]])));
  cntsfound=vecsort(mult(apply(z->z[2],res)),[1]);
  if(cntsfound!=cnts2, print("Counts do not match: ",cnts2," vs ",cntsfound);error("BAD"));
  apply(z->z[1], res);
}

\\ independent checks: degree, irreducible, disc supported on S, Galois group order p*d,
\\ polredabs canonical, pairwise distinct
sv_verify(res, p, S) =
{
  my(bad = 0, t = getabstime());
  for(i = 1, #res, my(f = res[i][1], d = res[i][2], D, G);
    if(poldegree(f) != p || !polisirreducible(f), bad++; print("  not degree p irreducible: ", f); next);
    D = nfdisc(f);
    if(#setminus(Set(factor(abs(D))[,1]~), Set(S)), bad++; print("  ramified outside S: ", f));
    /*if(p > 2 && p <= 11, G = galt(f); if(G[1] != p * d, bad++; print("  Galois group ", G, " != order ", p * d, ": ", f))); */
    if(polredabs(f) != f, bad++; print("  not polredabs: ", f)));
  if(#Set(apply(r -> r[1], res)) != #res, bad++; print("  duplicates"));
  printf("  verify: %d problems (%.2f s)\n", bad, (getabstime() - t) / 1000.);
  if(bad, error("verify failed");error("Verify failed"));
  bad;
}

default(realprecision,100);

onep(d,N,fl)={
  my(L=[]);
  forprime(pp=2,N,L=concat(L,solvfields(d,[pp],fl)));
  grh=if(fl==0, "", "grh");
  write(Str("Onepd", d,"B",N,grh),L);
  return(#L);
}

onec(d,N,fl)={
  forprime(pp=2,N,solvcount(d,[pp],fl));
  return(1);
}

twop(d,N,fl)={
  my(L=[]);
  forprime(pp=2,N,forprime(k=pp+1,N,L=concat(L,solvfields(d,[pp,k],fl))));
  grh=if(fl==0, "", "grh");
  write(Str("Twopd", d,"B",N,grh),L);
  return(#L);
}

threep(d,N,fl)={
  my(L=[]);
  forprime(pp=2,N,
    forprime(k=pp+1,N,
      forprime(ll=k+1,N,
        L=concat(L,solvfields(d,[pp,k,ll],fl)))));
  grh=if(fl==0, "", "grh");
  write(Str("Threepd", d,"B",N,grh),L);
  return(#L);
}

twoc(d,N,fl)={
  forprime(pp=2,N,forprime(k=pp+1,N,solvcount(d,[pp,k],fl)));
  return(1);
}

threec(d,N,fl)={
  forprime(pp=2,N,
    forprime(k=pp+1,N,
      forprime(ll=k+1,N,
        solvcount(d,[pp,k,ll],fl))));
  return(1);
}

