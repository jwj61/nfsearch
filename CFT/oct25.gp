\\ oct25.gp -- all octic fields F with Galois group 8T25 = 2^3:7 = AGL(1,8) (order 56)
\\ unramified outside a finite set S of primes, with a GRH-free proof of completeness.
\\
\\   oct25(S)        list of polredabs polynomials, complete without GRH (default)
\\   oct25(S, 0)     same, but trust bnfinit (GRH) -- faster
\\   oct25count(S)   just the (certified) number of fields, without constructing them
\\
\\ Method.  Let L be the Galois closure of F, G = Gal(L/Q) = V:C7, V = F_2^3 with C7
\\ acting irreducibly.  Then K = L^V is a cyclic septic field and L/K is a (C2)^3
\\ extension, Galois over Q, on which C7 acts irreducibly.  Conversely, for such K and L
\\ (H^2(C7,V) = 0, coprime orders) Gal(L/Q) = V:C7 and F = L^{C7} is an 8T25 field;
\\ all complements are conjugate (H^1 = 0), so L determines F up to isomorphism and vice
\\ versa (L is the Galois closure of F, K = L^V is determined by L).  Hence
\\
\\   #{F} = sum over cyclic septic K unramified outside S of the number of
\\          sigma-stable H < Cl_m(K) with Cl_m(K)/H = F_2^3 irreducible as C7-module,
\\
\\ where m = prod_{P | q in S} P^{a_P} * (all 7 real places), a_P = 1 for q odd and
\\ a_P = 2 e(P/2) + 1 for q = 2 (the conductor bound for quadratic extensions).
\\
\\ * K: conductor divides f0 = 49^[7 in S] * prod_{q in S, q = 1 mod 7} q (Kronecker-Weber
\\   + the conductor bound of Lemma 3 of grhfree_count.tex): unconditional.
\\ * Cl_m(K)/2 with its sigma-action: unconditional once bnfcertify(bnf,1) holds and
\\   the units, -1 and the alpha_j ((alpha_j) = g_j^{d_j}, d_j even) are independent in
\\   K^*/K^*2 (quadratic residue characters) -- Lemmas 5-8 of grhfree_count.tex with p = 2.
\\ * Cl_m/2 = M_1 + M_f + M_g (f, g the two irreducible cubics over F_2, 7 odd so
\\   semisimple); M_f = S_f^k and the quotients = S_f are counted by (8^k - 1)/7.
\\ Every returned F is checked (polgalois = 8T25, disc supported on S, pairwise
\\ non-isomorphic) and their number equals the certified count, so the list is complete.
\\
\\ Construction of F from H: bnrclassfield gives L = K(sqrt a1, sqrt a2, sqrt a3).  Any
\\ lift s of sigma to L has order 7 (the norm 1+sigma+..+sigma^6 kills V), so
\\ theta = sum_{i<7} s^i(sqrt a1) lies in F = L^<s>, and theta is not in Q (its components
\\ on the 7 nontrivial monomials sqrt(a_J) are nonzero), so Q(theta) = F (C7 is maximal).
\\ Its characteristic polynomial over K is computed in the multiquadratic algebra and
\\ has rational coefficients.

\\ ---------- cyclic septic fields of conductor dividing f0 ----------
oc_field(G, H, f) =
{
  if(f < 10^5, return(galoissubcyclo(G, H)));
  my(bnr = bnrinit(bnfinit(y, 1), [f, [1]]), M);
  M = Mat(vector(#G.gen, i, bnrisprincipal(bnr, G.gen[i], 0)));
  bnrclassfield(bnr, mathnf(matconcat([M * H, matdiagonal(bnr.cyc)])), 2);
}
oc_septics(S) =
{
  my(f0 = 1, G, L = List());
  foreach(S, q, if(q == 7, f0 *= 49, q % 7 == 1, f0 *= q));
  if(f0 == 1, return([]));
  G = znstar(f0, 1);
  foreach(subgrouplist(G.cyc, [7]), H, if(matdet(H) == 7,
    listput(L, polredabs(subst(oc_field(G, H, f0), x, y)))));
  Vec(L);
}

\\ ---------- certification (p = 2) ----------
oc_red(nf, x, P, modP) =
{
  if(type(x) != "t_MAT", x = Mat([x, 1]));
  for(i = 1, #x~, if(nfeltval(nf, x[i,1], P), return(0)));
  prod(i = 1, #x~, nfmodpr(nf, x[i,1], modP)^x[i,2]);
}
\\ proves the fundamental units, -1 and the alpha_j (d_j even) independent in K^*/K^*2,
\\ using quadratic residue characters at primes P of odd norm.  0 = gave up.
oc_psat2(bnf, maxprimes = 3000) =
{
  my(nf = bnf.nf, G = List(), U = bnfunits(bnf)[1], cyc = bnf.cyc, M = [;], q = 2, cnt = 0, n);
  for(i = 1, #U - 1, listput(G, U[i]));
  listput(G, bnf.tu[2]);                                   \\ w = 2 always here (K real)
  for(j = 1, #cyc, if(cyc[j] % 2 == 0,
    listput(G, bnfisprincipal(bnf, idealpow(nf, bnf.gen[j], cyc[j]), 4)[2])));
  n = #G;
  while(matrank(M) < n,
    if(cnt++ > maxprimes, return(0));
    q = nextprime(q + 1);
    if(nf.disc % q == 0, next);
    foreach(idealprimedec(nf, q), P,
      my(modP = nfmodprinit(nf, P), N = q^P[4], row = vector(n), ok = 1);
      for(i = 1, n, my(r = oc_red(nf, G[i], P, modP));
        if(r == 0, ok = 0; break);
        row[i] = (r^((N - 1) / 2) != 1));
      if(ok, M = matconcat([M; Mod(row, 2)]))));
  1;
}

\\ ---------- the ray class group of K and its sigma-action mod 2 ----------
\\ returns [bnr, sigma (as polmod), A2 = action on Cl_m/2 (columns), even indices]
oc_ray(Kpol, S, cert) =
{
  my(bnf = bnfinit(Kpol, 1), nf, m = 1, bnr, sig, A, ev);
  nf = bnf.nf;
  if(cert,
    if(!bnfcertify(bnf, 1), error("bnfcertify(,1) failed for ", Kpol));
    if(!oc_psat2(bnf), error("2-saturation not proved for ", Kpol)));
  foreach(S, q, foreach(idealprimedec(nf, q), P,
    m = idealmul(nf, m, idealpow(nf, P, if(q == 2, 2 * P[3] + 1, 1)))));
  bnr = bnrinit(bnf, [m, vector(7, i, 1)], 1);
  sig = nfgaloisconj(nf);
  sig = sig[select(s -> s != y, sig, 1)[1]];               \\ any generator of C7
  A = bnrgaloismatrix(bnr, sig);
  ev = select(d -> d % 2 == 0, bnr.cyc, 1);
  [bnr, Mod(sig, Kpol), matrix(#ev, #ev, i, j, Mod(A[ev[i], ev[j]], 2)), ev];
}

\\ the simple submodules W of ker f(A~) (dual of the S_f-isotypic part); for each, the
\\ subgroup H = W^perp of Cl_m (HNF), with Cl_m/H = F_2^3 irreducible
oc_subgroups(bnr, A2, ev, onlycount) =
{
  my(r = #ev, At = A2~, cyc = bnr.cyc, n = #cyc, res = List(), cnt = 0);
  if(r < 3, return(if(onlycount, 0, [])));
  foreach([x^3 + x + 1, x^3 + x^2 + 1], f,
    my(Kf = matker(subst(lift(f), x, At)), k = #Kf, seen = Map(), W);
    if(k == 0, next);
    if(k % 3, error("bug: dim ker f(A) = ", k));
    cnt += (8^(k/3) - 1) / 7;
    if(onlycount, next);
    forvec(c = vector(k, i, [0, 1]), if(c != vector(k),
      my(v = Kf * Mod(c, 2)~, B, key);
      if(mapisdefined(seen, lift(v)), next);
      B = matimage(matconcat([v, At * v, At^2 * v]));
      if(#B != 3, error("bug: orbit span"));
      forvec(e = vector(3, i, [0, 1]), if(e != [0,0,0], mapput(seen, lift(B * Mod(e, 2)~), 1)));
      \\ H = {x in Z^n : x_ev . w = 0 mod 2 for w in B}
      my(Hp = matker(B~), gens = List());
      for(j = 1, n, if(cyc[j] % 2, listput(gens, matid(n)[, j]), listput(gens, 2 * matid(n)[, j])));
      for(j = 1, #Hp, my(z = vectorv(n)); for(i = 1, r, z[ev[i]] = lift(Hp[i, j])); listput(gens, z));
      listput(res, mathnf(matconcat([Mat(Vec(gens)), matdiagonal(cyc)])))));
    if(#res != cnt && !onlycount, error("bug: subgroup count")));
  if(onlycount, cnt, Vec(res));
}

\\ ---------- the multiquadratic algebra K(sqrt a1, sqrt a2, sqrt a3) ----------
\\ an element is a vector of 8 polmods in y, index J+1, J a bitmask: sum u[J] sqrt(a_J)
oc_mul(u, v, a) =
{
  my(w = vector(8, i, 0 * a[1]));
  for(I = 0, 7, if(u[I+1] != 0, for(J = 0, 7, if(v[J+1] != 0,
    my(c = u[I+1] * v[J+1], B = bitand(I, J));
    for(j = 1, 3, if(bittest(B, j - 1), c *= a[j]));
    w[bitxor(I, J) + 1] += c))));
  w;
}
oc_sig(z, sig) = Mod(subst(lift(z), y, lift(sig)), sig.mod);
\\ image of u under the lift s with s(sqrt a_j) = img[j]
oc_lift(u, img, a, sig) =
{
  my(w = vector(8, i, 0 * a[1]), one = vector(8, i, 0 * a[1]));
  one[1] = 1 + 0 * a[1];
  for(J = 0, 7, if(u[J+1] != 0,
    my(t = one); t[1] = oc_sig(u[J+1], sig);
    for(j = 1, 3, if(bittest(J, j - 1), t = oc_mul(t, img[j], a)));
    w += t));
  w;
}
\\ sigma(a_j) = c^2 a_E: returns sqrt(sigma(a_j)) = c sqrt(a_E) as an algebra element
oc_sqrtimg(b, a, nf) =
{
  for(E = 0, 7,
    my(aE = 1 + 0 * a[1], r);
    for(j = 1, 3, if(bittest(E, j - 1), aE *= a[j]));
    r = nfroots(nf, x^2 - lift(b * aE));                  \\ sqrt(b aE) = c aE
    if(#r, my(u = vector(8, i, 0 * a[1])); u[E+1] = Mod(r[1], nf.pol) / aE; return(u)));
  error("bug: sigma(a_j) not in <a_1,a_2,a_3> K^*2");
}
\\ the octic field for H
oc_octic(bnr, H, sig) =
{
  my(nf = bnr.nf, Kp = nf.pol, P = bnrclassfield(bnr, H), a, img, e1, th, s, M, ch);
  if(#P != 3, error("bug: bnrclassfield gave ", #P, " polynomials"));
  \\ x^2 + b x + c: sqrt(b^2 - 4c) generates the same quadratic extension
  a = vector(3, j, my(p = P[j]); if(poldegree(p) != 2, error("bug: not quadratic"));
    Mod(lift(polcoef(p, 1))^2 - 4 * lift(polcoef(p, 0)), Kp));
  img = vector(3, j, oc_sqrtimg(oc_sig(a[j], sig), a, nf));
  e1 = vector(8, i, 0 * a[1]); e1[2] = 1 + 0 * a[1];
  th = e1; s = e1;
  for(i = 1, 6, s = oc_lift(s, img, a, sig); th += s);
  M = matrix(8, 8, i, j, my(b = vector(8, k, 0 * a[1])); b[j] = 1 + 0 * a[1]; oc_mul(th, b, a)[i]);
  ch = charpoly(M);
  ch = Pol(apply(c -> my(l = lift(c)); if(poldegree(l) > 0, error("bug: charpoly not over Q")); simplify(l), Vec(ch)));
  if(!polisirreducible(ch), error("bug: theta of degree < 8"));
  my(d = denominator(content(ch)));                          \\ monic in Z[x]: d^8 ch(x/d)
  polredabs(polredbest(d^8 * subst(ch, x, x / d)));
}

\\ ---------- verification of the output ----------
oc_verify(F, S) =
{
  if(polgalois(F)[1] != 56, error("not 8T25: ", F));
  my(D = nfdisc(F));
  foreach(factor(D)[, 1], q, if(!setsearch(S, q), error("ramified outside S: ", F)));
}

\\ ---------- main ----------
oct25count(S, cert = 1) =
{
  my(cnt = 0);
  S = Set(S);
  foreach(oc_septics(S), Kp, my(R = oc_ray(Kp, S, cert)); cnt += oc_subgroups(R[1], R[3], R[4], 1));
  cnt;
}
oct25(S, cert = 1) =
{
  my(t0 = getabstime(), ans = List(), cnt = 0, Ks);
  S = Set(S);
  Ks = oc_septics(S);
  foreach(Ks, Kp,
    my(R = oc_ray(Kp, S, cert), Hs = oc_subgroups(R[1], R[3], R[4], 0));
    cnt += #Hs;
    foreach(Hs, H, listput(ans, oc_octic(R[1], H, R[2]))));
  ans = Set(ans);
  foreach(ans, F, oc_verify(F, S));
  if(#ans != cnt, error("count mismatch: ", #ans, " fields, expected ", cnt));
  printf("oct25(%s): %d cyclic septic K, %d fields (%s, %d ms)\n", S, #Ks, #ans,
    if(cert, "GRH-free", "GRH"), getabstime() - t0);
  ans;
}
