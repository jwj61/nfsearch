/* Find quadratic lifts unramified outside a set
*/

default(primelimit, 10^8);
default(realprecision, 250);
/*bnf=bnfinit(x^5 - x^4 - 4*x^3 + 3*x^2 + 3*x - 1); */

ramps(f)=return(nfdiscfactors(f)[2][,1]~);

/* Shorthand */
q(pol) = return(newlars(ramps(pol),pol));

/* GRH-free certification of what degree-ell class field theory uses: Cl(K)/ell and E/ell.
   bnfcertify(bnf,1) proves that the factor base generates Cl(K); ql_psat proves that the
   fundamental units, zeta (if ell | w) and alpha_j with (alpha_j) = g_j^{d_j} (ell | d_j) are
   independent mod ell-th powers (ell-th power residue characters at degree-1 primes).  See
   grhfree_count.tex, Lemmas 5-8. */
ql_psatred(nf, x, P, modP) =
{
  if(type(x) != "t_MAT", x = Mat([x, 1]));
  for(i = 1, #x~, if(nfeltval(nf, x[i,1], P), return(0)));
  prod(i = 1, #x~, nfmodpr(nf, x[i,1], modP)^x[i,2]);
}
ql_psat(bnf, p, maxprimes = 2000) =
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
      for(i = 1, n, my(r = ql_psatred(nf, G[i], P, modP));
        if(r == 0, ok = 0; break);
        row[i] = fflog(r^((q - 1) / p), zp, p));
      if(ok, M = matconcat([M; Mod(row, p)]))));
  1;
}

/* Cyclic degree-ell extensions L of K = Q[x]/(pol) unramified outside S.
   S: rational primes and/or prime ideals of K (from nfinit(subst(pol,x,y))).
   arch: real places allowed to ramify ([-1]: all of them when ell = 2, none otherwise).
   flag = 0: just list fields (no polredabs)
   flag = 1: the same via rnfpolredabs(,,2);
   flag = 2: relative polynomials (rnfpolredabs), one per extension L/K, so L and its
             conjugates under Aut(K) are listed separately.
   wram > 0: only extensions whose conductor has exponent exactly wram at every prime of K
             above ell that is in S (wild part fixed; nothing else is restricted).
   Raises an error unless Cl(K)/ell and E/ell are certified without GRH. */
newlars(S, pol, ell=2, arch=[-1], flag=0, wram=0) =
{
  my(maxB, bnf, bnr, sgl, maxe, np, ans, wP = List());

  bnf = bnfinit(subst(pol,x,y),1);
  S = subst(S,x,y);
  if(bnfcertify(bnf, 1) != 1, error("bnfcertify(bnf,1) failed"));
  if(!ql_psat(bnf, ell), error("ell-saturation of the units failed"));
  maxB = 1;
  for(k=1,#S,
    my(ipd = if(type(S[k])=="t_INT", idealprimedec(bnf, S[k]), [S[k]]));
    for(j=1,#ipd,
      if(ipd[j].p == ell,
        maxe = idealpow(bnf, ipd[j], if(wram, wram, floor(ell*ipd[j].e/(ell-1))+1));
        listput(wP, ipd[j]);
        maxB = idealmul(bnf, maxB, maxe),
      /* else: tame, exponent 1 */
        maxB = idealmul(bnf, maxB, ipd[j]))));
  if(ell==2 && arch[1]==-1,
    arch = vector(polsturm(pol), h, 1);
    maxB = [maxB, arch],
    if(arch[1] != -1,
      maxB=[maxB, arch])
  );
  bnr = bnrinit(bnf,maxB,1);
  sgl = subgrouplist(bnr,[ell],1);       /* all subgroups of index ell */
  ans = [];
  for(k=1,#sgl,
    if(wram,
      my(c = bnrconductor(bnr, sgl[k])[1]);
      if([P | P <- wP, idealval(bnf, c, P) != wram] != [], next));
    np = if(flag == 2, rnfpolredabs(bnf, bnrclassfield(bnr, sgl[k], 1)),
         flag == 1, rnfpolredabs(bnf, bnrclassfield(bnr, sgl[k], 1), 2),
         bnrclassfield(bnr, sgl[k], 2));
    ans = concat(ans, np);
  );
  return(Set(ans));
}


/* List of conductors unramified outside of S */

/* This part is not working yet */
condlistS(bnf, S, ell=2) =
{
  local(plist,pbnd,ipd,plistn,ent,pl2,minn,spot);

/*  n= poldegree(bnf.pol); */
  /* list of [ideal, norm] */
  plist=List();
  listput(plist,[1,1],1);
  spot=1;
  for(k=1, #S,
    p=S[k];
    ipd=idealprimedec(bnf,p);
    for(j=1, #ipd,
      if(p==ell,
        for(kk=2, floor(p*(ipd[j].e)/(p-1)+1),
          ent= [p^((ipd[j].f)*kk), idealpow(bnf, ipd[j], kk)];
          spot++;
          listput(plist, ent,spot)
        );
      , /* else */
        ent=[p^(ipd[j].f), ipd[j]];
      );
    );
  );
  forprime(p=3, pbnd,
    print1(p, " -> ", pbnd , "(", length(plist),") ", );
    ipd=idealprimedec(bnf, p);
    print1(" ipd ", #ipd, " - ");
    for(j=1, #ipd,
      ent=[p^(ipd[j].f), ipd[j]];
      plistn=spot; k=1;
      print1(" ", ent[1]);
      while(k<=plistn && ent[1]*plist[k][1]<=pbnd,
        spot++;listput(plist, [plist[k][1]*ent[1], idealmul(bnf, plist[k][2],ent[2])], spot);
        k++;
      );
      listsort(plist);
    );
      print();
  );
  return(plist);
}

lifts(pol, B)=
{
  local(np,bnf,res, condl,arch,sgl,bnr,rnf,spot=0);
  bnf=bnfinit(subst(pol,x,y),1);
  if(bnfcertify(bnf) == 0, return(1/0));
  res=List();
  condl=condlist(bnf, B);
  print("Got cond list of length ", length(condl));
  arch=vector(bnf.r1,h,0);
  for(j=1,length(condl),
    bnr=bnrinit(bnf, [condl[j][2], arch], 1);
    sgl=subgrouplist(bnr, [2]);
    for(k=1,length(sgl),
      rnf=rnfkummer(bnr,sgl[k]);
      np=rnfpolredabs(bnf,rnf,2);
      spot++; listput(res,np, spot);
    );
  );
  return(res);
}

