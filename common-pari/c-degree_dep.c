#include <stdio.h>
#include "c-target.h"

#define min(x,y) (x<y) ? x : y

#if DEG==6
#define name "Sextic"
#else
#if DEG==5
#define name 	"Quintic"
#else
#if DEG==4
#define name "Quartic"
#else
#if DEG==3
#define name 	"Cubic"
#else
#if DEG==7
#if defined(IS_TR)
#define name 	"Septic-tr"
#else
#define name 	"Septic"
#endif
#else
#if DEG==8
#define name 	"Octic"
#else
#define name   "???"
#endif
#endif
#endif
#endif
#endif
#endif

void welcome()
{
  printf("\nJohn's search program for ");
  printf("%ss!\n\n", name);
#if DEG>8
  exit(1);
#endif
}

long
  pohst_thing()
{
  register int i, j;
  double    kpows[8], degroots[DEG+2], dff, newval;
  long      sgn_n, degover2;
  long long absAn;

  if(An>0) sgn_n=1; else sgn_n= -1;
  degover2 = DEG/2;
  absAn = mylabs(An);
  for(i=1; i<DEG; i++)
    pbounds[i] = -9999;
  for(i=1;i<=degover2; i++) kpows[i] = fpow(absAn, 2.0/i);

  for(i=1;i<=degover2; i++)  {    /* Use Newton's method to find some roots */
    degroots[i-1] =  0.9*fpow(Dmaxroot1/DEG, 0.5);

    while((fabs(dff = newt(degroots[i-1], i, kpows[i])) > 0.001) ||
          (fabs(dff) > 0.001*degroots[i-1]) || (degroots[i-1] < dff)) {
      if(degroots[i-1] < dff) degroots[i-1] = degroots[i-1]/2.0;
      else degroots[i-1] = degroots[i-1] - dff;
      
    }
  }

  for(i=1;i<=degover2; i++)  {    /* Use Newton's method to find more roots */
    degroots[i+degover2-1] = 2.0* fpow(Dmaxroot1/DEG, 0.5) - degroots[i-1];
    while(fabs(dff = newt(degroots[i+degover2-1], i, kpows[i])) > 0.001) {
      degroots[i+degover2-1] = degroots[i+degover2-1] - dff;
    }
  } /* Roots now found by Newton's method */
  
  for(i=1; i<3; i++) {
    for(j=0; j<2*degover2; j++) {
      newval = thef(-1*i, ((j%degover2)+1),degroots[j],absAn);
      if(newval > pbounds[i]) {
        pbounds[i] = newval;
      }
    }
  }
  for(i=3; i<DEG; i++) {
    for(j=0; j<2*degover2; j++) {
      newval = thef(i, ((j%degover2)+1),degroots[j],absAn);
      if(newval > pbounds[i]) {
        pbounds[i] = newval;
      }
    }
  }

  for(i=3; i<DEG; i++)
    pbounds[i] = min(bounds[i], pbounds[i]);
  return(sgn_n);
}

void
  makeoutfilename(char *s, long *j_p)
{
  /* int i;
     char t[10]; */

/*   sprintf(s, "Search_%1ld%c", plist[0], (int) j_p[0]+'a'); */
/*   for(i=1; i<totprimes; i++) { */
/*     sprintf(t, ".%1ld%c", plist[i], (int) j_p[i]+'a'); */
/*     strcat(s,t); */
/*   } */

/*   sprintf(s, "Search_%1ld-%1d", plist[0], (int) j_p[0]); */
/*   for(i=1; i<totprimes; i++) { */
/*     sprintf(t, ".%1ld-%1d", plist[i], (int) j_p[i]); */
/*     strcat(s,t); */
/*   } */
  if(ALTERNATING) {
	sprintf(s, "eSearch_%s.%1ld", outfilename, slop);
  } else {
    sprintf(s, "Search_%s.%1ld", outfilename, slop);
  }
}

/* Uses degree dependent macro "name" */
void
getlogfilename(char *t)
{
  if(ALTERNATING) {
    sprintf(t, "e%s.log", name);
  } else {
    sprintf(t, "%s.log", name);
  }
}

void
  spitsplittingmenu(void)
{
#if DEG==3
  printf("\na: 3  b: 21  *: user\n");
#else
#if DEG==4
  printf("\na: 4  b: 31  c: 22  d: 211  *: user\n");
#else
#if DEG==5
  printf("\na: 5  b: 41  c: 32  d: 311  e: 221  f: 2111  *: user\n");
#else
#if DEG==6
  printf(
         "\na: 6  b: 51  c: 42  d: 411  e: 33  f: 321  g: 3111  h: 222  i: 2211  j: 2  *: user\n");
#else
#if DEG==7
  printf(
         "\na: 7  b: 61  c: 52  d: 511  e: 43  f: 421  g: 4111  h: 331  i: 322\n j: 3211  k: 31111  l: 2221  m: 22111  n: 211111  *: user\n");
#else
#if DEG==8
  printf(
         "\na: 8  b: 71  c: 62  d: 611  e: 53  f: 521  g: 5111  h: 44  i: 431\n j: 422  k: 4211  l: 41111  m: 332  n: 3311  *: user\n");
#endif
#endif
#endif
#endif
#endif
#endif
  printf("\nYour selections: ");
  fflush(stdout);
}

void
  splittingmenustring(char *t, long j)
{
  switch (j) {
#if DEG==3
  case 0: strcpy(t, "3"); break;
  case 1: strcpy(t, "2"); break;
#else
#if DEG==4
  case 0: strcpy(t, "4"); break;
  case 1: strcpy(t, "3"); break;
  case 2: strcpy(t, "22"); break;
  case 3: strcpy(t, "2"); break;
#else
#if DEG==5
  case 0: strcpy(t, "5"); break;
  case 1: strcpy(t, "4"); break;
  case 2: strcpy(t, "32"); break;
  case 3: strcpy(t, "3"); break;
  case 4: strcpy(t, "22"); break;
  case 5: strcpy(t, "2"); break;
#else
#if DEG==6
  case 0: strcpy(t, "6"); break;
  case 1: strcpy(t, "5"); break;
  case 2: strcpy(t, "42"); break;
  case 3: strcpy(t, "4"); break;
  case 4: strcpy(t, "33"); break;
  case 5: strcpy(t, "32"); break;
  case 6: strcpy(t, "3"); break;
  case 7: strcpy(t, "222"); break;
  case 8: strcpy(t, "22"); break;
  case 9: strcpy(t, "2"); break;
#else
#if DEG==7
  case 0: strcpy(t, "7"); break;
  case 1: strcpy(t, "6"); break;
  case 2: strcpy(t, "52"); break;
  case 3: strcpy(t, "5"); break;
  case 4: strcpy(t, "43"); break;
  case 5: strcpy(t, "42"); break;
  case 6: strcpy(t, "4"); break;
  case 7: strcpy(t, "33"); break;
  case 8: strcpy(t, "322"); break;
  case 9: strcpy(t, "32"); break;
  case 10: strcpy(t, "3"); break;
  case 11: strcpy(t, "222"); break;
  case 12: strcpy(t, "22"); break;
  case 13: strcpy(t, "2"); break;
#else
#if DEG==8
  case 0: strcpy(t, "8"); break;
  case 1: strcpy(t, "7"); break;
  case 2: strcpy(t, "62"); break;
  case 3: strcpy(t, "6"); break;
  case 4: strcpy(t, "53"); break;
  case 5: strcpy(t, "52"); break;
  case 6: strcpy(t, "5"); break;
  case 7: strcpy(t, "44"); break;
  case 8: strcpy(t, "43"); break;
  case 9: strcpy(t, "422"); break;
  case 10: strcpy(t, "42"); break;
  case 11: strcpy(t, "4"); break;
  case 12: strcpy(t, "332"); break;
  case 13: strcpy(t, "33"); break;
#endif
#endif
#endif
#endif
#endif
#endif
  default: strcpy(t, "error");
  }
}

/* j contains the case and p is the prime */
int iswild(long j, long p){
  switch(j) {
#if DEG==3
  case 0: if(p==3) return 1;
    return 0;
  case 1: if(p==2) return 1;
    return 0;
#else
#if DEG==4
  case 0: if(p==2) return 1;
    return 0;
  case 1: if(p==3) return 1;
    return 0;
  case 2: if(p==2) return 1;
    return 0;
  case 3: if(p==2) return 1;
    return 0;
#else
#if DEG==5
  case 0: if(p==5) return 1;
    return 0;
  case 1: if(p==2) return 1;
    return 0;
  case 2: if((p==2) || (p==3)) return 1;
    return 0;
  case 3: if(p==3) return 1;
    return 0;
  case 4: if(p==2) return 1;
    return 0;
  case 5: if(p==2) return 1;
    return 0;
#else
#if DEG==6
  case 0: if((p==2) || (p==3)) return 1;
    return 0;
  case 1: if(p==5) return 1;
    return 0;
  case 2: if(p==2) return 1;
    return 0;
  case 3: if(p==2) return 1;
    return 0;
  case 4: if(p==3) return 1;
    return 0;
  case 5: if((p==2) || (p==3)) return 1;
    return 0;
  case 6: if(p==3) return 1;
    return 0;
  case 7: if(p==2) return 1;
    return 0;
  case 8: if(p==2) return 1;
    return 0;
  case 9: if(p==2) return 1;
    return 0;
#else
#if DEG==7
  case 0: if(p==7) return 1;
    return 0;
  case 1: if((p==2) || (p==3)) return 1;
    return 0;
  case 2: if((p==2) || (p==5)) return 1;
    return 0;
  case 3: if(p==5) return 1;
    return 0;
  case 4: if((p==2) || (p==3)) return 1;
    return 0;
  case 5: if(p==2) return 1;
    return 0;
  case 6: if(p==2) return 1;
    return 0;
  case 7: if(p==3) return 1;
    return 0;
  case 8: if((p==2) || (p==3)) return 1;
    return 0;
  case 9: if((p==2) || (p==3)) return 1;
    return 0;
  case 10: if(p==3) return 1;
    return 0;
  case 11: if(p==2) return 1;
    return 0;
  case 12: if(p==2) return 1;
    return 0;
  case 13: if(p==2) return 1;
    return 0;
#endif
#endif
#endif
#endif
#endif
  default:
    return 0;
  }
}
