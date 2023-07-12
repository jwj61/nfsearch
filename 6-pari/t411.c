#include "c-target.h"

#define pp plist[pnum]
#define bb pdata[pnum][7][0]
#define cc pdata[pnum][7][1]
#define aa pdata[pnum][7][2]

/* Assume p>3
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* Let's START *************************************************************/
int t411init(long pnum, long which, long what) {
  long long j, k;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][0][0] = 3;
    pdata[pnum][0][1] = 3;
    
    pflags[pnum][1] = 1;
    pflags[pnum][2] = -1;
    pflags[pnum][6] = 0;
    pflags[pnum][3] = -1;
    for(j=4; j<=5; j++)
      pflags[pnum][j] = 1;
    
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1;
    pdata[pnum][2][1] = 0;
    break;
  case 2:  /* c */
    bb = what; /* for future reference */
    break;
  case 6: /* 10*a^6 - 4*b*a^5 + c*a^4 = a^4*(10*a^2 - 4*b*a + c)*/
    pdata[pnum][6][0] = pp;
    pdata[pnum][6][1] = 0;
    aa = 0;
    cc = what;
    break;
  case 3: /* 20*a^3 - 10*b*a^2 + 4*c*a */
    pdata[pnum][3][0] = pp;
    j = truemod(20*aa*aa-10*bb*aa+4*cc, pp);
    pdata[pnum][3][1] = (aa*j) %  pp;
     /* 45*a^4 - 20*b*a^3 + 6*c*a^2 = a^2 (45*a^2 - 20*b*a + 6*c)  */
    pdata[pnum][4][0] = pp;
    j = (aa*aa) % pp;
    k = truemod(45*j-20*bb*aa+6*cc, pp);
    pdata[pnum][4][1] = (j*k) % pp;
     /* 36*a^5 - 15*b*a^4 + 4*c*a^3 = a^2(36*a^3 - 15*b*a^2 + 4*c*a)*/
    pdata[pnum][5][0] = pp;
    j = (aa*aa) % pp;
    k = truemod(36*aa*j - 15*bb*j + 4*cc*aa, pp);
    pdata[pnum][5][1] = (j*k) % pp;
    break;
  default:
    pdata[pnum][which][0] = pp;
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t411next(long pnum, long which) {
  long long  j, k, l;

  if(which==6) {
    ++(aa);
    if(aa==pp)
      return 0; /* Done with this case! */
    /* 20*a^3 - 10*b*a^2 + 4*c*a */
    j = (aa*aa) % pp;
    k = (j*j) % pp;
    l = truemod(10*j - 4*bb*aa + cc, pp);
    pdata[pnum][which][1] = (l*k) % pp;
    return 1;
  }
  /* in all other cases, there is only one case to consider */
  return 0;
}


/* 

411
    
[b, c,
  20*a^3 - 10*b*a^2 + 4*c*a,
  45*a^4 - 20*b*a^3 + 6*c*a^2,
  36*a^5 - 15*b*a^4 + 4*c*a^3,
  10*a^6 - 4*b*a^5 + c*a^4]


--------------------------------------

33

 [a, b, -5/27*a^3 + 2/3*b*a, -1/9*b*a^2 + 1/3*b^2,
   1/81*a^5 - 2/27*b*a^3 + 1/9*b^2*a,
   -1/729*a^6 + 1/81*b*a^4 - 1/27*b^2*a^2 + 1/27*b^3]
 
 --------------------------------------
321
 
[c, -6*a^2 + (-6*b + 3*c)*a + (-3*b^2 + 2*c*b),
  -8*a^3 + (-18*b + 3*c)*a^2 + (-12*b^2 + 6*c*b)*a + (-2*b^3 + c*b^2),
  -3*a^4 + (-18*b + c)*a^3 + (-18*b^2 + 6*c*b)*a^2 + (-6*b^3 + 3*c*b^2)*a,
  -6*b*a^4 + (-12*b^2 + 2*c*b)*a^3 + (-6*b^3 + 3*c*b^2)*a^2,
  -3*b^2*a^4 + (-2*b^3 + c*b^2)*a^3]

  --------------------------------------
3111

 [b, c, (10*a^6 - 6*b*a^5 + 3*c*a^4 + d)/a^3,
   (15*a^6 - 8*b*a^5 + 3*c*a^4 + 3*d)/a^2,
   (6*a^6 - 3*b*a^5 + c*a^4 + 3*d)/a, d]

 and

 [b, c, d, 0, 0, 0]
 
  --------------------------------------
222

 [a, b, -1/8*a^3 + 1/2*b*a + 2*c, 1/64*a^4 - 1/8*b*a^2 + c*a + 1/4*b^2,
   -1/4*c*a^2 + c*b, c^2]

  --------------------------------------
2211

  [c, d, 4*a^3 - 3*c*a^2 + (-6*b + 2*d)*a + 2*c*b,
    3*a^4 - 2*c*a^3 + d*a^2 - 2*c*b*a + (-3*b^2 + 2*d*b),
    6*b*a^3 - 4*c*b*a^2 + (-6*b^2 + 2*d*b)*a + c*b^2,
    3*b^2*a^2 - 2*c*b^2*a + (-2*b^3 + d*b^2)]

  --------------------------------------
21111

 From a,b,c,d,e

 [b, e, c, d, 0, 0]

 or
 
 [b, e, c,
   (-5*a^6 + 4*b*a^5 - 3*e*a^4 + 2*c*a^3 + d)/a^2,
   (-4*a^6 + 3*b*a^5 - 2*e*a^4 + c*a^3 + 2*d)/a,
   d]
*/
