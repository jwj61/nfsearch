#include "c-target.h"
#include <stdio.h>

#define pp plist[pnum]
#define bb pdata[pnum][8][0]
#define aa pdata[pnum][8][1]
#define cc pdata[pnum][8][2]
#define dd pdata[pnum][8][3]
#define ee pdata[pnum][8][4]
#define i2 pdata[pnum][8][5]
#define i4 pdata[pnum][8][6]

/* Assume p>3
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* Let's START *************************************************************/
int t2221init(long pnum, long which, long what) {
  long long j, b2, tmp;
  long long d2, d3, d4, a2,a3,i8,i16,i32,i64,ad2,ad3,ab,a4,a5;
  

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][0][0] = 3;
    pdata[pnum][0][1] = 5;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1;
    pdata[pnum][2][1] = 0;
    pdata[pnum][DEG][0] = 1;
    pdata[pnum][DEG][1] = 0;
    pdata[pnum][3][0] = pp;
    pdata[pnum][4][0] = pp;
    pdata[pnum][5][0] = pp;
    pdata[pnum][6][0] = pp;

    i2 = (pp+1)/2;
    i4 = (i2*i2) % pp;

    for(j=1; j<=DEG; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pflags[pnum][DEG] = CALL_INIT; 
    pflags[pnum][3] = CALL_BOTH;
    break;
    /*

if a == 0, loop on d (rest are 0)
if a != 0, any constant term is ok, but this determines d
    */
  case 2:  /* c */
    dd = what; /* for future reference */
    break;
  case DEG:
    bb = what;
    break;
  case 3:
    aa = what; /* for future reference */
    if(aa==0) {
      ee=0;
      cc=0;
      pdata[pnum][3][1] = (bb*dd) % pp;
      b2 = (bb*bb) % pp;
      pdata[pnum][4][1] = truemod(i4*b2,pp);
      tmp = (b2*dd) % pp;
      pdata[pnum][5][1] = truemod(i4*tmp,pp);
      pdata[pnum][6][1] = 0;
    } else {
      ee=pp;
      cc=1;
      d2 = (dd*dd) % pp;
      d3 = (dd*d2) % pp;
      d4 = (d2*d2) % pp;
      b2 = (bb*bb) % pp;
      a2 = (aa*aa) % pp;
      a3 = (a2*aa) % pp;
      a4 = (a3*aa) % pp;
      a5 = (a4*aa) % pp;
      i8 = (i4*i2) % pp;
      i16 = (i8*i2) % pp;
      i32 = (i16*i2) % pp;
      i64 = (i32*i2) % pp;
      ad2 = (aa*d2) % pp;
      ad3 = (aa*d3) % pp;
      tmp = (i2*bb - 6*a2) % pp;
      ab = (aa * bb) % pp;
      /*      -1/8*d^3 - 1/8*a*d^2 + (1/2*b - 3/8*a^2)*d + 1/2*a*b + 5/8*a^3 + 2 */
      pdata[pnum][3][1] = truemod(-i8*d3 - i2*ad2 + tmp*dd + 2*ab
                                  + 40*a3 + 1, pp);
      tmp = (-i8*bb-5*i2*a2) % pp;
      /*1/64*d^4 - 1/16*a*d^3 + (-1/8*b - 5/32*a^2)*d^2 + (1/4*a*b + (7/16*a^3 + 1))
       *d + (1/4*b^2 - 1/8*a^2*b + (-15/64*a^4 + a)) */
      pdata[pnum][4][1] = truemod(i64*d4 - i4*ad3 + tmp*d2
                                  +((ab + 28*a3 + i2) % pp)*dd
                                  + i4*b2 - 2*((a2*bb)%pp) - 60*a4+ 2*aa, pp);
      /*      1/64*a*d^4 + 1/16*a^2*d^3 + (-1/8*a*b + (-1/32*a^3 - 1/4))*d^2 + (-1/4*a^2*b + (-3/16*a^4 + 1/2*a))*d + (1/4*a*b^2 + (3/8*a^3 + 1)*b + (9/64*a^5 - 1/4*a^2)) */
      pdata[pnum][5][1] = truemod(i16*((aa*d4)%pp)
                                  + ((a2*d3)%pp)
                                  + ((-i2*ab -2*a3 - i8) %pp)*d2
                                  + ((-4*a2*bb -48*a4 + aa) % pp)*dd
                                  + ((aa*b2)%pp)
                                  + (24*a3 + i2)*bb
                                     + 144*a5 - 2*a2, pp);
      /* -1/4*a*d^2 - 1/2*a^2*d + (a*b + (3/4*a^3 + 1)) */
      pdata[pnum][6][1] = truemod(-i2*ad2 - 4*((a2*dd)%pp)
                                  +2* ab + 24*a3 + i4, pp);
    }
    break;
  default:
    /*    pdata[pnum][which][0] = pp; */
    fprintf(stderr,"Error in t2221 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t2221next(long pnum, long which) {
  long long a2, a3, invc, invc2, invc4, invc6, invc8, a4,a5, invc10;
  long long i8, d2, d3,b2,e3,e2,e4,e5,d4,c2;

  if(which != 3) {
    fprintf(stderr,"Error in t2221 next.\n");
    return 0;
  }

  if(aa==0) { /* Separate loop? */
    if(ee<pp-1) {
      ee++;
      /*
       [d, b, -5*e^3 + 6*d*e^2 + (-b - 2*d^2)*e + d*b,
       -15/4*e^4 + 4*d*e^3  + (-1/2*b - d^2)*e^2 + 1/4*b^2,
       -9/2*e^5 + 33/4*d*e^4 + (-3*b -  5*d^2)*e^3 + (7/2*d*b + d^3)*e^2 + (-1/2*b^2 - d^2*b)*e + 1/4*d*b^2,
  0, 0]
      */
      d2 = (dd*dd) % pp;
      d3 = (dd*d2) % pp;
      b2 = (bb*bb) % pp;
      e2 = (ee*ee) % pp;
      e3 = (e2*ee) % pp;
      e4 = (e3*ee) % pp;
      e5 = (e4*ee) % pp;
      pdata[pnum][3][1] = truemod(-5*e3 + 6*dd*e2 + (-bb - 2*d2)*ee + dd*bb,
                                  pp);
      pdata[pnum][4][1] = truemod(-15*i4*e4 + 4*dd*e3  +
                                  ((-i2*bb - d2) % pp)*e2 + i4*b2,pp);
      pdata[pnum][5][1] = truemod(-9*i2*e5 + 33*i4*((dd*e4) % pp) +
                                  (-3*bb -  5*d2)*e3 +
                                  ((7*i2*((dd*bb)%pp) + d3) %pp)*e2 +
                                  ((-i2*b2 - d2*bb) %pp)*ee +
                                  i4*((dd*b2) % pp),pp);
      return 1;
    }
  }
  cc++;
  if(cc==pp)
    return 0;
  /* Now, don't know about a, but a and e have been set */

  a2 = (aa*aa) % pp;
  a3 = (a2*aa) % pp;
  a4 = (a3*aa) % pp;
  a5 = (a4*aa) % pp;
  d2 = (dd*dd) % pp;
  d3 = (d2*dd) % pp;
  d4 = (d3*dd) % pp;
  b2 = (bb*bb) % pp;

  c2 = (cc*cc) % pp;
  invc = minv(cc, pp);
  invc2 = (invc*invc) % pp;
  invc4 = (invc2*invc2) % pp;
  invc6 = (invc4*invc2) % pp;
  invc8 = (invc6*invc2) % pp;
  invc10 = (invc8*invc2) % pp;
  i8 = (i2*i4) % pp;
    
  pdata[pnum][3][1] = truemod(40*((invc6*a3)%pp) - 6*dd*((invc4*a2)%pp)
                              + ((2*invc2*bb - i2*((d2*invc2)%pp)) %pp)*aa +
                              i2*((dd*bb)%pp) + cc- i8*d3, pp);
  pdata[pnum][4][1] = truemod(-60*invc8*a4 + 28*dd*((invc6*a3)%pp) +
                              ((-2*invc4*bb - 5*i2*((d2*invc4)%pp))%pp)*a2 +
                              ((dd*invc2)%pp)*((bb*aa)%pp)+
                              (((8*cc - d3)*i4)%pp)*((invc2*aa)%pp)
                              + i4*b2 - i8*d2*bb + i2*((dd*cc)%pp) +
                              i8*((i8*d4)%pp), pp);
  pdata[pnum][5][1] = truemod(144*invc10*a5 - 48*dd*((invc8*a4)%pp) +
                              ((24*invc6*bb - 2*d2*invc6)%pp)*a3 +
                              ((-4*dd*((invc4*bb)%pp) + (-2*cc + d3)*invc4)%pp)*a2 +
                              ((invc2*b2 - ((i2*d2)%pp)*((invc2*bb)%pp) +
  (((dd*cc + i8*((i2*d4)%pp))%pp)*invc2))%pp)*aa +
                              i2*((cc*bb)%pp) - i8*((d2*cc)%pp), pp);
  pdata[pnum][6][1] = truemod(24*((invc4*invc)%pp)*a3
                              - 4*((dd*invc2)%pp)*((invc*a2)%pp)
                              + ((2*invc*bb - i2*((d2*invc)%pp))%pp)*aa
                              + i4*c2, pp);
  
  return 1;
}


