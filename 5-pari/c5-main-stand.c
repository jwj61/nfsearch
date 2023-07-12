#include <stdio.h>
#ifdef ALPHA
#include <sys/time.h>
#else
/* #include <time.h> */
#include <sys/time.h>
#include <unistd.h>
#endif
#include "c-target.h"
#include "c-doout.h"


int main(int argc, char * argv [])
{
  struct        timeval time1, time2;
  struct        timezone tzn;
  FILE          *afile;
  unsigned long timetot;
  long          sgn_n;
  register long i, count[10];


  char      s[512];
  long long high[9], low[9], llval, newbnd, bigB;
  long      ltop;
  double    boundn, btmp, savelval, lval, newval;
  GEN       bigbgen;

  /* End of variable declarations */

  pari_init(2000000,500000);  /* take some memory for the stack */

  welcome();
  myinitprimes();
  /* getprimelist(argc, argv); */
  totprimes = 0;

  initmygens();
  /* getsearchcases(); */

/*  bigB = 1000;
  bigB = bigB*10; */
  bigB = atoi(argv[1]);
  bigbgen = stoi((long) bigB);
  verbosity = 3;

  /* makeoutfilename(s, j_p); */
  sprintf(s, "StdSearch_%Ld", bigB);

  /* Make a little file which says what case we are currently doing */
  afile = fopen("Current", "w");
  fprintf(afile, "Working on %s.\n", outfilename);
  fclose(afile);
  
  afile = fopen(s, "a");

  polcount = 0;
  postdisc = 0;
  t2count = 0;
  discfcount = 0;
  gettimeofday(&time1, &tzn);
  /* Loop through cases */

                  
  /* initcase(j_p); This will do the init 0 things */

  Dmaxroot = 1.0 * bigB/(1.0* DEG );
  Dmaxroot *= mygamma[DEG-1];
  Dmaxroot = t2tester fpow(Dmaxroot, 1.0/(DEG-1));
  Dmaxroot += 5./4.;
  printf("Max T2 for %Ld is %5.3f\n", bigB, Dmaxroot);

  for(i=0; i<7; i++) {
	arc[i] = 1;
	shift[i] = 0;
  }
  ltop = avma;

  /* ------------ Counter 1 ---------- */
  for(count[1] = 0;  2*count[1]  <= DEG; count[1]++) {
	A1 = count[1];
	myaffsi(A1, coefs[1]);
	
	/* This is |S_i| */
	for(i=3;i<=DEG; i++) {
	  bounds[i] =  fpow(Dmaxroot + A1*A1/DEG*1.0, i/2.0);
	  /* printf("\nBounds[%1d] = %6.2f ", i, bounds[i]); */
	}
                    
	boundn = fpow(Dmaxroot1/DEG, DEG/2.0);
	centers[2] = A1*A1*(DEG-1)/(2.0*DEG);
	/* This is generic */

	high[2] = (long long) S2c_upper(Dmaxroot, 2);
	low[2] = (long long) S2c_lower(- Dmaxroot, 2);
	
	/*  ------------  Counter 2 ---------- */
	
	for(count[2] = low[2]; count[2] <= high[2]; count[2]++) {
	  A2 = count[2];
	  myaffsi( A2, coefs[2]);
	  
	  /* This is generic */
      
	  high[DEG] = (long long) A2c_upper(boundn, DEG, arc[DEG]);
	  /* New adjustment for deg odd! */
	  if(A1==0)
		low[DEG] = 0;
	  else
		low[DEG] = (long long)
		  A2c_lower(-boundn, DEG, arc[DEG]);
	  
	  
	  centers[3] = A1*A2 - A1*A1*A1/3.0;
	  
	  if(verbosity>2) {
		printf("\nLoop at (%1ld,%1ld) ---> (%1d,",
			   count[1], count[2],
			   (int) ((DEG)/(2)));
		plint(high[2]);
		printf("), %1d range: ", DEG);
		plint(low[DEG]);
		printf(" to ");
		plint(high[DEG]);
		printf("\n");
		fflush(stdout);
	  }
                      
		
	  /*  ------------  Counter 5 ---------- */
	  
	  for(count[DEG] =(long) low[DEG];
		  count[DEG]<= high[DEG]; count[DEG]++) {
		An = count[DEG];
		myaffsi(An, coefs[DEG]);

/*
if(A1 == 0 && A2 == 70 && An == 1364) {
printf("testing ");plint(A1);printf(",");plint(A2);printf(",");plint(An);
printf("\n");
}
*/
		/*		printf("\nA5 = "); plint(An); */
		if(( An != 0) && (mylabs(An) <= boundn)) {
		  
                          /* It mostly does side effects, and captures
                             the sign of An (whoppee!!). */
		  
		  sgn_n = pohst_thing();
          
		  low[3] =  (long long) S2c_lower(-pbounds[3],3);
		  high[3] = (long long) S2c_upper(pbounds[3],3);
                          
		  
		  /* Inherited from An bound */
		  savelval =  -An*pbounds[2]/2.0;
		  if(An>0) {
			newbnd = (long long) ceil((savelval - shift[3])
									  /(1.0*arc[3]));
			if(low[3]<newbnd) low[3] = newbnd;
		  } else {
			newbnd = (long long) floor((savelval - shift[3])
									   /(1.0*arc[3]));
			if(high[3]>newbnd) high[3] = newbnd;
		  }
          
          
		  /*  ------------  Counter 3 ---------- */
		  
		  for(count[3] = low[3]; count[3]<= high[3];
			  count[3]++) {
			A3 = count[3];
			myaffsi(A3, coefs[3]);

			centers[4] = A1*A1*A1*A1/4.0 -A1*A1*A2 +
			  A2*A2/2.0 +A1*A3;
			
			if(((2.0*A3)/An+pbounds[2]) < 0) {
			  printf("Funny square root: %3.2f\n",
					 ((2.0*A3)/An+pbounds[2]));
			  printf("%3.2f is the bound, An = ", savelval);
			  plint(An);
			  printf(" and A3 = ");
			  plint(A3);
			  printf("\n");
			  printf("Is %3.2f plus %3.2f positive?\n",
					 (2.0 * A3)*(An),An*An*pbounds[2]);
			} else {
			  btmp = (Dmaxroot1- S2)*(Dmaxroot1-S2) -
				2.0*S3*S3/(S2+Dmaxroot1);
			  
			  if(btmp > pbounds[4])  btmp = pbounds[4];
			  low[4] =  (long long) S2c_lower(-pbounds[4],4);
			  high[4] = (long long) S2c_upper(btmp,4);
			  
			  /* |A4| <= lval */
			  lval = sqrt(((2.0*A3)/An+pbounds[2]))*An*sgn_n;
			  /* |A5| <= newval */
			  newval = pbounds[1]*mylabs(An);
			  if(lval > newval) lval = newval;
			  
			  llval = (long long) A2c_upper(lval,4,arc[4]);
			  if(llval< high[4]) high[4] = llval;
			  
			  llval = (long long) A2c_lower(-lval,4,arc[4]);
			  if(llval> low[4]) low[4] = llval;
			  
			  
			  
			  /*  ------------  Counter 4 ---------- */
			  
			  for(count[4] = low[4]; count[4]<= high[4];
				  count[4]++) {
				A4 = count[4];
				myaffsi(A4, coefs[4]);
                
				/*  ------------  The actual testing   ---------- */
                
                
				polcount++;
				/* take absolute value on the way */
				/* if(issquare(discsr(a))) {  */
                
			  postdisc++;
			  /* outbrute(td);printf(",\n"); 
				 printf("%l\n", gissquarerem(td, &sqrtdust)); */
				if(isirreducible(a)) {
				  td=absi( discf(a));
				  if(gcmp(td, bigbgen)<=0) {
				    discfcount++;
					pari_fprintf(afile, "%Ps\n", polredabs0(a,0));
				/*    DUMPIT */
				  } /* small disc*/
				} /* Irreducible */
/*				}  */  /* disc is a square */
                                    
                
				avma = ltop;
			  } /* Counter 4 */
			}
		  } /* Counter 3 */
		} /* Is A5 != 0 */
	  } /* Counter 5 */
	} /* Counter 2 */
  } /* Counter 1 */
  if(verbosity>2) {
	if(polcount!= 0) {
	  printf("\nSearched ");
	  plint(polcount);
	  printf(" polynomials for %1Ld", bigB);
	  printf(".  ");
	  printf("Poly disc check %1ld times, ", postdisc);
	  printf("Field disc check %1ld times.\n",
			 discfcount);
	} else printf("0");
    
  }
  gettimeofday(&time2, &tzn);
  timetot = time2.tv_sec - time1.tv_sec;
/*   fprintf(logfile, "(%Ld,",bigB); */
/*   fprintf(logfile,"), %1ld polynomials, %1ld seconds, \ */
/* %3.2f T2, %1ld past disc, %1ld past T2, %1ld discf\n", */
/*           polcount, timetot, Dmaxroot, postdisc, t2count, discfcount); */
/*   fflush(logfile); */
  fclose(afile);
  printf("\nDone.");
  printf("\n");
  return 0;
}


