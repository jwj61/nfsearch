#include <stdio.h>
#include <time.h>
#include "genpari.h"
#include "target.h"
#include "doout.h"

void main(int argc, char * argv [])
{
	struct        timeval time1, time2;
	struct        timezone tzn;
	FILE          *afile;
	unsigned long timetot;
	long          newmodu[10], divbyme[10], b1ore[10];
	long long     invertme[DEG+1], divmod[DEG+1], byme;
	long          sgn_n, j_p[10];
  register long i, count[10], ii;
  register long long b1;


	char      s[512];
	long long high[9], low[9], tmpl, llval, newbnd, specialn;
  long      ltop, vp;
	double    btmp, boundn;
	/*	double    newval, lval, savelval; */

	/* End of variable declarations */

  init(500000,500000);  /* take some memory for the stack */

	welcome();
	myinitprimes();
	getprimelist(argc, argv);
	adjustplistfordegree();	
	initmygens();
	getsearchcases();

	/* Loop through cases */

	for(j_p[0]=0; j_p[0]<=maxp[0]; ++j_p[0]) {
		for(j_p[1]=0; j_p[1]<=maxp[1]; ++j_p[1]) {
			for(j_p[2]=0; j_p[2]<=maxp[2]; ++j_p[2]) {
				for(j_p[3]=0; j_p[3]<=maxp[3]; ++j_p[3]) {
					for(j_p[4]=0; j_p[4]<=maxp[4]; ++j_p[4]) {
					for(j_p[5]=0; j_p[5]<=maxp[5]; ++j_p[5]) {
					for(j_p[6]=0; j_p[6]<=maxp[6]; ++j_p[6]) {
					for(j_p[7]=0; j_p[7]<=maxp[7]; ++j_p[7]) {
					for(j_p[8]=0; j_p[8]<=maxp[8]; ++j_p[8]) {

						printf("\nThis is case number %ld", j_p[0]+1);
						for(i=1; i<totprimes; i++)
							printf(" x %ld", j_p[i]+1);
						printf(" of %1ld", maxp[0]+1);
						for(i=1; i<totprimes; i++)
							printf(" x %ld", maxp[i]+1);
						printf("\n");

						initcase(j_p, divmod, invertme);
						showorefactors();

						makeoutfilename(s, j_p);
						afile = fopen(s, "w");


						gettimeofday(&time1, &tzn);
						affsi((long) 0, coefs[1]);
						boundn = pow(Dmaxroot/DEG, DEG/2.0);
						ltop = avma;

	/* arc[1] has been divided by 6 */

	for(count[1] = 0;  2*count[1]* arc[1] <= ell; count[1]++) {
		b1 = count[1] * arc[1];
		if(b1 != 0)
			for(i=0; i<totprimesplus; i++) 
				b1ore[i] = myval(b1, plist[i]);

		/* This is generic */
		shift[2] = ((long long) -DEG)*(DEG-1)*b1*b1/2;

		high[2] = (long long) S2c_upper(Dmaxroot, 2);
		low[2] = (long long) S2c_lower(- Dmaxroot, 2);

		/* A new invention, it is the effective multiplier for n */
		/* We remove some of the multipliers depending on the trace */

		specialn = getspecialn(b1);
		/*		printf("\nb1= %1ld, specialn = ",b1);
		plint(specialn);  */


/*  ------------  Counter 2 ---------- */

	for(count[2] = low[2]; count[2] <= high[2]; count[2]++) {
		A2 = arc[2]*count[2] + shift[2];
		myaffsi( A2, coefs[2]);

		/* This is generic */
		shift[DEG] = 	truemod(arc[2]*count[2]
												 	*intpowmod(b1,DEG-2, specialn) 
													 -(DEG-1)* intpowmod(b1,DEG, specialn), specialn);

		high[DEG] = (long long) A2c_upper(boundn, DEG, specialn);
		low[DEG] = (long long) A2c_lower(-boundn, DEG, specialn);


		shift[3] = ((long long) (DEG-2))*arc[2]*count[2]*b1 
			- b1*b1*b1*(DEG*(DEG-1)*(DEG-2))/3;

		centers[4] = (A2) * (A2)/2.0;

		printf("\nLoop at (%1ld,%1ld) ---> (%1d,", count[1],count[2],
				(int)	 (ell/(2*arc[1])));
		plint(high[2]);
		printf("), %1d range: ", DEG);
		plint(low[DEG]);
		printf(" to ");
		plint(high[DEG]);
		printf("\n");
		/*		printf("\n       specialn = %1ld\n", specialn); */
		fflush(stdout);
		
/*  ------------  Counter 4 ---------- */

	for(count[DEG] =(long) low[DEG]; count[DEG]<= high[DEG]; count[DEG]++) {
		An = ((long long int) count[DEG])*specialn + shift[DEG];
		myaffsi( An, coefs[DEG]);

		/* printf("\nA4 = "); plint(An); */
		if(( An != 0) && (mylabs(An) <= boundn)) {

			/* This is not needed in degree 3 */
			/* It mostly does side effects, and captures the sign of
				 An (whoppee!!). */
			sgn_n = pohst_thing();
			

			/* Let's talk about the shift.  We are trying to insure a certain
				 congruence for An modulo arc[DEG].  We only have forced the
				 congruence modulo specialn.  As we pick each coefficient after
				 that, we may inherit a congruence to keep An alive.  So we need
				 to compute that. */

			if(specialn == arc[DEG]) { /* No extra stuff happening */
				invertme[3] = 0;
				divmod[3] = 1;
			} else {
				for(i=0; i<totprimesplus; i++) {
					divbyme[i] = min(ore[i][4], b1ore[i]+ore[i][3]);
				}
				byme = vectoint(divbyme);
				
				if(((tmpl = (An + 3* intpow(b1,4)- arc[2]*count[2]*intpow(b1,2)))
						% byme) != 0) {
					printf("\nbyme="); plint(byme);
					printf(" shift[4]="); plint(shift[4]);
					printf(" arc[4]="); plint(arc[4]);
					printf(" tmpl="); plint(tmpl);
					printf(" computed:");
					plint(An +  3* intpow(b1,4)- arc[2]*count[2]*intpow(b1,2));
					printf(" An= ");
					plint(An); printf("\n");
					printf("\nI think there is a problem with b1, a2, An: ");
					plint(b1); printf(", ");
					plint(count[2]*arc[2]); printf(", ");
					plint(An); printf("\n"); exit(1);
				}
				tmpl = tmpl/byme;
				invertme[3] = (arc[3]*b1)/byme;
				divmod[3] = arc[4]/byme;
				/* Now we invert invertme */
				invertme[3] = truemod(((long long) tmpl)
															* gcdtrick(invertme[3] % divmod[3],
																				 divmod[3]), divmod[3]);
			}

			low[3] =  (long long) S2c_lower(-pbounds[3],3);
			high[3] = (long long) S2c_upper(pbounds[3],3);

		/*    printf("\nDivmod = %1d", divmod); */
		
		/* Save some values which get used by inner loops, and here */

/*  ------------  Counter 3 ---------- */


	for(count[3] = low[3]; count[3]<= high[3]; count[3]++) {
		A3 = arc[3] *(invertme[3] + count[3]*divmod[3]) + shift[3];
		myaffsi( A3, coefs[3]);


/*  ------------  The actual testing   ---------- */


				polcount++;
				gop1z(absi, discsr(a), td); /* take absolute value on the way */
				if(!gcmp0(td)) {
					for(i=0; i<totprimes; i++) {
						vp = pvaluation(td, gprimes[i], &tdtmp);
						gaffect(tdtmp, td);
					}
					if(gcmp1(td)) {
						if(gcmp1(gisirreducible(a))) {
							DUMPIT
						}
					} else {
						if(carrecomplet(td, &sqrtdust)) {
							postdisc++;
#ifdef TTWO
							if(gcmpsg((long) ceil(Dmaxroot), gnorml2(roots(a,3))) > 0) {
								t2count++;
#endif
								if(gcmp1(gisirreducible(a))) {
									gop1z(absi, discf(a), td);
									discfcount++;
									for(i=0; i<totprimes; i++) {
										vp = pvaluation(td, gprimes[i], &tdtmp);
										gaffect(tdtmp, td);
									}
									if(gcmp1(td)) {
										DUMPIT
									} /* Discf check */
								} /* Irreducible */
#ifdef TTWO
							} /* T2 check */
#endif
						}
					}
				}

				avma = ltop;
	} /* Counter 3 */
	/* } A6 satisfying a basic congruence */
		} /* Is A6 != 0 */
	} /* Counter 6 */
	} /* Counter 2 */
	} /* Counter 1 */
	printf("\nDone.  Searched ");
	plint(polcount);
	printf(" polynomials for %1ld^%1ld", plist[0], maxexpo[0]);
	for(i=1; i<totprimes; i++)
		printf(" %1ld^%1ld", plist[i], maxexpo[i]);
	printf(".\n");
	printf("Survived discrim check %1ld times.\n", postdisc);
	printf("Field discrim. checked %1ld times.\n\n", discfcount);

	gettimeofday(&time2, &tzn);
	timetot = time2.tv_sec - time1.tv_sec;
	fprintf(logfile, "(%1ld:%c",plist[0], (int) p_inps[0][j_p[0]]-1+'a');
	for(i=1; i<totprimes; i++)
		fprintf(logfile, ",%1ld:%c", plist[i], (int) p_inps[i][j_p[i]]-1+'a');
	fprintf(logfile,"), %1ld polynomials, %1ld seconds, %3.2f T2, %1ld past disc, %1ld past T2, %1ld discf\n",
				  polcount, timetot, Dmaxroot, postdisc, t2count, discfcount);
	fflush(logfile);
	fclose(afile);
					}
					}
					}
					}
					}
				}
			}
		}
	}
}

