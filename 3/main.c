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
	long long high[9], low[9], specialn /*, tmpl, llval,  newbnd */;
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

						makeoutfilename(s,j_p);
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

		high[2] = (long long) floor(((Dmaxroot/2.0) -shift[2])/arc[2]);
		low[2] = (long long) ceil(((- Dmaxroot/2.0) - shift[2])/arc[2]);

		/* A new invention, it is the effective multiplier for n */
		/* We remove some of the multipliers depending on the trace */

		/* Just for cubics */
		specialn = getspecialn(b1);


/*  ------------  Counter 2 ---------- */

	for(count[2] = low[2]; count[2] <= high[2]; count[2]++) {
		A2 = ((long long int) arc[2])*count[2] + shift[2];
		myaffsi( A2, coefs[2]);

		/* This is generic */
		shift[DEG] = 	(long long) truemod(((long long) arc[2])*count[2]
												 	*intpowmod(b1,DEG-2, specialn) 
													 -(DEG-1)* intpowmod(b1,DEG, specialn), specialn);

		high[DEG] = (long long) A2c_upper(boundn, DEG, specialn);
		low[DEG] = (long long) A2c_lower(-boundn, DEG, specialn);


		/*		shift[3] = ((long long) (DEG-2))*arc[2]*count[2]*b1 
			- b1*b1*b1*(DEG*(DEG-1)*(DEG-2))/3;

		centers[4] = (A2) * (A2)/2.0; */

		printf("\nLoop at (%1ld,%1ld) ---> (%1d,", count[1],count[2],
				(int)	 (ell/(2*arc[1])));
		plint(high[2]);
		printf("), %1d range: ", DEG);
		plint(low[DEG]);
		printf(" to ");
		plint(high[DEG]);
		printf("\n");
		fflush(stdout);
		/*		printf("\n       specialn = %1ld\n", specialn); */
		
/*  ------------  Counter 6 ---------- */

	for(count[DEG] =(long) low[DEG]; count[DEG]<= high[DEG]; count[DEG]++) {
		An = ((long long int) count[DEG])*specialn + shift[DEG];
		myaffsi( An, coefs[DEG]);

		/* printf("\nA6 = "); plint(A6); */
		if(( An != 0) && (mylabs(An) <= boundn)) {


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
										/*			fprintf(afile, "[x^3+");
										plintf(afile, A2);
										fprintf(afile, "*x+");
										plintf(afile, A3);
										fprintf(afile, ",%1ld], \\\n", b1);
										
										fflush(afile); */
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

