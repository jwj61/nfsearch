# nfsearch
Targeted number field searching

Programs were created for my own personal use.  The source code is full of artifacts of various features which may no longer be supported.

There are search programs for different degree fields, each in its own directory.

The find all primitive extensions of Q of the chosen degree with ramification patterns at primes as specified by the user.

Installation
============

In the directory you want, edit the Makefile to set paths to pari library and header files.  Also set the path for the Inputs directory (making it the full path to the Inputs directory in the same directory you are in).

To compile in degree 3, then run "make c3" (for degree 4, "make c4").

Running the program
===================

To search for degree 6 fields which are ramified at p and q (and only p and q), type
  c6 p q
You will be prompted for the ramification pattern at each prime.  Type the sequence of letters for the options you want.  (If the factorization of p in the desired field has invariants (e1,f1), ... (ek, fk), then the ramification pattern is the list of ei, each repeated fi times).

Note, if you want all primitive fields which are unramified outside a set, you may need to run the program multiple times since it does not allow you to pick "unramified" at any of the primes used when starting the program.

What is the output?
===================

  * When the program starts, it creates a file "Current" which describes the current search.  It is mainly of use when a script is used to do many runs of the program.
  * As it runs, it adds polynomials, one per line, to a file "Search..." where the name of the file contains the information about which search was run.
  * When the program finishes, it adds a line to a log file called "Cubic.log", "Quartic.log", etc., depending on the degree of the search.  It contains some information on the run time of the program, and a little diagnostic information on how many polynomials were examined, and how many survived various tests.
  * The polynomials in the output file Search... have been put through polredabs, so typically a polynomial will appear multiple times.  It does not check that the polynomials hit the target exactly, so the user may need to filter those out after the fact.

After running the program, the results can be easily read into a gp session with "readvec", and then take "Set" of the list to remove repeats.

Speed-ups
=========

  * Before the exact discriminant is computed, each polynomial goes through a cheap modular test (common-pari/c-qfilter.c).  For a target field, disc(f) = d_K * (index)^2 and the exponent of each p in S in d_K is fixed by the current case, so disc(f) has a known square class.  For 24 primes q = 1 mod 4 outside S, the program checks that disc(f) is in that class modulo q (whenever q does not divide disc(f)).  The answer depends only on the innermost coefficient modulo q, so it is cached per residue for each innermost loop.  Almost all polynomials are rejected this way, and the exact discriminant is computed only for the rest.  The run summary on the screen reports how many polynomials the filter tested and rejected.
  * Because the filter uses the discriminant exponents of the current case, fields with a different exponent at some p in S are no longer found incidentally.  They are found in their own case when that case is searched.  The filter is turned off automatically in slop mode, and it can be removed entirely by compiling with -DNO_QFILTER.
  * In even degree, when A1 = 0 the search takes A3 >= 0 (alpha -> -alpha fixes A1 and negates A3).
  * In degrees 6 and 7, a_{n-3} is also bounded through the reciprocal roots: |a_{n-3}| <= |a_n| C(n,3) (S_{-1}/n)^3, where S_{-1} is Pohst's bound on the sum of |1/alpha|.
