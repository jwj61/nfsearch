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
