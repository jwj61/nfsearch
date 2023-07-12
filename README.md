# nfsearch
Targeted number field searching

Programs were created for my own personal use.  The source code is full of artifacts of various features which may no longer be supported.

There are search programs for different degree fields, each in its own directory.

Installation
============

In the directory you want, edit the Makefile to set paths to pari library and header files.  Also set the path for the Inputs directory (making it the full path to the Inputs directory in the same directory you are in).

To compile in degree 3, then run "make c3" (for degree 4, "make c4").

Running the program
===================

To search for degree 6 fields which are ramified at p and q (and only p and q), type
  c6 p q
You will be prompted for the ramification pattern at each prime.  Type the sequence of letters for the options you want.  (If the factorization of p in the desired field has invariants (e1,f1), ... (ek, fk), then the ramification pattern is the list of ei, each repeated fi times).
