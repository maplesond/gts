
#Good Transcript Selector (GTS)



##Installation:

  - If you cloned the git repository you must first run "./autogen.sh" to create the configure and make files for your project.  Do not worry if this fails due to missing dependencies at this stage.  If you downloaded a source code distribution tarball then you can skip this step.
  - GTS depends on samtools.  Instead of having to install this separately we provide the source code bundled and customised for GTS here.  To compile samtools change into the samtools directory and simply type ```make```
  - Now, for a typical installation on a machine where you have root access type ```./configure; make; sudo make install;```

Type ```./configure --help``` for full details.

The Makefile for GTS can take several goals.  Full details of common make goals can be found in the INSTALL file.  Typically, the following options can optionally used by KAT:

  - ```make check``` - runs unit tests.  Requires boost unit test framework to be installed and available.
  - ```make dist``` - packages the installation into a tarballed distributable.
  - ```make distcheck``` - runs some sanity tests to ensure the tarballed distributable is likely to work.


##Operating Instructions:

After GTS has been installed, the `gts` executable should be available.

Typing `gts` or `gts --help` at the command line will present you with the GTS help message.

##Notes

We have found Transdecoder V2.0.1 to produce transcripts with CDS and 5' and 3' UTRs, when it's possible the 5' or 3' UTRs are actually partials and not labelled as such in the GFF type output (they are labelled in the attribute tags).  This information is currently ignored by GTS.  What this means is that stage 2 (the 5' and 3' UTR check) will end up with fewer transcripts being filtered (as more transcripts will have 5' and 3' UTRs marked up), however these partial UTRs should still get filtered out in stage 3 when we look at comparing coordinates with full lengther.  The end result is that the output will be correct but the filtering will be heavier in stage 3 rather than stage 2.

##Dependencies

Uses Boost.  Developed with Boost V1.52.  Please make sure this is properly installed
and configured before compiling or running GTS.


##Licensing:

GNU GPL V3.  See COPYING file for more details.


##Authors:

Daniel Mapleson
David Swarbreck

See AUTHORS file for more details.


##Acknowledgements:

Affiliation: The Genome Analysis Centre (TGAC)
Funding: The Biotechnology and Biological Sciences Research Council (BBSRC)
