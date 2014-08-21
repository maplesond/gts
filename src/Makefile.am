AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# Create gts binary
bin_PROGRAMS = gts \
	       gffids \
	       gfffilter \
	       gbfilter

gts_CXXFLAGS = -g3 @AM_CXXFLAGS@
gts_CPPFLAGS = -DCPLUSPLUS @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/src 
gts_LDFLAGS = @BOOST_LDFLAGS@ @AM_LDFLAGS@
gts_LDADD = @BOOST_LIBS@
gts_SOURCES =	filters/transcript_filter.hpp \
		filters/inconsistent_coords_filter.hpp \
		filters/multiple_orf_filter.hpp \
		filters/multiple_transcript_filter.hpp \
		filters/strand_filter.hpp \
		filters/transcript_filter.hpp \
		gff.hpp \
		fln.hpp \
		genbank.hpp \
		gts.cc
	
gffids_CXXFLAGS = -g3 @AM_CXXFLAGS@
gffids_CPPFLAGS = -DCPLUSPLUS @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/src 
gffids_LDFLAGS = @BOOST_LDFLAGS@ @AM_LDFLAGS@
gffids_LDADD = @BOOST_LIBS@
gffids_SOURCES =    gff.hpp \
		    gff_ids.cc

gfffilter_CXXFLAGS = -g3 @AM_CXXFLAGS@
gfffilter_CPPFLAGS = -DCPLUSPLUS @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/src 
gfffilter_LDFLAGS = @BOOST_LDFLAGS@ @AM_LDFLAGS@
gfffilter_LDADD = @BOOST_LIBS@
gfffilter_SOURCES = gff.hpp \
		gff_filter.cc
                
gbfilter_CXXFLAGS = -g3 @AM_CXXFLAGS@
gbfilter_CPPFLAGS = -DCPLUSPLUS @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/src 
gbfilter_LDFLAGS = @BOOST_LDFLAGS@ @AM_LDFLAGS@
gbfilter_LDADD = @BOOST_LIBS@
gbfilter_SOURCES = genbank.hpp \
		gb_filter.cc
	