AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# Create gts binary
bin_PROGRAMS = gts

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
		gts.cc
                
