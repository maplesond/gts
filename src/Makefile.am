AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# Create gts binary
bin_PROGRAMS = gts \
	       gffids \
	       gff_filter \
	       gtf_filter \
	       gbfilter \
	       fix_gtf

gts_CXXFLAGS = -g3 @AM_CXXFLAGS@
gts_CPPFLAGS = -DCPLUSPLUS @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/src 
gts_LDFLAGS = @BOOST_LDFLAGS@ @AM_LDFLAGS@
gts_LDADD = @BOOST_LIBS@
gts_SOURCES =	filters/transcript_filter.hpp \
		filters/inconsistent_coords_filter.hpp \
		filters/utr_filter.hpp \
		filters/multiple_orf_filter.hpp \
		filters/strand_filter.hpp \
		filters/overlap_filter.hpp \
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

gff_filter_CXXFLAGS = -g3 @AM_CXXFLAGS@
gff_filter_CPPFLAGS = -DCPLUSPLUS @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/src 
gff_filter_LDFLAGS = @BOOST_LDFLAGS@ @AM_LDFLAGS@
gff_filter_LDADD = @BOOST_LIBS@
gff_filter_SOURCES = gff.hpp \
		gff_filter.cc
	
gtf_filter_CXXFLAGS = -g3 @AM_CXXFLAGS@
gtf_filter_CPPFLAGS = -DCPLUSPLUS @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/src 
gtf_filter_LDFLAGS = @BOOST_LDFLAGS@ @AM_LDFLAGS@
gtf_filter_LDADD = @BOOST_LIBS@
gtf_filter_SOURCES = gff.hpp \
		gtf_filter.cc
                
gbfilter_CXXFLAGS = -g3 @AM_CXXFLAGS@
gbfilter_CPPFLAGS = -DCPLUSPLUS @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/src 
gbfilter_LDFLAGS = @BOOST_LDFLAGS@ @AM_LDFLAGS@
gbfilter_LDADD = @BOOST_LIBS@
gbfilter_SOURCES = genbank.hpp \
		gb_filter.cc
	
fix_gtf_CXXFLAGS = -g3 @AM_CXXFLAGS@
fix_gtf_CPPFLAGS = -DCPLUSPLUS @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/src 
fix_gtf_LDFLAGS = @BOOST_LDFLAGS@ @AM_LDFLAGS@
fix_gtf_LDADD = @BOOST_LIBS@
fix_gtf_SOURCES = gff.hpp \
		fix_gtf.cc
	