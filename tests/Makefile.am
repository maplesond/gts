AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

TESTS = check_gts
check_PROGRAMS = check_gts
check_gts_CPPFLAGS = @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/src -I$(top_srcdir)/src/.deps
#check_gts_CXXFLAGS = -g3 -Wall -Werror -Wl,-rpath,@BAMTOOLS_PATH@/lib @AM_CXXFLAGS@
check_gts_CXXFLAGS = -g3 @AM_CXXFLAGS@
check_gts_LDFLAGS = @BOOST_LDFLAGS@ -L$(top_srcdir)/src/.deps @AM_LDFLAGS@
check_gts_LDADD = -lboost_unit_test_framework @BOOST_LIBS@
check_gts_SOURCES =	check_genbank.cpp \
			check_gff.cpp \
			check_gts.cpp

			    
			    


