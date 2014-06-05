//  ********************************************************************
//  This file is part of GTS (Good Transcript Selector).
//
//  GTS is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Portculis is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Portculis.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <iostream>
#include <fstream>
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
using boost::to_upper_copy;
using boost::timer::auto_cpu_timer;
using boost::shared_ptr;
using boost::unordered_map;
namespace po = boost::program_options;
namespace bfs = boost::filesystem;

#include "gff.hpp"

namespace gts {

typedef boost::error_info<struct GTSError,string> GTSErrorInfo;
struct GTSException: virtual boost::exception, virtual std::exception { };

typedef std::vector< boost::shared_ptr<GFF> > GFFList;

class GFFData {
public:
    uint32_t mrna;
    uint32_t utr5;
    uint32_t utr3;
    shared_ptr<GFF> gff;
    
    GFFData():mrna(0),utr5(0),utr3(0) { }
    
    void incMrna() {
        mrna++;
    }
    
    void incUtr5() {
        utr5++;
    }
    
    void incUtr3() {
        utr3++;
    }
};


typedef boost::unordered_map<string, GFFData> IdMap;
typedef boost::unordered_map<string, uint32_t> IdCounter;


class GTS {
    
private:
    
    GFFList genomicGffs;
    IdMap transcripts;
    IdMap uniqGffs;
    IdMap uniqUtrGffs;
    
    double cds;
    bool include;
    bool verbose;
    
protected:
    
    void filterMultipleOrfs() {
    
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        
        cout << "Stage 1: Filter out transcripts with multiple ORFs" << endl;
        IdCounter utr5;
        IdCounter utr3;
        vector<string> orfs;
        
        BOOST_FOREACH(shared_ptr<GFF> gff, genomicGffs) {

            string id = gff->GetId();
            switch(gff->GetType()) {
            case MRNA:
                transcripts[id].incMrna();
                orfs.push_back(id);
                break;
            case UTR5:
                utr5[id]++;
                break;
            case UTR3:
                utr3[id]++;
                break;
            }
        }

        vector<string> unique;
        vector<string> uniqueUtr;
        BOOST_FOREACH(IdMap::value_type i, transcripts) {
            
            if (i.second.mrna == 1) { 
                uniqGffs[i.first] = i.second;
                if (utr3[i.first] > 0 && utr5[i.first] > 0) {
                    uniqUtrGffs[i.first] = i.second;
                }
            }
        }

        cout << " - # transcripts: " << transcripts.size() << endl
             << " - # ORFs: " << orfs.size() << endl
             << " - # transcripts with one ORF: " << uniqGffs.size() << endl
             << " - # transcripts with one ORF and 5' and 3' UTRs: " << uniqUtrGffs.size() << endl;

    }
    
public :
    
    GTS(string& genomicGffFile, double cds, bool include, bool verbose) :
        cds(cds), include(include), verbose(verbose)
    {
        // Check for all required inputs
        
        if (!bfs::exists(genomicGffFile) && ! bfs::symbolic_link_exists(genomicGffFile)) {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                    "Could not find specific genomic GFF file: ") + genomicGffFile));
        }
        
        /*if (!bfs::exists(transcriptGffFile) && ! bfs::symbolic_link_exists(transcriptGffFile)) {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                    "Could not find specific transcript GFF file: ") + transcriptGffFile));
        }
        
        if (!bfs::exists(flnResultsDir) && ! bfs::symbolic_link_exists(flnResultsDir)) {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                    "Could not find full lengther results directory: ") + flnResultsDir));
        }*/
            
        GFF::loadGFF(genomicGffFile, genomicGffs);        
        cout << "Found " << genomicGffs.size() << " GFF records in " << genomicGffFile << endl << endl;
    }

    virtual ~GTS() {}
    
    void filter() {
        filterMultipleOrfs();
    }
    


    
};

}




string helpHeader() {
    return string("\nGTS Help.\n\n") +
                  "GTS (Good transcript selector) is a tool to filter out all transcripts that we are not very confident are genuine\n\n" +
                  "Usage: gts [options] --ggff <genomic gff file> --tgff <transcript gff file> --fln_dir <full lengther results dir>\n\n" +
                 "\nAvailable options";
}


/**
 * Start point for portculis.
 */
int main(int argc, char *argv[]) {
    
    try {
        // Portculis args
        string genomicGffFile;
        string transcriptGffFile;
        string flnResultsDir;
        string outputPrefix;
        bool include;
        double cds;
                
        string cufflinksGtfFile;
        
        bool verbose;
        bool version;
        bool help;

        // Declare the supported options.
        po::options_description generic_options(helpHeader());
        generic_options.add_options()
                ("ggff,g", po::value<string>(&genomicGffFile), 
                    "Gff file containing the genomic coordinates for the transcript features.")
                ("tgff,t", po::value<string>(&transcriptGffFile),
                    "Gff file containing the transcript coordinates for the transcript features.")
                ("fln_dir,f", po::value<string>(&flnResultsDir),
                    "Full lengther results directory, containing the \"dbannotated.txt\" and \"new_coding.txt\" files.")
                ("output,o", po::value<string>(&outputPrefix)->default_value(string("gts_out")),
                    "The output prefix for all output files generated.")
                ("include,i", po::bool_switch(&include)->default_value(false), 
                    "Include transcripts with no full lengther homology hit.")
                ("cds,c", po::value<double>(&cds)->default_value(0.5), 
                    "Min percentage length of CDS relative to mRNA for hits with homology.  0.0 -> 1.0")
                ("version", po::bool_switch(&version)->default_value(false), "Print version string")
                ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
                ;

        // Combine non-positional options
        po::options_description cmdline_options;
        cmdline_options.add(generic_options);

        // Parse command line
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
        po::notify(vm);

        // Output help information the exit if requested
        if (argc == 1 || argc == 2 && help) {
            cout << generic_options << endl;
            return 1;
        }

        // Output version information then exit if requested
        if (version) {
#ifndef PACKAGE_NAME
#define PACKAGE_NAME "GTS"
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1.0"
#endif
            cout << PACKAGE_NAME << " V" << PACKAGE_VERSION << endl;
            return 0;
        }
        
        
        
        
        const string gffPassOut = outputPrefix + ".pass.gff3";
        const string gffFailOut = outputPrefix + ".pass.gff3";
        
        auto_cpu_timer timer(1, "Total wall time taken: %ws\n\n");
                
        gts::GTS gts(genomicGffFile, cds, include, verbose);
        gts.filter();
                
    } catch (boost::exception &e) { 
        std::cerr << boost::diagnostic_information(e); 
        return 4;
    } catch (exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 5;
    } catch (const char* msg) {
        cerr << "Error: " << msg << endl;
        return 6;
    } catch (...) {
        cerr << "Error: Exception of unknown type!" << endl;
        return 7;
    }

    return 0;
}


