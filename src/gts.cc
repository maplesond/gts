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
#include "fln.hpp"

namespace gts {

typedef boost::error_info<struct GTSError,string> GTSErrorInfo;
struct GTSException: virtual boost::exception, virtual std::exception { };

typedef std::vector< boost::shared_ptr<GFF> > GFFList;
typedef std::vector< boost::shared_ptr<DBAnnot> > FLNDBAnnotList;

typedef boost::unordered_map<string, shared_ptr<GFF> > GFFIdMap;
typedef boost::unordered_map<string, shared_ptr<DBAnnot> > DBAnnotIdMap;
typedef boost::unordered_map<string, uint32_t> IdCounter;

const int32_t POS_THRESHOLD = 10;
const int32_t LONG_CDS_LEN_THRESHOLD = 200;

class GTS {
    
private:
    
    GFFList genomicGffs;
    GFFList transdecoderCdsGffs;
    FLNDBAnnotList flnDbannots;
    FLNDBAnnotList flnNc;
    GFFIdMap genomicGffMap;
    GFFIdMap uniqGffs;
    GFFIdMap uniqUtr;
    GFFIdMap uniqCds;
    GFFIdMap consistentCds;
    DBAnnotIdMap uniqFlnCds;
    DBAnnotIdMap uniqFlnNcCds;
    
    vector<string> uniqUtrFlnCds;
    vector<string> uniqUtrFln;
    vector<string> uniqUtrNc;
    vector<string> stage2Pass;
    
    
    double cds;
    bool include;
    bool verbose;
    
protected:
    
    void updateTDVFLN(const string& id, shared_ptr<GFF> tdc, shared_ptr<DBAnnot> fln, bool longCds) {
        int32_t deltaStart = std::abs(tdc->GetStart() - fln->GetStart());
        int32_t deltaEnd = std::abs(tdc->GetEnd() - fln->GetEnd());
        int32_t tdcLen = tdc->GetEnd() - tdc->GetStart();

        if (deltaStart <= POS_THRESHOLD && deltaEnd <= POS_THRESHOLD && 
                (!longCds || (longCds && tdcLen >= LONG_CDS_LEN_THRESHOLD)) ) {

            uniqUtrFlnCds.push_back(id);

            const double seqFrac = (double)tdcLen / (double)fln->GetFastaLength();

            if (seqFrac >= cds) {
                consistentCds[id] = tdc;
            }
        }               
    }
    
    void filterInconsistent() {
        
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        
        cout << "Stage 2: Filter out CDS that are inconsistent between transdecoder and full lengther" << endl;
        
        BOOST_FOREACH(shared_ptr<GFF> gff, transdecoderCdsGffs) {
        
            if (gff->GetType() == CDS) {
            
                const string cdsid = gff->GetCdsid();
                if (uniqUtr.count(cdsid)) {
                    uniqCds[cdsid] = gff;
                }
            }
        }
        
        // Index fln cdss by id
        BOOST_FOREACH(shared_ptr<DBAnnot> db, flnDbannots) {
            uniqFlnCds[db->GetId()] = db;
        }
        
        // Index fln cdss by id
        BOOST_FOREACH(shared_ptr<DBAnnot> db, flnNc) {
            uniqFlnNcCds[db->GetId()] = db;
            uniqFlnCds[db->GetId()] = db;
        }
        
        BOOST_FOREACH(GFFIdMap::value_type i, uniqCds) {
            
            if (uniqFlnCds.count(i.first)) {
                
                shared_ptr<GFF> transdecoder = i.second;
                
                if (uniqFlnCds.count(i.first)) {
                    
                    shared_ptr<DBAnnot> fln = uniqFlnCds[i.first];
                    
                    if (fln->GetStatus() == COMPLETE) {                    
                        updateTDVFLN(i.first, transdecoder, fln, false);
                    }                
                }
                else if (include && uniqFlnNcCds.count(i.first)) {
                    updateTDVFLN(i.first, transdecoder, uniqFlnNcCds[i.first], false);
                }
            }
        }
        
        stage2Pass.reserve(uniqUtrFln.size() + uniqUtrNc.size());
        stage2Pass.insert(stage2Pass.end(), uniqUtrFln.begin(), uniqUtrFln.end());
        stage2Pass.insert(stage2Pass.end(), uniqUtrNc.begin(), uniqUtrNc.end());
        
        cout << " - # Transdecoder CDSs with IDs matching those from stage 1: " << uniqCds.size() << endl
             << " - # Full lengther CDSs found in dbannotated.txt: " << flnDbannots.size() << endl
             << " - # Full lengther CDSs found in new_coding.txt: " << flnNc.size() << endl
             << " - # Total distinct full lengther CDSs: " << uniqFlnCds.size() << endl            
             << " - # Transcripts with consistent fln and transdecoder CDS coordinates: " << uniqUtrFlnCds.size() << endl
             << " - # Transcripts with similarity passing CDS coordinate check: " << uniqUtrFln.size() << endl
             << " - # Transcripts with no similarity passing CDS coordinate check (will be 0 if --include wasn't used): " << uniqUtrNc.size() << endl
             << " - # Transcripts passing filter 2 (CDS coordinate check): " << stage2Pass.size() << endl;
                
    }
    
    void filterMultipleOrfs() {
    
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        
        cout << "Stage 1: Filter out transcripts with multiple ORFs and no 5' and 3' UTRs" << endl;
        IdCounter mrna;
        IdCounter utr5;
        IdCounter utr3;
        vector<string> orfs;
        
        BOOST_FOREACH(shared_ptr<GFF> gff, genomicGffs) {

            string id = gff->GetId();
            
            genomicGffMap[id] = gff;
            
            switch(gff->GetType()) {
            case MRNA:
                mrna[id]++;
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
        BOOST_FOREACH(IdCounter::value_type i, mrna) {
            
            if (i.second == 1) { 
                uniqGffs[i.first] = genomicGffMap[i.first];
                if (utr3[i.first] > 0 && utr5[i.first] > 0) {
                    uniqUtr[i.first] = genomicGffMap[i.first];
                }
            }
        }

        cout << " - # ORFs: " << orfs.size() << endl
             << " - # transcripts: " << mrna.size() << endl            
             << " - # transcripts with one ORF: " << uniqGffs.size() << endl
             << " - # transcripts with one ORF and both 5' and 3' UTRs: " << uniqUtr.size() << endl;

    }
    
public :
    
    GTS(const string& genomicGffFile, const string& transcriptGffFile, const string& flnDir, double cds, bool include, bool verbose) :
        cds(cds), include(include), verbose(verbose)
    {
        // Check for all required inputs
        
        if (!bfs::exists(genomicGffFile) && ! bfs::symbolic_link_exists(genomicGffFile)) {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                    "Could not find specific genomic GFF file: ") + genomicGffFile));
        }
        
        if (!bfs::exists(transcriptGffFile) && ! bfs::symbolic_link_exists(transcriptGffFile)) {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                    "Could not find specific transcript GFF file: ") + transcriptGffFile));
        }
        
        if (!bfs::exists(flnDir) && ! bfs::symbolic_link_exists(flnDir)) {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                    "Could not find full lengther results directory: ") + flnDir));
        }
        
        string dbAnnotFile = flnDir + "/dbannotated.txt";
        if (!bfs::exists(dbAnnotFile)) {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                    "Could not find full lengther dbannotated.txt file at: ") + dbAnnotFile));
        }
        
        string ncFile = flnDir + "/new_coding.txt";
        if (!bfs::exists(ncFile)) {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                    "Could not find full lengther new_coding.txt file at: ") + ncFile));
        }
        
        auto_cpu_timer timer(1, "Total load time: %ws\n\n");
        
        cout << "Loading inputs" << endl
             << "--------------" << endl << endl;
        
        cout << "Loading Genomic GFF file" << endl;
        GFF::load(genomicGffFile, genomicGffs);

        cout << "Loading Transcript GFF file" << endl;
        GFF::load(transcriptGffFile, transdecoderCdsGffs);

        cout << "Loading Full Lengther DB Annot file" << endl;
        DBAnnot::load(dbAnnotFile, flnDbannots);
        
        cout << "Loading Full Lengther New Coding file" << endl;
        DBAnnot::load(ncFile, flnNc);
        
    }

    virtual ~GTS() {}
    
    void filter() {
        
        auto_cpu_timer timer(1, "Total filter time: %ws\n\n");
        
        cout << "Filtering transcripts" << endl
             << "---------------------" << endl << endl;
        
        filterMultipleOrfs();
        filterInconsistent();
    }
    
    void saveGB(const string& gbFile) {
        
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
                    "Include transcripts with no full lengther homology hit, providing it has a full lengther new_coding hit.")
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
                
        // Select good transcripts
        gts::GTS gts(genomicGffFile, transcriptGffFile, flnResultsDir, cds, include, verbose);
        gts.filter();
        
        // Create GB from GFF objects (output only single transcript per locus)
        gts.saveGB("out");
                
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


