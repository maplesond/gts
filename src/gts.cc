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
//  along with GTS.  If not, see <http://www.gnu.org/licenses/>.
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
#include "filters/transcript_filter.hpp"
#include "filters/multiple_orf_filter.hpp"
#include "filters/inconsistent_coords_filter.hpp"
#include "filters/multiple_transcript_filter.hpp"
#include "filters/strand_filter.hpp"

namespace gts {

typedef boost::error_info<struct GTSError,string> GTSErrorInfo;
struct GTSException: virtual boost::exception, virtual std::exception { };

typedef std::vector< boost::shared_ptr<DBAnnot> > FLNDBAnnotList;


class GTS {
    
private:
    
    
    shared_ptr<GFFList> genomicGffs;
    GFFList transdecoderCdsGffs;
    GFFList gtfs;
    FLNDBAnnotList flnDbannots;
    FLNDBAnnotList flnNc;
    
    Maps maps;
    
    const string genomicGffFile;
    const string transcriptGffFile;
    const string flnDir;
    
    string outputPrefix;
    string gtfsFile;
    double cds;
    bool include;
    bool outputAllStages;
    bool verbose;
    
protected:
    
    
    
    void load() {
        
        auto_cpu_timer timer(1, "Total load time: %ws\n\n");
        
        cout << "Loading inputs" << endl
             << "--------------" << endl << endl;
        
        if (!bfs::exists(genomicGffFile) && ! bfs::symbolic_link_exists(genomicGffFile)) {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                    "Could not find specific genomic GFF file: ") + genomicGffFile));
        }
        
        if (!bfs::exists(transcriptGffFile) && ! bfs::symbolic_link_exists(transcriptGffFile)) {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                    "Could not find specific transcript GFF file: ") + transcriptGffFile));
        }
        
        if (!gtfsFile.empty()) {
            if (!bfs::exists(gtfsFile) && ! bfs::symbolic_link_exists(gtfsFile)) {
                BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                        "Could not find specific GTF file: ") + gtfsFile));
            }
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
        
        cout << "Loading Genomic GFF file" << endl;
        GFF::load(GFF3, genomicGffFile, *genomicGffs);
        
        cout << "Loading Transcript GFF file" << endl;
        GFF::load(GFF3, transcriptGffFile, transdecoderCdsGffs);
        
        cout << "Loading Cufflinks GTF file" << endl;
        GFF::load(GTF, gtfsFile, gtfs);
        
        cout << "Loading Full Lengther DB Annot file" << endl;
        DBAnnot::load(dbAnnotFile, flnDbannots);
        
        cout << "Loading Full Lengther New Coding file" << endl;
        DBAnnot::load(ncFile, flnNc);

    }
    
    void createMaps() {
        
        auto_cpu_timer timer(1, "Total indexing time: %ws\n\n");
        
        cout << "Creating indices" << endl
             << "----------------" << endl << endl;
        
        // Index genomic GFFs by Id
        BOOST_FOREACH(shared_ptr<GFF> gff, *genomicGffs) {
            if (gff->GetType() == MRNA) {
                maps.genomicGffMap[gff->GetRootId()] = gff;
            }
        }
        
        cout << "Indexed " << maps.genomicGffMap.size() << " mRNAs from genomic GFF file keyed to Root ID" << endl;
        
        // Index transdecoder CDSes
        BOOST_FOREACH(shared_ptr<GFF> gff, transdecoderCdsGffs) {        
            if (gff->GetType() == CDS) {            
                maps.transdecoderCdsGffMap[gff->GetSeqId()] = gff;
            }
        }
        
        cout << "Indexed " << maps.transdecoderCdsGffMap.size() << " CDSes from transcript GFF file keyed to Root ID" << endl;
        
        if (!gtfsFile.empty()) {
            // Index cufflinks transcripts
            BOOST_FOREACH(shared_ptr<GFF> gff, gtfs) {        
                if (gff->GetType() == TRANSCRIPT) {            
                    maps.gtfMap[gff->GetRootTranscriptId()] = gff;
                }
            }
        
            cout << "Indexed " << maps.gtfMap.size() << " cufflinks transcripts" << endl;
        }
        
        // Index fln cdss by id
        BOOST_FOREACH(shared_ptr<DBAnnot> db, flnDbannots) {
            maps.allDistinctFlnCds[db->GetId()] = db;
            if (db->GetStatus() == COMPLETE) {
                maps.uniqFlnCds[db->GetId()] = db;                
            }
        }
        
        cout << "Indexed " << maps.uniqFlnCds.size() << " complete and known transcripts from full lengther" << endl;
        
        // Index fln cdss by id
        BOOST_FOREACH(shared_ptr<DBAnnot> db, flnNc) {
            maps.uniqFlnNcCds[db->GetId()] = db;
            maps.allDistinctFlnCds[db->GetId()] = db;
        }
        
        cout << "Indexed " << maps.uniqFlnNcCds.size() << " full lengther new coding transcripts" << endl;
                
        cout << "Indexed " << maps.allDistinctFlnCds.size() << " total full lengther transcripts" << endl;
    }
    
    void filter(std::vector< shared_ptr<GFFList> >& stages) {
        
        auto_cpu_timer timer(1, "Total filtering time: %ws\n\n");
        
        cout << "Filtering transcripts" << endl
             << "---------------------" << endl << endl;
        
        std::vector< shared_ptr<TranscriptFilter> > filters;
        
        filters.push_back(shared_ptr<TranscriptFilter>(new MultipleOrfFilter()));
        filters.push_back(shared_ptr<TranscriptFilter>(new InconsistentCoordsFilter(include, cds)));
        filters.push_back(shared_ptr<TranscriptFilter>(new MultipleTranscriptFilter()));
        filters.push_back(shared_ptr<TranscriptFilter>(new StrandFilter()));
                        
        stages.push_back(genomicGffs);
        
        for(int i = 0; i < filters.size(); i++) {
            
            stages.push_back(shared_ptr<GFFList>(new GFFList()));
            shared_ptr<GFFList> in = stages[i];
            shared_ptr<GFFList> out = stages[i+1];
            
            cout << "Executing filter " << i+1 << " of " << filters.size() << endl
                 << "Name: " << filters[i]->getName() << endl
                 << "Description: " << filters[i]->getDescription() << endl
                 << "Filter input contains " << in->size() << " GFF records" << endl;
        
            
            // Do the filtering for this stage
            filters[i]->filter(*in, maps, *out);            
        
            // Record how many entries have been filtered
            size_t diff = in->size() - out->size();
            
            cout << "Report: " << endl
                 << filters[i]->getReport() << endl;
            
            cout << "Filtered out " << diff << " GFF records" << endl
                 << "Output contains " << out->size() << " GFF records" << endl 
                 << "Filter " << i+1 << " of " << filters.size() << " completed" << endl << endl;
            
            // Output filtered GFF for this stage if requested
            if (outputAllStages) {
                
                setSourceForOutput(*out, "gts");
                
                std::stringstream ss;
                ss << outputPrefix << ".stage." << (i+1) << ".gff3";
                const string stageOut = ss.str();
                gts::GFF::save(stageOut, *out);
            }
            
            cout << "--------------------------------------" << endl << endl;           
            
            
        }
        
    }
    
    void setSourceForOutput(GFFList& gffs, string source) {
        
        BOOST_FOREACH(shared_ptr<GFF> gff, gffs) {
            gff->SetSource(source);
        }
    }
    
    void output(GFFList& gffs) {        
        
        auto_cpu_timer timer(1, "Total writing time: %ws\n\n");
        
        // Save passed output
        std::stringstream ssp;
        ssp << outputPrefix << ".pass.gff3";
        const string passOut = ssp.str(); 
        std::ofstream pass(passOut.c_str());
        
        // Save passed output
        std::stringstream ssf;
        ssf << outputPrefix << ".fail.gff3";
        const string failOut = ssf.str();        
        std::ofstream fail(failOut.c_str());
        
        cout << "--------------------------------------" << endl << endl
             << "Re-processing: " << genomicGffFile << endl
             << "Splitting file based on transcripts that passed all the filters" << endl
             << "Passed: " << passOut << endl
             << "Failed: " << failOut << endl << endl;
        
        // Set the output source
        setSourceForOutput(gffs, "gts");
        
        // Sort the output
        std::sort(gffs.begin(), gffs.end(), GFFOrdering());        
        
        // Create list of failed transcripts
        GFFIdMap passed;
        BOOST_FOREACH(shared_ptr<GFF> gff, gffs) {
            passed[gff->GetRootId()] = gff;
        } 
        
        std::ifstream in(genomicGffFile.c_str());
        std::string line; 
        while (std::getline(in, line)) {            
            boost::trim(line);
            if (!line.empty()) {
                shared_ptr<GFF> gff = GFF::parse(GFF3, line);
                
                gff->SetSource("gts");
                
                if (passed.count(gff->GetRootId()) > 0) {
                    gff->write(pass);
                }
                else {
                    gff->write(fail);
                }
            }
        }
        in.close();        
        
        pass.close();
        fail.close();
    }
    
public :
    
    GTS(const string& genomicGffFile, const string& transcriptGffFile, const string& flnDir) :
        genomicGffFile(genomicGffFile), transcriptGffFile(transcriptGffFile), flnDir(flnDir),
                outputPrefix("gts_out"), gtfsFile(""), cds(0.5), include(false), outputAllStages(false), verbose(false)
    {
        genomicGffs = shared_ptr<GFFList>(new GFFList());
    }
    
    string getOutputPrefix() const
    {
        return outputPrefix;
    }

    void setOutputPrefix(string outputPrefix)
    {
        this->outputPrefix = outputPrefix;
    }

    string getGTFFile() const
    {
        return gtfsFile;
    }

    void setGTFFile(string gtfFile)
    {
        this->gtfsFile = gtfFile;
    }

        
    double getCds() const
    {
        return cds;
    }

    void setCds(double cds)
    {
        this->cds = cds;
    }

    bool isInclude() const
    {
        return include;
    }

    void setInclude(bool include)
    {
        this->include = include;
    }

    bool isOutputAllStages() const
    {
        return outputAllStages;
    }

    void setOutputAllStages(bool outputAllStages)
    {
        this->outputAllStages = outputAllStages;
    }

    bool isVerbose() const
    {
        return verbose;
    }

    void setVerbose(bool verbose)
    {
        this->verbose = verbose;
    }


    virtual ~GTS() {}
    
    
    void execute() {
        
        auto_cpu_timer timer(1, "Total execution time: %ws\n");
        
        // Load the input data
        load();
        
        // Create maps to assist filters
        createMaps();
        
        // Do the filtering
        std::vector< shared_ptr<GFFList> > stages;
        filter(stages);
        
        // Output consolidated GFFs
        output(*(stages[stages.size() - 1]));
        
        cout << "--------------------------------------" << endl << endl;  
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
        string gtfsFile;
        string flnResultsDir;
        string outputPrefix;
        bool include;
        double cds;
        bool outputAllStages;
                
        string cufflinksGtfFile;
        
        bool verbose;
        bool version;
        bool help;

        // Declare the supported options.
        po::options_description generic_options(helpHeader());
        generic_options.add_options()
                ("ggff,g", po::value<string>(&genomicGffFile), 
                    "GFF3 file containing the genomic coordinates for the transcript features.")
                ("tgff,t", po::value<string>(&transcriptGffFile),
                    "GFF3 file containing the transcript coordinates for the transcript features.")
                ("gtf", po::value<string>(&gtfsFile),
                    "GTF file containing transcripts.")
                ("fln_dir,f", po::value<string>(&flnResultsDir),
                    "Full lengther results directory, containing the \"dbannotated.txt\" and \"new_coding.txt\" files.")
                ("output,o", po::value<string>(&outputPrefix)->default_value(string("gts_out")),
                    "The output prefix for all output files generated.")
                ("include", po::bool_switch(&include)->default_value(false), 
                    "Include transcripts with no full lengther homology hit, providing it has a full lengther new_coding hit.")
                ("cds", po::value<double>(&cds)->default_value(0.5), 
                    "Min percentage length of CDS relative to mRNA for hits with homology.  0.0 -> 1.0")
                ("all,a", po::bool_switch(&outputAllStages)->default_value(false), 
                    "Whether or not to output GFF entries filtered at each stage.")
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
        
        
        // Select good transcripts
        gts::GTS gts(genomicGffFile, transcriptGffFile, flnResultsDir);
        gts.setOutputPrefix(outputPrefix);
        gts.setGTFFile(gtfsFile);
        gts.setCds(cds);
        gts.setInclude(include);
        gts.setOutputAllStages(outputAllStages);
        gts.setVerbose(verbose);
        
        gts.execute();        
                
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


