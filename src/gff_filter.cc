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
#include <vector>
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::vector;

#include <boost/unordered_set.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered/unordered_map.hpp>
using boost::timer::auto_cpu_timer;
using boost::shared_ptr;
using boost::unordered_set;
namespace po = boost::program_options;
namespace bfs = boost::filesystem;

#include "gff.hpp"
using gts::gff::GFF;
using gts::gff::GFF3;
using gts::gff::GFFErrorInfo;
using gts::gff::GFFException;

typedef std::vector<shared_ptr<GFF> > GFFList;

string helpHeader() {
    return string("\nGFF Filter Help.\n\n") +
                  "The gfffilter tool is used to filter out listed entries from the provided GFF file.  The tool will automatically try to determine parent and child relationships between the entries and filter those as well.\n\n" +
                  "Usage: gfffilter [options] -i <gff file> -l <list file> -o <gff file>\n\n" +
                 "\nAvailable options";
}

/**
 * Create gene map
 * @param input
 * @param gffMap
 */
void indexGenes(GFFList& input, unordered_set<string>& transcriptSet, unordered_set<string>& geneSet) {
    
    auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
    
    BOOST_FOREACH(shared_ptr<GFF> gff, input) {
        if (gff->GetType() == gts::gff::MRNA && transcriptSet.count(gff->GetId()) > 0) {
            geneSet.insert(gff->GetParentId());
        }
    }
    
    cout << " - Found " << geneSet.size() << " genes to exclude from GFF." << endl;
}


void loadEntries(string& path, unordered_set<string>& transcriptSet) {
    
    auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
    
    cout << " - Loading entry set: " << path << endl;

    std::ifstream file(path.c_str());
    std::string line; 
    while (std::getline(file, line)) {            
        boost::trim(line);
        if (!line.empty()) {
            transcriptSet.insert(line);
        }
    }
    file.close();
        
    cout << " - Loaded " << transcriptSet.size() << " entries from file." << endl;
}

void filter(GFFList& input, unordered_set<string>& geneSet, unordered_set<string>& transcriptSet, GFFList& output) {
    
    auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
    
    uint32_t keepGenes = 0;
    uint32_t keepTranscripts = 0;
    uint32_t totalGenes = 0;
    uint32_t totalTranscripts = 0;
    
    // Index genomic GFFs by Id
    BOOST_FOREACH(shared_ptr<GFF> gff, input) {

        if (gff->GetType() == gts::gff::GENE) {
            totalGenes++;                
        }
        else if (gff->GetType() == gts::gff::MRNA) {
            totalTranscripts++;
        }
        
        if (transcriptSet.count(gff->GetId()) > 0 ||      // Covers mRNA
            transcriptSet.count(gff->GetParentId()) > 0 ||  // Covers child entries of mRNA
            geneSet.count(gff->GetId()) > 0) {            // Covers genes
            // Do nothing
        }
        else {
            output.push_back(gff);
            
            if (gff->GetType() == gts::gff::GENE) {
                keepGenes++;                
            }
            else if (gff->GetType() == gts::gff::MRNA) {
                keepTranscripts++;
            }
        }
    }
    cout << " - Keeping " << output.size() << " out of " << input.size() << " GFF records" << endl;
    cout << " - Keeping " << keepGenes << " out of " << totalGenes << " genes" << endl;
    cout << " - Keeping " << keepTranscripts << " out of " << totalTranscripts << " transcripts" << endl;
}


/**
 * Start point for portculis.
 */
int main(int argc, char *argv[]) {
    
    try {
        // Portculis args
        string inputFile;
        string listFile;
        string outputFile;
        
        bool version;
        bool help;

        // Declare the supported options.
        po::options_description generic_options(helpHeader());
        generic_options.add_options()
                ("input,i", po::value<string>(&inputFile), 
                    "The input GFF file to filter")
                ("list,l", po::value<string>(&listFile),
                    "The list of GFF ids to filter")
                ("output,o", po::value<string>(&outputFile),
                    "The tab separated output file which will contain IDs")
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
        
        auto_cpu_timer timer(1, "\nTotal wall time taken: %ws\n\n");
        
        cout << endl << "Extracting mRNA IDs and associated parent gene IDs from GFF file" << endl << endl;
        
        cout << "Loading all GFF records file" << endl;
        vector<shared_ptr<GFF> > gffs;
        GFF::load(GFF3, inputFile, gffs);
        
        cout << "Loading transcript IDs (mRNAs) to filter" << endl;
        unordered_set<string> transcriptsToExclude;
        loadEntries(listFile, transcriptsToExclude);
        
        cout<< "Finding set of all genes to exclude GFF" << endl;
        unordered_set<string> genesToExclude;
        indexGenes(gffs, transcriptsToExclude, genesToExclude);
                
        cout << "Filtering listed entries from GFF" << endl;
        vector<shared_ptr<GFF> > filtered;
        filter(gffs, genesToExclude, transcriptsToExclude, filtered);
        
        cout << "Writing filtered GFF to disk" << endl;        
        gts::gff::GFF::save(outputFile, filtered);            
        
        cout << "Completed" << endl;
                
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


