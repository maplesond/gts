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

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <boost/shared_ptr.hpp>
using boost::timer::auto_cpu_timer;
using boost::shared_ptr;
namespace po = boost::program_options;
namespace bfs = boost::filesystem;

#include "gff.hpp"
using gts::gff::GFF;
using gts::gff::GFF3;
using gts::gff::GFFErrorInfo;
using gts::gff::GFFException;

typedef std::vector<shared_ptr<GFF> > GFFList;


string helpHeader() {
    return string("\nGFF ID Extractor Help.\n\n") +
                  "The gffids tool is used to extract mRNA IDs, and their associated parent gene IDs to a tab separated file\n\n" +
                  "Usage: gffids [options] -i <gff file> -o <tsv file>\n\n" +
                 "\nAvailable options";
}

void output(string& outputFile, GFFList& mRNAs) {
    
    auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
    
    cout << " - Saving to: " << outputFile << endl;
    
    std::ofstream file(outputFile.c_str());
        
    BOOST_FOREACH(shared_ptr<GFF> gff, mRNAs) {
        bool done = false;
        
        string id = gff->GetId();
        string parent = gff->GetParentId();
        
        if (id.empty()) {
            BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                        "mRNA does not contain ID entry")));
        }
        
        if (parent.empty()) {
            BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                        "mRNA does not contain Parent entry")));
        }
        
        file << parent << "\t" << id << std::endl;        
    }
    
    file.close();
}


/**
 * Start point for portculis.
 */
int main(int argc, char *argv[]) {
    
    try {
        // Portculis args
        string inputFile;
        string outputFile;
        
        bool version;
        bool help;

        // Declare the supported options.
        po::options_description generic_options(helpHeader());
        generic_options.add_options()
                ("input,i", po::value<string>(&inputFile), 
                    "The input GFF file to extract IDs from")
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
        
        cout << "Loading mRNA entries from GFF file" << endl;
        vector<shared_ptr<GFF> > gffs;
        GFF::load(GFF3, inputFile, gffs, gts::gff::MRNA);
        
        cout << "Writing IDs to output" << endl;        
        output(outputFile, gffs);        
        
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


