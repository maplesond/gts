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
#include "genbank.hpp"
using gts::gff::GFF;
using gts::gff::GFF3;
using gts::gb::Genbank;
using gts::gb::Feature;
using gts::gb::Property;

typedef boost::unordered_map<string, shared_ptr<GFF> > GFFIdMap;

string helpHeader() {
    return string("\nGenbank Filter Help.\n\n") +
                  "The genbank filter tool is used to filter a genbank file based on transcripts found in a provided GFF file\n\n" +
                  "Usage: gbfilter [options] -gb <genbank file> -g <gff file> -o <output file>\n\n" +
                 "\nAvailable options";
}


/**
 * Start point for portculis.
 */
int main(int argc, char *argv[]) {
    
    try {
        // Portculis args
        string passGffFile;
        string genbankFile;
        string outputFile;
        
        bool version;
        bool help;

        // Declare the supported options.
        po::options_description generic_options(helpHeader());
        generic_options.add_options()
                ("pass_gff,p", po::value<string>(&passGffFile), 
                    "GFF file containing transcripts that should be kept in the genbank file")
                ("genbank,b", po::value<string>(&genbankFile),
                    "The genbank file to filter")
                ("out,o", po::value<string>(&outputFile),
                    "The output genbank file")
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
        
        
        // Load genbank file to filter
        cout << "Loading genbank file" << endl;
        vector<shared_ptr<Genbank> > genbank;
        gts::gb::Genbank::load(genbankFile, genbank);
        
        // Index genomic GFFs by Id
        cout << "Filtering unsuitable genbank records" << endl;
        vector<shared_ptr<Genbank> > genbankFiltered;
        BOOST_FOREACH(shared_ptr<Genbank> gb, genbank) {
            
            uint16_t nbMrnas = 0;
            uint16_t nbCDS = 0;
            BOOST_FOREACH(shared_ptr<Feature> f, gb->getFeatures()->features) {
                if (boost::iequals(f->type, "mRNA")) {
                    nbMrnas++;
                }
                else if (boost::iequals(f->type, "CDS")) {
                    nbCDS++;
                }
            }
            
            if (nbMrnas == 1 && nbCDS == 1) {
                genbankFiltered.push_back(gb);                
            }
        }
        cout << " = Keeping " << genbankFiltered.size() << " out of " << genbank.size() << " genbank records" << endl << endl;
        
        cout << "Loading GFF file" << endl;
        vector<shared_ptr<GFF> > gffs;
        GFF::load(GFF3, passGffFile, gffs, gts::gff::MRNA);
        
        // Index genomic GFFs by Id
        cout << "Indexing GFF file" << endl;
        GFFIdMap gffAsmblMap;
        GFFIdMap gffMrnaMap;
        GFFIdMap gffIdMap;
        BOOST_FOREACH(shared_ptr<GFF> gff, gffs) {
            vector<string> idElements;
            string id(gff->GetId());
            boost::split( idElements, id, boost::is_any_of("|"), boost::token_compress_on );                

            gffAsmblMap[idElements[0]] = gff;
            gffMrnaMap[idElements[0]] = gff;
            gffIdMap[gff->GetId()] = gff;
        }
        cout << " = done" << endl << endl;
        
        // Keeping only 
        cout << "Cross checking genbank records with GFF.  Keeping matches." << endl;
        vector<shared_ptr<Genbank> > genbankFiltered2;
        BOOST_FOREACH(shared_ptr<Genbank> gb, genbank) {
            bool done = false;
            BOOST_FOREACH(shared_ptr<Feature> f, gb->getFeatures()->features) {
                if (boost::iequals(f->type, "CDS")) {
                   BOOST_FOREACH(shared_ptr<Property> p, f->properties) {
                       if (boost::iequals(p->name, "gene")) {                           
                           if (gffIdMap.count(p->value) != 0) {
                               genbankFiltered2.push_back(gb);
                               done = true;
                               break;
                           }
                       }                                              
                   }
                   
                   if (done) {
                       break;
                   }
                }
            }
        }
        cout << " - Keeping " << genbankFiltered2.size() << " out of " << genbankFiltered.size() << " genbank records" << endl;
        
                
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


