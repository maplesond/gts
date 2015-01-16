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
#include <boost/unordered_set.hpp>
using boost::to_upper_copy;
using boost::timer::auto_cpu_timer;
using boost::shared_ptr;
using boost::unordered_map;
using boost::unordered_set;
namespace po = boost::program_options;
namespace bfs = boost::filesystem;

#include "gff.hpp"
#include "genbank.hpp"
#include "fln.hpp"
#include "filters/transcript_filter.hpp"
#include "filters/multiple_orf_filter.hpp"
#include "filters/inconsistent_coords_filter.hpp"
#include "filters/strand_filter.hpp"
#include "filters/overlap_filter.hpp"
using gts::gff::GFF;
using gts::gff::GFFModel;

const double DEFAULT_CDS_LEN_RATIO = 0.4;
const double DEFAULT_CDNA_LEN_RATIO = 0.5;
const uint32_t DEFAULT_WINDOW_SIZE = 1000;

namespace gts {

typedef boost::error_info<struct GTSError,string> GTSErrorInfo;
struct GTSException: virtual boost::exception, virtual std::exception { };

typedef std::vector< boost::shared_ptr<DBAnnot> > FLNDBAnnotList;


class GTS {
    
private:    
    
    GFFModelPtr genomicGffModel;
    GFFModelPtr alignmentGffModel;
    GFFModelPtr genomicGffModelFixed;
    GFFModelPtr alignmentGffModelFixed;
    FLNDBAnnotList flnDbannots;
    FLNDBAnnotList flnNc;
    
    Maps maps;
    
    const string genomicGffFile;
    const string transcriptGffFile;
    const string flnDir;
    
    string outputPrefix;
    string gtfsFile;
    double cdsLenRatio;
    double cdnaLenRatio;
    bool include;
    uint32_t windowSize;
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
        this->genomicGffModel = GFFModel::load(genomicGffFile);
        
        cout << endl << "Loading Cluster Alignment GFF file" << endl;
        this->alignmentGffModel = GFFModel::load(transcriptGffFile);
        
        cout << endl <<"Loading GTF file" << endl;
        GFFListPtr gtfs = make_shared<GFFList>();
        GFF::load(GTF, gtfsFile, *gtfs);
        uint32_t nbGtfTranscripts = 0;
        BOOST_FOREACH(GFFPtr gff, *gtfs) {        
            if (gff->GetType() == TRANSCRIPT) {
                const string transcriptRootId = gff->GetRootTranscriptId();
                maps.gtfMap[transcriptRootId] = gff;
                nbGtfTranscripts++;
            }
        }
        cout << " = Indexed " << maps.gtfMap.size() << " distinct GTF transcripts" << endl;
        
        cout << endl << "Loading Full Lengther DB Annot file" << endl;
        DBAnnot::load(dbAnnotFile, flnDbannots);
        
        cout << "Loading Full Lengther New Coding file" << endl;
        DBAnnot::load(ncFile, flnNc);

    }
    
    void createMaps() {
        
        auto_cpu_timer timer(1, "Total indexing time: %ws\n\n");
        
        cout << "Creating indices" << endl
             << "----------------" << endl << endl;
        
        // Index transdecoder CDSes
        GFFListPtr cdses = alignmentGffModel->getAllOfType(CDS);
        BOOST_FOREACH(GFFPtr cds, *cdses) {        
            maps.transdecoderCdsGffMap[cds->GetId()] = cds;            
        }
        
        cout << "Indexed " << maps.transdecoderCdsGffMap.size() << " CDSes from transcript GFF file keyed to Target ID" << endl;
        
        // Index transdecoder CDNAs
        GFFListPtr exons = alignmentGffModel->getAllOfType(EXON);
        BOOST_FOREACH(GFFPtr cdna, *exons) {        
            maps.transdecoderCDNAGffMap[cdna->GetId()] = cdna;            
        }
        
        cout << "Indexed " << maps.transdecoderCDNAGffMap.size() << " CDNAs (exons) from transcript GFF file keyed to Target ID" << endl;
        
        
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
    
    /**
     * Probably input from transdecoder contains transcripts that have a unique 
     * gene assigned.  What we need though is each gene to contain multiple transcripts
     * if those transcripts live at the same locus.
     * @param genes
     */
    GFFModelPtr resolveGeneModel(GFFModel& genes, bool genomicCoords) {
        
        const size_t nbGenes = genes.getNbGenes();
        const size_t nbTranscripts = genes.getTotalNbTranscripts();
        
        GFFModelPtr newGeneModel = make_shared<GFFModel>();
        
        if (nbGenes < nbTranscripts) {
            cout << " - Gene count and transcript count are already different.  Skipping step" << endl;
            
            BOOST_FOREACH(GFFPtr gene, *genes.getGeneList()) {
                newGeneModel->addGene(gene);
            }
        }
        else if (nbGenes > nbTranscripts) {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                        "Corrupt gene model.  Gene model contains more genes than transcripts.")));
        }
        else if (genes.getNbGenes() >= 2) {
            
            cout << " - Combining transcripts" << endl;
            
            GFFIdMap newGeneMap;
            
            // Change gene name to whatever was found in the GTF file
            BOOST_FOREACH(GFFPtr thisGene, *genes.getGeneList()) {
                
                GFFListPtr transcripts = thisGene->GetChildList();
                 
                const string rootId = thisGene->GetRootId();                
                
                string gtfGeneId = maps.gtfMap[rootId]->GetGeneId();
                string gtfTranscriptId = maps.gtfMap[rootId]->GetTranscriptId();
                
                // Check to see if this gene cane be merged with the last one
                if (newGeneMap.count(gtfGeneId)) {
                
                    GFFPtr lastGene = newGeneMap[gtfGeneId];
                    
                    // Sanity check
                    if (genomicCoords && !boost::iequals(thisGene->GetSeqId(), lastGene->GetSeqId())) {
                        
                        BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                            "Genes with same Id are not on the same target sequence: Gene Id: ") 
                                + thisGene->GetId() + " - Target Seq Id: " + thisGene->GetSeqId()));                        
                    }
                    
                    // Move transcripts across to last Gene
                    BOOST_FOREACH(GFFPtr transcript, *transcripts) {

                        // Update this transcript's parent gene id.
                        transcript->SetParentId(gtfGeneId);
                        transcript->SetAlias(gtfTranscriptId);
                        
                        // Fix CDS Ids
                        fixCDSIds(transcript);
                        
                        // Add this transcript to the existing gene
                        lastGene->addChild(transcript);
                        
                        // Extend start coord if necessary
                        if (transcript->GetStart() < lastGene->GetStart()) {
                            lastGene->SetStart(transcript->GetStart());
                        }
                        
                        // Extend end coord if necessary 
                        if (transcript->GetEnd() > lastGene->GetEnd()) {
                            lastGene->SetEnd(transcript->GetEnd());
                        }
                    }
                }
                else {                                    
                    
                    if (gtfGeneId.empty()) {
                        BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                            "Could not find transcript assembly id in GTF file: ") + rootId));
                    }
                    
                    // Set the new gene id
                    thisGene->SetId(gtfGeneId);
                    
                    newGeneMap[gtfGeneId] = thisGene;
                    
                    // Update the transcripts' parent gene id.
                    BOOST_FOREACH(GFFPtr transcript, *transcripts) {
                        transcript->SetParentId(gtfGeneId);
                        transcript->SetAlias(gtfTranscriptId);
                        fixCDSIds(transcript);
                    }                                       
                }
            }
            
            BOOST_FOREACH(GFFIdMap::value_type& i, newGeneMap) {
                newGeneModel->addGene(i.second);
            }            
        }
        else {
            BOOST_THROW_EXCEPTION(GTSException() << GTSErrorInfo(string(
                        "Must have at least 2 or more genes to rebuild gene model")));
        }
        
        // We should have a proper gene model now, so now just run through and update
        // attributes so they have sensible values
        BOOST_FOREACH(GFFPtr g, *(newGeneModel->getGeneList())) {
            g->SetName(g->GetId());
            
            BOOST_FOREACH(GFFPtr t, *(g->GetChildList())) {
                
                t->SetName(t->GetId());
                t->SetParentId(g->GetId());
                t->SetNote(g->GetId());                
            }
        }
        
        return newGeneModel;        
    }
    
    /**
     * Resolve CDS ids (convert from "cds.xxxx" to xxxx.cds<n>)
     * @param transcript The transcript containing CDSs to fix
     */
    void fixCDSIds(GFFPtr transcript) {
        GFFListPtr cdses = transcript->GetAllOfType(gts::gff::CDS);                
        uint16_t cdsIndex = 1;
        BOOST_FOREACH(GFFPtr cds, *cdses) {

            const string cdsid = cds->GetId();

            if (cdsid.find("cds.") == 0) {

                stringstream ss;
                ss << cdsid.substr(4) << ".cds" << cdsIndex++;

                cds->SetId(ss.str());
            }                    
        }
    }
    
    void resolveGenes() {
        
        auto_cpu_timer timer(1, "Total resolving time: %ws\n\n");
        
        cout << "Resolving gene model using GTF gene names" << endl
             << "-----------------------------------------" << endl << endl
             << "Processing genomic gene model" << endl;
        
        this->genomicGffModelFixed = resolveGeneModel(*(this->genomicGffModel), true);
        
        cout << " = Gene model contains " << this->genomicGffModelFixed->getNbGenes() 
                << " genes and " << this->genomicGffModelFixed->getTotalNbTranscripts() 
                << " transcripts" << endl << endl;
        
        cout << "Processing cluster alignment gene model" << endl;        
        
        this->alignmentGffModelFixed = resolveGeneModel(*(this->alignmentGffModel), false);        
        
        cout << " = Gene model contains " << this->alignmentGffModelFixed->getNbGenes() 
                << " genes and " << this->alignmentGffModelFixed->getTotalNbTranscripts() 
                << " transcripts" << endl << endl;
    }
    
    void filter(std::vector< shared_ptr<GFFModel> >& stages) {
        
        auto_cpu_timer timer(1, "Total filtering time: %ws\n\n");
        
        cout << "Filtering genomic gene model" << endl
             << "----------------------------" << endl << endl;
        
        vector< shared_ptr<TranscriptFilter> > filters;
        
        filters.push_back(make_shared<MultipleOrfFilter>());
        filters.push_back(make_shared<InconsistentCoordsFilter>(include, cdsLenRatio, cdnaLenRatio));
        filters.push_back(make_shared<StrandFilter>());
        filters.push_back(make_shared<OverlapFilter>(windowSize, genomicGffModelFixed));
                        
        stages.push_back(genomicGffModelFixed);
        
        for(int i = 0; i < filters.size(); i++) {
            
            stages.push_back(make_shared<GFFModel>());
            shared_ptr<GFFModel> in = stages[i];
            shared_ptr<GFFModel> out = stages[i+1];
            
            cout << "Executing filter " << i+1 << " of " << filters.size() << endl
                 << "Name: " << filters[i]->getName() << endl
                 << "Description: " << filters[i]->getDescription() << endl
                 << "Filter input contains " << in->getNbGenes() << " genes and " << in->getTotalNbTranscripts() << " transcripts" << endl;
                    
            // Do the filtering for this stage
            filters[i]->filter(*in, maps, *out);            
        
            // Record how many entries have been filtered
            size_t geneDiff = in->getNbGenes() - out->getNbGenes();
            size_t transcriptDiff = in->getTotalNbTranscripts() - out->getTotalNbTranscripts();
            
            cout << "Report: " << endl
                 << filters[i]->getReport() << endl;
            
            cout << "Filtered out " << geneDiff << " genes and " << transcriptDiff << " transcripts" << endl
                 << "Output contains " << out->getNbGenes() << " genes and " << out->getTotalNbTranscripts() << " transcripts" << endl 
                 << "Filter " << i+1 << " of " << filters.size() << " completed" << endl << endl;
            
            // Output filtered GFF for this stage if requested
            if (outputAllStages) {
                
                std::stringstream ss;
                ss << outputPrefix << ".stage." << (i+1) << ".gff3";
                const string stageOut = ss.str();
                out->save(stageOut, true, string("gts"));
            }
            
            cout << "--------------------------------------" << endl << endl;
        }
    }
    
    
    void output(GFFModel& goodGeneModel) {        
        
        auto_cpu_timer timer(1, "Total writing time: %ws\n\n");
        
        // Save passed output
        std::stringstream ssp;
        ssp << outputPrefix << ".pass.gff3";
        const string passOut = ssp.str(); 
        
        // Save passed output
        std::stringstream ssf;
        ssf << outputPrefix << ".fail.gff3";
        const string failOut = ssf.str();        
        std::ofstream fail(failOut.c_str());
        
        cout << "--------------------------------------" << endl << endl
             << "Saving final output" << endl
             << "-------------------" << endl << endl
             << "Re-processing: " << genomicGffFile << endl
             << "Splitting file based on transcripts that passed all the filters" << endl << endl;
        
        // Output passed results
        goodGeneModel.save(passOut, true, "gts");
                
        uint32_t failGeneCount = 0;
        uint32_t failTranscriptCount = 0;
        
        // Output failed results
        BOOST_FOREACH(shared_ptr<GFF> gene, *(genomicGffModel->getGeneList())) {
            
            string id = gene->GetId();
            
            if (goodGeneModel.containsGene(id)) {
               
                GFFList failedTranscripts;
                
                BOOST_FOREACH(shared_ptr<GFF> transcript, *(gene->GetChildList())) {
                    
                    if (!goodGeneModel.containsTranscript(transcript->GetId())) {
                        failedTranscripts.push_back(transcript);
                    }
                }
                
                if (!failedTranscripts.empty()) {
                    
                    // Write out just the gene contents
                    gene->write(fail, "gts", false);
                    failGeneCount++;
                
                    // Now only write out the failed transcripts
                    BOOST_FOREACH(shared_ptr<GFF> failedTranscript, failedTranscripts) {
                        failedTranscript->write(fail, "gts", true);
                        failTranscriptCount++;
                    }
                }                
            }
            else {
                // The gene wasn't in the passed gene model, so just write out the complete gene and child entries
                gene->write(fail, "gts", true);
                failGeneCount++;
                failTranscriptCount += gene->GetChildList()->size();
            }
            
            // Write a gap between the genes
            fail << endl;            
        }
        
        fail.close();
        
        cout << "Processed " << genomicGffModel->getNbGenes() << " genes and " << genomicGffModel->getTotalNbTranscripts() << " transcripts" << endl
             << "Sent " << goodGeneModel.getNbGenes() << " genes and " << goodGeneModel.getTotalNbTranscripts() << " transcripts to " << passOut << endl
             << "Sent " << failGeneCount << " genes and " << failTranscriptCount << " transcripts to " << failOut << endl << endl
             << "NOTE: the sum of passed and failed gene counts may exceed the number of processed genes due to multi-transcript genes." << endl << endl;
    }
    
public :
    
    GTS(const string& genomicGffFile, const string& transcriptGffFile, const string& flnDir) :
        genomicGffFile(genomicGffFile), transcriptGffFile(transcriptGffFile), flnDir(flnDir),
                outputPrefix("gts_out"), gtfsFile(""), 
                cdsLenRatio(DEFAULT_CDS_LEN_RATIO), cdnaLenRatio(DEFAULT_CDNA_LEN_RATIO), 
                include(false), windowSize(DEFAULT_WINDOW_SIZE), outputAllStages(false), verbose(false)
    {
        genomicGffModel = make_shared<GFFModel>();
        alignmentGffModel = make_shared<GFFModel>();
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

        
    double getCdsLenRatio() const
    {
        return cdsLenRatio;
    }

    void setCdsLenRatio(double cdsLenRatio)
    {
        this->cdsLenRatio = cdsLenRatio;
    }
    
    double getCdnaLenRatio() const
    {
        return cdnaLenRatio;
    }

    void setCdnaLenRatio(double cdnaLenRatio)
    {
        this->cdnaLenRatio = cdnaLenRatio;
    }


    bool isInclude() const
    {
        return include;
    }

    void setInclude(bool include)
    {
        this->include = include;
    }
    
    uint32_t getWindowSize() const
    {
        return windowSize;
    }

    void setWindowSize(uint32_t windowSize)
    {
        this->windowSize = windowSize;
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
        
        // Resolve genes so that we have one gene and one or more transcripts
        resolveGenes();
        
        // Create maps to assist filters
        createMaps();
        
        // Do the filtering
        vector<GFFModelPtr> stages;
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
        uint32_t windowSize;
        double cdsLenRatio;
        double cdnaLenRatio;
        bool outputAllStages;
                
        string cufflinksGtfFile;
        
        bool verbose;
        bool version;
        bool help;

        // Declare the supported options.
        po::options_description generic_options(helpHeader());
        generic_options.add_options()
                ("genomic_gff,g", po::value<string>(&genomicGffFile), 
                    "Transdecoder GFF3 file containing the genomic coordinates for the transcript features.")
                ("transcript_gff,t", po::value<string>(&transcriptGffFile),
                    "Transdecoder GFF3 file containing the transcript coordinates for the transcript features.")
                ("gtf", po::value<string>(&gtfsFile),
                    "GTF file containing transcripts.")
                ("fln_dir,f", po::value<string>(&flnResultsDir),
                    "Full lengther results directory, containing the \"dbannotated.txt\" and \"new_coding.txt\" files.")
                ("output,o", po::value<string>(&outputPrefix)->default_value(string("gts_out")),
                    "The output prefix for all output files generated.")
                ("include_putative,i", po::bool_switch(&include)->default_value(false), 
                    "Include putative transcripts, i.e. transcripts with a full lengther new_coding hit.")
                ("window_size,w", po::value<uint32_t>(&windowSize)->default_value(DEFAULT_WINDOW_SIZE), 
                    "The gap to enforce between genes.")
                ("cds_ratio", po::value<double>(&cdsLenRatio)->default_value(DEFAULT_CDS_LEN_RATIO), 
                    "Min percentage length of CDS content relative to full length transcripts for hits with homology.  0.0 -> 1.0")
                ("cdna_ratio", po::value<double>(&cdnaLenRatio)->default_value(DEFAULT_CDNA_LEN_RATIO), 
                    "Min percentage length of cDNA (exon) content relative to full length transcripts for hits with homology.  0.0 -> 1.0")
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
        gts.setCdsLenRatio(cdsLenRatio);
        gts.setCdnaLenRatio(cdnaLenRatio);
        gts.setInclude(include);
        gts.setWindowSize(windowSize);
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


