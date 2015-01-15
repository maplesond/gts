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

#pragma once

#include <sstream>
using std::stringstream;        

#include "../gff.hpp"
#include "transcript_filter.hpp"
using gts::gff::GFFModel;


namespace gts {

typedef boost::unordered_map<string, uint32_t> IdCounter;
    
class MultipleOrfFilter : public TranscriptFilter {
      
public:
    
    MultipleOrfFilter() : TranscriptFilter() {}
    
    ~MultipleOrfFilter() {}
    
    string getName() {
        return string("Multiple ORF Filter");
    }
    
    string getDescription() {
        return string("Keeps transcripts with a single ORF (in transdecoder terms, this means 1 transcript per gene/locus) and at least one 5' and 3' UTR");
    }
    
    
protected:    
    
    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        uint32_t inGenes = in.getGeneList()->size();
        uint32_t inOrfs = 0;
        uint32_t in5UTRs = 0;
        uint32_t in3UTRs = 0;
        
        uint32_t outOrfs = 0;
        uint32_t out5UTRs = 0;
        uint32_t out3UTRs = 0;
        
        // Count duplicate mRNAs, 5' UTRs and 3' UTRs
        BOOST_FOREACH(GFFPtr gene, *(in.getGeneList())) {
            
            // Gene should only have mRNAs (transcripts) as children
            uint32_t geneTranscripts = gene->GetNbChildren();
            
            if (geneTranscripts == 0) {
                cerr << "Found gene with no transcripts: " << gene->GetId() << endl;
            }
            // Only interested in genes that have a single transcript
            else {
            
                uint32_t geneOrfs = 0;
                uint32_t gene5UTRs = 0;
                uint32_t gene3UTRs = 0;

                GFFList goodTranscripts;
                uint32_t goodOrfs = 0;
                uint32_t good5UTRs = 0;
                uint32_t good3UTRs = 0;

                BOOST_FOREACH(GFFPtr transcript, *(gene->GetChildList())) {

                    // Gene should only have mRNAs (transcripts) as children
                    uint32_t transcriptOrfs = 0;
                    uint32_t transcript5UTRs = 0;
                    uint32_t transcript3UTRs = 0;

                    GFFListPtr children = transcript->GetChildList();

                    if (children) {
                        BOOST_FOREACH(GFFPtr gff, *children) {

                            switch(gff->GetType()) {
                                case CDS:
                                    transcriptOrfs++;
                                    geneOrfs++;
                                    inOrfs++;
                                    break;
                                case UTR5:
                                    transcript5UTRs++;
                                    gene5UTRs++;
                                    in5UTRs++;
                                    break;
                                case UTR3:
                                    transcript3UTRs++;
                                    gene3UTRs++;
                                    in3UTRs++;
                                    break;                            
                            }
                        } 
                    }
                    else {
                        BOOST_THROW_EXCEPTION(TranscriptFilterException() << TranscriptFilterErrorInfo(string(
                                "Invalid GFF.  Found a transcript with no children: ") + transcript->GetId()));
                    }

                    // We only want transcripts with at least 1 5' and 3' UTR
                    if (transcript5UTRs >= 1 && transcript3UTRs >= 1) {
                        goodTranscripts.push_back(transcript);
                        goodOrfs += transcriptOrfs;
                        good5UTRs += transcript5UTRs;
                        good3UTRs += transcript3UTRs;
                    }
                }

                // We only want genes with 1 transcript and at least 1 5' and 3' UTR
                if (geneTranscripts == 1 && goodTranscripts.size() >= 1) {                

                    outOrfs += goodOrfs;
                    out5UTRs += good5UTRs;
                    out3UTRs += good3UTRs;

                    // Copy gene without child info
                    GFFPtr newGene = make_shared<GFF>(*gene);

                    BOOST_FOREACH(GFFPtr goodTranscript, goodTranscripts) {
                        newGene->addChild(goodTranscript);
                    }

                    out.addGene(newGene);
                }
            }
        }

        
        stringstream ss;
        
        ss << " - # Genes: " << out.getNbGenes() << " / " << in.getNbGenes() << endl
           << " - # Transcripts (mRNA): " << out.getTotalNbTranscripts() << " / " << in.getTotalNbTranscripts() << endl     
           << " - # ORFs (CDS): " << outOrfs << " / " << inOrfs << endl
           << " - # 5' UTRs: " << out5UTRs << " / " << in5UTRs << endl
           << " - # 3' UTRs: " << out3UTRs << " / " << in3UTRs << endl;
        
        report = ss.str();
    }
};

}