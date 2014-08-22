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
        return string("Filters out transcripts with multiple ORFs and no 5' and 3' UTRs");
    }
    
    
protected:    
    
    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        uint32_t inGenes = in.getGeneList()->size();
        uint32_t inTranscripts = 0;
        uint32_t inOrfs = 0;
        uint32_t in5UTRs = 0;
        uint32_t in3UTRs = 0;
        
        uint32_t outGenes = 0;
        uint32_t outTranscripts = 0;
        uint32_t outOrfs = 0;
        uint32_t out5UTRs = 0;
        uint32_t out3UTRs = 0;
        
        // Count duplicate mRNAs, 5' UTRs and 3' UTRs
        BOOST_FOREACH(shared_ptr<GFF> gene, *(in.getGeneList())) {
            
            // Gene should only have mRNAs (transcripts) as children
            uint32_t geneTranscripts = gene->GetNbChildren();
            
            inTranscripts += geneTranscripts;
            
            uint32_t geneOrfs = 0;
            uint32_t gene5UTRs = 0;
            uint32_t gene3UTRs = 0;
            
            vector<shared_ptr<GFF> > goodTranscripts;
            uint32_t goodOrfs = 0;
            uint32_t good5UTRs = 0;
            uint32_t good3UTRs = 0;
            
            BOOST_FOREACH(shared_ptr<GFF> transcript, *gene->GetChildList()) {
                
                // Gene should only have mRNAs (transcripts) as children
                uint32_t transcriptOrfs = 0;
                uint32_t transcript5UTRs = 0;
                uint32_t transcript3UTRs = 0;
                
                BOOST_FOREACH(shared_ptr<GFF> gff, *transcript->GetChildList()) {
                    
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
                
                // We only want transcripts with 1 CDS (ORF) and at least 1 5' and 3' UTR
                if (transcriptOrfs == 1 && transcript5UTRs >= 1 && transcript3UTRs >= 1) {
                    goodTranscripts.push_back(transcript);
                    goodOrfs += transcriptOrfs;
                    good5UTRs += transcript5UTRs;
                    good3UTRs += transcript3UTRs;
                }
            }
            
            // We only want genes with 1 or more good transcript
            if (goodTranscripts.size() >= 1) {                
               
                outGenes++;
                outTranscripts++;
                outOrfs += goodOrfs;
                out5UTRs += good5UTRs;
                out3UTRs += good3UTRs;
                
                // Copy gene without child info
                shared_ptr<GFF> newGene = shared_ptr<GFF>(new GFF(*gene));
                
                BOOST_FOREACH(shared_ptr<GFF> goodTranscript, goodTranscripts) {
                    newGene->addChild(goodTranscript);
                }
                
                out.addGene(newGene);
            }
        }

        
        stringstream ss;
        
        ss << " - # Genes: " << outGenes << " / " << inGenes << endl
           << " - # Transcripts (mRNA): " << outTranscripts << " / " << inTranscripts << endl     
           << " - # ORFs (CDS): " << outOrfs << " / " << inOrfs << endl;
        
        report = ss.str();
    }
};

}