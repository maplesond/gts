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

#include "transcript_filter.hpp"

using std::stringstream;        

namespace gts {

typedef boost::unordered_map<string, uint32_t> IdCounter;
    
class StrandFilter : public TranscriptFilter {
      
public:
    
    StrandFilter() : TranscriptFilter() {}
    
    ~StrandFilter() {}
    
    string getName() {
        return string("Strand Filter");
    }
    
    string getDescription() {
        return string("Filters out genomic transcripts which have an inconsistent strand when compared to the GTF file");
    }
    
    
protected:    
    
    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        BOOST_FOREACH(GFFPtr gene, *(in.getGeneList())) {
            
            GFFList goodTranscripts;
            const char gffGeneStrand = gene->GetStrand();            
            
            if (gffGeneStrand != '.') {
                
                GFFListPtr transcripts = gene->GetChildList();
                BOOST_FOREACH(GFFPtr transcript, *transcripts) {

                    const char gffTranscriptStrand = transcript->GetStrand();
                    const char gtfTranscriptStrand = maps.gtfMap[transcript->GetRootId()]->GetStrand();

                    if (gffTranscriptStrand == gffGeneStrand &&
                        (gffTranscriptStrand == gtfTranscriptStrand || gtfTranscriptStrand == '.')) {

                        goodTranscripts.push_back(transcript);
                    } 
                }

                // We only want genes with 1 or more good transcript
                if (goodTranscripts.size() >= 1) {                

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
           << " - # Transcripts: " << out.getTotalNbTranscripts() << " / " << in.getTotalNbTranscripts() << endl;
        
        report = ss.str();
    }
    
};
}
       
