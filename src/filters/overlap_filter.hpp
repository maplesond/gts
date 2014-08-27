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
#include <stdlib.h>

#include "transcript_filter.hpp"

using std::stringstream;        

namespace gts {

typedef boost::unordered_map<string, uint32_t> IdCounter;
    
class OverlapFilter : public TranscriptFilter {
   
private:
    uint32_t windowSize;
    
public:
    
    OverlapFilter(uint32_t windowSize) : TranscriptFilter() {
        this->windowSize = windowSize;
    }
    
    ~OverlapFilter() {}
    
    string getName() {
        return string("Overlap Filter");
    }
    
    string getDescription() {
        return string("Filters out genomic transcripts which overlap with each other or are within a given window.");
    }
    
    
protected:    
    
    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        GFFPtr lastOkGene;
        
        bool first = true;
        bool ok = true;

        stringstream ss;
        
        ss << " - Window Size: " << this->windowSize << endl;
        
        BOOST_FOREACH(shared_ptr<GFF> gene1, *(in.getGeneList())) {
            
            uint16_t overlapCount = 0;
            
            BOOST_FOREACH(shared_ptr<GFF> gene2, *(in.getGeneList())) {
                
                if (gene1 != gene2) {
                
                    if (boost::equals(gene1->GetSeqId(), gene2->GetSeqId())) {

                        int32_t startDiff = abs(gene1->GetStart() - gene2->GetStart());
                        int32_t endDiff = abs(gene1->GetEnd() - gene2->GetEnd());
                        int32_t rangeDiff = gene1->GetStart() - gene2->GetEnd();
                        if (startDiff <= this->windowSize || endDiff <= this->windowSize || rangeDiff <= this->windowSize) {
                            overlapCount++;
                        }
                    }                    
                }
            }

            if (overlapCount == 0) {
                out.addGene(gene1);
            }
        }
                
        ss << " - # Genes: " << out.getNbGenes() << " / " << in.getNbGenes() << endl
           << " - # Transcripts: " << out.getTotalNbTranscripts() << " / " << in.getNbGenes() << endl;
        
        report = ss.str();
    }
    
};
}
       
