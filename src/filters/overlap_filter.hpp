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

class OverlapFilter : public TranscriptFilter {
   
private:
    uint32_t windowSize;
    GFFModelPtr fullModel;
    
public:
    
    OverlapFilter(uint32_t windowSize, GFFModelPtr fullModel) : TranscriptFilter() {
        this->windowSize = windowSize;
        this->fullModel = fullModel;
    }
    
    ~OverlapFilter() {}
    
    string getName() {
        return string("Overlap Filter");
    }
    
    string getDescription() {
        return string("Filters out genomic transcripts which overlap with each other or are within a given window.  Checks the genes that have passed all previous filters against the genes present in the original model.");
    }
    
    
protected:    
    
    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        GFFPtr lastOkGene;
        
        bool first = true;
        bool ok = true;

        stringstream ss;
        
        ss << " - Window Size: " << this->windowSize << endl;
        
        BOOST_FOREACH(GFFPtr gene1, *(in.getGeneList())) {
            
            uint16_t overlapCount = 0;
            
            BOOST_FOREACH(GFFPtr gene2, *(fullModel->getGeneList())) {
                
                if (gene1 != gene2 && 
                    !boost::equals(gene1->GetId(), gene2->GetId()) &&
                    boost::equals(gene1->GetSeqId(), gene2->GetSeqId())) {
                
                    // This is probably overkill but let's test everything just to be sure!
                    int32_t startDiff = abs(gene1->GetStart() - gene2->GetStart());
                    int32_t endDiff = abs(gene1->GetEnd() - gene2->GetEnd());
                    int32_t rangeDiff1 = abs(gene1->GetStart() - gene2->GetEnd());
                    int32_t rangeDiff2 = abs(gene1->GetEnd() - gene2->GetStart());
                    bool overlapping = ((gene1->GetStart() < gene2->GetStart() && gene1->GetEnd() > gene2->GetStart()) || 
                                        (gene2->GetStart() < gene1->GetStart() && gene2->GetEnd() > gene1->GetStart()));

                    if (overlapping ||
                        startDiff <= this->windowSize || 
                        endDiff <= this->windowSize || 
                        rangeDiff1 <= this->windowSize ||
                        rangeDiff2 <= this->windowSize ) {

                        overlapCount++;
                    }                                      
                }
            }

            if (overlapCount == 0) {
                out.addGene(gene1);
            }
        }
                
        ss << " - Checking " << in.getNbGenes() << " passed genes against the " << fullModel->getNbGenes() << " genes present in original model" << endl
           << " - # Genes: " << out.getNbGenes() << " / " << in.getNbGenes() << endl
           << " - # Transcripts: " << out.getTotalNbTranscripts() << " / " << in.getNbGenes() << endl;
        
        report = ss.str();
    }
    
};
}
       
