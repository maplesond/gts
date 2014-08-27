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
    
class MultipleTranscriptFilter : public TranscriptFilter {
      
public:
    
    MultipleTranscriptFilter() : TranscriptFilter() {}
    
    ~MultipleTranscriptFilter() {}
    
    string getName() {
        return string("Multiple transcript Filter");
    }
    
    string getDescription() {
        return string("Selects the longest transcript per gene");
    }
    
    
protected:    
    
    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        BOOST_FOREACH(shared_ptr<GFF> gene, *(in.getGeneList())) {
            
            const string id = gene->GetRootId();            
            
            size_t nbTranscripts = gene->GetChildList()->size();
            
            if (nbTranscripts <= 0) {
                // Maybe should throw an error here... but for now we just ignore
                // Hopefully this should never happen
            }
            else if (nbTranscripts == 1) {
                // Just one transcript for this gene, so pass on.
                out.addGene(gene);
            }
            else {
                
                // More than one transcript so we need to pick one
                GFFPtr longest;
                int32_t longestLength = 0;
                
                // Go through each transcript and 
                BOOST_FOREACH(GFFPtr transcript, *(gene->GetChildList())) {
                    
                    // Sum the CDS entries for this transcript
                    int32_t cdsLength = 0;
                    BOOST_FOREACH(GFFPtr child, *(transcript->GetChildList())) {
                        if (child->GetType() == CDS) {
                            cdsLength += child->GetLength();
                        }
                    }
                    
                    if (cdsLength > longestLength) {
                        longestLength = cdsLength;
                        longest = transcript;                        
                    }
                }
                
                // Copy the gene GFF (without children) and then just tag on the 
                // existing transcript as a child
                shared_ptr<GFF> newGene = make_shared<GFF>(*gene);
                newGene->addChild(longest);
                out.addGene(newGene);                
            }                
        }
                
        stringstream ss;
        
        ss << " - # Genes: " << out.getNbGenes() << " / " << in.getNbGenes() << endl;
        ss << " - # Transcripts: " << out.getTotalNbTranscripts() << " / " << in.getTotalNbTranscripts() << endl;
        
        report = ss.str();  
    }
};
}