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

typedef boost::unordered_map<string, uint32_t> GeneId2LenMap;
typedef boost::unordered_map<string, string> TranscriptId2GeneIdMap;

    
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
        
        TranscriptId2GeneIdMap geneList;
        GeneId2LenMap geneCdnaLen;                
        GFFIdMap geneMap;
        
        BOOST_FOREACH(shared_ptr<GFF> gff, *(in.getGeneList())) {
            
            const string id = gff->GetRootId();            
            const string geneId = geneList[id];
            
            if (maps.allDistinctFlnCds.count(id) > 0 &&
                (geneCdnaLen.count(geneId) == 0 || 
                    (geneCdnaLen.count(geneId) > 0 && 
                    maps.allDistinctFlnCds[id]->GetFastaLength() > geneCdnaLen[geneId]))) {
                
                geneCdnaLen[geneId]
                maps.allDistinctFlnCds[id]->GetFastaLength();
                                
                if ()
                shared_ptr<GFF> newGene = shared_ptr<GFF>(new GFF(*gff));
                
                newGene->addChild(transcript);
                
                out.addGene(newGene);
            }                
        }
        
        BOOST_FOREACH(GFFIdMap::value_type i, geneMap) {
            
        }
        
        
        stringstream ss;
        
        ss << " - # Genes: " << out.getNbGenes() << " / " << in.getNbGenes() << endl;
        ss << " - # Transcripts: " << out.getTotalNbTranscripts() << " / " << in.getTotalNbTranscripts() << endl;
        
        report = ss.str();  
    }
};
}