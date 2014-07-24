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
    
    void filterInternal(GFFList& in, Maps& maps, GFFList& out) {
        
        IdCounter mrna, utr5, utr3;
        uint32_t orfs = 0;
        
        // Count duplicate mRNAs, 5' UTRs and 3'UTRs
        BOOST_FOREACH(shared_ptr<GFF> gff, in) {

            switch(gff->GetType()) {
            case MRNA:
                mrna[gff->GetRootId()]++;
                orfs++;
                break;
            case UTR5:
                utr5[gff->GetRootId()]++;
                break;
            case UTR3:
                utr3[gff->GetRootId()]++;
                break;
            }
        }

        uint32_t uniqueTranscripts = 0;
        BOOST_FOREACH(IdCounter::value_type i, mrna) {
            
            const string id = i.first;
            
            if (i.second == 1) {                 
                shared_ptr<GFF> uniqMrna = maps.genomicGffMap[id];
                uniqueTranscripts++;
                if (utr3[id] > 0 && utr5[id] > 0) {
                    //uniqUtr[i.first] = uniqMrna;
                    out.push_back(uniqMrna);                
                }
            }
        }

        stringstream ss;
        
        ss << " - # ORFs: " << orfs << endl
           << " - # transcripts: " << mrna.size() << endl            
           << " - # transcripts with only one ORF: " << uniqueTranscripts << endl
           << " - # transcripts with only one ORF and at least one 5' and 3' UTR: " << out.size() << endl;
        
        report = ss.str();
    }
};

}