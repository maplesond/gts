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

#include "transcript_filter.hpp"


namespace gts {

const int32_t POS_THRESHOLD = 10;
const int32_t LONG_CDS_LEN_THRESHOLD = 200;
    
class InconsistentCoordsFilter : public TranscriptFilter {
     
private:
    bool include;
    double cdsFrac;
    double cdnaFrac;
    
public:
    
    InconsistentCoordsFilter(bool include, double cdsFrac, double cdnaFrac) : 
        TranscriptFilter(), include(include), cdsFrac(cdsFrac), cdnaFrac(cdnaFrac) {}
    
    ~InconsistentCoordsFilter() {}
    
    string getName() {
        return string("Inconsistent Transcript Filter");
    }
    
    string getDescription() {
        return string("Filters out CDSes that are inconsistent between transdecoder and full lengther");
    }
    
    double getCdsFrac() const {
        return cdsFrac;
    }

    void setCdsFrac(double cdsFrac) {
        this->cdsFrac = cdsFrac;
    }
    
    double getCdnaFrac() const {
        return cdnaFrac;
    }

    void setCdnaFrac(double cdnaFrac) {
        this->cdnaFrac = cdnaFrac;
    }

    bool isInclude() const {
        return include;
    }

    void setInclude(bool include) {
        this->include = include;
    }

   
protected:    
    
    /**
     * Filters out transcripts with inconsistent coordinates.  Ensures consistency between:
     * 1 - input transcripts, 2 - transdecoder, 3 - full lengther
     * @param in Unique Transcripts that have both 5' and 3' UTRs
     * @param out 
     */
    void filterInternal(GFFList& in, Maps& maps, GFFList& out) {
        
        // Create a map of the input transcripts
        GFFIdMap uniqueTranscripts;            
        BOOST_FOREACH(shared_ptr<GFF> gff, in) {
           uniqueTranscripts[gff->GetRootId()] = gff;
        }
        
        // Create a map of transdecoder CDSes that are also found in the input transcripts
        GFFIdMap uniqCds;
        BOOST_FOREACH(GFFIdMap::value_type i, maps.transdecoderCdsGffMap) {        
            if (uniqueTranscripts.count(i.second->GetSeqId())) {
                uniqCds[i.first] = i.second;
            }
        }
        
        uint32_t matchingCdsFlnIds = 0;
        uint32_t flnCompleteConsistent = 0;
        uint32_t flnNewConsistent = 0;
        uint32_t similarTranscripts = 0;
        uint32_t notSimilarTranscripts = 0;
        
        BOOST_FOREACH(GFFIdMap::value_type i, uniqCds) {
            
            const string id = i.second->GetSeqId();
            shared_ptr<GFF> gff = i.second;
            int32_t tdcLen = gff->GetEnd() - gff->GetStart();
            
            bool consistent = false;
            bool longEnough = false;
            
            if (maps.uniqFlnCds.count(id)) {
                matchingCdsFlnIds++;
                consistent = isTDCAndFLNConsistent(gff, maps.uniqFlnCds[id], gts::POS_THRESHOLD, 2);

                if (consistent) {
                    flnCompleteConsistent++;
                    longEnough = isSeqLongEnough(tdcLen, maps.uniqFlnCds[id]->GetFastaLength(), cdsFrac);
                    if (longEnough) {
                        similarTranscripts++;
                    }
                }
            }            
            else if (include && maps.uniqFlnNcCds.count(id)) {
                matchingCdsFlnIds++;
                consistent = isTDCAndFLNConsistent(gff, maps.uniqFlnNcCds[id], gts::POS_THRESHOLD, gts::POS_THRESHOLD) &&
                        tdcLen >= gts::LONG_CDS_LEN_THRESHOLD;
                
                if (consistent) {
                    flnCompleteConsistent++;
                    longEnough = isSeqLongEnough(tdcLen, maps.uniqFlnNcCds[id]->GetFastaLength(), 0.5);
                    if (longEnough) {
                        notSimilarTranscripts++;
                    }
                }
            }
            
            if (consistent && longEnough) {
                out.push_back(uniqueTranscripts[i.second->GetSeqId()]);
            }
        }
        
        stringstream ss;
        
        ss << " - Including consistent full lengther new coding hits: " << std::boolalpha << include << endl
           << " - Min required ratio of transdecoder to full lengther length: " << cdsFrac << endl
           << " - # Transdecoder CDSs with IDs matching those from input transcripts: " << uniqCds.size() << endl
           << " - # Transdecoder CDSs with IDs matching Full Lengther transcripts: " << matchingCdsFlnIds << endl
           << " - # Transcripts with consistent transdecoder CDS and Full Lengther coordinates: " << flnCompleteConsistent << endl
           << " - # Consistent and long transcripts with similarity to Complete Full Lengther transcripts: " << similarTranscripts << endl           
           << " - # Consistent and long transcripts with no similarity (will be 0 if --include wasn't used): " << notSimilarTranscripts << endl;           
       
        report = ss.str();
    }

    bool isTDCAndFLNConsistent(shared_ptr<GFF> cds, shared_ptr<DBAnnot> fln, int32_t startThreshold, int32_t endThreshold) {
        
        int32_t deltaStart = std::abs(cds->GetStart() - fln->GetOrfStart());
        int32_t deltaEnd = std::abs(cds->GetEnd() - fln->GetOrfEnd());        

        return deltaStart <= startThreshold && deltaEnd <= endThreshold;
    }

    bool isSeqLongEnough(int32_t seqLen, int32_t fullTranscriptLen, double threshold) {        
        return ((double)seqLen / (double)fullTranscriptLen) >= threshold;
    }
};

}
    
