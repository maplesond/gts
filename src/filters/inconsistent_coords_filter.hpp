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
    double cds;
    
public:
    
    InconsistentCoordsFilter(bool include, double cds) : TranscriptFilter(), include(include), cds(cds) {}
    
    ~InconsistentCoordsFilter() {}
    
    string getName() {
        return string("Inconsistent Transcript Filter");
    }
    
    string getDescription() {
        return string("Filters out CDSes that are inconsistent between transdecoder and full lengther");
    }
    
    double getCds() const {
        return cds;
    }

    void setCds(double cds) {
        this->cds = cds;
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
        
        uint32_t flnCompleteConsistent = 0;
        uint32_t flnNewConsistent = 0;
        BOOST_FOREACH(GFFIdMap::value_type i, uniqCds) {
            
            const string id = i.second->GetSeqId();
            shared_ptr<GFF> gff = i.second;
            
            bool consistent = false;
            // If 
            if (maps.uniqFlnCds.count(id)) {
                consistent = isTDCAndFLNConsistent(id, gff, maps.uniqFlnCds[id], false);

                if (consistent) {
                    flnCompleteConsistent++;
                }
            }
            else if (include && maps.uniqFlnNcCds.count(id)) {
                consistent = isTDCAndFLNConsistent(id, gff, maps.uniqFlnNcCds[id], false);
                
                if (consistent) {
                    flnNewConsistent++;
                }
            }
            
            if (consistent) {
                out.push_back(uniqueTranscripts[i.second->GetSeqId()]);
            }
        }
        
        stringstream ss;
        
        ss << " - Including consistent full lengther new coding hits: " << std::boolalpha << include << endl
           << " - Min required ratio of transdecoder to full lengther length: " << cds << endl
           << " - # Transdecoder CDSs with IDs matching those from input transcripts: " << uniqCds.size() << endl
           << " - # Transcripts with consistent and complete Full lengther annotations: " << flnCompleteConsistent << endl
           << " - # Transcripts with no complete full lengther annotation but with new coding entry (will be 0 if --include wasn't used): " << flnNewConsistent << endl;           
       
        report = ss.str();
    }

    bool isTDCAndFLNConsistent(const string& id, shared_ptr<GFF> tdc, shared_ptr<DBAnnot> fln, bool longCds) {
        int32_t deltaStart = std::abs(tdc->GetStart() - fln->GetStart());
        int32_t deltaEnd = std::abs(tdc->GetEnd() - fln->GetEnd());
        int32_t tdcLen = tdc->GetEnd() - tdc->GetStart();

        if (deltaStart <= gts::POS_THRESHOLD && deltaEnd <= gts::POS_THRESHOLD && 
                (!longCds || (longCds && tdcLen >= gts::LONG_CDS_LEN_THRESHOLD)) ) {

            const double seqFrac = (double)tdcLen / (double)fln->GetFastaLength();

            if (seqFrac >= cds) {
                return true;
            }
        }
        
        return false;
    }    
};

}
    
