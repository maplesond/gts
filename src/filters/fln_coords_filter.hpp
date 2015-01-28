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

#include "../gff.hpp"
#include "transcript_filter.hpp"

using gts::gff::GffType;


namespace gts {

const int32_t POS_THRESHOLD = 10;
const int32_t LONG_CDS_LEN_THRESHOLD = 200;
    
class FlnCoordsFilter : public TranscriptFilter {
     
private:
    bool include;
    
public:
    
    FlnCoordsFilter(bool include) : 
        TranscriptFilter(), include(include) {}
    
    ~FlnCoordsFilter() {}
    
    string getName() {
        return string("Inconsistent Transcript Filter");
    }
    
    string getDescription() {
        return string("Filters out transcripts whose CDS is inconsistent with full lengther coordinates");
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
     * cluster aligned transdecoder transcripts and full lengther
     * @param in Unique Transcripts that have both 5' and 3' UTRs
     * @param out 
     */
    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        uint32_t nbConsistent = 0;
        uint32_t nbFlnDbAnnotConsistent = 0;
        uint32_t nbFlnNewCodingConsistent = 0;
        
        BOOST_FOREACH(GFFPtr gene, *in.getGeneList()) {        
            
            GFFList goodTranscripts;
            
            BOOST_FOREACH(GFFPtr transcript, *gene->GetChildList()) {
                
                const string transcriptId = transcript->GetId();
                const string rootId = transcript->GetRootId();
                
                int32_t cdsLength = transcript->GetLengthOfAllTypes(CDS);

                if (cdsLength == 0) {
                    BOOST_THROW_EXCEPTION(TranscriptFilterException() << TranscriptFilterErrorInfo(string(
                            "No CDSes found for this transcript: " + transcript->GetId())));
                }
                
                int32_t cdsStartOffset = getCdsStartOffset(transcript);
                int32_t cdsEndOffset = cdsStartOffset + cdsLength - 1;
                
                bool consistent = false;
                    
                // Check: this transcript's CDS is consistent with full lengther (dbannotated)
                if (maps.uniqFlnCds.count(rootId)) {

                    if (consistent = isTDCAndFLNConsistent(cdsStartOffset, cdsEndOffset, maps.uniqFlnCds[rootId], gts::POS_THRESHOLD, 2)) {
                        nbFlnDbAnnotConsistent++;
                    }
                }
                // Check: If requested by the user, and we didn't find anything is fln dbannotated, check if this transcript's 
                // CDS is consistent with full lengther (new_coding)
                else if (include && maps.uniqFlnNcCds.count(rootId)) {
                    if (consistent = (isTDCAndFLNConsistent(cdsStartOffset, cdsEndOffset, maps.uniqFlnNcCds[rootId], gts::POS_THRESHOLD, gts::POS_THRESHOLD) &&
                            cdsLength >= gts::LONG_CDS_LEN_THRESHOLD)) {
                        nbFlnNewCodingConsistent++;
                    }
                }

                if (consistent) {
                    goodTranscripts.push_back(transcript);
                    nbConsistent++;
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
        
        stringstream ss;
        
        ss << " - Including consistent full lengther new coding hits: " << std::boolalpha << include << endl
           << " ------------" << endl
           << " - # Transcripts with CDS consistent with Full Lengther coordinates: " << nbConsistent << endl
           << "   - # From DBAnnotated file: " << nbFlnDbAnnotConsistent << endl           
           << "   - # From NewCoding file (will be 0 if not requested): " << nbFlnNewCodingConsistent << endl
           << " -----------" << endl
           << " - # Genes: " << out.getNbGenes() << " / " << in.getNbGenes() << endl
           << " - # Transcripts (mRNA): " << out.getTotalNbTranscripts() << " / " << in.getTotalNbTranscripts() << endl;
        
        report = ss.str();
    }

    bool isTDCAndFLNConsistent(int32_t cdsStartOffset, int32_t cdsEnd, shared_ptr<DBAnnot> fln, int32_t startThreshold, int32_t endThreshold) {
        
        int32_t deltaStart = std::abs(cdsStartOffset - fln->GetOrfStart());
        int32_t deltaEnd = std::abs(cdsEnd - fln->GetOrfEnd());

        //cout << "TranscriptID: " << fln->GetId() << endl;
        //cout << "Start\t" << cdsStartOffset << "\t" << fln->GetOrfStart() << "\t" << deltaStart << endl;
        //cout << "End\t" << cdsEnd << "\t" << fln->GetOrfEnd() << "\t" << deltaEnd << endl;

        return deltaStart <= startThreshold && deltaEnd <= endThreshold;
    }
    
    int32_t getCdsStartOffset(GFFPtr transcript) {
        
        int32_t offset = 0;
        
        GFFListPtr cdses = transcript->GetAllOfType(CDS);
        GFFListPtr exons = transcript->GetAllOfType(EXON);
        
        if (!cdses) {
           BOOST_THROW_EXCEPTION(TranscriptFilterException() << TranscriptFilterErrorInfo(string(
                    "No CDS found for this transcript: " + transcript->GetId())));
        }
                
        if (!exons) {
           BOOST_THROW_EXCEPTION(TranscriptFilterException() << TranscriptFilterErrorInfo(string(
                    "No Exons found for this transcript: " + transcript->GetId())));
        }
         
        GFFPtr cds = cdses->at(0);        
        
        int32_t cDnaLen = 0;
        for(size_t i = 0; i < exons->size(); i++) {
            
            GFFPtr exon = exons->at(i);
            
            if (exon->GetStart() < cds->GetStart()) {
                cDnaLen += exon->GetLength();
            }
            else {
                break;
            }
        } 
        
        // Add one to get into the 1-based coord system of full lengther
        return cDnaLen - cds->GetLength() + 1;
    }
};

}
    
