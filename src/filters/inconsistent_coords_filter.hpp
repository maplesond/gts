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

#include "gff.hpp"
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
     * cluster aligned transdecoder transcripts and full lengther
     * @param in Unique Transcripts that have both 5' and 3' UTRs
     * @param out 
     */
    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        uint32_t totalGenomicCDSes = 0;
        uint32_t totalMatchedCACDSes = 0;
        uint32_t matchingCdsFlnIds = 0;
        uint32_t flnCompleteConsistent = 0;
        uint32_t flnNewConsistent = 0;
        uint32_t similarTranscripts = 0;
        uint32_t notSimilarTranscripts = 0;
        
        BOOST_FOREACH(GFFPtr gene, *in.getGeneList()) {        
            
            GFFList goodTranscripts;
            
            BOOST_FOREACH(GFFPtr transcript, *gene->GetChildList()) {
                
                const string transcriptId = transcript->GetId();
                const string rootId = transcript->GetRootId();
                
                GFFListPtr genomicCDSes = transcript->GetAllOfType(gts::gff::CDS);
                
                if (genomicCDSes) {
                    
                    totalGenomicCDSes += genomicCDSes->size();
                    
                    // This loop just counts how may genomic CDSes we also find in the cluster aligned file
                    BOOST_FOREACH(GFFPtr genomicCDS, *genomicCDSes) {

                        const string cdsid = genomicCDS->GetId();

                        // Check 1: Can we find the genomic CDS in the cluster aligned GFF
                        if (maps.transdecoderCdsGffMap.count(cdsid)) {
                            totalMatchedCACDSes ++;
                        }
                    }
                    
                    BOOST_FOREACH(GFFPtr genomicCDS, *genomicCDSes) {

                        const string cdsid = genomicCDS->GetId();

                        // Check 1: Can we find the genomic CDS in the cluster aligned GFF
                        if (maps.transdecoderCdsGffMap.count(cdsid)) {

                            GFFPtr cacds = maps.transdecoderCdsGffMap[cdsid];

                            // Just a sanity check... make sure the transcript Ids in the 
                            // cluster aligned GFF and the genomic GFF for this CDS match
                            const string caTranscriptId = cacds->GetParentId();                        
                            if (!boost::equals(transcriptId, caTranscriptId)) {
                                BOOST_THROW_EXCEPTION(TranscriptFilterException() << TranscriptFilterErrorInfo(string(
                                    "Incompatible GFFs.  The genomic transcriptId and cluster aligned transcript Id for cluster CDS are not consistent: ") + cdsid));
                            }

                            // Get the cluster aligned CDS length
                            int32_t tdcLen = cacds->GetEnd() - cacds->GetStart();

                            bool consistent = false;
                            bool longEnough = false;

                            // Check 2: this cluster in full lengther to see if we can find something similar
                            if (maps.uniqFlnCds.count(rootId)) {
                                matchingCdsFlnIds++;
                                consistent = isTDCAndFLNConsistent(cacds, maps.uniqFlnCds[rootId], gts::POS_THRESHOLD, 2);

                                if (consistent) {
                                    flnCompleteConsistent++;
                                    longEnough = isSeqLongEnough(tdcLen, maps.uniqFlnCds[rootId]->GetFastaLength(), cdsFrac);
                                    if (longEnough) {
                                        similarTranscripts++;
                                    }
                                }
                            }            
                            else if (include && maps.uniqFlnNcCds.count(rootId)) {
                                matchingCdsFlnIds++;
                                consistent = isTDCAndFLNConsistent(cacds, maps.uniqFlnNcCds[rootId], gts::POS_THRESHOLD, gts::POS_THRESHOLD) &&
                                        tdcLen >= gts::LONG_CDS_LEN_THRESHOLD;

                                if (consistent) {
                                    flnCompleteConsistent++;
                                    longEnough = isSeqLongEnough(tdcLen, maps.uniqFlnNcCds[rootId]->GetFastaLength(), 0.5);
                                    if (longEnough) {
                                        notSimilarTranscripts++;
                                    }
                                }
                            }

                            if (consistent && longEnough) {
                                goodTranscripts.push_back(transcript);
                                break;  // Found a good enough CDS in this transcripts so skip out of the genomic CDS loop
                            }
                        }                    
                    }
                }
                else {
                    
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
           << " - Min required ratio of transdecoder to full lengther length: " << cdsFrac << endl
           << " - # Transdecoder cluster aligned CDSs with IDs also found in found in genomic GFF: " << totalMatchedCACDSes << " / " << totalGenomicCDSes << endl
           << " - # Transdecoder cluster aligned CDSs with IDs matching Full Lengther transcripts: " << matchingCdsFlnIds << " / " << totalMatchedCACDSes << endl
           << " - # Transcripts with consistent transdecoder CDS and Full Lengther coordinates: " << flnCompleteConsistent << " / " << matchingCdsFlnIds << endl
           << " - # Consistent and long transcripts with similarity to Complete Full Lengther transcripts: " << similarTranscripts << " / " << flnCompleteConsistent << endl           
           << " - # Consistent and long transcripts with no similarity (will be 0 if --include wasn't used): " << notSimilarTranscripts << " / " << flnCompleteConsistent << endl << endl           
           << " - # Genes: " << out.getNbGenes() << " / " << in.getNbGenes() << endl
           << " - # Transcripts (mRNA): " << out.getTotalNbTranscripts() << " / " << in.getTotalNbTranscripts() << endl;
        
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
    
