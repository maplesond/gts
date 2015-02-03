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

class Cds2CDnaFilter : public TranscriptFilter {
     
private:
    double cdsFrac;
    
public:
    
    Cds2CDnaFilter(double cdsFrac) : 
        TranscriptFilter(), cdsFrac(cdsFrac) {}
    
    ~Cds2CDnaFilter() {}
    
    string getName() {
        return string("CDS 2 cDNA Length Ratio Transcript Filter");
    }
    
    string getDescription() {
        return string("Filters out transcripts whose CDS 2 cDNA length ratio is below threshold");
    }
    
    double getCdsFrac() const {
        return cdsFrac;
    }

    void setCdsFrac(double cdsFrac) {
        this->cdsFrac = cdsFrac;
    }
    
   
protected:    
    
    /**
     * Filters out transcripts with inconsistent coordinates.  Ensures consistency between:
     * cluster aligned transdecoder transcripts and full lengther
     * @param in Unique Transcripts that have both 5' and 3' UTRs
     * @param out 
     */
    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        uint32_t nbLongEnough = 0;
        
        BOOST_FOREACH(GFFPtr gene, *in.getGeneList()) {        
            
            GFFList goodTranscripts;
            
            BOOST_FOREACH(GFFPtr transcript, *gene->GetChildList()) {
                
                const string transcriptId = transcript->GetId();
                const string rootId = transcript->GetRootId();
                
                int32_t cdsLength = transcript->GetLengthOfAllTypes(CDS);
                int32_t cdnaLength = transcript->GetLengthOfAllTypes(EXON);

                if (cdnaLength == 0) {
                    BOOST_THROW_EXCEPTION(TranscriptFilterException() << TranscriptFilterErrorInfo(string(
                            "No Exons found for this transcript: " + transcript->GetId())));
                }

                if (cdsLength == 0) {
                    BOOST_THROW_EXCEPTION(TranscriptFilterException() << TranscriptFilterErrorInfo(string(
                            "No CDSes found for this transcript: " + transcript->GetId())));
                }
                
                // First check to see if CDS is long enough as a proportion of cDNA length
                bool longEnough = isSeqLongEnough(cdsLength, cdnaLength, cdsFrac);
                
                if (longEnough) {
                    goodTranscripts.push_back(transcript);
                    nbLongEnough++;
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
        
        ss << " - Min required ratio of transdecoder CDS to cDNA length: " << cdsFrac << endl
           << " -----------" << endl
           << " - # Genes: " << out.getNbGenes() << " / " << in.getNbGenes() << endl
           << " - # Transcripts (mRNA): " << out.getTotalNbTranscripts() << " / " << in.getTotalNbTranscripts() << endl;
        
        report = ss.str();
    }

    bool isSeqLongEnough(int32_t cdsLength, int32_t cdnaLength, double threshold) {        
        return ((double)cdsLength / (double)cdnaLength) >= threshold;
    }
};

}
    
