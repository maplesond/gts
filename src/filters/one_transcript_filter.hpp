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

#include <algorithm>
#include <sstream>
using std::stringstream;

#include <boost/unordered_map.hpp>
using boost::unordered_map;

#include "../gff.hpp"
#include "transcript_filter.hpp"
using gts::gff::GFFModel;
using gts::gff::GffType;


namespace gts {

typedef boost::unordered_map<string, uint32_t> IdCounter;
    
class OneTranscriptFilter : public TranscriptFilter {
      
public:

    OneTranscriptFilter() : TranscriptFilter() {}
    
    ~OneTranscriptFilter() {}
    
    string getName() {
        return string("One Transcript Per Gene Filter");
    }
    
    string getDescription() {
        return string("Selects the longest ORF transcript per gene");
    }
    
    
protected:    
    
    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        uint32_t inGenes = in.getGeneList()->size();


        // Identify multiple ORFs and remove them
        BOOST_FOREACH(GFFPtr gene, *(in.getGeneList())) {

            int32_t maxCdsLength = 0;
            GFFPtr longestTranscript;
            
            BOOST_FOREACH(GFFPtr transcript, *(gene->GetChildList())) {

                int32_t cdsLength = transcript->GetLengthOfAllTypes(CDS);
                
                if (cdsLength > maxCdsLength) {
                    maxCdsLength = cdsLength;
                    longestTranscript = transcript;
                }
            }
            
            if (longestTranscript) {

                // Copy gene without child info
                GFFPtr newGene = make_shared<GFF>(*gene);

                newGene->addChild(longestTranscript);
                
                out.addGene(newGene);
            }
            else {
                BOOST_THROW_EXCEPTION(TranscriptFilterException() << TranscriptFilterErrorInfo(string(
                            "No Transcript found for this gene: " + gene->GetId())));
            }
        }

        
        stringstream ss;
        
        ss << " - # Genes: " << out.getNbGenes() << " / " << in.getNbGenes() << endl
           << " - # Transcripts (mRNA): " << out.getTotalNbTranscripts() << " / " << in.getTotalNbTranscripts() << endl;
        
        report = ss.str();
    }
};

}

