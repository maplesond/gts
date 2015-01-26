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
using std::stringstream;

#include <boost/unordered_map.hpp>
using boost::unordered_map;

#include "../gff.hpp"
#include "transcript_filter.hpp"
using gts::gff::GFFModel;
using gts::gff::GffType;


namespace gts {

typedef boost::unordered_map<string, uint32_t> IdCounter;
    
class UTRFilter : public TranscriptFilter {
      
public:

    UTRFilter() : TranscriptFilter() {}
    
    ~UTRFilter() {}
    
    string getName() {
        return string("UTR Filter");
    }
    
    string getDescription() {
        return string("Requires at least one 5' and 3' UTR for each transcript");
    }
    
    
protected:    
    
    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        uint32_t inGenes = in.getGeneList()->size();


        // Idenitfy multiple ORFs and remove them
        BOOST_FOREACH(GFFPtr gene, *(in.getGeneList())) {

            GFFList goodTranscripts;

            BOOST_FOREACH(GFFPtr transcript, *(gene->GetChildList())) {

                GFFListPtr utr5List = transcript->GetAllOfType(UTR5);
                GFFListPtr utr3List = transcript->GetAllOfType(UTR3);

                if (!utr5List->empty() && !utr3List->empty()) {
                    goodTranscripts.push_back(transcript);
                }
            }

            // We only want genes with at least 1 5' and 3' UTR
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
        
        ss << " - # Genes: " << out.getNbGenes() << " / " << in.getNbGenes() << endl
           << " - # Transcripts (mRNA): " << out.getTotalNbTranscripts() << " / " << in.getTotalNbTranscripts() << endl;
        
        report = ss.str();
    }
};

}

