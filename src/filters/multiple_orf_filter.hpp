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
        return string("Keeps transcripts with a single ORF (in transdecoder terms, this means 1 transcript per gene/locus) and at least one 5' and 3' UTR");
    }
    
    
protected:

    bool isExonListSame(GFFList& exonList1, GFFList& exonList2) {

        if (exonList1.size() != exonList2.size()) {

            return false;
        }

        // Assume exon list is in coord order
        for(size_t i = 0; i < exonList1.size(); i++) {
            GFFPtr e1 = exonList1[i];
            GFFPtr e2 = exonList2[i];

            if (e1->GetStart() != e2->GetStart() ||
                    e1->GetEnd() != e2->GetEnd()) {
                return false;
            }
        }

        return true;
    }

    /**
    * Adds all transcripts that are share ORFs to the provided set.
    */
    void addMORFs(GFFList& transcripts, unordered_set<GFFPtr>& morfs) {

        for(size_t i = 0; i < transcripts.size(); i++) {

            GFFPtr ti = transcripts[i];

            vector<GFFPtr> multiOrfTranscripts;

            for (size_t j = 0; j < transcripts.size(); j++) {

                if (i != j) {

                    GFFPtr tj = transcripts[j];

                    if (ti->GetStart() == tj->GetStart() &&
                            ti->GetEnd() == tj->GetEnd()) {

                        cout << "Same transcript: " << ti->GetId() << endl;

                        if (isExonListSame(*(ti->GetAllOfType(EXON)), *(tj->GetAllOfType(EXON)))) {
                            morfs.insert(ti);
                            morfs.insert(tj);
                        }
                    }
                }
            }
        }
    }

    void filterInternal(GFFModel& in, Maps& maps, GFFModel& out) {
        
        // Idenitfy multiple ORFs and remove them
        BOOST_FOREACH(GFFPtr gene, *(in.getGeneList())) {

            // Gene should only have mRNAs (transcripts) as children
            uint32_t geneTranscripts = gene->GetNbChildren();

            GFFList goodTranscripts;
            unordered_set<GFFPtr> morfs;

            // Identify all those transcripts which share open reading frames
            addMORFs(*(gene->GetChildList()), morfs);

            if (morfs.size() > 0) {
                cout << morfs.size() << endl;
            }

            BOOST_FOREACH(GFFPtr transcript, *(gene->GetChildList())) {

                // Only add the transcript if it doesn't share an open reading frame with another transcript in this gene
                if (morfs.find(transcript) == morfs.end()) {
                    goodTranscripts.push_back(transcript);
                }
            }

            // We only want genes with 1 transcript and at least 1 5' and 3' UTR
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