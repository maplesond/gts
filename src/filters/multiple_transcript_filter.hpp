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
    
    void filterInternal(GFFList& in, Maps& maps, GFFList& out) {
        
        TranscriptId2GeneIdMap geneList; 
        GFFIdMap geneMap;            
        
        BOOST_FOREACH(GFFIdMap::value_type i, maps.gtfMap) {
            geneList[i.second->GetRootTranscriptId()] = i.second->GetGeneId();
        }
        
        GeneId2LenMap geneCdnaLen;
        
        BOOST_FOREACH(shared_ptr<GFF> gff, in) {
            
            const string id = gff->GetRootId();
            
            const string geneId = geneList[id];
            
            if (geneCdnaLen.count(geneId) == 0 || 
                    (geneCdnaLen.count(geneId) > 0 && 
                        maps.allDistinctFlnCds[id]->GetFastaLength() > geneCdnaLen[geneId])) {
                
                   maps.allDistinctFlnCds[id]->GetFastaLength();
                   geneMap[geneId] = gff;
            }                
        }
        
        BOOST_FOREACH(GFFIdMap::value_type i, geneMap) {
            out.push_back(i.second);
        }
        
        stringstream ss;
        
        ss << " - # Genes: " << geneMap.size() << endl;
        
        report = ss.str();
        
        /*
my @lines = <CUFF>;
my %strand_cuff;
my %gene_list = ();
my %gene_list_exclude = ();
for (@lines){
        chomp;
    my @a = split(/\t/, $_);

    (my $id) = $a[8] =~ /.*; transcript_id "(.*?)"/;
    (my $gene) = $a[8] =~ /gene_id "(.*?)"/;
    $strand_cuff{$id} = $a[6];
    $gene_list{$id} = $gene;
    #print "$gene\n";

}

foreach my $b (@transcript) {
(my $gene) = ($gene_list{$b});
$gene_list_exclude{$b} = $gene;

}


foreach my $b (@merged_pass_cds) {
my @fln = ();
(my $gene) = ($gene_list{$b});
#print "$gene\n";

@fln = split(/\-/,$merged_fln_cds_ns_cds{$b});
#print "$fln[3]\n";

#$gene_list_exclude{$b} = $gene;

if ($fln[3] > $gene_cdna_length{$gene}) {
$gene_cdna_length{$gene} = $fln[3];
$gene_transcript{$gene} = $b;

}


}

foreach my $values (values %gene_transcript) {
    push(@transcripts_pass,$values);
}


my $size = @transcripts_pass;
print OUTPUTFILELOG "#Transcripts passing filter 3 (selecting 1 transcript per gene) $size\n";*/
    }
};
}