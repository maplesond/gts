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
    
class StrandFilter : public TranscriptFilter {
      
public:
    
    StrandFilter() : TranscriptFilter() {}
    
    ~StrandFilter() {}
    
    string getName() {
        return string("Strand Filter");
    }
    
    string getDescription() {
        return string("Filters out transcripts which have an inconsistent strand");
    }
    
    
protected:    
    
    void filterInternal(GFFList& in, Maps& maps, GFFList& out) {
    /*    
        my %strand_decoder;
my @transcripts_pass_all;


for (@gff_pass){
    my @a = split;
    my $id;
    my $g;

    if ($a[2] =~ m/mRNA/i) {
    ($id) = $a[8] =~ /ID=(.*?)\|/;
    $strand_decoder{$id} = $a[6];
    }
}

foreach my $key (keys %strand_decoder) {

if ($strand_decoder{$key} eq $strand_cuff{$key} || $strand_cuff{$key} eq '.') {

push(@transcripts_pass_all,$key);

}
}

for (@lines_gff){
    my @a = split;
    my $id;
    if ($a[2] =~ m/gene|mRNA|exon|five_prime_UTR|three_prime_UTR/i) {
    ($id) = $a[8] =~ /ID=(.*?)\|/;
    }
    if ($a[2] =~ m/CDS/i) {
    ($id) = $a[8] =~ /ID=cds\.(.*?)\|/;
    }

    foreach my $b (@transcripts_pass_all) {
    if ($b eq $id) {
    my $c = join("\t", @a);
    push (@gff_pass_all,$c);

    }
    }
}

my $size = @transcripts_pass_all;
print OUTPUTFILELOG "#Transcripts passing filter 4 (final strand check) $size\n";
*/
    }
    
};
}
       
