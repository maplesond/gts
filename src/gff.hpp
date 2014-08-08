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

#include <string>
#include <vector>
using std::string;
using std::vector;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
using boost::lexical_cast;
using boost::shared_ptr;

namespace gts {
namespace gff {

typedef boost::error_info<struct GFFError,string> GFFErrorInfo;
struct GFFException: virtual boost::exception, virtual std::exception { };

enum FileFormat {
    GFF2,
    GFF3,
    GTF
};

static FileFormat fileFormatFromString(string& s) {
    
    if (s.compare("GFF2") == 0) {
        return GFF2;
    }
    else if (s.compare("GFF3") == 0) {
        return GFF3;
    }
    else if (s.compare("GTF") == 0) {
        return GTF;
    }
    else {
        BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
            "Could not recognise GFF style file format: ") + s));
    }
}

enum GffType {
    GENE,
    MRNA,
    UTR5,
    UTR3,
    CDS,
    TRANSCRIPT,
    EXON,
    OTHER
};

static GffType gffTypeFromString(string& s) {
    
    if (boost::iequals(s, "gene")) {
        return GENE;
    }
    else if (boost::iequals(s, "mRNA")) {
        return MRNA;
    }
    else if (boost::iequals(s, "five_prime_utr")) {
        return UTR5;
    }
    else if (boost::iequals(s, "three_prime_utr")) {
        return UTR3;
    }
    else if (boost::iequals(s, "cds")) {
        return CDS;
    }
    else if (boost::iequals(s, "transcript")) {
        return TRANSCRIPT;
    }
    else if (boost::iequals(s, "exon")) {
        return EXON;
    }
    else {
        return OTHER;
    }
}

static string gffTypeToString(GffType type) {
    
    switch(type) {
        case GENE:
            return "gene";
        case MRNA:
            return "mRNA";
        case UTR5:
            return "five_prime_utr";
        case UTR3:
            return "three_prime_utr";
        case CDS:
            return "cds";
        case TRANSCRIPT:
            return "transcript";
        case EXON:
            return "exon";
        case OTHER:
            return ".";
    }    
}

class GFF {
    
private:
    
    // Enum indicating which format this file is: (GFF2, GFF3, GTF)
    FileFormat fileFormat;
    
    // GFF2, GFF3 and GTF features
    string seqId;
    string source;    
    GffType type;
    int32_t start;
    int32_t end;
    double score;
    char strand;
    int8_t phase;
    
    // GFF3 attributes
    string id;
    string cdsid;
    string name;
    string alias;
    string parent;
    string target;
    string gap;
    bool circular;
    
    // GTF attributes
    string geneId;
    string transcriptId;
    
    // Cufflinks GTF attributes
    uint16_t exonNumber;
    double fpkm;
    double frac;
    double confLo;
    double confHigh;    
    double coverage;
    
public:

    GFF(FileFormat fileFormat) : 
        fileFormat(fileFormat),
        seqId(""), 
        source(""), 
        type(OTHER), 
        start(0), end(0),
        score(-1.0),
        strand('.'),
        phase(-1) {
    }

    virtual ~GFF() {}
    
    int8_t GetPhase() const {
        return phase;
    }

    void SetPhase(int8_t phase) {
        this->phase = phase;
    }

    double GetScore() const {
        return score;
    }

    void SetScore(double score) {
        this->score = score;
    }
    
    char GetStrand() const {
        return strand;
    }

    void SetStrand(char strand) {
        this->strand = strand;
    }

    string GetSeqId() const {
        return seqId;
    }

    void SetSeqId(string seqId) {
        this->seqId = seqId;
    }

    string GetSource() const {
        return source;
    }

    void SetSource(string source) {
        this->source = source;
    }

    
    double GetCoverage() const {
        return coverage;
    }

    void SetCoverage(double coverage) {
        this->coverage = coverage;
    }

    string GetId() const {
        return id;
    }

    void SetId(string id) {
        this->id = id;
    }
    
    string GetRootId() const {
        
        vector<string> idElements;
        boost::split( idElements, id, boost::is_any_of("|"), boost::token_compress_on );
        size_t pos = idElements[0].find("cds");
        return pos != std::string::npos ? idElements[0].substr(pos+4) : idElements[0];        
    }

    string GetName() const {
        return name;
    }

    void SetName(string name) {
        this->name = name;
    }

    string GetParent() const {
        return parent;
    }

    void SetParent(string parent) {
        this->parent = parent;
    }

    GffType GetType() const {
        return type;
    }

    void SetType(GffType type) {
        this->type = type;
    }
    
    int32_t GetEnd() const {
        return end;
    }

    void SetEnd(int32_t end) {
        this->end = end;
    }
    
    int32_t GetLength() const {
        return std::abs(end - start);
    }

    int32_t GetStart() const {
        return start;
    }

    void SetStart(int32_t start) {
        this->start = start;
    }
    
    string GetAlias() const {
        return alias;
    }

    void SetAlias(string alias) {
        this->alias = alias;
    }

    bool IsCircular() const {
        return circular;
    }

    void SetCircular(bool circular) {
        this->circular = circular;
    }

    FileFormat GetFileFormat() const {
        return fileFormat;
    }

    void SetFileFormat(FileFormat fileFormat) {
        this->fileFormat = fileFormat;
    }

    string GetGap() const {
        return gap;
    }

    void SetGap(string gap) {
        this->gap = gap;
    }

    string GetGeneId() const {
        return geneId;
    }

    void SetGeneId(string geneId) {
        this->geneId = geneId;
    }

    string GetTarget() const {
        return target;
    }

    void SetTarget(string target) {
        this->target = target;
    }

    string GetTranscriptId() const {
        return transcriptId;
    }
    
    string GetRootTranscriptId() const {
        vector<string> idElements;
        boost::split( idElements, transcriptId, boost::is_any_of("|"), boost::token_compress_on );
        
        return idElements.size() == 1 ? idElements[0] :
                    idElements.size() == 2 ? idElements[1] :
                        "";
    }

    void SetTranscriptId(string transcriptId) {
        this->transcriptId = transcriptId;
    }
    
    double GetConfHigh() const {
        return confHigh;
    }

    void SetConfHigh(double confHigh) {
        this->confHigh = confHigh;
    }

    double GetConfLo() const {
        return confLo;
    }

    void SetConfLo(double confLo) {
        this->confLo = confLo;
    }

    uint16_t GetExonNumber() const {
        return exonNumber;
    }

    void SetExonNumber(uint16_t exonNumber) {
        this->exonNumber = exonNumber;
    }

    double GetFpkm() const {
        return fpkm;
    }

    void SetFpkm(double fpkm) {
        this->fpkm = fpkm;
    }

    double GetFrac() const {
        return frac;
    }

    void SetFrac(double frac) {
        this->frac = frac;
    }


    void writeGFF3Attribs(std::ostream& out) {
        
        vector<string> elems;
        
        elems.push_back(string("ID=") + id);
        
        if (!name.empty()) {
            elems.push_back(string("Name=") + name);
        }
        
        if (!alias.empty()) {
            elems.push_back(string("Alias=") + alias);
        }
        
        if (!parent.empty()) {
            elems.push_back(string("Parent=") + parent);
        }
        
        if (!target.empty()) {
            elems.push_back(string("Target=") + target);
        }
        
        if (!gap.empty()) {
            elems.push_back(string("Gap=") + gap);
        }

        out << boost::algorithm::join(elems, ";");        
    }
    
    void writeGTFAttribs(std::ostream& out) {
        
        vector<string> elems;
        
        elems.push_back(string("gene_id \"") + geneId + "\"");
        elems.push_back(string("transcript_id \"") + transcriptId + "\"");
        
        if (exonNumber >= 0) {
            elems.push_back(string("exon_number \"") + boost::lexical_cast<string>(exonNumber) + "\"");
        }
        
        if (fpkm >= 0) {
            elems.push_back(string("FPKM \"") + boost::lexical_cast<string>(fpkm) + "\"");
        }
        
        if (frac >= 0.0) {
            elems.push_back(string("frac \"") + boost::lexical_cast<string>(frac) + "\"");
        }
        
        if (confLo >= 0.0) {
            elems.push_back(string("conf_lo \"") + boost::lexical_cast<string>(confLo) + "\"");
        }
        
        if (confHigh >= 0.0) {
            elems.push_back(string("conf_hi \"") + boost::lexical_cast<string>(confHigh) + "\"");
        }
        
        if (coverage >= 0.0) {
            elems.push_back(string("cov \"") + boost::lexical_cast<string>(coverage) + "\"");
        }

        out << boost::algorithm::join(elems, ";");
    }

    void write(std::ostream& out) {
        
        std::stringstream ss;
        ss << "ID=" << id << ";";
        
        if (!parent.empty()) {
            ss << "Parent=" << parent << ";";
        }
        
        out << seqId << "\t" 
            << source << "\t"
            << gffTypeToString(type) << "\t"
            << start << "\t"
            << end << "\t"
            << (score == -1.0 ? "." : boost::lexical_cast<std::string>(score)) << "\t"
            << strand << "\t"
            << (phase == -1 ? "." : boost::lexical_cast<std::string>(phase)) << "\t";
            
        if (fileFormat == GFF3) {
            writeGFF3Attribs(out);
        }
        else if (fileFormat == GTF) {
            writeGTFAttribs(out);
        }
        
        out << endl;
    }
    
    
    static shared_ptr<GFF> parse(FileFormat fileFormat, const string& line) {
        vector<string> parts;
        boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_on );

        if (parts.size() != 9) {
            BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                "Could not parse GFF line due to incorrect number of columns. Expected 9 columns: ") + line));
        }
        
        shared_ptr<GFF> gff = shared_ptr<GFF>(new GFF(fileFormat));
        
        gff->SetSeqId(parts[0]);
        gff->SetSource(parts[1]);
        gff->SetType(gffTypeFromString(parts[2]));
        gff->SetStart(lexical_cast<int32_t>(parts[3]));
        gff->SetEnd(lexical_cast<int32_t>(parts[4]));
        gff->SetScore(parts[5][0] == '.' ? -1.0 : lexical_cast<double>(parts[5]));
        gff->SetStrand(parts[6][0]);
        gff->SetPhase(parts[7][0] == '.' ? -1 : lexical_cast<int8_t>(parts[7]));
        
        vector<string> attrParts;
        boost::split( attrParts, parts[8], boost::is_any_of(";"), boost::token_compress_on );
        
        BOOST_FOREACH(string attr, attrParts) {
            
            boost::trim(attr);
            
            if (!attr.empty()) {
                vector<string> attrElements;
                
                if (fileFormat == GFF3) {
                    boost::split( attrElements, attr, boost::is_any_of("="), boost::token_compress_on );

                    string key = attrElements[0];
                    string val = attrElements[1];
                    if (boost::iequals(key, "ID")) {
                        gff->SetId(val);
                    }
                    else if (boost::iequals(key, "Name")) {
                        gff->SetName(val);
                    }
                    else if (boost::iequals(key, "Parent")) {
                        gff->SetParent(val);
                    }
                    else if (boost::iequals(key, "Alias")) {
                        gff->SetAlias(val);
                    }
                    else if (boost::iequals(key, "Target")) {
                        gff->SetTarget(val);
                    }
                    else if (boost::iequals(key, "Gap")) {
                        gff->SetGap(val);
                    }
                }
                else if(fileFormat == GTF || fileFormat == GFF2) {
                    boost::split( attrElements, attr, boost::is_any_of(" "), boost::token_compress_on );

                    string key = attrElements[0];
                    string val = attrElements[1].substr(1, attrElements[1].size()-2);
                    
                    if (boost::iequals(key, "gene_id")) {
                        gff->SetGeneId(val);
                    }
                    else if (boost::iequals(key, "transcript_id")) {
                        gff->SetTranscriptId(val);
                    }
                    else if (boost::iequals(key, "exon_number")) {
                        gff->SetExonNumber(lexical_cast<uint16_t>(val));
                    }
                    else if (boost::iequals(key, "FPKM")) {
                        gff->SetFpkm(lexical_cast<double>(val));
                    }
                    else if (boost::iequals(key, "frac")) {
                        gff->SetFrac(lexical_cast<double>(val));
                    }
                    else if (boost::iequals(key, "conf_lo")) {
                        gff->SetConfLo(lexical_cast<double>(val));
                    }
                    else if (boost::iequals(key, "conf_hi")) {
                        gff->SetConfHigh(lexical_cast<double>(val));
                    }
                    else if (boost::iequals(key, "coverage")) {
                        gff->SetCoverage(lexical_cast<double>(val));
                    }                 
                }
                else {
                    BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                        "Could not recognise GFF style file format.")));
                }        
            }            
        }
        
        return gff;
    }
    
    
    static void load(FileFormat fileFormat, const string& path, std::vector< boost::shared_ptr<GFF> >& gffs) {
    
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        cout << " - Loading GFF: " << path << endl;
        
        std::ifstream file(path.c_str());
        std::string line; 
        while (std::getline(file, line)) {            
            boost::trim(line);
            if (!line.empty()) {
                gffs.push_back(parse(fileFormat, line));
            }
        }
        file.close();
        
        cout << " - Found " << gffs.size() << " GFF records." << endl;
    }
    
    static void save(const string& path, std::vector< boost::shared_ptr<GFF> >& gffs) {
        
        auto_cpu_timer timer(1, "- Wall time taken: %ws\n\n");
        cout << "Saving GFF: " << path << endl;
        
        std::ofstream file(path.c_str());
        BOOST_FOREACH(shared_ptr<GFF> gff, gffs) {
            gff->write(file);
        }
        file.close();
    }
};

struct GFFOrdering {
    inline bool operator ()(const shared_ptr<GFF>& a, const shared_ptr<GFF>& b) {
        
        int seqId = a->GetSeqId().compare(b->GetSeqId());
        if (seqId != 0) {
            return seqId < 0;            
        }
        else {
            int sDiff = a->GetStart() - b->GetStart();
            if (sDiff != 0) {
                return a->GetStart() < b->GetStart();
            }
            else {
                return a->GetEnd() < b->GetEnd();
            }
        }        
    }
};

}
}