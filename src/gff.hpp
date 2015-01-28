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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::ostream;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/enable_shared_from_this.hpp>
using boost::lexical_cast;
using boost::timer::auto_cpu_timer;
using boost::make_shared;
using boost::shared_ptr;
using boost::unordered_map;

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
    MIRNA,
    PROTEIN,
    UTR5,
    UTR3,
    CDS,
    TRANSCRIPT,
    EXON,
    TSS,
    TTS,
    OTHER,
    ANY
};

static GffType gffTypeFromString(string& s) {
    
    if (boost::iequals(s, "gene")) {
        return GENE;
    }
    else if (boost::iequals(s, "mRNA")) {
        return MRNA;
    }
    else if (boost::iequals(s, "miRNA")) {
        return MIRNA;
    }
    else if (boost::iequals(s, "protein")) {
        return PROTEIN;
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
    else if (boost::iequals(s, "tss")) {
        return TSS;
    }
    else if (boost::iequals(s, "tts")) {
        return TTS;
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
        case MIRNA:
            return "miRNA";
        case PROTEIN:
            return "protein";
        case UTR5:
            return "five_prime_utr";
        case UTR3:
            return "three_prime_utr";
        case CDS:
            return "CDS";
        case TRANSCRIPT:
            return "transcript";
        case EXON:
            return "exon";
        case TSS:
            return "tss";
        case TTS:
            return "tts";
        case OTHER:
            return ".";
    }    
}

class GFF;
class GFFModel;

typedef boost::shared_ptr<gts::gff::GFF> GFFPtr;
typedef boost::shared_ptr<gts::gff::GFFModel> GFFModelPtr;

typedef std::vector<GFFPtr> GFFList;
typedef boost::unordered_map<string, GFFPtr> GFFIdMap;

typedef boost::shared_ptr<GFFList> GFFListPtr;
typedef boost::shared_ptr<GFFIdMap> GFFIdMapPtr;

class GFF : public boost::enable_shared_from_this<GFF> {
    
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
    string note;
    string parentId;
    string target;
    string gap;
    bool circular;
    string derivesFrom;
    string index;
    
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
    
    // Children
    GFFPtr parent;
    GFFIdMapPtr childMap;
    GFFListPtr childList;
    
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
            
        childMap = make_shared<GFFIdMap>();
        childList = make_shared<GFFList>();
        
    }
        
    GFF(const GFF& gff) : 
        fileFormat(gff.fileFormat),
        seqId(gff.seqId), 
        source(gff.source), 
        type(gff.type), 
        start(gff.start), end(gff.end),
        score(gff.score),
        strand(gff.strand),
        phase(gff.phase) {
        
        id = gff.id;
        cdsid = gff.cdsid;
        name = gff.name;
        alias = gff.alias;
        note = gff.note;
        parentId = gff.parentId;
        target = gff.target;
        gap = gff.gap;
        circular = gff.circular;
        derivesFrom = gff.derivesFrom;
        index = gff.index;
        
        childMap = make_shared<GFFIdMap>();
        childList = make_shared<GFFList>();
        
    }

    virtual ~GFF() {

        removeLinks();
    }
    
    void removeLinks() {
        
        if (childList && !childList->empty()) {
            BOOST_FOREACH(GFFPtr gff, *childList) {
                if (gff) {
                    gff->removeLinks();                
                    gff.reset();
                }
            }
            childList.reset();
        }
        
        if (childMap && !childMap->empty()) {
           BOOST_FOREACH(GFFIdMap::value_type& i, *childMap) {
                if (i.second) {
                    i.second->removeLinks();                
                    i.second.reset();
                }
            }
            childList.reset(); 
        }
        
        
    }
        
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

    string GetParentId() const {
        return parentId;
    }

    void SetParentId(string parent) {
        this->parentId = parent;
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
        return std::abs(end - start) + 1;
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
    
    string GetNote() const {
        return note;
    }

    void SetNote(string note) {
        this->note = note;
    }


    bool IsCircular() const {
        return circular;
    }

    void SetCircular(bool circular) {
        this->circular = circular;
    }
    
    string GetDerivesFrom() const {
        return derivesFrom;
    }

    void SetDerivesFrom(string derivesFrom) {
        this->derivesFrom = derivesFrom;
    }
    
    string GetIndex() const {
        return index;
    }

    void SetIndex(string index) {
        this->index = index;
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

    void addChild(GFFPtr gff) {        
        addChild(gff, false);
    }

    void addChild(GFFPtr gff, bool noMap) {
        
        string childId = gff->id;
        
        if (!noMap) {            
            if (this->childMap->count(childId)) {
                BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                    "Invalid GFF: Already seen this GFF in child map: ") + childId));
            }
            (*(this->childMap))[childId] = gff;        
        }
        
        gff->parent = shared_from_this();
        this->childList->push_back(gff);
    }
    
    size_t GetNbChildren() {
        return this->childList ? this->childList->size() : 0;
    }
    
    GFFListPtr GetChildList() const {
        return childList;
    }
    
    GFFListPtr GetAllChildren() const {
        
        GFFListPtr gffs = make_shared<GFFList>();
        
        if (childList) {
            BOOST_FOREACH(GFFPtr child, *(childList)) {            

                gffs->push_back(child);

                child->GetAllChildren(*gffs);
            }
        }
        
        return gffs;
    }
    
    
    void GetAllChildren(GFFList& gffs) const {
        
        if (childList) {
            BOOST_FOREACH(GFFPtr child, *(childList)) {            

                gffs.push_back(child);

                child->GetAllChildren(gffs);
            }
        }        
    }
    
    GFFListPtr GetAllOfType(GffType type) const {
        
        GFFListPtr gffs = make_shared<GFFList>();
        
        if (childList) {
            BOOST_FOREACH(GFFPtr child, *(childList)) {            

                if (child->GetType() == type) {
                    gffs->push_back(child);
                }

                child->GetAllOfType(type, *gffs);
            }
        }
        
        return gffs;
    }
    
    int32_t GetLengthOfAllTypes(GffType type) const {
        
        GFFListPtr features = this->GetAllOfType(type);
        
        int32_t len = 0;
        
        if (features) {
            // This loop just counts how may genomic CDSes we also find in the cluster aligned file
            BOOST_FOREACH(GFFPtr feature, *features) {
                len += feature->GetLength();
            }
        }
        
        return len;
    }
    
    
    void GetAllOfType(GffType type, GFFList& gffs) const {
        
        if (childList) {
            BOOST_FOREACH(GFFPtr child, *(childList)) {            

                if (child->GetType() == type) {
                    gffs.push_back(child);
                }

                child->GetAllOfType(type, gffs);
            }
        }        
    }

    GFFIdMapPtr GetChildMap() const {
        return childMap;
    }


    GFFPtr GetParent() {
        return this->parent;
    }



    void writeGFF3Attribs(ostream& out) {
        
        vector<string> elems;
        
        elems.push_back(string("ID=") + id);
        
        if (!parentId.empty()) {
            elems.push_back(string("Parent=") + parentId);
        }
        
        if (!name.empty()) {
            elems.push_back(string("Name=") + name);
        }
        
        if (!note.empty()) {
            elems.push_back(string("Note=") + note);
        }
        
        if (!alias.empty()) {
            elems.push_back(string("Alias=") + alias);
        }
                        
        if (!target.empty()) {
            elems.push_back(string("Target=") + target);
        }
        
        if (!gap.empty()) {
            elems.push_back(string("Gap=") + gap);
        }
        
        if (!derivesFrom.empty()) {
            elems.push_back(string("Derives_from=") + derivesFrom);
        }
        
        if (!index.empty()) {
            elems.push_back(string("Index=") + index);
        }

        out << boost::algorithm::join(elems, ";");        
    }
    
    void writeGTFAttribs(ostream& out) {
        
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

    // Sort by: seqid / start / end / type
    struct GFFOrdering {
        inline bool operator ()(const GFFPtr& a, const GFFPtr& b) {

            int seqId = a->GetSeqId().compare(b->GetSeqId());
            if (seqId != 0) {
                return seqId < 0;
            }
            else {
                int32_t sDiff = a->GetStart() - b->GetStart();
                if (sDiff != 0) {
                    return a->GetStart() < b->GetStart();
                }
                else {
                    int32_t eDiff = a->GetEnd() - b->GetEnd();

                    if (eDiff != 0) {
                        // We actually want the record with the largest end point to come
                        // first... i.e. gene before exon
                        return a->GetEnd() > b->GetEnd();
                    }
                    else {
                        return a->GetType() < b->GetType();
                    }
                }
            }
        }
    };
    
    void write(ostream& out) {
        write(out, false);
    }
    
    void write(ostream& out, bool writeChildren) {
        write(out, this->source, writeChildren);
    }

    void write(ostream& out, string newSource, bool writeChildren) {
        
        out << seqId << "\t" 
            << newSource << "\t"
            << gffTypeToString(type) << "\t"
            << start << "\t"
            << end << "\t"
            << (score == -1.0 ? "." : lexical_cast<string>(score)) << "\t"
            << strand << "\t"
            << (phase == -1 ? "." : lexical_cast<string>(phase)) << "\t";
            
        if (fileFormat == GFF3) {
            writeGFF3Attribs(out);
        }
        else if (fileFormat == GTF) {
            writeGTFAttribs(out);
        }
        
        out << endl;
            
        if (writeChildren) {

            std::sort(this->childList->begin(), this->childList->end(), GFF::GFFOrdering());

            BOOST_FOREACH(GFFPtr child, *(this->childList)) {
                child->write(out, newSource, true);
            }
        }
    }
    
    
    static GFFPtr parse(FileFormat fileFormat, const string& line) {
        vector<string> parts;
        boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_on );

        if (parts.size() != 9) {
            BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                "Could not parse GFF line due to incorrect number of columns. Expected 9 columns: ") + line));
        }
        
        GFFPtr gff = make_shared<GFF>(fileFormat);
        
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
                        gff->SetParentId(val);
                    }
                    else if (boost::iequals(key, "Alias")) {
                        gff->SetAlias(val);
                    }
                    else if (boost::iequals(key, "Note")) {
                        gff->SetNote(val);
                    }
                    else if (boost::iequals(key, "Target")) {
                        gff->SetTarget(val);
                    }
                    else if (boost::iequals(key, "Gap")) {
                        gff->SetGap(val);
                    }
                    else if (boost::iequals(key, "Derives_from")) {
                        gff->SetDerivesFrom(val);
                    }
                    else if (boost::iequals(key, "Index")) {
                        gff->SetIndex(val);
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
    
    
    static void load(FileFormat fileFormat, const string& path, GFFList& gffs) {
    
        load(fileFormat, path, gffs, ANY);
    }
    
    static void load(FileFormat fileFormat, const string& path, GFFList& gffs, GffType filter) {
    
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n");
        cout << " - Loading GFF: " << path << endl;
        
        if (filter != ANY) {
            cout << " - Keeping only : " << gffTypeToString(filter) << endl;
        }
        
        
        uint32_t totalCount = 0;
        ifstream file(path.c_str());
        string line; 
        while (std::getline(file, line)) {            
            boost::trim(line);
            if (!line.empty()) {
                
                GFFPtr gff = parse(fileFormat, line);
                totalCount++;
                
                if (filter == ANY || gff->GetType() == filter) {
                    gffs.push_back(gff);
                }
            }
        }
        file.close();
        
        cout << " - Loaded " << gffs.size() << " out of " << totalCount << " GFF records." << endl;
    }
    
    static void save(const string& path, GFFList& gffs) {
        save(path, gffs, string(""));
    }
    
    static void save(const string& path, GFFList& gffs, string source) {
        
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        cout << " - Saving to: " << path << endl;
        
        ofstream file(path.c_str());
        BOOST_FOREACH(GFFPtr gff, gffs) {
            
            if (source.empty()) {
                gff->write(file);
            }
            else {
                gff->write(file, source, false);
            }
        }
        file.close();
    }    
    
};





class GFFModel {

private:
    GFFListPtr geneList;
    GFFIdMap geneMap;
    GFFIdMap transcriptMap;

public:

    GFFModel() {
        geneList = make_shared<GFFList>();
    }
    
    virtual ~GFFModel() {
        
        BOOST_FOREACH(GFFPtr gff, *geneList) {
            if (gff) {
                gff->removeLinks();
                gff.reset();
            }
        }
    }
    
    GFFPtr getGeneByIndex(uint32_t index) {
        return (*geneList)[index];
    }
    
    bool containsGene(string id) {
        return geneMap.count(id);
    }

    GFFPtr getGeneById(string& id) {
        return geneMap[id];
    }
    
    bool containsTranscript(string id) {
        return transcriptMap.count(id);
    }
    
    GFFPtr getTranscriptById(string id) {
        return transcriptMap[id];
    }
    
    size_t getTotalNbTranscripts() {
        return transcriptMap.size();
    }
    
    size_t getNbTranscripts(string& id) {
        return geneMap[id]->GetNbChildren();
    }
    
    size_t getNbGenes() {
        return geneList->size();
    }
    
    GFFListPtr getGeneList() {
        return geneList;
    }
    
    GFFListPtr getFullList() {
        
        GFFListPtr gffs = make_shared<GFFList>();
        
        BOOST_FOREACH(GFFPtr gene, *(this->geneList)) {
            
            gffs->push_back(gene);            
            gene->GetAllChildren(*gffs);            
        }
        
        return gffs;        
    }
    
    GFFListPtr getAllOfType(GffType type) {
        
        GFFListPtr gffs = make_shared<GFFList>();
        
        GFFListPtr genes = this->geneList;
        
        if (genes) {
            BOOST_FOREACH(GFFPtr gene, *genes) {
                gffs->push_back(gene);            
                gene->GetAllOfType(type, *gffs);            
            }
        }
        
        return gffs;        
    }
        
    void addGene(GFFPtr gff) {
        
        string id = gff->GetId();
        
        if (gff->GetType() != GENE) {
            BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                "The GFF provided does not represent a gene: ") + id));
        }
        
        if (this->geneMap.count(id)) {
            BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                "Invalid GFF: Already loaded gene with this Id: ") + id));
        }
        else {
            this->geneMap[id] = gff;
            this->geneList->push_back(gff);
            
            BOOST_FOREACH(GFFPtr transcript, *(gff->GetChildList())) {
                this->transcriptMap[transcript->GetId()] = transcript;
            }
        }
    }
    
    void rebuildGeneMap() {
       
        this->geneMap.clear();
        
        BOOST_FOREACH(GFFPtr gene, *(this->geneList)) {
            this->geneMap[gene->GetId()] = gene;
        }
    }
    
    static GFFModelPtr load(const string& path) {
        
        GFFModelPtr geneModel = make_shared<GFFModel>();
        
        GFFList gffs;
        GFF::load(GFF3, path, gffs);
        
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n");
        cout << " - Linking GFF records to create gene model" << endl;
        
        BOOST_FOREACH(GFFPtr gff, gffs) {
            
            string id = gff->GetId();
            
            if (gff->GetType() == GENE) {
                geneModel->addGene(gff);
            }
            else if (gff->GetType() == MRNA || gff->GetType() == MIRNA) {
                
                string parent = gff->GetParentId();
            
                // We assume the gene is already present... should be in most GFFs
                if (geneModel->geneMap.count(parent) > 0) {                    
                    geneModel->geneMap[parent]->addChild(gff);
                    geneModel->transcriptMap[id] = gff;
                }
                else {
                    BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                        "Invalid GFF: Could not find parent gene for mRNA: ") + id));
                }
            }
            else if (gff->GetType() == PROTEIN) {
                
                cerr << "Ignoring protein: " << id << endl;
                
                /*string derivesFrom = gff->GetDerivesFrom();
                
                // We assume the gene is already present... should be in most GFFs
                if (geneModel->transcriptMap.count(derivesFrom) > 0) {                    
                    geneModel->transcriptMap[derivesFrom]->addChild(gff, true);
                }
                else {
                    BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                        "Invalid GFF: Could not find parent gene for mRNA: ") + id));
                }*/
            }
            else {
            
                string parent = gff->GetParentId();
            
                vector<string> parents;
                boost::split( parents, parent, boost::is_any_of(","), boost::token_compress_off );
                
                vector<string> filteredParents;
                
                BOOST_FOREACH(string p, parents) {
                    
                    if (p.find("-Protein") == std::string::npos) {
                        filteredParents.push_back(p);
                    }
                }
                
                if (filteredParents.size() > 1) {
                    cerr << "Ignoring GFF entry: id-" << id << "; type-" << gff->GetType() << endl;
                }
                else {                
                
                    if (geneModel->transcriptMap.count(filteredParents[0]) > 0) {
                        geneModel->transcriptMap[filteredParents[0]]->addChild(gff, true);                        
                    }
                    else {
                        BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                            "Invalid GFF: Could not find parent transcript for GFF entry: ") + id));                    
                    }
                }
            }
        }
        
        cout << " - Found " << geneModel->getNbGenes() << " genes and " << geneModel->getTotalNbTranscripts() << " transcripts" << endl;
        
        return geneModel;
    }
    
    void save(const string path) {
        save(path, false, string(""));
    }
    
    void save(const string path, const bool sort) {
        save(path, sort, string(""));
    }
    
    void save(const string path, const bool sort, const string source) {
        
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        cout << " - Saving to: " << path << endl;
        
        if (this->geneList) {
            if (this->geneList->empty()) {
                cerr << "No GFFs to save!" << endl;
            }
            else {
            
                const string s = source.empty() ? this->geneList->at(0)->GetSource() : source;

                if (sort) {
                    cout << " - Sorting GFF records" << endl;
                    std::sort(this->geneList->begin(), this->geneList->end(), GFF::GFFOrdering());
                }

                ofstream file(path.c_str());

                BOOST_FOREACH(GFFPtr gene, *(this->geneList)) {

                    gene->write(file, s, true);

                    // Separate genes with an extra line
                    file << endl;

                }
                file.close();                            
            }
        }
    }
    
};

}
}