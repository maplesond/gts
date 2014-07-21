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

typedef boost::error_info<struct GFFError,string> GFFErrorInfo;
struct GFFException: virtual boost::exception, virtual std::exception { };

enum GffType {
    MRNA,
    UTR5,
    UTR3,
    CDS,
    TRANSCRIPT,
    OTHER
};

static GffType gffTypeFromString(string& s) {
    
    if (s.compare("mRNA") == 0) {
        return MRNA;
    }
    else if (s.compare("five_prime_utr") == 0) {
        return UTR5;
    }
    else if (s.compare("three_prime_utr") == 0) {
        return UTR3;
    }
    else if (s.compare("CDS") == 0) {
        return CDS;
    }
    else if (s.compare("TRANSCRIPT") == 0) {
        return TRANSCRIPT;
    }
    else {
        return OTHER;
    }
}

class GFF {
    
private:
    
    string target;
    string tool;
    
    GffType type;
    int32_t start;
    int32_t end;
    
    char strand;
    
    string id;
    string cdsid;
    string name;
    string parent;
    double coverage;
    string identity;
    
public:

    GFF() {
    }

    virtual ~GFF() {}
    
    char GetStrand() const {
        return strand;
    }

    void SetStrand(char strand) {
        this->strand = strand;
    }

    string GetTarget() const {
        return target;
    }

    void SetTarget(string target) {
        this->target = target;
    }

    string GetTool() const {
        return tool;
    }

    void SetTool(string tool) {
        this->tool = tool;
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
    
    string GetCdsid() const {
        return cdsid;
    }

    void SetCdsid(string cdsid) {
        this->cdsid = cdsid;
    }


    string GetIdentity() const {
        return identity;
    }

    void SetIdentity(string identity) {
        this->identity = identity;
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

    int32_t GetStart() const {
        return start;
    }

    void SetStart(int32_t start) {
        this->start = start;
    }

    void write(std::ostream& out) {
        
        std::stringstream ss;
        ss << "ID=" << id << ";";
        
        if (!parent.empty()) {
            ss << "Parent=" << parent << ";";
        }
        
        out << target << "\t" 
            << tool << "\t"
            << type << "\t"
            << start << "\t"
            << end << "\t"
            << "." << "\t"
            << strand << "\t"
            << "." << "\t"
            << ss.str() << endl;
    }
    
    
    static shared_ptr<GFF> parse(const string& line) {
        vector<string> parts;
        boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_on );

        if (parts.size() != 9) {
            BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                "Could not parse GFF line due to incorrect number of columns. Expected 9 columns: ") + line));
        }
        
        shared_ptr<GFF> gff = shared_ptr<GFF>(new GFF());
        
        gff->SetTarget(parts[0]);
        gff->SetTool(parts[1]);
        gff->SetType(gffTypeFromString(parts[2]));
        gff->SetStart(lexical_cast<int32_t>(parts[3]));
        gff->SetEnd(lexical_cast<int32_t>(parts[4]));
        //gff->SetEnd(parts[5][0] == '.' ? -1 : lexical_cast<int32_t>(parts[5]));
        gff->SetStrand(parts[6][0]);
        
        vector<string> attrParts;
        boost::split( attrParts, parts[8], boost::is_any_of(";"), boost::token_compress_on );
        
        BOOST_FOREACH(string attr, attrParts) {
            vector<string> attrElements;
            boost::split( attrElements, attr, boost::is_any_of("="), boost::token_compress_on );
            
            string* key = &attrElements[0];
            string val = attrElements[1];
            if (key->compare("ID") == 0) {
                vector<string> idElements;
                boost::split( idElements, attr, boost::is_any_of("|"), boost::token_compress_on );
                size_t pos = idElements[0].find("cds");
                if (pos != std::string::npos) {
                    gff->SetCdsid(idElements[0].substr(pos+4));
                }
                
                gff->SetId(val);
            }
            else if (key->compare("Name") == 0) {
                gff->SetName(val);
            }
            else if (key->compare("Parent") == 0) {
                gff->SetParent(val);
            }
            else if (key->compare("coverage") == 0) {
                gff->SetCoverage(lexical_cast<double>(val));
            }
            else if (key->compare("identity") == 0) {
                gff->SetIdentity(val);
            }
        
        }
        
        return gff;

    }
    
    
    typedef std::vector< boost::shared_ptr<GFF> > GFFList;
typedef boost::unordered_map<string, shared_ptr<GFF> > GFFIdMap;

    static void load(const string& path, std::vector< boost::shared_ptr<GFF> >& gffs) {
    
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        cout << " - Loading GFF: " << path << endl;
        
        std::ifstream file(path.c_str());
        std::string line; 
        while (std::getline(file, line)) {            
            boost::trim(line);
            if (!line.empty()) {
                gffs.push_back(parse(line));
            }
        }
        file.close();
        
        cout << " - Found " << gffs.size() << " GFF records." << endl;
    }
    
    static void save(const string& path, std::vector< boost::shared_ptr<GFF> >& gffs) {
        
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        cout << " - Saving GFF: " << path << endl;
        
        std::ofstream file(path.c_str());
        BOOST_FOREACH(shared_ptr<GFF> gff, gffs) {
            gff->write(file);
        }
        file.close();
    }
};

}