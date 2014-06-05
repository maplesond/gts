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
//  along with Portculis.  If not, see <http://www.gnu.org/licenses/>.
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
    else {
        return OTHER;
    }
}

class GFF {
    
private:
    
    GffType type;
    
    
    
    string id;
    string name;
    string parent;
    double coverage;
    string identity;
    
public:

    GFF() {
    }

    
    GFF(GffType type, double coverage, string id, string identity, string name, string parent) :
        type(type), coverage(coverage), id(id), identity(identity), name(name), parent(parent)  {
    }

    virtual ~GFF() {}
    
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

    
    static shared_ptr<GFF> parse(const string& line) {
        vector<string> parts;
        boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_on );

        if (parts.size() != 9) {
            BOOST_THROW_EXCEPTION(GFFException() << GFFErrorInfo(string(
                "Could not parse GFF line due to incorrect number of columns. Expected 9 columns: ") + line));
        }
        
        shared_ptr<GFF> gff = shared_ptr<GFF>(new GFF());
        
        gff->SetType(gffTypeFromString(parts[2]));
        
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
                
                gff->SetId(idElements[0]);
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
    
    static void loadGFF(const string& path, vector< shared_ptr<GFF> >& gffs) {
    
        std::ifstream file(path.c_str());
        std::string line; 
        while (std::getline(file, line)) {            
            gffs.push_back(parse(line));
        }
        file.close();
    }
    
    
};

}