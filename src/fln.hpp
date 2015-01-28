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

    
typedef boost::error_info<struct FLNError,string> FLNErrorInfo;
struct FLNException: virtual boost::exception, virtual std::exception { };

enum FLNStatus {
    INTERNAL,
    COMPLETE,
    PUTATIVE_COMPLETE,
    CTERM,
    NTERM,
    PUTATIVE_CTERM,
    PUTATIVE_NTERM,
    MISASSEMBLED,
    CODING,
    PUTATIVE_CODING,
    UNKNOWN,
    DB_OTHER
};


static FLNStatus flnStatusFromString(string& s) {
    
    if (boost::iequals(s, "Internal")) {
        return INTERNAL;
    }
    else if (boost::iequals(s, "Complete")) {
        return COMPLETE;
    }
    else if (boost::iequals(s, "Putative Complete")) {
        return PUTATIVE_COMPLETE;
    }
    else if (boost::iequals(s, "C-terminus")) {
        return CTERM;
    }
    else if (boost::iequals(s, "N-terminus")) {
        return NTERM;
    }
    else if (boost::iequals(s, "Putative C-terminus")) {
        return NTERM;
    }
    else if (boost::iequals(s, "Misassembled")) {
        return MISASSEMBLED;
    }
    else if (boost::iequals(s, "coding")) {
        return CODING;
    }
    else if (boost::iequals(s, "putative_coding")) {
        return PUTATIVE_CODING;
    }
    else if (boost::iequals(s, "unknown")) {
        return UNKNOWN;
    }
    else {
        return DB_OTHER;
    }
}

class DBAnnot {
private:
    string id;
    int32_t fastaLength;
    FLNStatus status;
    int32_t orfStart;
    int32_t orfEnd;
    int32_t sStart;
    int32_t sEnd;
    
public:
    DBAnnot() : id(""), fastaLength(0), status(DB_OTHER), orfStart(-1), orfEnd(-1), sStart(-1), sEnd(-1) {};
    
    virtual ~DBAnnot() {};
    
    int32_t GetOrfEnd() const {
        return orfEnd;
    }

    void SetOrfEnd(int32_t end) {
        this->orfEnd = end;
    }

    int32_t GetFastaLength() const {
        return fastaLength;
    }

    void SetFastaLength(int32_t fastaLength) {
        this->fastaLength = fastaLength;
    }

    string GetId() const {
        return id;
    }

    void SetId(string id) {
        this->id = id;
    }

    int32_t GetOrfStart() const {
        return orfStart;
    }

    void SetOrfStart(int32_t start) {
        this->orfStart = start;
    }

    FLNStatus GetStatus() const {
        return status;
    }

    void SetStatus(FLNStatus status) {
        this->status = status;
    }
    
    int32_t GetSEnd() const {
        return sEnd;
    }

    void SetSEnd(int32_t sEnd) {
        this->sEnd = sEnd;
    }

    int32_t GetSStart() const {
        return sStart;
    }

    void SetSStart(int32_t sStart) {
        this->sStart = sStart;
    }


    
    static shared_ptr<DBAnnot> parse(const string& line) {
        vector<string> parts;
        boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_off );

        if (parts.size() < 8 || parts.size() > 18) {
            BOOST_THROW_EXCEPTION(FLNException() << FLNErrorInfo(string(
                "Could not parse GFF line due to incorrect number of columns. Expected at least 8 columns.  Found ") + 
                    lexical_cast<string>(parts.size()) + " columns.  Line: " + line));
        }
        
        shared_ptr<DBAnnot> db = make_shared<DBAnnot>();
        
        db->SetId(parts[0]);
        db->SetFastaLength(lexical_cast<int32_t>(parts[1]));
        db->SetStatus(flnStatusFromString(parts[4]));
        
        if (db->GetStatus() != MISASSEMBLED && parts.size() >= 14) {
            db->SetOrfStart(parts[12].empty() ? -1 : lexical_cast<int32_t>(parts[12]));
            db->SetOrfEnd(parts[13].empty() ? -1 : lexical_cast<int32_t>(parts[13]));
            
            if (parts.size() >= 16) {
                int32_t ss = parts[14].empty() ? -1 : lexical_cast<int32_t>(parts[14]);
                int32_t se = parts[15].empty() ? -1 : lexical_cast<int32_t>(parts[15]);

                db->SetSStart(ss < se ? ss : se);
                db->SetSEnd(ss < se ? se : ss);
            }
        }
                  
        return db;

    }
    
    
    
    static void load(const string& path, vector< shared_ptr<DBAnnot> >& dbannots) {
    
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        cout << " - Loading FLN DBAnnot: " << path << endl;
        
        std::ifstream file(path.c_str());
        std::string line;
        std::getline(file, line); // Ignore the header line
        while (std::getline(file, line)) {            
            boost::trim(line);
            if (!line.empty()) {
                dbannots.push_back(parse(line));
            }
        }
        file.close();
        
        cout << " - Found " << dbannots.size() << " DB Annot records." << endl;
    }
    
};

class NonCoding {
    
};

}