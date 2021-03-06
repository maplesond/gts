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
#include <string>
#include <vector>
#include <fstream>
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/timer/timer.hpp>
#include <boost/unordered_map.hpp>
using boost::lexical_cast;
using boost::make_shared;
using boost::shared_ptr;
using boost::timer::auto_cpu_timer;
using boost::unordered_map;

namespace gts{
namespace gb {

typedef boost::error_info<struct GenbankError,string> GenbankErrorInfo;
struct GenbankException: virtual boost::exception, virtual std::exception { };


enum BlockType {
    LOCUS,    
    DEFINITION,
    ACCESSION,    
    VERSION,
    KEYWORDS,
    SOURCE,
    REFERENCE,
    FEATURES,
    BASE_COUNT,
    ORIGIN,
    END_RECORD,
    UNKNOWN_BLOCKTYPE
};

static BlockType blockTypeFromString(string & s) {
    if (boost::equals(s, "LOCUS")) {
        return LOCUS;
    }
    else if (boost::equals(s, "DEFINITION")) {
        return DEFINITION;
    }
    else if (boost::equals(s, "ACCESSION")) {
        return ACCESSION;
    }
    else if (boost::equals(s, "VERSION")) {
        return VERSION;
    }
    else if (boost::equals(s, "KEYWORDS")) {
        return KEYWORDS;
    }
    else if (boost::equals(s, "SOURCE")) {
        return SOURCE;
    }
    else if (boost::equals(s, "REFERENCE")) {
        return REFERENCE;
    }
    else if (boost::equals(s, "FEATURES")) {
        return FEATURES;
    }
    else if (boost::equals(s, "BASE COUNT")) {
        return BASE_COUNT;
    }
    else if (boost::equals(s, "BASE")) {
        return BASE_COUNT;
    }
    else if (boost::equals(s, "ORIGIN")) {
        return ORIGIN;
    }
    else if (boost::equals(s, "//")) {
        return END_RECORD;
    }     
    else {
        return UNKNOWN_BLOCKTYPE;
    }
}

static string blockTypeToString(BlockType type) {
    
    switch(type) {
        case LOCUS:
            return "LOCUS";
        case DEFINITION:
            return "DEFINITION";
        case ACCESSION:
            return "ACCESSION";
        case VERSION:
            return "VERSION";
        case KEYWORDS:
            return "KEYWORDS";
        case SOURCE:
            return "SOURCE";
        case REFERENCE:
            return "REFERENCE";
        case FEATURES:
            return "FEATURES";
        case BASE_COUNT:
            return "BASE COUNT";
        case ORIGIN:
            return "ORIGIN";
        case END_RECORD:
            return "//";
    }

    return "";    
}

enum FeaturesType {
    SOURCE_FEATURE,
    CDS,
    GENE,
    UNKNOWN_FEATURE
};

static FeaturesType featureFromString(string& s) {
    
    if (boost::iequals(s, "gene")) {
        return GENE;
    }
    else if (boost::iequals(s, "source")) {
        return SOURCE_FEATURE;
    }
    else if (boost::iequals(s, "CDS")) {
        return CDS;
    }    
    else {
        return UNKNOWN_FEATURE;
    }
}

static string featureToString(FeaturesType type) {
    
    switch(type) {
        case GENE:
            return "gene";
        case CDS:
            return "CDS";
        case SOURCE_FEATURE:
            return "source";        
    }

    return "";    
}

struct Property {    
    string name;
    string value;
    
    void write(std::ostream& out) {
        out << "                     " << "/" << name << "=" << value << endl;
    }
    
    static shared_ptr<Property> parse(vector<string>& lines) {
        
        shared_ptr<Property> p = make_shared<Property>();        
        
        string line = boost::trim_copy(lines[0]);
        
        vector<string> strVec;
        boost::algorithm::split(strVec, line, boost::is_any_of("="), boost::algorithm::token_compress_on);
        
        p->name = strVec[0].substr(1);
        
        std::stringstream ss;
        ss << strVec[1];
        
        for(size_t i = 1; i < lines.size(); i++) {
            ss << boost::trim_copy(lines[i]);            
        }

        string val = ss.str();
        
        p->value = (val[0] == '\"' && val[val.length() - 1] == '\"') ? 
            val.substr(1, val.length() - 2) :
            val;
        
        return p;
    }
};

struct Block {
    BlockType name;
    vector<string> lines;
    
    void write(std::ostream& out) {
        BOOST_FOREACH(string s, lines) {
            out << s << endl;
        }
    }
};

typedef unordered_map<string, shared_ptr<Property> > PropertyMap;

struct Feature {
    string type;
    string location;
    PropertyMap properties;
    
    void write(std::ostream& out) {
        std::stringstream space;
        
        for(int i = 0; i < 16 - type.length(); i++) {
            space << " ";
        }
        
        out << "     " << type << space.str() << location << endl;
        
        BOOST_FOREACH(PropertyMap::value_type& p, properties) {
            p.second->write(out);
        }
    }
    
    static shared_ptr<Feature> parse(vector<string>& lines) {
        
        /*cout << endl << "Feature lines: " << endl;
        BOOST_FOREACH(string s, lines) {
            cout << s << endl;
        }
        cout << endl;*/
        
        shared_ptr<Feature> f = make_shared<Feature>();        
                
        vector<string> strVec;
        boost::algorithm::split(strVec, lines[0], boost::is_any_of("\t "), boost::algorithm::token_compress_on);
        
        f->type = strVec[0];
        
        std::ostringstream loc;
        loc << strVec[1];
        
        size_t i = 1;
        for(i = 1; i < lines.size(); i++) {
            string line = boost::trim_copy(lines[i]);
            
            if (line[0] != '/') {
                loc << line;                
            }
            else {
                break;
            }
        }
        
        f->location = loc.str();
        
        bool first = true;
        vector<string> property;
        
        for(i; i < lines.size(); i++) {
            
            string line = boost::trim_copy(lines[i]);
            
            if (line[0] == '/' && !first) {                
                shared_ptr<Property> p = Property::parse(property);
                f->properties[p->name] = p;
                property.clear();
                first = false;
            }            
            
            property.push_back(line);            
        }
        
        if (!property.empty()) {
            shared_ptr<Property> p = Property::parse(property);
            f->properties[p->name] = p;
        }
        
        return f;
    }
};

typedef unordered_map<string, shared_ptr<Feature> > FeatureMap;
typedef vector<shared_ptr<Feature> > FeatureList;

struct Features {
    
    string header;
    FeatureMap featureMap;
    FeatureList featureList;
      
    void write(std::ostream& out) {
        
        out << header << endl;
        
        BOOST_FOREACH(shared_ptr<Feature> f, featureList) {
            f->write(out);
        }         
    }
    
    bool noFeatures() {
        return featureMap.size() == 0;
    }
    
    static shared_ptr<Features> parse(Block& block) {
        
            
        /*cout << endl << "Features block: " << endl;
        BOOST_FOREACH(string s, block.lines) {
            cout << s << endl;
        }
        cout << endl;*/
        
        shared_ptr<Features> newFeature = make_shared<Features>();        
        
        vector<string> strVec;
        boost::algorithm::split(strVec, block.lines[0], boost::is_any_of("\t "), boost::algorithm::token_compress_on);
        
        newFeature->header = strVec[1];
        
        vector<string> feature;
        bool first = true;
            
        for(size_t i = 1; i < block.lines.size(); i++) {
            
            vector<string> strVec2;
            string line = boost::trim_copy(block.lines[i]);
            std::replace(line.begin(), line.end(), '\t', ' ');

            // Look for a large gap not starting at the beginning of the string... 
            // this is probably a feature header
            if (line.find("    ") != string::npos && !first) {                
                shared_ptr<Feature> f = Feature::parse(feature);
                newFeature->featureList.push_back(f);
                                    
                if (newFeature->featureMap.count(f->type)) {
                    //std::cerr << "Found " << f->type << " before" << endl;
                }
                else {
                    newFeature->featureMap[f->type] = f;
                    feature.clear();
                    first = true;
                }
            }
            else {
                first = false;                
            }
            
            feature.push_back(line);            
        }
        
        if (!feature.empty()) {
            shared_ptr<Feature> f = Feature::parse(feature);                
            newFeature->featureList.push_back(f);
            
            if (newFeature->featureMap.count(f->type)) {
                //std::cerr << "Found " << f->type << " before" << endl;
            }
            else {
                newFeature->featureMap[f->type] = f;                
            }
        }
                
        return newFeature;
    }
};

class Genbank {
    
private:
    
    vector<shared_ptr<Block> > blocks;
    
    shared_ptr<Block> locus;
    shared_ptr<Block> definition;
    shared_ptr<Block> accession;
    shared_ptr<Block> version;
    shared_ptr<Block> keywords;
    shared_ptr<Block> baseCount;
    shared_ptr<Block> source;
    vector<shared_ptr<Block> > references;
    shared_ptr<Features> features;
    shared_ptr<Block> origin;
    
    
public:

    Genbank() {
    }

    virtual ~Genbank() {}
    

    shared_ptr<Block> getAccession() const {
        return accession;
    }

    void setAccession(shared_ptr<Block> accession) {
        this->accession = accession;
    }

    shared_ptr<Block> getBaseCount() const {
        return baseCount;
    }

    void setBaseCount(shared_ptr<Block> baseCount) {
        this->baseCount = baseCount;
    }

    vector<shared_ptr<Block> > getBlocks() const {
        return blocks;
    }

    void setBlocks(vector<shared_ptr<Block> > blocks) {
        this->blocks = blocks;
    }

    shared_ptr<Block> getDefinition() const {
        return definition;
    }

    void setDefinition(shared_ptr<Block> definition) {
        this->definition = definition;
    }
    
    /*shared_ptr<Feature> getFeature(uint16_t index) const {
        return features->features[index];
    }*/

    shared_ptr<Features> getFeatures() const {
        return features;
    }

    void setFeatures(shared_ptr<Features> features) {
        this->features = features;
    }

    shared_ptr<Block> getKeywords() const {
        return keywords;
    }

    void setKeywords(shared_ptr<Block> keywords) {
        this->keywords = keywords;
    }

    shared_ptr<Block> getLocus() const {
        return locus;
    }

    void setLocus(shared_ptr<Block> locus) {
        this->locus = locus;
    }

    shared_ptr<Block> getOrigin() const {
        return origin;
    }

    void setOrigin(shared_ptr<Block> origin) {
        this->origin = origin;
    }

    vector<shared_ptr<Block> > getReferences() const {
        return references;
    }

    void setReferences(vector<shared_ptr<Block> > references) {
        this->references = references;
    }

    shared_ptr<Block> getSource() const {
        return source;
    }

    void setSource(shared_ptr<Block> source) {
        this->source = source;
    }

    shared_ptr<Block> getVersion() const {
        return version;
    }

    void setVersion(shared_ptr<Block> version) {
        this->version = version;
    }

    
    void write(std::ostream& out) {
        
        // Just write the blocks
        BOOST_FOREACH(shared_ptr<Block> b, blocks) {
            b->write(out);
        }                
    }
    
private:
    
    static bool readBlock(ifstream& in, string& currentLine, Block& block) {
        
        std::istringstream iss(currentLine);
        string word;
        iss >> word;
         
        BlockType bt = blockTypeFromString(word);
        
        if (bt == UNKNOWN_BLOCKTYPE) {
            BOOST_THROW_EXCEPTION(GenbankException() << GenbankErrorInfo(string(
                        "Unknown block type detected: ") + word));
        }
        
        block.name = bt;        
        block.lines.push_back(currentLine);
        
        if (bt != END_RECORD) {
            
            string line;
            while (getline(in, line)) {

                if (line.empty() || line[0] == ' ' || line[0] == '\t') {
                    // Add this line to the block
                    block.lines.push_back(line); 
                }
                else {
                
                    std::istringstream iss(line);
                    string word;
                    iss >> word;

                    BlockType bt = blockTypeFromString(word);

                    if (bt != UNKNOWN_BLOCKTYPE) {
                        // We got to another block (or something unknown) so override 
                        // the current line and end
                        currentLine = line;
                        return true;
                    }                    
                }
            }
            
            // End of file
            return false;
        }
        
        return true;
    }
    
       
    
    
public:
   
    static shared_ptr<Genbank> readRecord(ifstream& in) {
        
        shared_ptr<Genbank> gb = make_shared<Genbank>();
        
        bool done = false;
        
        string line;
        
        if (!getline(in, line)) {
            // Reached end of file
            return shared_ptr<Genbank>();
        }
        
                    
        while (!done) {

            shared_ptr<Block> b = make_shared<Block>();
            done = !readBlock(in, line, *b);
            
            // Exit loop if this is an end record marker
            if (b->name == END_RECORD) {
                done = true;
            }
            else {
            
                switch (b->name) {
                    case LOCUS:
                        gb->locus = b;
                        break;
                    case DEFINITION:
                        gb->definition = b;
                        break;
                    case ACCESSION:
                        gb->accession = b;
                        break;
                    case VERSION:
                        gb->version = b;
                        break;
                    case KEYWORDS:
                        gb->keywords = b;
                        break;
                    case REFERENCE:
                        gb->references.push_back(b);
                        break;
                    case FEATURES:
                        gb->features = Features::parse(*b);
                        break;
                    case BASE_COUNT:
                        gb->baseCount = b;
                        break;
                    case ORIGIN:
                        gb->origin = b;
                        break;
                }
                
                gb->blocks.push_back(b);                
            }
            
        }
        
        return gb;
    }
    
    
    static void load(const string& path, std::vector< boost::shared_ptr<Genbank> >& genbank) {
    
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        cout << " - Loading Genbank: " << path << endl;
        
        ifstream file(path.c_str());        
        
        while (file.good()) {            
            shared_ptr<Genbank> gb = readRecord(file);
            if (gb != NULL) {                
                genbank.push_back(gb);                
            }
        }
        file.close();
        
        cout << " - Loaded " << genbank.size() << " genbank records." << endl;
    }
    
    static void save(const string& path, std::vector< boost::shared_ptr<Genbank> >& genbank) {
        
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        cout << " - Output file: " << path << endl;
        
        bool first = true;
        ofstream file(path.c_str());
        BOOST_FOREACH(shared_ptr<Genbank> gb, genbank) {
            
            if (first) {
                first = false;
            }
            else {
                file << "//" << endl;
            }
            gb->write(file);
        }
        file.close();
    }
};

}
}