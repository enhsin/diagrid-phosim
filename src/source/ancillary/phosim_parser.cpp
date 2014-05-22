///
/// @package phosim
/// @file phosim_parser.cpp
/// @brief Parse configuration files formatted as key-value pairs.
/// Provide PhosimPar class to provide implicit type conversions.
///
/// @brief Created by:
/// @author J. Chiang
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <fstream>
#include <utility>
#include <iostream>
#include <iomanip> 

#include "ancillary/phosim_parser.h"

namespace ancillary {

    PhosimParser::PhosimParser() {}

    PhosimParser::PhosimParser(std::istream & in_stream) {
        read_stream(in_stream);
    }

    PhosimParser::PhosimParser(const std::string & parfile) {
        read_file(parfile);
    }

    void PhosimParser::read_stream(std::istream & in_stream) {
        // Each line is parsed as a key-value pair assuming space-delimited
        // fields.  The first field is the key and the remaining fields are
        // re-concatenated as a single string.
        std::string line;
        while (std::getline(in_stream, line, '\n')) {
            std::vector<std::string> tokens;
            stringTokenize(line, " ", tokens);
            if (tokens.size() >= 2) {
                // Reject lines with 0 or 1 tokens.
                std::string key(tokens[0]);
                std::string value(tokens[1]);
                for (size_t i(2); i < tokens.size(); i++) {
                    value += " " + tokens[i];
                }
                if (has_key(key)) {
                    m_data[key].push_back(PhosimPar(value));
                } else {
                    std::vector<PhosimPar> values;
                    values.push_back(PhosimPar(value));
                    m_data.insert(std::make_pair(tokens[0], values));
                }
            }
        }
    }
    // ********************************************************************

    void PhosimParser::readStream(std::istream& inStream) 
    // ***********************************************************************
    // This method is just like read_stream except it ignores comments
    // ***********************************************************************
    // Each line is parsed as a key-value pair assuming space-delimited
    // fields.  The first field is the key and the remaining fields are
    // re-concatenated as a single string.
    // ************************************************************************
    // Note we also have ignore comments (anything beyond a # sign , may be 
    // whole lines
    {
        // ***********************************************************************
        // Each line is parsed as a key-value pair assuming white space delimited
        // fields.  (This allows for spaces or tab to seperate values. Both are used
        // in lsst files.) 
        // For Commands: the first field is the key and the remaining fields are
        // re-concatenated as a single string.
        // ***********************************************************************
        std::string line;

        while (std::getline(inStream, line, '\n')) {
            // ********************************************************************
            // seach the line for a # sign It deniotes all following chars are a 
            // comment. Delete the # and all followin chars from the string. 
            // *********************************************************************
            std::vector<std::string> tokens;
            std::string value;

            bool goodLine=stringTokenize(line,tokens);
            if(!goodLine){
                continue;
            }

            int numTokens=tokens.size();

            size_t j=1;
            std::string key=tokens[0];     //Normal single key, single value ocmmand
            if (numTokens > 1) {
                value=tokens[j];
            } else {
                value=tokens[0];
            }
            j++;
            for (size_t i=j; i <(size_t) numTokens; i++) {
                value += " " + tokens[i];
            }
            if (has_key(key)) {
                m_data[key].push_back(PhosimPar(value));
            } else {
                std::vector<PhosimPar> values;
                values.push_back(PhosimPar(value));
                m_data.insert(std::make_pair(key, values));
            }
        }
        return;
    }
    // ************************************************************************

    void PhosimParser::readCommandStream(std::istream& inStream, int& numExtra) 
    // ***********************************************************************
    // This version parses a stream that may have 2 types of commands:
    // A. Commands with a command verb and a single value
    // B. Commands that are "control" commands that have body and zernike
    //    specs. These have multiple values. We use a concatination of the first
    //    (Ex: M1 or L1 or F1E )with the tag (Ex: psi or theta or zdis or z11 ) 
    //    elements on the line (with a space between) as the string "key" in
    //    the m_data map.
    // ************************************************************************
    // Note we also have ignore comments (anything beyond a # sign , may be 
    // whole lines
    {
        // ***********************************************************************
        // Each line is parsed as a key-value pair assuming white space delimited
        // fields.  (This allows for spaces or tab to seperate values. Both are used
        // in lsst files.) 
        // For Commands: the first field is the key and the remaining fields are
        // re-concatenated as a single string.
        // For control statements: The first 2 are tags, concatinated with a 
        // seperating space as the key and all the remaining fields are
        // re-concatenated as a single string.
        // ***********************************************************************
        std::string line;
 
        numExtra=0;
        while (std::getline(inStream, line, '\n')) {
            // ********************************************************************
            // seach the line for a # sign It denotes all following chars are a 
            // comment. Delete the # and all following chars from the string. 
            // *********************************************************************
            std::vector<std::string> tokens;

            bool goodLine=stringTokenize(line,tokens);
            if(!goodLine){
                continue;
            }

            int numTokens=tokens.size();
    
            size_t j=1;
            std::string key=tokens[0];     //Normal single key, single value ocmmand
            if ( numTokens > 2 ){      //If a "control" comand make the key a
                key=key+" "+tokens[1];   //contatination(seperated by a space) of the
                j=2;                     //first 2 strings. Rest re-concatenated as a 
                numExtra++;
            }                          //single string (space delimited).
    
            std::string value=tokens[j];
            j++;
            for (size_t i=j; i <(size_t) numTokens; i++) {
                value += " " + tokens[i];
            }
            if (has_key(key)) {
                m_data[key].push_back(PhosimPar(value));
            } else {
                std::vector<PhosimPar> values;
                values.push_back(PhosimPar(value));
                m_data.insert(std::make_pair(key, values));
            }
        }
        return;
    }
    // ************************************************************************

    void PhosimParser::readSegmentation(std::istream& in_stream) {
        // Read segmentation.txt
        std::string line;
        while (std::getline(in_stream, line, '\n')) {
            std::vector<std::string> tokens;
            bool goodLine = stringTokenize(line,tokens);
            if(!goodLine){
                continue;
            }
            std::string chipID(tokens[0]);
            std::string value(tokens[1]);
            for (size_t i(2); i < tokens.size(); i++) {
                value += " " + tokens[i];
            }
            std::vector<PhosimPar> values;
            values.push_back(PhosimPar(value));
            m_data.insert(std::make_pair(chipID, values));
            int numAmplifiers = std::atoi(tokens[1].c_str());
            for (int j(0); j < numAmplifiers; j++) {
                std::string line2;
                std::getline(in_stream, line2, '\n');
                std::vector<std::string> tokens2;
                bool goodLine = stringTokenize(line2,tokens2);
                if(!goodLine){
                    j--;
                    continue;
                }
                std::string value2(tokens2[0]);
                for (size_t i(1); i < tokens2.size(); i++) {
                    value2 += " " + tokens2[i];
                }
                m_data[chipID].push_back(PhosimPar(value2));
            }
        }
    }

    void PhosimParser::read_file(const std::string & parfile) {
        std::ifstream in_stream(parfile.c_str());
        readStream(in_stream);
    }

    bool PhosimParser::has_key(const std::string & key) const {
        std::map<std::string, std::vector<PhosimPar> >::const_iterator 
            it(m_data.find(key));
        if (it == m_data.end()) {
            return false;
        }
        return true;
    }

    const PhosimPar & PhosimParser::operator[](const std::string & key) const {
        return getVector(key).back();
    }

    const std::vector<PhosimPar> & 
    PhosimParser::getVector(const std::string & key) const {
        std::map<std::string, std::vector<PhosimPar> >::const_iterator 
            it(m_data.find(key));
        if (it == m_data.end()) {
            std::ostringstream message;
            message << "PhosimParser::operator[]: " 
                    << "Requested key " << key << " not found.";
            throw PhosimParserException(message.str());
        }
        return it->second;
    }

    void PhosimParser::getNameVector(const std::string & key,
                                     std::vector<std::string> & names) const {
        const std::vector<PhosimPar> & values(getVector(key));
        if (names.size() < values.size()) names.resize(values.size());
        for (size_t i(0); i < values.size(); i++) {
            std::vector<std::string> tokens;
            stringTokenize(values[i], " ", tokens);
            names[std::atoi(tokens[0].c_str())] = tokens[1];
        }
    }

    void PhosimParser::getNameVector(const std::string & key,
                                     std::vector<int> & names) const {
        const std::vector<PhosimPar> & values(getVector(key));
        if (names.size() < values.size()) names.resize(values.size());
        for (size_t i(0); i < values.size(); i++) {
            std::vector<std::string> tokens;
            stringTokenize(values[i], " ", tokens);
            names[std::atoi(tokens[0].c_str())] = atoi(tokens[1].c_str());
        }
    }

    void PhosimParser::getNameVector(const std::string & key,
                                     std::vector<float> & names) const {
        const std::vector<PhosimPar> & values(getVector(key));
        if (names.size() < values.size()) names.resize(values.size());
        for (size_t i(0); i < values.size(); i++) {
            std::vector<std::string> tokens;
            stringTokenize(values[i], " ", tokens);
            names[std::atoi(tokens[0].c_str())] = atof(tokens[1].c_str());
        }
    }

    void PhosimParser::getNameVector(const std::string & key,
                                     std::vector<double> & names) const {
        const std::vector<PhosimPar> & values(getVector(key));
        if (names.size() < values.size()) names.resize(values.size());
        for (size_t i(0); i < values.size(); i++) {
            std::vector<std::string> tokens;
            stringTokenize(values[i], " ", tokens);
            names[std::atoi(tokens[0].c_str())] = atof(tokens[1].c_str());
        }
    }

    void PhosimParser::getNameMatrix(const std::string & key,
                                     std::vector<std::vector<double> > & names) const {
        const std::vector<PhosimPar> & values(getVector(key));
        for (size_t i(0); i < values.size(); i++) {
            std::vector<std::string> tokens;
            stringTokenize(values[i], " ", tokens);
            names[std::atoi(tokens[0].c_str())][std::atoi(tokens[1].c_str())] = atof(tokens[2].c_str());
        }
    }

    void PhosimParser::getNameList(const std::string & key,
                                   std::vector<std::string> & names) const {
        const std::vector<PhosimPar> & values(getVector(key));
        names.resize(values.size());
        for (size_t i(0); i < values.size(); i++) {
            names[i].assign(values[i]);
        }
    }

    int PhosimParser::getNameListSize(const std::string & key) {
        const std::vector<PhosimPar> & values(getVector(key));
        return values.size();
    }

    int PhosimParser::getNameListSize(const std::string & key, 
                                      int index) {
        const std::vector<PhosimPar> & values(getVector(key));
        int num=0;
        for (size_t i(0); i < values.size(); i++) {
            std::vector<std::string> tokens;
            stringTokenize(values[i], " ", tokens);
            if (std::atoi(tokens[0].c_str()) == index) num++;
        }
        return num;
    }

    int PhosimParser::getNameListSize(const std::string & key,
                                      int index1, int index2) {
        const std::vector<PhosimPar> & values(getVector(key));
        int num=0;
        for (size_t i(0); i < values.size(); i++) {
            std::vector<std::string> tokens;
            stringTokenize(values[i], " ", tokens);
            if (std::atoi(tokens[0].c_str()) == index1 && std::atoi(tokens[1].c_str()) == index2) num++;
        }
        return num;
    }


    int PhosimParser::integer(const std::string & key) const {
        return std::atoi(static_cast<std::string>(operator[](key)).c_str());
    }

    long PhosimParser::long_int(const std::string & key) const {
        return std::atoi(static_cast<std::string>(operator[](key)).c_str());
    }

    float PhosimParser::floating(const std::string & key) const {
        return std::atof(static_cast<std::string>(operator[](key)).c_str());
    }

    double PhosimParser::dfloat(const std::string & key) const {
        return std::atof(static_cast<std::string>(operator[](key)).c_str());
    }

    void PhosimParser::set(const std::string & key, const std::string & value) {
        if (!has_key(key)) {
            std::vector<PhosimPar> values;
            values.push_back(PhosimPar(value));
            m_data.insert(std::make_pair(key, values));
        } else {
            m_data[key].push_back(PhosimPar(value));
        }
    }

    void PhosimParser::get(const std::string & key, std::string & value) const {
        try {
            value = static_cast<std::string>(operator[](key));
        } catch (PhosimParserException & eObj) {
            // do nothing
        }
    }

    void PhosimParser::get(const std::string & key, int & value) const {
        try {
            value = integer(key);
        } catch (PhosimParserException & eObj) {
            // do nothing
        }
    }

    void PhosimParser::get(const std::string & key, long & value) const {
        try {
            value = long_int(key);
        } catch (PhosimParserException & eObj) {
            // do nothing
        }
    }

    void PhosimParser::get(const std::string & key, float & value) const {
        try {
            value = floating(key);
        } catch (PhosimParserException & eObj) {
            // do nothing
        }
    }

    void PhosimParser::get(const std::string & key, double & value) const {
        try {
            value = dfloat(key);
        } catch (PhosimParserException & eObj) {
            // do nothing
        }
    }

    void PhosimParser::stringTokenize(std::string input, 
                                      const std::string & delimiters,
                                      std::vector<std::string> & tokens) {
        std::string::size_type j;
        while ( (j = input.find_first_of(delimiters)) != std::string::npos ) {
            if (j != 0) tokens.push_back(input.substr(0, j));
            input = input.substr(j+1);
        }
        tokens.push_back(input);
        if (tokens.back() == "") {
            tokens.pop_back();
        }
    }
    // *************************************************************************

    bool PhosimParser::stringTokenize(std::string line, 
                                      std::vector<std::string>& tokens) 
    // **********************************************************************
    // Get all the tokens from a line using white space as a delimiter.
    // **********************************************************************
    {
        // ********************************************************************
        // seach the line for a # sign It deniotes all following chars are a 
        // comment. Delete the # and all followin chars from the string. 
        // *********************************************************************
        std::string::size_type idx=line.find("#");
        if(idx!=std::string::npos){           //ok there is a comment in the line
            if(idx!=0){
                line=line.substr(0,idx-1);
            }
            else{
                return false;
            }
        }
        std::istringstream iss(line);
        std::string value;
        while(iss>>value){     //This uses white space (see above) to 
            tokens.push_back(value);  //seperate values
        }
        int numTokens=tokens.size();
        if(numTokens==0){         //blank line
            return false;
        }
        return true;
    }
    // **********************************************************************

    bool PhosimParser::getKey(int keyLocation, std::string& thisKey)
    // **********************************************************
    // Get key at location n returns false if beyon m_data.end()
    // **********************************************************
    {
        // *********************************************************
        //Check we haven't requested to go beyiond the end of m_data
        // *********************************************************
        if(keyLocation>=(int)m_data.size() || keyLocation<0 ){
            return false;
        }

        pos=m_data.begin();
        advance(pos,keyLocation);
        thisKey=pos->first;
        return true;
    }
    // ************************************************************
    int PhosimParser::removeKey(std::string & key)
    {
        int numRemoved=m_data.erase(key);
        return numRemoved;
    }
  
} // namespace ancillary 
