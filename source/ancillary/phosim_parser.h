///
/// @package phosim
/// @file PhosimParser.h
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

#ifndef ancillary_phosim_parser_h
#define ancillary_phosim_parser_h

#include <cstdlib>
#include <istream>
#include <sstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace ancillary {

    class PhosimPar {
    public:
    PhosimPar() : m_value("") {}

        explicit PhosimPar(const std::string & value) : m_value(value) {}
        operator std::string() const {
            return m_value;
        }
        operator int() const {
            return std::atoi(m_value.c_str());
        }
        operator long() const {
            return std::atoi(m_value.c_str());
        }
        operator float() const {
            return std::atof(m_value.c_str());
        }
        operator double() const {
            return std::atof(m_value.c_str());
        }
        bool operator==(const std::string & other) const {
            return other==m_value;
        }
        bool operator==(int other) const {
            return other==std::atoi(m_value.c_str());
        }
        const char * c_str() const {
            return m_value.c_str();
        }
    private:
        std::string m_value;
    };

    class PhosimParser {

    public:

        PhosimParser();

        PhosimParser(std::istream & in_stream);

        PhosimParser(const std::string & parfile);

        /// The typical use case for this in client code would be
        /// read_stream(std::cin).
        void read_stream(std::istream & in_stream);

        void readStream(std::istream & in_stream);//ignores comments

        //And a read especially for the control.txt file
        void readCommandStream(std::istream& inStream, int& numExtra);

        //for segmentation.txt
        void readSegmentation(std::istream & in_stream);

        /// This allows from several files to be read in succession,
        /// obviating the need to concatenate them into a single input
        /// file.
        void read_file(const std::string & parfile);

        bool has_key(const std::string & key) const;

        const PhosimPar & operator[](const std::string & key) const;

        const std::vector<PhosimPar> & getVector(const std::string & key) const;

        void getNameVector(const std::string & key,
                           std::vector<std::string> & names) const;
        void getNameVector(const std::string & key,
                           std::vector<int> & names) const;
        void getNameVector(const std::string & key,
                           std::vector<float> & names) const;
        void getNameVector(const std::string & key,
                           std::vector<double> & names) const;
        void getNameMatrix(const std::string & key,
                           std::vector<std::vector<double> > & names) const;
        void getNameList(const std::string & key,
                         std::vector<std::string> & names) const;

        int getNameListSize(const std::string & key);
        int getNameListSize(const std::string & key, int index);
        int getNameListSize(const std::string & key, int index1, int index2);

        void set(const std::string & key, const std::string & value);

        /// If the requested key does not exist, the get functions leave
        /// value unchanged.
        void get(const std::string & key, std::string & value) const;
        void get(const std::string & key, int & value) const;
        void get(const std::string & key, long & value) const;
        void get(const std::string & key, float & value) const;
        void get(const std::string & key, double & value) const;

        size_t getNumberKeys(){return m_data.size();};
        bool getKey(int keyLocation, std::string& thisKey);

        bool stringTokenize(std::string input,
                            std::vector<std::string>& tokens);
        static void stringTokenize(std::string input,
                                   const std::string & delimiters,
                                   std::vector<std::string> & tokens);
        int removeKey(std::string & key);

    private:

        std::map<std::string, std::vector<PhosimPar> > m_data;
        std::map<std::string, std::vector<PhosimPar> >::iterator pos;

        int integer(const std::string & key) const;
        long long_int(const std::string & key) const;
        float floating(const std::string & key) const;
        double dfloat(const std::string & key) const;

    };

    class PhosimParserException : public std::runtime_error {
    public:
    PhosimParserException(const std::string & what) 
        : std::runtime_error(what) {}
    };

} // namespace ancillary

typedef std::map<std::string, ancillary::PhosimPar> PhosimParMap;

#endif // ancillary_phosim_parser_h
