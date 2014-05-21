///
/// @package phosim
/// @file test_phosim_parser.cpp
/// @brief Unit tests for PhosimParser class.
///
/// @brief Created by:
/// @author J. Chiang
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <cassert>
#include <cmath>

#include <sstream>
#include <iostream>

#include "ancillary/phosim_parser.h"

#define ASSERT_EQUALS(X, Y) assert(fabs( (X - Y)/Y ) < 1e-5)

void test_stringTokenize() {
    /// Testing the general version.
    std::string test_string("1  \t2 \t 3");
    std::vector<std::string> tokens;
    ancillary::PhosimParser::stringTokenize(test_string, " \t", tokens);
    assert(tokens[0] == "1");
    //    if (tokens[0] == "1") std::cout << "PhosimParser_test1               Pass\n";
    //    else std::cout << "PhosimParser_test1               Fail\n";
    assert(tokens[1] == "2");
    assert(tokens[2] == "3");

    /// Testing the redundant version that requires an object and which only
    /// works with whitespace.
    ancillary::PhosimParser parser;
    tokens.clear();
    parser.stringTokenize(test_string, tokens);
    assert(tokens[0] == "1");
    assert(tokens[1] == "2");
    assert(tokens[2] == "3");
}

void test_PhosimParser() {
    ancillary::PhosimParser parser(std::cin);
    int my_integer;
    float my_float;
    assert(parser["my_string"] == "foobar");
    parser.get("my_integer", my_integer);
    assert(342 == my_integer);
    parser.get("my_float", my_float);
    ASSERT_EQUALS(2.1, my_float);

    ancillary::PhosimParser file_parser("test.pars");
    assert(file_parser["my_string"] == "foobar");
    file_parser.get("my_integer", my_integer);
    assert(342 == my_integer);
    file_parser.get("my_float", my_float);
    ASSERT_EQUALS(2.1, my_float);

    ancillary::PhosimPar par("1.20");

    ASSERT_EQUALS(float(par), 1.2);
    assert(par == "1.20");

    double x = par;
    std::string y = par;
    std::ostringstream message;
    message << x << "  " << y;
    assert(message.str() == "1.2  1.20");

    const std::vector<ancillary::PhosimPar> & pars(parser.getVector("chipid"));
    assert(pars[0] == "0 R01_S00");
    assert(pars[1] == "1 R01_S01");
    assert(pars[2] == "2 R01_S02");
    assert(pars[3] == "3 R01_S10");
    assert(pars.size() == 9);

    std::vector<std::string> chipid;
    parser.getNameVector("chipid", chipid);
    assert(chipid[0] == "R01_S00");
    assert(chipid[1] == "R01_S01");
    assert(chipid[2] == "R01_S02");
    assert(chipid[3] == "R01_S10");
}

int main(int iargc, char * argv[]) {
    test_stringTokenize();
    test_PhosimParser();
    return 0;
}
