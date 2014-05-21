///
/// @package phosim
/// @file test_ancillary.cpp
/// @brief Unit tests for ancillary subpackage classes using Boost.Test.
///
/// @brief Created by:
/// @author J. Chiang
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <cmath>

#include <iomanip>
#include <iostream>
#include <sstream>

#include "raytrace/basic_types.h"
#include "raytrace/rng_mwc.h"

#include "ancillary/fits_utils.h"
#include "ancillary/phosim_parser.h"

#define BOOST_TEST_MODULE ancillary package tests
#include "boost/test/included/unit_test.hpp"
#include "boost/test/results_collector.hpp"
using namespace boost::unit_test;

void log_results(const results_collector_t & results_collector, 
                 std::ostream & outstream=std::cout) {
    const test_results & my_results = 
        results_collector.results(framework::current_test_case().p_id);
    std::ostringstream message;
    message << std::setw(33) << std::setiosflags(std::ios::left) 
            << framework::current_test_case().p_name;
    if (my_results.passed()) {
        message << "Pass";
    } else {
        message << "Fail";
    }
    BOOST_TEST_MESSAGE(message.str());
    //    outstream << message.str() << std::endl;
}

BOOST_AUTO_TEST_SUITE(PhosimParser_tests)

BOOST_AUTO_TEST_CASE(test_stringTokenize) {
    /// Testing the general version.
    std::string test_string("1  \t2 \t 3");
    std::vector<std::string> tokens;
    ancillary::PhosimParser::stringTokenize(test_string, " \t", tokens);
    BOOST_CHECK_EQUAL(tokens[0], "1");
    BOOST_CHECK_EQUAL(tokens[1], "2");
    BOOST_CHECK_EQUAL(tokens[2], "3");

    /// Testing the version that requires an object and which only works
    /// with whitespace.
    ancillary::PhosimParser parser;
    tokens.clear();
    parser.stringTokenize(test_string, tokens);
    BOOST_CHECK_EQUAL(tokens[0], "1");
    BOOST_CHECK_EQUAL(tokens[1], "2");
    BOOST_CHECK_EQUAL(tokens[2], "3");

    log_results(results_collector);
}

BOOST_AUTO_TEST_CASE(test_PhosimParser) {
    ancillary::PhosimParser parser("test.pars");
    int my_integer;
    float my_float;
    BOOST_CHECK(parser["my_string"] == "foobar");
    parser.get("my_integer", my_integer);
    BOOST_CHECK_EQUAL(342, my_integer);
    parser.get("my_float", my_float);
    BOOST_CHECK_CLOSE(2.1, my_float, 1e-5);

    ancillary::PhosimParser file_parser("test.pars");
    BOOST_CHECK(file_parser["my_string"] == "foobar");
    file_parser.get("my_integer", my_integer);
    BOOST_CHECK_EQUAL(342, my_integer);
    file_parser.get("my_float", my_float);
    BOOST_CHECK_CLOSE(2.1, my_float, 1e-5);

    ancillary::PhosimPar par("1.20");

    BOOST_CHECK_EQUAL(static_cast<double>(par), 1.2);
    BOOST_CHECK(par == "1.20");

    double x = par;
    std::string y = par;
    std::ostringstream message;
    message << x << "  " << y;
    BOOST_CHECK_EQUAL(message.str(), "1.2  1.20");

    const std::vector<ancillary::PhosimPar> & pars(parser.getVector("chipid"));
    BOOST_CHECK(pars[0] == "0 R01_S00");
    BOOST_CHECK(pars[1] == "1 R01_S01");
    BOOST_CHECK(pars[2] == "2 R01_S02");
    BOOST_CHECK(pars[3] == "3 R01_S10");
    BOOST_CHECK_EQUAL(pars.size(), 9);

    std::vector<std::string> chipid;
    parser.getNameVector("chipid", chipid);
    BOOST_CHECK_EQUAL(chipid[0], "R01_S00");
    BOOST_CHECK_EQUAL(chipid[1], "R01_S01");
    BOOST_CHECK_EQUAL(chipid[2], "R01_S02");
    BOOST_CHECK_EQUAL(chipid[3], "R01_S10");

    log_results(results_collector);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(fits_utils_tests)

BOOST_AUTO_TEST_CASE(test_FitsImage) {
    long nx(90);
    long ny(100);
    std::vector<float> output_image;
    for (long i(0); i < nx; i++) {
        for (long j(0); j < ny; j++) {
            output_image.push_back(RandomNumbers::RngFloat());
        }
    }
    std::string image_file("test_image.fits");
    fits_utils::FitsImage<float> output(output_image, nx, ny, "!" + image_file);
    output.write();
    output.close();

    long naxes[2];
    std::vector<float> input_image;
    fits_utils::FitsImage<float> input(input_image, naxes, image_file);
   
    BOOST_CHECK_EQUAL(naxes[0], nx);
    BOOST_CHECK_EQUAL(naxes[1], ny);
    for (long i(0); i < nx*ny; i++) {
        BOOST_CHECK_EQUAL(output_image[i], input_image[i]);
    }
    log_results(results_collector);
}

BOOST_AUTO_TEST_SUITE_END()
