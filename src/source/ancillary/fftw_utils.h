///
/// @package phosim
/// @file fftw_utils.h
/// @brief Helper functions for FFTW.
///
/// @brief Created by
/// @author J. Chiang (SLAC)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef fftw_util_h
#define fftw_util_h

#include <stdexcept>

#include "fftw3.h"

namespace fftw_utils {

    class DoubleArray {

    public:

    DoubleArray(long nx, long ny) 
        : m_nx(nx), m_ny(ny), m_data(double_array(nx*ny)) {}

        ~DoubleArray() {
            fftw_free(m_data);
        }

    DoubleArray(const DoubleArray & other) 
        : m_nx(other.m_nx), m_ny(other.m_ny), m_data(double_array(m_nx*m_ny)) {
            copy_data(other);
        }

        DoubleArray & operator=(const DoubleArray & rhs) {
            if (this != &rhs) {
                m_nx = rhs.m_nx;
                m_ny = rhs.m_ny;
                m_data = double_array(m_nx*m_ny);
                copy_data(rhs);
            }
            return *this;
        }

        double * operator()() {
            return m_data;
        }

        double & operator[](long indx) {
            return m_data[indx];
        }

    private:

        long m_nx;
        long m_ny;
        double * m_data;

        inline double * double_array(int size) {
            double * out 
                = reinterpret_cast<double *>(fftw_malloc(sizeof(double)*size));
            if (out == NULL) {
                throw std::runtime_error("fftw_malloc double failed");
            }
            return out;
        }
        inline void copy_data(const DoubleArray & other) {
            for (long i(0); i < m_nx*m_ny; i++) {
                m_data[i] = other.m_data[i];
            }
        }
    };

    class ComplexArray {

    public:

    ComplexArray(long nx, long ny) 
        : m_nx(nx), m_ny(ny), m_data(complex_array(nx*ny)) {}

        ~ComplexArray() {
            fftw_free(m_data);
        }

    ComplexArray(const ComplexArray & other) 
        : m_nx(other.m_nx), m_ny(other.m_ny), m_data(complex_array(m_nx*m_ny)) {
            copy_data(other);
        }

        ComplexArray & operator=(const ComplexArray & rhs) {
            if (this != &rhs) {
                m_nx = rhs.m_nx;
                m_ny = rhs.m_ny;
                m_data = complex_array(m_nx*m_ny);
                copy_data(rhs);
            }
            return *this;
        }

        fftw_complex * operator()() {
            return m_data;
        }

        fftw_complex & operator[](long indx) {
            return m_data[indx];
        }

        void inverse_fft(ComplexArray & other) {
            fftw_plan pb(fftw_plan_dft_2d(m_nx, m_ny, m_data, other.m_data, 
                                          FFTW_BACKWARD, FFTW_ESTIMATE));
            fftw_execute(pb);
            fftw_destroy_plan(pb);
        }

    private:

        long m_nx;
        long m_ny;
        fftw_complex * m_data;

        inline fftw_complex * complex_array(int size) {
            fftw_complex * out = reinterpret_cast<fftw_complex *>
                (fftw_malloc(sizeof(fftw_complex)*size));
            if (out == NULL) {
                throw std::runtime_error("fftw_malloc fftw_complex failed");
            }
            return out;
        }
        inline void copy_data(const ComplexArray & other) {
            for (long i(0); i < m_nx*m_ny; i++) {
                m_data[i][0] = other.m_data[i][0];
                m_data[i][1] = other.m_data[i][1];
            }
        }
    };

} // namespace fftw_util

#endif // fftw_util_h
