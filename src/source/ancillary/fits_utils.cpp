///
/// @package phosim
/// @file fits_utils.cpp
/// @brief Classes and functions to handle FITS data.
///
/// @brief Created by
/// @author J. Chiang (SLAC)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <stdexcept>

#include "fitsio2.h"

#include "ancillary/fits_utils.h"

namespace fits_utils {

    void write_fits_image(const std::vector<double> & image,
                          long nx, long ny,
                          const std::string & outfile,
                          const keyword_map_t & keywords) {
        FitsImage<double> foo(image, nx, ny, outfile);

        keyword_map_t::const_iterator it(keywords.begin());
        for ( ; it !=  keywords.end(); ++it) {
            foo.write_keyword(it->first, it->second);
        }
        foo.write();
    }

    void write_fits_image(const std::vector<float> & image, 
                          long nx, long ny,
                          const std::string & outfile,
                          const keyword_map_t & keywords) {
        FitsImage<float> foo(image, nx, ny, outfile);

        keyword_map_t::const_iterator it(keywords.begin());
        for ( ; it !=  keywords.end(); ++it) {
            foo.write_keyword(it->first, it->second);
        }
        foo.write();
    }

    void write_fits_image_cphead(const std::vector<unsigned short> & image,
                                 long nx, long ny,
                                 const std::string & outfile,
                                 const std::string & infile,
                                 float *crpix, int *fflag,
                                 const keyword_map & keywords) {

        FitsImage<unsigned short> foo(image, nx, ny, outfile, infile, crpix, 16);

        keyword_map::const_iterator it(keywords.begin());
        for ( ; it !=  keywords.end(); ++it) {
            foo.write_keyword(it->first, it->second);
        }
        // rotate and flip
        if (fflag[0] == 1) {
            foo.flip(fflag[1], fflag[2]);
        }
        foo.write();
    }

    void read_fits_image(std::vector<float> & image,
                         long *naxes,
                         const std::string & infile) {
        FitsImage<float> foo(image, naxes, infile);
    }

    template<class T>
    FitsImage<T>::FitsImage(const std::vector<T> & image, long nx, long ny,
                            const std::string & outfile)
        : m_image(image), m_filename(outfile), m_fptr(0) {
        int status(0);
        fits_create_file(&m_fptr, outfile.c_str(), &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::FitsImage: cfitsio error");
        }
        m_naxes[0] = nx;
        m_naxes[1] = ny;
        create_image();
    }

    template<class T>
    FitsImage<T>::FitsImage(const std::vector<T> & image, long nx, long ny,
                            const std::string & outfile,
                            const std::string & infile,
                            float *dcrpix, long bitpix) 
        : m_image(image), m_filename(outfile), m_fptr(0) {
        int status(0);
        fits_create_file(&m_fptr, outfile.c_str(), &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::FitsImage: cfitsio error");
        }

        fitsfile * m_fptr_in;
        fits_open_file(&m_fptr_in, infile.c_str(), READONLY, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::FitsImage: cfitsio error");
        }

        fits_copy_header(m_fptr_in, m_fptr, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::FitsImage: cfitsio error");
        }
   
        int keynum[] = {10, 13};
        float crpix[2];
        char name[4096], value[4096], comment[4096];
        for (int i(0); i < 2; i++) {
            fits_read_keyn(m_fptr, keynum[i], name, value, comment, &status);
            if (status != 0) {
                fits_report_error(stderr, status);
                throw std::runtime_error("FitsImage::FitsImage: cfitsio error");
            }
            crpix[i] = atof(value) + dcrpix[i];
        }

        fits_close_file(m_fptr_in, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::FitsImage: cfitsio error");
        }

        m_naxes[0] = nx;
        m_naxes[1] = ny;
        for (int i=0; i < 2; i++) {
            std::ostringstream ss;
            ss << "NAXIS" << (i+1);
            fits_update_key(m_fptr, TLONG, ss.str().c_str(), &m_naxes[i], NULL, &status);
            if (status != 0) {
                fits_report_error(stderr, status);
                throw std::runtime_error("FitsImage::FitsImage: cfitsio error");
            }
        }
        for (int i=0; i < 2; i++) {
            std::ostringstream ss;
            ss << "CRPIX" << (i+1);
            fits_update_key(m_fptr, TFLOAT, ss.str().c_str(), &crpix[i], NULL, &status);
            if (status != 0) {
                fits_report_error(stderr, status);
                throw std::runtime_error("FitsImage::FitsImage: cfitsio error");
            }
        }
        fits_update_key(m_fptr, TLONG, "BITPIX", &bitpix, NULL, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::FitsImage: cfitsio error");
        }
    }


    template <class T>
    void FitsImage<T>::flip(int serialread, int parallelread) {
        int status(0);
        // CRPIX1, CRVAL1, CRPIX2, CRVAL2, CD1_1, CD1_2, CD2_1, CD2_2
        int keynum[] = {10, 11, 13, 14, 15, 16, 17, 18}; 
        float keys[8], keys_orig[8];
        std::string keywords[8];
        char name[4096], value[4096], comment[4096];
   
        for (int i(0); i < 8; i++) {
            fits_read_keyn(m_fptr, keynum[i], name, value, comment, &status);
            if (status != 0) {
                fits_report_error(stderr, status);
                throw std::runtime_error("FitsImage::flip: cfitsio error");
            }
            keys_orig[i] = atof(value);
            std::ostringstream ss;
            ss << name;
            keywords[i] = ss.str();
        }
      
        /* x --> y */
        keys[0] = keys_orig[2];
        keys[1] = keys_orig[1];
        keys[2] = keys_orig[0];
        keys[3] = keys_orig[3];
        keys[4] = keys_orig[5];
        keys[5] = keys_orig[4];
        keys[6] = keys_orig[7];
        keys[7] = keys_orig[6];
      
        /* flip x */
        if (serialread == 1) {
            keys[0] = -keys[0] + m_naxes[0] - 1;
            keys[4] = -keys[4];
            keys[6] = -keys[6];
        }
        /* flip y */
        if (parallelread == 1) {
            keys[2] = -keys[2] + m_naxes[1] - 1;
            keys[5] = -keys[5];
            keys[7] = -keys[7];
        }
      
        for (int i(0); i < 8; i++) {
            fits_update_key(m_fptr, TFLOAT, keywords[i].c_str(), &keys[i],
                            NULL, &status);
            if (status != 0) {
                fits_report_error(stderr, status);
                throw std::runtime_error("FitsImage::flip: cfitsio error");
            }
        }
    }

    template<class T>
    FitsImage<T>::FitsImage(std::vector<T> & image, long *naxes, 
                            const std::string & infile)
        : m_image(image), m_filename(infile), m_fptr(0) {
        int status(0);
        int nfound(0);
        fits_open_file(&m_fptr, infile.c_str(), READONLY, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::FitsImage: cfitsio error");
        }
        fits_read_keys_lng(m_fptr, (char*)"NAXIS", 1, 2, m_naxes, &nfound, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::FitsImage: cfitsio error");
        }
        naxes[0] = m_naxes[0];
        naxes[1] = m_naxes[1];
        read_image(image); 
    }

    template<class T>
    FitsImage<T>::~FitsImage() {
        if (m_fptr) {
            close();
        }
    }

    template<class T>
    void FitsImage<T>::write_keyword(const std::string & keyword, float value) {
        int status(0);
        /* Convert keyword to char* since earlier versions of cfitsio cannot handle
           const char* */
        std::vector<char> v_tmp(keyword.begin(), keyword.end());
        v_tmp.push_back(0);
        fits_write_key(m_fptr, TFLOAT, &v_tmp[0], &value, "", &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::write_keyword: cfitsio error");
        }
    }

    template<class T>
    void FitsImage<T>::write_keyword(const std::string & keyword,
                                     const keyProperties & keyProperties) {
        int status(0);
        /* Convert keyword to char* since earlier versions of cfitsio cannot handle
           const char* */
        std::vector<char> v_tmp(keyword.begin(), keyword.end());
        v_tmp.push_back(0);
        if (keyProperties.datatype == "TLONG") {
            long value = std::atoi(keyProperties.value.c_str());
            fits_write_key(m_fptr, TLONG, &v_tmp[0], &value,
                           keyProperties.comment.c_str(), &status);
        } else if (keyProperties.datatype == "TINT") {
            int value = std::atoi(keyProperties.value.c_str());
            fits_write_key(m_fptr, TINT, &v_tmp[0], &value,
                           keyProperties.comment.c_str(), &status);
        } else if (keyProperties.datatype == "TFLOAT") {
            float value = std::atof(keyProperties.value.c_str());
            fits_write_key(m_fptr, TFLOAT, &v_tmp[0], &value,
                           keyProperties.comment.c_str(), &status);
        } else if (keyProperties.datatype == "TDOUBLE") {
            double value = std::atof(keyProperties.value.c_str());
            fits_write_key(m_fptr, TDOUBLE, &v_tmp[0], &value,
                           keyProperties.comment.c_str(), &status);
        } else if (keyProperties.datatype == "TSTRING") {
            fits_write_key(m_fptr, TSTRING, &v_tmp[0],
                           (char *)keyProperties.value.c_str(),
                           keyProperties.comment.c_str(), &status);
        } else {
            throw std::runtime_error("FitsImage::write_keyword: "
                                     "wrong keyword datatype");
        }
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::write_keyword: cfitsio error");
        }
    }

    template<class T>
    void FitsImage<T>::write() {
        write_image();
    }

    template<class T>
    void FitsImage<T>::close() {
        int status(0);
        fits_close_file(m_fptr, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::close: cfitsio error");
        }
        m_fptr = 0;
    }

    template<>
    void FitsImage<double>::create_image() {
        int status(0);
        fits_create_img(m_fptr, DOUBLE_IMG, 2, m_naxes, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::create_image: cfitsio error");
        }
    }

    template<>
    void FitsImage<float>::create_image() {
        int status(0);
        fits_create_img(m_fptr, FLOAT_IMG, 2, m_naxes, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::create_image: cfitsio error");
        }
    }

    template<>
    void FitsImage<unsigned short>::create_image() {
        int status(0);
        fits_create_img(m_fptr, SHORT_IMG, 2, m_naxes, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::create_image: cfitsio error");
        }
    }

    template<>
    void FitsImage<double>::write_image() {
        int status(0);
        long fele(1);
        long imgsize(m_naxes[0]*m_naxes[1]);
        fits_write_img(m_fptr, TDOUBLE, fele, imgsize, 
                       const_cast<double *>(&m_image[0]), &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::write_image: cfitsio error");
        }
    }

    template<>
    void FitsImage<float>::write_image() {
        int status(0);
        long fele(1);
        long imgsize(m_naxes[0]*m_naxes[1]);
        fits_write_img(m_fptr, TFLOAT, fele, imgsize, 
                       const_cast<float *>(&m_image[0]), &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::write_image: cfitsio error");
        }
    }

    template<>
    void FitsImage<unsigned short>::write_image() {
        int status(0);
        long fele(1);
        long imgsize(m_naxes[0]*m_naxes[1]);
        fits_write_img(m_fptr, TUSHORT, fele, imgsize,
                       const_cast<unsigned short *>(&m_image[0]), &status);
        if (status != 0 && status != 412) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::write_image: cfitsio error");
        }
    }

    template<>
    void FitsImage<float>::read_image(std::vector<float> & image) {
        int status(0);
        long fele(1);
        long imgsize(m_naxes[0]*m_naxes[1]);
        int anynull(0);
        float nullval(0);
        image.resize(imgsize, 0);
        fits_read_img(m_fptr, TFLOAT, fele, imgsize, &nullval, &image[0],
                      &anynull, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            throw std::runtime_error("FitsImage::read_image: cfitsio error");
        }
    }

    std::string filename(const std::string & prefix, 
                         const std::string & basename) {
        if (prefix.empty()) {
            return basename;
        }
        std::ostringstream my_filename;
        my_filename << "!" << prefix << "_" << basename;
        return my_filename.str();
    }

} // namespace fits_utils
