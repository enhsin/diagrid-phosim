///
/// @package phosim
/// @file fits_utils.h
/// @brief Classes and functions to handle FITS data.
///
/// @brief Created by
/// @author J. Chiang (SLAC)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef fits_utils_h
#define fits_utils_h

#include "fitsio.h"

#include <cstdlib>

#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace fits_utils {

    class keyProperties {
    public:

    keyProperties(const std::string & dtype, const std::string & v,
                  const std::string & com)
        : datatype(dtype), value(v), comment(com) {
        }
    keyProperties(const std::string & dtype, double v,
                  const std::string & com)
        : datatype(dtype), value(""), comment(com) {
            std::ostringstream ss;
            ss << v;
            value=ss.str();
        }
    keyProperties(const std::string & dtype, float v,
                  const std::string & com)
        : datatype(dtype), value(""), comment(com) {
            std::ostringstream ss;
            ss << v;
            value=ss.str();
        }
    keyProperties(const std::string & dtype, long v,
                  const std::string & com)
        : datatype(dtype), value(""), comment(com) {
            std::ostringstream ss;
            ss << v;
            value=ss.str();
        }
    keyProperties(const std::string & dtype, int v,
                  const std::string & com)
        : datatype(dtype), value(""), comment(com) {
            std::ostringstream ss;
            ss << v;
            value=ss.str();
        }

        std::string datatype;
        std::string value;
        std::string comment;

    };

    /**
     * @class FitsImage
     * @brief A class to handle basic 2D FITS images.  Currently, will
     * only write FITS PrimaryHDU-only files, with minimal keyword
     * handling.  Template specializations only for double, float and
     * unsigned short images, so far.
     */

    template<class T>
        class FitsImage {
    public:

        FitsImage(const std::vector<T> & image, long nx, long ny,
                  const std::string & outfile);

        FitsImage(const std::vector<T> & image, long nx, long ny,
                  const std::string & outfile, 
                  const std::string & infile, 
                  float *dcrpix, long bitpix);

        FitsImage(std::vector<T> & image, long *naxes, const std::string & infile);

        ~FitsImage();

        void flip(int serialread, int parallelread);

        void write_keyword(const std::string & keyword, float value);

        void write_keyword(const std::string & keyword,
                           const keyProperties & keyProperties);

        void write();

        void close();

        /// @return Pixel values as a const vector<T>.
        const std::vector<T> & image() const {
            return m_image;
        }

        /// @return Number of pixels in the x-direction.
        long nx() const {
            return m_naxes[0];
        }
   
        /// @return Number of pixels in the y-direction.
        long ny() const {
            return m_naxes[1];
        }

    private:
        const std::vector<T> & m_image;
        std::string m_filename;
        fitsfile * m_fptr;
        long m_naxes[2];

        void create_image();
        void write_image();
        void read_image(std::vector<T> & image);
    };

    typedef std::map<std::string, float> keyword_map_t;
    typedef std::map<std::string, keyProperties> keyword_map;

    void write_fits_image(const std::vector<double> & in, long nx, long ny,
                          const std::string & outfile,
                          const keyword_map_t & keywords=keyword_map_t());

    void write_fits_image(const std::vector<float> & in, long nx, long ny,
                          const std::string & outfile,
                          const keyword_map_t & keywords=keyword_map_t());

    void write_fits_image_cphead(const std::vector<unsigned short> & in, 
                                 long nx, long ny,
                                 const std::string & outfile,
                                 const std::string & infile,
                                 float *dcrpix, int *fflag,
                                 const keyword_map & keywords=keyword_map());

    void read_fits_image(std::vector<float> & image, long *naxes, 
                         const std::string & infile);

    /// Implement naming convention used in C version of turb2d.cpp.
    std::string filename(const std::string & prefix, 
                         const std::string & basename);

} // namespace fits_utils

#endif // fits_utils_h
