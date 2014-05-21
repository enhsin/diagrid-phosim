///
/// @package phosim
/// @file atmosphere_creator.h
/// @brief atmosphere parameter distribution calculator
///
/// @brief Created by:
/// @author Mallory Young (Purdue)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
/// @author En-Hsin Peng (Purdue)
/// @author J. Chiang (SLAC)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef atmosphere_AtmosphereCreator_h
#define atmosphere_AtmosphereCreator_h

#include <string>
#include <vector>

namespace atmosphere {

    /// @class AtmosphereCreator

    class AtmosphereCreator {

    public:

        AtmosphereCreator(int numlevel, float groundlevel, const std::string & datadir, const std::string & instrdir);

        void run(float monthnum, float constrainseeing, const std::string & outputfilename, const std::vector<int> & cloudscreen, long seed);

        const std::vector<float> & altitudes() const {
            return m_altitudes;
        }

        const std::vector<float> & osests() const {
            return m_osests;
        }

    private:

        int m_numlevel;
        float m_groundlevel;
        std::string m_datadir;
        std::string m_instrdir;
        std::vector<float> m_osests;
        std::vector<float> m_altitudes;
        std::vector<float> m_jests;

        void ccalc();
        void magcalc(float monthnum, float altitude, float & magest, float & direst);
        float outerscale(float altitude) const;
    };

}
#endif // atmosphere_AtmosphereCreator_h
