#ifdef WIN32 //For microsoft's C++ compiler
    #include <climits>

    typedef __int64             int64;
    typedef unsigned __int64   uint64;
    typedef __int32             int32;
    typedef unsigned __int32   uint32;
    typedef __int16             int16;
    typedef unsigned __int16   uint16;
    typedef __int8              int8 ;
    typedef unsigned __int8    uint8 ;

    const uint32 uint32_MAX=_UI32_MAX;
#elif __GNUC__  //for g++

    #define __STDC_LIMIT_MACROS
    #include <stdint.h>

    typedef int64_t  int64;
    typedef uint64_t uint64;
    typedef int32_t  int32;
    typedef uint32_t uint32;
    typedef int16_t  int16;
    typedef uint16_t uint16;
    typedef int8_t   int8;
    typedef uint8_t  uint8;

     const uint32 uint32_MAX = 0xFFFFFFFFU;

#else // for all other compilers, standard-compilant includes
    #include <cstdint.h>

    typedef std:: int64_t  int64;
    typedef std::uint64_t uint64;
    typedef std:: int32_t  int32;
    typedef std::uint32_t uint32;
    typedef std:: int16_t  int16;
    typedef std::uint16_t uint16;
    typedef std:: int8_t   int8;
    typedef std::uint8_t  uint8;

    const uint32 uint32_MAX=std::uint32_MAX;
#endif
