#ifndef GSWTEOS10_SHARED_EXPORTS_H
#define GSWTEOS10_SHARED_EXPORTS_H

#ifdef GSWTEOS10_STATIC
#  define GSWTEOS10_SHARED_EXPORTS
#else
#  ifndef GSWTEOS10_SHARED_EXPORTS
#    ifdef _WIN32
#       ifdef GSWTEOS_10_EXPORTS
            /* We are building this library */
#           define GSWTEOS10_SHARED_EXPORTS __declspec(dllexport)
#       else
            /* We are using this library */
#           define GSWTEOS10_SHARED_EXPORTS __declspec(dllimport)
#       endif
#    else
#       define GSWTEOS10_SHARED_EXPORTS
#    endif
#  endif
#endif

#endif /* GSWTEOS10_SHARED_EXPORTS_H */
