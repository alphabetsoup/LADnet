#ifndef _GLOBAL_H_
#define _GLOBAL_H_
//#define ARMA_USE_ATLAS
//#define ARMA_USE_WRAPPER
#undef ARMA_USE_WRAPPER
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#ifdef _M_X64
#define ARMA_64BIT_WORD // 64bit only
#endif
#define USE_SCHEMA_DATA
//#define ARMA_ATLAS_INCLUDE_DIR "/usr/include/atlas/"
//#include "/usr/include/atlas/cblas.h"
//#include "/usr/include/atlas/clapack.h"
#include "armadillo"
#include <stdexcept>
#include <assert.h>
#endif
