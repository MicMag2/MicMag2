#ifndef MAGNUM_CONFIG_H
#define MAGNUM_CONFIG_H

// Enable CVODE support?
/* #undef HAVE_CVODE */

// Enable CUDA support?
#define HAVE_CUDA
// Additionally, compile 64-bit CUDA code?
#define HAVE_CUDA_64

//#define ENABLE_CUDA_THREAD_SYNCHRONIZE 1
#define HAVE_FFTW_THREADS 1

#endif
