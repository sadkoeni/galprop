#ifndef _mlc_Malloc_h_
#define _mlc_Malloc_h_

#include <stdlib.h>

#define ALIGN_BYTES 64L
#define FLOATS_IN_ALIGN_BYTES 16L
#define ALIGNED     __attribute__((aligned(ALIGN_BYTES)))

#define ALIGN_BYTES_4K 4096L
#define FLOATS_IN_ALIGN_BYTES_4K 1024L
#define ALIGNED_4K     __attribute__((aligned(ALIGN_BYTES_4K)))

namespace mlc {

  bool Aligned(const void* ptr);

  void* MallocAligned(const size_t bytes);

  void FreeAligned(void* ptr);

  bool Aligned4K(const void* ptr);

  void* MallocAligned4K(const size_t bytes);

  void FreeAligned4K(void* ptr);

};

#endif
