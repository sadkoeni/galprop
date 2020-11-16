#include <Malloc.h>

namespace mlc {

  bool Aligned(const void* ptr)
  {
    // Returns true if the pointer is pointing to a location aligned to alignment-byte boundary, or false otherwise
    return (((size_t)ptr)%ALIGN_BYTES == 0);
  }

  void* MallocAligned(const size_t bytes) {
    void* ptr;
    const size_t ptrSize = sizeof(ptr);
    // Returns a pointer to a location aligned on an `alignment'-byte boundary
    // with at least `bytes' bytes allocated.
    // `holder' is a pointer that should be free()ed to release memory.
    ptr = malloc(bytes+ALIGN_BYTES-1+ptrSize);
    if (ptr==NULL) return NULL;

    const void* offsPtr = (void*)((char*)ptr + ptrSize); // Offset by ptrSize to make room for storing the original pointer

    const size_t ptrAddr = (size_t)(offsPtr); // Numerical value  for futher operations
    const size_t ptrMisalignment = ALIGN_BYTES - (ptrAddr % ALIGN_BYTES); // How far away from the nearest aligned address on the right
    const size_t ptrMoveBy = ptrMisalignment % ALIGN_BYTES; // If ptrMisalignment==ALIGN_BYTES, then ptrMoveBy=0.
    
    void* alignedAddress = (void*) ( (char*)offsPtr + ptrMoveBy ); // The nearest aligned address
    void** saveAddress = (void**) ( (char*)(alignedAddress) - ptrSize ); // Where the original pointer is to be stored
    *saveAddress = (void*)ptr; // Storing the original pointer for deletion
    return alignedAddress;
  }

  void FreeAligned(void* ptr) {
     if (ptr == NULL) return;
    void** savePtrLoc = (void**)((char*)ptr - sizeof(void*)); // Memory address where the original pointer is stored
    void* originalPtr = *savePtrLoc; // Original pointer returned by malloc()
    free(originalPtr);
  }

  bool Aligned4K(const void* ptr)
  {
    // Returns true if the pointer is pointing to a location aligned to alignment-byte boundary, or false otherwise
    return (((size_t)ptr)%ALIGN_BYTES_4K == 0);
  }

  void* MallocAligned4K(const size_t bytes) {
    void* ptr;
    const size_t ptrSize = sizeof(ptr);
    // Returns a pointer to a location aligned on an `alignment'-byte boundary
    // with at least `bytes' bytes allocated.
    // `holder' is a pointer that should be free()ed to release memory.
    ptr = malloc(bytes+ALIGN_BYTES_4K-1+ptrSize);
    if (ptr==NULL) return NULL;

    const void* offsPtr = (void*)((char*)ptr + ptrSize); // Offset by ptrSize to make room for storing the original pointer

    const size_t ptrAddr = (size_t)(offsPtr); // Numerical value  for futher operations
    const size_t ptrMisalignment = ALIGN_BYTES_4K - (ptrAddr % ALIGN_BYTES_4K); // How far away from the nearest aligned address on the right
    const size_t ptrMoveBy = ptrMisalignment % ALIGN_BYTES_4K; // If ptrMisalignment==ALIGN_BYTES, then ptrMoveBy=0.
    
    void* alignedAddress = (void*) ( (char*)offsPtr + ptrMoveBy ); // The nearest aligned address
    void** saveAddress = (void**) ( (char*)(alignedAddress) - ptrSize ); // Where the original pointer is to be stored
    *saveAddress = (void*)ptr; // Storing the original pointer for deletion
    return alignedAddress;
  }

  void FreeAligned4K(void* ptr) {
     if (ptr == NULL) return;
    void** savePtrLoc = (void**)((char*)ptr - sizeof(void*)); // Memory address where the original pointer is stored
    void* originalPtr = *savePtrLoc; // Original pointer returned by malloc()
    free(originalPtr);
  }


};
