/*
  This is memory.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "memory.h"
#include <limits.h>
#include <cstring> // for |memcpy|, |memset|
#include <cassert>

#include "error.h"

namespace memory {
  using namespace error;
};

namespace {
  using namespace memory;
  const Ulong MEMORY_MAX = ULONG_MAX;
  const Ulong ABYTES = sizeof(Align);
  const Ulong ARENA_BITS = 16;
};

/****************************************************************************

  This module contains the memory allocator for this program. Since we are
  doing a lot of dynamic memory allocation, and use many variable-sized
  structures, both small and large, efficient memory allocation seems to
  be a crucial performance issue.

  We try to take an approach which is as simple as possible, for the sake
  of efficiency, while trying to be not too wasteful. First of all, memory
  is gotten from the system, via calloc, in chunks of 2^N.a bytes, where
  a = ABYTES will ensure proper alignment for all the types used
  in this program (here we make the crucial assumption that all builtin
  types have a size which is a power of two.) The only exception is when
  a memory request does not fit in such a chunk; then we allocate a block
  of size 2^m.a for large enough m. These chunks, once gotten,
  are never given back; however, we try to reuse them insofar as possible.

  Only three basic operations are provided : alloc, which returns a pointer
  to a new memory block; grow, which resizes a block to a larger block;
  and free, which frees a block. All blocks handed out in this way have
  sizes of the form 2^m.a, 0 <= m < BITS(Ulong). On average this should
  lead to a memory efficiency of 75%, which is not worse than what I was
  doing earlier, anyway.

 ****************************************************************************/

namespace memory {

Arena& arena()

{
  static Arena a(ARENA_BITS);
  return a;
}

};

/****************************************************************************

        Chapter I -- The Arena class.

 ****************************************************************************/

namespace memory {

Arena::Arena(Ulong bsBits)
  : d_bsBits(bsBits), d_count(0)
{
  for (unsigned i=0; i<BITS(Ulong); ++i)
  {
    d_list[i] = nullptr;
    d_used[i] = 0;
    d_allocated[i] = 0;
  }
}

Arena::~Arena()

/*
  Nothing to do here! All memory is allocated in fixed-size arrays.
*/

{}

/*
  Provide a possibly new block of size |2^b| (this is only called when
  none are directly available, i.e., when |d_list[b]==nullptr|). We look
  first whether a free block of larger size is available; if so, we
  split that one up, leaving two blocks of size |2^b| in |d_list[b]|
  If no such free block is found, we request through |calloc| a block of
  size |2^b| if that is large enough, or otherwise one of the minimal size
  |2^d_bsBits|, wich is then split up as before.
  (All sizes are in units of |ABYTES|, the effective word size for us.)

  Note that this is the only place where we can run out of memory.
  In that case, the error OUT_OF_MEMORY is set, and the corresponding
  error function run. In most cases, this will terminate the program.

  There is no return value, but if no |OUT_OF_MEMORY| error occurs,
  |d_list[b]| will point to a list of at least one block of size |2^b|.
  The caller must test that |ERRNO==0| before using the provided block.

  NOTE : as this function will be heavily used, it should be rewritten
  using a bitmap of available blocks.
*/
void Arena::newBlock(unsigned b)
{
  for (unsigned j = b+1; j < BITS(Ulong); ++j) {
    if (d_list[j]!=nullptr) /* split this block up */
      {
	Align *ptr = reinterpret_cast<Align*> (d_list[j]);
	d_list[j] = d_list[j]->next;
	d_allocated[j]--;
	for (unsigned i = b; i < j; ++i) {
	  assert(d_list[i]==nullptr); // tested before we came here
	  d_list[i] = reinterpret_cast<MemBlock*> (ptr + (1L<<i));
	  assert(d_list[i]->next==nullptr); // since free blocks are 0-filled
	  d_allocated[i]++;
	}
	d_list[b]->next = reinterpret_cast<MemBlock*> (ptr); // link @0 after @1
	d_list[b]->next->next = nullptr; // this one element was not 0-filled
	d_allocated[b]++;
	return;
      }
  }

  /* If we get here we need more memory from the system.
     When |Error| is called, it may or may not invoke |exit(0)|, depending
     on |CATCH_MEMORY_OVERFLOW|; the caller must therefore be aware and
     always test that |ERRNO==0| before using |d_list[b]|
 */

  if (b >= d_bsBits) { /* get block directly */
    if (d_count > MEMORY_MAX-(1L<<b)) // total memory would exceed 2^64 words
    {
      Error(OUT_OF_MEMORY); // actually out of compiled-in limit
      return;
    }
    d_list[b] = static_cast<MemBlock *> (calloc(1L<<b,ABYTES));
    if (d_list[b] == nullptr) // the system failed to honor the request
    {
      Error(OUT_OF_MEMORY);
      return;
    }
    d_count += 1L<<b; // record the number of words bought from system
    d_allocated[b]++;
    return;
  }

  // now we are asking for a very small block, and it needs fresh allocation
  if (d_count > MEMORY_MAX-(1L<<d_bsBits)) {
    Error(OUT_OF_MEMORY);
    return;
  }

  Align *ptr = static_cast<Align *> (calloc(1L<<d_bsBits,ABYTES));
  if (ptr == 0) {
    Error(OUT_OF_MEMORY);
    return;
  }

  d_count += 1L<<d_bsBits;

  for (unsigned j = b; j < d_bsBits; ++j) {
    d_list[j] = reinterpret_cast<MemBlock *> (ptr + (1L<<j));
    d_allocated[j]++;
  }

  d_list[b]->next = reinterpret_cast<MemBlock *> (ptr);
  d_allocated[b]++;

  return;
}

/*
  Return a pointer to a block of |2^m.ABYTES| bytes, where |m| is the
  smallest integer such that |2^m.ABYTES >= n|.

  It is assumed that |ABYTES| is a power of 2.

  The memory is zero-initialized.
*/
void* Arena::alloc(size_t n)
{
  using namespace constants; // for |lastBit|
  if (n == 0)
    return nullptr; // silly request, silly reply

  /* compute size of block */

  unsigned b = 0; // default to |2^0==1| word
  if (n > ABYTES) // but if more is asked for, request chunck of |2^b| words
    b = lastBit(n-1)-lastBit(ABYTES)+1; // with |b| large enough
  // the conditional above ensures we never set |b<0|, or call |lastBit(0)|

  if (d_list[b] == nullptr) { /* need to make a new block */
    newBlock(b);
    if (ERRNO!=0)
      return nullptr; // kick the error bucket down the road
  }

  /* take block off from list */

  MemBlock *block = d_list[b];
  assert(block != nullptr); // because we called |newBlock| and tested |ERRNO|

  d_list[b] = d_list[b]->next; // pop off now used first block from the list
  block->next = nullptr; // and clear up that pointer so whole block is zero
  d_used[b]++;

  return static_cast<void*> (block);
}

/*
  Return the size of the actual memory allocation provided on a request
  of |n| nodes of size |m|, in units of |m| (so result is at least |n|)
*/
Ulong Arena::allocSize(Ulong n, Ulong m)
{
  return byteSize(n,m)/m;
}

/*
  Return the actual number of bytes of the memory allocation (as opposed
  to allocSize, which rounds the allocation to the largest multiple of m.)
*/
Ulong Arena::byteSize(Ulong n, Ulong m)
{
  using namespace constants; // for |lastBit| and friends
  if (n == 0)
    return 0;
  if (n*m <= ABYTES)
    return ABYTES;
  return (1 << (lastBit(n*m-1)-lastBit(ABYTES)+1))*ABYTES;
}

/*
  Resize |ptr| to size new_size. This involves getting the larger block,
  copying the contents of ptr to it, and freeing ptr; we never try to
  fuse smaller adjacent blocks together.

  Returns 0 and sets the error MEMORY_WARNING in case of overflow, if
  CATCH_MEMORY_OVERFLOW is set.

  NOTE : equivalent to alloc if old_size = 0.
*/

void *memory::Arena::realloc(void *ptr, size_t old_size, size_t new_size)
{
  assert(old_size<new_size); // this is only for growing, not shrinking
  void *new_ptr = alloc(new_size);
  if (new_ptr==nullptr) /* memory overflow occurred */
    return nullptr; // kick the error bucket down the road
  if (old_size>0) // move contents of smaller old block into new block
  {
    std::memcpy(new_ptr,ptr,old_size);
    free(ptr,old_size);
  }

  return new_ptr;
}

/*
  Return the memory block allocated to |ptr| to the free list. In order to
  know to which list the pointer should be pushed (it will in fact be
  prepended), we need to pass the size |n| to which |ptr| was allocated.
  Actually |n| can be the size that was asked for, actual size is power of 2
*/
void Arena::free(void *ptr, size_t n)
{
  using namespace constants; // for |lastBit|
  if (ptr == 0)
    return;
  if (n == 0)
    return;

  unsigned b = 0;
  if (n > ABYTES)
    b = lastBit(n-1)-lastBit(ABYTES)+1;

  memset(ptr,0,(1L<<b)*ABYTES); // erase full contents of allocted block
  MemBlock *block = static_cast<MemBlock*>(ptr); // prepare to link in free list
  block->next = d_list[b]; // push to list (uses one pointer insize the block)
  d_list[b] = block; // set the list to start with this freed block
  d_used[b]--;
}

// Print information about given out memory, to |file|
void Arena::print(FILE *file) const
{
  fprintf(file,"%-10s%10s/%-10s\n","size : 2^","used","allocated"); // header

  Ulong used_count = 0;

  for (unsigned j = 0; j < BITS(Ulong); ++j) {
    fprintf(file,"%3u%7s%10lu/%-10lu\n",j,"",d_used[j],d_allocated[j]);
    used_count += (1L<<j)*d_used[j];
  }

  fprintf(file,"\n");
  fprintf(file,"total : %10lu/%-10lu %lu-byte units used/allocated\n",
	  used_count,static_cast<Ulong>(d_count),ABYTES);
}

};


namespace memory {

void pause()

{
  ;
}

};
