// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#pragma once

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#if defined(__APPLE__)
#else
#include <malloc.h>
#endif
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "bridge.h"
#include "pbbs_strings.h"

namespace gbbs_io {

template <class E>
struct pairFirstCmp {
  bool operator()(std::pair<uintE, E> a, std::pair<uintE, E> b) {
    return a.first < b.first;
  }
};

template <class E>
struct getFirst {
  uintE operator()(std::pair<uintE, E> a) { return a.first; }
};

// returns a pointer and a length
inline std::pair<char*, size_t> mmapStringFromFile(const char* filename) {
  struct stat sb;
  int fd = open(filename, O_RDONLY);
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1) {
    perror("fstat");
    exit(-1);
  }
  if (!S_ISREG(sb.st_mode)) {
    perror("not a file\n");
    exit(-1);
  }
  char* p =
      static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (p == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  if (close(fd) == -1) {
    perror("close");
    exit(-1);
  }
  size_t n = sb.st_size;
  return std::make_pair(p, n);
}

void unmmap(char* bytes, size_t bytes_size) {
  if (bytes) {
    if (munmap(bytes, bytes_size) == -1) {
      perror("munmap");
      exit(-1);
    }
  }
}

inline sequence<char> readStringFromFile(char* fileName) {
  std::ifstream file(fileName, std::ios::in | std::ios::binary | std::ios::ate);
  if (!file.is_open()) {
    debug(std::cout << "Unable to open file: " << fileName << "\n";);
    abort();
  }
  uint64_t end = file.tellg();
  file.seekg(0, std::ios::beg);
  uint64_t n = end - file.tellg();
  auto bytes = sequence<char>(n); // n+1?
  file.read(bytes.begin(), n);
  file.close();
  return bytes;
}

std::tuple<char*, size_t> read_o_direct(char* fname) {
  /* read using O_DIRECT, which bypasses caches. */
  int fd;
#if defined(__APPLE__)
  if ((fd = open(fname, O_RDONLY)) != -1) {
#else
  if ((fd = open(fname, O_RDONLY | O_DIRECT)) != -1) {
#endif
    debug(std::cout << "input opened!"
              << "\n";);
  } else {
    std::cout << "can't open input file!";
  }
  //    posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);

  size_t fsize = lseek(fd, 0, SEEK_END);
  lseek(fd, 0, 0);

  /* allocate properly memaligned buffer for bytes */
#if defined(__APPLE__)
  char* bytes = NULL;
  posix_memalign((void**)&bytes, 4096 * 2, fsize + 4096);
#else
  char* bytes = (char*)memalign(4096 * 2, fsize + 4096);
#endif
  debug(std::cout << "fsize = " << fsize << "\n";);

  size_t sz = 0;

  size_t pgsize = getpagesize();
  debug(std::cout << "pgsize = " << pgsize << "\n";);

  size_t read_size = 1024 * 1024 * 1024;
  if (sz + read_size > fsize) {
    size_t k = std::ceil((fsize - sz) / pgsize);
    read_size = std::max(k * pgsize, pgsize);
    debug(std::cout << "set read size to: " << read_size << " " << (fsize - sz)
              << " bytes left"
              << "\n";);
  }

  while (sz + read_size < fsize) {
    void* buf = bytes + sz;
    debug(std::cout << "reading: " << read_size << "\n";);
    sz += read(fd, buf, read_size);
    debug(std::cout << "read: " << sz << " bytes"
              << "\n";);
    if (sz + read_size > fsize) {
      size_t k = std::ceil((fsize - sz) / pgsize);
      read_size = std::max(k * pgsize, pgsize);
      debug(std::cout << "set read size to: " << read_size << " " << (fsize - sz)
                << " bytes left"
                << "\n";);
    }
  }
  if (sz < fsize) {
    debug(std::cout << "last read: rem = " << (fsize - sz) << "\n";);
    void* buf = bytes + sz;
    sz += read(fd, buf, pgsize);
    debug(std::cout << "read " << sz << " bytes "
              << "\n";);
  }
  close(fd);
  return std::make_tuple(bytes, fsize);
}

} // namespace gbbs_io
