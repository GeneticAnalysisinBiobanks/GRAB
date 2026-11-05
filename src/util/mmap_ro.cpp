// mmap_ro.cpp — read-only memory-mapped file, POSIX / Windows backends.

#include "util/mmap_ro.hpp"

#include <stdexcept>

#if defined(_WIN32)
#  include <windows.h>
#else
#  include <fcntl.h>
#  include <sys/mman.h>
#  include <sys/stat.h>
#  include <unistd.h>
   // MAP_POPULATE pre-faults pages on Linux; macOS / BSD have no
   // equivalent so we fall back to plain MAP_PRIVATE.
#  ifndef MAP_POPULATE
#    define MAP_POPULATE 0
#  endif
#endif

namespace util {

MmapReadOnly::MmapReadOnly(const std::string &path) {
#if defined(_WIN32)
    HANDLE hFile = ::CreateFileA(path.c_str(), GENERIC_READ, FILE_SHARE_READ,
                                 nullptr, OPEN_EXISTING,
                                 FILE_ATTRIBUTE_NORMAL, nullptr);
    if (hFile == INVALID_HANDLE_VALUE)
        throw std::runtime_error("Cannot open file: " + path);

    LARGE_INTEGER li;
    if (!::GetFileSizeEx(hFile, &li)) {
        ::CloseHandle(hFile);
        throw std::runtime_error("Cannot stat file: " + path);
    }
    const std::size_t fileSize = static_cast<std::size_t>(li.QuadPart);
    if (fileSize == 0) {
        ::CloseHandle(hFile);
        throw std::runtime_error("Empty file: " + path);
    }

    HANDLE hMap = ::CreateFileMappingA(hFile, nullptr, PAGE_READONLY, 0, 0,
                                       nullptr);
    if (!hMap) {
        ::CloseHandle(hFile);
        throw std::runtime_error("Cannot create file mapping: " + path);
    }

    void *mapped = ::MapViewOfFile(hMap, FILE_MAP_READ, 0, 0, 0);
    if (!mapped) {
        ::CloseHandle(hMap);
        ::CloseHandle(hFile);
        throw std::runtime_error("Cannot map view of file: " + path);
    }

    addr_  = mapped;
    size_  = fileSize;
    hFile_ = hFile;
    hMap_  = hMap;
#else
    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0)
        throw std::runtime_error("Cannot open file: " + path);

    struct stat st;
    if (::fstat(fd, &st) != 0) {
        ::close(fd);
        throw std::runtime_error("Cannot stat file: " + path);
    }
    const std::size_t fileSize = static_cast<std::size_t>(st.st_size);
    if (fileSize == 0) {
        ::close(fd);
        throw std::runtime_error("Empty file: " + path);
    }

    void *mapped = ::mmap(nullptr, fileSize, PROT_READ,
                          MAP_PRIVATE | MAP_POPULATE, fd, 0);
    ::close(fd);
    if (mapped == MAP_FAILED)
        throw std::runtime_error("Cannot mmap file: " + path);

    ::madvise(mapped, fileSize, MADV_SEQUENTIAL);

    addr_ = mapped;
    size_ = fileSize;
#endif
}

MmapReadOnly::~MmapReadOnly() {
#if defined(_WIN32)
    if (addr_)  ::UnmapViewOfFile(addr_);
    if (hMap_)  ::CloseHandle(static_cast<HANDLE>(hMap_));
    if (hFile_) ::CloseHandle(static_cast<HANDLE>(hFile_));
#else
    if (addr_) ::munmap(addr_, size_);
#endif
}

} // namespace util
