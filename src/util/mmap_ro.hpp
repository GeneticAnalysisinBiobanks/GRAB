// mmap_ro.hpp — read-only memory-mapped file with RAII cleanup.
//
// Portable wrapper that picks POSIX mmap on Linux/macOS and
// CreateFileMapping/MapViewOfFile on Windows (MinGW), exposing a
// pointer-and-length view of the file for zero-copy parsing.
//
// The header intentionally hides <windows.h> from translation units
// that include it; platform handles are stored as opaque `void *`.

#pragma once

#include <cstddef>
#include <string>

namespace util {

class MmapReadOnly {
public:
    explicit MmapReadOnly(const std::string &path);
    ~MmapReadOnly();

    MmapReadOnly(const MmapReadOnly &)            = delete;
    MmapReadOnly &operator=(const MmapReadOnly &) = delete;

    const char *data() const noexcept { return static_cast<const char *>(addr_); }
    std::size_t size() const noexcept { return size_; }

private:
    void       *addr_  = nullptr;
    std::size_t size_  = 0;
    // Windows: HANDLE for the file and the mapping (both `void *`).
    // POSIX: unused, remain nullptr.
    void       *hFile_ = nullptr;
    void       *hMap_  = nullptr;
};

} // namespace util
