// logging.hpp — Thread-safe timestamped stderr logger
//
// Usage:  infoMsg("processed %d markers", n);
// Output: [INFO] 2026-03-28 20:05:56 processed 236 markers
#pragma once

#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <mutex>

// Function-local static ensures the mutex is initialized exactly once
// even across translation units (C++17 magic statics are thread-safe).
inline void infoMsg(
    const char *fmt,
    ...
) {
    static std::mutex s_logMutex;
    char buf[512];
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    std::lock_guard<std::mutex> lk(s_logMutex);
    auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    char ts[20];
    std::strftime(ts, sizeof(ts), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    std::fprintf(stderr, "[INFO] %s %s\n", ts, buf);
}
