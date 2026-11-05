#pragma once

#include <cstdlib>
#include <filesystem>

namespace test_output {

inline std::filesystem::path directory(const char *name) {
    const char *runDirectory = std::getenv("GRAB_TEST_OUTPUT_DIR");
    const auto root = runDirectory && *runDirectory
        ? std::filesystem::path(runDirectory)
        : std::filesystem::path("tests_output");
    const auto path = root / name;
    std::filesystem::create_directories(path);
    return path;
}

} // namespace test_output
