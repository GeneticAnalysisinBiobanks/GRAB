#pragma once

#include <string>

namespace util {

// Strip leading zeros only from entirely numeric chromosome IDs.
inline void stripChromosomeZeroPadding(std::string &chromosome) {
    if (chromosome.size() < 2 || chromosome[0] != '0' ||
        chromosome.find_first_not_of("0123456789") != std::string::npos)
        return;
    const auto first = chromosome.find_first_not_of('0');
    chromosome.erase(0, first == std::string::npos ? chromosome.size() - 1 : first);
}

} // namespace util
