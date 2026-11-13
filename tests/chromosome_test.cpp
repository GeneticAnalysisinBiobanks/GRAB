#include "util/chromosome.hpp"
#include "tinytest.hpp"

TEST(numeric_chromosomes_ignore_zero_padding) {
    for (int chr = 1; chr <= 22; ++chr) {
        const auto expected = std::to_string(chr);
        std::string actual = "00" + expected;
        util::stripChromosomeZeroPadding(actual);
        CHECK(actual == expected);
    }
    std::string zero = "000";
    util::stripChromosomeZeroPadding(zero);
    CHECK(zero == "0");
}

TEST(chromosome_normalization_preserves_other_labels) {
    for (const char *label : {"", "0", "1", "22", "X", "XY", "PAR1", "0X", "chr01"}) {
        std::string actual = label;
        util::stripChromosomeZeroPadding(actual);
        CHECK(actual == label);
    }
}

TINYTEST_MAIN
