#pragma once

#include <cstddef>
#include <string>
#include <vector>

// Select written columns without changing the method's internal result layout.
struct OutputOptions {
    bool spaStatus = false;
    bool z = false;
};

enum class ResultColumnKind { Always, Z, SpaStatus };

struct ResultColumn {
    std::string name;
    ResultColumnKind kind = ResultColumnKind::Always;
};

using ResultColumns = std::vector<ResultColumn>;

// Build once per output schema; share the immutable selection across workers.
// The header and missing-value suffix include their leading tabs.
struct OutputLayout {
    std::vector<size_t> indices;
    std::string header;
    std::string naSuffix;

    OutputLayout(const ResultColumns &columns, const OutputOptions &options) {
        indices.reserve(columns.size());
        for (size_t i = 0; i < columns.size(); ++i) {
            const auto &column = columns[i];
            if (column.kind == ResultColumnKind::Z && !options.z) continue;
            if (column.kind == ResultColumnKind::SpaStatus && !options.spaStatus) continue;
            indices.push_back(i);
            header += '\t';
            header += column.name;
            naSuffix += "\tNA";
        }
    }
};
