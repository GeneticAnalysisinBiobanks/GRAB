// main.cpp — GRAB CLI entry point
//
// Canonical main.cpp: stays at src/ root, contains only the entry point.
// All parsing, help, and dispatch logic lives in src/cli/.

#include "cli/cli.hpp"

int main(
    int argc,
    char *argv[]
) {
    return cli::run(argc, argv);
}
