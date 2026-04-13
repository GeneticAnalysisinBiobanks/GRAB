#!/usr/bin/bash

find src -type f -name '*pp' -exec clang-format -i {} +
