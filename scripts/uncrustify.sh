#!/usr/bin/bash

find src -type f -name '*pp' -exec uncrustify -c uncrustify.cfg --replace --no-backup {} \;
