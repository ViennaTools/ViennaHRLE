#!/bin/bash

find `pwd` -iname "*.hpp" -o -iname "*.cpp" | xargs clang-format -i
