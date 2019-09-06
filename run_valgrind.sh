#!/bin/bash

# This script uses valgrind to check whether this project is working properly or
# not.  valgrind will perform a full leak analysis on the 'main' binary, and
# report any possible mismanagement of memory.

valgrind --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --verbose \
         --log-file=valgrind-out.txt \
         ./main_tests

echo -e "\n" 
echo -e "*******************************************************"
echo -e " Valgrind output was saved to file: 'valgrind-out.txt'"
echo -e "*******************************************************\n"

cat ./valgrind-out.txt
