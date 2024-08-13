#!/bin/bash

NUMTHREADS=$(grep -c ^processor /proc/cpuinfo)

for i in $(seq $((NUMTHREADS - 1))); do
    echo $i &
done
