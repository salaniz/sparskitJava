#!/usr/bin/env bash
# simple script to download and untar all files

mkdir -p files
cd files
grep -v '^#' ../files.txt | xargs -I % sh -c 'curl  http://www.cise.ufl.edu/research/sparse/MM/%.tar.gz | tar xvz'
