#!/bin/bash

## -------------------------------------------------------------
## This is an example integrate-script for Herwig-Parallel.
## Put it in the same folder as the infile you would like to run
## Herwig/Matchbox with and then execute 'herwig-parallel-read'.
## -------------------------------------------------------------

@HOSTNAME@

echo ''
echo '---------------------------------------------'
echo 'J O B   D A T A   A N D   H O S T   U S A G E'
echo '---------------------------------------------'
echo 'start time:           '$(date)
echo 'hostname:             '$(hostname)
echo 'operating system:     '$(uname -o)
echo 'kernel:               '$(uname -srv)
echo 'architecture:'
echo ' machine:             '$(uname -m)
echo ' processor:           '$(uname -p)
echo ' hardware platform:   '$(uname -i)
echo ''
echo 'memory & swap:'
free -g
echo ''
echo 'w:'
w
echo '---------------------------------------------'
echo ''

source ~/.bashrc

echo ''
echo '-----------------------------------------------------------'
echo 'H E R W I G / M A T C H B O X   I N T E G R A T E   S T E P'
echo '-----------------------------------------------------------'
echo ''
time Herwig++ integrate @RUNFILE@ --jobid=@JOBID@ @SETUPFILE@
