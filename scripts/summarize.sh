#!/bin/bash
awk '
  FNR==1 {sec=""}
  /^=== sequential/      {sec="sequential"}
  /^=== mpi efficient/   {sec="efficient"}
  /^=== mpi inefficient/ {sec="inefficient"}
  /^Matrix dimension/    {n=$3; dets=""}
  /^Product seconds/     {prod=$3}
  /^Seconds/             {dets=dets "  " $2}
  /^Identity verified/   {ok=$3}
  /^Total seconds/ {printf "n=%-4s %-12s product %10s | detA detB detC:%s | identity %s\n", n, sec, prod, dets, ok}
  /CANCELLED/ {print FILENAME ": " $0}
' "$@"
