#!/bin/zsh

R_SCRIPT=/mnt/f/Projects/NE_Ontario_Bird_Models/R/bash_run_brt.R
FILENAME=/mnt/f/Projects/NE_Ontario_Bird_Models/specieslist.txt

echo "Starting processing of $R_SCRIPT"

while IFS= read -r line; do
  echo "--- Starting BRT models for $line---"
  Rscript --vanilla $R_SCRIPT $line

done < $FILENAME