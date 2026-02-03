#!/bin/bash
# ============================================================
# Description: Time series from Landsat and MODIS
# Script: run_tss.sh
# Purpose: Run the tss tool with configurable parameters.
# Author: Eduardo Jimenez Hernandez <eduardojh@arizona.edu>
# Date: 2026-02-02
# ============================================================

Program="./bin/tss"

OutputRes="750"
Sensor1="OLI08"
Sensor2="MOD"
Tile="h09v05"

CWD="/mnt/error/RESAMPLED/"
DirLandsat="${CWD}LANDSAT/${Tile}/"
DirMODIS="${CWD}MODIS/${Tile}/"
InputDirs="${DirLandsat},${DirMODIS}"
OutputDir="${CWD}TIMES_SERIES/"

Datasets="NDVI,EVI,EVI2"
Patterns=".A2013"
OutputFile="TS_${Sensor1}_${Sensor2}_${OutputRes}"

# === CONFIGURE BOUNDS ===
if [ "$Tile" == "h09v05" ]; then
  if [ "$OutputRes" == "250" ]; then
      Bounds="0,0,649,771"
  elif [ "$OutputRes" == "500" ]; then
      Bounds="0,0,324,385"
  elif [ "$OutputRes" == "750" ]; then
      Bounds="0,0,216,257"
  elif [ "$OutputRes" == "1500" ]; then
      Bounds="0,0,108,128"
  else
      echo "Error: Either unsupported resolution '$OutputRes'. Please use 250, 500, 750, or 1500."
      exit 1
  fi
elif [ "$Tile" == "h11v08" ]; then
  if [ "$OutputRes" == "250" ]; then
      Bounds="0,0,651,646"
  elif [ "$OutputRes" == "500" ]; then
      Bounds="0,0,325,323"
  elif [ "$OutputRes" == "750" ]; then
      Bounds="0,0,216,215"
  elif [ "$OutputRes" == "1500" ]; then
      Bounds="0,0,108,107"
  else
      echo "Error: Either unsupported resolution '$OutputRes'. Please use 250, 500, 750, or 1500."
      exit 1
  fi
elif [ "$Tile" == "h11v11" ]; then
  if [ "$OutputRes" == "250" ]; then
      Bounds="0,0,640,704"
  elif [ "$OutputRes" == "500" ]; then
      Bounds="0,0,320,351"
  elif [ "$OutputRes" == "750" ]; then
      Bounds="0,0,213,234"
  elif [ "$OutputRes" == "1500" ]; then
      Bounds="0,0,106,116"
  else
      echo "Error: Either unsupported resolution '$OutputRes'. Please use 250, 500, 750, or 1500."
      exit 1
  fi
elif [ "$Tile" == "h17v07" ]; then
  if [ "$OutputRes" == "250" ]; then
      Bounds="0,0,648,666"
  elif [ "$OutputRes" == "500" ]; then
      Bounds="0,0,323,333"
  elif [ "$OutputRes" == "750" ]; then
      Bounds="0,0,215,222"
  elif [ "$OutputRes" == "1500" ]; then
      Bounds="0,0,107,111"
  else
      echo "Error: Either unsupported resolution '$OutputRes'. Please use 250, 500, 750, or 1500."
      exit 1
  fi
else
  echo "       Or unsupported tile '$Tile'. Please use h09v05, h11v08, h11v11, or h17v07."
  exit 1
fi

# === DISPLAY SETTINGS ===
echo "------------------------------------------------------------"
echo " Running tss with the following parameters:"
echo " Program:       $Program"
echo " Input dirs:    $InputDirs"
echo " Output dir:    $OutputDir"
echo " Datasets:      $Datasets"
echo " Pattern:       $Patterns"
echo " Output file:   $OutputFile"
echo " Bounds:        $Bounds"
echo "------------------------------------------------------------"
echo

# === EXECUTION ===
"$Program" \
  --idir "$InputDirs" \
  --odir "$OutputDir" \
  --datasets "$Datasets" \
  --pattern "$Patterns" \
  --name "$OutputFile" \
  --bounds "$Bounds"
  # --debug

# === CHECK EXIT STATUS ===
exit_code=$?
if [ $exit_code -eq 0 ]; then
  echo "tss completed SUCCESSFULLY!"
else
  echo "tss encountered an error (exit code $exit_code)."
  exit $exit_code
fi
