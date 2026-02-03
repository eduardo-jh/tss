# tss
Time series from HDF files and multiple sources

```
Usage: tss [OPTIONS]
Extraction of time series from HDF files in multiple dirs.
Options:
  -f, --format FORMAT  Specify HDF format: hdf4 or hdf5 (default: hdf4)
  -i, --idir DIR       Input directories to scan HDF files (Landsat,MODIS)
  -o, --odir DIR       Output directory to write CSV time series files
  -d, --datasets LIST  List of datasets to read (e.g. "NDVI,EVI (calc)" or 0,1,2)
  -p, --pattern STR    Pattern to filter files to process (e.g. '.A2008')
  -n, --name STR       Basename for the CSV time series files
  -b, --bounds LIST    Row,col to start and row,col to finish extraction (e.g. 0,0,100,100)
  -s, --slack          Matches dates to -1 and +1 days of DOY
  -l, --lower VALUE    Exclude values less than VALUE from extraction (default: -10000)
  -u, --upper VALUE    Exclude values greater than VALUE from extraction (default: 10000)
  -v, --version        Show version information
  -h, --help           Show this help message
```

# Author

Copyright (C) 2025-2026 Eduardo Jimenez Hernandez <eduardojh@arizona.edu>.

# License

License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
here is NO WARRANTY, to the extent permitted by law.