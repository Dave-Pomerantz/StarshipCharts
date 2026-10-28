# StarshipCharts

Weekly sets of NOAA's electronic charts and gridded depths for Penobscot Bay to Narragansett Bay, for the Starship app's own devices.

Each set is a release. Its `set.tsv` lists every file of the set: the kind (chart or grid), NOAA's name and version, the release that holds the file, its size and SHA-256. A file unchanged since an earlier set stays in the earlier release, so every set can be read again later without copying unchanged files.

- **Charts** are NOAA's own zips, unaltered, as published at charts.noaa.gov, with NOAA's user agreement inside. A chart enters a set when NOAA's catalogue lists it as active and its bounds meet the region in `chartset/region.tsv`; it changes when NOAA's edition or update number changes, and a cancelled chart leaves the next set.
- **Grids** are converted from NOAA's S-102 Bathymetric Surface files (bucket `noaa-s102-pds`, edition 3.0.0). The conversion keeps every cell's depth, NOAA's stated uncertainty and the survey behind it, and adds coarser levels that keep each block's shallowest depth. The layout is described at the top of `chartset/grid_file.py`.

The grids are not NOAA's product and are not for navigation. NOAA publishes S-102 as a test product, to be used with the official charts, and the converted files are derived from it by this repository's code. NOAA's charts and grids are works of the United States government.

`.github/workflows/chart-set.yml` builds a set each Friday and on request.
