ECCO-V algorithm for application on CloudSat Cloud Profiling Radar (CPR) data.

The details of the algorithm are described in

Romatschke, U., Dixon, M., 2022: Vertically Resolved Convective/StratiformEcho Type Identification and ConvectivityRetrieval for Vertically Pointing Radars. EarthArXiv, https://doi.org/10.31223/X54S77.

CloudSat data are available at https://www.cloudsat.cira.colostate.edu/. We use radar reflectivity from the 2B-GEOPROF product (Marchand et al. 2008), and the convective/stratiform flag and freezing level from the 2C-PRECIP-COLUMN product (Haynes, 2018).

convStrat_cloudSat.m runs ECCO-V on a specific example.

convStrat_cloudSat_compare.m runs ECCO-V on a specific example and compares it with the CloudSat echo type flag.

convStrat_compare_stats.m runs a statistical comparison between ECCO-V and the CloudSat echo type flag.

All other files are functions.
