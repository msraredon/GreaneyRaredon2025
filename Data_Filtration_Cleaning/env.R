#### ENVIRONMENT
#### Code for publication defining conda environment & R session info
#### Used in filter.R, clean.R, explore.R, focus.R, figure.R




#### Setting up conda environment
conda create -n englung -c conda-forge \
                        -c bioconda \
                        r-seurat=4* \
                        r-ggplot2 \
                        r-viridis \
                        r-rcolorbrewer \
                        r-velocyto.r \
                        bioconductor-complexheatmap




#### conda list
# packages in environment at .../englung:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       1_gnu    conda-forge
_r-mutex                  1.0.1               anacondar_1    conda-forge
binutils_impl_linux-64    2.36.1               h193b22a_2    conda-forge
binutils_linux-64         2.36                 hf3e587d_9    conda-forge
bioconductor-batchelor    1.10.0            r41h619a076_1    bioconda
bioconductor-beachmat     2.10.0            r41h619a076_1    bioconda
bioconductor-biobase      2.54.0            r41h5c21468_1    bioconda
bioconductor-biocgenerics 0.40.0            r41hdfd78af_0    bioconda
bioconductor-biocneighbors 1.12.0            r41h619a076_1    bioconda
bioconductor-biocparallel 1.28.3            r41h619a076_0    bioconda
bioconductor-biocsingular 1.10.0            r41h619a076_1    bioconda
bioconductor-biocviews    1.62.0            r41hdfd78af_0    bioconda
bioconductor-complexheatmap 2.10.0            r41hdfd78af_0    bioconda
bioconductor-delayedarray 0.20.0            r41h5c21468_1    bioconda
bioconductor-delayedmatrixstats 1.16.0            r41hdfd78af_0    bioconda
bioconductor-genomeinfodb 1.30.0            r41hdfd78af_0    bioconda
bioconductor-genomeinfodbdata 1.2.7             r41hdfd78af_1    bioconda
bioconductor-genomicranges 1.46.1            r41h5c21468_0    bioconda
bioconductor-graph        1.72.0            r41h5c21468_1    bioconda
bioconductor-hsmmsinglecell 1.14.0            r41hdfd78af_0    bioconda
bioconductor-iranges      2.28.0            r41h5c21468_1    bioconda
bioconductor-limma        3.50.1            r41h5c21468_0    bioconda
bioconductor-matrixgenerics 1.6.0             r41hdfd78af_0    bioconda
bioconductor-pcamethods   1.86.0            r41h619a076_1    bioconda
bioconductor-rbgl         1.70.0            r41h619a076_1    bioconda
bioconductor-residualmatrix 1.4.0             r41hdfd78af_0    bioconda
bioconductor-s4vectors    0.32.3            r41h5c21468_0    bioconda
bioconductor-scaledmatrix 1.2.0             r41hdfd78af_0    bioconda
bioconductor-scuttle      1.4.0             r41h619a076_1    bioconda
bioconductor-singlecellexperiment 1.16.0            r41hdfd78af_0    bioconda
bioconductor-sparsematrixstats 1.6.0             r41h619a076_1    bioconda
bioconductor-summarizedexperiment 1.24.0            r41hdfd78af_0    bioconda
bioconductor-xvector      0.34.0            r41h5c21468_1    bioconda
bioconductor-zlibbioc     1.40.0            r41h5c21468_1    bioconda
blosc                     1.21.0               h9c3ff4c_0    conda-forge
boost-cpp                 1.74.0               h6cacc03_7    conda-forge
bwidget                   1.9.14               ha770c72_1    conda-forge
bzip2                     1.0.8                h7f98852_4    conda-forge
c-ares                    1.18.1               h7f98852_0    conda-forge
ca-certificates           2022.9.24            ha878542_0    conda-forge
cairo                     1.16.0            ha12eb4b_1010    conda-forge
cfitsio                   4.0.0                h9a35b8e_0    conda-forge
curl                      7.82.0               h2283fc2_0    conda-forge
expat                     2.4.8                h27087fc_0    conda-forge
font-ttf-dejavu-sans-mono 2.37                 hab24e00_0    conda-forge
font-ttf-inconsolata      3.000                h77eed37_0    conda-forge
font-ttf-source-code-pro  2.038                h77eed37_0    conda-forge
font-ttf-ubuntu           0.83                 hab24e00_0    conda-forge
fontconfig                2.13.96              h8e229c2_2    conda-forge
fonts-conda-ecosystem     1                             0    conda-forge
fonts-conda-forge         1                             0    conda-forge
freetype                  2.10.4               h0708190_1    conda-forge
freexl                    1.0.6                h7f98852_0    conda-forge
fribidi                   1.0.10               h36c2ea0_0    conda-forge
gcc_impl_linux-64         10.3.0              hf2f2afa_14    conda-forge
gcc_linux-64              10.3.0               hc39de41_9    conda-forge
geos                      3.10.2               h9c3ff4c_0    conda-forge
geotiff                   1.7.0                h6593c0a_6    conda-forge
gettext                   0.19.8.1          h73d1719_1008    conda-forge
gfortran_impl_linux-64    10.3.0              h73f4979_14    conda-forge
gfortran_linux-64         10.3.0               hb09a455_9    conda-forge
giflib                    5.2.1                h36c2ea0_2    conda-forge
gmp                       6.2.1                h58526e2_0    conda-forge
graphite2                 1.3.13            h58526e2_1001    conda-forge
gsl                       2.7                  he838d99_0    conda-forge
gxx_impl_linux-64         10.3.0              hf2f2afa_14    conda-forge
gxx_linux-64              10.3.0               h2593f52_9    conda-forge
harfbuzz                  4.2.0                h40b6f09_0    conda-forge
hdf4                      4.2.15               h10796ff_3    conda-forge
hdf5                      1.12.1          nompi_h4df4325_104    conda-forge
icu                       69.1                 h9c3ff4c_0    conda-forge
jbig                      2.1               h7f98852_2003    conda-forge
jpeg                      9e                   h7f98852_0    conda-forge
json-c                    0.15                 h98cffda_0    conda-forge
kealib                    1.4.14               h87e4c3c_3    conda-forge
kernel-headers_linux-64   2.6.32              he073ed8_15    conda-forge
keyutils                  1.6.1                h166bdaf_0    conda-forge
krb5                      1.19.3               h08a2579_0    conda-forge
lcms2                     2.12                 hddcbb42_0    conda-forge
ld_impl_linux-64          2.36.1               hea4e1c9_2    conda-forge
lerc                      3.0                  h9c3ff4c_0    conda-forge
libblas                   3.9.0           13_linux64_openblas    conda-forge
libcblas                  3.9.0           13_linux64_openblas    conda-forge
libcurl                   7.82.0               h2283fc2_0    conda-forge
libdap4                   3.20.6               hd7c4107_2    conda-forge
libdeflate                1.10                 h7f98852_0    conda-forge
libedit                   3.1.20191231         he28a2e2_2    conda-forge
libev                     4.33                 h516909a_1    conda-forge
libffi                    3.4.2                h7f98852_5    conda-forge
libgcc-devel_linux-64     10.3.0              he6cfe16_14    conda-forge
libgcc-ng                 11.2.0              h1d223b6_14    conda-forge
libgdal                   3.4.1                h4ae554a_2    conda-forge
libgfortran-ng            11.2.0              h69a702a_14    conda-forge
libgfortran5              11.2.0              h5c6108e_14    conda-forge
libglib                   2.70.2               h174f98d_4    conda-forge
libgomp                   11.2.0              h1d223b6_14    conda-forge
libiconv                  1.16                 h516909a_0    conda-forge
libkml                    1.3.0             h238a007_1014    conda-forge
liblapack                 3.9.0           13_linux64_openblas    conda-forge
libnetcdf                 4.8.1           nompi_hb3fd0d9_101    conda-forge
libnghttp2                1.47.0               he49606f_0    conda-forge
libnsl                    2.0.0                h7f98852_0    conda-forge
libopenblas               0.3.18          pthreads_h8fe5266_0    conda-forge
libpng                    1.6.37               h21135ba_2    conda-forge
libpq                     14.2                 h676c864_0    conda-forge
librttopo                 1.1.0                hf69c175_9    conda-forge
libsanitizer              10.3.0              h26c7422_14    conda-forge
libspatialite             5.0.1               h0e567f8_14    conda-forge
libssh2                   1.10.0               ha35d2d1_2    conda-forge
libstdcxx-devel_linux-64  10.3.0              he6cfe16_14    conda-forge
libstdcxx-ng              11.2.0              he4da1e4_14    conda-forge
libtiff                   4.3.0                h542a066_3    conda-forge
libuuid                   2.32.1            h7f98852_1000    conda-forge
libwebp-base              1.2.2                h7f98852_1    conda-forge
libxcb                    1.13              h7f98852_1004    conda-forge
libxml2                   2.9.12               h885dcf4_1    conda-forge
libzip                    1.8.0                h1c5bbd1_1    conda-forge
libzlib                   1.2.11            h166bdaf_1014    conda-forge
lz4-c                     1.9.3                h9c3ff4c_1    conda-forge
make                      4.3                  hd18ef5c_1    conda-forge
ncurses                   6.3                  h9c3ff4c_0    conda-forge
nspr                      4.32                 h9c3ff4c_1    conda-forge
nss                       3.77                 h2350873_0    conda-forge
openjpeg                  2.4.0                hb52868f_1    conda-forge
openssl                   3.0.2                h166bdaf_1    conda-forge
pandoc                    2.17.1.1             ha770c72_0    conda-forge
pango                     1.50.6               hbd2fdc8_0    conda-forge
pcre                      8.45                 h9c3ff4c_0    conda-forge
pcre2                     10.37                h032f7d1_0    conda-forge
pixman                    0.40.0               h36c2ea0_0    conda-forge
poppler                   21.11.0              ha39eefc_0    conda-forge
poppler-data              0.4.11               hd8ed1ab_0    conda-forge
postgresql                14.2                 hce44dc1_0    conda-forge
proj                      8.2.1                h277dcde_0    conda-forge
pthread-stubs             0.4               h36c2ea0_1001    conda-forge
r-abind                   1.4_5           r41hc72bb7e_1003    conda-forge
r-askpass                 1.1               r41hcfec24a_2    conda-forge
r-assertthat              0.2.1             r41hc72bb7e_2    conda-forge
r-backports               1.4.1             r41hcfec24a_0    conda-forge
r-base                    4.1.3                hd930d0e_0    conda-forge
r-base64enc               0.1_3           r41hcfec24a_1004    conda-forge
r-bh                      1.78.0_0          r41hc72bb7e_0    conda-forge
r-biocmanager             1.30.16           r41hc72bb7e_0    conda-forge
r-bit                     4.0.4             r41hcfec24a_0    conda-forge
r-bit64                   4.0.5             r41hcfec24a_0    conda-forge
r-bitops                  1.0_7             r41hcfec24a_0    conda-forge
r-boot                    1.3_28            r41hc72bb7e_0    conda-forge
r-brio                    1.1.3             r41hcfec24a_0    conda-forge
r-bslib                   0.3.1             r41hc72bb7e_0    conda-forge
r-cachem                  1.0.6             r41hcfec24a_0    conda-forge
r-callr                   3.7.0             r41hc72bb7e_0    conda-forge
r-catools                 1.18.2            r41h03ef668_0    conda-forge
r-circlize                0.4.14            r41hc72bb7e_0    conda-forge
r-class                   7.3_20            r41hcfec24a_0    conda-forge
r-classint                0.4_3             r41h859d828_2    conda-forge
r-cli                     3.2.0             r41h03ef668_0    conda-forge
r-clue                    0.3_60            r41hcfec24a_0    conda-forge
r-cluster                 2.1.3             r41h8da6f51_0    conda-forge
r-codetools               0.2_18            r41hc72bb7e_0    conda-forge
r-colorspace              2.0_3             r41h06615bd_0    conda-forge
r-combinat                0.0_8           r41hc72bb7e_1003    conda-forge
r-commonmark              1.8.0             r41h06615bd_0    conda-forge
r-cowplot                 1.1.1             r41hc72bb7e_0    conda-forge
r-cpp11                   0.4.2             r41hc72bb7e_0    conda-forge
r-crayon                  1.5.1             r41hc72bb7e_0    conda-forge
r-crosstalk               1.2.0             r41hc72bb7e_0    conda-forge
r-curl                    4.3.2             r41hcfec24a_0    conda-forge
r-data.table              1.14.2            r41hcfec24a_0    conda-forge
r-dbi                     1.1.2             r41hc72bb7e_0    conda-forge
r-ddrtree                 0.1.5           r41h9c3ff4c_1004    conda-forge
r-deldir                  1.0_6             r41h859d828_0    conda-forge
r-densityclust            0.3             r41h03ef668_1006    conda-forge
r-desc                    1.4.1             r41hc72bb7e_0    conda-forge
r-diffobj                 0.3.5             r41hcfec24a_0    conda-forge
r-digest                  0.6.29            r41h03ef668_0    conda-forge
r-docopt                  0.7.1             r41hc72bb7e_1    conda-forge
r-doparallel              1.0.17            r41hc72bb7e_0    conda-forge
r-dplyr                   1.0.8             r41h7525677_0    conda-forge
r-dqrng                   0.3.0             r41h03ef668_0    conda-forge
r-e1071                   1.7_9             r41h03ef668_0    conda-forge
r-ellipsis                0.3.2             r41hcfec24a_0    conda-forge
r-evaluate                0.15              r41hc72bb7e_0    conda-forge
r-fansi                   1.0.3             r41h06615bd_0    conda-forge
r-farver                  2.1.0             r41h03ef668_0    conda-forge
r-fastica                 1.2_3             r41he454529_0    conda-forge
r-fastmap                 1.1.0             r41h03ef668_0    conda-forge
r-fitdistrplus            1.1_8             r41hc72bb7e_0    conda-forge
r-fnn                     1.1.3             r41h03ef668_2    conda-forge
r-fontawesome             0.2.2             r41hc72bb7e_0    conda-forge
r-foreach                 1.5.2             r41hc72bb7e_0    conda-forge
r-formatr                 1.12              r41hc72bb7e_0    conda-forge
r-fs                      1.5.2             r41h03ef668_0    conda-forge
r-furrr                   0.2.3             r41hc72bb7e_0    conda-forge
r-futile.logger           1.4.3           r41hc72bb7e_1003    conda-forge
r-futile.options          1.0.1           r41hc72bb7e_1002    conda-forge
r-future                  1.24.0            r41hc72bb7e_0    conda-forge
r-future.apply            1.8.1             r41hc72bb7e_0    conda-forge
r-generics                0.1.2             r41hc72bb7e_0    conda-forge
r-getoptlong              1.0.5             r41hc72bb7e_0    conda-forge
r-ggplot2                 3.3.5             r41hc72bb7e_0    conda-forge
r-ggrepel                 0.9.1             r41h03ef668_0    conda-forge
r-ggridges                0.5.3             r41hc72bb7e_0    conda-forge
r-globaloptions           0.1.2             r41ha770c72_0    conda-forge
r-globals                 0.14.0            r41hc72bb7e_0    conda-forge
r-glue                    1.6.2             r41h06615bd_0    conda-forge
r-goftest                 1.2_3             r41hcfec24a_0    conda-forge
r-gplots                  3.1.1             r41hc72bb7e_0    conda-forge
r-gridextra               2.3             r41hc72bb7e_1003    conda-forge
r-grr                     0.9.5           r41h03ef668_1004    conda-forge
r-gtable                  0.3.0             r41hc72bb7e_3    conda-forge
r-gtools                  3.9.2             r41hcfec24a_0    conda-forge
r-hdf5r                   1.3.5             r41h3416e65_0    conda-forge
r-here                    1.0.1             r41hc72bb7e_0    conda-forge
r-hexbin                  1.28.2            r41h8da6f51_0    conda-forge
r-highr                   0.9               r41hc72bb7e_0    conda-forge
r-htmltools               0.5.2             r41h03ef668_0    conda-forge
r-htmlwidgets             1.5.4             r41hc72bb7e_0    conda-forge
r-httpuv                  1.6.5             r41h03ef668_0    conda-forge
r-httr                    1.4.2             r41hc72bb7e_0    conda-forge
r-hunspell                3.0.1             r41h03ef668_0    conda-forge
r-ica                     1.0_2             r41hc72bb7e_2    conda-forge
r-igraph                  1.3.0             r41hf10d5bd_0    conda-forge
r-irlba                   2.3.5             r41he454529_0    conda-forge
r-isoband                 0.2.5             r41h03ef668_0    conda-forge
r-iterators               1.0.14            r41hc72bb7e_0    conda-forge
r-jquerylib               0.1.4             r41hc72bb7e_0    conda-forge
r-jsonlite                1.8.0             r41h06615bd_0    conda-forge
r-kernsmooth              2.23_20           r41h742201e_0    conda-forge
r-knitr                   1.38              r41hc72bb7e_0    conda-forge
r-labeling                0.4.2             r41hc72bb7e_1    conda-forge
r-lambda.r                1.2.4             r41hc72bb7e_1    conda-forge
r-later                   1.2.0             r41h03ef668_0    conda-forge
r-lattice                 0.20_45           r41hcfec24a_0    conda-forge
r-lazyeval                0.2.2             r41hcfec24a_2    conda-forge
r-leiden                  0.3.9             r41hc72bb7e_0    conda-forge
r-leidenbase              0.1.3             r41h1aed7a7_2    bioconda
r-lifecycle               1.0.1             r41hc72bb7e_0    conda-forge
r-listenv                 0.8.0             r41hc72bb7e_1    conda-forge
r-lmtest                  0.9_40            r41h8da6f51_0    conda-forge
r-lobstr                  1.1.1             r41h03ef668_1    conda-forge
r-lsei                    1.3_0             r41h92ddd45_1    conda-forge
r-magrittr                2.0.3             r41h06615bd_0    conda-forge
r-mass                    7.3_56            r41h06615bd_0    conda-forge
r-matrix                  1.4_1             r41h0154571_0    conda-forge
r-matrix.utils            0.9.8             r41hc72bb7e_1    conda-forge
r-matrixstats             0.61.0            r41hcfec24a_0    conda-forge
r-mgcv                    1.8_40            r41h0154571_0    conda-forge
r-mime                    0.12              r41hcfec24a_0    conda-forge
r-miniui                  0.1.1.1         r41hc72bb7e_1002    conda-forge
r-modeldata               0.1.1             r41hc72bb7e_0    conda-forge
r-munsell                 0.5.0           r41hc72bb7e_1004    conda-forge
r-nlme                    3.1_157           r41h8da6f51_0    conda-forge
r-npsurv                  0.5_0             r41hc72bb7e_0    conda-forge
r-openssl                 2.0.0             r41h1f3e0c5_0    conda-forge
r-parallelly              1.30.0            r41hc72bb7e_0    conda-forge
r-patchwork               1.1.1             r41hc72bb7e_0    conda-forge
r-pbapply                 1.5_0             r41hc72bb7e_0    conda-forge
r-pbmcapply               1.5.0             r41hcfec24a_1    conda-forge
r-pheatmap                1.0.12            r41hc72bb7e_2    conda-forge
r-pillar                  1.7.0             r41hc72bb7e_0    conda-forge
r-pkgconfig               2.0.3             r41hc72bb7e_1    conda-forge
r-pkgload                 1.2.4             r41h03ef668_0    conda-forge
r-plotly                  4.10.0            r41hc72bb7e_0    conda-forge
r-plyr                    1.8.7             r41h7525677_0    conda-forge
r-png                     0.1_7           r41hcfec24a_1004    conda-forge
r-polyclip                1.10_0            r41h03ef668_2    conda-forge
r-praise                  1.0.0           r41hc72bb7e_1005    conda-forge
r-processx                3.5.3             r41h06615bd_0    conda-forge
r-promises                1.2.0.1           r41h03ef668_0    conda-forge
r-proxy                   0.4_26            r41hcfec24a_0    conda-forge
r-pryr                    0.1.5             r41h03ef668_0    conda-forge
r-ps                      1.6.0             r41hcfec24a_0    conda-forge
r-pscl                    1.5.5             r41h7f98852_1    conda-forge
r-purrr                   0.3.4             r41hcfec24a_1    conda-forge
r-qlcmatrix               0.9.7           r41hc72bb7e_1002    conda-forge
r-r.methodss3             1.8.2             r41hc72bb7e_1    conda-forge
r-r.oo                    1.25.0            r41hc72bb7e_1    conda-forge
r-r.utils                 2.12.0            r41hc72bb7e_1    conda-forge
r-r6                      2.5.1             r41hc72bb7e_0    conda-forge
r-rann                    2.6.1             r41h03ef668_2    conda-forge
r-rappdirs                0.3.3             r41hcfec24a_0    conda-forge
r-raster                  3.5_15            r41h7525677_1    conda-forge
r-rcolorbrewer            1.1_3             r41h785f33e_0    conda-forge
r-rcpp                    1.0.8.3           r41h7525677_0    conda-forge
r-rcppannoy               0.0.19            r41h03ef668_0    conda-forge
r-rcpparmadillo           0.11.0.0.0        r41h43535f1_0    conda-forge
r-rcppeigen               0.3.3.9.1         r41h306847c_0    conda-forge
r-rcpphnsw                0.3.0             r41h03ef668_0    conda-forge
r-rcppparallel            5.1.5             r41h03ef668_0    conda-forge
r-rcppprogress            0.4.2             r41hc72bb7e_1    conda-forge
r-rcpptoml                0.1.7             r41h03ef668_1    conda-forge
r-rcurl                   1.98_1.6          r41hcfec24a_0    conda-forge
r-rematch2                2.1.2             r41hc72bb7e_1    conda-forge
r-reshape2                1.4.4             r41h03ef668_1    conda-forge
r-reticulate              1.24              r41h03ef668_0    conda-forge
r-rhpcblasctl             0.21_247.1        r41hcfec24a_0    conda-forge
r-rjson                   0.2.21            r41h03ef668_0    conda-forge
r-rlang                   1.0.2             r41h7525677_0    conda-forge
r-rmarkdown               2.13              r41hc72bb7e_1    conda-forge
r-rocr                    1.0_11            r41hc72bb7e_1    conda-forge
r-rpart                   4.1.16            r41hcfec24a_0    conda-forge
r-rprojroot               2.0.3             r41hc72bb7e_0    conda-forge
r-rsample                 0.1.1             r41hc72bb7e_0    conda-forge
r-rspectra                0.16_0            r41h306847c_4    conda-forge
r-rstudioapi              0.13              r41hc72bb7e_0    conda-forge
r-rsvd                    1.0.5             r41hc72bb7e_0    conda-forge
r-rtsne                   0.15              r41h6dc32e9_3    conda-forge
r-runit                   0.4.32          r41hc72bb7e_1002    conda-forge
r-s2                      1.0.7             r41h09fa465_1    conda-forge
r-sass                    0.4.1             r41h7525677_0    conda-forge
r-scales                  1.1.1             r41hc72bb7e_0    conda-forge
r-scattermore             0.8               r41hcfec24a_0    conda-forge
r-sctransform             0.3.3             r41hf361f1a_1    conda-forge
r-seurat                  4.1.0             r41h03ef668_0    conda-forge
r-seuratobject            4.0.4             r41h03ef668_0    conda-forge
r-sf                      1.0_6             r41h7edceaa_0    conda-forge
r-shape                   1.4.6             r41ha770c72_0    conda-forge
r-shiny                   1.7.1             r41h785f33e_0    conda-forge
r-sitmo                   2.0.2             r41h03ef668_0    conda-forge
r-slam                    0.1_50            r41hb699f27_1    conda-forge
r-slider                  0.2.2             r41hcfec24a_0    conda-forge
r-snow                    0.4_4             r41hc72bb7e_0    conda-forge
r-sourcetools             0.1.7           r41h03ef668_1002    conda-forge
r-sp                      1.4_6             r41hcfec24a_0    conda-forge
r-sparsesvd               0.2               r41hcfec24a_1    conda-forge
r-spatstat                2.3_4             r41hc72bb7e_0    conda-forge
r-spatstat.core           2.4_2             r41h7525677_0    conda-forge
r-spatstat.data           2.1_4             r41hc72bb7e_0    conda-forge
r-spatstat.geom           2.4_0             r41h06615bd_0    conda-forge
r-spatstat.linnet         2.3_2             r41hcfec24a_0    conda-forge
r-spatstat.random         2.2_0             r41h7525677_0    conda-forge
r-spatstat.sparse         2.1_0             r41hcfec24a_0    conda-forge
r-spatstat.utils          2.3_0             r41h7f98852_0    conda-forge
r-spdata                  2.0.1             r41hc72bb7e_0    conda-forge
r-spdep                   1.2_3             r41h06615bd_0    conda-forge
r-speedglm                0.3_4             r41hc72bb7e_0    conda-forge
r-spelling                2.2               r41hc72bb7e_0    conda-forge
r-stringi                 1.7.6             r41h337692f_1    conda-forge
r-stringr                 1.4.0             r41hc72bb7e_2    conda-forge
r-survival                3.3_1             r41h06615bd_0    conda-forge
r-sys                     3.4               r41hcfec24a_0    conda-forge
r-tensor                  1.5             r41hc72bb7e_1003    conda-forge
r-terra                   1.5_21            r41hd904c4b_0    conda-forge
r-testthat                3.1.3             r41h7525677_0    conda-forge
r-tibble                  3.1.6             r41hcfec24a_0    conda-forge
r-tidyr                   1.2.0             r41h03ef668_0    conda-forge
r-tidyselect              1.1.2             r41hc72bb7e_0    conda-forge
r-tinytex                 0.38              r41hc72bb7e_0    conda-forge
r-units                   0.8_0             r41h03ef668_0    conda-forge
r-utf8                    1.2.2             r41hcfec24a_0    conda-forge
r-uwot                    0.1.11            r41h03ef668_0    conda-forge
r-vctrs                   0.4.0             r41h7525677_0    conda-forge
r-velocyto.r              0.6               r41h46c59ee_4    bioconda
r-vgam                    1.1_6             r41h859d828_0    conda-forge
r-viridis                 0.6.2             r41hc72bb7e_0    conda-forge
r-viridislite             0.4.0             r41hc72bb7e_0    conda-forge
r-waldo                   0.4.0             r41hc72bb7e_0    conda-forge
r-warp                    0.2.0             r41hcfec24a_1    conda-forge
r-withr                   2.5.0             r41hc72bb7e_0    conda-forge
r-wk                      0.6.0             r41h03ef668_0    conda-forge
r-xfun                    0.30              r41h7525677_0    conda-forge
r-xml                     3.99_0.9          r41h06615bd_0    conda-forge
r-xml2                    1.3.3             r41h03ef668_0    conda-forge
r-xtable                  1.8_4             r41hc72bb7e_3    conda-forge
r-yaml                    2.3.5             r41h06615bd_0    conda-forge
r-zoo                     1.8_9             r41h06615bd_1    conda-forge
readline                  8.1                  h46c0cb4_0    conda-forge
sed                       4.8                  he412f7d_0    conda-forge
sqlite                    3.37.1               h4ff8645_0    conda-forge
sysroot_linux-64          2.12                he073ed8_15    conda-forge
tiledb                    2.6.4                h3f4058f_0    conda-forge
tk                        8.6.12               h27826a3_0    conda-forge
tktable                   2.10                 hb7b940f_3    conda-forge
tzcode                    2022a                h166bdaf_0    conda-forge
tzdata                    2022a                h191b570_0    conda-forge
udunits2                  2.2.28               hc3e0081_0    conda-forge
xerces-c                  3.2.3                h8ce2273_4    conda-forge
xorg-kbproto              1.0.7             h7f98852_1002    conda-forge
xorg-libice               1.0.10               h7f98852_0    conda-forge
xorg-libsm                1.2.3             hd9c2040_1000    conda-forge
xorg-libx11               1.7.2                h7f98852_0    conda-forge
xorg-libxau               1.0.9                h7f98852_0    conda-forge
xorg-libxdmcp             1.1.3                h7f98852_0    conda-forge
xorg-libxext              1.3.4                h7f98852_1    conda-forge
xorg-libxrender           0.9.10            h7f98852_1003    conda-forge
xorg-libxt                1.2.1                h7f98852_2    conda-forge
xorg-renderproto          0.11.1            h7f98852_1002    conda-forge
xorg-xextproto            7.3.0             h7f98852_1002    conda-forge
xorg-xproto               7.0.31            h7f98852_1007    conda-forge
xz                        5.2.5                h516909a_1    conda-forge
zlib                      1.2.11            h166bdaf_1014    conda-forge
zstd                      1.5.2                ha95c52a_0    conda-forge




#### sessionInfo()
R version 4.1.3 (2022-03-10)

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] slingshot_2.2.1             TrajectoryUtils_1.2.0      
 [3] SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
 [5] Biobase_2.54.0              GenomicRanges_1.46.1       
 [7] GenomeInfoDb_1.30.1         IRanges_2.28.0             
 [9] S4Vectors_0.32.4            BiocGenerics_0.40.0        
[11] MatrixGenerics_1.6.0        matrixStats_0.62.0         
[13] princurve_2.1.6             readxl_1.4.1               
[15] patchwork_1.1.2             gridExtra_2.3              
[17] sva_3.42.0                  BiocParallel_1.28.3        
[19] genefilter_1.76.0           mgcv_1.8-40                
[21] nlme_3.1-157                cowplot_1.1.1              
[23] ComplexHeatmap_2.10.0       circlize_0.4.15            
[25] dplyr_1.0.9                 scales_1.2.1               
[27] RColorBrewer_1.1-3          viridis_0.6.2              
[29] viridisLite_0.4.1           ggplot2_3.3.6              
[31] sp_1.5-0                    SeuratObject_4.1.0         
[33] Seurat_4.1.1               

loaded via a namespace (and not attached):
  [1] plyr_1.8.7             igraph_1.3.1           lazyeval_0.2.2        
  [4] splines_4.1.3          listenv_0.8.0          scattermore_0.8       
  [7] digest_0.6.30          foreach_1.5.2          htmltools_0.5.3       
 [10] fansi_1.0.3            magrittr_2.0.3         memoise_2.0.1         
 [13] tensor_1.5             cluster_2.1.3          doParallel_1.0.17     
 [16] ROCR_1.0-11            limma_3.50.3           Biostrings_2.62.0     
 [19] globals_0.16.1         annotate_1.72.0        spatstat.sparse_2.1-1 
 [22] colorspace_2.0-3       blob_1.2.3             ggrepel_0.9.1         
 [25] RCurl_1.98-1.7         crayon_1.5.2           jsonlite_1.8.3        
 [28] progressr_0.11.0       spatstat.data_2.2-0    survival_3.3-1        
 [31] zoo_1.8-10             iterators_1.0.14       glue_1.6.2            
 [34] polyclip_1.10-4        gtable_0.3.1           zlibbioc_1.40.0       
 [37] XVector_0.34.0         leiden_0.4.3           DelayedArray_0.20.0   
 [40] GetoptLong_1.0.5       future.apply_1.9.1     shape_1.4.6           
 [43] abind_1.4-5            edgeR_3.36.0           DBI_1.1.2             
 [46] spatstat.random_2.2-0  miniUI_0.1.1.1         Rcpp_1.0.9            
 [49] xtable_1.8-4           clue_0.3-61            reticulate_1.25       
 [52] spatstat.core_2.4-4    bit_4.0.5              htmlwidgets_1.5.4     
 [55] httr_1.4.4             ellipsis_0.3.2         ica_1.0-3             
 [58] pkgconfig_2.0.3        XML_3.99-0.10          uwot_0.1.14           
 [61] deldir_1.0-6           locfit_1.5-9.5         utf8_1.2.2            
 [64] tidyselect_1.2.0       rlang_1.0.6            reshape2_1.4.4        
 [67] later_1.3.0            AnnotationDbi_1.56.2   cellranger_1.1.0      
 [70] munsell_0.5.0          tools_4.1.3            cachem_1.0.6          
 [73] cli_3.3.0              generics_0.1.3         RSQLite_2.2.14        
 [76] ggridges_0.5.4         stringr_1.4.1          fastmap_1.1.1         
 [79] goftest_1.2-3          bit64_4.0.5            fitdistrplus_1.1-8    
 [82] purrr_0.3.5            RANN_2.6.1             KEGGREST_1.34.0       
 [85] pbapply_1.5-0          future_1.28.0          mime_0.12             
 [88] compiler_4.1.3         plotly_4.10.0          png_0.1-7             
 [91] spatstat.utils_2.3-1   tibble_3.1.7           stringi_1.7.8         
 [94] rgeos_0.5-9            lattice_0.20-45        Matrix_1.5-1          
 [97] vctrs_0.4.1            pillar_1.8.1           lifecycle_1.0.3       
[100] spatstat.geom_2.4-0    lmtest_0.9-40          GlobalOptions_0.1.2   
[103] RcppAnnoy_0.0.19       bitops_1.0-7           data.table_1.14.4     
[106] irlba_2.3.5.1          httpuv_1.6.6           R6_2.5.1              
[109] promises_1.2.0.1       KernSmooth_2.23-20     parallelly_1.32.1     
[112] codetools_0.2-18       MASS_7.3-57            rjson_0.2.21          
[115] withr_2.5.0            sctransform_0.3.5      GenomeInfoDbData_1.2.7
[118] parallel_4.1.3         rpart_4.1.16           tidyr_1.2.1           
[121] Rtsne_0.16             shiny_1.7.2

