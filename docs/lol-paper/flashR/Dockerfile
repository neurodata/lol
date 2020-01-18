FROM flashx/flashx:latest

RUN Rscript -e "install.packages(c('oro.nifti', 'ramify', 'profmem', 'tidyverse'), repo='http://cran.r-project.org')"
RUN mkdir /data
COPY MNI152_T1_1mm_brain_mask.nii.gz /data/
RUN Rscript -e "install.packages(c('https://cran.r-project.org/src/contrib/Archive/mvtnorm/mvtnorm_1.0-8.tar.gz'))"
RUN Rscript -e "install.packages(c('robust'), repo='http://cran.r-project.org')"
RUN Rscript -e "require(devtools); install_github('neurodata/lol')"
RUN Rscript -e "install.packages('randomForest', repos='http://cran.r-project.org')"
