FROM flashx/flashx:latest

RUN Rscript -e "install.packages(c('oro.nifti', 'ramify', 'profmem', 'tidyverse'), repo='http://cran.r-project.org')"
RUN mkdir /data
COPY MNI152_T1_1mm_brain_mask.nii.gz /data/
