# Figure Reproduction

Note: Code for Figures 1, 2, 4, and Supplementary Figure 7 can be found in the corresponding [MATLAB package](https://github.com/neurodata/lol-paper#reproduction-instructions).

### Figure 3

Figure 3 can be reproduced (from scratch) using the [simulations driver](https://github.com/neurodata/lol/blob/master/docs/lol-paper/sim-figure/sims_driver.R), with dependencies suitably installed, from terminal:

```
Rscript sims_driver.R
```

This script, in full, took us about 500 core-hours.

this will produce a file in the directory, `data/sims/lol_sims_lda.rds`, or you can use our version of the outputs, found [here](https://github.com/neurodata/lol/blob/master/docs/lol-paper/data/sims/lol_sims_lda.rds). This can then be used in the Rmarkdown [Figure 3](https://github.com/neurodata/lol/blob/master/docs/lol-paper/sim-figure/sims_driver.R) to produce Figure 3.

### Figure 5

To reproduce Figure 5, the simplest approach, in our opinion, is to leverage the [FlashLol docker container](https://hub.docker.com/r/neurodata/flashlol). This docker container contains the necessary dependencies to use both [FlashX](http://flashx.io/), a package for rapid numerical computation (optionally) using semi-external memory, and [the LOL FlashX implementation](https://github.com/neurodata/lol/blob/master/docs/lol-paper/flashlol-figure/flashLol.R), which also implements PCA, CCA, RRLDA, and Random Projections.

To build the docker container, navigate to the appropriate directory, and build like normal:

```
cd flashR/
docker build neurodata/flashlol:0.0.3 .
```

Optionally, pull the existing docker container:

```
docker pull neurodata/flashlol:0.0.3
```

Next, make a new directory somewhere on your machine; let's assume it's called `<data>`. Create a sub-directory, `<data>/dwi`.

one can navigate to the [neurodata.io/mri](https://neurodata.io/mri/) cloud, and download the "Diffusion MRI" >> "Aligned Images" for all of the individuals with Diffusion MRI connectomes. Put in a directory, called `<data>/<Dataset>/dwi`, for each of the datasets employed.

Next, download the coresponding phenotypic `.csv` file for each dataset, also from [neurodata.io/mri](https://neurodata.io/mri/). Place this at `<data>/<Dataset>/<Dataset>.csv`.

Finally, for each dataset, one can run the file using the docker container:

```
docker run -ti --entrypoint /bin/bash -v <data>/<Dataset>:/brains -v <path>/<to>/<lol>/<repository>/:/lol neurodata/flashlol:0.0.3
cd /lol/docs/lol-paper/flashlol-figure/
Rscript flashlol_corr_driver.R
```
which will create a file at `<data>/<Dataset>/Dataset-<Dataset>_flashlol.rds`. A procedure to set this up on Amazon AMIs can be found [here](https://github.com/neurodata/lol/blob/master/docs/lol-paper/scratch/flashlol.md). 

