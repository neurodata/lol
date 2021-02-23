# Figure Reproduction

Note: Code for Figures 1, 2, 4, and Supplementary Figure 7 can be found in the corresponding [MATLAB package](https://github.com/neurodata/lol-paper#reproduction-instructions).

### Figure 3

Figure 3 can be reproduced (from scratch) using the [simulations driver](https://github.com/neurodata/lol/blob/master/docs/lol-paper/sim-figure/sims_driver.R), with dependencies suitably installed, from terminal:

```
Rscript sims_driver.R
```

This script, in full, took us about 500 core-hours.

this will produce a file in the directory, `data/sims/lol_sims_lda.rds`, or you can use our version of the outputs, found [here](https://github.com/neurodata/lol/blob/master/docs/lol-paper/data/sims/lol_sims_lda.rds). This can then be used in the Rmarkdown [Figure 3](https://github.com/neurodata/lol/blob/master/docs/lol-paper/sim-figure/sims_driver.R) to produce Figure 3.
