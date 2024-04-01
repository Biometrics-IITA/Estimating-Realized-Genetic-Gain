
## Introduction

Estimating genetic gain is crucial in plant breeding as it allows
breeders to quantify the progress in developing new varieties with
desirable traits. By measuring the increase in performance over
generations, breeders can evaluate the effectiveness of their breeding
strategies and adjust breeding strategies and logistics as necessary for
increased efficiency.

## The data

The data [data](./data/cassava_uyt_2016_2021.csv) are from trials
regularly conducted at IITA experiment stations and other locations
within West Africa through collaborative efforts with partners, from
2016 to 2021. This forms an integral component of the IITA cassava
breeding strategy, aimed at channeling elite products into the national
agricultural research system.

## The method

The realized genetic gain [estimation](realized%20genetic%20gain.R) is
done in three steps:

- Step 1: Linear mixed model (LMM) is used for single-trial analysis to
  estimate heritability for each environment. Environments with a
  heritability lower than 0.2 are discarded.

- Step 2: Done using a two-stages approach:

  - LMM is used to estimate BLUEs and corresponding weights for each
    environment.
  - Combined LMM is performed to estimate combined BLUEs and
    corresponding weights.

- Step 3: A weighted-linear regression of the combined BLUEs on the year
  of origin

## Acknowledgement

The IITA West Africa cassava breeding program, in collaboration with
NARS partners, has generated these data across 84 environments for 6
years.
