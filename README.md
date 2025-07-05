# phain-tf
A repository and [webpage](____________) accompanying the paper *_________PHAIN_TF_______________*.

> ****Abstract****

## Contents of the repository

The paper compares the recent methods abbreviated DPAI and Janssen-TF with a newly proposed method U-PHAIN-TF.
- DPAI codes are not a part of this repository but are available [here](https://github.com/fmiotello/dpai).
- Janssen_TF codes are not a part of this repository but are available [here](https://github.com/rajmic/spectrogram-inpainting).
- Matlab codes of our method are available in the `U-PHAIN-TF` folder.
- For reproducibility reasons, the codes are set to read the input (uncorrupted) audio files from the `dataset` folder.
- The spectrogram masks used in our experiments are read from the `spectrogram_masks` folder.
- ...
 <!--- Regarding the experiment using the IRMAS dataset, the folder `audio-irmas` includes a list of the files used in our experiment and a Matlab script which crops the files to a length of 5 seconds and subsamples them to 16 kHz. The original files can be downloaded [here](https://www.upf.edu/web/mtg/irmas). -->

...
to do
...


## Dependencies
The Matlab codes for U-PHAIN-TF use the [LTFAT](https://ltfat.org/) and the Signal Processing Toolbox. To compute the perceptually-motivated evaluation, we have used the [PEMO-Q package](https://uol.de/en/mediphysics/downloads/pemo-q) (version 1.4.1).
