# phain-tf
A repository and [webpage](____________) accompanying the paper *_________PHAIN_TF_______________*.

> ****Abstract****

## Contents of the repository

The paper compares the recent methods abbreviated DPAI and JanssenTF with a newly proposed method U-PHAIN-TF.
- DPAI and JanssenTF codes are not a part of this repository but are available at [Git &mdash DPAI](https://github.com/fmiotello/dpai) and [Git &mdash JanssenTF](https://github.com/rajmic/spectrogram-inpainting).
- Matlab codes of our method are available in the `U-PHAIN-TF` folder.
- For reproducibility reasons, the codes are set to read the input (uncorrupted) audio files from the `dataset` folder.
 - `DPAI_originals` is the default dataset, it contains audio files used in DPAI paper
 - `IRMAS_originals` contains a subset of the IRMAS dataset create using the information available at [Git &mdash JanssenTF](https://github.com/rajmic/spectrogram-inpainting). The original IRMAS dataset can be downloaded [here](https://www.upf.edu/web/mtg/irmas)
- The spectrogram masks used in our experiments are read from the `spectrogram_masks` folder.
- ...

...
to do
...


## Dependencies
The Matlab codes for U-PHAIN-TF use the [LTFAT](https://ltfat.org/) and the Signal Processing Toolbox. To compute the perceptually-motivated evaluation, we have used the [PEMO-Q package](https://uol.de/en/mediphysics/downloads/pemo-q) (version 1.4.1).
