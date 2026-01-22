# U-PHAIN-TF
A repository and [webpage](https://sedemto.github.io/phain-tf/) accompanying the TASLP research paper *Audio Inpainting in Time-Frequency Domain
with Phase-Aware Prior*.

> The so-called audio inpainting problem in the time domain refers to estimating missing segments of samples within a signal. Over the years, several methods have been developed for such type of audio inpainting. In contrast to this case, a time-frequency variant of inpainting appeared in the literature, where the challenge is to reconstruct missing spectrogram columns with reliable information. We propose a method to address this time-frequency audio inpainting problem. Our approach is based on the recently introduced phase-aware signal prior that exploits an estimate of the instantaneous frequency. An optimization problem is formulated and solved using the generalized Chambolle–Pock algorithm. The proposed method is evaluated both objectively and subjectively against other time-frequency inpainting methods, specifically a deep-prior neural network and the autoregression-based approach known as Janssen-TF. Our proposed approach surpassed these methods in the objective evaluation as well as in the conducted listening test. Moreover, this outcome is achieved with a substantially reduced computational requirement compared to alternative methods.

## Contents of the repository

The paper compares the recent methods DPAI and JanssenTF with a newly proposed method U-PHAIN-TF.
- DPAI and JanssenTF codes are not a part of this repository but are available at [DPAI](https://github.com/fmiotello/dpai) and [JanssenTF](https://github.com/rajmic/spectrogram-inpainting).
- Matlab codes of our method are available in the `U-PHAIN-TF` folder.
- For reproducibility reasons, the codes are set to read the input (uncorrupted) audio files from the `dataset` folder.
 - `DPAI_originals` is the default dataset, it contains audio files used in DPAI paper.
 - `IRMAS_five_seconds` contains a subset of the IRMAS dataset created using the information available at [JanssenTF](https://github.com/rajmic/spectrogram-inpainting). The original IRMAS dataset can be downloaded [here](https://www.upf.edu/web/mtg/irmas).
- The spectrogram masks used in our experiments are read from the `spectrogram_masks` folder.

Note that to exactly reproduce the SNR and ODG results, first, .wav signals need to be created from the reconstructed signals.

## Dependencies
The Matlab codes for U-PHAIN-TF use the [LTFAT](https://ltfat.org/) and the Signal Processing Toolbox. We used Matlab R2025a in our experiments.

To compute the perceptually-motivated evaluation, we have used the [PEMO-Q package](https://uol.de/en/mediphysics/downloads/pemo-q) (version 1.4.1).
