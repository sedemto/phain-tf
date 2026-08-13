# U-PHAIN-TF
A repository and [webpage](https://sedemto.github.io/phain-tf/) accompanying the preprint research paper *Audio Inpainting in Time-Frequency Domain
with Phase-Aware Prior*.

> We address the problem of time-frequency audio inpainting, where the goal is to fill missing spectrogram portions with consistent information. Despite recent advances, existing approaches still face limitations in both reconstruction quality and computational efficiency. To bridge this gap, we propose a method that utilizes a phase-aware signal prior which exploits estimates of the instantaneous frequency. An optimization problem is formulated and solved using the generalized Chambolle–Pock algorithm. The proposed method is evaluated against other time-frequency inpainting methods, specifically a deep-prior audio inpainting neural network, the autoregression-based approach known as Janssen-TF, and a sparsity-driven baseline. For short gap durations, the proposed approach achieves superior SNR, while performing comparably to Janssen-TF on larger gaps. In terms of perceptual quality (both objective and subjective), the proposed method consistently outperforms existing methods across all gap lengths. In addition, the reconstructions are obtained with a substantially reduced computational cost compared to alternative methods.

The submitted manuscript is available at [arXiv](https://arxiv.org/abs/2601.18535).
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
