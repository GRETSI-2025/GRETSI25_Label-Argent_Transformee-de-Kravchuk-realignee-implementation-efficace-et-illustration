# A novel aligned Kravchuk transform based on Spin Spherical Harmonics Transform


This project contains the Python code associated to the paper

> **Pascal, B.**, Flamant, J. & Bardenet, R. (2025). ``*Transformée de Kravchuk réalignée, implémentation efficace et illustrations sur quelques signaux élémentaires*". Submitted.  
>  [hal-]()

## Project description

Recently, a novel discrete generalized time-frequency transform has been introduced, namely the *Kravchuk transform*, which represents a discrete complex signal by a function on the unit sphere.

Three demonstration notebooks are provided:

- [`demo-kravchuk-bat`](kravchuk-bat-signal.ipynb)

> This toolbox provides a stable implementation of this novel *Kravchuk* transform and the code to reproduce `Figures 1, 2 and 6` of the paper ``*A covariant, discrete time-frequency representation tailored for zero-based signal detection*", comparing the standard and the *Kravchuk* spectrograms of noisy chirps, with a peculiar focus on the zeros.  
> A demonstration is given in the notebook [`kravchuk-spectrogram-and-zeros`](/demos/kravchuk-spectrogram-and-zeros.ipynb).

A novel efficient methodology relying on the functional statistics of the point process formed by the zeros of the Kravchuk spectrogram for detecting the presence of some signal is implemented.

> The detection procedure based on the functional statistics of the zeros of the Kravchuk spectrogram is implemented.
> For sake of comparison, we provide also an implementation of the counterpart strategy relying on the zeros of the Short-Time Fourier transform developed in the paper ``*On the zeros of the spectrogram of white noise*" by Bardenet R., Flamant, J. & Chainais, P. (2021) Applied and Computational Harmonic Analysis.
>
> The interested reader can then reproduce `Figures 7, 8 and 9` of the paper ``*A covariant, discrete time-frequency representation tailored for zero-based signal detection*".  
> A demonstration is given in the notebook [`detection-test-Kravchuk-zeros`](/demos/detection-test-Kravchuk-zeros.ipynb).

## Dependencies

The following Python libraries are necessary:
- `matplotlib`
- `numpy`
- `scipy`
- `pyssht`

## References

[1] Pascal, B., & Bardenet, R. (2022). A covariant, discrete time-frequency representation tailored for zero-based signal detection. *IEEE Transactions on Signal Processing*, 70, 2950–2961.

[2] Pascal, B., & Bardenet, R. (2022, September). Une famille de représentations covariantes de signaux discrets et son application à la détection de signaux à partir de leurs zéros. *GRETSI’22 XXVIIIème Colloque Francophone De Traitement Du Signal Et Des Images*.