# A novel aligned Kravchuk transform based on Spin Spherical Harmonics Transform


This project contains the Python code associated to the paper

> **Pascal, B.**, Flamant, J. & Bardenet, R. (2025). ``*Transformée de Kravchuk réalignée, implémentation efficace et illustrations sur quelques signaux élémentaires*". Submitted.  
>  [hal-]()

## Project description

Recently, a novel discrete generalized time-frequency transform has been introduced, namely the *Kravchuk transform*, which represents a discrete complex signal by a function on the unit sphere.

Three demonstration notebooks are provided:

- [`kravchuk-elementary-signals`](notebooks/kravchuk-elementary-signals.ipynb)

> *Illustrates the novel aligned Kravchuk transform on several elementary signals, namely a Dirac delta, a pure sine, a chirp and a coherent state, reproducing `Figure 4` of the paper.*
> 
- [`kravchuk-bat-signal`](notebooks/kravchuk-bat-signal.ipynb)

> *Loads a real recording of bat call used for echolocation and compare its Kravchuk and Fourier spectrograms reproducing `Figure 5` of the paper.*
>
> - [`kravchuk-aligned-signal`](notebooks/kravchuk-aligned-signal.ipynb)

> *Generate an elementary signal and compare its original and novel aligned Kravchuk spectrograms reproducing `Figure 1` of the paper. Also enable to rotate the signals.*

## Dependencies

The following Python libraries are necessary:
- `matplotlib`
- `numpy`
- `scipy`
- `pyssht`

## References

[1] **Pascal, B.**, & Bardenet, R. (2022). A covariant, discrete time-frequency representation tailored for zero-based signal detection. *IEEE Transactions on Signal Processing*, 70, 2950–2961. [hal-03553433](https://hal.archives-ouvertes.fr/hal-03553433/document)

[2] **Pascal, B.**, & Bardenet, R. (2022, September). Une famille de représentations covariantes de signaux discrets et son application à la détection de signaux à partir de leurs zéros. *GRETSI’22 XXVIIIème Colloque Francophone De Traitement Du Signal Et Des Images*. [hal-03614725](https://hal.archives-ouvertes.fr/hal-03614725/document)