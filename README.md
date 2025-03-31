[![Poetry](https://github.com/bpascal-fr/from-kravchuk-to-ssht/actions/workflows/install-with-poetry.yml/badge.svg?branch=main)](https://github.com/bpascal-fr/from-kravchuk-to-ssht/actions/workflows/install-with-poetry.yml)
[![Requirements.txt](https://github.com/bpascal-fr/from-kravchuk-to-ssht/actions/workflows/install-with-requirements.yml/badge.svg?branch=main)](https://github.com/bpascal-fr/from-kravchuk-to-ssht/actions/workflows/install-with-requirements.yml)

# A novel aligned Kravchuk transform based on Spin Spherical Harmonics Transform


This project contains the Python code associated to the paper

> **Pascal, B.**, Flamant, J. & Bardenet, R. (2025). ``*Transformée de Kravchuk réalignée, implémentation efficace et illustrations sur quelques signaux élémentaires*". Submitted.  
>  [hal-]()

## Project description

Recently, a novel discrete generalized time-frequency transform has been introduced, namely the *Kravchuk transform*, which represents a discrete complex signal by a function on the unit sphere.

Three demonstration notebooks are provided:

- [`kravchuk-elementary-signals`](notebooks/kravchuk-elementary-signals.ipynb)

> *Illustrates the novel aligned Kravchuk transform on several elementary signals, namely a Dirac delta, a pure sine, a chirp and a coherent state, reproducing `Figure 4` of the paper.*
 
- [`kravchuk-bat-signal`](notebooks/kravchuk-bat-signal.ipynb)

> *Loads a real recording of bat call used for echolocation and compare its Kravchuk and Fourier spectrograms reproducing `Figure 5` of the paper.*

- [`kravchuk-aligned-signal`](notebooks/kravchuk-aligned-transform.ipynb)

> *Generate an elementary signal and compare its original and novel aligned Kravchuk spectrograms reproducing `Figure 1` of the paper if the signal is chosen to be the **chirp** with `snr = Inf`.
> Compare the computational load of the original and the novel aligned Kravchuk transform as reported in `Figure 3` (top).
> Evaluated the precision of the implementation of the aligned Kravchuk transform and of its inverse as reported in `Figure 3` (bottom).
> Exemplify the rotation of a signal.*

## Installation Instructions

### Using Poetry (recommended)

Poetry is a dependency manager for Python that simplifies package management. To install dependencies using Poetry, follow these steps:

1. Ensure [Poetry](https://python-poetry.org/docs/) (v2.0+) is installed. If not, install it using:
   ```sh
   curl -sSL https://install.python-poetry.org | python3 -
   ```
   Or, if you already have Poetry, update it:
   ```sh
   poetry self update
   ```

2. Navigate to your project directory:
   ```sh
   cd /path/to/your/project
   ```

3. Install dependencies from `pyproject.toml`:
   ```sh
   poetry install
   ```
   The installation will automatically create a virtual environment ``.venv/``in the same folder, that you can use for running files in ``notebooks/``.

### Using requirements.txt

If you do want to use Poetry, installation is also possible using the proposed `requirements.txt` file.
Follow these steps:

1. Ensure you have Python and pip installed. You can check by running:
   ```sh
   python3 --version
   pip --version
   ```

2. Navigate to your project directory:
   ```sh
   cd /path/to/your/project
   ```

3. Create and activate a virtual environment (recommended):
   ```sh
   python3 -m venv venv
   source venv/bin/activate  # On Windows use: venv\Scripts\activate
   ```

4. Install dependencies from `requirements.txt`:
   ```sh
   pip install -r requirements.txt
   ```

You're now ready to compute Kravchuk transforms! 🚀


## References

[1] **Pascal, B.**, & Bardenet, R. (2022). A covariant, discrete time-frequency representation tailored for zero-based signal detection. *IEEE Transactions on Signal Processing*, 70, 2950–2961. [hal-03553433](https://hal.archives-ouvertes.fr/hal-03553433/document)

[2] **Pascal, B.**, & Bardenet, R. (2022, September). Une famille de représentations covariantes de signaux discrets et son application à la détection de signaux à partir de leurs zéros. *GRETSI’22 XXVIIIème Colloque Francophone De Traitement Du Signal Et Des Images*. [hal-03614725](https://hal.archives-ouvertes.fr/hal-03614725/document)

[3] **Pascal, B.**, Flamant J., & Bardenet, R. (2025). Transformée de Kravchuk réalignée, implémentation efficace et illustration sur signaux élémentaires et réels. *Submitted*. [hal-](https://hal.archives-ouvertes.fr/)
