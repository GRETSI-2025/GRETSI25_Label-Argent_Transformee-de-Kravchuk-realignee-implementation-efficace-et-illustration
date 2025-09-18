[![](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![](https://img.shields.io/badge/poetry-v2.+-blue.svg)](https://python-poetry.org/)
[![Poetry](https://github.com/bpascal-fr/from-kravchuk-to-ssht/actions/workflows/install-with-poetry.yml/badge.svg?branch=main)](https://github.com/bpascal-fr/from-kravchuk-to-ssht/actions/workflows/install-with-poetry.yml)
[![Requirements.txt](https://github.com/bpascal-fr/from-kravchuk-to-ssht/actions/workflows/install-with-requirements.yml/badge.svg?branch=main)](https://github.com/bpascal-fr/from-kravchuk-to-ssht/actions/workflows/install-with-requirements.yml)

# A novel aligned Kravchuk transform based on Spin Spherical Harmonics Transform

<hr>

**_Dépôt labelisé dans le cadre du [Label Reproductible du GRESTI'25](https://gretsi.fr/colloque2025/recherche-reproductible/)_**

| Label décerné | Auteur | Rapporteur | Éléments reproduits | Liens |
|:-------------:|:------:|:----------:|:-------------------:|:------|
| ![](label_or.png) | Barbara PASCAL<br>[@bpascal-fr](https://github.com/bpascal-fr) | Mathilde DUPOUY<br>[@MathildeDupouy](https://github.com/MathildeDupouy) |  Figures 1, 3, 4 et 5 | 📌&nbsp;[Dépôt&nbsp;original](https://github.com/bpascal-fr/from-kravchuk-to-ssht)<br>⚙️&nbsp;[Issue](https://github.com/GRETSI-2025/Label-Reproductible/issues/24)<br>📝&nbsp;[Rapport](https://github.com/GRETSI-2025/Label-Reproductible/tree/main/rapports/Rapport_issue_24) |

<hr>

This project contains the Python code associated to the paper

> **Pascal, B.**, Flamant, J. & Bardenet, R. (2025). ``*Transformée de Kravchuk réalignée, implémentation efficace et illustration sur signaux élémentaires et réels*". Submitted.  [[pdf]](paper/2025-kravchuk-ssht-hal.pdf)  
>  [hal-05013793](https://hal.science/hal-05013793)

## Project description

Recently, a novel discrete generalized time-frequency transform has been introduced, namely the *Kravchuk transform*, which represents a discrete complex signal by a function on the unit sphere.

Three demonstration notebooks are provided:

- [`kravchuk-elementary-signals`](notebooks/kravchuk-elementary-signals.ipynb)

> *Illustrates the novel aligned Kravchuk transform on several elementary signals, namely a Dirac delta, a pure sine, a chirp and a coherent state, reproducing `Figure 4` of the paper.*
 
- [`kravchuk-bat-signal`](notebooks/kravchuk-bat-signal.ipynb)

> *Loads a real recording of bat call used for echolocation and compare its Kravchuk and Fourier spectrograms reproducing `Figure 5` of the paper.*

- [`kravchuk-aligned-transform`](notebooks/kravchuk-aligned-transform.ipynb)

> *Generate an elementary signal and compare its original and novel aligned Kravchuk spectrograms reproducing `Figure 1` of the paper if the signal is chosen to be the **chirp** with `snr = Inf`.  
> Compare the computational load of the original and the novel aligned Kravchuk transform on one example.  
> Evaluated the precision of the implementation of the aligned Kravchuk transform and of its inverse on one example.  
> Exemplify the rotation of a signal.*

- [`kravchuk-time-precision`](notebooks/kravchuk-time-precision.ipynb)

> *Compare the computational load of the original and the novel aligned Kravchuk transform for different sizes of signal across several realizations to reproduce  `Figure 3` (top).  
> Evaluated the signal-to-error ratio of the implementation of the aligned Kravchuk transform and of its inverse for different sizes of signal across several realizations to reproduce `Figure 3` (bottom).*

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

### Font handling and rendering with LaTeX

The displays of signals and Kravchuk spectrogram leverages the ability of [`matplotlib`](https://matplotlib.org/) to render text using LaTeX. Instructions are provided in the documentation on [matplotlib backends](https://matplotlib.org/stable/install/dependencies.html#optional-dependencies) and [text rendering with LaTeX](https://matplotlib.org/stable/users/explain/text/usetex.html).

You're now ready to compute Kravchuk transforms! 🚀


## References

[1] **Pascal, B.**, & Bardenet, R. (2022). A covariant, discrete time-frequency representation tailored for zero-based signal detection. *IEEE Transactions on Signal Processing*, 70, 2950–2961. [hal-03553433](https://hal.archives-ouvertes.fr/hal-03553433/document)

[2] **Pascal, B.**, & Bardenet, R. (2022, September). Une famille de représentations covariantes de signaux discrets et son application à la détection de signaux à partir de leurs zéros. *GRETSI’22 XXVIIIème Colloque Francophone De Traitement Du Signal Et Des Images*. [hal-03614725](https://hal.archives-ouvertes.fr/hal-03614725/document)

[3] **Pascal, B.**, Flamant J., & Bardenet, R. (2025). Transformée de Kravchuk réalignée, implémentation efficace et illustration sur signaux élémentaires et réels. *Submitted*. [hal-05013793](https://hal.science/hal-05013793)
