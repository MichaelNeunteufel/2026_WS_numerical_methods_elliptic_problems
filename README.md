# Numerical Methods for Elliptic Problems

Public course material for the lecture *Numerical Methods for Elliptic Problems*, winter term 2026 at Johannes Kepler University Linz (JKU Linz).

## Use the notebooks

- Run in the browser with JupyterLite: [open the course index](https://michaelneunteufel.github.io/2026_WS_numerical_methods_elliptic_problems/lab?path=index.ipynb)
- Run locally: install NGSolve and JupyterLab, then start `python -m jupyter lab` in a clone of this repository.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ngsolve jupyterlab anywidget
python -m jupyter lab
```

On Windows, activate the environment with `.venv\Scripts\Activate.ps1`.

## Repository layout

- `content/index.ipynb`: entry point for JupyterLite
- `content/notebooks/`: lecture notebooks
- `content/exercises/`: exercise sheets and notebooks
- `content/script/`: lecture notes
