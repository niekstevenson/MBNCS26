# BayesFlow / uv Setup

This folder contains a tiny `pyproject.toml` for creating a Python environment with:

- [BayesFlow](https://github.com/bayesflow-org/bayesflow)
- [JAX](https://github.com/jax-ml/jax)
- [HSSM](https://github.com/lnccbrown/HSSM) and [PyMC](https://www.pymc.io/)
- `ipykernel`, so this environment can be used as a Jupyter notebook kernel in VS Code


The dependency manager used here is [uv](https://docs.astral.sh/uv/). Think of uv as a fast tool that creates a clean Python environment, installs packages into it, and remembers the exact versions it chose.

## 1. Install uv

Open a terminal.

On Windows PowerShell:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

On macOS or Linux:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

After installing, close and reopen the terminal. This helps your terminal find the new `uv` command.

Check that uv is available:

```bash
uv --version
```

If that prints a version number, you are good.

## 2. Go to this tutorial folder

From the root of this repository, move into this environment tutorial folder:

```bash
cd Day5-SBI
```

You should see a file called `pyproject.toml`. That file is the recipe `uv` reads.

Important: you do not run the TOML file directly. You run `uv` commands in the same folder as the TOML file.

## 3. Create the environment and install the packages

Run:

```bash
uv sync
```

uv will:

- read `pyproject.toml`
- choose a compatible Python version
- create a local `.venv` folder
- install `BayesFlow`, `HSSM`, `PyMC`, `JAX`, and `ipykernel`
- create a `uv.lock` file with the exact resolved package versions

If uv says it cannot find the right Python version, install Python 3.12 through uv:

```bash
uv python install 3.12
uv sync --python 3.12
```

## 4. Run Python inside the uv environment

Use `uv run` whenever you want to run Python with these packages available:

```bash
uv run python
```

That opens a Python prompt. You can then try:

```python
import bayesflow
```

BayesFlow will print an informative message if it finds the DL backend (i.e., JAX).

To leave the Python prompt:

```python
exit()
```

## 5. Use notebooks in VS Code

Install the VS Code Python and Jupyter extensions if you do not already have them.

Open this repository folder in VS Code. The workspace setting in `.vscode/settings.json` points VS Code at:

```text
Day5-SBI/.venv
```

After you run `uv sync`, VS Code should be able to find that environment.

For a notebook:

1. Open or create an `.ipynb` file.
2. Click `Select Kernel` in the top right of the notebook.
3. Choose `Python Environments`.
4. Select the `.venv` under `Day5-SBI`.

## 6. CPU vs GPU

The default `pyproject.toml` is for CPU installs. This is the easiest setup and is the right starting point for most people.

If you want NVIDIA GPU support, check the platform-specific install notes first:

- [JAX installation guide](https://docs.jax.dev/en/latest/installation.html)

GPU installs depend on your operating system, CUDA version, and driver version, so do not guess here.

## 7. Common fixes

If `uv` is not found, close and reopen your terminal.

If package resolution fails, try using Python 3.12:

```bash
uv python install 3.12
uv sync --python 3.12
```

If VS Code does not show the environment, make sure you ran `uv sync` from `Day5-SBI`, then reload the VS Code window and use `Python: Select Interpreter` from the Command Palette.

If a notebook asks to install a kernel, it usually means `ipykernel` is missing from the selected environment. This setup includes it, so run `uv sync` again from `tutorials`.

If you get GPU-related JAX errors, switch back to the CPU setup first and confirm the basic install works.

## 8. What the TOML file means

The important part of `pyproject.toml` is its `dependencies = [...]` list — read it there rather than here, so the two cannot drift apart. It names the packages `uv sync` installs into the local tutorial environment.
