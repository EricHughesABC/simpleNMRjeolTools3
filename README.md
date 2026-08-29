# simpleNMRjeolTools

JEOL/JASON client for [simpleNMR](https://github.com/EricHughesABC/simpleNMRtools) —
reads a `.jjh5` NMR file, lets the user assign spectrum types via a shared
dialog, builds the JSON payload the simpleNMR server expects, submits it,
and opens the result in a PyQt viewer.

## Installing

```bash
pip install -e .
```

This installs a `simplenmr-jeol` command (via `[project.scripts]`) and
pulls in [`simpleNMRbuilder`](https://github.com/EricHughesABC/simpleNMRbuilder)
— the shared JSON-contract/dialog/submission library also used by the
Bruker converter — as a real dependency.

## Configuring JASON's External Tools

Point JASON's **Program** field at the installed console script's fixed
path, e.g.:

```
/path/to/your/conda/env/bin/simplenmr-jeol
```

Leave **Arguments** empty — the script finds its input automatically
(JASON drops `input.jjh5` in the run folder it creates and sets that as
the working directory; `find_input_jjh5()` picks it up from there).

This replaces the previous setup, which pointed JASON's Program field
directly at the raw Python interpreter with a hardcoded path to
`simpleNMRjeolTools_v8.py` as an argument — that path broke every time
the source moved. The console script's path is stable across code
changes; only a `pip install -e .` reinstall (needed if the package's
own dependencies change) would move it.

## Structure

```
src/simplenmrjeol/
    __init__.py             # exports jeolData
    __main__.py             # enables `python -m simplenmrjeol`
    json_converter.py       # jeolData class: reads .jjh5, builds the JSON
                             #   payload (mirrors simpleNMRbrukerTools'
                             #   core/json_converter.py naming/role)
    jason_simpleNMR_cli.py  # the runnable program: commandline()/main(),
                             #   JASON environment quirks, dialog flow,
                             #   server submission, opening the viewer
tests/
    test_builder_integration.py
```

`archive/` holds everything from before this repo was turned into a
proper installable package — see `archive/README.md` for what's there
and why it wasn't just deleted.

## Development

```bash
pip install -e ".[test]"
pytest
```
