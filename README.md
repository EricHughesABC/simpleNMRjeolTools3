# simpleNMRjeolTools

JEOL/JASON client for [simpleNMR](https://github.com/EricHughesABC/simpleNMRtools) —
reads a `.jjh5` NMR file, lets the user assign spectrum types via a shared
dialog, builds the JSON payload the simpleNMR server expects, submits it,
and opens the result in a PyQt viewer.

## Installing

```bash
pip install lxml
pip install -i https://test.pypi.org/simple/ beautifuljason==1.3.0b1
pip install -e .
```

The first two lines are needed while `beautifuljason` 1.3.0 (which adds
`Molecule.set_assignments()`, used to write assignments back to JASON) is
only on TestPyPI. Install `lxml` from normal PyPI first: with
`-i https://test.pypi.org/simple/`, pip looks for *every* missing
dependency on TestPyPI, which doesn't host `lxml`. Once 1.3.0 is on PyPI,
`pip install -e .` alone is enough.

This installs a `simplenmr-jeol` command (via `[project.scripts]`) and
pulls in [`simpleNMRbuilder`](https://github.com/EricHughesABC/simpleNMRbuilder)
— the shared JSON-contract/dialog/submission library also used by the
Bruker converter — as a real dependency.

## Configuring JASON's External Tools

Create the tool with these settings:

- **Output → Result Type:** Document
- **Output → Output Mode:** Modify the input file
- **Output → After the Tool Finishes:** Replace the current document
- **Command → Arguments** (one per line):

### step 1
<img width="521" height="470" alt="image" src="https://github.com/user-attachments/assets/62fd2c7c-e417-42e0-a352-9e45095be616" />

### step 2
<img width="522" height="472" alt="image" src="https://github.com/user-attachments/assets/40a80891-b0cf-40a7-879b-55a30b37192d" />

### step 3
<img width="525" height="467" alt="image" src="https://github.com/user-attachments/assets/b5c87618-9c56-4faa-9452-78077c09a727" />

### step 4
<img width="520" height="469" alt="image" src="https://github.com/user-attachments/assets/7bcbffcf-4080-4009-ba50-353ac047fadf" />

### step 5
<img width="521" height="563" alt="image" src="https://github.com/user-attachments/assets/873ccef5-cb8d-4a12-9fcc-8121e5b1fc24" />

### Notes
  - **in step 1**, set the path to your local installation of python. On Mac Os with Mx chipsets the python needs to be a native installation for the Mx processor.
  - **in step 3**, if you do not want to run the local simpleNMR server then change the address from https://127.0.0.1:5000 to https://simplenmr.pythonanywhere.com/


  ```
  /path/to/your/conda/env/bin/simplenmr-jeol
  <input>
  ```

JASON creates a run folder for each execution, saves the current document
there as `input.jjh5`, and passes that path in place of `<input>`. (If no
path is given, `find_input_jjh5()` looks for `input.jjh5` in the working
directory, which JASON also sets to the run folder.) Add `--backup` as a
third line to keep a timestamped copy of `input.jjh5` before it is
modified; it is off by default because the files can be large.

### What happens in a run

1. The spectrum-assignment dialog runs, the data is submitted to the
   simpleNMR server, and the result opens in the viewer (a separate
   process; this program waits for it to close).
2. Pressing **Export** in the viewer saves `html/input_py.json` in the run
   folder.
3. When the viewer closes, `assignment_writer` checks the export against
   `input.jjh5` (atoms must match atom-for-atom; every ppm value must
   match a peak exactly) and shows a summary. On **Yes**, the assignments
   are written into `input.jjh5` with `Molecule.set_assignments()`.
4. The program exits and JASON reloads `input.jjh5`, replacing the
   current document.

If the viewer is closed without exporting, or the summary is declined,
`input.jjh5` is left unchanged.

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
                             #   server submission, opening the viewer,
                             #   writing assignments back after it closes
    assignment_writer.py    # matches the viewer's export to the .jjh5 and
                             #   writes it with Molecule.set_assignments()
tests/
    test_builder_integration.py
    test_assignment_writer.py   # mock JASON document from data/alpha_ionone_py.json
```

`archive/` holds everything from before this repo was turned into a
proper installable package — see `archive/README.md` for what's there
and why it wasn't just deleted.

## Development

```bash
pip install -e ".[test]"
pytest
```
