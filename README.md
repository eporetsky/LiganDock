# Ligandock

Ligandock is a flexible and automated molecular docking toolkit built to streamline ligand and receptor preparation, grid box definition, and batch docking workflows. It interfaces with AutoDock Vina, Meeko, RDKit, OpenBabel, and PyKVFinder for state-of-the-art small molecule docking workflows.

## Key Features
- **Automated ligand/receptor prep**: RDKit, Meeko, and OpenBabel integration.
- **Flexible grid box definition:**
  - Center by marker ligand (PDB/SDF),
  - Center by select residues,
  - Use PyKVFinder-defined pocket.
- **Batch processing**: Dock one or more ligands across single or multiple protein structures.
- **Simple modern command line interface (CLI)** (via [Click](https://click.palletsprojects.com/)), ready for workflows and scripting.
- **Beautiful output:** Docked poses, energy tables, and grid box visualizations.

## Installation

```bash
git clone https://github.com/eporetsky/LiganDock.git
cd LiganDock/
conda env create -f environment.yml
conda activate ligandock
```

## Command Line Usage

### General Example

Dock one or more ligands (ligand names should be exactly 3 letters) against one protein or a directory of proteins:

```bash
ligandock dock \
    --smiles "SMILES_STRING1" --ligand_name "LG1" \
    --smiles "SMILES_STRING2" --ligand_name "LG2" \
    --pdb "path/to/myprotein.pdb" \
    --output_folder "results" \
    --box_size 20 20 20
```

For grid box definition, use **one** of these options:
- `--marker_file MARKER.SDF`
- `--res_box CHAIN:RESID (A213,A245,...) ...`
- `--pykvfinder_box CHAIN+RESID (A213)`

### For folders
If you want to dock against all `.pdb` files in a folder:
```bash
ligandock dock ... --pdb myreceptors/
```

## Example Workflows: Docking GPP into (+)-Limonene Synthase

These workflows demonstrate docking geranyl diphosphate (GPP) into the active site of (+)-Limonene Synthase (PDB: 5UV1). Three different methods are shown for defining the center of the binding box:

**Preparation:**
- The original ligand is removed from the clean structure (`5UV1.clean.pdb`).
- The marker file (`0FV.sdf`) contains the coordinates of the crystallographic ligand as reference for grid placement.

### Example 1: Using a Marker Ligand File
Use a user-provided SDF ligand file (e.g., from crystal structure, DiffDock) to center the binding box:

```bash
ligandock dock \
    --smiles "CC(=CCC/C(=C/COP(=O)(O)OP(=O)(O)O)/C)C" \
    --ligand_name "GPP" \
    --pdb "example/5UV1.clean.pdb" \
    --marker_file "example/0FV.sdf" \
    --output_folder "example_lig" \
    --box_size 20 20 20
```

### Example 2: Using PyKVFinder Pocket Detection
Define a box around a binding pocket adjacent to a user-provided residue (binding pocket identified using PyKVFinder):

```bash
ligandock dock \
    --smiles "CC(=CCC/C(=C/COP(=O)(O)OP(=O)(O)O)/C)C" \
    --ligand_name GPP \
    --pdb example/5UV1.clean.pdb \
    --output_folder example_pocket \
    --pykvfinder_box A343
```

### Example 3: Using Residue Centers
Center the box on any number of user-provided residues (e.g., on either side of a known binding pocket):

```bash
ligandock dock \
    --smiles "CC(=CCC/C(=C/COP(=O)(O)OP(=O)(O)O)/C)C" \
    --ligand_name GPP \
    --pdb example/5UV1.clean.pdb \
    --output_folder example_res \
    --res_box A343 A488
```

**Output:**
- Docked structures in `example/` (poses in `.pdbqt`)
- Grid box and visualization files
- Tabular docking energies

**Visualization Example:**
- The image below (see `example/5UV1.clean.GPP.png`) shows an overlap of:
  - The **original ligand** (from 0FV, orange)
  - The **docked GPP ligand** (cyan)
  - The clean PDB structure (transparent cartoon)

![Docking visualization](example/5UV1.clean.GPP.png)

**Caption:**
> Top docked GPP pose (cyan) as predicted by `ligandock`, overlaying the position of the original crystallographic ligand (2-Fluorogeranyl Diphosphate; orange sticks) in the (+)-Limonene Synthase (PDB 5UV1, cartoon).


## Options Reference
- `--smiles [SMILES]` : Ligand SMILES (repeat for more than one)
- `--ligand_name [NAME]` : Ligand names (repeat, positional match to --smiles)
- `--pdb [FILE_or_DIR]`: Receptor PDB file or folder containing PDBs
- `--output_folder [FOLDER]` : Where output files go
- `--box_size X Y Z` : Grid box size (Angstrom)
- `--marker_file` / `--res_box` / `--pykvfinder_box` : Box centering options
- [`more options --run 'ligandock dock --help'`]

## Output Files
- Docked ligand poses (`.pdbqt`)
- Energy tables (`.tsv`)
- For --pykvfinder_box:
  - Grid box visualization PDBs
  - PyKVFinder pocket PDBs

## Contributing & Issues
While best efforts have been made to ensure accuracy and adherence to standard practices, users are encouraged to evaluate the software for their specific use cases. If you encounter issues or wish to contribute improvements, please open an issue or submit a pull request. 