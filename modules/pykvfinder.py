# conda activate lig

import os
import glob
import pandas as pd
import pyKVFinder
import warnings
from openbabel import pybel
import argparse
import numpy as np
warnings.filterwarnings('ignore')

from openbabel import pybel
def sdf_to_pdb(sdf_file, pdb_file):
    # Load the SDF file
    molecules = pybel.readfile("sdf", sdf_file)

    # Convert and save each molecule as PDB
    with open(pdb_file, "w") as outfile:
        for molecule in molecules:
            pdb_data = molecule.write("pdb")
            outfile.write(pdb_data)

residues_df = []

os.makedirs("pykvfinder", exist_ok=True)

#for fl in glob.glob("structures_aln/*.pdb"):
fl = "backup/H26B_NB.pdb"
file_name = os.path.basename(fl)

print("Starting to process:", file_name)

out_pdb = os.path.join("pykvfinder", file_name)

def pykvfinder_pdb(chain_id, residue_number, input_pdb, output_pdb):
    """Generate a PDB file for the specified residue using PyKVFinder and save to output_folder."""
    mol = next(pybel.readfile("pdb", input_pdb))
    residue_atoms = [atom for atom in mol if atom.OBAtom.GetResidue().GetNum() == residue_number and atom.OBAtom.GetResidue().GetChain() == chain_id]
    residue_mol = pybel.Molecule(pybel.ob.OBMol())
    for atom in residue_atoms:
        residue_mol.OBMol.AddAtom(atom.OBAtom)
    residue_mol.write("pdb", "tmp/pybel.pdb", overwrite=True)
    print(f"Generated PDB file for residue {chain_id}{residue_number} at {output_pdb}")

    results = pyKVFinder.run_workflow(input_pdb, ligand="tmp/pybel.pdb",
                volume_cutoff=100, probe_in=1.4, probe_out=5, step=0.5, ligand_cutoff=10, removal_distance=1.4)
    results.export_all(output=output_pdb)
    return output_pdb

def extract_coords_from_pdb(pdb_file):
    """
    Extract all atom coordinates from a PDB file, even if the file is not standard.
    Returns a numpy array of shape (N, 3).
    """
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except Exception:
                    continue
    return np.array(coords)

def add_box_padding(box_size, padding):
    """Add padding to box size in each dimension."""
    return box_size + padding

def calculate_docking_box(pdb_file, padding):
    """Return min_coords, max_coords, centroid, box_size for all atoms in pdb_file."""
    try:
        coords = extract_coords_from_pdb(pdb_file)
        if coords.shape[0] == 0:
            raise ValueError("No atoms found in PDB file.")
        min_coords = coords.min(axis=0)
        max_coords = coords.max(axis=0)
        centroid = (min_coords + max_coords) / 2
        box_size = max_coords - min_coords
        
        # Add padding to box size
        box_size = add_box_padding(box_size, padding)
        
        # Add debug output
        print(f"Found {coords.shape[0]} atoms in PDB file")
        print(f"Min coordinates: {min_coords}")
        print(f"Max coordinates: {max_coords}")
        print(f"Box center: {centroid}")
        print(f"Box size (with {padding}Å padding): {box_size}")
        
        return min_coords, max_coords, centroid, box_size
    except Exception as e:
        print(f"Error calculating docking box: {str(e)}")
        return None, None, None, None

def generate_box_pdb(centroid, box_size, output_pdb):
    """Generate a PDB file that visualizes the docking box using dummy atoms at the corners, centered at centroid, with given box_size (axis-aligned)."""
    try:
        # Calculate half-sizes
        half_size = np.array(box_size) / 2.0
        # 8 corners of the box
        corners = [
            centroid + np.array([ -half_size[0], -half_size[1], -half_size[2] ]),
            centroid + np.array([ -half_size[0], -half_size[1],  half_size[2] ]),
            centroid + np.array([ -half_size[0],  half_size[1], -half_size[2] ]),
            centroid + np.array([ -half_size[0],  half_size[1],  half_size[2] ]),
            centroid + np.array([  half_size[0], -half_size[1], -half_size[2] ]),
            centroid + np.array([  half_size[0], -half_size[1],  half_size[2] ]),
            centroid + np.array([  half_size[0],  half_size[1], -half_size[2] ]),
            centroid + np.array([  half_size[0],  half_size[1],  half_size[2] ]),
        ]
        with open(output_pdb, 'w') as f:
            for i, corner in enumerate(corners, start=1):
                f.write(f"HETATM{i:5d}  DUM DUM A   1    {corner[0]:8.3f}{corner[1]:8.3f}{corner[2]:8.3f}  1.00  0.00           D\n")
        print(f"Generated docking box PDB file at {output_pdb}")
    except Exception as e:
        print(f"Error generating docking box PDB: {str(e)}")