# dock.py
"""Docking logic for ligandock 'dock' command."""

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="prody")
# ligandock/lib/python3.9/site-packages/prody/utilities/misctools.py:424: UserWarning: pkg_resources is deprecated as an API

import pandas as pd
import rdkit
import meeko
from meeko import PDBQTWriterLegacy
from vina import Vina
import os
import glob
import shutil
from openbabel import openbabel
import numpy as np
import subprocess

from ligandock.utils import pykvfinder_pdb, calculate_docking_box, generate_box_pdb

def prepare_ligand(smiles_string):
    lig = rdkit.Chem.MolFromSmiles(smiles_string)
    protonated_lig = rdkit.Chem.AddHs(lig)
    rdkit.Chem.AllChem.EmbedMolecule(protonated_lig)
    meeko_prep = meeko.MoleculePreparation()
    molecule_setups = meeko_prep.prepare(protonated_lig)
    writer = PDBQTWriterLegacy()
    result = writer.write_string(molecule_setups[0])
    if isinstance(result, tuple):
        return result[0]
    else:
        return result

def convert_to_sdf(input_file, output_file):
    obconversion = openbabel.OBConversion()
    obconversion.SetInAndOutFormats("pdb", "sdf")
    mol = openbabel.OBMol()
    if obconversion.ReadFile(mol, input_file):
        obconversion.WriteFile(mol, output_file)
        return True
    return False

def get_residue_centroid(pdb_file, residue_specs):
    residue_specs = residue_specs.split(',')
    try:
        all_coords = []
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        for res_spec in residue_specs:
            chain_id, res_num = res_spec[0], int(res_spec[1:])
            res_coords = []
            for line in lines:
                if line.startswith(('ATOM', 'HETATM')):
                    pdb_chain = line[21].strip()
                    pdb_res_num = int(line[22:26].strip())
                    if pdb_res_num == res_num:
                        if chain_id is None or pdb_chain == chain_id:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            res_coords.append([x, y, z])
            if not res_coords:
                print(f"Warning: No atoms found for residue {res_spec}")
                continue
            res_centroid = np.mean(res_coords, axis=0)
            all_coords.append(res_centroid)
            print(f"  Residue {res_spec} centroid: {res_centroid}")
        if not all_coords:
            print("Error: No valid residues found")
            return None
        geometric_center = np.mean(all_coords, axis=0)
        return geometric_center.tolist()
    except Exception as e:
        print(f"Error processing residues: {str(e)}")
        return None

def get_marker_centroid(marker_file):
    try:
        file_ext = os.path.splitext(marker_file)[1].lower()
        if file_ext == '.pdb':
            temp_sdf = "tmp/temp_marker.sdf"
            if not convert_to_sdf(marker_file, temp_sdf):
                print(f"Warning: Could not convert PDB file to SDF: {marker_file}")
                return None
            marker_file = temp_sdf
        marker = next(rdkit.Chem.SDMolSupplier(marker_file))
        if marker is None:
            print(f"Warning: Could not read marker file: {marker_file}")
            return None
        centroid = rdkit.Chem.rdMolTransforms.ComputeCentroid(marker.GetConformer())
        return [centroid.x, centroid.y, centroid.z]
    except Exception as e:
        print(f"Error processing marker file: {str(e)}")
        return None
    finally:
        if file_ext == '.pdb' and os.path.exists("tmp/temp_marker.sdf"):
            os.remove("tmp/temp_marker.sdf")

def prepare_receptor(prot_fl):
    try:
        print("  Trying mk_prepare_receptor.py receptor preparation...")
        print(f"mk_prepare_receptor.py --read_pdb {prot_fl} -o tmp/receptor -p -a")
        os.system(f"mk_prepare_receptor.py --read_pdb {prot_fl} -o tmp/receptor -p -a")
        if os.path.exists("tmp/receptor.pdbqt"):
            print("  mk_prepare_receptor.py receptor preparation successful")
            return True
        return False
    except Exception as e:
        print(f"  Error in receptor preparation: {str(e)}")
        return False

def prepare_flexible_receptor(prot_fl, flex_res_list, box_center, box_size, output_prefix="tmp/receptor"):
    try:
        flex_args = " ".join([f"-f {res}" for res in flex_res_list])
        cmd = (
            f"mk_prepare_receptor.py --read_pdb {prot_fl} -o {output_prefix} -p"
            f"--box_size {box_size[0]} {box_size[1]} {box_size[2]} "
            f"--box_center {box_center[0]} {box_center[1]} {box_center[2]} "
            f"{flex_args}"
        )
        print(f"  Running: {cmd}")
        subprocess.run(cmd, shell=True)
        rigid_pdbqt = f"{output_prefix}_rigid.pdbqt"
        flex_pdbqt = f"{output_prefix}_flex.pdbqt"
        box_txt = f"{output_prefix}.box.txt"
        if os.path.exists(rigid_pdbqt) and os.path.exists(flex_pdbqt) and os.path.exists(box_txt):
            return rigid_pdbqt, flex_pdbqt, box_txt
        else:
            print("  Flexible receptor preparation failed.")
            return None, None, None
    except Exception as e:
        print(f"  Error in flexible receptor preparation: {str(e)}")
        return None, None, None

def cleanup_temp_files():
    if os.path.exists("tmp"):
        shutil.rmtree("tmp")
    os.makedirs("tmp", exist_ok=True)

def prepare_ligand_with_scrubber(smiles_string, ligand_name, tmp_dir="tmp"):
    os.makedirs(tmp_dir, exist_ok=True)
    sdf_file = os.path.join(tmp_dir, f"{ligand_name}.sdf")
    scrub_cmd = f'scrub.py "{smiles_string}" -o {sdf_file} --ff "mmff94" --skip_tautomers --skip_acidbase --ph 7'
    print(f"Running: {scrub_cmd}")
    subprocess.run(scrub_cmd, shell=True, check=True)
    pdbqt_prefix = os.path.join(tmp_dir, f"{ligand_name}_protomer")
    mk_lig_cmd = f"mk_prepare_ligand.py -i {sdf_file} --multimol_prefix {pdbqt_prefix}"
    print(f"Running: {mk_lig_cmd}")
    subprocess.run(mk_lig_cmd, shell=True, check=True)
    pdbqt_file = os.path.join(tmp_dir, f"{ligand_name}_protomer-1.pdbqt")
    if os.path.exists(pdbqt_file):
        return pdbqt_file
    raise FileNotFoundError(f"Expected PDBQT file not found: {pdbqt_file}")

def run_docking(**kwargs):
    # Map kwargs as needed for compatibility
    smiles = kwargs.get('smiles')
    pdb_path = kwargs.get('pdb')
    output_folder = kwargs.get('output_folder')
    ligand_name = kwargs.get('ligand_name') or []
    exhaustiveness = kwargs.get('exhaustiveness')
    n_poses = kwargs.get('n_poses')
    box_size = tuple(kwargs.get('box_size') or (15.0, 15.0, 15.0))
    marker_file = kwargs.get('marker_file')
    res_box = kwargs.get('res_box') or []
    pykvfinder_box = kwargs.get('pykvfinder_box')
    flex_res = kwargs.get('flex_res') or []

    cleanup_temp_files()
    if ligand_name and len(ligand_name) != len(smiles):
        print(f"Error: Number of ligand names ({len(ligand_name)}) must match number of SMILES ({len(smiles)})")
        return
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs("tmp", exist_ok=True)
    # Prepare ligands
    ligand_data = []
    for i, smi in enumerate(smiles):
        name = ligand_name[i] if i < len(ligand_name) else f"lig{i+1}"
        ligand_file = prepare_ligand_with_scrubber(smi, name, tmp_dir="tmp")
        ligand_data.append({'name': name, 'file': ligand_file, 'smiles': smi})
    ligand_suffix = ".".join(ligand_name) if ligand_name else "multi_ligand"
    complete_files = set(glob.glob(os.path.join(output_folder, "*.pdbqt")))
    energies_df = pd.DataFrame()

    # Determine if pdb_path is file or dir
    pdb_files = []
    if os.path.isfile(pdb_path):
        if pdb_path.lower().endswith('.pdb'):
            pdb_files = [pdb_path]
        else:
            print(f"Error: --pdb must be a .pdb file or folder containing .pdb files")
            return
    elif os.path.isdir(pdb_path):
        pdb_files = sorted(glob.glob(os.path.join(pdb_path, "*.pdb")))
        if not pdb_files:
            print(f"Error: Directory {pdb_path} contains no .pdb files")
            return
    else:
        print(f"Error: --pdb must be an existing file or directory")
        return

    for prot_fl in pdb_files:
        try:
            prot_id = os.path.basename(prot_fl)
            print(f"Processing {prot_id}")
            if pykvfinder_box:
                chain_id = pykvfinder_box[0]
                residue_number = int(pykvfinder_box[1:])
                output_name = prot_id.replace(".pdb", ".pykvf.pdb")
                pykvf_pdb = os.path.join(output_folder, output_name)
                pykvfinder_pdb(chain_id, residue_number, prot_fl, pykvf_pdb)
                min_coords, max_coords, centroid, box_size = calculate_docking_box(pykvf_pdb, padding=0)
                if min_coords is None or max_coords is None:
                    print(f"  Error: Could not calculate docking box for {prot_id}. Skipping.")
                    continue
                box_pdb = os.path.join(output_folder, f"{prot_id}.box.pdb")
                generate_box_pdb(centroid, box_size, box_pdb)
            elif marker_file:
                centroid = get_marker_centroid(marker_file)
                if centroid is None:
                    print(f"  Error: Could not process marker file for {prot_id}. Skipping.")
                    continue
            else: # res_box
                centroid = get_residue_centroid(prot_fl, res_box)
                if centroid is None:
                    print(f"  Error: Could not process residues for {prot_id}. Skipping.")
                    continue
            print(f"  Docking box center for {prot_id}: {centroid}")
            print(f"  Docking box size for {prot_id}: {box_size}")
            output_name = os.path.join(output_folder, prot_id.replace(".pdb", f".{ligand_suffix}.pdbqt"))
            if output_name in complete_files:
                print(f"Skipping {prot_id} - already processed")
                continue
            if flex_res:
                rigid_pdbqt, flex_pdbqt, box_txt = prepare_flexible_receptor(
                    prot_fl, flex_res, centroid, box_size, output_prefix="tmp/receptor"
                )
                if not (rigid_pdbqt and flex_pdbqt and box_txt):
                    print(f"  Failed to prepare flexible receptor for {prot_id}")
                    continue
                receptor_file = rigid_pdbqt
                flex_file = flex_pdbqt
                config_file = box_txt
            else:
                if not prepare_receptor(prot_fl):
                    print(f"  Failed to prepare receptor for {prot_id}")
                    continue
                receptor_file = "tmp/receptor.pdbqt"
                flex_file = None
                config_file = None
            ligand_files = [lig['file'] for lig in ligand_data]
            v = Vina(sf_name='vina', verbosity=0)
            v.set_receptor(receptor_file)
            if flex_file:
                v.set_flex(flex_file)
            print(f"Setting ligands: {ligand_files}")
            v.set_ligand_from_file(ligand_files)
            if flex_file and config_file:
                v.compute_vina_maps(config=config_file)
            else:
                v.compute_vina_maps(center=centroid, box_size=box_size)
            v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
            v.write_poses(output_name, overwrite=True, n_poses=n_poses)
            energies = v.energies(n_poses=n_poses, energy_range=3.0)
            energy_data = []
            for i, energy_row in enumerate(energies):
                energy_dict = {
                    'protein_id': prot_id.replace('.pdb', ''),
                    'ligand_names': ligand_suffix,
                    'pose': i + 1,
                    'total_energy': energy_row[0],
                    'inter_energy': energy_row[1],
                    'intra_energy': energy_row[2],
                    'torsional_energy': energy_row[3],
                    'intra_best_pose': energy_row[4] if len(energy_row) > 4 else None
                }
                energy_data.append(energy_dict)
            protein_energies_df = pd.DataFrame(energy_data)
            energies_df = pd.concat([energies_df, protein_energies_df], ignore_index=True)
            print(f"  Completed docking for {prot_id}")
        except Exception as e:
            print(f"Error processing {prot_id}: {str(e)}")
        finally:
            pass
    cleanup_temp_files()
    if not energies_df.empty:
        energy_output_file = os.path.join(output_folder, prot_id.replace(".pdb", f".{ligand_suffix}.tsv"))
        energies_df.to_csv(energy_output_file, sep='\t', index=False)
        print(f"\nEnergy results saved to: {energy_output_file}")
        print("\nSummary statistics:")
        print(f"Total proteins processed: {energies_df['protein_id'].nunique()}")
        print(f"Total poses generated: {len(energies_df)}")
        print(f"Best total energy: {energies_df['total_energy'].min():.3f} kcal/mol")
        print(f"Mean total energy: {energies_df['total_energy'].mean():.3f} kcal/mol")
    else:
        print("\nNo energy data collected.")
