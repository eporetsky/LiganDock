import pandas as pd
import rdkit
import meeko
from meeko import PDBQTWriterLegacy
import scipy as sc
from vina import Vina
import os
import glob
import argparse
import shutil
from openbabel import openbabel
import numpy as np
from modules.pykvfinder import pykvfinder_pdb, calculate_docking_box, generate_box_pdb
import subprocess

def prepare_ligand(smiles_string):
    """Prepare ligand from SMILES string and return PDBQT string."""
    lig = rdkit.Chem.MolFromSmiles(smiles_string)
    protonated_lig = rdkit.Chem.AddHs(lig)
    rdkit.Chem.AllChem.EmbedMolecule(protonated_lig)
    
    meeko_prep = meeko.MoleculePreparation()
    molecule_setups = meeko_prep.prepare(protonated_lig)
    
    # Use the new API
    writer = PDBQTWriterLegacy()
    result = writer.write_string(molecule_setups[0])
    
    # Handle case where result might be a tuple
    if isinstance(result, tuple):
        return result[0]  # Take first element if it's a tuple
    else:
        return result

def convert_to_sdf(input_file, output_file):
    """Convert PDB to SDF format using OpenBabel."""
    obconversion = openbabel.OBConversion()
    obconversion.SetInAndOutFormats("pdb", "sdf")
    mol = openbabel.OBMol()
    
    if obconversion.ReadFile(mol, input_file):
        obconversion.WriteFile(mol, output_file)
        return True
    return False

def get_residue_centroid(pdb_file, residue_specs):
    """Get centroid coordinates from specified residues in PDB file."""
    try:
        all_coords = []
        
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        for res_spec in residue_specs:
            # Parse residue specification (e.g., "A:123" or "123" for any chain)
            if ':' in res_spec:
                chain_id, res_num = res_spec.split(':')
                res_num = int(res_num)
            else:
                chain_id = None
                res_num = int(res_spec)
            
            # Find atoms for this residue
            res_coords = []
            for line in lines:
                if line.startswith(('ATOM', 'HETATM')):
                    pdb_chain = line[21].strip()
                    pdb_res_num = int(line[22:26].strip())
                    
                    # Check if this atom belongs to our target residue
                    if pdb_res_num == res_num:
                        if chain_id is None or pdb_chain == chain_id:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            res_coords.append([x, y, z])
            
            if not res_coords:
                print(f"Warning: No atoms found for residue {res_spec}")
                continue
            
            # Add residue centroid to all coordinates
            res_centroid = np.mean(res_coords, axis=0)
            all_coords.append(res_centroid)
            print(f"  Residue {res_spec} centroid: {res_centroid}")
        
        if not all_coords:
            print("Error: No valid residues found")
            return None
        
        # Calculate geometric mean of all residue centroids
        geometric_center = np.mean(all_coords, axis=0)
        return geometric_center.tolist()
        
    except Exception as e:
        print(f"Error processing residues: {str(e)}")
        return None

def get_marker_centroid(marker_file):
    """Get centroid coordinates from marker file (PDB or SDF)."""
    try:
        # Check file extension
        file_ext = os.path.splitext(marker_file)[1].lower()
        
        # If PDB, convert to SDF first
        if file_ext == '.pdb':
            temp_sdf = "tmp/temp_marker.sdf"
            if not convert_to_sdf(marker_file, temp_sdf):
                print(f"Warning: Could not convert PDB file to SDF: {marker_file}")
                return None
            marker_file = temp_sdf
        
        # Read the SDF file
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
        # Clean up temporary SDF file if it was created
        if file_ext == '.pdb' and os.path.exists("tmp/temp_marker.sdf"):
            os.remove("tmp/temp_marker.sdf")

def prepare_receptor(prot_fl, use_fallback=False):
    """Prepare receptor using Meeko only."""
    try:
        # Use mk_prepare_receptor.py directly on the PDB file
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
    """
    Prepare rigid and flexible receptor PDBQT files using Meeko's mk_prepare_receptor.py.
    Returns (rigid_pdbqt, flex_pdbqt, box_txt) if successful, else (None, None, None).
    """
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
    """Remove temporary files."""
    if os.path.exists("tmp"):
        shutil.rmtree("tmp")
    os.makedirs("tmp", exist_ok=True)

def prepare_ligand_with_scrubber(smiles_string, ligand_name, tmp_dir="tmp"):
    """
    Prepare ligand using scrub.py and mk_prepare_ligand.py from a SMILES string.
    Returns the path to the generated PDBQT file (always {ligand_name}_protomer-1.pdbqt).
    """
    os.makedirs(tmp_dir, exist_ok=True)
    sdf_file = os.path.join(tmp_dir, f"{ligand_name}.sdf")
    # 1. Use scrub.py to generate SDF
    scrub_cmd = f'scrub.py "{smiles_string}" -o {sdf_file} --ff "mmff94" --skip_tautomers --skip_acidbase --ph 7'
    print(f"Running: {scrub_cmd}")
    subprocess.run(scrub_cmd, shell=True, check=True)
    # 2. Use mk_prepare_ligand.py to generate PDBQT
    pdbqt_prefix = os.path.join(tmp_dir, f"{ligand_name}_protomer")
    mk_lig_cmd = f"mk_prepare_ligand.py -i {sdf_file} --multimol_prefix {pdbqt_prefix}"
    print(f"Running: {mk_lig_cmd}")
    subprocess.run(mk_lig_cmd, shell=True, check=True)
    # Always use the -1 protomer file
    pdbqt_file = os.path.join(tmp_dir, f"{ligand_name}_protomer-1.pdbqt")
    if os.path.exists(pdbqt_file):
        return pdbqt_file
    raise FileNotFoundError(f"Expected PDBQT file not found: {pdbqt_file}")

def main():
    parser = argparse.ArgumentParser(description='Run AutoDock Vina docking with SMILES ligand and PDB proteins.')
    parser.add_argument('--smiles', required=True, nargs='+', help='SMILES string(s) of the ligand(s)')
    parser.add_argument('--pdb_folder', required=True, help='Folder containing PDB files')
    
    # Create mutually exclusive group for box center definition
    box_group = parser.add_mutually_exclusive_group(required=True)
    box_group.add_argument('--marker_file', help='PDB or SDF file containing marker molecule for docking box center')
    box_group.add_argument('--res_box', nargs='+', help='Residue specifications for docking box center (e.g., "A:123" "B:456" or "123" "456")')
    box_group.add_argument('--pykvfinder_box', help='Residue specification for PyKVFinder box (e.g., A202)')
    
    parser.add_argument('--output_folder', required=True, help='Output folder for docking results')
    parser.add_argument('--ligand_name', nargs='*', help='Name(s) of the ligand(s) (optional, will be added to output filenames)')
    parser.add_argument('--exhaustiveness', type=int, default=64, help='Docking exhaustiveness (default: 1000)')
    parser.add_argument('--n_poses', type=int, default=20, help='Number of poses to generate (default: 20)')
    parser.add_argument('--box_size', type=float, nargs=3, default=[15, 15, 15], 
                      help='Docking box size in Angstroms (default: 12 12 12)')
    parser.add_argument('--use_fallback', action='store_true', 
                      help='Skip Meeko and use prepare_receptor4.py directly')
    parser.add_argument('--flex_res', nargs='+', help='Flexible residue specifications (e.g., "A:315" "B:123")')
    
    args = parser.parse_args()
    
    cleanup_temp_files()

    # Validate ligand names if provided
    if args.ligand_name and len(args.ligand_name) != len(args.smiles):
        print(f"Error: Number of ligand names ({len(args.ligand_name)}) must match number of SMILES ({len(args.smiles)})")
        return
    
    # Create output folder and tmp folder
    os.makedirs(args.output_folder, exist_ok=True)
    os.makedirs("tmp", exist_ok=True)
    
    # Prepare ligands using new workflow
    ligand_data = []
    for i, smiles in enumerate(args.smiles):
        ligand_name = args.ligand_name[i] if args.ligand_name else f"lig{i+1}"
        ligand_file = prepare_ligand_with_scrubber(smiles, ligand_name, tmp_dir="tmp")
        ligand_data.append({
            'name': ligand_name,
            'file': ligand_file,
            'smiles': smiles
        })
    
    # Create ligand names string for output
    if args.ligand_name:
        ligand_suffix = ".".join(args.ligand_name)
    else:
        ligand_suffix = "multi_ligand"
    
    # Get list of completed files
    complete_files = set(glob.glob(os.path.join(args.output_folder, "*.pdbqt")))
    
    # Initialize DataFrame for energy results
    energies_df = pd.DataFrame()
    
    # Process each PDB file
    for prot_fl in glob.glob(os.path.join(args.pdb_folder, "*.pdb")):
        try:
            prot_id = os.path.basename(prot_fl)
            print(f"Processing {prot_id}")
            
            # Clean tmp folder for this iteration
            
            # Determine docking box center and size
            if args.pykvfinder_box:
                chain_id = args.pykvfinder_box[0]
                residue_number = int(args.pykvfinder_box[1:])
                pykvf_pdb = os.path.join(args.output_folder, f"{prot_id}.pykvf.pdb")
                pykvfinder_pdb(chain_id, residue_number, prot_fl, pykvf_pdb)
                min_coords, max_coords, centroid, box_size = calculate_docking_box(pykvf_pdb, padding=0)
                if min_coords is None or max_coords is None:
                    print(f"  Error: Could not calculate docking box for {prot_id}. Skipping.")
                    continue
                args.box_size = box_size
                box_pdb = os.path.join(args.output_folder, f"{prot_id}.box.pdb")
                generate_box_pdb(centroid, box_size, box_pdb)
            elif args.marker_file:
                centroid = get_marker_centroid(args.marker_file)
                box_size = args.box_size
                if centroid is None:
                    print(f"  Error: Could not process marker file for {prot_id}. Skipping.")
                    continue
            else:  # args.res_box
                centroid = get_residue_centroid(prot_fl, args.res_box)
                box_size = args.box_size
                if centroid is None:
                    print(f"  Error: Could not process residues for {prot_id}. Skipping.")
                    continue
            
            print(f"  Docking box center for {prot_id}: {centroid}")
            print(f"  Docking box size for {prot_id}: {box_size}")
        
            #break

            # Create output filename
            output_name = os.path.join(args.output_folder, prot_id.replace(".pdb", f".{ligand_suffix}.pdbqt"))
                
            if output_name in complete_files:
                print(f"Skipping {prot_id} - already processed")
                continue
            
            # If flexible residues are specified, prepare flexible receptor
            if args.flex_res:
                rigid_pdbqt, flex_pdbqt, box_txt = prepare_flexible_receptor(
                    prot_fl, args.flex_res, centroid, box_size, output_prefix="tmp/receptor"
                )
                if not (rigid_pdbqt and flex_pdbqt and box_txt):
                    print(f"  Failed to prepare flexible receptor for {prot_id}")
                    continue
                receptor_file = rigid_pdbqt
                flex_file = flex_pdbqt
                config_file = box_txt
            else:
                # Fallback to original receptor preparation
                if not prepare_receptor(prot_fl):
                    print(f"  Failed to prepare receptor for {prot_id}")
                    continue
                receptor_file = "tmp/receptor.pdbqt"
                flex_file = None
                config_file = None
            
            # Recreate ligand files in clean tmp folder
            ligand_files = [ligand['file'] for ligand in ligand_data]
            
            # Run docking with multiple ligands
            v = Vina(sf_name='vina', verbosity=0)
            v.set_receptor(receptor_file)
            if flex_file:
                v.set_flex(flex_file)
            print(f"Setting ligands: {ligand_files}")
            v.set_ligand_from_file(ligand_files)
            if flex_file and config_file:
                v.compute_vina_maps(config=config_file)
            else:
                v.compute_vina_maps(center=centroid, box_size=args.box_size)
            v.dock(exhaustiveness=args.exhaustiveness, n_poses=args.n_poses)
            v.write_poses(output_name, overwrite=True, n_poses=args.n_poses)
            
            # Extract energy data
            energies = v.energies(n_poses=args.n_poses, energy_range=3.0)
            
            # Convert to DataFrame and add metadata
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
            
            # Add to main DataFrame
            protein_energies_df = pd.DataFrame(energy_data)
            energies_df = pd.concat([energies_df, protein_energies_df], ignore_index=True)
            
            print(f"  Completed docking for {prot_id}")
            
        except Exception as e:
            print(f"Error processing {prot_id}: {str(e)}")
        finally:
            pass  # Keep tmp files for next iteration
    
    # Clean up tmp folder at the end
    cleanup_temp_files()
    
    # Save energy results
    if not energies_df.empty:
        energy_output_file = os.path.join(args.output_folder, prot_id.replace(".pdb", f".{ligand_suffix}.tsv"))
        energies_df.to_csv(energy_output_file, sep='\t', index=False)
        print(f"\nEnergy results saved to: {energy_output_file}")
        
        # Print summary statistics
        print("\nSummary statistics:")
        print(f"Total proteins processed: {energies_df['protein_id'].nunique()}")
        print(f"Total poses generated: {len(energies_df)}")
        print(f"Best total energy: {energies_df['total_energy'].min():.3f} kcal/mol")
        print(f"Mean total energy: {energies_df['total_energy'].mean():.3f} kcal/mol")
    else:
        print("\nNo energy data collected.")

if __name__ == "__main__":
    main() 