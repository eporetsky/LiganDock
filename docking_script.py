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
    """Prepare receptor using Meeko or fallback to prepare_receptor4.py."""
    try:
        # First try pdb2pqr + Meeko
        if not use_fallback:
            print("  Trying Meeko receptor preparation...")
            os.system(f"pdb2pqr --ffout AMBER --keep-chain {prot_fl} tmp/receptor.pdbqr 2>/dev/null")
            if os.path.exists("tmp/receptor.pdbqr"):  # pdb2pqr succeeded
                os.system("mk_prepare_receptor.py --pdb tmp/receptor.pdbqr -o tmp/receptor.pdbqt --skip_gpf 2>/dev/null")
                if os.path.exists("tmp/receptor.pdbqt"):
                    print("  Meeko receptor preparation successful")
                    return True
        
        # Fallback: try prepare_receptor4.py directly on PDB
        print("  Trying prepare_receptor4.py fallback...")
        os.system(f"prepare_receptor4.py -r {prot_fl} -o tmp/receptor.pdbqt 2>/dev/null")
        if os.path.exists("tmp/receptor.pdbqt"):
            print("  prepare_receptor4.py successful")
            return True
            
        # If both fail, try OpenBabel conversion + prepare_receptor4.py
        print("  Trying OpenBabel + prepare_receptor4.py...")
        # Convert PDB to MOL2 using OpenBabel (sometimes cleans up formatting)
        obconversion = openbabel.OBConversion()
        obconversion.SetInAndOutFormats("pdb", "mol2")
        mol = openbabel.OBMol()
        if obconversion.ReadFile(mol, prot_fl):
            obconversion.WriteFile(mol, "tmp/temp_receptor.mol2")
            # Convert back to PDB
            obconversion.SetInAndOutFormats("mol2", "pdb")
            mol2 = openbabel.OBMol()
            if obconversion.ReadFile(mol2, "tmp/temp_receptor.mol2"):
                obconversion.WriteFile(mol2, "tmp/temp_receptor.pdb")
                os.system("prepare_receptor4.py -r tmp/temp_receptor.pdb -o tmp/receptor.pdbqt 2>/dev/null")
                if os.path.exists("tmp/receptor.pdbqt"):
                    print("  OpenBabel + prepare_receptor4.py successful")
                    return True
        
        return False
        
    except Exception as e:
        print(f"  Error in receptor preparation: {str(e)}")
        return False

def cleanup_temp_files():
    """Remove temporary files."""
    if os.path.exists("tmp"):
        shutil.rmtree("tmp")
    os.makedirs("tmp", exist_ok=True)

def main():
    parser = argparse.ArgumentParser(description='Run AutoDock Vina docking with SMILES ligand and PDB proteins.')
    parser.add_argument('--smiles', required=True, nargs='+', help='SMILES string(s) of the ligand(s)')
    parser.add_argument('--pdb_folder', required=True, help='Folder containing PDB files')
    parser.add_argument('--marker_file', required=True, help='PDB or SDF file containing marker molecule for docking box center')
    parser.add_argument('--output_folder', required=True, help='Output folder for docking results')
    parser.add_argument('--ligand_name', nargs='*', help='Name(s) of the ligand(s) (optional, will be added to output filenames)')
    parser.add_argument('--exhaustiveness', type=int, default=16, help='Docking exhaustiveness (default: 1000)')
    parser.add_argument('--n_poses', type=int, default=20, help='Number of poses to generate (default: 20)')
    parser.add_argument('--box_size', type=float, nargs=3, default=[15, 15, 15], 
                      help='Docking box size in Angstroms (default: 12 12 12)')
    parser.add_argument('--use_fallback', action='store_true', 
                      help='Skip Meeko and use prepare_receptor4.py directly')
    
    args = parser.parse_args()
    
    # Validate ligand names if provided
    if args.ligand_name and len(args.ligand_name) != len(args.smiles):
        print(f"Error: Number of ligand names ({len(args.ligand_name)}) must match number of SMILES ({len(args.smiles)})")
        return
    
    # Create output folder and tmp folder
    os.makedirs(args.output_folder, exist_ok=True)
    os.makedirs("tmp", exist_ok=True)
    
    # Get marker centroid
    centroid = get_marker_centroid(args.marker_file)
    if centroid is None:
        print("Error: Could not process marker file. Exiting.")
        return
    
    print(f"Docking box center: {centroid}")
    
    # Prepare ligands
    ligand_data = []
    for i, smiles in enumerate(args.smiles):
        lig_pdbqt = prepare_ligand(smiles)
        ligand_name = args.ligand_name[i] if args.ligand_name else f"lig{i+1}"
        ligand_file = f"tmp/ligand_{ligand_name}.pdbqt"
        
        with open(ligand_file, 'w') as f:
            f.write(lig_pdbqt)
        
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
            cleanup_temp_files()
            
            # Create output filename
            output_name = os.path.join(args.output_folder, prot_id.replace(".pdb", f".{ligand_suffix}.pdbqt"))
                
            if output_name in complete_files:
                print(f"Skipping {prot_id} - already processed")
                continue
            
            # Prepare receptor with fallback options
            if not prepare_receptor(prot_fl, args.use_fallback):
                print(f"  Failed to prepare receptor for {prot_id}")
                continue
            
            # Recreate ligand files in clean tmp folder
            ligand_files = []
            for ligand in ligand_data:
                lig_pdbqt = prepare_ligand(ligand['smiles'])
                ligand_file = f"tmp/ligand_{ligand['name']}.pdbqt"
                
                with open(ligand_file, 'w') as f:
                    f.write(lig_pdbqt)
                
                ligand_files.append(ligand_file)
            
            # Run docking with multiple ligands
            v = Vina(sf_name='vina', verbosity=0)
            v.set_receptor("tmp/receptor.pdbqt")
            
            # Set multiple ligands
            #for ligand_file in ligand_files:
            v.set_ligand_from_file(ligand_files)
            
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