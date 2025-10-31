import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="prody")

import click
from .dock import run_docking

@click.group()
def cli():
    """Ligandock CLI entry point"""
    pass

@cli.command()
@click.option('--smiles', multiple=True, required=True, help='SMILES string(s) of the ligand(s)')
@click.option('--pdb', required=True, type=click.Path(exists=True), help='PDB file or folder with PDB files')
@click.option('--output_folder', required=True, type=click.Path(file_okay=False), help='Output folder for docking results')
@click.option('--ligand_name', multiple=True, help='Name(s) of the ligand(s) (optional, will be added to output filenames)')
@click.option('--exhaustiveness', default=32, show_default=True, type=int, help='Docking exhaustiveness')
@click.option('--n_poses', default=20, show_default=True, type=int, help='Number of poses to generate')
@click.option('--box_size', nargs=3, type=float, default=(15.0, 15.0, 15.0), show_default=True, help='Docking box size in Angstroms')
@click.option('--marker_file', type=click.Path(), help='PDB or SDF file containing marker molecule for docking box center')
@click.option('--res_box', help='Residue specifications for docking box center (e.g., "A123,A456")')
@click.option('--pykvfinder_box', help='Residue specification for PyKVFinder box (e.g., A202)')
@click.option('--flex_res', multiple=True, help='Flexible residue specifications (e.g., "A:315" "B:123")')
def dock(**kwargs):
    """Dock ligands against proteins."""
    run_docking(**kwargs)

@cli.command()
@click.option('--input', required=True, type=click.Path(exists=True), help='Input PDBQT file')
@click.option('--pose', required=True, type=int, help='Pose number to extract (was --model)')
@click.option('--ligand_names', required=True, help='Comma-separated ligand names in the order they appear')
@click.option('--output', required=True, type=click.Path(file_okay=False), help='Output folder (will be created if not exists)')
@click.option('--pdb', required=False, type=click.Path(exists=True), help='Receptor PDB file to merge with ligands')
def export(input, pose, ligand_names, output, pdb):
    """Extract a pose from multi-ligand PDBQT, export PDBQT/SDF, optionally merge with receptor."""
    from ligandock import export as export_module
    import os
    os.makedirs(output, exist_ok=True)
    pdbqt_files, sdf_files = export_module.extract_model_and_rename_ligands(input, pose, ligand_names, output)
    merged_sdf = os.path.join(output, 'merged.sdf')
    export_module.merge_sdf_files(sdf_files, merged_sdf)
    # Save merged PDB if --pdb is provided
    if pdb:
        basename = os.path.splitext(os.path.basename(input))[0]
        merged_pdb_path = os.path.join(output, f'{basename}.merged.pdb')
        export_module.merge_receptor_and_ligands_to_pdb(pdb, pdbqt_files, merged_pdb_path)
        print(f"Merged PDB written to {merged_pdb_path}")
    print("Generated files:")
    for fpath in pdbqt_files + sdf_files + [merged_sdf]:
        print(fpath)
