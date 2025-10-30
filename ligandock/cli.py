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
