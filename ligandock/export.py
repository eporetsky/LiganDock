import os
import re
import subprocess
import math

def extract_model_and_rename_ligands(input_pdbqt, pose_num, ligand_names, output_folder):
    with open(input_pdbqt, 'r') as f:
        lines = f.readlines()
    # Find the MODEL (POSE) X block
    model_start = None
    model_end = None
    for i, line in enumerate(lines):
        if line.startswith(f"MODEL {pose_num}"):
            model_start = i
        elif model_start is not None and line.startswith("MODEL") and not line.startswith(f"MODEL {pose_num}"):
            model_end = i
            break
    if model_start is None:
        raise ValueError(f"MODEL {pose_num} not found in {input_pdbqt}")
    if model_end is None:
        model_end = len(lines)
    model_lines = lines[model_start:model_end]
    # Find ligands (ROOT...REMARK SMILES)
    ligand_blocks = []
    block = []
    in_ligand = False
    for line in model_lines:
        if line.startswith("ROOT"):
            in_ligand = True
            block = [line]
        elif in_ligand:
            block.append(line)
            if line.startswith("REMARK SMILES"):
                ligand_blocks.append(block)
                in_ligand = False
                block = []
        else:
            continue
    # If the last ligand block didn't end with REMARK SMILES, add it
    if in_ligand and block:
        ligand_blocks.append(block)
    ligand_names_list = [name.strip() for name in ligand_names.split(",")]
    if len(ligand_blocks) != len(ligand_names_list):
        raise ValueError(f"Number of ligands in model ({len(ligand_blocks)}) does not match number of ligand names ({len(ligand_names_list)})")
    pdbqt_files = []
    sdf_files = []
    for block, lig_name in zip(ligand_blocks, ligand_names_list):
        # Rename UNL to ligand name in all relevant lines
        new_block = []
        for line in block:
            if (line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("ROOT") or line.startswith("BRANCH") or line.startswith("ENDROOT") or line.startswith("ENDBRANCH")) and " UNL " in line:
                new_block.append(line.replace(" UNL ", f" {lig_name} "))
            else:
                new_block.append(line)
        pdbqt_file = os.path.join(output_folder, f"{lig_name}.pdbqt")
        with open(pdbqt_file, 'w') as f:
            f.write(f"MODEL {pose_num}\n")
            for l in new_block:
                f.write(l)
        fix_pdbqt_element_columns(pdbqt_file)
        pdbqt_files.append(pdbqt_file)
        sdf_file = os.path.join(output_folder, f"{lig_name}.sdf")
        print(f"obabel {pdbqt_file} -O {sdf_file}")
        subprocess.run(f"obabel {pdbqt_file} -O {sdf_file}", shell=True, check=True)
        sdf_files.append(sdf_file)
    return pdbqt_files, sdf_files

def fix_pdbqt_element_columns(pdbqt_file):
    import re
    lines = []
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16].strip()
                element = re.sub('[^A-Za-z]', '', atom_name)
                if len(element) > 1 and element[0] == 'H':
                    element = 'H'
                elif len(element) > 1:
                    element = element[:2].title()
                elif element:
                    element = element[0].upper()
                else:
                    element = 'X'
                line = line[:76] + f"{element:>2}" + line[78:]
            lines.append(line)
    with open(pdbqt_file, 'w') as f:
        f.writelines(lines)

def merge_sdf_files(sdf_files, merged_sdf):
    with open(merged_sdf, 'w') as outfile:
        for sdf_file in sdf_files:
            with open(sdf_file, 'r') as infile:
                outfile.write(infile.read())

def pdbqt_ligand_to_pdb_lines(pdbqt_file, start_atom_num=1, chain_id='L'):
    pdb_lines = []
    atom_num = start_atom_num
    res_num = 1
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16]
                res_name = line[17:20]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                occ = line[54:60] if len(line) >= 60 else '  1.00'
                temp = line[60:66] if len(line) >= 66 else ' 20.00'
                element = line[76:78].strip() if len(line) >= 78 else atom_name.strip()[0]
                pdb_line = f"HETATM{atom_num:5d} {atom_name:<4} {res_name:>3} {chain_id}{res_num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occ:>6}{temp:>6}          {element:>2}\n"
                pdb_lines.append(pdb_line)
                atom_num += 1
    return pdb_lines, atom_num-1

def merge_receptor_and_ligands_to_pdb(receptor_pdb, ligand_pdbqt_files, output_pdb):
    # Read receptor atoms up to and including the first TER line
    receptor_lines = []
    with open(receptor_pdb, 'r') as f:
        for line in f:
            if line.startswith('TER'):
                receptor_lines.append('TER\n')
                break
            receptor_lines.append(line)
    # Find last atom number in receptor
    last_atom_num = 0
    for line in receptor_lines:
        if line.startswith(('ATOM', 'HETATM')):
            try:
                last_atom_num = max(last_atom_num, int(line[6:11]))
            except Exception:
                continue
    # Add ligands, continuing atom numbers
    merged_lines = receptor_lines[:]
    chain_id = 'L'
    for i, lig_file in enumerate(ligand_pdbqt_files):
        lig_lines, lig_last_atom = pdbqt_ligand_to_pdb_lines(lig_file, start_atom_num=last_atom_num+1, chain_id=chr(ord('L')+i))
        merged_lines.extend(lig_lines)
        last_atom_num = lig_last_atom  # Continue atom numbering for next ligand
    merged_lines.append('END\n')
    with open(output_pdb, 'w') as f:
        f.writelines(merged_lines)
