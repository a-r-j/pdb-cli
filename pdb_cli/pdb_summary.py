import os
from collections import Counter

from Bio.PDB import MMCIFParser, PDBParser


def summarize(
    pdb_file,
    chains=False,
    sequence=False,
    composition=False,
    residues=False,
    atoms=False,
    unique_residues=False,
    bbox=False,
    ligands=False,
    ptm=False,
    missing_residues=False,
    bfactors=False,
    binding_sites=False,
    binding_radius=5.0,
    interface=False,
    interface_radius=5.0,
):
    # Detect file format and use appropriate parser
    file_extension = os.path.splitext(pdb_file)[1].lower()
    if file_extension in [".cif", ".mmcif"]:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)
    result = {}

    if chains:
        result["Chains"] = get_chains_info(structure)
    if sequence:
        result["Sequence(s)"] = get_sequences_with_chains(structure)
    if composition:
        result["Residue composition"] = get_residue_composition(structure)
    if residues:
        result["Residues"] = get_all_residues(structure)
    if atoms:
        result["Number of atoms"] = count_atoms(structure)
    if unique_residues:
        result["Number of unique residues"] = count_unique_residues(structure)
    if bbox:
        result["Bounding box"] = get_bounding_box(structure)
    if ligands:
        result["Non-protein ligands"] = get_non_protein_ligands(structure)
    if ptm:
        result["Post-translational modifications"] = (
            get_post_translational_modifications(structure)
        )
    if missing_residues:
        result["Missing residues"] = get_missing_residues(structure)
    if bfactors:
        result["B-factors"] = get_bfactor_stats(structure)
    if binding_sites:
        result["Substrate binding sites"] = get_binding_sites(structure, binding_radius)
    if interface:
        result["Protein-protein interfaces"] = get_protein_interfaces(
            structure, interface_radius
        )
    return result


def get_chains_info(structure):
    chains = []
    chain_lengths = {}
    for model in structure:
        for chain in model:
            chains.append(chain.id)
            chain_lengths[f"Chain {chain.id}"] = len(list(chain))
    return {"Count": len(chains), "Chain IDs": chains, "Chain lengths": chain_lengths}


def get_sequences_with_chains(structure):
    from Bio.PDB.Polypeptide import PPBuilder

    ppb = PPBuilder()
    seqs = {}
    for model in structure:
        for chain in model:
            polypeptides = ppb.build_peptides(chain)
            if polypeptides:
                seqs[f"Chain {chain.id}"] = str(polypeptides[0].get_sequence())
            else:
                seqs[f"Chain {chain.id}"] = "N/A"
    return seqs if seqs else "N/A"


def get_residue_composition(structure):
    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                residues.append(residue.get_resname())
    composition = dict(Counter(residues))
    # Sort by count (descending) and then by name
    return dict(sorted(composition.items(), key=lambda x: (-x[1], x[0])))


def get_all_residues(structure):
    """Get all residues in the format chain:3letter_code:residue_number"""
    all_residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Get the author residue number (residue.id[1]) and 3-letter code
                residue_info = f"{chain.id}:{residue.get_resname()}:{residue.id[1]}"
                all_residues.append(residue_info)

    # Sort by chain, then by residue number
    all_residues.sort(key=lambda x: (x.split(":")[0], int(x.split(":")[2])))
    return all_residues


def count_atoms(structure):
    return sum(1 for _ in structure.get_atoms())


def count_unique_residues(structure):
    unique = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                unique.add((chain.id, residue.get_resname(), residue.id[1]))
    return len(unique)


def get_bounding_box(structure):
    coords = [atom.get_coord() for atom in structure.get_atoms()]
    if not coords:
        return "N/A"
    import numpy as np

    coords = np.array(coords)
    min_corner = coords.min(axis=0)
    max_corner = coords.max(axis=0)
    return {
        "Min coordinates (x,y,z)": tuple(min_corner.round(3)),
        "Max coordinates (x,y,z)": tuple(max_corner.round(3)),
        "Size (x,y,z)": tuple((max_corner - min_corner).round(3)),
    }


def get_non_protein_ligands(structure):
    """Identify non-protein ligands (small molecules, cofactors, etc.)"""
    # Standard amino acid residues (3-letter codes)
    protein_residues = {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "SEC",
        "PYL",  # Selenocysteine and Pyrrolysine
        "HID",
        "HIE",
        "HIP",  # Histidine variants
        "ASH",
        "GLH",  # Protonated forms
        "LYN",
        "CYM",
        "CYX",  # Other variants
    }

    ligands = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                # Skip water molecules
                if resname in ["HOH", "WAT", "TIP3", "SOL"]:
                    continue
                # Check if it's not a protein residue
                if resname not in protein_residues:
                    if resname not in ligands:
                        ligands[resname] = 0
                    ligands[resname] += 1

    if ligands:
        # Sort by count (descending) and then by name
        return dict(sorted(ligands.items(), key=lambda x: (-x[1], x[0])))
    else:
        return "No non-protein ligands found"


def get_post_translational_modifications(structure):
    """Identify post-translational modifications"""
    # Common PTM residue names
    ptm_residues = {
        "SEP": "Phosphoserine",
        "TPO": "Phosphothreonine",
        "PTR": "Phosphotyrosine",
        "HYP": "Hydroxyproline",
        "HIP": "Histidine variants",
        "CYX": "Disulfide-linked cysteine",
        "CSO": "S-hydroxycysteine",
        "CSD": "S-mercaptocysteine",
        "CSX": "Cysteine sulfenic acid",
        "CME": "S-carboxymethylcysteine",
        "CMH": "S-carboxymethylcysteine",
        "CSR": "S-carboxymethylcysteine",
        "CSS": "S-carboxymethylcysteine",
        "CSW": "S-carboxymethylcysteine",
        "CSY": "S-carboxymethylcysteine",
        "CSZ": "S-carboxymethylcysteine",
        "MSE": "Selenomethionine",
        "SEC": "Selenocysteine",
        "PYL": "Pyrrolysine",
        "KCX": "Lysine with unknown modification",
        "LLP": "Lysine with unknown modification",
        "MLY": "Methylated lysine",
        "MLZ": "Methylated lysine",
        "M3L": "Trimethylated lysine",
        "FME": "N-formylmethionine",
        "ACE": "N-acetylated residue",
        "NME": "N-methylated residue",
    }

    ptms = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in ptm_residues:
                    if resname not in ptms:
                        ptms[resname] = []
                    ptms[resname].append(f"{chain.id}:{residue.id[1]}")

    if ptms:
        # Sort by residue name
        return dict(sorted(ptms.items()))
    else:
        return "No post-translational modifications detected"


def get_missing_residues(structure):
    """Identify missing residues in the structure"""
    missing = {}
    for model in structure:
        for chain in model:
            residues = list(chain)
            if len(residues) < 2:
                continue

            # Sort residues by residue number
            residues.sort(key=lambda r: r.id[1])

            missing_in_chain = []
            for i in range(len(residues) - 1):
                current_num = residues[i].id[1]
                next_num = residues[i + 1].id[1]

                # Check for gaps in residue numbering
                if next_num - current_num > 1:
                    missing_range = list(range(current_num + 1, next_num))
                    missing_in_chain.extend(missing_range)

            if missing_in_chain:
                missing[f"Chain {chain.id}"] = missing_in_chain

    if missing:
        return missing
    else:
        return "No missing residues detected"


def get_bfactor_stats(structure):
    """Calculate B-factor statistics"""
    bfactors = []
    for atom in structure.get_atoms():
        bfactors.append(atom.get_bfactor())

    if not bfactors:
        return "No B-factor data available"

    import numpy as np

    bfactors = np.array(bfactors)

    return {
        "Mean B-factor": round(float(np.mean(bfactors)), 2),
        "Median B-factor": round(float(np.median(bfactors)), 2),
        "Min B-factor": round(float(np.min(bfactors)), 2),
        "Max B-factor": round(float(np.max(bfactors)), 2),
        "Std deviation": round(float(np.std(bfactors)), 2),
        "High B-factor atoms (>50)": int(np.sum(bfactors > 50)),
        "Very high B-factor atoms (>80)": int(np.sum(bfactors > 80)),
    }


def get_binding_sites(structure, radius=5.0):
    """Identify potential substrate binding sites near ligands"""
    # Get non-protein ligands
    protein_residues = {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "SEC",
        "PYL",
        "HID",
        "HIE",
        "HIP",
        "ASH",
        "GLH",
        "LYN",
        "CYM",
        "CYX",
    }

    # Find ligand atoms and organize by ligand type
    ligand_atoms_by_type = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname not in protein_residues and resname not in [
                    "HOH",
                    "WAT",
                    "TIP3",
                    "SOL",
                ]:
                    if resname not in ligand_atoms_by_type:
                        ligand_atoms_by_type[resname] = []
                    for atom in residue:
                        ligand_atoms_by_type[resname].append(atom)

    if not ligand_atoms_by_type:
        return "No ligands found for binding site analysis"

    # Find protein residues within radius of each ligand type
    binding_sites_by_ligand = {}

    for ligand_type, ligand_atoms in ligand_atoms_by_type.items():
        binding_residues = {}

        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname()
                    if resname in protein_residues:
                        for atom in residue:
                            min_distance = float("inf")
                            for ligand_atom in ligand_atoms:
                                distance = atom - ligand_atom
                                if distance <= radius:
                                    min_distance = min(min_distance, distance)

                            if min_distance != float("inf"):
                                key = f"{chain.id}:{residue.get_resname()}:{residue.id[1]}"
                                if key not in binding_residues:
                                    binding_residues[key] = {
                                        "Chain": chain.id,
                                        "Residue": resname,
                                        "Number": residue.id[1],
                                        "Min distance": min_distance,
                                    }

        if binding_residues:
            # Sort by distance
            sorted_residues = sorted(
                binding_residues.items(), key=lambda x: x[1]["Min distance"]
            )

            formatted_residues = {}
            for key, info in sorted_residues:
                formatted_residues[key] = f"Distance: {info['Min distance']:.2f}Å"

            binding_sites_by_ligand[ligand_type] = formatted_residues

    if binding_sites_by_ligand:
        # Create a formatted string with proper indentation
        output_lines = []
        output_lines.append(f"Binding radius used: {radius}Å")
        output_lines.append(f"Total ligands found: {len(binding_sites_by_ligand)}")
        output_lines.append("Ligand details:")

        for ligand_type, residues in binding_sites_by_ligand.items():
            output_lines.append(f"  Ligand: {ligand_type}")
            output_lines.append(f"    Number of interacting residues: {len(residues)}")
            output_lines.append("    Residues:")
            for residue_key, distance_info in residues.items():
                output_lines.append(f"      {residue_key}: {distance_info}")

        return "\n".join(output_lines)
    else:
        return f"No protein residues found within {radius}Å of ligands"


def get_protein_interfaces(structure, radius=5.0):
    """Identify protein-protein interfaces between different chains"""
    # Standard amino acid residues (3-letter codes)
    protein_residues = {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "SEC",
        "PYL",
        "HID",
        "HIE",
        "HIP",
        "ASH",
        "GLH",
        "LYN",
        "CYM",
        "CYX",
    }

    # Get all chains
    chains = []
    for model in structure:
        for chain in model:
            chains.append(chain)

    if len(chains) < 2:
        return "No protein-protein interfaces found (fewer than 2 chains)"

    # Find interfaces between different chains
    interfaces = {}
    interface_count = 0

    for i, chain1 in enumerate(chains):
        for j, chain2 in enumerate(chains[i + 1 :], i + 1):
            interface_residues = {}

            # Check residues from chain1 against chain2
            for residue1 in chain1:
                resname1 = residue1.get_resname()
                if resname1 in protein_residues:
                    for atom1 in residue1:
                        min_distance = float("inf")
                        closest_residue = None

                        for residue2 in chain2:
                            resname2 = residue2.get_resname()
                            if resname2 in protein_residues:
                                for atom2 in residue2:
                                    distance = atom1 - atom2
                                    if distance <= radius and distance < min_distance:
                                        min_distance = distance
                                        closest_residue = residue2

                        if min_distance != float("inf") and closest_residue is not None:
                            key1 = f"{chain1.id}:{resname1}:{residue1.id[1]}"
                            if key1 not in interface_residues:
                                interface_residues[key1] = {
                                    "Chain": chain1.id,
                                    "Residue": resname1,
                                    "Number": residue1.id[1],
                                    "Min distance": min_distance,
                                    "Interacting with": f"{chain2.id}:{closest_residue.get_resname()}:{closest_residue.id[1]}",
                                }

            # Check residues from chain2 against chain1
            for residue2 in chain2:
                resname2 = residue2.get_resname()
                if resname2 in protein_residues:
                    for atom2 in residue2:
                        min_distance = float("inf")
                        closest_residue = None

                        for residue1 in chain1:
                            resname1 = residue1.get_resname()
                            if resname1 in protein_residues:
                                for atom1 in residue1:
                                    distance = atom2 - atom1
                                    if distance <= radius and distance < min_distance:
                                        min_distance = distance
                                        closest_residue = residue1

                        if min_distance != float("inf") and closest_residue is not None:
                            key2 = f"{chain2.id}:{resname2}:{residue2.id[1]}"
                            if key2 not in interface_residues:
                                interface_residues[key2] = {
                                    "Chain": chain2.id,
                                    "Residue": resname2,
                                    "Number": residue2.id[1],
                                    "Min distance": min_distance,
                                    "Interacting with": f"{chain1.id}:{closest_residue.get_resname()}:{closest_residue.id[1]}",
                                }

            if interface_residues:
                interface_count += 1
                interface_name = f"Interface {chain1.id}-{chain2.id}"
                interfaces[interface_name] = interface_residues

    if interfaces:
        # Create a formatted string with proper indentation
        output_lines = []
        output_lines.append(f"Interface radius used: {radius}Å")
        output_lines.append(f"Total interfaces found: {interface_count}")
        output_lines.append("Interface details:")

        for interface_name, residues in interfaces.items():
            output_lines.append(f"  {interface_name}")
            output_lines.append(f"    Number of interface residues: {len(residues)}")
            output_lines.append("    Residues:")
            # Sort by distance (closest first)
            sorted_residues = sorted(
                residues.items(), key=lambda x: x[1]["Min distance"]
            )
            for residue_key, residue_info in sorted_residues:
                output_lines.append(
                    f"      {residue_key} -> {residue_info['Interacting with']} ({residue_info['Min distance']:.2f}Å)"
                )

        return "\n".join(output_lines)
    else:
        return f"No protein-protein interfaces found within {radius}Å"
