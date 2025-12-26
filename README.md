# pdb-cli

A Python CLI utility to extract and summarize useful information from PDB (Protein Data Bank) and mmCIF files.

## Features
- Number of chains and chain names
- Sequence extraction with chain names
- Residue composition
- Number of atoms
- Number of unique residues
- Bounding box calculation
- All features can be toggled with command-line flags
- Optional colored output using [`rich`](https://github.com/Textualize/rich)

## Installation

```bash
pip install .
```

Or, for development:

```bash
git clone https://github.com/a-r-j/pdb-cli.git
cd pdb-cli
pip install -e .
```

## Usage

```bash
pdb-cli download 4hhb
pdb-cli myfile.pdb [flags]
pdb-cli myfile.cif [flags]
```

```
📊 Summary for 4hhb.pdb:
====================================
🔹 Chains: 
  Count: 4
  Chain IDs: ['A', 'B', 'C', 'D']
  Chain lengths: {'Chain A': 198, 'Chain B': 205, 'Chain C': 201, 'Chain D': 197}
--------------------------------------------------
🔹 Sequence(s): 
  Chain A: VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR
  Chain B: VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH
  Chain C: VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR
  Chain D: VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH
====================================
```

### Example

```bash
pdb-cli summary 4hhb.pdb --chains --sequence --composition --residues --atoms --unique-residues --bbox --ligands --ptm --missing-residues --bfactors --binding-sites --interface --color
pdb-cli summary 4hhb.cif --chains --sequence --composition --residues --atoms --unique-residues --bbox --ligands --ptm --missing-residues --bfactors --binding-sites --interface --color
```
```

📊 Summary for 4hhb.pdb:
====================================
🔹 Chains: 
  Count: 4
  Chain IDs: ['A', 'B', 'C', 'D']
  Chain lengths: {'Chain A': 198, 'Chain B': 205, 'Chain C': 201, 'Chain D': 197}
--------------------------------------------------
🔹 Sequence(s): 
  Chain A: VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR
  Chain B: VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH
  Chain C: VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR
  Chain D: VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH
--------------------------------------------------
🔹 Residue composition: 
  HOH: 221
  ALA: 72
  LEU: 72
  VAL: 62
  LYS: 44
  GLY: 40
  HIS: 38
  SER: 32
  THR: 32
  ASP: 30
  PHE: 30
  PRO: 28
  GLU: 24
  ASN: 20
  ARG: 12
  TYR: 12
  GLN: 8
  CYS: 6
  MET: 6
  TRP: 6
  HEM: 4
  PO4: 2
--------------------------------------------------
🔹 Number of atoms: 4779
--------------------------------------------------
🔹 Number of unique residues: 801
--------------------------------------------------
🔹 Bounding box: 
  Min coordinates (x,y,z): (np.float32(-34.578), np.float32(-27.563), np.float32(-30.367))
  Max coordinates (x,y,z): (np.float32(35.602), np.float32(26.492), np.float32(30.048))
  Size (x,y,z): (np.float32(70.18), np.float32(54.055), np.float32(60.415))
--------------------------------------------------
🔹 Non-protein ligands: 
  HEM: 4
  PO4: 2
--------------------------------------------------
🔹 Post-translational modifications: No post-translational modifications detected
--------------------------------------------------
🔹 Missing residues: No missing residues detected
--------------------------------------------------
🔹 B-factors: 
  Mean B-factor: 24.8
  Median B-factor: 20.52
  Min B-factor: 4.91
  Max B-factor: 80.12
  Std deviation: 13.2
  High B-factor atoms (>50): 271
  Very high B-factor atoms (>80): 16
====================================
```
### Flags
- `--chains`           Show number of chains and chain names
- `--sequence`         Show sequence(s) with chain names
- `--composition`      Show residue composition
- `--residues`         Show all residues (chain:code:number)
- `--atoms`            Show number of atoms
- `--unique-residues`  Show number of unique residues
- `--bbox`             Show bounding box
- `--ligands`          Show non-protein ligands (small molecules, cofactors)
- `--ptm`              Show post-translational modifications
- `--missing-residues` Show missing residues in the structure
- `--bfactors`         Show B-factor statistics (quality indicator)
- `--binding-sites`    Show substrate binding sites near ligands
- `--binding-radius`   Radius (Å) for binding site analysis (default: 5.0)
- `--interface`        Show protein-protein interfaces between chains
- `--interface-radius` Radius (Å) for interface analysis (default: 5.0)
- `--color`            Enable colored output (requires `rich`)

## Downloading Structure Files

You can download structure files from various databases using the CLI:

### Download PDB or mmCIF files from RCSB PDB

```bash
# Download PDB file
pdb-cli download 1abc

# Download mmCIF file
pdb-cli download 1abc --cif

# Specify output directory
pdb-cli download 1abc --outdir ./structures
```

### Download AlphaFold predictions from AlphaFold Database

```bash
# Download AlphaFold prediction as PDB file
pdb-cli alphafold P12345

# Download AlphaFold prediction as mmCIF file
pdb-cli alphafold P12345 --cif

# Specify output directory
pdb-cli alphafold P12345 --outdir ./predictions
pdb-cli alphafold P12345 --cif --outdir ./predictions
```

The AlphaFold files are downloaded with the following formats:
- **PDB**: `AF-{uniprot_id}-F1-model_v6.pdb`
- **mmCIF**: `AF-{uniprot_id}-F1-model_v6.cif`

You can then use any downloaded file with the summary/analysis features as usual.

## License
MIT 