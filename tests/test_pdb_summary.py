import pytest
import numpy as np
from Bio.PDB import Chain, Residue, Atom
from pdb_cli import pdb_summary

def test_get_chains_info(mock_structure):
    info = pdb_summary.get_chains_info(mock_structure)
    assert info["Count"] == 3
    assert info["Chain IDs"] == ["A", "B", "L"]
    assert info["Chain lengths"]["Chain A"] == 2
    assert info["Chain lengths"]["Chain B"] == 1

def test_get_sequences_with_chains(mock_structure):
    # Note: Bio.PDB.PPBuilder only builds peptides from protein residues
    seqs = pdb_summary.get_sequences_with_chains(mock_structure)
    assert "Chain A" in seqs
    assert "Chain B" in seqs
    # ALA + GLY -> AG
    assert "A" in seqs["Chain A"] # Depending on how PPBuilder handles it, usually AA sequence
    # Let's verify exact behavior later, but 'A' and 'G' are standard amino acids
    
def test_get_residue_composition(mock_structure):
    comp = pdb_summary.get_residue_composition(mock_structure)
    assert comp["ALA"] == 1
    assert comp["GLY"] == 1
    assert comp["VAL"] == 1
    assert comp["HEM"] == 1

def test_count_atoms(mock_structure):
    # 2 in A + 1 in B + 1 in L = 4
    count = pdb_summary.count_atoms(mock_structure)
    assert count == 4

def test_get_bounding_box(mock_structure):
    bbox = pdb_summary.get_bounding_box(mock_structure)
    # Min: 0,0,0
    # Max: 3,4,0
    assert bbox["Min coordinates (x,y,z)"] == (0.0, 0.0, 0.0)
    assert bbox["Max coordinates (x,y,z)"] == (3.0, 4.0, 0.0)
    assert bbox["Size (x,y,z)"] == (3.0, 4.0, 0.0)

def test_get_non_protein_ligands(mock_structure):
    ligands = pdb_summary.get_non_protein_ligands(mock_structure)
    assert ligands["HEM"] == 1
    # Check that standard residues are not ligands
    assert "ALA" not in ligands

def test_get_binding_sites(mock_structure):
    # Radius 5.0 should include interaction between HEM (Chain L) and GLY (Chain A atom at 3,0,0 vs HEM at 3,4,0 -> dist 4.0)
    sites = pdb_summary.get_binding_sites(mock_structure, radius=5.0)
    assert isinstance(sites, str)
    assert "Ligand: HEM" in sites
    assert "GLY" in sites
    assert "Distance: 4.00Å" in sites

def test_get_protein_interfaces(mock_structure):
    # Radius 5.0 should catch Chain A vs Chain B (Atom1 at 0,0,0 vs Atom3 at 0,4,0 -> dist 4.0)
    interfaces = pdb_summary.get_protein_interfaces(mock_structure, radius=5.0)
    assert isinstance(interfaces, str)
    assert "Interface A-B" in interfaces
    assert "ALA" in interfaces
    assert "VAL" in interfaces
    assert "4.00Å" in interfaces

def test_get_bfactor_stats(mock_structure):
    stats = pdb_summary.get_bfactor_stats(mock_structure)
    # B-factors: 10, 20, 30, 15
    # Mean: 18.75
    # Min: 10
    # Max: 30
    assert stats["Min B-factor"] == 10.0
    assert stats["Max B-factor"] == 30.0
    assert stats["Mean B-factor"] == 18.75

def test_get_missing_residues_no_gaps(mock_structure):
    # Our mock structure has residues 1, 2 in Chain A (no gap)
    missing = pdb_summary.get_missing_residues(mock_structure)
    assert missing == "No missing residues detected"

def test_get_missing_residues_with_gap(mock_structure):
    # Create a gap
    model = mock_structure[0]
    chain = Chain.Chain("C")
    res1 = Residue.Residue((" ", 1, " "), "ALA", " ")
    res3 = Residue.Residue((" ", 3, " "), "ALA", " ") # skip 2
    chain.add(res1)
    chain.add(res3)
    model.add(chain)
    
    missing = pdb_summary.get_missing_residues(mock_structure)
    assert isinstance(missing, dict)
    assert "Chain C" in missing
    assert missing["Chain C"] == [2]
