import pytest
import numpy as np
from Bio.PDB import Structure, Model, Chain, Residue, Atom

@pytest.fixture
def mock_structure():
    """Create a mock Bio.PDB structure for testing."""
    # Structure -> Model -> Chain -> Residue -> Atom
    structure = Structure.Structure("test_struct")
    model = Model.Model(0)
    structure.add(model)

    # Chain A: Protein
    chain_a = Chain.Chain("A")
    model.add(chain_a)
    
    # Residue 1: ALA
    res1 = Residue.Residue((" ", 1, " "), "ALA", " ")
    atom1 = Atom.Atom("CA", np.array([0.0, 0.0, 0.0], dtype=np.float32), 10.0, 1.0, " ", "CA", 1, "C")
    res1.add(atom1)
    chain_a.add(res1)
    
    # Residue 2: GLY (nearby)
    res2 = Residue.Residue((" ", 2, " "), "GLY", " ")
    atom2 = Atom.Atom("CA", np.array([3.0, 0.0, 0.0], dtype=np.float32), 20.0, 1.0, " ", "CA", 2, "C")
    res2.add(atom2)
    chain_a.add(res2)

    # Chain B: Protein (Interface with A)
    chain_b = Chain.Chain("B")
    model.add(chain_b)
    
    # Residue 1: VAL (close to Chain A res1)
    res3 = Residue.Residue((" ", 1, " "), "VAL", " ")
    atom3 = Atom.Atom("CA", np.array([0.0, 4.0, 0.0], dtype=np.float32), 30.0, 1.0, " ", "CA", 3, "C") # distance 4.0 from atom1
    res3.add(atom3)
    chain_b.add(res3)

    # Chain L: Ligand
    chain_l = Chain.Chain("L")
    model.add(chain_l)
    
    # Ligand: HEM (close to Chain A res2)
    ligand = Residue.Residue(("H_HEM", 1, " "), "HEM", " ")
    lig_atom = Atom.Atom("FE", np.array([3.0, 4.0, 0.0], dtype=np.float32), 15.0, 1.0, " ", "FE", 4, "FE") # distance 4.0 from atom2
    ligand.add(lig_atom)
    chain_l.add(ligand)

    return structure

@pytest.fixture
def mock_pdb_file(tmp_path):
    """Create a temporary PDB file."""
    p = tmp_path / "test.pdb"
    p.write_text("HEADER    TEST STRUCTURE\nATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C  \n")
    return str(p)
