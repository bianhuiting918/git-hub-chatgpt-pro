import pathlib
import sys

FLOW = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FLOW / "scripts"))

from prepare_nylc_nalpha_h2_ash306 import drop_gro_atom, export_chain_pdb, extract_molecule_itp, replace_gro_slice, rewrite_chain_itp

ITP = """[ moleculetype ]
Protein_chain_H 3

[ atoms ]
; nr type resnr residue atom cgnr charge mass
; residue 1 THR rtp NTHR q +1.0
1 N3 1 THR N 1 0.18120 14.01000
2 H 1 THR H1 2 0.19340 1.00800
3 H 1 THR H2 3 0.19340 1.00800
4 H 1 THR H3 4 0.19340 1.00800
5 CT 1 THR CA 5 0.00340 12.01000
6 C 1 THR C 6 0.61630 12.01000 ; qtot 1

[ bonds ]
1 2 1
1 3 1
1 4 1
1 5 1
5 6 1

[ angles ]
2 1 4 1
2 1 5 1

[ dihedrals ]
2 1 5 6 9
4 1 5 6 9

[ cmap ]
2 1 5 6 5 1
"""


def test_rewrite_chain_itp_removes_h3_and_exactly_neutralizes_nalpha():
    rewritten, audit = rewrite_chain_itp(ITP, remove_index1=4, nalpha_charge=-0.62540)
    assert " THR H3 " not in rewritten
    assert audit["removed_atom"] == {"index1": 4, "name": "H3", "charge_e": 0.1934}
    assert abs(audit["charge_delta_e"] + 1.0) < 1e-8
    assert audit["new_atom_count"] == 5
    assert audit["dropped_interaction_counts"] == {"bonds": 1, "angles": 1, "dihedrals": 1}
    assert "1 4 1" in rewritten
    assert "1 5 1" not in rewritten
    assert "4 1 5 6" not in rewritten
    assert "NTHR_NH2_PROXY q 0.0" in rewritten
    assert "qtot 0" in rewritten


def test_drop_gro_atom_preserves_box_and_renumbers_atoms():
    gro = """fixture
5
    1THR      N    1   0.000   0.000   0.000
    1THR     H1    2   0.010   0.000   0.000
    1THR     H2    3   0.000   0.010   0.000
    1THR     H3    4   0.000   0.000   0.010
    1THR     CA    5   0.100   0.000   0.000
   9.00000   9.00000   9.00000
"""
    rewritten, audit = drop_gro_atom(gro, remove_index1=4, expected_name="H3")
    lines = rewritten.splitlines()
    assert lines[1].strip() == "4"
    assert "H3" not in rewritten
    assert lines[-1] == "   9.00000   9.00000   9.00000"
    assert int(lines[-2][15:20]) == 4
    assert audit["old_atom_count"] == 5
    assert audit["new_atom_count"] == 4


def test_extract_molecule_itp_strips_standalone_forcefield_water_and_system():
    standalone = """; header
#include \"amber99sb-ildn.ff/forcefield.itp\"

[ moleculetype ]
Protein_chain_H 3
[ atoms ]
1 N3 267 THR N 1 0.0 14.01
#ifdef POSRES
#include \"posre_chainH.itp\"
#endif

; Include water topology
#include \"amber99sb-ildn.ff/tip3p.itp\"
[ system ]
Protein
[ molecules ]
Protein_chain_H 1
"""
    extracted = extract_molecule_itp(standalone)
    assert extracted.startswith("[ moleculetype ]")
    assert "posre_chainH.itp" in extracted
    assert "forcefield.itp" not in extracted
    assert "tip3p.itp" not in extracted
    assert "[ system ]" not in extracted

def test_replace_gro_slice_requires_same_chain_length_and_preserves_suffix():
    full = """full
5
    1AAA      A    1   0.000   0.000   0.000
    2BBB      B    2   0.100   0.000   0.000
    2BBB      C    3   0.200   0.000   0.000
    3LIG      X    4   0.300   0.000   0.000
    3LIG      Y    5   0.400   0.000   0.000
   4.00000   4.00000   4.00000
"""
    chain = """chain
2
  267THR      N    1   1.100   0.000   0.000
  267THR     CA    2   1.200   0.000   0.000
   4.00000   4.00000   4.00000
"""
    merged, audit = replace_gro_slice(full, chain, start_index1=2, end_index1=3)
    lines = merged.splitlines()
    assert lines[1].strip() == "5"
    assert "267THR" in lines[3]
    assert "3LIG" in lines[5]
    assert int(lines[5][15:20]) == 4
    assert int(lines[6][15:20]) == 5
    assert audit["replaced_atom_count"] == 2
    assert audit["full_atom_count"] == 5


def test_export_chain_pdb_marks_only_asp306_as_ash():
    full = """full
5
  267THR      N    1   0.000   0.000   0.000
  306ASP     CG    2   0.100   0.200   0.300
  306ASP    OD1    3   0.200   0.200   0.300
  306ASP    OD2    4   0.100   0.300   0.300
  307ALA     CA    5   0.500   0.500   0.500
   4.00000   4.00000   4.00000
"""
    pdb, audit = export_chain_pdb(full, start_index1=1, end_index1=5, ash_resid=306)
    assert pdb.count(" ASH H 306") == 3
    assert " THR H 267" in pdb
    assert " ALA H 307" in pdb
    assert audit["chain_atom_count"] == 5
    assert audit["renamed_ash306_atom_count"] == 3
