import json, csv, subprocess, sys
sys.dont_write_bytecode = True
from pathlib import Path
from build_iccg_lg1_viewer import carbon_rings_from_pdb

ROOT=Path(__file__).parent
PDB=ROOT/'ICCG_active_cartoon_fixed_sites_with_LG1.pdb'
TSV=ROOT/'ICCG_active_fixed_sites.tsv'
AUDIT=ROOT/'ICCG_active_fixed_sites_audit.json'

def test_outputs_contract():
    assert PDB.exists() and TSV.exists() and AUDIT.exists()
    lines=PDB.read_text().splitlines()
    assert sum(1 for l in lines if l.startswith(('HELIX','SHEET'))) == 21
    assert not any('HEMT' in l or (l.startswith(('ATOM','HETATM')) and l[22:26].strip()=='856') for l in lines)
    atom_lines=[l for l in lines if l.startswith('ATOM')]
    residues={(l[21], int(l[22:26])) for l in atom_lines}
    assert len(residues)==258
    assert any(l.startswith('ATOM') and l[12:16].strip()=='OG' and l[17:20].strip()=='SER' and int(l[22:26])==165 for l in lines)
    lig=[l for l in lines if l.startswith('HETATM') and l[17:20].strip()=='LG1' and l[21]=='L' and int(l[22:26])==301]
    assert len(lig)==54
    assert sum(1 for l in lig if not l[76:78].strip().upper().startswith('H'))==32
    assert len(carbon_rings_from_pdb(PDB))==2
    conect=[l for l in lines if l.startswith('CONECT')]
    assert conect
    pairs=set()
    for l in conect:
        a=int(l[6:11]);
        for i in range(11,len(l),5):
            s=l[i:i+5].strip()
            if s: pairs.add((a,int(s)))
    assert all((b,a) in pairs for a,b in pairs)
    audit=json.loads(AUDIT.read_text())
    assert audit['geometry_status'] in {'PASS','PROVISIONAL_DISPLAY_ONLY_LG1_CLASH'}
    for k in ['worst_contact','min_nonbonded_protein_lg1_heavy_distance','max_vdw_overlap','alignment_anchor_residue_count','alignment_backbone_atom_count','backbone_fit_rmsd']:
        assert audit.get(k) not in (None,'',[])
    rows=list(csv.DictReader(TSV.open(), delimiter='\t'))
    fixed={int(r['resi']) for r in rows if r['fixed']=='1'}
    expected=set(audit['catalytic_residues']) | set(audit['new_LG1_5A_pocket_residues']) | set(audit['high_conservation_residues']) | set(audit['native_cysteine_residues'])
    assert fixed==expected==set(audit['fixed_residues'])
    assert audit['counts']['protein_residues']==258
    assert audit['counts']['lg1_atoms']==54
    assert audit['counts']['lg1_heavy_atoms']==32
    assert audit['counts']['fixed_residues']==len(fixed)
    assert audit['counts']['tsv_rows']==len(rows)

def test_anchor_fit_is_proper_kabsch_optimum():
    from build_iccg_lg1_viewer import compute_model, det3, rmsd_after_transform
    model = compute_model()
    audit = model['audit']
    rotation = model['rotation']
    assert abs(det3(rotation) - 1.0) < 1e-6
    recomputed = rmsd_after_transform(model['fit_reference_points'], model['fit_iccg_points'], rotation, model['translation'])
    assert abs(recomputed - audit['backbone_fit_rmsd']) < 1e-6
    assert audit['backbone_fit_rmsd'] < 2.0
    assert audit['alignment_key_mappings']['159'] == 164
    assert audit['alignment_key_mappings']['160'] == 165
    assert audit['alignment_key_mappings']['206'] == 210
    assert audit['alignment_key_mappings']['237'] == 242
    assert '238' in audit['alignment_key_mappings']



def test_b_factors_distances_counts_and_directory_contract():
    from build_iccg_lg1_viewer import compute_model, norm, vsub
    model = compute_model()
    audit = model['audit']
    lines=PDB.read_text().splitlines()
    fixed=set(audit['fixed_residues'])
    for l in lines:
        if l.startswith('ATOM'):
            resi=int(l[22:26]); b=float(l[60:66])
            assert abs(b-(0.0 if resi in fixed else 99.0)) < 0.01
        if l.startswith('HETATM') and l[17:20].strip()=='LG1':
            assert abs(float(l[60:66])-50.0) < 0.01
    orig=model['ligand_original']; trans=model['ligand_transformed']
    for i in range(len(orig)):
        for j in range(i+1,len(orig)):
            assert abs(norm(vsub(orig[i]['coord'],orig[j]['coord'])) - norm(vsub(trans[i]['coord'],trans[j]['coord']))) < 1e-6
    assert audit['counts']['protein_residues'] == len({int(l[22:26]) for l in lines if l.startswith('ATOM')})
    assert audit['counts']['lg1_atoms'] == sum(1 for l in lines if l.startswith('HETATM') and l[17:20].strip()=='LG1')
    assert audit['counts']['fixed_residues'] == sum(1 for r in csv.DictReader(TSV.open(), delimiter='\t') if r['fixed']=='1')
    allowed={'TASK.md','inputs','build_iccg_lg1_viewer.py','test_iccg_lg1_viewer.py','ICCG_active_cartoon_fixed_sites_with_LG1.pdb','ICCG_active_fixed_sites.tsv','ICCG_active_fixed_sites_audit.json','RUNBOOK.md','RUN_HISTORY.tsv'}
    assert {p.name for p in ROOT.iterdir()} <= allowed

if __name__ == '__main__':
    test_outputs_contract()
    test_anchor_fit_is_proper_kabsch_optimum()
    test_b_factors_distances_counts_and_directory_contract()
    print('all tests: PASS')
