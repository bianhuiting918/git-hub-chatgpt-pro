#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np

MOTIF = [40,49,56]

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--status-tsv', required=True)
    ap.add_argument('--reference-pdb', required=True)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--summary-json', required=True)
    ap.add_argument('--accepted-tsv', required=True)
    ap.add_argument('--motif-ca-rmsd-cutoff', type=float, default=1.0)
    ap.add_argument('--accepted-per-bin', type=int, default=10)
    return ap.parse_args()

def read_tsv(path: Path):
    if not path.exists():
        return []
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write_tsv(path: Path, rows: list[dict], fields: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def parse_fixed_positions(text: str) -> list[int]:
    out = []
    for tok in text.split():
        tok = tok.strip()
        if tok.startswith('A') and tok[1:].isdigit():
            out.append(int(tok[1:]))
    return sorted(set(out))

def parse_ca(path: Path):
    ca = {}
    if not path.exists():
        return ca
    for line in path.open():
        if not line.startswith('ATOM'):
            continue
        if line[12:16].strip() != 'CA':
            continue
        chain = line[21].strip() or 'A'
        if chain != 'A':
            continue
        resid = int(line[22:26])
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except Exception:
            parts = line.split(); xyz = (float(parts[6]), float(parts[7]), float(parts[8]))
        ca[resid] = xyz
    return ca

def rmsd(points_a, points_b):
    ref = []
    cand = []
    missing = []
    for resid in MOTIF:
        if resid not in points_a or resid not in points_b:
            missing.append(f'A{resid}')
            continue
        ref.append(points_a[resid])
        cand.append(points_b[resid])
    if missing or not ref:
        return '', ' '.join(missing)
    a = np.asarray(ref, dtype=float)
    b = np.asarray(cand, dtype=float)
    ac = a - a.mean(axis=0)
    bc = b - b.mean(axis=0)
    h = bc.T @ ac
    u, s, vt = np.linalg.svd(h)
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt
    b_aligned = bc @ r
    diff = ac - b_aligned
    return round(float(np.sqrt((diff * diff).sum() / len(ref))), 4), ''

def main() -> int:
    args = parse_args()
    manifest_rows = read_tsv(Path(args.manifest))
    status_rows = read_tsv(Path(args.status_tsv))
    latest_status = {}
    for row in status_rows:
        latest_status[row['sample_id']] = row
    ref_ca = parse_ca(Path(args.reference_pdb))
    out_rows = []
    for row in manifest_rows:
        sample = row['sample_id']
        status = latest_status.get(sample, {})
        pred = Path(row['postseq_predicted_pdb'])
        fixed = parse_fixed_positions(row.get('fixed_positions',''))
        seq = row['sequence']
        fixed_mut = []
        # Identity reference is encoded by generation: fixed positions were never mutated.
        # This check verifies the generated sequence still matches its own manifest fixed mask.
        for mut in row.get('mutations','').split(';'):
            if not mut:
                continue
            pos = int(mut.split(':',1)[0][1:])
            if pos in fixed:
                fixed_mut.append(mut)
        candidate_ca = parse_ca(pred) if pred.exists() and pred.stat().st_size > 1000 else {}
        motif_rmsd, missing = rmsd(ref_ca, candidate_ca)
        if status.get('status') in {'OK','SKIP_EXISTS'} and not fixed_mut and motif_rmsd != '' and motif_rmsd <= args.motif_ca_rmsd_cutoff:
            gate = 'PASS'
            reasons = ''
        else:
            gate = 'FAIL' if status.get('status') in {'OK','SKIP_EXISTS'} else 'NOT_EVALUATED'
            reasons_list = []
            if status.get('status') not in {'OK','SKIP_EXISTS'}:
                reasons_list.append('esmfold_not_ok')
            if fixed_mut:
                reasons_list.append('fixed_position_mutation')
            if missing:
                reasons_list.append('missing_motif_ca')
            if motif_rmsd != '' and motif_rmsd > args.motif_ca_rmsd_cutoff:
                reasons_list.append('motif_ca_rmsd_gt_1A')
            reasons = ';'.join(reasons_list)
        out = dict(row)
        out.update({
            'esmfold_status': status.get('status','NOT_EVALUATED'),
            'esmfold_bytes': status.get('bytes',''),
            'motif_ca_rmsd_A': motif_rmsd,
            'missing_motif_ca': missing,
            'fixed_mutation_count': len(fixed_mut),
            'fixed_mutations': ';'.join(fixed_mut),
            'routeB_active_pocket_only_gate': gate,
            'gate_fail_reasons': reasons,
            'global_ca_rmsd_role': 'NOT_USED_HARD_GATE',
            'sidechain_ligand_geometry_role': 'NOT_USED_AT_SEQUENCE_STAGE_ESMFOLD_APO',
        })
        out_rows.append(out)
    fields = list(out_rows[0].keys()) if out_rows else []
    write_tsv(Path(args.out_tsv), out_rows, fields)
    accepted = [r for r in out_rows if r['routeB_active_pocket_only_gate'] == 'PASS']
    accepted.sort(key=lambda r: (int(r['bin']), float(r['motif_ca_rmsd_A'])))
    selected = []
    by_bin = {}
    status_counts = {}
    for r in out_rows:
        b = r['bin']; g = r['routeB_active_pocket_only_gate']
        by_bin.setdefault(b, {'evaluated':0,'pass':0,'selected_for_target':0})
        by_bin[b]['evaluated'] += 1
        if g == 'PASS': by_bin[b]['pass'] += 1
        status_counts[g] = status_counts.get(g,0)+1
    selected_count = {}
    for r in accepted:
        b = r['bin']
        if selected_count.get(b,0) >= args.accepted_per_bin:
            continue
        selected.append(r)
        selected_count[b] = selected_count.get(b,0) + 1
    for b,n in selected_count.items():
        by_bin[b]['selected_for_target'] = n
    write_tsv(Path(args.accepted_tsv), selected, fields)
    summary = {
        'manifest': args.manifest,
        'status_tsv': args.status_tsv,
        'reference_pdb': args.reference_pdb,
        'out_tsv': args.out_tsv,
        'accepted_tsv': args.accepted_tsv,
        'evaluated_universe_rows': len(out_rows),
        'gate_definition': 'PASS if ESMFold OK/SKIP_EXISTS, fixed positions unchanged, motif CA RMSD over A40/A49/A56 <= 1.0 A',
        'global_ca_rmsd_role': 'not a hard gate',
        'sidechain_ligand_geometry_role': 'diagnostic only until ligand-aware holo stage',
        'counts': status_counts,
        'by_bin': by_bin,
        'target_per_bin': args.accepted_per_bin,
        'bins_meeting_target': sorted([b for b,v in by_bin.items() if v.get('selected_for_target',0) >= args.accepted_per_bin], key=int),
        'target_complete': all(by_bin.get(str(b),{}).get('selected_for_target',0) >= args.accepted_per_bin for b in [90,80,70,60,50]),
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
