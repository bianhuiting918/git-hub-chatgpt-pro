#!/bin/bash
set -euo pipefail
if [[ $# -ne 4 ]]; then
  echo "usage: $0 PRMTOP RST7 OUTDIR LABEL" >&2
  exit 2
fi
PRMTOP=$1
RST7=$2
OUTDIR=$3
LABEL=$4
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
CPPTRAJ=$(command -v cpptraj)
test -s "$PRMTOP"
test -s "$RST7"
mkdir -p "$OUTDIR"
FULL="$OUTDIR/${LABEL}.full_explicit_water.pdb"
FOCUS="$OUTDIR/${LABEL}.active_site_8A.pdb"
cat > "$OUTDIR/${LABEL}.full.cpptraj.in" <<EOF
parm $PRMTOP
trajin $RST7 1 1
trajout $FULL pdb
run
EOF
cat > "$OUTDIR/${LABEL}.focus.cpptraj.in" <<EOF
parm $PRMTOP
reference $RST7
trajin $RST7 1 1
strip !(:L2<:8.0)
trajout $FOCUS pdb
run
EOF
"$CPPTRAJ" -i "$OUTDIR/${LABEL}.full.cpptraj.in" > "$OUTDIR/${LABEL}.full.cpptraj.log" 2>&1
"$CPPTRAJ" -i "$OUTDIR/${LABEL}.focus.cpptraj.in" > "$OUTDIR/${LABEL}.focus.cpptraj.log" 2>&1
test -s "$FULL"
test -s "$FOCUS"
sha256sum "$PRMTOP" "$RST7" "$FULL" "$FOCUS" > "$OUTDIR/${LABEL}.sha256.tsv"
{
  printf 'label\t%s\n' "$LABEL"
  printf 'source_prmtop\t%s\n' "$PRMTOP"
  printf 'source_rst7\t%s\n' "$RST7"
  printf 'focus_selection\tcomplete residues within 8 A of :L2\n'
  printf 'solvent\tfull PDB retains explicit water; focus PDB retains selected nearby residues/waters\n'
  printf 'scientific_status\tattack-only completed r01; 2.834598 A; not TS/PMF\n'
} > "$OUTDIR/${LABEL}.metadata.tsv"
stamp=$(date -Is)
printf '%s\texport_attack_seed_pdb\tPASS_TECHNICAL\tlabel=%s;source=%s;full=%s;focus=%s\n' "$stamp" "$LABEL" "$RST7" "$FULL" "$FOCUS" >> "$TASK/run_history.tsv"
printf '%s\n' "$FULL" "$FOCUS"
