# NylC A1 representative NAC frame extraction design

Date: 2026-07-27

## Goal

Promote one auditable full explicit-water coordinate set from the selected
A1 replica `nac_evt25_time1462ps / seed26723`. The nominated time is 354 ps,
inside the 350--358 ps longest continuous NAC event and with the lowest
instantaneous classical-MM potential in that event.

## Scientific boundaries

- Microstate: Thr267 O-gamma-minus / N-alpha-H3+.
- Joint NAC: PA66-L2 reactive carbonyl C to Thr267 OG1 <=0.35 nm and carbonyl
  O-C-OG1 angle 95--115 degrees.
- Gate group: NylC residues 261--266; Thr267 is excluded.
- The full explicit-water/ion system is preserved.
- A selected classical-MM NAC coordinate is not a QM/MM GS optimization, TS,
  PMF, proton-transfer path, or activation barrier.

## Selected approach

GROMACS-native extraction plus an independent audit. This preserves trajectory
and PBC semantics while freezing provenance and preventing promotion of a
wrong-time or wrong-microstate frame.

## Data flow

1. Resolve the selected replica through its immutable A1 audit manifest.
2. Verify SHA256 of source TPR/XTC and confirm the requested time exists.
3. Extract only the 354 ps System frame to `source.tmp.gro` on SCNet.
4. Independently recompute atom identities, joint NAC geometry, gate opening,
   minimum heavy-atom contacts, atom count, box and finite coordinates.
5. Verify A1 proton inventory: Thr267 OG1 has no HG1 and N-alpha has three
   hydrogens; do not infer a proton pathway.
6. Run `grompp -maxwarn 0` with the existing A1 topology as a topology
   readiness check without dynamics.
7. Promote an immutable GRO and viewer PDB only after all gates pass. Write a
   manifest, SHA table and remote run-history terminal status.

## Outputs and acceptance

Coordinates stay on SCNet. Only compact JSON/TSV audit and documentation enter
GitHub. No trajectory, large topology, checkpoint, credential or secret is
committed. Technical PASS requires locked source hashes/time, complete finite
coordinates, correct A1 proton inventory, acceptable contacts and grompp
success. Scientific coordinate PASS additionally requires the joint NAC.
