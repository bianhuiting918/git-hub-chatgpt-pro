# NylC A1 nine-replica NAC audit design

Date: 2026-07-27

## Scope

Audit the fully unrestrained 1 ns NPT window from three A1 starting conformations
(`nac_evt18_time1206ps`, `nac_evt25_time1462ps`, and
`nac_evt08_time1086ps`) with velocity seeds 26711, 26723, and 26737. A1 is
Thr267 O-gamma deprotonated and N-alpha-H3+. This audit is classical-MM
screening evidence only.

## Authoritative definitions

- NAC attack distance: PA66-L2 reactive carbonyl carbon to Thr267 OG1, at most
  0.35 nm.
- NAC attack angle: the same substrate carbonyl oxygen-carbon-Thr267 OG1 angle,
  95 to 115 degrees.
- A frame passes only when distance and angle pass jointly.
- Gate opening uses NylC residues 261-266; Thr267 is excluded.
- Residence continuity uses the actual trajectory sampling interval and never
  combines separate distance and angle occupancies.
- Only `npt300free` is scientific evidence. Restrained stages are excluded.

## Data flow

For each replica, generate periodicity-corrected primitive time series from
`run.tpr` and `run.xtc`: NAC distance, NAC angle, gate opening, ligand
retention/contact, and proton-network geometry. Extract temperature, pressure,
volume, and potential energy from `run.edr`. Scan `run.log` for
LINCS/SETTLE/NaN/FATAL.

Audit each replica independently, preserving `NOT_EVALUATED_*` for missing or
technically corrupt inputs. Merge all nine audit JSON files into one frozen
denominator table and rank only replicas that are numerically valid,
thermodynamically stable, bound, and contain NAC frames.

## Scientific outputs

Report per replica: frame denominator, NAC occupancy, longest continuous NAC,
distance/angle distributions, gate-opening distribution, potential-energy
summary, temperature/pressure/volume stability, ligand retention, and numerical
issues. Candidate-level summaries retain all three seeds and do not hide failed
replicas. A low-energy frame may be selected only from an actually NAC-positive
unrestrained frame; potential energy is a within-Hamiltonian screening metric,
not a QM/MM barrier.

## Failure handling

One failed replica does not stop the other eight. Technical completion is
separate from scientific status. No QM/MM Step1/Step2/PMF submission occurs
until this audit produces an eligible, stable NAC structure.
