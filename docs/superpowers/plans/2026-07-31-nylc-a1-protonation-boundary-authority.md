# NylC A1 protonation/boundary authority implementation

1. Add a failing contract test for one 4-task array:
   - tasks 0-1: neutral Thr267, current 146-QM mask.
   - tasks 2-3: activated Thr267, expanded mask including complete Thr268 (atoms 8964-8977).
   - all tasks: four inherited 250-cycle minimization blocks, no positional or reaction-coordinate restraints.
2. Implement one driver that reuses the validated MM-frame/PBC reconstruction, prepares the neutral HG1 coordinate or expanded QM mask, validates the expected QM contracts and audits numerical/heavy-skeleton/proton state after every block.
3. Add one runner and one Slurm array file (0-3, eight MPI ranks per task).
4. Deploy one immutable SCNet snapshot, verify targeted tests, Python syntax, Bash syntax, source/prmtop hashes and dry-run contracts, then submit exactly once.
5. Report technical and chemical-integrity denominators separately. Do not start steering, release, shooting, PMF, NEB or string calculations.
