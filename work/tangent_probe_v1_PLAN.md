# Tangent sphere accessibility pilot (2026-09-09)
User approved carbonyl-carbon tangent sphere endpoint plus external approach path.
- No density/support term in exposure score. Original parent coordinate/graph provenance retained.
- Sphere radii 0.5,1,1.5,2,3,4,5,6,8,10 A; 1024 Fibonacci directions in chemical C=O/leaving-atom frame. Refinement 2048 on 3 cases/material.
- Endpoint at carbon center+(1.7+r)*unit_direction; no deletion of target oxygen/N/O, neighbors, target C or explicit H. RAD inherited C1.7 N1.55 O1.52 H1.2 S1.8. Exact tangency allowed, tolerance1e-7 A.
- Continuous straight approach: swept-sphere segments of length<=2 A checked analytically against periodic neighbor atoms; free exterior is top/bottom half-space above/below ALL atom vdW extents plus probe radius. Ray <=120 A; unresolved reported separately and never counted as pass. Straight-ray restriction is not proof no curved path exists.
- First run six analytic tests RED then GREEN.
- Enumerate all model-chain C(=O)-O-C / C(=O)-N-C bonds, no support filter. Initial endpoint64-direction r0.5 screen used ONLY to select diverse pilot cases, not final exposure classification. Chain caps/loose-tail topology status not assessed by this pilot.
- Pilot up to12 sites/material spanning initial endpoint fractions0,(0,.05],(.05,.15],>.15; maintain all-pilot and original candidate denominators. Compute A(r), solid angles, connected opening components, sampled Rmax and censored flags. Score is log-radius average of A(r), empirical geometric summary, not energetic or enzyme activity.
- Verify arithmetic, path <=endpoint, 1024/2048 changes, replicate monotonic nesting for exact fixed directions, and empty/unresolved states. Write versioned outputs and append RUN_LOG/RUNBOOK.
- Do not restart abandoned material consensus fits or support50/30 tasks. No new donor scan or full enzyme fitting here.
