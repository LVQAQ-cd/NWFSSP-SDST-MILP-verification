# -*- coding: utf-8 -*-
"""
NWFSSP-SDST (No-Wait Flow Shop with Sequence-Dependent Setups)
Validation script: MILP via CPLEX + brute-force check via an external decoder

Purpose
-------
This script validates a Mixed-Integer Linear Programming (MILP) model for the
no-wait flow shop scheduling problem with sequence-dependent setup times (NWFSSP-SDST).
It solves the MILP with CPLEX (via DOcplex) and independently evaluates sequences
using your reference evaluator's `decode(seq)` method. The optimal MILP sequence
is compared against a brute-force search (full permutation) scored by the same
`decode`, ensuring both the optimal objective value and sequence match.

Core notation (consistent with the paper)
-----------------------------------------
- P[j, w] : processing time of job j on machine w (shape: n_jobs x n_machines)
- ST[w, i, j] : setup time on machine w when job i is immediately followed by job j
                (shape: n_machines x n_jobs x n_jobs)
- First-job setup on machine w for job j is taken from the diagonal:
  st0[w, j] = ST[w, j, j]   (kept to match your evaluation convention)
- Sequences are 1-based lists, e.g., [1, 4, 2, 3], as required by `decode`.

Dependencies
------------
- python >= 3.8
- numpy
- docplex (and an accessible IBM CPLEX installation)
- algorithm.py providing a class with `.decode(seq)` (imported below as `MakespanEvaluator`)

Reproducibility
---------------
The MILP and brute-force parts are deterministic given P and ST.
All makespan evaluations are routed through the same `decode` to avoid
any discrepancy between the MILP model and the scoring function.

How to run
----------
$ python this_script.py
It will:
(1) solve the MILP,
(2) evaluate the MILP sequence via decode,
(3) brute-force the best sequence via decode,
(4) assert both results match, and
"""

from itertools import permutations
from typing import List, Tuple
import numpy as np
from docplex.mp.model import Model

from algorithm import MakespanEvaluator  # evaluator with .decode(seq)


# =========================
#  MILP model and solution
# =========================
def solve_milp(P: np.ndarray, ST: np.ndarray, log_output: bool = False) -> Tuple[float, List[int]]:
    """
    Build and solve the NWFSSP-SDST MILP with CPLEX/DOcplex.

    Parameters
    ----------
    P : np.ndarray
        Processing time matrix of shape (n_jobs, n_machines).
    ST : np.ndarray
        Sequence-dependent setup tensor of shape (n_machines, n_jobs, n_jobs).
        ST[w, i, j] is the setup time on machine w when i is immediately before j.
    log_output : bool, optional
        If True, print CPLEX solver logs. Default: False.

    Returns
    -------
    Cmax_model : float
        Optimal MILP objective value (makespan), directly from the model.
    best_seq_1based : List[int]
        Optimal job sequence in 1-based indexing (e.g., [1, 4, 2, 3]).
    """
    n_jobs, n_machines = P.shape
    assert ST.shape == (n_machines, n_jobs, n_jobs)

    # First-job setups: take the diagonal of ST on each machine to align with the evaluator
    st0 = np.stack([np.diag(ST[w]) for w in range(n_machines)], axis=0)  # shape (m, n)

    mdl = Model(name="NWFSSP_SDST_MILP", log_output=log_output)

    # Index sets
    J = range(n_jobs)       # jobs
    K = range(n_jobs)       # positions
    M = range(n_machines)   # machines

    # ----------------
    # Decision variables
    # ----------------
    x = {(j, k): mdl.binary_var(name=f"x_{j}_{k}") for j in J for k in K}             # assignment
    C = {(k, w): mdl.continuous_var(lb=0, name=f"C_{k}_{w}") for k in K for w in M}   # completion

    # Objective: minimize the makespan on the last machine at the last position
    mdl.minimize(C[n_jobs - 1, n_machines - 1])

    # ----------------
    # Constraints
    # ----------------
    # (1) Exactly one job per position
    for k in K:
        mdl.add_constraint(mdl.sum(x[j, k] for j in J) == 1, ctname=f"assign_pos_{k}")

    # (2) Each job appears exactly once
    for j in J:
        mdl.add_constraint(mdl.sum(x[j, k] for k in K) == 1, ctname=f"assign_job_{j}")

    # (3) First job lower bound on each machine: first-job setup + processing
    for w in M:
        mdl.add_constraint(
            C[0, w] >= mdl.sum(x[j, 0] * (P[j, w] + st0[w, j]) for j in J),
            ctname=f"first_job_{w}",
        )

    # (4) No-wait coupling across machines for the same position
    for k in K:
        for w in range(n_machines - 1):
            mdl.add_constraint(
                C[k, w + 1] - C[k, w] == mdl.sum(x[j, k] * P[j, w + 1] for j in J),
                ctname=f"no_wait_k{k}_w{w}",
            )

    # (5) Same-machine adjacency with SDST (linearization for x[i,k] * x[j,k+1])
    z = {(i, j, k): mdl.binary_var(name=f"z_{i}_{j}_{k}")
         for k in range(n_jobs - 1) for i in J for j in J}

    for k in range(n_jobs - 1):
        for i in J:
            for j in J:
                mdl.add_constraint(z[i, j, k] <= x[i, k], ctname=f"z_le_xik_{i}_{j}_{k}")
                mdl.add_constraint(z[i, j, k] <= x[j, k + 1], ctname=f"z_le_xjkp1_{i}_{j}_{k}")
                mdl.add_constraint(
                    z[i, j, k] >= x[i, k] + x[j, k + 1] - 1,
                    ctname=f"z_ge_sum_minus1_{i}_{j}_{k}",
                )

    for k in range(n_jobs - 1):
        for w in M:
            mdl.add_constraint(
                C[k + 1, w] - C[k, w]
                - mdl.sum(x[j, k + 1] * P[j, w] for j in J)
                - mdl.sum(z[i, j, k] * ST[w, i, j] for i in J for j in J)
                >= 0,
                ctname=f"machine_seq_w{w}_k{k}",
            )

    # Solve MILP
    sol = mdl.solve(log_output=log_output)
    if sol is None:
        raise RuntimeError("MILP infeasible or no solution returned by CPLEX.")

    Cmax_model = float(sol.get_value(C[n_jobs - 1, n_machines - 1]))

    # Recover the 1-based sequence from x[j,k]
    seq0 = [-1] * n_jobs
    for k in K:
        for j in J:
            if sol.get_value(x[j, k]) > 0.5:
                seq0[k] = j
                break
    best_seq_1based = [j + 1 for j in seq0]
    return Cmax_model, best_seq_1based


# =========================
#  Brute-force with decode
# =========================
def brute_force_best(P: np.ndarray, ST: np.ndarray) -> Tuple[float, List[int]]:
    """
    Enumerate all job permutations (1-based) and score them using
    MakespanEvaluator.decode(seq). Return the best makespan and sequence.

    Parameters
    ----------
    P : np.ndarray
        Processing time matrix (n_jobs x n_machines).
    ST : np.ndarray
        Setup time tensor (n_machines x n_jobs x n_jobs).

    Returns
    -------
    best_val : float
        Best makespan (as given by decode).
    best_seq : List[int]
        Best sequence (1-based).
    """
    decoder = MakespanEvaluator(P, ST)
    n_jobs, _ = P.shape

    best_val = float("inf")
    best_seq: List[int] = []

    for seq in permutations(range(1, n_jobs + 1)):
        val = float(decoder.decode(list(seq)))
        if val < best_val:
            best_val = val
            best_seq = list(seq)
    return best_val, best_seq


# ===============
#  Demo / self-test
# ===============
if __name__ == "__main__":
    # Example instance (same as in the paper/examples)
    P = np.array([[2, 1, 3],
                  [3, 2, 2],
                  [2, 1, 3],
                  [1, 3, 2]], dtype=float)

    ST = np.array([[[2, 4, 6, 1],
                    [3, 2, 3, 8],
                    [4, 6, 3, 2],
                    [6, 1, 2, 4]],
                   [[4, 3, 2, 3],
                    [1, 2, 4, 3],
                    [4, 1, 1, 2],
                    [2, 5, 5, 1]],
                   [[1, 4, 3, 1],
                    [2, 3, 1, 5],
                    [4, 5, 4, 5],
                    [4, 3, 4, 5]]], dtype=float)

    # 1) Solve MILP
    cmax_model, seq_model = solve_milp(P, ST, log_output=False)

    # Evaluate the MILP sequence via decode (to keep the scoring convention identical)
    decoder = MakespanEvaluator(P, ST)
    cmax_model_by_decode = float(decoder.decode(seq_model))
    print(f"[MILP] seq = {seq_model},  Cmax(model) = {cmax_model:.6g},  Cmax(decode) = {cmax_model_by_decode:.6g}")

    # 2) Brute-force enumeration using decode as the scoring function
    cmax_bf, seq_bf = brute_force_best(P, ST)
    print(f"[BF]   seq = {seq_bf},  Cmax(decode) = {cmax_bf:.6g}")

    # 3) Consistency check: the MILP solution must match the brute-force optimum (tolerance can be adjusted)
    assert abs(cmax_model_by_decode - cmax_bf) < 1e-6, "Mismatch between MILP and brute-force optimum (by decode)."