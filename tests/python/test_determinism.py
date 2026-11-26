"""Bit-for-bit replay."""

import rbsim


def run_hashes(n_steps, broadphase=rbsim.Broadphase.SWEEP_AND_PRUNE):
    w = rbsim.pile(100, broadphase)
    hashes = []
    for _ in range(n_steps):
        w.step(1.0 / 60.0)
        hashes.append(w.state_hash())
    return hashes


def test_python_run_reproduces_the_golden_hash():
    # Same constant the C++ test checks: stepping from Python must not
    # perturb a single bit.
    assert run_hashes(600)[-1] == rbsim.PILE100_HASH600


def test_repeat_runs_match_every_step():
    assert run_hashes(300) == run_hashes(300)
