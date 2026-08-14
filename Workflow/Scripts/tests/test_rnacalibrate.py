from __future__ import annotations
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

import pytest
import rnacalibrate as rc

SCRIPTS = str(Path(__file__).resolve().parent.parent)


def _fasta(tmp_path, text):
    p = tmp_path / "t.fna"
    p.write_text(text)
    return p


def test_length_stats_population_stdev(tmp_path):
    f = _fasta(tmp_path, ">a\nAAAA\n>b\nAAAAAA\n")
    stats = rc.compute_target_length_stats(f)
    assert stats == {"count": 2, "mean": 5.0, "std": 1.0}


def test_length_stats_single_record_has_zero_stdev(tmp_path):
    f = _fasta(tmp_path, ">a\nAAAA\n")
    assert rc.compute_target_length_stats(f)["std"] == 0.0


def test_length_stats_rejects_empty_file(tmp_path):
    with pytest.raises(RuntimeError, match="No FASTA records"):
        rc.compute_target_length_stats(_fasta(tmp_path, ""))


def test_build_length_arg_rounds_to_ints():
    assert rc.build_length_arg({"count": 3, "mean": 986.4, "std": 792.2}) == {
        "mean": 986, "std": 792, "value": "986,792",
    }


LENGTH_ARG = {"mean": 986, "std": 792, "value": "986,792"}


def test_build_command_minimal_omits_optional_flags():
    assert rc.build_command("RNAcalibrate", "q.fa", "t.fna", 10000, 50000, LENGTH_ARG) == [
        "RNAcalibrate", "-k", "10000", "-q", "q.fa", "-t", "t.fna",
        "-m", "50000", "-l", "986,792",
    ]


def test_build_command_appends_every_optional_flag_in_order():
    assert rc.build_command(
        "RNAcalibrate", "q.fa", "t.fna", 10000, 50000, LENGTH_ARG,
        randomize_targets=True, max_internal_loop=4, max_bulge_loop=5, seed="2,7",
    ) == [
        "RNAcalibrate", "-k", "10000", "-q", "q.fa", "-t", "t.fna",
        "-m", "50000", "-l", "986,792", "-u", "4", "-v", "5", "-f", "2,7", "-s",
    ]


def test_build_command_seed_is_forced_helix_not_rng():
    # Guards the naming trap: `seed` here means RNAhybrid's -f, never an RNG seed.
    cmd = rc.build_command("RNAcalibrate", "q.fa", "t.fna", 1, 1, LENGTH_ARG, seed="2,7")
    assert "-f" in cmd and "2,7" in cmd


def test_parse_output_returns_one_entry_per_row():
    parsed = rc.parse_rnacalibrate_output(
        "hsa-let-7a-5p 64 2.118843 0.166056\nhsa-miR-21-5p 63 2.025677 0.163574\n"
    )
    assert parsed == {"per_query": [
        {"query": "hsa-let-7a-5p", "sample_size": 64, "xi": 2.118843, "theta": 0.166056},
        {"query": "hsa-miR-21-5p", "sample_size": 63, "xi": 2.025677, "theta": 0.163574},
    ]}


def test_parse_output_rejects_nan():
    with pytest.raises(RuntimeError, match="NaN"):
        rc.parse_rnacalibrate_output("hsa-miR-2861 105 -nan -nan\n")


def test_parse_output_rejects_wrong_column_count():
    with pytest.raises(RuntimeError, match="expected 4 columns"):
        rc.parse_rnacalibrate_output("hsa-let-7a-5p 64 2.118843\n")


def test_parse_output_rejects_empty():
    with pytest.raises(RuntimeError, match="did not produce any calibration rows"):
        rc.parse_rnacalibrate_output("\n\n")


def test_derive_seed_is_stable_against_known_values():
    # Hardcoded on purpose: these must not drift between runs, machines, or releases.
    assert rc.derive_seed(1, "UGAGGUAGUAGGUUGUAUAGUU") == 2093734617
    assert rc.derive_seed(1, "GGGGCCUGGCGGUGGGCGG") == 1920725615
    assert rc.derive_seed(1, "ACGU") == 600126125


def test_derive_seed_differs_per_sequence_and_per_base():
    assert rc.derive_seed(1, "ACGU") != rc.derive_seed(1, "ACGA")
    assert rc.derive_seed(1, "ACGU") != rc.derive_seed(2, "ACGU")


def test_derive_seed_stays_in_faketime_safe_range():
    for seq in ("ACGU", "UGAGGUAGUAGGUUGUAUAGUU", "G" * 200):
        value = rc.derive_seed(2**31 - 1, seq)
        assert 0 <= value < 2**31


def test_derive_seed_ignores_pythonhashseed():
    # The regression guard for using built-in hash(), which is salted per process.
    code = (
        f"import sys; sys.path.insert(0, {SCRIPTS!r}); "
        "import rnacalibrate; print(rnacalibrate.derive_seed(1, 'ACGU'))"
    )
    outputs = set()
    for hashseed in ("0", "1", "12345"):
        env = {**os.environ, "PYTHONHASHSEED": hashseed}
        done = subprocess.run([sys.executable, "-c", code], capture_output=True,
                              text=True, env=env, check=True)
        outputs.add(done.stdout.strip())
    assert outputs == {"600126125"}


def test_faketime_env_uses_ld_preload_on_linux():
    env = rc.build_faketime_env(42, "/opt/lib/libfaketime.so.1", "linux", {"PATH": "/bin"})
    assert env["LD_PRELOAD"] == "/opt/lib/libfaketime.so.1"
    assert env["FAKETIME_FMT"] == "%s"
    assert env["FAKETIME"] == "@42"
    assert env["PATH"] == "/bin"           # base env is preserved
    assert "DYLD_INSERT_LIBRARIES" not in env


def test_faketime_env_uses_dyld_on_darwin():
    env = rc.build_faketime_env(42, "/opt/lib/libfaketime.1.dylib", "darwin", {})
    assert env["DYLD_INSERT_LIBRARIES"] == "/opt/lib/libfaketime.1.dylib"
    assert env["DYLD_FORCE_FLAT_NAMESPACE"] == "1"
    assert env["FAKETIME"] == "@42"
    assert "LD_PRELOAD" not in env


def test_faketime_env_does_not_mutate_the_base_env():
    base = {"PATH": "/bin"}
    rc.build_faketime_env(1, "/x.so", "linux", base)
    assert base == {"PATH": "/bin"}


def test_faketime_env_rejects_unknown_platform():
    with pytest.raises(RuntimeError, match="Unsupported platform"):
        rc.build_faketime_env(1, "/x.so", "win32", {})


def test_faketime_env_seed_is_an_absolute_frozen_timestamp():
    # The '@' prefix freezes the clock; without it libfaketime treats the value as an offset.
    assert rc.build_faketime_env(0, "/x.so", "linux", {})["FAKETIME"] == "@0"


def test_resolve_faketime_library_reports_the_env_file_when_missing(monkeypatch):
    monkeypatch.setattr(rc, "which_required",
                        lambda *a: (_ for _ in ()).throw(RuntimeError("not found")))
    with pytest.raises(RuntimeError, match="rnahybrid.yaml"):
        rc.resolve_faketime_library("linux")


def test_resolve_faketime_library_finds_the_object_beside_the_wrapper(tmp_path, monkeypatch):
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "faketime").write_text("")
    (tmp_path / "lib" / "faketime").mkdir(parents=True)
    (tmp_path / "lib" / "faketime" / "libfaketime.so.1").write_text("")
    monkeypatch.setattr(rc, "which_required", lambda *a: str(tmp_path / "bin" / "faketime"))
    found = rc.resolve_faketime_library("linux")
    # .resolve() on both sides: macOS symlinks /tmp -> /private/tmp.
    assert found == (tmp_path / "lib" / "faketime" / "libfaketime.so.1").resolve()


def test_resolve_faketime_library_errors_when_wrapper_exists_but_object_does_not(tmp_path, monkeypatch):
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "faketime").write_text("")
    monkeypatch.setattr(rc, "which_required", lambda *a: str(tmp_path / "bin" / "faketime"))
    with pytest.raises(RuntimeError, match="libfaketime.so.1"):
        rc.resolve_faketime_library("linux")


@pytest.mark.parametrize("first, repeat, other, verdict", [
    ("A", "A", "B", "verified"),
    ("A", "B", "C", "blocked"),
    # Two identical runs alone would look like success; the third run is what proves
    # the seed actually controls the output.
    ("A", "A", "A", "blind"),
    ("A", "B", "A", "blocked"),     # blocked outranks blind
])
def test_classify_probe(first, repeat, other, verdict):
    assert rc.classify_probe(first, repeat, other) == verdict


@pytest.mark.parametrize("rng_seed, randomize_targets, expected", [
    (1, True, True),
    (None, True, False),
    # `-s` is the only consumer of the RNG. Without it RNAcalibrate is already
    # byte-reproducible, so pinning is inert — and probing would report a false "blind"
    # and fail a valid config.
    (1, False, False),
    (None, False, False),
])
def test_pinning_enabled(rng_seed, randomize_targets, expected):
    assert rc.pinning_enabled(rng_seed, randomize_targets) is expected


def test_query_key_takes_the_first_whitespace_token():
    assert rc.query_key("hsa-miR-21-5p") == "hsa-miR-21-5p"
    assert rc.query_key("hsa-miR-21-5p some description") == "hsa-miR-21-5p"


def test_iter_query_records_yields_header_and_sequence(tmp_path):
    q = tmp_path / "q.fa"
    q.write_text(">hsa-let-7a-5p\nUGAGGUAG\nUAGGUU\n>hsa-miR-21-5p\nUAGCUUAU\n")
    assert list(rc.iter_query_records(q)) == [
        ("hsa-let-7a-5p", "UGAGGUAGUAGGUU"),
        ("hsa-miR-21-5p", "UAGCUUAU"),
    ]


def test_iter_query_records_rejects_duplicate_ids(tmp_path):
    # Today RNAcalibrate emits two rows and load_per_query_distributions silently
    # keeps the last. split_query_per_miRNA already guards this; match it.
    q = tmp_path / "q.fa"
    q.write_text(">hsa-let-7a-5p\nACGU\n>hsa-let-7a-5p\nUGCA\n")
    with pytest.raises(RuntimeError, match="Duplicate query ID"):
        list(rc.iter_query_records(q))


def test_iter_query_records_rejects_an_empty_file(tmp_path):
    q = tmp_path / "q.fa"
    q.write_text("")
    with pytest.raises(RuntimeError, match="No FASTA records"):
        list(rc.iter_query_records(q))


def test_expect_single_row_returns_the_entry():
    parsed = {"per_query": [{"query": "m", "sample_size": 1, "xi": 0.1, "theta": 0.2}]}
    assert rc.expect_single_row(parsed, "m")["xi"] == 0.1


def test_expect_single_row_rejects_multiple_rows():
    parsed = {"per_query": [{"query": "a"}, {"query": "b"}]}
    with pytest.raises(RuntimeError, match="exactly 1 calibration row"):
        rc.expect_single_row(parsed, "a")


def test_sha256_file_matches_a_known_digest(tmp_path):
    f = tmp_path / "x.txt"
    f.write_bytes(b"ACGU\n")
    assert rc.sha256_file(f) == (
        "d57db2001ce2b4da9b363b703c19d7f7a94772c0c895388d06812de225366798"
    )


@pytest.mark.parametrize("variant", ["acgu", "ACGT", "acgt", " ACGU \n", "\tACGU"])
def test_derive_seed_and_sequence_sha256_share_one_normalisation(variant):
    # The contract is that both consumers hash the same string. Asserted without restating
    # derive_seed's mixing arithmetic, which would only ever agree with itself.
    assert rc.derive_seed(1, variant) == rc.derive_seed(1, "ACGU")
    assert rc.sequence_sha256(variant) == rc.sequence_sha256("ACGU")


def test_count_fasta_records(tmp_path):
    f = tmp_path / "x.fa"
    f.write_text(">a\nACGU\n>b\nUGCA\n>c\nAAAA\n")
    assert rc.count_fasta_records(f) == 3


@pytest.mark.parametrize("rng_seed", [1, None])
def test_build_inputs_block_shape(tmp_path, rng_seed):
    q = tmp_path / "q.fa"
    q.write_text(">m\nACGU\n")
    t = tmp_path / "t.fna"
    t.write_text(">g\nAAAA\n>h\nCCCC\n")
    block = rc.build_inputs_block(
        query=q, target=t, k=10000, max_target_length=50000,
        length_arg={"mean": 986, "std": 792, "value": "986,792"},
        forced_helix=None, max_internal_loop=None, max_bulge_loop=None, rng_seed=rng_seed,
    )
    assert block["query"]["n_records"] == 1
    assert block["target"]["n_records"] == 2
    assert block["query"]["sha256"] == rc.sha256_file(q)
    assert block["target"]["sha256"] == rc.sha256_file(t)
    assert block["params"] == {
        "k": 10000, "max_target_length": 50000, "length_arg": "986,792",
        "forced_helix": None, "max_internal_loop": None, "max_bulge_loop": None,
    }
    assert block["rng_seed"] == rng_seed


# ---- verify_faketime: probe sequencing and the sleep between the two same-seed probes ----
# subprocess.run and time.sleep are faked, so these run instantly despite the real 1.1s sleep.


def _completed(stdout):
    return SimpleNamespace(stdout=stdout)


def test_verify_faketime_sleeps_between_first_and_second_probe_then_verifies(monkeypatch):
    calls = []

    def fake_run(command, **kwargs):
        calls.append(("run", kwargs["env"]["FAKETIME"]))
        return _completed(kwargs["env"]["FAKETIME"])

    monkeypatch.setattr(rc.subprocess, "run", fake_run)
    monkeypatch.setattr(rc.time, "sleep", lambda seconds: calls.append(("sleep", seconds)))

    rc.verify_faketime("RNAcalibrate", "q.fa", "t.fna", 50000, LENGTH_ARG, 1,
                       "/opt/lib/libfaketime.so.1", "linux")

    # Same seed twice (A, B) around the sleep, then a different seed (C) with no sleep
    # before it — matching classify_probe's A==B, A!=C contract.
    assert [c[0] for c in calls] == ["run", "sleep", "run", "run"]
    assert calls[1][1] >= 1.0


def test_verify_faketime_retries_at_larger_k_when_blind_at_k5(monkeypatch):
    # k=5 can't discriminate (every probe returns the same output regardless of seed);
    # k=200 (PROBE_K_RETRY) does. Worst case is 3 probes + 1 sleep per k level.
    calls = []

    def fake_run(command, **kwargs):
        probe_k = command[2]
        calls.append(("run", probe_k))
        if probe_k == str(rc.PROBE_K):
            return _completed("indistinguishable")
        return _completed(kwargs["env"]["FAKETIME"])

    monkeypatch.setattr(rc.subprocess, "run", fake_run)
    monkeypatch.setattr(rc.time, "sleep", lambda seconds: calls.append(("sleep", seconds)))

    rc.verify_faketime("RNAcalibrate", "q.fa", "t.fna", 50000, LENGTH_ARG, 1,
                       "/opt/lib/libfaketime.so.1", "linux")

    assert [c[0] for c in calls] == [
        "run", "sleep", "run", "run", "run", "sleep", "run", "run",
    ]


# ---- error messages carry the full command, not just command[:2] ----

def test_verify_faketime_probe_error_includes_the_full_command(monkeypatch):
    def fake_run(command, **kwargs):
        raise subprocess.CalledProcessError(2, command, output="", stderr="blocked")

    monkeypatch.setattr(rc.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError) as excinfo:
        rc.verify_faketime("RNAcalibrate", "q.fa", "t.fna", 50000, LENGTH_ARG, 1,
                           "/opt/lib/libfaketime.so.1", "linux")

    message = str(excinfo.value)
    # The old bug rendered only command[:2] ("RNAcalibrate -k"); -t/-m only appear once
    # the full command is joined.
    assert "-t t.fna" in message
    assert "-m 50000" in message


def test_run_rnacalibrate_error_includes_the_full_command(tmp_path, monkeypatch):
    q = tmp_path / "q.fa"; q.write_text(">m\nACGU\n")
    t = tmp_path / "t.fna"; t.write_text(">g\nAAAA\n")

    monkeypatch.setattr(rc, "which_required", lambda *a: "RNAcalibrate")

    def fake_run(command, **kwargs):
        raise subprocess.CalledProcessError(1, command, output="", stderr="boom")

    monkeypatch.setattr(rc.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError) as excinfo:
        rc.run_rnacalibrate(q, t, tmp_path / "out.json", k=10, max_target_length=100)

    message = str(excinfo.value)
    assert "RNAcalibrate failed for miRNA m" in message
    assert "-t" in message and "-m 100" in message


@pytest.mark.parametrize("randomize_targets, rng_seed, warns", [
    (False, 42, True),          # set but inert: `-s` is never passed, so nothing to pin
    (False, None, False),
    (True, 42, False),          # the normal pinned case
])
def test_run_rnacalibrate_warns_only_when_the_seed_is_inert(
        tmp_path, monkeypatch, capsys, randomize_targets, rng_seed, warns):
    q = tmp_path / "q.fa"; q.write_text(">m\nACGU\n")
    t = tmp_path / "t.fna"; t.write_text(">g\nAAAA\n")

    monkeypatch.setattr(rc, "which_required", lambda *a: "RNAcalibrate")
    monkeypatch.setattr(rc, "resolve_faketime_library", lambda *a: "/opt/lib/libfaketime.so.1")
    monkeypatch.setattr(rc, "verify_faketime", lambda *a, **k: None)  # probed above
    monkeypatch.setattr(rc.subprocess, "run", lambda *a, **k: _completed("m 1 1.0 1.0\n"))

    rc.run_rnacalibrate(q, t, tmp_path / "out.json", k=10, max_target_length=100,
                        randomize_targets=randomize_targets, rng_seed=rng_seed)

    err = capsys.readouterr().err
    if warns:
        assert f"rng_seed={rng_seed}" in err and "randomize_targets" in err
    else:
        assert err == ""
