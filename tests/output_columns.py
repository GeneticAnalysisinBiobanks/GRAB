#!/usr/bin/env python3
"""Run all four output-flag combinations and compare retained cells exactly.

Run from the repository root:
  python3 tests/output_columns.py --work-dir /path/to/scratch \
      [--reference-dir /path/to/pre-change-baseline] [--reference-bin /path/to/old/grab2]

Requires the baseline script's developer tools (Python, plink2, zstd).
The reference directory must contain a full-column pre-change baseline run.
"""

import argparse
import concurrent.futures
import gzip
import itertools
import math
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


MODES = {
    "full": ["--output-z", "--output-spa-status"],
    "default": [],
    "z": ["--output-z"],
    "status": ["--output-spa-status"],
}
TABLE = re.compile(r"\.(SPACox|SPAGRM|SAGELD|GALLOP|SPAGxE|SPAGxEmix|SPAmix|"
                   r"SPAmixPlus|SPAsqr|WtCoxG|LEAF|LocalP)(?:\.gz|\.zst)?$")


def read_lines(path):
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as stream:
            return stream.read().splitlines()
    if path.suffix == ".zst":
        return subprocess.check_output(["zstd", "-dc", str(path)], text=True).splitlines()
    return path.read_text().splitlines()


def tables(directory):
    return {p.relative_to(directory) for p in directory.rglob("*") if TABLE.search(p.name)}


def compare(full, selected, mode):
    paths = tables(full)
    if paths != tables(selected):
        raise AssertionError(f"{mode}: result-file sets differ")
    cells = 0
    for relative in sorted(paths):
        original = read_lines(full / relative)
        actual = read_lines(selected / relative)
        header = original[0].split("\t")
        indices = []
        for i, name in enumerate(header):
            status = re.search(r"(?:^|_)SPA_STATUS(?:_|$)", name) is not None
            z = (name == "Z" or name.startswith("Z_")) and not name.startswith("Z_Norm")
            if status and mode not in ("full", "status"):
                continue
            if z and mode not in ("full", "z"):
                continue
            indices.append(i)
        for line_no, (before, after) in enumerate(itertools.zip_longest(original, actual), 1):
            if before is None or after is None:
                raise AssertionError(f"{mode}/{relative}: row count differs")
            fields = before.split("\t")
            if len(fields) != len(header):
                raise AssertionError(f"{full / relative}:{line_no}: malformed reference row")
            expected = "\t".join(fields[i] for i in indices)
            if after != expected:
                raise AssertionError(f"{mode}/{relative}:{line_no}: selected columns or values differ")
            cells += len(indices)
    if not paths:
        raise AssertionError(f"No result tables in {full}")
    print(f"PASS {mode}: {len(paths)} tables, {cells} exact cells", flush=True)


def run_baseline(repo, binary, work, mode):
    destination = work / mode
    destination.mkdir()
    wrapper = work / f"grab2-{mode}"
    wrapper.write_text(
        f"#!{sys.executable}\nimport os, sys\n"
        "args = sys.argv[1:]\n"
        f"if '--method' in args: args += {MODES[mode]!r}\n"
        f"os.execv({str(binary)!r}, [{str(binary)!r}, *args])\n")
    wrapper.chmod(0o755)
    script = (repo / "examples/baseline.sh").read_text()
    script = script.replace("OUT_DIR=examples_output/baseline", "OUT_DIR=" + shlex.quote(str(destination)))
    script = script.replace("build/grab2", shlex.quote(str(wrapper)))
    script = script.replace("plink2 --pfile", "plink2 --threads 2 --pfile")
    script_path = work / f"baseline-{mode}.sh"
    script_path.write_text(script)
    with (work / f"baseline-{mode}.log").open("w") as log:
        subprocess.run(["bash", str(script_path)], cwd=repo, stdout=log, stderr=log, check=True)
    print(f"PASS baseline {mode}", flush=True)


def local_fixture(repo, work):
    fixture = work / "local-fixture"
    fixture.mkdir()
    # The bundled BIM/FAM fixtures contain literal backslash-t separators.
    # Decode copies in the test directory; preserve the repository fixtures.
    for suffix in ("fam", "bim"):
        source = repo / f"examples/hgdp_lanc.{suffix}"
        (fixture / f"input.{suffix}").write_text(source.read_text().replace("\\t", "\t"))
    shutil.copyfile(repo / "examples/hgdp_lanc.lanc", fixture / "input.lanc")
    with (fixture / "resid.txt").open("w") as stream:
        stream.write("FID\tIID\tR1\tR2\n")
        for i, line in enumerate((fixture / "input.fam").read_text().splitlines()):
            fid, iid = line.split()[:2]
            stream.write(f"{fid}\t{iid}\t{math.sin(i + 0.2):.17g}\t{math.cos(i + 0.3):.17g}\n")
    # Empty off-diagonal kinship; the method still computes diagonal variance.
    (fixture / "input.phi").write_text("idx1\tidx2\n")
    return ["--method", "SPAmixLocalPlus", "--lanc", str(fixture / "input"),
            "--admix-phi", str(fixture / "input.phi"), "--pheno", str(fixture / "resid.txt"),
            "--resid-name", "R1,R2", "--threads", "2"]


def run_local(binary, work, mode, common, flags):
    destination = work / f"local-{mode}"
    destination.mkdir()
    with (work / f"local-{mode}.log").open("w") as log:
        for label, extra in (("scan", ["--compression", "gz"]),
                             ("failed_qc", ["--mac", "1000000000", "--compression", "zst"])):
            subprocess.run([str(binary), *common, *extra, *flags,
                            "--out", str(destination / label)],
                           stdout=log, stderr=log, check=True)


def check_cli(binary):
    for flag in ("--output-z", "--output-spa-status"):
        help_result = subprocess.run([str(binary), "--help", flag[2:]], capture_output=True, text=True)
        assert help_result.returncode == 0 and flag in help_result.stderr
        duplicate = subprocess.run([str(binary), flag, flag], capture_output=True, text=True)
        assert duplicate.returncode != 0 and "specified more than once" in duplicate.stderr
    removed = subprocess.run([str(binary), "--output-z-norm"], capture_output=True, text=True)
    assert removed.returncode != 0 and "unknown option" in removed.stderr
    print("PASS CLI help, duplicate rejection, and unsupported flag rejection", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=Path("build/grab2"))
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path)
    parser.add_argument("--reference-bin", type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    binary, work = args.binary.resolve(), args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    check_cli(binary)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        jobs = [pool.submit(run_baseline, repo, binary, work, mode) for mode in MODES]
        for job in jobs:
            job.result()
    if args.reference_dir:
        compare(args.reference_dir.resolve(), work / "full", "full")
    for mode in MODES:
        compare(work / "full", work / mode, mode)
    common = local_fixture(repo, work)
    for mode, flags in MODES.items():
        run_local(binary, work, mode, common, flags)
    if args.reference_bin:
        run_local(args.reference_bin.resolve(), work, "before", common, [])
        compare(work / "local-before", work / "local-full", "full")
    for mode in MODES:
        compare(work / "local-full", work / f"local-{mode}", mode)


if __name__ == "__main__":
    main()
