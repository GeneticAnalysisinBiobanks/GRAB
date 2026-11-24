#!/usr/bin/env python3
"""Check SPAsqr defaults and overrides; pass --work-dir for disposable fixtures."""
import argparse
import csv
from pathlib import Path
import random
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="build/grab2")
    parser.add_argument("--work-dir", required=True)
    args = parser.parse_args()
    work = Path(args.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=False)
    rng = random.Random(71)
    n = 1000
    ids = [f"s{i}" for i in range(n)]
    y = [rng.gauss(3, 2) for _ in ids]
    (work / "pheno").write_text("IID L R\n" + "".join(
        f"{iid} {value} {value}\n" for iid, value in zip(ids, y)))
    (work / "ldak").write_text("FID IID Chr1 Chr2\n" + "".join(
        f"0 {iid} 0 0\n" for iid in ids))
    (work / "regenie").write_text("FID_IID " + " ".join(ids) + "\n" + "".join(
        f"{chrom} " + " ".join(["0"] * n) + "\n" for chrom in (1, 2)))
    (work / "pred").write_text(f"L {work / 'ldak'}\nR {work / 'regenie'}\n")
    with (work / "geno.vcf").open("w") as stream:
        stream.write('##fileformat=VCFv4.2\n##contig=<ID=1>\n##contig=<ID=2>\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                     '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
                     + "\t".join(ids) + "\n")
        for i, count in enumerate((10, 20, 30, 10)):
            g = ["0/1"] * count + ["0/0" if i < 3 else "1/1"] * (n - count)
            rng.shuffle(g)
            stream.write(f"{1 + i // 2}\t{i + 1}\tv{i}\tA\tC\t.\tPASS\t.\tGT\t"
                         + "\t".join(g) + "\n")
    common = [args.binary, "--method", "SPAsqr", "--vcf", str(work / "geno.vcf"),
              "--pheno", str(work / "pheno"), "--pheno-name", "L,R",
              "--spasqr-taus", "0.3,0.7", "--mac", "0", "--maf", "0",
              "--output-z", "--output-spa-status"]

    def run(name, extra=(), mode="wald"):
        result = subprocess.run(common + ["--out", str(work / name), "--spasqr-mode", mode]
                                + list(extra), text=True, capture_output=True)
        log = result.stdout + result.stderr
        (work / f"{name}.log").write_text(log)
        assert result.returncode == 0, log
        tables = []
        for pheno in ("L", "R"):
            with (work / f"{name}.{pheno}.SPAsqr").open() as stream:
                tables.append(list(csv.DictReader(stream, delimiter="\t")))
        return tables, log

    base, log = run("default")
    assert log.count("skipped heterogeneity test for 2 markers with MAF < 0.01") == 2
    assert "[WARN]" in log and log.count("pheno-transform = int") == 2
    enabled, log = run("enabled", ["--spasqr-het-min-maf", "0"])
    assert "skipped heterogeneity test" not in log
    het = {"BETA_HOM", "SE_HOM", "HET_Q", "LOG10P_HET"}
    for default, all_tests in zip(base, enabled):
        assert len(default) == 4
        for a, b in zip(default, all_tests):
            assert {k: v for k, v in a.items() if k not in het} == {
                k: v for k, v in b.items() if k not in het}
            assert all(b[k] != "NA" for k in het), b
            if a["ID"] in ("v0", "v3"):
                assert all(a[k] == "NA" for k in het), a
            else:  # Includes MAF exactly 0.01.
                assert all(a[k] == b[k] for k in het), a
    raised, log = run("raised", ["--spasqr-het-min-maf", "0.02"])
    assert log.count("skipped heterogeneity test for 4 markers") == 2
    assert all(row[k] == "NA" for table in raised for row in table for k in het)
    threaded, _ = run("threaded", ["--threads", "2", "--chunk-ksnp", "0.5"])
    assert threaded == base
    for value in ("-0.1", "0.51", "nan", "inf"):
        result = subprocess.run(common + ["--spasqr-het-min-maf", value], capture_output=True)
        assert result.returncode != 0, value

    pred = ["--pred-list", str(work / "pred")]
    for mode in ("score", "wald"):
        automatic, log = run(f"{mode}-auto", pred, mode)
        assert "[L]: pheno-transform = raw" in log
        assert "[R]: pheno-transform = standardize" in log
        assert "overriding" not in log
        if mode == "score":
            assert "bound=0.65" in log
            _, bound_log = run("bound", pred + ["--spasqr-outlier-abs-bound", "0.55"], mode)
            assert "bound=0.55" in bound_log
        else:
            assert log.count("skipped heterogeneity test for 2 markers") == 2
        for transform in ("raw", "int", "standardize"):
            explicit, log = run(f"{mode}-{transform}", pred + ["--pheno-transform", transform], mode)
            assert f"[L]: pheno-transform = {transform}" in log
            assert "[R]: pheno-transform = standardize" in log
            assert ("overriding" in log) == (transform != "standardize")
            assert explicit[1] == automatic[1]
            if transform == "raw":
                assert explicit == automatic
    print("PASS: SPAsqr thresholds, counts, transforms, overrides and thread reproducibility")


if __name__ == "__main__":
    main()
