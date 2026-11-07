#!/usr/bin/env python3
"""Download and validate an MGnify genome catalogue in the layout required by
MetaPilot ("Add MAG database").

    <OUT>/<catalogue>/<version>/
        genomes-all_metadata.tsv
        original_db/<genomeID>.faa
        eggNOG/<genomeID>_eggNOG.tsv

Cross-platform (Windows, macOS, Linux). Python >= 3.8, standard library only.
Re-runnable: files already complete are skipped.

Author: Primo Veneroso. License: MIT.
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes"
ID_RE = re.compile(r"^MGYG[0-9]+(\.[0-9]+)?$")
FASTA_SEQ_RE = re.compile(r"^[A-Za-z*. -]+$")
USER_AGENT = "mgnify-metapilot-db/1.0 (Primo Veneroso)"


# ---------------------------------------------------------------- download
def fetch(url: str, dest: Path, attempts: int = 6, show_bar: bool = False) -> bool:
    """Download url to dest via dest.part; never leaves a partial final file."""
    part = dest.with_name(dest.name + ".part")
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as resp, open(part, "wb") as out:
                total = int(resp.headers.get("Content-Length") or 0)
                done = 0
                while True:
                    chunk = resp.read(1 << 16)
                    if not chunk:
                        break
                    out.write(chunk)
                    done += len(chunk)
                    if show_bar and total and sys.stdout.isatty():
                        pct = done * 100 // total
                        print(f"\r  {pct:3d}% ({done >> 20}/{total >> 20} MB)", end="", flush=True)
            if show_bar and sys.stdout.isatty():
                print()
            if part.stat().st_size > 0:
                os.replace(part, dest)
                return True
        except urllib.error.HTTPError as e:
            if e.code == 404:
                part.unlink(missing_ok=True)
                return False
        except (urllib.error.URLError, OSError, TimeoutError):
            pass
        part.unlink(missing_ok=True)
        time.sleep(attempt * 5)
    return False


def unpack_gzip(gz: Path, final: Path) -> bool:
    """Verify and decompress gz into final; on failure leave no final file."""
    tmp = final.with_name(final.name + ".part")
    try:
        if not gz.is_file() or gz.stat().st_size == 0:
            return False
        with gzip.open(gz, "rb") as src, open(tmp, "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 20)
        if tmp.stat().st_size == 0:
            tmp.unlink(missing_ok=True)
            return False
        os.replace(tmp, final)
        gz.unlink(missing_ok=True)
        return True
    except (OSError, EOFError, gzip.BadGzipFile):
        tmp.unlink(missing_ok=True)
        return False


def fetch_gz(url: str, final: Path) -> bool:
    if final.is_file() and final.stat().st_size > 0:
        return True
    gz = final.with_name(final.name + ".gz")
    for p in (gz, final.with_name(final.name + ".part")):
        p.unlink(missing_ok=True)
    if not fetch(url, gz):
        return False
    if unpack_gzip(gz, final):
        return True
    gz.unlink(missing_ok=True)
    return False


def download_genome(base_url: str, gid: str, orig: Path, egg: Path) -> list[str]:
    url = f"{base_url}/species_catalogue/{gid[:11]}/{gid}/genome"
    failed = []
    if not fetch_gz(f"{url}/{gid}.faa.gz", orig / f"{gid}.faa"):
        failed.append(f"{gid} faa")
    if not fetch_gz(f"{url}/{gid}_eggNOG.tsv.gz", egg / f"{gid}_eggNOG.tsv"):
        failed.append(f"{gid} eggNOG")
    return failed


def progress(done: int, total: int, start: float, errors: int, final: bool = False) -> None:
    width = 26
    filled = done * width // total if total else width
    elapsed = max(time.time() - start, 1e-6)
    rate = done / elapsed
    eta = int((total - done) / rate) if rate > 0 and done < total else 0
    line = (f"[{'#' * filled}{'-' * (width - filled)}] {done * 100 / total:5.1f}% "
            f"({done}/{total}) | {rate:.1f} gen/s | ETA {eta // 60:02d}:{eta % 60:02d} | err:{errors}")
    if sys.stdout.isatty():
        print("\r\033[K" + line, end="\n" if final else "", flush=True)
    elif final:
        print(line)


# ---------------------------------------------------------------- metadata
def read_species(meta: Path) -> list[str]:
    with open(meta, encoding="utf-8", newline="") as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for needed in ("Genome", "Lineage", "Species_rep"):
            if needed not in col:
                sys.exit(f"ERROR: column '{needed}' missing in {meta}")
        ids = set()
        for line in fh:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) > col["Species_rep"]:
                sid = fields[col["Species_rep"]].strip()
                if ID_RE.match(sid):
                    ids.add(sid)
    return sorted(ids)


# ---------------------------------------------------------------- validation
def check_fasta(path: Path) -> list[str]:
    gid = path.name[:-4]
    problems, seen_header, has_seq = [], False, False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                seen_header = True
                h = line[1:].strip()
                if not h.startswith(gid + "_"):
                    problems.append(f"{path.name}: header genomeID mismatch: {h}")
                    break
                rest = h[len(gid) + 1:]
                sp = rest.find(" ")
                if sp <= 0 or not rest[sp + 1:].strip():
                    problems.append(f"{path.name}: missing proteinID or protein name: {h}")
                    break
            elif line.strip():
                if not seen_header:
                    problems.append(f"{path.name}:{n}: sequence before first header")
                    break
                if not FASTA_SEQ_RE.match(line):
                    problems.append(f"{path.name}:{n}: unexpected characters in sequence")
                    break
                has_seq = True
    if not problems and not has_seq:
        problems.append(f"{path.name}: no sequences")
    return problems


def check_eggnog(path: Path) -> list[str]:
    gid = path.name[: -len("_eggNOG.tsv")]
    rows = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        first = fh.readline().rstrip("\r\n")
        hdr_nf = len(first.split("\t"))
        if hdr_nf < 2:
            return [f"{path.name}: header not tab-separated"]
        for n, line in enumerate(fh, 2):
            line = line.rstrip("\r\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            rows += 1
            fields = line.split("\t")
            if len(fields) != hdr_nf:
                return [f"{path.name}:{n}: {len(fields)} columns, header has {hdr_nf}"]
            if not fields[0].startswith(gid + "_"):
                return [f"{path.name}:{n}: query ID does not start with {gid}_"]
    return [] if rows else [f"{path.name}: no data rows"]


def validate(db: Path, meta: Path, expected: list[str], limit: int, jobs: int,
             failed: list[str]) -> int:
    orig, egg = db / "original_db", db / "eggNOG"
    errors = 0

    def err(msg: str) -> None:
        nonlocal errors
        errors += 1
        print("✗ " + msg)

    # Rule 6: Genome/Lineage filled on every row
    with open(meta, encoding="utf-8", newline="") as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        gi, li = header.index("Genome"), header.index("Lineage")
        empty = sum(1 for line in fh
                    if (f := line.rstrip("\r\n").split("\t")) and
                    (len(f) <= max(gi, li) or not f[gi].strip() or not f[li].strip()))
    if empty:
        err(f"Rule 6: {empty} metadata rows with empty Genome or Lineage")

    faa = sorted(p.name[:-4] for p in orig.glob("*.faa") if p.is_file())
    eggs = sorted(p.name[: -len("_eggNOG.tsv")] for p in egg.glob("*_eggNOG.tsv") if p.is_file())
    print(f"original_db: {len(faa)} files   eggNOG: {len(eggs)} files   expected: {len(expected)}")
    if len(faa) != len(eggs):
        err("Rule 1: file counts differ")

    unexpected = [p for p in orig.iterdir() if p.is_file() and not p.name.endswith(".faa")]
    unexpected += [p for p in egg.iterdir() if p.is_file() and not p.name.endswith("_eggNOG.tsv")]
    if unexpected:
        err(f"Rule 2: {len(unexpected)} unexpected or leftover (.part/.gz) files, e.g. {unexpected[0]}")
    empty_files = [p for d in (orig, egg) for p in d.iterdir() if p.is_file() and p.stat().st_size == 0]
    if empty_files:
        err(f"{len(empty_files)} empty (0 byte) files, e.g. {empty_files[0]}")

    sf, se, sx = set(faa), set(eggs), set(expected)
    if sf ^ se:
        err(f"Rule 3: {len(sf ^ se)} genomes without a pair")
    if limit == 0:
        for name, have in (("FASTA", sf), ("eggNOG", se)):
            missing, extra = sorted(sx - have), sorted(have - sx)
            if missing or extra:
                err(f"{name} IDs differ from catalogue: {len(missing)} missing "
                    f"{missing[:5]}, {len(extra)} unexpected {extra[:5]}")
    else:
        miss = len(sx - sf) + len(sx - se)
        if miss:
            err(f"{miss} files missing for the genomes selected with -n")

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        bad_faa = [m for r in ex.map(check_fasta, (orig / f"{g}.faa" for g in faa)) for m in r]
        bad_egg = [m for r in ex.map(check_eggnog, (egg / f"{g}_eggNOG.tsv" for g in eggs)) for m in r]
    if bad_faa:
        err(f"Rule 4: {len(bad_faa)} FASTA problems, e.g. {bad_faa[0]}")
    if bad_egg:
        err(f"Rule 5: {len(bad_egg)} eggNOG problems, e.g. {bad_egg[0]}")
    if failed:
        err(f"{len(failed)} downloads failed (see failed.txt). Re-run to retry.")
    return errors


# ---------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", default="MGnify_MAG_DB", help="output root (default: %(default)s)")
    ap.add_argument("-c", "--catalogue", default="human-gut", help="MGnify catalogue (default: %(default)s)")
    ap.add_argument("-v", "--version", default="v3.0", help="catalogue version (default: %(default)s)")
    ap.add_argument("-j", "--jobs", type=int, default=8, help="parallel downloads (default: %(default)s)")
    ap.add_argument("-n", "--limit", type=int, default=0, help="only first N species (quick test)")
    ap.add_argument("-V", "--validate-only", action="store_true", help="validate local folder, no download")
    args = ap.parse_args()
    if args.jobs < 1:
        ap.error("--jobs must be >= 1")

    base_url = f"{BASE}/{args.catalogue}/{args.version}"
    db = Path(args.out) / args.catalogue / args.version
    orig, egg, work = db / "original_db", db / "eggNOG", db / ".download_work"
    for d in (orig, egg, work):
        d.mkdir(parents=True, exist_ok=True)
    meta = db / "genomes-all_metadata.tsv"

    print(f"Catalogue: {base_url}\nFolder   : {db}\nParallel : {args.jobs}")

    if not (meta.is_file() and meta.stat().st_size > 0):
        if args.validate_only:
            sys.exit(f"ERROR: metadata missing: {meta}")
        print("Downloading genomes-all_metadata.tsv (~110 MB) ...")
        if not fetch(f"{base_url}/genomes-all_metadata.tsv", meta, show_bar=True):
            sys.exit("ERROR: metadata download failed")

    species = read_species(meta)
    if args.limit > 0:
        species = species[: args.limit]
    if not species:
        sys.exit("ERROR: no species found in metadata")
    print(f"Representative species: {len(species)}")

    failed: list[str] = []
    if not args.validate_only:
        print("Downloading FASTA and eggNOG ...")
        start, done, lock = time.time(), 0, threading.Lock()
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futures = [ex.submit(download_genome, base_url, g, orig, egg) for g in species]
            for fut in as_completed(futures):
                with lock:
                    failed.extend(fut.result())
                    done += 1
                    progress(done, len(species), start, len(failed))
        progress(done, len(species), start, len(failed), final=True)
        (work / "failed.txt").write_text("".join(f + "\n" for f in failed), encoding="utf-8")
    else:
        print("Validate-only mode (-V): download skipped.")

    print("\n=== Validation ===")
    errors = validate(db, meta, species, args.limit, args.jobs, failed)
    print()
    if args.limit > 0:
        if errors == 0:
            print(f"✓ PARTIAL TEST COMPLETED ({len(species)} genomes). Full catalogue not verified.")
            return 0
        print(f"Partial test finished with {errors} problem(s).")
        return 1
    if errors == 0:
        print(f"✓ Database ready: {db}\n  Add this folder in MetaPilot as MAG database; "
              "Marker4.fasta will be generated by MetaPilot.")
        return 0
    print(f"Finished with {errors} problem(s). Fix them and re-run (completed files are skipped).")
    return 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nInterrupted. Re-run to resume.")
        sys.exit(130)
