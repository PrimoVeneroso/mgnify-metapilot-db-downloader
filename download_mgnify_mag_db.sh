#!/usr/bin/env bash
# Download an MGnify genome catalogue (default: human-gut v3.0) in the layout
# required by MetaPilot "Add MAG database":
#
#   <OUT>/<catalogue>/<version>/
#       genomes-all_metadata.tsv          (columns Genome and Lineage)
#       original_db/<genomeID>.faa        (protein FASTA, decompressed)
#       eggNOG/<genomeID>_eggNOG.tsv      (eggNOG-mapper output, decompressed)
#
# Since v3.0 EBI files are .gz: they are decompressed because the MetaPilot
# validator only accepts plain .faa/.fasta and .tsv/.txt files.
# Marker4.fasta is NOT created here: MetaPilot generates it when the database is added.
#
# Usage: ./download_mgnify_mag_db.sh [-o DIR] [-c catalogue] [-v version] [-j parallel_downloads]
#        -n N download only the first N species (quick test)
#        -V   validate the local folder only, no download
# Re-runnable: complete files are skipped, only missing ones are retried.
#
# Author: Primo Veneroso. License: MIT.

set -euo pipefail

OUT_ROOT="./MGnify_MAG_DB"
CATALOGUE="human-gut"
VERSION="v3.0"
JOBS=8
LIMIT=0
VALIDATE_ONLY=0

usage() {
    sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

while getopts ":o:c:v:j:n:Vh" opt; do
    case "$opt" in
        o) OUT_ROOT="$OPTARG" ;;
        c) CATALOGUE="$OPTARG" ;;
        v) VERSION="$OPTARG" ;;
        j) JOBS="$OPTARG" ;;
        n) LIMIT="$OPTARG" ;;
        V) VALIDATE_ONLY=1 ;;
        h) usage 0 ;;
        *) usage 2 ;;
    esac
done

for tool in curl gzip awk xargs find sort; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: missing '$tool'" >&2; exit 2; }
done
[[ "$JOBS" =~ ^[0-9]+$ && "$JOBS" -ge 1 ]] || { echo "ERROR: -j must be an integer >= 1" >&2; exit 2; }

BASE_URL="https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/${CATALOGUE}/${VERSION}"
DB_DIR="${OUT_ROOT%/}/${CATALOGUE}/${VERSION}"
ORIG_DIR="$DB_DIR/original_db"
EGG_DIR="$DB_DIR/eggNOG"
WORK_DIR="$DB_DIR/.download_work"
META="$DB_DIR/genomes-all_metadata.tsv"
FAILED="$WORK_DIR/failed.txt"
PROGRESS_LOG="$WORK_DIR/progress.log"

mkdir -p "$ORIG_DIR" "$EGG_DIR" "$WORK_DIR"
: > "$FAILED"
: > "$PROGRESS_LOG"

echo "Catalogue: $BASE_URL"
echo "Folder   : $DB_DIR"
echo "Parallel : $JOBS"

# Download with retries and increasing back-off. Writes to .part first,
# so an interrupted file is never mistaken for a complete one.
fetch() {
    local url="$1" dest="$2" show_bar="${3:-0}" attempt
    local -a mode=(--silent --show-error)
    if [[ "$show_bar" == "1" && -t 1 ]]; then
        mode=(--progress-bar)
    fi
    for attempt in 1 2 3 4 5 6; do
        if curl --fail "${mode[@]}" --location \
                --connect-timeout 30 --speed-limit 1024 --speed-time 120 \
                --output "$dest.part" "$url"; then
            if [[ -s "$dest.part" ]]; then
                mv -f "$dest.part" "$dest"
                return 0
            fi
        fi
        rm -f "$dest.part"
        sleep $((attempt * 5))
    done
    return 1
}

# Verify, decompress and install a GZIP archive.
# On failure it never leaves a final file that looks valid.
unpack_gzip() {
    local gz="$1"
    local final="$2"
    local tmp="${2}.part"

    [[ -s "$gz" ]] || return 1

    if ! gzip -t -- "$gz" 2>/dev/null; then
        return 1
    fi

    if ! gzip -dc -- "$gz" > "$tmp"; then
        rm -f -- "$tmp"
        return 1
    fi

    if [[ ! -s "$tmp" ]]; then
        rm -f -- "$tmp"
        return 1
    fi

    if ! mv -f -- "$tmp" "$final"; then
        rm -f -- "$tmp"
        return 1
    fi

    rm -f -- "$gz"
    return 0
}

# Download and decompress a GZIP file.
fetch_gz() {
    local url="$1"
    local final="$2"
    local gz="${2}.gz"

    [[ -s "$final" ]] && return 0

    rm -f -- "$final.part" "$gz" "$gz.part"

    if ! fetch "$url" "$gz"; then
        return 1
    fi

    if unpack_gzip "$gz" "$final"; then
        return 0
    fi

    rm -f -- "$gz" "$gz.part" "$final.part"
    return 1
}

# One representative genome: protein FASTA + eggNOG annotation.
download_genome() {
    local id="$1"
    local parent="${id:0:11}"   # MGYG000000001 -> MGYG0000000
    local base="$BASE_URL/species_catalogue/$parent/$id/genome"
    fetch_gz "$base/$id.faa.gz"        "$ORIG_DIR/$id.faa"        || printf '%s faa\n' "$id" >> "$FAILED"
    fetch_gz "$base/${id}_eggNOG.tsv.gz" "$EGG_DIR/${id}_eggNOG.tsv" || printf '%s eggNOG\n' "$id" >> "$FAILED"
}

# Download a batch of genomes over a single TCP/TLS connection (HTTP keep-alive);
# any file that fails in the batch is retried individually.
download_batch() {
    local id parent base faa_final egg_final faa_gz egg_gz
    local -a curl_args=()
    for id in "$@"; do
        parent="${id:0:11}"
        base="$BASE_URL/species_catalogue/$parent/$id/genome"
        faa_final="$ORIG_DIR/$id.faa"
        egg_final="$EGG_DIR/${id}_eggNOG.tsv"
        faa_gz="$faa_final.gz.part"
        egg_gz="$egg_final.gz.part"
        [[ ! -s "$faa_final" ]] && curl_args+=(-o "$faa_gz" "$base/$id.faa.gz")
        [[ ! -s "$egg_final" ]] && curl_args+=(-o "$egg_gz" "$base/${id}_eggNOG.tsv.gz")
    done

    if [[ "${#curl_args[@]}" -gt 0 ]]; then
        curl --fail --silent --show-error --location \
            --connect-timeout 30 --speed-limit 1024 --speed-time 120 \
            --retry 3 --retry-delay 2 "${curl_args[@]}" || true
    fi

    for id in "$@"; do
        parent="${id:0:11}"
        base="$BASE_URL/species_catalogue/$parent/$id/genome"
        faa_final="$ORIG_DIR/$id.faa"
        egg_final="$EGG_DIR/${id}_eggNOG.tsv"
        faa_gz="$faa_final.gz.part"
        egg_gz="$egg_final.gz.part"

        if [[ ! -s "$faa_final" ]]; then
            if ! unpack_gzip "$faa_gz" "$faa_final"; then
                rm -f -- "$faa_gz" "$faa_final.part"
                if ! fetch_gz "$base/$id.faa.gz" "$faa_final"; then
                    printf '%s faa\n' "$id" >> "$FAILED"
                fi
            fi
        fi

        if [[ ! -s "$egg_final" ]]; then
            if ! unpack_gzip "$egg_gz" "$egg_final"; then
                rm -f -- "$egg_gz" "$egg_final.part"
                if ! fetch_gz "$base/${id}_eggNOG.tsv.gz" "$egg_final"; then
                    printf '%s eggNOG\n' "$id" >> "$FAILED"
                fi
            fi
        fi

        echo 1 >> "$PROGRESS_LOG"
    done
}
export -f fetch unpack_gzip fetch_gz download_genome download_batch
export BASE_URL ORIG_DIR EGG_DIR FAILED PROGRESS_LOG

# Draw the progress bar (uses \r on interactive terminals).
render_progress() {
    local done_count="$1" total="$2" start_ts="$3" active="$4" err_count="$5" final_newline="${6:-0}"
    if (( done_count > total )); then done_count="$total"; fi
    local width=26
    local filled=$(( total > 0 ? done_count * width / total : width ))
    local empty=$(( width - filled ))
    local bar_fill bar_empty
    printf -v bar_fill '%*s' "$filled" ''
    printf -v bar_empty '%*s' "$empty" ''
    bar_fill="${bar_fill// /#}"
    bar_empty="${bar_empty// /-}"

    local now elapsed pct_int pct_dec rate_str eta_str
    now=$(date +%s)
    elapsed=$(( now - start_ts ))
    if (( elapsed < 1 )); then elapsed=1; fi
    pct_int=$(( total > 0 ? done_count * 100 / total : 100 ))
    pct_dec=$(( total > 0 ? (done_count * 1000 / total) % 10 : 0 ))

    if (( done_count > 0 && done_count < total )); then
        local rem=$(( total - done_count ))
        local eta_sec=$(( rem * elapsed / done_count ))
        printf -v eta_str "%02d:%02d" $(( eta_sec / 60 )) $(( eta_sec % 60 ))
    elif (( done_count >= total )); then
        eta_str="00:00"
    else
        eta_str="--:--"
    fi
    rate_str=$(LC_NUMERIC=C awk -v d="$done_count" -v e="$elapsed" 'BEGIN { printf "%.1f", d / e }')

    if [[ -t 1 ]]; then
        printf "\r\033[K[%s%s] %3d.%d%% (%d/%d) | %s gen/s | ETA %s | active:%d | err:%d" \
            "$bar_fill" "$bar_empty" "$pct_int" "$pct_dec" "$done_count" "$total" \
            "$rate_str" "$eta_str" "$active" "$err_count"
        if [[ "$final_newline" == "1" ]]; then printf "\n"; fi
    elif [[ "$final_newline" == "1" ]]; then
        printf "[%s%s] %3d.%d%% (%d/%d) | %s gen/s | ETA %s | err:%d\n" \
            "$bar_fill" "$bar_empty" "$pct_int" "$pct_dec" "$done_count" "$total" \
            "$rate_str" "$eta_str" "$err_count"
    fi
    return 0
}

monitor_progress() {
    set +e
    local total="$1" start_ts="$2" last_log=0
    while true; do
        local done_count active err_count now
        done_count=$(wc -l < "$PROGRESS_LOG" 2>/dev/null || echo 0)
        active=$(find "$ORIG_DIR" "$EGG_DIR" -maxdepth 1 -name '*.part' 2>/dev/null | wc -l)
        err_count=$(wc -l < "$FAILED" 2>/dev/null || echo 0)
        render_progress "$done_count" "$total" "$start_ts" "$active" "$err_count" 0
        now=$(date +%s)
        if [[ ! -t 1 ]] && (( now - last_log >= 15 )); then
            render_progress "$done_count" "$total" "$start_ts" "$active" "$err_count" 1
            last_log="$now"
        fi
        LC_NUMERIC=C sleep 0.5
    done
}

# 1) Taxonomic metadata (~110 MB). The species list is read from here,
#    avoiding crawling thousands of server folders over HTTP.
if [[ ! -s "$META" ]]; then
    if [[ "$VALIDATE_ONLY" -eq 1 ]]; then
        echo "ERROR: metadata missing or empty: $META" >&2
        exit 1
    fi
    echo "Download genomes-all_metadata.tsv (~110 MB) ..."
    fetch "$BASE_URL/genomes-all_metadata.tsv" "$META" 1 || { echo "ERROR: metadata download failed" >&2; exit 1; }
fi

SPECIES_LIST="$WORK_DIR/species_reps.txt"
awk -F'\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) { gsub(/\r/, "", $i); col[$i] = i }
        if (!("Genome" in col) || !("Lineage" in col) || !("Species_rep" in col)) {
            print "ERROR: columns Genome/Lineage/Species_rep missing" > "/dev/stderr"
            exit 3
        }
        next
    }
    {
        id = $(col["Species_rep"])
        gsub(/\r/, "", id)
        if (id ~ /^MGYG[0-9]+(\.[0-9]+)?$/) print id
    }
' "$META" | LC_ALL=C sort -u > "$SPECIES_LIST"
if [[ "$LIMIT" -gt 0 ]]; then
    head -n "$LIMIT" "$SPECIES_LIST" > "$SPECIES_LIST.tmp" && mv -f "$SPECIES_LIST.tmp" "$SPECIES_LIST"
fi

TOTAL=$(wc -l < "$SPECIES_LIST")
[[ "$TOTAL" -gt 0 ]] || { echo "ERROR: no species found in metadata" >&2; exit 1; }
echo "Representative species: $TOTAL"

# 2) Parallel download in batches of 5 genomes per keep-alive connection,
#    with a live progress bar (skipped with -V).
if [[ "$VALIDATE_ONLY" -eq 0 ]]; then
    echo "Downloading FASTA and eggNOG ..."
    START_TS=$(date +%s)
    monitor_progress "$TOTAL" "$START_TS" &
    MONITOR_PID=$!
    trap 'kill "$MONITOR_PID" 2>/dev/null || true' EXIT

    xargs -P "$JOBS" -n 5 bash -c 'download_batch "$@"' _ < "$SPECIES_LIST"

    kill "$MONITOR_PID" 2>/dev/null || true
    wait "$MONITOR_PID" 2>/dev/null || true
    trap - EXIT
    DONE_FINAL=$(wc -l < "$PROGRESS_LOG" 2>/dev/null || echo "$TOTAL")
    ERR_FINAL=$(wc -l < "$FAILED" 2>/dev/null || echo 0)
    render_progress "$DONE_FINAL" "$TOTAL" "$START_TS" 0 "$ERR_FINAL" 1
else
    echo "Validate-only mode (-V): download skipped."
fi

# 3) Full validation against the MetaPilot import rules.
echo
echo "=== Validation ==="
ERRORS=0

# Rule 6: Genome and Lineage filled on every row of genomes-all_metadata.tsv.
if ! awk -F'\t' '
    NR == 1 { for (i = 1; i <= NF; i++) { gsub(/\r/, "", $i); c[$i] = i }; next }
    {
        gsub(/\r/, "", $c["Genome"]); gsub(/\r/, "", $c["Lineage"])
        if ($c["Genome"] == "" || $c["Lineage"] == "") bad++
    }
    END { exit (bad > 0 ? 1 : 0) }
' "$META"; then
    echo "✗ Rule 6: some rows of genomes-all_metadata.tsv have empty Genome or Lineage"
    ERRORS=$((ERRORS + 1))
fi

N_FAA=$(find "$ORIG_DIR" -maxdepth 1 -name '*.faa' -type f | wc -l)
N_EGG=$(find "$EGG_DIR" -maxdepth 1 -name '*_eggNOG.tsv' -type f | wc -l)
echo "original_db: $N_FAA files   eggNOG: $N_EGG files   expected: $TOTAL"
[[ "$N_FAA" -eq "$N_EGG" ]] || { echo "✗ Rule 1: file counts differ"; ERRORS=$((ERRORS + 1)); }

# Unexpected extensions or leftover .part/.gz files in the final folders.
UNEXPECTED_ORIG=$(find "$ORIG_DIR" -maxdepth 1 -type f ! -name '*.faa' | wc -l)
UNEXPECTED_EGG=$(find "$EGG_DIR" -maxdepth 1 -type f ! -name '*_eggNOG.tsv' | wc -l)
if [[ "$UNEXPECTED_ORIG" -gt 0 || "$UNEXPECTED_EGG" -gt 0 ]]; then
    echo "✗ Rule 2: unexpected or leftover .part/.gz files ($UNEXPECTED_ORIG in original_db, $UNEXPECTED_EGG in eggNOG)"
    find "$ORIG_DIR" -maxdepth 1 -type f ! -name '*.faa' | head -n 10
    find "$EGG_DIR" -maxdepth 1 -type f ! -name '*_eggNOG.tsv' | head -n 10
    ERRORS=$((ERRORS + 1))
fi

# Empty (0 byte) files.
EMPTY_FILES=$(find "$ORIG_DIR" "$EGG_DIR" -maxdepth 1 -type f -size 0 | wc -l)
if [[ "$EMPTY_FILES" -gt 0 ]]; then
    echo "✗ Empty (0 byte) files found: $EMPTY_FILES"
    find "$ORIG_DIR" "$EGG_DIR" -maxdepth 1 -type f -size 0 | head -n 10
    ERRORS=$((ERRORS + 1))
fi

ACTUAL_FAA="$WORK_DIR/actual_faa.txt"
ACTUAL_EGG="$WORK_DIR/actual_egg.txt"
find "$ORIG_DIR" -maxdepth 1 -type f -name '*.faa' -printf '%f\n' |
    sed 's/\.faa$//' | LC_ALL=C sort -u > "$ACTUAL_FAA"
find "$EGG_DIR" -maxdepth 1 -type f -name '*_eggNOG.tsv' -printf '%f\n' |
    sed 's/_eggNOG\.tsv$//' | LC_ALL=C sort -u > "$ACTUAL_EGG"

# Rule 3: genomeID pairing.
UNPAIRED=$(LC_ALL=C comm -3 "$ACTUAL_FAA" "$ACTUAL_EGG" | wc -l)
[[ "$UNPAIRED" -eq 0 ]] || { echo "✗ Rule 3: $UNPAIRED genomes without a pair"; ERRORS=$((ERRORS + 1)); }

# Compare expected IDs with the IDs on disk.
if [[ "$LIMIT" -eq 0 ]]; then
    if ! cmp -s "$SPECIES_LIST" "$ACTUAL_FAA"; then
        echo "✗ ERROR: FASTA IDs do not match the catalogue."
        echo "  Genomes without FASTA:"
        LC_ALL=C comm -23 "$SPECIES_LIST" "$ACTUAL_FAA" | head -n 20
        echo "  Unexpected FASTA:"
        LC_ALL=C comm -13 "$SPECIES_LIST" "$ACTUAL_FAA" | head -n 20
        ERRORS=$((ERRORS + 1))
    fi
    if ! cmp -s "$SPECIES_LIST" "$ACTUAL_EGG"; then
        echo "✗ ERROR: eggNOG IDs do not match the catalogue."
        echo "  Genomes without eggNOG:"
        LC_ALL=C comm -23 "$SPECIES_LIST" "$ACTUAL_EGG" | head -n 20
        echo "  Unexpected eggNOG files:"
        LC_ALL=C comm -13 "$SPECIES_LIST" "$ACTUAL_EGG" | head -n 20
        ERRORS=$((ERRORS + 1))
    fi
else
    MISSING_LIMIT_FAA=$(LC_ALL=C comm -23 "$SPECIES_LIST" "$ACTUAL_FAA" | wc -l)
    MISSING_LIMIT_EGG=$(LC_ALL=C comm -23 "$SPECIES_LIST" "$ACTUAL_EGG" | wc -l)
    if [[ "$MISSING_LIMIT_FAA" -gt 0 || "$MISSING_LIMIT_EGG" -gt 0 ]]; then
        echo "✗ ERROR: files missing for the genomes selected with -n ($MISSING_LIMIT_FAA FASTA, $MISSING_LIMIT_EGG eggNOG)."
        ERRORS=$((ERRORS + 1))
    fi
fi

# Rule 4: >genomeID_proteinID name headers, sequences present, valid FASTA characters.
BAD_HEADERS="$WORK_DIR/bad_headers.txt"
find "$ORIG_DIR" -maxdepth 1 -name '*.faa' -type f -size +0c -print0 | xargs -0 -r -P "$JOBS" awk '
    FNR == 1 {
        if (prev != "" && !has_sequence) print prev ": no valid FASTA sequence"
        prev = FILENAME
        g = FILENAME; sub(/.*\//, "", g); sub(/\.faa$/, "", g)
        seen_header = 0; has_sequence = 0
    }
    /^>/ {
        seen_header = 1
        h = substr($0, 2); sub(/^[ \t]+/, "", h); sub(/[ \t\r]+$/, "", h)
        if (index(h, g "_") != 1) { print FILENAME ": genomeID mismatch: " h; next }
        rest = substr(h, length(g) + 2); sp = index(rest, " ")
        if (sp <= 1) { print FILENAME ": missing proteinID or name: " h; next }
        name = substr(rest, sp + 1); gsub(/^[ \t]+|[ \t]+$/, "", name)
        if (name == "") print FILENAME ": missing protein name: " h
        next
    }
    /^[ \t\r]*$/ { next }
    {
        sub(/\r$/, "", $0)
        if (!seen_header) print FILENAME ": sequence before first header at line " FNR
        if ($0 !~ /^[A-Za-z*. -]+$/) print FILENAME ": unexpected characters in sequence at line " FNR
        else if (seen_header) has_sequence = 1
    }
    END {
        if (prev != "" && !has_sequence) print prev ": no valid FASTA sequence"
    }
' > "$BAD_HEADERS"
if [[ -s "$BAD_HEADERS" ]]; then
    echo "✗ Rule 4: $(wc -l < "$BAD_HEADERS") FASTA problems (see $BAD_HEADERS)"
    ERRORS=$((ERRORS + 1))
fi

# Rule 5: eggNOG tab-separated, non-empty, rows consistent with the header, queries start with <genomeID>_.
BAD_EGG_FILE="$WORK_DIR/bad_eggnog.txt"
find "$EGG_DIR" -maxdepth 1 -name '*_eggNOG.tsv' -type f -size +0c -print0 | xargs -0 -r -P "$JOBS" awk -F'\t' '
    FNR == 1 {
        if (prev != "" && data_rows == 0 && !bad_file) print prev ": no data rows besides headers"
        prev = FILENAME
        g = FILENAME; sub(/.*\//, "", g); sub(/_eggNOG\.tsv$/, "", g)
        hdr_nf = NF; data_rows = 0; bad_file = 0
        if (NF < 2) { print FILENAME ": header not tab-separated (columns=" NF ")"; bad_file = 1 }
        next
    }
    bad_file { next }
    /^[ \t]*#/ || /^[ \t\r]*$/ { next }
    {
        data_rows++
        if (NF != hdr_nf) {
            print FILENAME ":" FNR ": " NF " columns, header has " hdr_nf
            bad_file = 1
        } else if (index($1, g "_") != 1) {
            print FILENAME ":" FNR ": query ID does not start with " g "_ (" $1 ")"
            bad_file = 1
        }
    }
    END {
        if (prev != "" && data_rows == 0 && !bad_file) print prev ": no data rows besides headers"
    }
' > "$BAD_EGG_FILE"
if [[ -s "$BAD_EGG_FILE" ]]; then
    echo "✗ Rule 5: $(wc -l < "$BAD_EGG_FILE") invalid eggNOG files (see $BAD_EGG_FILE)"
    ERRORS=$((ERRORS + 1))
fi

if [[ -s "$FAILED" ]]; then
    echo "✗ Failed downloads: $(wc -l < "$FAILED") (see $FAILED). Re-run the script to retry."
    ERRORS=$((ERRORS + 1))
fi

echo
# -n means a partial test, not a complete database.
if [[ "$LIMIT" -gt 0 ]]; then
    if [[ "$ERRORS" -eq 0 ]]; then
        echo "✓ PARTIAL TEST COMPLETED."
        echo "  Selected genomes: $TOTAL"
        echo "  FASTA files: $N_FAA"
        echo "  eggNOG files: $N_EGG"
        echo "  The full catalogue was not verified."
        exit 0
    fi
    echo "Partial test finished with problems: $ERRORS"
    exit 1
fi

# Only a complete download can be declared ready.
if [[ "$ERRORS" -eq 0 && "$N_FAA" -eq "$TOTAL" ]]; then
    echo "✓ Database ready: $DB_DIR"
    echo "  In MetaPilot, add this folder as MAG database."
    echo "  Marker4.fasta will be generated by MetaPilot."
    exit 0
fi

echo "Finished with problems: $ERRORS"
echo "See $FAILED and the validation reports in $WORK_DIR."
echo "Re-run the script after fixing the problems."
exit 1
