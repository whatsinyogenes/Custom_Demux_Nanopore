#!/usr/bin/env bash

set -Eeuo pipefail


###############################################################################
# nanopore custom demultiplexing pipeline
#
# Workflow:
#
#   0. Validate inputs
#   1. Validate sequencing manifest
#   2. Reduce index manifest to UDPs actually sequenced
#   3. Generate Dorado TOMLs + FASTAs
#   4. Dorado i7 forward
#   5. Dorado i7 reverse
#   6. Dorado i5 forward
#   7. Dorado i5 reverse
#   8. Reconstruct sample FASTQs
#   9. Fastplong trimming × N
#  10. minimap2 -ax sr
#  11. samtools sort/index
#
###############################################################################


###############################################################################
# Defaults
###############################################################################

DORADO="dorado"
FASTPLONG="fastplong"
MINIMAP2="minimap2"
SAMTOOLS="samtools"
PYTHON="python3"

THREADS=4

FASTPLONG_PASSES=3
DISTANCE_THRESHOLD="0.10"
TRIMMING_EXTENSION=0
MIN_LENGTH=50

#change minimap preset if for longreads etc 
MINIMAP2_PRESET="sr"

# Pipeline controls
RESUME=true
FORCE=false
KEEP_INTERMEDIATE=true

SKIP_DEMUX=false
SKIP_TRIMMING=false
SKIP_ALIGNMENT=false


###############################################################################
# Dorado scoring defaults
###############################################################################

MAX_BARCODE_PENALTY=11
BARCODE_END_PROXIMITY=1000
MIN_BARCODE_PENALTY_DIST=3
MIN_SEPARATION_ONLY_DIST=6

FLANK_LEFT_PAD=5
FLANK_RIGHT_PAD=10

FRONT_BARCODE_WINDOW=1000
REAR_BARCODE_WINDOW=1000

MIDSTRAND_FLANK_SCORE="0.90"

REAR_ONLY_BARCODES="false"


###############################################################################
# Barcode flank sequences
###############################################################################

I7_FRONT="CAAGCAGAAGACGGCATACGAGAT"
I7_REAR="GTCTCGTGGGCTCGG"

I5_FRONT="AATGATACGGCGACCACCGAGATCTACAC"
I5_REAR="GTCGGCAGCGTC"


###############################################################################
# Variables populated by arguments
###############################################################################

INPUT=""
INDEXES=""
SAMPLES_MANIFEST=""
REFERENCE=""
ADAPTERS=""
OUTPUT=""


###############################################################################
# Usage
###############################################################################

usage() {
    cat <<'EOF'

Nano_process_pipeline

Usage:

  bash Nano_process_pipeline \
      --input FASTQ_OR_DIRECTORY \
      --indexes INDEX_CSV \
      --samples SEQUENCING_MANIFEST_TSV \
      --reference REFERENCE_FASTA \
      --adapters FASTPLONG_ADAPTER_FASTA \
      --output OUTPUT_DIRECTORY


Required arguments:

  --input PATH
        FASTQ file or directory containing FASTQ/FASTQ.GZ files.

  --indexes CSV
        Master index manifest CSV containing all UDP barcode sequences.

  --samples TSV
        Sequencing manifest specifying which UDPs were actually used.

        Format:

            sample_name    index_pair
            patient1       UDP0001
            patient2       UDP0003
            patient3       UDP0007

  --reference FASTA
        Reference genome FASTA for minimap2.

  --adapters FASTA
        Adapter FASTA for Fastplong.

  --output DIRECTORY
        Output directory.


General options:

  --threads N
        Number of threads.
        Default: 4

  --dorado PATH
        Dorado executable.
        Default: dorado

  --fastplong PATH
        Fastplong executable.
        Default: fastplong

  --minimap2 PATH
        minimap2 executable.
        Default: minimap2

  --samtools PATH
        samtools executable.
        Default: samtools

  --python PATH
        Python executable.
        Default: python3


Fastplong:

  --fastplong-passes N
        Number of Fastplong passes.
        Default: 3

  --distance-threshold VALUE
        Fastplong distance threshold.
        Default: 0.10

  --trimming-extension N
        Fastplong trimming extension.
        Default: 0

  --min-length N
        Minimum retained read length.
        Default: 50


Alignment:

  --minimap2-preset PRESET
        minimap2 preset.
        Default: sr


Dorado scoring:

  --max-barcode-penalty N
  --barcode-end-proximity N
  --min-barcode-penalty-dist N
  --min-separation-only-dist N
  --flank-left-pad N
  --flank-right-pad N
  --front-barcode-window N
  --rear-barcode-window N
  --midstrand-flank-score VALUE
  --rear-only-barcodes true|false


Barcode flanks:

  --i7-front SEQUENCE
  --i7-rear SEQUENCE
  --i5-front SEQUENCE
  --i5-rear SEQUENCE


Pipeline control:

  --resume
        Resume completed stages.
        Default: enabled

  --no-resume
        Do not use .done markers.

  --force
        Re-run stages even if .done markers exist.
        Existing stage output is cleaned before rerunning.

  --keep-intermediate
        Keep intermediate files.
        Default: enabled

  --no-keep-intermediate
        Remove intermediate files when possible.

  --skip-demux
        Skip Dorado demultiplexing.

  --skip-trimming
        Skip Fastplong.

  --skip-alignment
        Skip minimap2/samtools.


Help:

  -h
  --help


Example:

  bash Nano_process_pipeline \
      --input fastqs/ \
      --indexes tagmentation_index_manifest_v2.csv \
      --samples sequencing_manifest.tsv \
      --reference hg38.fa \
      --adapters adaptors_for_fastplong.fasta \
      --output Nano_processed/

EOF
}


###############################################################################
# Error handling
###############################################################################

die() {
    echo
    echo "ERROR: $*" >&2
    exit 1
}


trap 'echo; echo "ERROR: Pipeline failed at line $LINENO." >&2' ERR


###############################################################################
# Argument parser
###############################################################################

while [[ $# -gt 0 ]]; do

    case "$1" in

        --input)
            INPUT="$2"
            shift 2
            ;;

        --indexes)
            INDEXES="$2"
            shift 2
            ;;

        --samples)
            SAMPLES_MANIFEST="$2"
            shift 2
            ;;

        --reference)
            REFERENCE="$2"
            shift 2
            ;;

        --adapters)
            ADAPTERS="$2"
            shift 2
            ;;


        --output)
            OUTPUT="$2"
            shift 2
            ;;


        --threads)
            THREADS="$2"
            shift 2
            ;;

        --dorado)
            DORADO="$2"
            shift 2
            ;;

        --fastplong)
            FASTPLONG="$2"
            shift 2
            ;;

        --minimap2)
            MINIMAP2="$2"
            shift 2
            ;;

        --samtools)
            SAMTOOLS="$2"
            shift 2
            ;;

        --python)
            PYTHON="$2"
            shift 2
            ;;


        --fastplong-passes)
            FASTPLONG_PASSES="$2"
            shift 2
            ;;

        --distance-threshold)
            DISTANCE_THRESHOLD="$2"
            shift 2
            ;;

        --trimming-extension)
            TRIMMING_EXTENSION="$2"
            shift 2
            ;;

        --min-length)
            MIN_LENGTH="$2"
            shift 2
            ;;


        --minimap2-preset)
            MINIMAP2_PRESET="$2"
            shift 2
            ;;


        --max-barcode-penalty)
            MAX_BARCODE_PENALTY="$2"
            shift 2
            ;;

        --barcode-end-proximity)
            BARCODE_END_PROXIMITY="$2"
            shift 2
            ;;

        --min-barcode-penalty-dist)
            MIN_BARCODE_PENALTY_DIST="$2"
            shift 2
            ;;

        --min-separation-only-dist)
            MIN_SEPARATION_ONLY_DIST="$2"
            shift 2
            ;;

        --flank-left-pad)
            FLANK_LEFT_PAD="$2"
            shift 2
            ;;

        --flank-right-pad)
            FLANK_RIGHT_PAD="$2"
            shift 2
            ;;

        --front-barcode-window)
            FRONT_BARCODE_WINDOW="$2"
            shift 2
            ;;

        --rear-barcode-window)
            REAR_BARCODE_WINDOW="$2"
            shift 2
            ;;

        --midstrand-flank-score)
            MIDSTRAND_FLANK_SCORE="$2"
            shift 2
            ;;

        --rear-only-barcodes)
            REAR_ONLY_BARCODES="$2"
            shift 2
            ;;


        --i7-front)
            I7_FRONT="$2"
            shift 2
            ;;

        --i7-rear)
            I7_REAR="$2"
            shift 2
            ;;

        --i5-front)
            I5_FRONT="$2"
            shift 2
            ;;

        --i5-rear)
            I5_REAR="$2"
            shift 2
            ;;


        --resume)
            RESUME=true
            shift
            ;;

        --no-resume)
            RESUME=false
            shift
            ;;

        --force)
            FORCE=true
            shift
            ;;

        --keep-intermediate)
            KEEP_INTERMEDIATE=true
            shift
            ;;

        --no-keep-intermediate)
            KEEP_INTERMEDIATE=false
            shift
            ;;

        --skip-demux)
            SKIP_DEMUX=true
            shift
            ;;

        --skip-trimming)
            SKIP_TRIMMING=true
            shift
            ;;

        --skip-alignment)
            SKIP_ALIGNMENT=true
            shift
            ;;


        -h|--help)
            usage
            exit 0
            ;;

        *)
            die "Unknown argument: $1. Use --help."
            ;;

    esac

done


###############################################################################
# Validate required arguments
###############################################################################

[[ -n "$INPUT" ]] \
    || die "Missing --input"

[[ -n "$INDEXES" ]] \
    || die "Missing --indexes"

[[ -n "$SAMPLES_MANIFEST" ]] \
    || die "Missing --samples"

[[ -n "$REFERENCE" ]] \
    || die "Missing --reference"

[[ -n "$ADAPTERS" ]] \
    || die "Missing --adapters"

[[ -n "$OUTPUT" ]] \
    || die "Missing --output"


[[ -e "$INPUT" ]] \
    || die "Input does not exist: $INPUT"

[[ -f "$INDEXES" ]] \
    || die "Index CSV does not exist: $INDEXES"

[[ -f "$SAMPLES_MANIFEST" ]] \
    || die "Sequencing manifest does not exist: $SAMPLES_MANIFEST"

[[ -f "$REFERENCE" ]] \
    || die "Reference does not exist: $REFERENCE"

[[ -f "$ADAPTERS" ]] \
    || die "Adapter FASTA does not exist: $ADAPTERS"


[[ "$THREADS" =~ ^[0-9]+$ ]] \
    || die "--threads must be an integer"

[[ "$FASTPLONG_PASSES" =~ ^[0-9]+$ ]] \
    || die "--fastplong-passes must be an integer"

[[ "$MIN_LENGTH" =~ ^[0-9]+$ ]] \
    || die "--min-length must be an integer"


###############################################################################
# Normalize paths
###############################################################################

INPUT="$(realpath "$INPUT")"
INDEXES="$(realpath "$INDEXES")"
SAMPLES_MANIFEST="$(realpath "$SAMPLES_MANIFEST")"
REFERENCE="$(realpath "$REFERENCE")"
ADAPTERS="$(realpath "$ADAPTERS")"
OUTPUT="$(realpath -m "$OUTPUT")"


###############################################################################
# Output directories
###############################################################################

GENERATED="${OUTPUT}/generated"
TOMLS="${GENERATED}/tomls"
FASTAS="${GENERATED}/fastas"

DEMUXED="${OUTPUT}/demuxed"
SAMPLES="${OUTPUT}/samples"
TRIMMED="${OUTPUT}/trimmed"
ALIGNED="${OUTPUT}/aligned"
LOGS="${OUTPUT}/logs"
DONE="${OUTPUT}/.done"

REDUCED_INDEXES="${GENERATED}/sequenced_indexes.csv"
SAMPLE_MAP="${GENERATED}/sample_map.tsv"

REFERENCE_MMI="${OUTPUT}/reference.mmi"

SUMMARY="${OUTPUT}/summary.tsv"


###############################################################################
# Create directories
###############################################################################

mkdir -p \
    "$OUTPUT" \
    "$GENERATED" \
    "$TOMLS" \
    "$FASTAS" \
    "$DEMUXED" \
    "$SAMPLES" \
    "$TRIMMED" \
    "$ALIGNED" \
    "$LOGS" \
    "$DONE"


###############################################################################
# Check executables
###############################################################################

command -v "$DORADO" >/dev/null 2>&1 \
    || die "Dorado not found: $DORADO"

command -v "$FASTPLONG" >/dev/null 2>&1 \
    || die "Fastplong not found: $FASTPLONG"

command -v "$MINIMAP2" >/dev/null 2>&1 \
    || die "minimap2 not found: $MINIMAP2"

command -v "$SAMTOOLS" >/dev/null 2>&1 \
    || die "samtools not found: $SAMTOOLS"

command -v "$PYTHON" >/dev/null 2>&1 \
    || die "Python not found: $PYTHON"


###############################################################################
# Logging
###############################################################################

LOGFILE="${LOGS}/pipeline.log"

exec > >(tee -a "$LOGFILE") 2>&1


###############################################################################
# Utility functions
###############################################################################

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}


log() {
    echo "[$(timestamp)] $*"
}


stage_done() {

    local marker="$1"

    if [[ "$FORCE" == true ]]; then
        return 1
    fi

    if [[ "$RESUME" == true && -f "$marker" ]]; then
        return 0
    fi

    return 1
}


mark_done() {

    local marker="$1"

    mkdir -p "$(dirname "$marker")"
    touch "$marker"
}


###############################################################################
# Find FASTQ files
###############################################################################

find_fastqs() {

    local source="$1"

    if [[ -f "$source" ]]; then

        case "$source" in
            *.fastq|*.fq|*.fastq.gz|*.fq.gz)
                printf '%s\n' "$source"
                ;;
            *)
                die "Input file does not look like FASTQ: $source"
                ;;
        esac

        return
    fi


    find "$source" -type f \
        \( \
            -name '*.fastq' \
            -o -name '*.fq' \
            -o -name '*.fastq.gz' \
            -o -name '*.fq.gz' \
        \) \
        -print0 |
        sort -z |
        while IFS= read -r -d '' f; do
            printf '%s\n' "$f"
        done
}


mapfile -t FASTQS < <(find_fastqs "$INPUT")


[[ "${#FASTQS[@]}" -gt 0 ]] \
    || die "No FASTQ files found in $INPUT"


log "Found ${#FASTQS[@]} FASTQ file(s)."


###############################################################################
# Validate index manifest
###############################################################################

validate_index_manifest() {

    local csv="$1"

    "$PYTHON" - "$csv" <<'PY'
import csv
import re
import sys

path = sys.argv[1]

required = [
    "Index Name",
    "i7 Bases in Adapter",
    "i7 Bases for Sample Sheet",
    "i5 Bases in Adapter",
    "i5 Bases for Sample Sheet in Forward Orientation",
    "I5 Bases for Sample Sheet in Reverse Complement Orientation",
]

dna = re.compile(r"^[ACGT]+$", re.IGNORECASE)

with open(path, newline="") as fh:
    reader = csv.DictReader(fh)

    if reader.fieldnames is None:
        raise SystemExit("CSV has no header.")

    missing = [x for x in required if x not in reader.fieldnames]

    if missing:
        raise SystemExit(
            "CSV is missing required columns:\n  " +
            "\n  ".join(missing)
        )

    names = set()
    rows = 0

    for line_number, row in enumerate(reader, start=2):

        rows += 1

        name = row["Index Name"].strip()

        if not name:
            raise SystemExit(
                f"Empty Index Name on CSV line {line_number}"
            )

        if name in names:
            raise SystemExit(
                f"Duplicate Index Name '{name}' on CSV line {line_number}"
            )

        names.add(name)

        for column in required[1:]:

            seq = row[column].strip()

            if not seq:
                raise SystemExit(
                    f"Empty sequence in '{column}' "
                    f"on CSV line {line_number}"
                )

            if not dna.fullmatch(seq):
                raise SystemExit(
                    f"Invalid DNA sequence in '{column}' "
                    f"on CSV line {line_number}: {seq}"
                )

if rows == 0:
    raise SystemExit("CSV contains no samples.")

print(f"Validated {rows} index/sample rows.")
PY
}


###############################################################################
# Validate sequencing manifest and create reduced index manifest
###############################################################################

prepare_sequencing_manifest() {

    local marker="${DONE}/sequencing_manifest.done"

    if stage_done "$marker" && \
       [[ -f "$REDUCED_INDEXES" ]] && \
       [[ -f "$SAMPLE_MAP" ]]; then

        log "Sequencing manifest already processed. Skipping."
        return
    fi


    log "Validating sequencing manifest:"
    log "  ${SAMPLES_MANIFEST}"


    "$PYTHON" - \
        "$INDEXES" \
        "$SAMPLES_MANIFEST" \
        "$REDUCED_INDEXES" \
        <<'PY'

import csv
import sys

index_csv = sys.argv[1]
sample_tsv = sys.argv[2]
reduced_csv = sys.argv[3]


###############################################################################
# Read master index manifest
###############################################################################

with open(index_csv, newline="") as fh:
    index_rows = list(csv.DictReader(fh))

if not index_rows:
    raise SystemExit("Master index CSV contains no rows.")


index_lookup = {}

for row in index_rows:

    udp = row["Index Name"].strip()

    if udp in index_lookup:
        raise SystemExit(
            f"Duplicate UDP in master index CSV: {udp}"
        )

    index_lookup[udp] = row


###############################################################################
# Read sequencing manifest
###############################################################################

with open(sample_tsv, newline="") as fh:

    reader = csv.DictReader(fh, delimiter="\t")

    if reader.fieldnames is None:
        raise SystemExit(
            "Sequencing manifest has no header."
        )

    required = [
        "sample_name",
        "index_pair",
    ]

    missing = [
        x for x in required
        if x not in reader.fieldnames
    ]

    if missing:
        raise SystemExit(
            "Sequencing manifest is missing required columns:\n  " +
            "\n  ".join(missing)
        )

    sample_rows = list(reader)


if not sample_rows:
    raise SystemExit(
        "Sequencing manifest contains no samples."
    )


###############################################################################
# Validate sample names and UDPs
###############################################################################

sample_names = set()
udp_names = set()

selected_rows = []

for line_number, row in enumerate(sample_rows, start=2):

    sample = row["sample_name"].strip()
    udp = row["index_pair"].strip()

    if not sample:
        raise SystemExit(
            f"Empty sample_name on sequencing manifest line "
            f"{line_number}"
        )

    if not udp:
        raise SystemExit(
            f"Empty index_pair on sequencing manifest line "
            f"{line_number}"
        )

    if sample in sample_names:
        raise SystemExit(
            f"Duplicate sample_name '{sample}' "
            f"on sequencing manifest line {line_number}"
        )

    if udp in udp_names:
        raise SystemExit(
            f"UDP '{udp}' is assigned more than once "
            f"in sequencing manifest. "
            f"Each UDP must map to one sample."
        )

    if udp not in index_lookup:
        raise SystemExit(
            f"UDP '{udp}' from sequencing manifest "
            f"is not present in master index CSV."
        )

    sample_names.add(sample)
    udp_names.add(udp)

    selected_rows.append(index_lookup[udp])


###############################################################################
# Write reduced index CSV
###############################################################################

with open(reduced_csv, "w", newline="") as fh:

    fieldnames = index_rows[0].keys()

    writer = csv.DictWriter(
        fh,
        fieldnames=fieldnames
    )

    writer.writeheader()

    writer.writerows(selected_rows)


###############################################################################
# Report
###############################################################################

print()
print("SEQUENCING MANIFEST")
print("===================")
print(f"Samples specified: {len(sample_rows)}")
print(f"UDP indexes used:  {len(selected_rows)}")
print()

for row in sample_rows:
    print(
        f"  {row['sample_name'].strip()}"
        f"\t{row['index_pair'].strip()}"
    )

print()
print(f"Reduced index manifest: {reduced_csv}")
print()

PY


    ###########################################################################
    # Create sample map
    #
    # This maps the actual sample names to their generated Dorado barcode IDs.
    ###########################################################################

    "$PYTHON" - \
        "$SAMPLES_MANIFEST" \
        "$SAMPLE_MAP" \
        <<'PY'

import csv
import sys

manifest = sys.argv[1]
output = sys.argv[2]

with open(manifest, newline="") as fh:
    rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )

with open(output, "w") as fh:

    fh.write(
        "sample_name\tindex_pair\ti7_F\ti7_REV\ti5_F\ti5_REV\n"
    )

    for number, row in enumerate(rows, start=1):

        sample = row["sample_name"].strip()
        udp = row["index_pair"].strip()

        i7_id = f"i7{number:02d}"
        i5_id = f"i5{number:02d}"

        fh.write(
            f"{sample}\t"
            f"{udp}\t"
            f"{i7_id}\t"
            f"{i7_id}\t"
            f"{i5_id}\t"
            f"{i5_id}\n"
        )

PY


    mark_done "$marker"

    log "Sequencing manifest processed."
    log "Reduced index manifest:"
    log "  ${REDUCED_INDEXES}"
}


###############################################################################
# Generate Dorado TOML
###############################################################################

generate_toml() {

    local name="$1"
    local output="$2"
    local n="$3"
    local i7front="$4"
    local i7rear="$5"

    cat > "$output" <<EOF
[arrangement]
name = "${name}"
kit = "${name}"

mask1_front = "${i7front}"
mask1_rear = "${i7rear}"

barcode1_pattern = "i7%02i"

first_index = 1
last_index = ${n}

rear_only_barcodes = ${REAR_ONLY_BARCODES}

[scoring]
max_barcode_penalty = ${MAX_BARCODE_PENALTY}
barcode_end_proximity = ${BARCODE_END_PROXIMITY}
min_barcode_penalty_dist = ${MIN_BARCODE_PENALTY_DIST}
min_separation_only_dist = ${MIN_SEPARATION_ONLY_DIST}

flank_left_pad = ${FLANK_LEFT_PAD}
flank_right_pad = ${FLANK_RIGHT_PAD}

front_barcode_window = ${FRONT_BARCODE_WINDOW}
rear_barcode_window = ${REAR_BARCODE_WINDOW}

midstrand_flank_score = ${MIDSTRAND_FLANK_SCORE}
EOF
}


generate_toml_i5() {

    local name="$1"
    local output="$2"
    local n="$3"

    cat > "$output" <<EOF
[arrangement]
name = "${name}"
kit = "${name}"

mask1_front = "${I5_FRONT}"
mask1_rear = "${I5_REAR}"

barcode1_pattern = "i5%02i"

first_index = 1
last_index = ${n}

rear_only_barcodes = ${REAR_ONLY_BARCODES}

[scoring]
max_barcode_penalty = ${MAX_BARCODE_PENALTY}
barcode_end_proximity = ${BARCODE_END_PROXIMITY}
min_barcode_penalty_dist = ${MIN_BARCODE_PENALTY_DIST}
min_separation_only_dist = ${MIN_SEPARATION_ONLY_DIST}

flank_left_pad = ${FLANK_LEFT_PAD}
flank_right_pad = ${FLANK_RIGHT_PAD}

front_barcode_window = ${FRONT_BARCODE_WINDOW}
rear_barcode_window = ${REAR_BARCODE_WINDOW}

midstrand_flank_score = ${MIDSTRAND_FLANK_SCORE}
EOF
}


###############################################################################
# Generate FASTAs + TOMLs from reduced manifest
###############################################################################

generate_barcode_files() {

    local manifest="$1"
    local destination_tomls="$2"
    local destination_fastas="$3"
    local marker="$4"

    if stage_done "$marker"; then
        log "Barcode generation already complete. Skipping."
        return
    fi

    log "Generating Dorado barcode files from:"
    log "  $manifest"


    rm -f \
        "$destination_fastas/i7_F.fasta" \
        "$destination_fastas/i7_REV.fasta" \
        "$destination_fastas/i5_F.fasta" \
        "$destination_fastas/i5_REV.fasta" \
        "$destination_tomls/i7_F.toml" \
        "$destination_tomls/i7_REV.toml" \
        "$destination_tomls/i5_F.toml" \
        "$destination_tomls/i5_REV.toml"


    local n

    n="$(
        "$PYTHON" - "$manifest" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="") as fh:
    print(sum(1 for _ in csv.DictReader(fh)))
PY
    )"


    "$PYTHON" - \
        "$manifest" \
        "$destination_fastas" \
        <<'PY'

import csv
import os
import sys

manifest = sys.argv[1]
fastas = sys.argv[2]

os.makedirs(fastas, exist_ok=True)

with open(manifest, newline="") as fh:
    rows = list(csv.DictReader(fh))

if not rows:
    raise SystemExit(
        "Reduced index manifest contains no samples."
    )


with open(os.path.join(fastas, "i7_F.fasta"), "w") as i7f, \
     open(os.path.join(fastas, "i7_REV.fasta"), "w") as i7r, \
     open(os.path.join(fastas, "i5_F.fasta"), "w") as i5f, \
     open(os.path.join(fastas, "i5_REV.fasta"), "w") as i5r:

    for number, row in enumerate(rows, start=1):

        i7_id = f"i7{number:02d}"
        i5_id = f"i5{number:02d}"

        i7_adapter = row[
            "i7 Bases in Adapter"
        ].strip()

        i7_sample = row[
            "i7 Bases for Sample Sheet"
        ].strip()

        i5_forward = row[
            "i5 Bases for Sample Sheet in Forward Orientation"
        ].strip()

        i5_reverse = row[
            "I5 Bases for Sample Sheet in Reverse Complement Orientation"
        ].strip()

        i7f.write(
            f">{i7_id}\n"
            f"{i7_adapter}\n"
        )

        i7r.write(
            f">{i7_id}\n"
            f"{i7_sample}\n"
        )

        i5f.write(
            f">{i5_id}\n"
            f"{i5_forward}\n"
        )

        i5r.write(
            f">{i5_id}\n"
            f"{i5_reverse}\n"
        )


print(
    f"Generated barcode FASTAs for {len(rows)} sample(s)."
)

PY


    generate_toml \
        "Nano_I7_F" \
        "${destination_tomls}/i7_F.toml" \
        "$n" \
        "$I7_FRONT" \
        "$I7_REAR"


    generate_toml \
        "Nano_I7_REV" \
        "${destination_tomls}/i7_REV.toml" \
        "$n" \
        "$I7_FRONT" \
        "$I7_REAR"


    generate_toml_i5 \
        "Nano_I5_F" \
        "${destination_tomls}/i5_F.toml" \
        "$n"


    generate_toml_i5 \
        "Nano_I5_REV" \
        "${destination_tomls}/i5_REV.toml" \
        "$n"


    mark_done "$marker"

    log "Barcode generation complete."
}


###############################################################################
# Run one Dorado pass
###############################################################################

run_dorado_pass() {

    local pass_name="$1"
    local toml="$2"
    local fasta="$3"
    local input_dir="$4"
    local output_dir="$5"
    local marker="$6"

    if stage_done "$marker"; then
        log "Dorado ${pass_name} already complete. Skipping."
        return
    fi


    if [[ "$FORCE" == true ]]; then
        rm -rf "$output_dir"
    fi


    mkdir -p "$output_dir"


    log "Running Dorado ${pass_name}."


    "$DORADO" demux \
        --kit-name "$pass_name" \
        --barcode-arrangement "$toml" \
        --barcode-sequences "$fasta" \
        --emit-fastq \
        --output-dir "$output_dir" \
        "$input_dir" \
        > "${LOGS}/dorado_${pass_name}.log" \
        2>&1


    mark_done "$marker"

    log "Dorado ${pass_name} complete."
}


###############################################################################
# Find unclassified FASTQs
###############################################################################

collect_unclassified() {

    local source_dir="$1"
    local destination_dir="$2"

    rm -rf "$destination_dir"
    mkdir -p "$destination_dir"

    local count=0


    while IFS= read -r -d '' file; do

        local base
        base="$(basename "$file")"

        cp -L "$file" "${destination_dir}/${base}"

        count=$((count + 1))

    done < <(
        find "$source_dir" -type f \
            \( \
                -name '*_unclassified.fastq' \
                -o -name '*_unclassified.fastq.gz' \
                -o -name '*_unclassified.fq' \
                -o -name '*_unclassified.fq.gz' \
            \) \
            -print0
    )


    log "Collected ${count} unclassified FASTQ file(s) from ${source_dir}."
}


###############################################################################
# Prepare Dorado input
###############################################################################

prepare_dorado_input() {

    local source_dir="$1"

    local count

    count="$(
        find "$source_dir" -type f \
            \( \
                -name '*.fastq' \
                -o -name '*.fq' \
                -o -name '*.fastq.gz' \
                -o -name '*.fq.gz' \
            \) |
        wc -l
    )"


    [[ "$count" -gt 0 ]] \
        || die "No FASTQ files available for Dorado input: $source_dir"
}


###############################################################################
# Four-pass production demultiplexing
###############################################################################

run_production_demux() {

    local marker="${DONE}/demux_complete.done"


    if [[ "$SKIP_DEMUX" == true ]]; then
        log "Dorado demultiplexing skipped."
        return
    fi


    if stage_done "$marker"; then
        log "Production demultiplexing already complete. Skipping."
        return
    fi


    log "Starting production four-pass demultiplexing."


    ###########################################################################
    # Stage 1: i7 forward
    #
    # Input: ORIGINAL FASTQs
    ###########################################################################

    local pass1="${DEMUXED}/01_i7_F"


    run_dorado_pass \
        "Nano_I7_F" \
        "${TOMLS}/i7_F.toml" \
        "${FASTAS}/i7_F.fasta" \
        "$PRODUCTION_INPUT" \
        "$pass1" \
        "${DONE}/demux_01_i7_F.done"


    ###########################################################################
    # Stage 2: i7 reverse
    #
    # Input: ONLY i7_F unclassified reads
    ###########################################################################

    local input2="${DEMUXED}/_unclassified_after_i7_F"


    collect_unclassified \
        "$pass1" \
        "$input2"


    prepare_dorado_input "$input2"


    local pass2="${DEMUXED}/02_i7_REV"


    run_dorado_pass \
        "Nano_I7_REV" \
        "${TOMLS}/i7_REV.toml" \
        "${FASTAS}/i7_REV.fasta" \
        "$input2" \
        "$pass2" \
        "${DONE}/demux_02_i7_REV.done"


    ###########################################################################
    # Stage 3: i5 forward
    #
    # Input: ONLY i7_REV unclassified reads
    ###########################################################################

    local input3="${DEMUXED}/_unclassified_after_i7_REV"


    collect_unclassified \
        "$pass2" \
        "$input3"


    prepare_dorado_input "$input3"


    local pass3="${DEMUXED}/03_i5_F"


    run_dorado_pass \
        "Nano_I5_F" \
        "${TOMLS}/i5_F.toml" \
        "${FASTAS}/i5_F.fasta" \
        "$input3" \
        "$pass3" \
        "${DONE}/demux_03_i5_F.done"


    ###########################################################################
    # Stage 4: i5 reverse
    #
    # Input: ONLY i5_F unclassified reads
    ###########################################################################

    local input4="${DEMUXED}/_unclassified_after_i5_F"


    collect_unclassified \
        "$pass3" \
        "$input4"


    prepare_dorado_input "$input4"


    local pass4="${DEMUXED}/04_i5_REV"


    run_dorado_pass \
        "Nano_I5_REV" \
        "${TOMLS}/i5_REV.toml" \
        "${FASTAS}/i5_REV.fasta" \
        "$input4" \
        "$pass4" \
        "${DONE}/demux_04_i5_REV.done"


    mark_done "$marker"

    log "Production four-pass demultiplexing complete."
}


###############################################################################
# Reconstruct samples from four Dorado passes
###############################################################################

reconstruct_samples() {

    local marker="${DONE}/sample_reconstruction.done"

    if stage_done "$marker"; then
        log "Sample reconstruction already complete. Skipping."
        return
    fi

    if [[ "$SKIP_DEMUX" == true ]]; then
        die "Cannot reconstruct samples because demultiplexing was skipped."
    fi

    log "Reconstructing sample FASTQs."

    shopt -s nullglob

    mapfile -t SAMPLE_LINES < <(
        tail -n +2 "$SAMPLE_MAP"
    )

    for line in "${SAMPLE_LINES[@]}"; do

        IFS=$'\t' read -r \
            sample \
            original \
            i7_f \
            i7_rev \
            i5_f \
            i5_rev \
            <<< "$line"

        local output="${SAMPLES}/${sample}.fastq"

        if [[ "$FORCE" == true ]]; then
            rm -f "$output"
        fi

        if [[ "$RESUME" == true && \
              "$FORCE" == false && \
              -s "$output" ]]; then

            log "  ${sample}: already reconstructed. Skipping."
            continue
        fi

        local i7_f_barcode="barcode7${i7_f#i7}"
	local i7_rev_barcode="barcode7${i7_rev#i7}"

	local i5_f_barcode="barcode5${i5_f#i5}"
	local i5_rev_barcode="barcode5${i5_rev#i5}"

        log "  ${sample}: ${original}"
        log "    01_i7_F   -> ${i7_f_barcode}"
        log "    02_i7_REV -> ${i7_rev_barcode}"
        log "    03_i5_F   -> ${i5_f_barcode}"
        log "    04_i5_REV -> ${i5_rev_barcode}"

        local found_files=0

        #######################################################################
        # 01_i7_F
        #######################################################################

        for f in \
            "$DEMUXED/01_i7_F"/*_"${i7_f_barcode}".fastq
        do
            log "    Adding: $(basename "$f")"
            cat "$f" >> "$output"
            found_files=$((found_files + 1))
        done


        #######################################################################
        # 02_i7_REV
        #######################################################################

        for f in \
            "$DEMUXED/02_i7_REV"/*_"${i7_rev_barcode}".fastq
        do
            log "    Adding: $(basename "$f")"
            cat "$f" >> "$output"
            found_files=$((found_files + 1))
        done


        #######################################################################
        # 03_i5_F
        #######################################################################

        for f in \
            "$DEMUXED/03_i5_F"/*_"${i5_f_barcode}".fastq
        do
            log "    Adding: $(basename "$f")"
            cat "$f" >> "$output"
            found_files=$((found_files + 1))
        done


        #######################################################################
        # 04_i5_REV
        #######################################################################

        for f in \
            "$DEMUXED/04_i5_REV"/*_"${i5_rev_barcode}".fastq
        do
            log "    Adding: $(basename "$f")"
            cat "$f" >> "$output"
            found_files=$((found_files + 1))
        done


        #######################################################################
        # Check reconstruction
        #######################################################################

        if [[ "$found_files" -eq 0 ]]; then
            log "  ${sample}: no classified reads"
            : > "$output"
            continue
        fi


        #######################################################################
        # Validate FASTQ structure
        #######################################################################

        local line_count
        line_count=$(wc -l < "$output")

        if (( line_count % 4 != 0 )); then
            die \
                "Invalid reconstructed FASTQ for ${sample}: " \
                "${line_count} lines is not divisible by 4."
        fi

        local read_count=$((line_count / 4))

        log "  ${sample}: ${read_count} reads from ${found_files} FASTQ file(s)"

    done

    shopt -u nullglob

    mark_done "$marker"

    log "Sample reconstruction complete."
}


###############################################################################
# Fastplong
###############################################################################

run_fastplong() {

    if [[ "$SKIP_TRIMMING" == true ]]; then
        log "Fastplong trimming skipped."
        return
    fi


    log "Starting Fastplong trimming."


    mapfile -t SAMPLE_FASTQS < <(
        find "$SAMPLES" \
            -maxdepth 1 \
            -type f \
            -name '*.fastq' |
        sort
    )


    for input in "${SAMPLE_FASTQS[@]}"; do

        local filename
        filename="$(basename "$input" .fastq)"

        local sample_dir="${TRIMMED}/${filename}"

        mkdir -p "$sample_dir"


        local previous="$input"


        for ((pass=1; pass<=FASTPLONG_PASSES; pass++)); do

            local output="${sample_dir}/${filename}.trim${pass}.fastq"

            local marker="${DONE}/fastplong_${filename}_pass${pass}.done"


            if stage_done "$marker"; then

                log "  ${filename}: Fastplong pass ${pass} already complete."

                previous="$output"

                continue

            fi


            log "  ${filename}: Fastplong pass ${pass}/${FASTPLONG_PASSES}"


            "$FASTPLONG" \
                -i "$previous" \
                -o "$output" \
                -a "$ADAPTERS" \
                --distance_threshold "$DISTANCE_THRESHOLD" \
                --trimming_extension "$TRIMMING_EXTENSION" \
                --thread "$THREADS" \
                -l "$MIN_LENGTH" \
                > "${LOGS}/fastplong_${filename}_pass${pass}.log" \
                2>&1


            [[ -s "$output" ]] \
                || die "Fastplong produced empty output: $output"


            mark_done "$marker"

            previous="$output"

        done

    done


    log "Fastplong trimming complete."
}


###############################################################################
# Build minimap2 index
###############################################################################

build_minimap2_index() {

    local marker="${DONE}/minimap2_index.done"


    if [[ "$SKIP_ALIGNMENT" == true ]]; then
        return
    fi


    if stage_done "$marker" && \
       [[ -f "$REFERENCE_MMI" ]]; then

        log "minimap2 reference index already exists. Skipping."

        return

    fi


    if [[ "$FORCE" == true ]]; then
        rm -f "$REFERENCE_MMI"
    fi


    log "Building minimap2 reference index."


    "$MINIMAP2" \
        -d "$REFERENCE_MMI" \
        "$REFERENCE" \
        > "${LOGS}/minimap2_index.log" \
        2>&1


    mark_done "$marker"

    log "minimap2 index complete."
}


###############################################################################
# Alignment
###############################################################################

run_alignment() {

    if [[ "$SKIP_ALIGNMENT" == true ]]; then
        log "Alignment skipped."
        return
    fi


    build_minimap2_index


    log "Starting minimap2 alignment."


    local sample
    local input


    ###########################################################################
    # Select inputs
    ###########################################################################

    if [[ "$SKIP_TRIMMING" == true ]]; then

        mapfile -t ALIGN_INPUTS < <(
            find "$SAMPLES" \
                -maxdepth 1 \
                -type f \
                -name '*.fastq' |
            sort
        )

    else

        mapfile -t ALIGN_INPUTS < <(
            find "$TRIMMED" \
                -type f \
                -name "*.trim${FASTPLONG_PASSES}.fastq" |
            sort
        )

    fi


    ###########################################################################
    # Align each sample
    ###########################################################################

    for input in "${ALIGN_INPUTS[@]}"; do

        if [[ "$SKIP_TRIMMING" == true ]]; then

            sample="$(basename "$input" .fastq)"

        else

            sample="$(
                basename \
                    "$input" \
                    ".trim${FASTPLONG_PASSES}.fastq"
            )"

        fi


        local bam="${ALIGNED}/${sample}.sorted.bam"
        local marker="${DONE}/alignment_${sample}.done"


        if stage_done "$marker" && \
           [[ -s "$bam" ]]; then

            log "  ${sample}: alignment already complete."

            continue

        fi


        log "  ${sample}: minimap2 -ax ${MINIMAP2_PRESET}"


        "$MINIMAP2" \
            -ax "$MINIMAP2_PRESET" \
            -t "$THREADS" \
            "$REFERENCE_MMI" \
            "$input" \
            2> "${LOGS}/minimap2_${sample}.log" |
        "$SAMTOOLS" sort \
            -@ "$THREADS" \
            -o "$bam" \
            -


        [[ -s "$bam" ]] \
            || die "samtools produced empty BAM: $bam"


        "$SAMTOOLS" index \
            -@ "$THREADS" \
            "$bam"


        mark_done "$marker"

    done


    log "Alignment complete."
}


###############################################################################
# Summary
###############################################################################

generate_summary() {

    log "Generating summary."


    {
        printf "sample\treads_raw"

        for ((pass=1; pass<=FASTPLONG_PASSES; pass++)); do
            printf "\treads_trim%d" "$pass"
        done

        printf "\treads_aligned\n"


        mapfile -t samples < <(
            find "$SAMPLES" \
                -maxdepth 1 \
                -type f \
                -name '*.fastq' |
            sort
        )


        for raw in "${samples[@]}"; do

            local sample
            sample="$(basename "$raw" .fastq)"


            local raw_reads

            raw_reads="$(
                awk 'END {print NR/4}' "$raw"
            )"


            printf "%s\t%s" \
                "$sample" \
                "$raw_reads"


            if [[ "$SKIP_TRIMMING" == false ]]; then

                for ((pass=1; pass<=FASTPLONG_PASSES; pass++)); do

                    local trim="${TRIMMED}/${sample}/${sample}.trim${pass}.fastq"


                    if [[ -f "$trim" ]]; then

                        local n

                        n="$(
                            awk 'END {print NR/4}' "$trim"
                        )"

                        printf "\t%s" "$n"

                    else

                        printf "\tNA"

                    fi

                done

            else

                for ((pass=1; pass<=FASTPLONG_PASSES; pass++)); do
                    printf "\tNA"
                done

            fi


            if [[ "$SKIP_ALIGNMENT" == false ]]; then

                local bam="${ALIGNED}/${sample}.sorted.bam"


                if [[ -f "$bam" ]]; then

                    local aligned

                    aligned="$(
                        "$SAMTOOLS" view -c -F 4 "$bam"
                    )"

                    printf "\t%s" "$aligned"

                else

                    printf "\tNA"

                fi

            else

                printf "\tNA"

            fi


            printf "\n"

        done

    } > "$SUMMARY"


    log "Summary:"
    log "  ${SUMMARY}"
}


###############################################################################
# Main pipeline
###############################################################################

log ""
log "============================================================"
log "Nano Demux Pipeline"
log "============================================================"
log ""
log "Input:"
log "  ${INPUT}"
log ""
log "Master index manifest:"
log "  ${INDEXES}"
log ""
log "Sequencing manifest:"
log "  ${SAMPLES_MANIFEST}"
log ""
log "Reference:"
log "  ${REFERENCE}"
log ""
log "Output:"
log "  ${OUTPUT}"
log ""


###############################################################################
# Validate master index CSV
###############################################################################

log "Validating master index manifest:"
validate_index_manifest "$INDEXES"


###############################################################################
# Reduce master index CSV using sequencing manifest
###############################################################################

prepare_sequencing_manifest


###############################################################################
# Production index manifest
###############################################################################

PRODUCTION_MANIFEST="$REDUCED_INDEXES"


###############################################################################
# Generate production barcode files
###############################################################################

generate_barcode_files \
    "$PRODUCTION_MANIFEST" \
    "$TOMLS" \
    "$FASTAS" \
    "${DONE}/production_barcode_generation.done"


###############################################################################
# Prepare production input directory
###############################################################################

PRODUCTION_INPUT="${OUTPUT}/input"


if [[ "$FORCE" == true || ! -d "$PRODUCTION_INPUT" ]]; then

    rm -rf "$PRODUCTION_INPUT"

    mkdir -p "$PRODUCTION_INPUT"


    for fq in "${FASTQS[@]}"; do

        ln -sf \
            "$fq" \
            "${PRODUCTION_INPUT}/$(basename "$fq")"

    done

fi


###############################################################################
# Production demultiplexing
###############################################################################

run_production_demux


###############################################################################
# Reconstruct samples
###############################################################################

reconstruct_samples


###############################################################################
# Fastplong
###############################################################################

run_fastplong


###############################################################################
# Alignment
###############################################################################

run_alignment


###############################################################################
# Summary
###############################################################################

generate_summary


###############################################################################
# Optional cleanup
###############################################################################

if [[ "$KEEP_INTERMEDIATE" == false ]]; then

    log "Removing intermediate files."

    rm -rf "$PRODUCTION_INPUT"

fi


###############################################################################
# Final report
###############################################################################

log ""
log "============================================================"
log "Nano Demux PIPELINE COMPLETE"
log "============================================================"
log ""
log "Output:"
log "  ${OUTPUT}"
log ""
log "Sequencing manifest:"
log "  ${SAMPLES_MANIFEST}"
log ""
log "Reduced index manifest:"
log "  ${REDUCED_INDEXES}"
log ""
log "Sample map:"
log "  ${SAMPLE_MAP}"
log ""
log "Generated barcode files:"
log "  ${FASTAS}"
log "  ${TOMLS}"
log ""
log "Sample FASTQs:"
log "  ${SAMPLES}"
log ""
log "Trimmed FASTQs:"
log "  ${TRIMMED}"
log ""
log "Aligned BAMs:"
log "  ${ALIGNED}"
log ""
log "Summary:"
log "  ${SUMMARY}"
log ""

log "Pipeline finished successfully."
