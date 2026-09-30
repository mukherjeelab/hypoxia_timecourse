#!/usr/bin/env bash
# RNAplfold accessibility for the MCF7-SIX1 dominant-transcript set, run LOCALLY in parallel.
#
#   bash cluster_scripts/viennarna/02d_run_rnaplfold_mcf7six1.sh [n_parallel]
#
# Parameters are matched PER REGION to the MDA-MB-231 run, so the accessibility values are
# comparable across cell lines. They are not the same for all three regions - see below.
#
# THE 5'UTR USES DIFFERENT WINDOW PARAMETERS FROM THE CDS AND 3'UTR.
# Recovered by regenerating stored MDA-MB-231 files and requiring a byte-identical match
# (02e_validate_rnaplfold_params.R), over 60 transcripts per region:
#
#   region   parameters              stored files reproduced
#   5'UTR    -W  80 -L  40 -u 30     yes
#   CDS      -W 150 -L 100 -u 30     60/60
#   3'UTR    -W 150 -L 100 -u 30     60/60
#
# Using -W 150 -L 100 on the 5'UTR reproduced only 12 of 60 - and those 12 are short 5'UTRs
# where the window exceeds the sequence, so -W has no effect. On the rest the values differed by
# up to 0.94 on a 0-1 scale, i.e. entirely different numbers. The 5'UTR run predates the
# CDS/3'UTR run and its parameters were never written down; CLAUDE.md documents "-W 150 -L 100"
# for the struct_* features generally, which is right for two regions out of three.
#
# CONSEQUENCE BEYOND THIS SCRIPT: tco_struct_accessibility_utr5_* is not on the same footing as
# the CDS and 3'UTR columns. A smaller window means less long-range pairing is considered, so
# 5'UTR accessibility is systematically higher than it would be at -W 150. Cross-region
# comparisons of those features are comparing two different calculations.
#
# -u 30 IS REQUIRED AND THE EARLIER SCRIPTS OMIT IT.
# 02b_run_rnaplfold_cds_utr3.bsub and 02c_run_rnaplfold_teleman.sh both call
# "RNAplfold -W 150 -L 100" with no -u, then count *_lunp files. On ViennaRNA 2.7.2 that
# combination writes only _dp.ps and produces NO _lunp at all, so those scripts as committed
# cannot have generated the stored output. The stored files' own header runs "l=1" through 30,
# which is what -u 30 means, so 30 is read back off the data rather than guessed. Regenerating
# two stored CDS files with these flags reproduced them BYTE-IDENTICALLY, which is the check
# that licenses generating the rest.
#
# Output goes into the SAME directories as the 231 run. An _lunp file is a property of the
# sequence, so a transcript dominant in both cell lines needs one copy; 01d_ already excluded
# those from the FASTAs. Nothing is ever overwritten - the move step asserts that.

set -euo pipefail
NPAR="${1:-6}"
CELL="${2:-mcf7six1}"   # mcf7six1 (default) | hela | dhx29
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
FA_DIR="${REPO}/accessories/plfold_output/${CELL}_fasta"
PLF="${REPO}/accessories/plfold_output/rnaplfold_output"
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/plf6_XXXXXX")"
trap 'rm -rf "${SCRATCH}"' EXIT

# macOS ships bash 3.2, which has no associative arrays: "declare -A" there makes [utr5]=... an
# arithmetic expression and dies under set -u with "utr5: unbound variable". A function keeps
# this runnable with /bin/bash on this machine as well as with a newer bash on the cluster.
region_dir() {
  case "$1" in
    utr5) echo rnaplfold_output ;;
    cds)  echo cds_rnaplfold_output ;;
    utr3) echo utr3_rnaplfold_output ;;
    *)    echo "unknown region: $1" >&2; exit 1 ;;
  esac
}

# Per-region parameters. Not a single global setting - see the note above.
region_flags() {
  case "$1" in
    utr5) echo "-W 80 -L 40 -u 30" ;;
    cds)  echo "-W 150 -L 100 -u 30" ;;
    utr3) echo "-W 150 -L 100 -u 30" ;;
    *)    echo "unknown region: $1" >&2; exit 1 ;;
  esac
}

for REGION in utr5 cds utr3; do
  FA="${FA_DIR}/${CELL}_${REGION}_sequences.fa"
  RDIR="$(region_dir "${REGION}")"
  FLAGS="$(region_flags "${REGION}")"
  DEST="${PLF}/${RDIR}"
  [[ -s "${FA}" ]] || { echo "skip ${REGION}: ${FA} missing or empty"; continue; }
  NSEQ=$(grep -c '^>' "${FA}" || true)
  echo "[$(date +%H:%M:%S)] ${REGION}: ${NSEQ} sequences, ${NPAR} workers, ${FLAGS} -> ${RDIR}"

  # Split by record, round-robin, so every worker gets a mix of short and long sequences
  # rather than one worker inheriting all the long 3'UTRs.
  awk -v n="${NPAR}" -v d="${SCRATCH}" -v r="${REGION}" '
    /^>/ { i = (i % n) + 1 } { print >> (d "/" r "_part" i ".fa") }' "${FA}"

  for p in $(seq 1 "${NPAR}"); do
    PART="${SCRATCH}/${REGION}_part${p}.fa"
    [[ -s "${PART}" ]] || continue
    W="${SCRATCH}/${REGION}_w${p}"; mkdir -p "${W}"
    # shellcheck disable=SC2086  -- FLAGS must word-split into separate arguments
    ( cd "${W}" && RNAplfold ${FLAGS} < "${PART}" > /dev/null 2>&1 ) &
  done
  wait

  moved=0; clashed=0
  for p in $(seq 1 "${NPAR}"); do
    W="${SCRATCH}/${REGION}_w${p}"
    [[ -d "${W}" ]] || continue
    for f in "${W}"/*_lunp; do
      [[ -e "${f}" ]] || continue
      b="$(basename "${f}")"
      if [[ -e "${DEST}/${b}" ]]; then clashed=$((clashed+1)); continue; fi
      mv "${f}" "${DEST}/${b}"; moved=$((moved+1))
    done
    rm -rf "${W}"
  done
  echo "[$(date +%H:%M:%S)] ${REGION}: moved ${moved} new _lunp; ${clashed} already existed (left untouched)"
  if [[ "${clashed}" -gt 0 ]]; then
    echo "  NOTE: a pre-existing file was never overwritten. 01d_ should have excluded these," >&2
    echo "        so a nonzero count means the FASTA and the directory disagree - check 01d_." >&2
  fi
  # find, not a glob: these directories pass 15,000 files and "ls dir/*_lunp" then exceeds
  # ARG_MAX and silently reports 0, which reads exactly like data loss.
  echo "[$(date +%H:%M:%S)] ${REGION}: ${DEST} now holds $(find "${DEST}" -name '*_lunp' | wc -l | tr -d ' ') files"
done
echo "All RNAplfold regions complete."
