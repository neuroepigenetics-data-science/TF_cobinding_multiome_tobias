#!/usr/bin/env bash
# ============================================================================
# sync_to_share.sh
# ----------------------------------------------------------------------------
# Copy this project to the neuroepigenetik SMB share so it can be picked up on
# another machine, and stage Claude Code's memory alongside it.
#
# MEASURED FACTS ABOUT THE SHARE (2026-09-03), which drive every choice below:
#
#   mount        //gdurmaz@138.245.214.48/neuroepigenetik -> /Volumes/neuroepigenetik
#   free space   395 Ti   (project is 134 G -- space is not a constraint)
#   throughput   ~7.5 MB/s measured on a 200 MB write
#   project      134 G in 43,307 files across 17,913 directories
#
#   => THE FULL COPY IS AN OVERNIGHT JOB, roughly 5 h at streaming speed and
#      realistically longer, because SMB pays a round trip per file and 43k
#      files is a lot of round trips. Hence the two phases: ESSENTIALS first,
#      so the project is usable on the other machine within a minute, then
#      BULK, which can run unattended and be interrupted freely.
#
# WHY THESE rsync FLAGS. macOS does not ship real rsync; /usr/bin/rsync is
# openrsync advertising "2.6.9 compatible". Tested against this share:
#     -rt --partial --no-perms --exclude=PATTERN     all work
#     --no-owner, --no-group, --modify-window 2      REJECTED
#     --exclude .DS_Store (space form)               REJECTED
# openrsync needs --opt=value; the space-separated form is a parse error.
# Perms are dropped because SMB cannot represent POSIX ownership anyway.
# `brew install rsync` would give rsync 3.x and --info=progress2, but is not
# required -- --progress works here.
#
# RESUMABLE. Re-running skips files already present with matching size+mtime,
# so an interrupted bulk phase is restarted by just running it again.
# NOTHING IS EVER DELETED AT THE DESTINATION -- --delete is deliberately unused.
#
# *** DO NOT RUN THE PIPELINE OFF THIS SHARE. *** At 7.5 MB/s, ATACorrect --
# which streams whole BAMs repeatedly during bias estimation -- would take days
# instead of hours. Treat the share as transport and archive; copy to local
# disk on the far machine before computing.
#
# USAGE
#   pipeline/sync_to_share.sh essentials        # dry run, shows what would move
#   pipeline/sync_to_share.sh essentials --go   # ~155 MB, under a minute
#   pipeline/sync_to_share.sh tfcomb --go       # ~11 G, what Tier 3 needs
#   pipeline/sync_to_share.sh bulk --go         # ~133 G, overnight
#   pipeline/sync_to_share.sh all --go
# ============================================================================
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
DEST="${DEST:-/Volumes/neuroepigenetik/TF_cobinding_multiome_tobias}"
MEMSRC="$HOME/.claude/projects/-Users-gorkemdurmaz-projects-TF-cobinding-multiome-TF-cobinding-multiome-tobias"
MEMDEST="$(dirname "$DEST")/claude_memory_TF_cobinding_multiome_tobias"

PHASE="${1:-}"
GO="${2:-}"
[[ "$PHASE" =~ ^(essentials|tfcomb|bulk|all)$ ]] || { sed -n '/^# USAGE/,/^# ===/p' "$0"; exit 1; }

RS=(rsync -rt --partial --no-perms --progress
    --exclude=.DS_Store --exclude=__pycache__ --exclude=.ipynb_checkpoints)
[[ "$GO" == "--go" ]] || { RS+=(--dry-run); echo ">>> DRY RUN -- add --go to actually copy"; }

# The bulk directories. Everything else counts as essentials.
BULK=(bam tobias qc logs)

mountpoint_check() {
  [[ -d /Volumes/neuroepigenetik ]] || { echo "ERROR: share not mounted at /Volumes/neuroepigenetik"; exit 1; }
}

sync_essentials() {
  echo "=== ESSENTIALS: code, notes, criteria, figures, meta (~155 MB) ==="
  local ex=(); for d in "${BULK[@]}"; do ex+=("--exclude=$d"); done
  "${RS[@]}" "${ex[@]}" "$SRC/" "$DEST/"

  echo
  echo "=== CLAUDE MEMORY -> $MEMDEST ==="
  # Staged OUTSIDE the project tree so nothing pollutes the git repo.
  # install_claude_memory.sh on the far machine puts it in the right place.
  "${RS[@]}" "$MEMSRC/" "$MEMDEST/"
  [[ "$GO" == "--go" ]] && cp "$SRC/pipeline/install_claude_memory.sh" "$MEMDEST/" || true
}

sync_tfcomb() {
  # ------------------------------------------------------------------------
  # The subset needed to run TF-COMB (Tier 3) and to remake the footprint
  # figures -- ~11 G of the 133 G, chosen by what each consumer actually reads.
  #
  # INCLUDED
  #   *_bound.bed        10,444 files, 5.0 G. TF-COMB's market-basket input.
  #                      ALL 746 motifs, not just Tier 3's four pairs: the
  #                      co-occurrence of a pair is ranked against the full
  #                      background, so a partial motif set changes the answer.
  #   *_results.txt      per-group BINDetect summary, 30 M with the below
  #   *_distances.txt    the 746x746 motif-similarity matrix -- the file that
  #                      corrected three claims on 2026-08-26 (SOX9/SOX10 are
  #                      separable at 0.999, voiding Tier 3 limitation 2)
  #   footprints/*.bw    6.1 G in 16 files. PlotAggregate input for the
  #                      Extended Data 6b analogue.
  #
  # DELIBERATELY EXCLUDED, and why -- these are not oversights:
  #   *_overview.txt     18 G. Per-site scores; nothing downstream reads them.
  #   *_all.bed          16 G. Every motif occurrence bound or not; PlotAggregate
  #                      does use _all.bed, so RE-ADD THIS if the aggregate
  #                      figures need remaking from scratch rather than viewing.
  #   *_unbound.bed      22 G. The complement of _bound; never an input.
  #   bam/               35 G. Only needed to re-run ATACorrect.
  #   atacorrect/        26 G. Only needed to re-run ScoreBigwig.
  #   bindetect_dars/     5 G. The DAR secondary analysis, concluded 08-20.
  #
  # WHY A FILE LIST AND NOT --include/--exclude: tested 2026-09-03, openrsync
  # IGNORES include/exclude ordering and copies everything anyway (4 files
  # where 2 were asked for). --files-from is exact. Do not "simplify" this.
  # ------------------------------------------------------------------------
  echo "=== TF-COMB ESSENTIALS: bound BEDs + tables + footprint bigwigs (~11 G) ==="
  local list; list="$(mktemp)"
  cd "$SRC"
  {
    find tobias/bindetect -name '*_bound.bed'
    find tobias/bindetect -maxdepth 2 \( -name '*_results.txt' -o -name '*_distances.txt' \
                                       -o -name '*_results.xlsx' -o -name 'bindetect_*.html' \)
    find tobias/footprints -name '*.bw'
  } | sed 's|^\./||' | sort -u > "$list"
  echo "    $(wc -l < "$list" | tr -d ' ') files listed"
  "${RS[@]}" --files-from="$list" "$SRC/" "$DEST/"
  rm -f "$list"
}

sync_bulk() {
  echo "=== BULK: ${BULK[*]} (~133 G, hours) ==="
  for d in "${BULK[@]}"; do
    [[ -d "$SRC/$d" ]] || continue
    echo "--- $d ($(du -sh "$SRC/$d" | cut -f1)) ---"
    "${RS[@]}" "$SRC/$d/" "$DEST/$d/"
  done
}

mountpoint_check
case "$PHASE" in
  essentials) sync_essentials ;;
  tfcomb)     sync_essentials; sync_tfcomb ;;
  bulk)       sync_bulk ;;
  all)        sync_essentials; sync_bulk ;;
esac

echo
echo "Done. Verify counts with:"
echo "  find '$SRC' -type f | wc -l ; find '$DEST' -type f | wc -l"
