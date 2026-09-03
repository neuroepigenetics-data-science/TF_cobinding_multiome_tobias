#!/usr/bin/env bash
# ============================================================================
# install_claude_memory.sh
# ----------------------------------------------------------------------------
# Put Claude Code's project memory in the right place on a NEW machine.
# Run this ON THE NEW MACHINE, from inside the project directory, BEFORE
# starting Claude Code there.
#
# THE PROBLEM THIS SOLVES. Claude Code does not keep memory in the project.
# It keeps it in ~/.claude/projects/<slug>/, where <slug> is the project's
# ABSOLUTE PATH with '/', '_' and '.' each replaced by '-'. So:
#
#   /Users/gorkemdurmaz/projects/TF_cobinding_multiome/TF_cobinding_multiome_tobias
#     -> -Users-gorkemdurmaz-projects-TF-cobinding-multiome-TF-cobinding-multiome-tobias
#
#   /Volumes/neuroepigenetik/TF_cobinding_multiome_tobias
#     -> -Volumes-neuroepigenetik-TF-cobinding-multiome-tobias
#
# BECAUSE THE SLUG ENCODES THE PATH, memory copied under the old machine's
# slug is invisible if the project lives anywhere else. This script derives
# the slug from the CURRENT directory, so it is correct wherever you put the
# project -- same-path Mac, different-path Mac, or a Linux box.
#
# WHAT IT INSTALLS
#   memory/*.md   the 21 memory files -- the project's accumulated context
#   MEMORY.md     the index Claude loads every session; without it the
#                 individual memory files are never surfaced
#   *.jsonl       past session transcripts (optional, --with-transcripts),
#                 which is what /resume reads
#
# USAGE
#   cd /path/to/TF_cobinding_multiome_tobias
#   /path/to/claude_memory_bundle/install_claude_memory.sh [BUNDLE_DIR] [--with-transcripts]
#
# BUNDLE_DIR defaults to this script's own directory, so running it straight
# out of the staged bundle needs no arguments.
# ============================================================================
set -euo pipefail

BUNDLE="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
[[ "$BUNDLE" == --* ]] && BUNDLE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WANT_TRANSCRIPTS=false
for a in "$@"; do [[ "$a" == "--with-transcripts" ]] && WANT_TRANSCRIPTS=true; done

PROJECT_DIR="$(pwd -P)"
SLUG="$(printf '%s' "$PROJECT_DIR" | sed 's/[/_.]/-/g')"
TARGET="$HOME/.claude/projects/$SLUG"

echo "project : $PROJECT_DIR"
echo "slug    : $SLUG"
echo "target  : $TARGET"
echo "bundle  : $BUNDLE"
echo

# Sanity: refuse to install memory against a directory that is not the project.
if [[ ! -f "$PROJECT_DIR/VALIDATION_CRITERIA.md" ]]; then
  echo "ERROR: no VALIDATION_CRITERIA.md here -- is this really the project root?"
  echo "       cd into the project directory and re-run."
  exit 1
fi
[[ -d "$BUNDLE/memory" ]] || { echo "ERROR: no memory/ in bundle: $BUNDLE"; exit 1; }

mkdir -p "$TARGET/memory"
cp "$BUNDLE"/memory/*.md "$TARGET/memory/"
cp "$BUNDLE"/MEMORY.md   "$TARGET/"
echo "installed $(ls "$TARGET"/memory/*.md | wc -l | tr -d ' ') memory files + MEMORY.md"

if $WANT_TRANSCRIPTS; then
  # Transcripts are named by session UUID and are path-independent, so they
  # can be copied across unchanged.
  shopt -s nullglob
  t=("$BUNDLE"/*.jsonl)
  if (( ${#t[@]} )); then cp "${t[@]}" "$TARGET/"; echo "installed ${#t[@]} session transcripts"; fi
fi

echo
echo "Done. Start Claude Code from $PROJECT_DIR and ask it to read its memory."
