#!/usr/bin/env bash
#
# Usage:
#   ./parallel_tailboxbg.sh <ALIGN_DIR> <GENBANK_DIR> <SCRIPT_PY> <OUTFILE>
#
# 1) We create (or truncate) 'progress.log' and display it with dialog --tailboxbg.
# 2) We run a single find -> parallel pipeline, writing output rows to $OUTFILE
#    and sending progress messages to 'progress.log'.
# 3) Tailbox updates in real time.
#

###############################################################################
# 0) Parse arguments
###############################################################################
ALIGN_DIR="$1"    # Directory containing alignment .tsv files
GENBANK_DIR="$2"  # Directory with .gbff
SCRIPT_PY="$3"    # Python coverage script
OUTFILE="$4"      # Single combined TSV result

if [[ -z "$ALIGN_DIR" || -z "$GENBANK_DIR" || -z "$SCRIPT_PY" || -z "$OUTFILE" ]]; then
  echo "Usage: $0 <ALIGN_DIR> <GENBANK_DIR> <SCRIPT_PY> <OUTFILE>"
  exit 1
fi

# Make sure the directory for OUTFILE exists
mkdir -p "$(dirname "$OUTFILE")"


N=$( find $GENBANK_DIR -type f -name "*.gbff" | wc -l )
THEO_MAX=$(( N*(N-1)/2 ))
###############################################################################
# 1) Prepare a log file to display with dialog --tailboxbg
###############################################################################
LOGFILE="/tmp/progress.log"
mkfifo "$LOGFILE"  # Clear or create

# Start tailbox in background to display $LOGFILE in real-time
dialog --keep-window --begin 10 70 --gauge 'Expected Files: %s'<${THEO_MAX} 7 70 --and-widget --title "Coverage Progress" --keep-window --begin 18 55 --tailboxbg $LOGFILE 18 155 &
DIALOG_PID=$!

###############################################################################
# 2) Write just one header row to $OUTFILE
###############################################################################
cat <<EOF > "$OUTFILE"
QUERY_ID	QUERY_LEN	QUERY_COORDS	QUERY_COV(%)	QUERY_NUM_GENES	REF_ID	REF_LEN	REF_COORDS	REF_COV(%)	REF_NUM_GENES	IDENTITY(%)
EOF

# Let the user see that in the log right away
echo "[INFO] Wrote header to $OUTFILE" > "$LOGFILE"

###############################################################################
# 3) Single-pass pipeline: find -> parallel
###############################################################################
#    - We do NOT run 'find' twice.
#    - For each file found, we run Python coverage script appending to $OUTFILE.
#    - We also tee all stderr lines into the log so you can watch them live.
#
#    If your Python script prints progress to stderr, you'll see it in the dialog.
###############################################################################
echo "[INFO] Starting alignment processing at $(date)" > "$LOGFILE"

# We can do one big group so we capture ALL output in the log
{
  # Each line from 'find' is a .tsv alignment file
  find "$ALIGN_DIR" -type f -name "*.tsv" | pv -n -s $THEO_MAX | parallel --jobs 4 --no-notice '
    echo "[INFO] Processing file: {}" >&2
    python "'"$SCRIPT_PY"'" \
      --annotations_dir "'"$GENBANK_DIR"'" \
      --alignment_file {} \
      --no_header \
    >> "'"$OUTFILE"'"
    echo "[INFO] Done with file: {}" >&2
  '
} 2>$LOGFILE

echo "[INFO] Processing finished at $(date)" > "$LOGFILE"

###############################################################################
# 4) Cleanup
###############################################################################
# Kill the background dialog
sleep 1
kill $DIALOG_PID 2>/dev/null || true

dialog --msgbox "All alignments processed!\n\nResults in:\n$OUTFILE" 8 50
clear
