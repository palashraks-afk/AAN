#!/bin/sh
# Make ~/AAN_mac: a copy of the repo with macOS paths, so ~/AAN stays clean for commits (Windows paths).
set -e
mkdir -p "$HOME/AAN_mac"
rsync -a --delete --exclude .git "$HOME/AAN/" "$HOME/AAN_mac/"
sed -i '' -e "s#D:/AAN_data#$HOME/AAN_data#g" -e "s#C:/AAN_ref#$HOME/AAN_ref#g" -e "s#D:/AAN#$HOME/AAN_mac#g" -e 's#magma\.exe#magma#g' "$HOME"/AAN_mac/scripts/*.py
echo synced
