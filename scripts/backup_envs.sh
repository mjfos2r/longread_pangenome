#!/bin/bash

# timestamp output directory for the envs.
BACKUP_DIR="conda_backups_$(date +%Y%m%d_%H%M)"
mkdir -p "$BACKUP_DIR"

# Save list of all present environments
conda env list > "$BACKUP_DIR/env_list.txt"

# Export each environment
while read -r line; do
    # Skip empty lines and comments
    [[ -z "$line" ]] && continue
    [[ $line =~ ^# ]] && continue

    # Extract environment name (first word of the line)
    env_name=$(echo "$line" | awk '{print $1}')

    # Skip if empty
    [[ -z "$env_name" ]] && continue

    echo "Backing up environment: $env_name"

    # Export environment YAML (WITHOUT BUILDS SO YOU CAN IMPORT QUICK N EZ)
    conda env export -n "$env_name" --no-builds > "$BACKUP_DIR/${env_name}_environment.yml"

    # Export explicit package list
    conda list -n "$env_name" > "$BACKUP_DIR/${env_name}_package_list.txt"

    # Export pip packages
    source activate "$env_name" && pip freeze > "$BACKUP_DIR/${env_name}_pip_packages.txt"

done < "$BACKUP_DIR/env_list.txt"

echo "Backup completed in $BACKUP_DIR"