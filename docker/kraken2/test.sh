#!/bin/bash
set -x  # Enable command echo
exec 1> >(tee -a /dev/stdout)  # Force output to stdout

echo "####################"
echo "# Starting tests   #"
echo "####################"

echo -e "\n# Testing Kraken2:"
echo "Kraken2 version:"
kraken2 --version

echo -e "\n# Testing KrakenTools:"
if command -v kreport2krona.py > /dev/null 2>&1; then
    echo "✓ kreport2krona.py exists"
    echo "Location: $(which kreport2krona.py)"
else 
    echo "✗ kreport2krona.py not found"
    exit 1
fi

echo -e "\n# Testing KronaTools:"
if command -v ktImportText > /dev/null 2>&1; then
    echo "✓ ktImportText exists"
    echo "Location: $(which ktImportText)"
else 
    echo "✗ ktImportText not found"
    exit 1
fi

echo -e "\n####################"
echo "# Tests complete   #"
echo "####################"
exit 0