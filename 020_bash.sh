#!/bin/bash

# ===== AUTO RUN NAME =====
base_dir="020a_params"
mkdir -p "$base_dir"

# find next available A number
last=$(ls "$base_dir" 2>/dev/null | grep -E '^A[0-9]+' | sed 's/A//' | sort -n | tail -1)

if [ -z "$last" ]; then
    next=1
else
    next=$((last + 1))
fi

run_name="A${next}"

if [ -z "$run_name" ]; then
    echo "Usage: ./run.sh A1"
    exit 1
fi

# ===== FOLDER SETUP =====
base_dir="020a_params"
folder="$base_dir/$run_name"

mkdir -p "$folder"

echo "Running case: $run_name"
echo "Using folder: $folder"

# ===== COPY FILES =====
cp params.inc "$folder"/
cp 020_Hernandez_complete.f "$folder"/

# copy processing scripts
cp 016_converter.py "$folder"/
cp data_sheet_trimmer.sh "$folder"/
cp 016b_cutterx_converter.py "$folder"/
cp 016b_cuttery_converter.py "$folder"/
cp 019_gnu "$folder"/
cp Push.py "$folder"/
cp Completion.py "$folder"/

# ===== ENTER FOLDER =====
cd "$folder" || exit

# ===== COMPILE & RUN =====
gfortran -I/usr/include 020_Hernandez_complete.f -lfftw3 -lm
./a.out "$run_name"

# ===== POST-PROCESS =====
python3 016_converter.py
bash data_sheet_trimmer.sh
python3 016b_cutterx_converter.py
python3 016b_cuttery_converter.py

python3 Push.py
python3 Completion.py

# ===== EXIT =====
cd - || exit

echo "Run completed for $run_name"
