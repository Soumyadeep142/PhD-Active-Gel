#!/bin/bash

# ===== AUTO RUN NAME =====
base_dir="021_params"
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
base_dir="021_params"
folder="$base_dir/$run_name"

mkdir -p "$folder"

echo "Running case: $run_name"
echo "Using folder: $folder"

# ===== COPY FILES =====
cp params.inc "$folder"/
cp 021_Hernandez_rotated_ellipsoid.f "$folder"/

# copy processing scripts
cp 016_converter.py "$folder"/
cp 21_data_sheet_trimmer.sh "$folder"/
cp 021_cutter_ellipse_converter.py "$folder"/
cp 021_cutter_circle_converter.py "$folder"/
cp Push.py "$folder"/
cp Completion.py "$folder"/

# ===== ENTER FOLDER =====
cd "$folder" || exit

# ===== COMPILE & RUN =====
gfortran -I/usr/include 021_Hernandez_rotated_ellipsoid.f -lfftw3 -lm
./a.out "$run_name"

# ===== POST-PROCESS =====
python3 016_converter.py
bash 21_data_sheet_trimmer.sh
python3 021_cutter_ellipse_converter.py
python3 021_cutter_circle_converter.py

python3 Push.py
python3 Completion.py

# ===== EXIT =====
cd - || exit

echo "Run completed for $run_name"
