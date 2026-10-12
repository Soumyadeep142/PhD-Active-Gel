#!/bin/bash

# ===== BASE DIRECTORY =====
base_dir="024_params"
mkdir -p "$base_dir"

MAX_JOBS=$(nproc)

# ===== FIND NEXT AVAILAsLE A NUMBER =====
last=$(ls "$base_dir" 2>/dev/null \
| grep -E '^A[0-9]+' \
| sed -E 's/^A([0-9]+).*/\1/' \
| sort -n \
| tail -1)

if [ -z "$last" ]; then
    next=1
else
    next=$((last + 1))
fi

run_name="A${next}"

echo "Base run name: $run_name"

# ============================================================
# PARALLEL LOOP : a_x = 10 -> 21
# ============================================================

for sig_VALUE in $(seq 0.0 0.02 0.3)
do
while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]
do
    wait -n
done
(
    # format sigma safely for folder names
    sig_TAG=$(printf "%.2f" $sig_VALUE | sed 's/\./p/g')

    local_run="${run_name}_sig_${sig_TAG}"

    folder="$base_dir/$local_run"

    mkdir -p "$folder"

    # ===== COPY FILES =====
    cp params.inc "$folder"/
    cp 024_Hernandez_updated_surface_tension.f "$folder"/

    cp 016_converter.py "$folder"/
    cp 21_data_sheet_trimmer.sh "$folder"/
    cp 021_cutter_ellipse_converter.py "$folder"/
    cp 021_cutter_circle_converter.py "$folder"/
    cp Push.py "$folder"/
    cp Completion.py "$folder"/

    # ===== ENTER FOLDER =====
    cd "$folder" || exit

    # ===== MODIFY a_x =====
#    sed -i "s/\(a_x[[:space:]]*=[[:space:]]*\).*/\1${AX_VALUE}.d0/" params.inc
# ===== MODIFY sig =====
    sed -i "s/\(sig[[:space:]]*=[[:space:]]*\).*/\1${sig_VALUE}d0/" params.inc
    # ===== COMPILE =====
    gfortran -I/usr/include 024_Hernandez_updated_surface_tension.f -lfftw3 -lm

    # ===== RUN =====
    ./a.out "$local_run"

    # ===== POST PROCESS =====
    python3 016_converter.py

    bash 21_data_sheet_trimmer.sh

    python3 021_cutter_ellipse_converter.py

    python3 021_cutter_circle_converter.py


) &
sleep 0.2

done

# ===== WAIT FOR ALL PARALLEL JOBS =====
wait

python3 Push.py
python3 Completion.py
