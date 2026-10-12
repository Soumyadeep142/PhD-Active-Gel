#!/bin/bash

# 1. Compile and run Fortran
gfortran -I/usr/include 019_Hernandez_c_everywhere.f -lfftw3 -lm
./a.out

# 2. Find the latest params folder
folder=$(ls -td 019_params/* | head -1)

echo "Using folder: $folder"

# 3. Copy processing scripts into folder
cp 016_converter.py "$folder"/
cp data_sheet_trimmer.sh "$folder"/
cp 016b_cutterx_converter.py "$folder"/
cp 016b_cuttery_converter.py "$folder"/
cp 016b_cuttery_converter.py "$folder"/
cp 019_gnu "$folder"/
# 4. Enter folder
cd "$folder" || exit

# 5. Run processing
python3 016_converter.py
bash data_sheet_trimmer.sh
python3 016b_cutterx_converter.py
python3 016b_cuttery_converter.py

# 6. Come back out
cd - || exit

# 7. Continue remaining steps
python3 Push.py
python3 Completion.py
