gfortran -I/usr/include 016_Hernandez_conc_curve_phi.f -lfftw3 -lm
./a.out 
python3 016_converter.py
bash data_sheet_trimmer.sh
python3 016b_cutterx_converter.py
python3 016b_cuttery_converter.py
python3 Completion.py
