
mkdir -p ./build;
cd ./build && cmake .. ; make ; ./main &
sleep(2)
mv output/ ; cd ..
python3 ../common/plot_output.py