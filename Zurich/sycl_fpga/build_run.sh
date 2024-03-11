
mkdir -p ./build;
cd ./build && cmake .. ; make ;mv ./$1 ../../
cd ../../ && ./$1