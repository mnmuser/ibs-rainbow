make clean
make all
./rainbow-gen-masterkey mpk.txt msk.txt
./rainbow-gen-usersk msk.txt David usk.txt
./rainbow-gen-userpk mpk.txt David upk.txt
./rainbow-sign usk.txt rainbow-sign.c | tee signature.txt
./rainbow-verify upk.txt signature.txt rainbow-sign.c