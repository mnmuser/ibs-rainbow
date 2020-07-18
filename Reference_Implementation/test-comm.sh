./rainbow-gen-masterkey mpk.txt msk.txt
./rainbow-gen-usersk msk.txt abcdef usk.txt
./rainbow-gen-userpk mpk.txt abcdef upk.txt
./rainbow-sign usk.txt rainbow-sign.c | tee signature.txt
./rainbow-verify upk.txt signature.txt rainbow-sign.c
