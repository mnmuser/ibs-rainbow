./rainbow-genkey mpk.txt msk.txt
./rainbow-sign sk.txt rainbow-sign.c | tee signature.txt
./rainbow-verify pk.txt signature.txt rainbow-sign.c
