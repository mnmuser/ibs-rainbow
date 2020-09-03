# ID-Rainbow 
FOR RESEARCH PURPOSES ONLY

based on the code of NIST-candidate Rainbow (Round 2) with the adapted variables after the successful Band-Separation-Attack

only Reference_Implementation

## How-To

### Compiling
```
make all
```

### Generate Master-Keys
```
./rainbow-gen-masterkey mpk.txt msk.txt
```

### Generate a User-Secret-Key
```
./rainbow-gen-usersk msk.txt <YOUR_ID> usk.txt
```
### Generate a User-Public-Key
```
./rainbow-gen-userpk mpk.txt <YOUR_ID> upk.txt
```

### Sign a message
(you need a message-file to sign (e.g. message.txt))
```
./rainbow-sign usk.txt message.txt| tee signature.txt
```

### Verify a signature
```
./rainbow-genkey upk.txt signature.txt message.txt
```

## Choose method

At the moment you can choose security parameters and ID-length (_ID) in api.h and rainbow_config.h
