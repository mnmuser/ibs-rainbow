# ID-Rainbow

based on the code of NIST-candidate Rainbow (Round 2)

working currently only on Reference_Implementation


## How-To

You can use the compiled artifacts in this repository to try it on your linux machine (tested on Ubuntu 19.10).
Download the artifacts and run them in your shell:

### Generate Keys
```
./rainbow-genkey pk.txt sk.txt
```

### Sign a message
(you need a message-file to sign (e.g. message.txt))
```
./rainbow-sign sk.txt message.txt| tee signature.txt
```

### Verify a signature
```
./rainbow-genkey pk.txt signature.txt message.txt
```

## Choose method

At the moment you can choose security level and method in api.h and rainbow_config.h
