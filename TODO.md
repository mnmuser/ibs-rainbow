# TODOs for MA

## needed for MA

- integrate ID to hash function ✔️
- construct Master Secret Key ✔
    - keep sk-generating
    - sk -> msk
    - (add S' and T' (based on ID?))
- construct user secret key 
    - multiply sk-parts with ID-hash
    - _Problem: Platzhalter für ID in pk_
    - Bilde endlichen Körper (GF16) auf ID-Hash ab?
    - (build real key with S' and T'...)
    - (base S and T on ID)
- construct master pk
    - keep everything
- construct user public key
    - multiply pk-parts with ID-hash (like sk)
- check verifying
- check if this is working

msk = S, S', F, T', T   
usk = S * S', S'⁻1 * F * T'⁻1 , T' * T

## optional things

- get Level I, III and V in one directory (set compile flags in makefile) (almost no differences between directories) ✔
- get Classic, Cyclic and Cyclic_compressed in one file (set compile flags in makefile) ✔
- same for AVX2 and SSE3?