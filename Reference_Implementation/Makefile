CC=    gcc
LD=    gcc

DFLAGS= -D "aligned_alloc(x,y)=malloc(y)"
CFLAGS= -O2 -std=c11 -Wall -Wextra -fno-omit-frame-pointer $(DFLAGS)
INCPATH= -I/usr/local/include -I/opt/local/include -I/usr/include -I$(PROJ_DIR)
LDFLAGS=
LIBPATH= -L/usr/local/lib -L/opt/local/lib -L/usr/lib
LIBS=    -lcrypto -lm


PROJ_DIR = rainbow_intern

ifeq ($(PROJ_DIR),avx2)
CFLAGS += -mavx2
CXXFLAGS += -mavx2
endif



SRCS = $(wildcard $(PROJ_DIR)/*.c)
SRCS_O = $(SRCS:.c=.o)
SRCS_O_ND = $(subst $(PROJ_DIR)/,,$(SRCS_O))

OBJ = $(SRCS_O_ND)

EXE= rainbow-gen-masterkey rainbow-gen-userpk rainbow-gen-usersk rainbow-sign rainbow-verify


CSRC= $(wildcard *.c)


ifdef DEBUG_FAST # normal debugging with sanitizers needs 3 minutes, so I decided for a fast variant
        CFLAGS=  -D_DEBUG_ -g -O2 -march=native -mtune=native -std=c11 -Wall -Wextra
        CXXFLAGS= -D_DEBUG_ -g -O2 -march=native -mtune=native -Wall -Wextra
endif

ifdef DEBUG
        CFLAGS=  -D_DEBUG_ -g -Og -std=c11 -Wall -Wextra  -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer -mavx2 #c99(sic!)
        CXXFLAGS= -D_DEBUG_ -g -Og -Wall -Wextra -fsanitize=address -fsanitize=undefined -fno-exceptions -fno-rtti -nostdinc++ -mavx2
        LDFLAGS= -fsanitize=address -fsanitize=undefined
endif

# genkey-benchmark (gcc10):
#-> O3 and march=native 12 seconds
#-> O3 and march=znver2 12 seconds -> Really new instruction set for Zen 2 (is used by native)
#-> O2 and march=native 18 seconds
#-> O2 w/o march 		17 seconds
#-> O3 w/o march		13 seconds
#-> Debug O3 march=native 3:10 minutes
#-> DEBUG O3 march=native w/o sanitizers 14 seconds
#-> O3 optimizes parts of the keys out.. -> back to O2

ifeq ($(shell uname -s),Darwin) # for macOS you need special openssl (homebrew)
DFLAGS=
LIBPATH= -L/usr/local/lib -L/usr/lib
LDFLAGS+= -L/usr/local/opt/openssl@1.1/lib
INCPATH+= -I/usr/local/opt/openssl@1.1/include/
endif

ifdef GPROF
        CFLAGS += -pg
        XXFLAGS += -pg
        LDFLAGS += -pg
endif

.PHONY: all tests tables clean

all: $(OBJ) $(EXE)


%-test: $(OBJ) %-test.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

%-benchmark: $(OBJ) %-benchmark.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

rainbow-gen-masterkey: $(OBJ) rainbow-gen-masterkey.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

rainbow-gen-userpk: $(OBJ) rainbow-gen-userpk.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

rainbow-gen-usersk: $(OBJ) rainbow-gen-usersk.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

rainbow-sign: $(OBJ) rainbow-sign.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

rainbow-verify: $(OBJ) rainbow-verify.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<

%.o: $(PROJ_DIR)/%.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<

clean_without_keys:
	rm *.o  rainbow-gen-masterkey rainbow-gen-usersk rainbow-gen-userpk rainbow-sign rainbow-verify signature.txt || true;

clean:
	rm *.o  rainbow-gen-masterkey rainbow-gen-usersk rainbow-gen-userpk rainbow-sign rainbow-verify mpk.txt msk.txt usk.txt upk.txt signature.txt || true;
