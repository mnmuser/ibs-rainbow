CC=    gcc
LD=    gcc

DFLAGS=-D "aligned_alloc(x,y)=malloc(y)"
CFLAGS= -O2 -std=c11 -Wall -Wextra -fno-omit-frame-pointer $(DFLAGS)
INCPATH= -I/usr/local/include -I/opt/local/include -I/usr/include -I$(PROJ_DIR)
LDFLAGS=
LIBPATH= -L/usr/local/lib -L/opt/local/lib -L/usr/lib
LIBS=    -lcrypto

ifeq ($(shell uname -s),Darwin) # for macOS you need special openssl (homebrew)
DFLAGS=
LDFLAGS+= -L/usr/local/opt/openssl@1.1/lib
INCPATH+= -I/usr/local/opt/openssl@1.1/include/
endif

PROJ_DIR = rainbow_intern

ifeq ($(PROJ_DIR),avx2)
CFLAGS += -mavx2
CXXFLAGS += -mavx2
endif



SRCS = $(wildcard $(PROJ_DIR)/*.c)
SRCS_O = $(SRCS:.c=.o)
SRCS_O_ND = $(subst $(PROJ_DIR)/,,$(SRCS_O))

OBJ = $(SRCS_O_ND)

EXE= rainbow-genkey rainbow-sign rainbow-verify PQCgenKAT_sign


CSRC= $(wildcard *.c)


ifdef DEBUG
        CFLAGS=  -D_DEBUG_ -g -O1 -std=c11 -Wall -Wextra -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer -mavx2 #c99(sic!)
        CXXFLAGS= -D_DEBUG_ -g -O1 -Wall -Wextra -fsanitize=address -fsanitize=undefined -fno-exceptions -fno-rtti -nostdinc++ -mavx2
        LDFLAGS= -fsanitize=address -fsanitize=undefined
endif

ifdef GPROF
	CFLAGS += -pg
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

.PHONY: all tests tables clean

all: $(OBJ) $(EXE)


%-test: $(OBJ) %-test.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

%-benchmark: $(OBJ) %-benchmark.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

rainbow-genkey: $(OBJ) rainbow-genkey.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

rainbow-sign: $(OBJ) rainbow-sign.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

rainbow-verify: $(OBJ) rainbow-verify.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

PQCgenKAT_sign: $(OBJ) PQCgenKAT_sign.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<

%.o: $(PROJ_DIR)/%.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<

clean:
	rm *.o  rainbow-genkey rainbow-sign rainbow-verify PQCgenKAT_sign pk.txt sk.txt signature.txt || true;
