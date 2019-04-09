CFLAGS = -Wall -O3
LDFLAGS = 
CC = g++
LD = g++
EXEC = eos.x
LIBS =
OBJS = atmdatb.o \
	   molecb.o \
	   partf.o \
	   saha.o \
	   leeabun.o \
	   Wittmann.o

######### Suffix rules ########################################
.SUFFIXES :    .o .cc

.cc.o:
	$(CC) $(CFLAGS) -c $<

all: 
	$(MAKE) $(OBJS)
	$(LD) -o $(EXEC) $(OBJS) $(LIBS) $(LDFLAGS)

clean:
	rm -f $(EXEC) $(OBJS)
