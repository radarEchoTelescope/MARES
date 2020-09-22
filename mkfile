TARGET = return_field
DIR= ../IceRayTracing/namespace/woROOT/

CC = g++
CFLAGS = -I$(DIR) #-Wall -O -g
# -O performs optimizations, -g enables debugging
LDLIBS= -lgsl

all: $(TARGET) #clean

OBJECTS= return_field.o $(DIR)IceRayTracing.o macro_scatter.o  antenna.o cascade.o


# For a given compilind order  xxxx: yyyy zzzz ttttt
#
# $@ is the left hand side of the order = xxxx
# $^ is the right hand side = yyyy zzzz tttt
# $< is the first term of the right hand side = yyyy
#
# Use this to avoid including headers in compilation lines.
# If you do, you will produce a yyyy.h.gsh file which is okay but might cause trouble.

$(DIR)IceRayTracing.o: $(DIR)IceRayTracing.cc $(DIR)IceRayTracing.hh
	$(CC) $(CFLAGS) -c $@ $<

antenna.o: antenna.cc antenna.hh $(DIR)IceRayTracing.hh
	$(CC) $(CFLAGS) -c $<

cascade.o: cascade.cc cascade.hh
	$(CC) $(CFLAGS) -c $<

macro_scatter.o: macro_scatter.cc macro_scatter.hh antenna.hh cascade.hh
	$(CC) $(CFLAGS) -c $<

return_field.o: return_field.cpp
	$(CC) $(CFLAGS) -c $<

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

.PHONY : clean
clean:
	rm -f *.o #$(DIR)*.o
