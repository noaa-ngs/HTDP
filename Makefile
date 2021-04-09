RM = rm -f

CC = CC 
FC = f90

CFLAGS = -C   
INCLUDES = 
LIBS = 
.SUFFIXES : .o .f 
# rule to create .o files from .f files
.f.o:
	$(RM) $*.o
	$(FC) -c  $(CFLAGS) $(INCLUDES) $*.f

# rule to create executables from .o files
.o:
	$(RM) $@
	$(FC) $(CFLAGS) $(INCLUDES) -o $@ $<  $(LIBS) 
                         
all: htdp

htdp.o:          htdp.f
initvl.o:        initvl.f
initeq.o:        initeq.f
initbd.o:        initbd.f
initps.o:        initps.f

htdp:   htdp.o initvl.o initeq.o initbd.o initps.o 
	$(FC) $(CFLAGS) $(INCLUDES) -o $@ $?  $(LIBS) 

