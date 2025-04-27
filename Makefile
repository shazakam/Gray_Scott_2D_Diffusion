F95 = mpif90
OPTS = -O3 -Wall
LIBS=-L${BLASDIR} ${BLASLIB} -lpthread

OBJS = header.o create_matrix.o initial.o matmult.o rhs.o timestepping.o grayscott.o save_fields.o sparsegather_opt.o
OBJS_TEST =  header.o create_matrix.o initial.o matmult.o rhs.o timestepping.o save_fields.o sparsegather_opt.o testing.o

grayscott: $(OBJS)
	$(F95) $(OPTS) $(OBJS) -o grayscott $(LIBS)

testing: $(OBJS_TEST)
	$(F95) $(OPTS) $(OBJS_TEST) -o testing $(LIBS)

%.o:  %.f90 header.f90
	$(F95) $(OPTS) -c $<

clean:
	rm -rf grayscott testing *.o *.mod core.* solution.*

