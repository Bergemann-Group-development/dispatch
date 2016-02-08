#FC = ifort
#FOPT =  -axAVX -xSSE4.1

FC = gfortran
FOPT =  -O3 -ftree-vectorize

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FOPT) $< -c

default: test_transfer test_transfer2

test_transfer2: test_transfer2.x
	@echo ""
	./test_transfer2.x
TEST_TRANSFER2 = test_transfer_sin2.o transfer_feautrier_fast2.o
test_transfer2.x: $(TEST_TRANSFER2)
	$(FC) $(TEST_TRANSFER2) -o test_transfer2.x

test_transfer: test_transfer.x
	@echo ""
	./test_transfer.x
TEST_TRANSFER = test_transfer.o transfer_feautrier.o
test_transfer.x: $(TEST_TRANSFER)
	$(FC) $(TEST_TRANSFER) -o test_transfer.x
