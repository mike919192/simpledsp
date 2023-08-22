all: run_test


PICKY_FLAGS:=-Wall -Wextra

INC_FLAGS:= -Iinclude

GCC_FLAGS:=$(PICKY_FLAGS) -march=x86-64-v3 -O2

CATCH2_SRC:= test/catch2/catch_amalgamated.cpp
CATCH2_o:= catch_amalgamated.o


.PHONY: clean
clean:
	rm $(CATCH2_o)
	rm testFFT.out

catch_amalgamated.o: test/catch2/catch_amalgamated.cpp test/catch2/catch_amalgamated.hpp
	g++ $(INC_FLAGS) $(GCC_FLAGS) -c -Iinclude $(CATCH2_SRC)

testFFT.out: test/testFFT.cpp include/sdsp/fft.h catch_amalgamated.o
	g++ $(INC_FLAGS) $(GCC_FLAGS) -o testFFT.out -Iinclude test/testFFT.cpp $(CATCH2_o)

run_test: testFFT.out
	./testFFT.out


