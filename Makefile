PROG = SpiralMeanField
CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++11 -march=native -I ~/Eigen_3.2.2 $(DEBUG)
LIBS = -llapack
OBJS = FinalSuperblock.o FreeFunctions.o Lanczos.o main.o $(PROG).o TheBlock.o
COMMONHS1 = GlobalHamiltonianParameters.h main.h
COMMONHS2 = $(COMMONHS1) Hamiltonian.h TheBlock.h FinalSuperblock.h
light = rm -f *.cpp~ *.h~ Makefile~
git = rm -f $(PROG) ./Output/*
deep = $(git) *.o

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LIBS) -o $(PROG) $(OBJS)

FinalSuperblock.o: $(COMMONHS2)

FreeFunctions.o: $(COMMONHS2) GlobalPrecisionParameters.h

Lanczos.o: $(COMMONHS1) Hamiltonian.h TheBlock.h GlobalPrecisionParameters.h

main.o: $(COMMONHS2) FreeFunctions.h GlobalPrecisionParameters.h ObservableOps.h

$(PROG).o: $(COMMONHS1) Hamiltonian.h TheBlock.h

TheBlock.o: $(COMMONHS2) GlobalPrecisionParameters.h

lightclean:
	$(light)

gitclean:
	$(git)

deepclean:
	$(deep)

clean:
	$(light)
	$(deep)

upload:
	scp -r *.cpp *.h Makefile $(OTHERS) knot.cnsi.ucsb.edu:~/$(DEST)

download:
	scp knot.cnsi.ucsb.edu:~/$(SOURCE)/Output/* Cluster
