matlab_path = /usr/local/MATLAB/R2022b

sources = main.cpp matlab_util.cpp JobQueue.cpp util.cpp sdp_branch_and_bound.cpp ThreadPool.cpp
objects = $(sources:.cpp=.o)

all: bb

bb: $(objects)
	g++ -L $(matlab_path)/extern/bin/glnxa64/ -L $(gurobi_path)/lib/ -o $@ $^ -larmadillo -pthread -llapack -lopenblas -lMatlabDataArray -lMatlabEngine -Wl,-rpath,$(matlab_path)/extern/bin/glnxa64/

%.o: %.cpp
		g++ -g -Wall -O2 -o $@ -c $< -I $(matlab_path)/extern/include/ -I $(gurobi_path)/include/


.PHONY: clean
		
clean:
		rm -f bb *.o
