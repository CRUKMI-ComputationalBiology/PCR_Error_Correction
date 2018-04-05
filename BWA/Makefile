
OBJS =  BWA.o argProcessor.o
EXECUTABLE = BWA 

CXX=mpicxx
CXXFLAG=-Wall -O3 
COMPILE=$(CXX) $(CXXFLAG) -c

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS) 
	$(CXX) -o $(EXECUTABLE) $(OBJS)

%.o: %.cpp
	$(COMPILE) -o $@ $<

clean:
	-rm $(OBJS) $(EXECUTABLE)

