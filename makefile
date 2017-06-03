CXX=g++
CXXFLAGS=-g -Wall -O2
SOURCES=NetSim.cpp
BOOST_INC=../boost_1_60_0/
BOOST_LIB=../boost_1_60_0/stage/lib
EXECUTABLE=test
all:$(EXECUTABLE)

$(EXECUTABLE): $(SOURCES) 
	$(CXX) -I $(BOOST_INC) -L $(BOOST_LIB) -lboost_random -lboost_system $(CXXFLAGS) $(SOURCES) -o $(EXECUTABLE) 

clean:
	rm $(EXECUTABLE)

