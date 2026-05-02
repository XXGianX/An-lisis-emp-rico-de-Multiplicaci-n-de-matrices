CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra
TARGET = matrix_analysis
SOURCES = main.cpp matriz.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	mkdir -p data
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJECTS)
	rm -rf data
