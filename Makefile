CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra
TARGET   = benchmark
SOURCES  = main.cpp matriz.cpp

all: $(TARGET)

$(TARGET):
	mkdir -p data figures
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
	rm -rf data/ figures/

.PHONY: all run clean
