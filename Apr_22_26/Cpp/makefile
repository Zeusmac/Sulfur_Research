# =========================
# Compiler + Architecture
# =========================
CXX = g++
ARCH = -arch x86_64

# =========================
# ROOT configuration
# =========================
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lSpectrum

# =========================
# Project structure
# =========================
INCLUDES = -Iheaders

SRCS = GammaDB.cpp \
       PeakTracker.cpp \
       PeakFitter.cpp \
       FitGammaAnalyzer.cpp \
       MainPeakfit.cpp

OBJS = $(SRCS:.cpp=.o)

TARGET = gamma_fit

# =========================
# Compile flags
# =========================
CXXFLAGS = $(ARCH) -O2 -Wall -std=c++17 $(ROOTCFLAGS) $(INCLUDES)

# =========================
# Link flags
# =========================
LDFLAGS = $(ARCH) $(ROOTLIBS)

# =========================
# Build rules
# =========================
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# =========================
# Utilities
# =========================
clean:
	rm -f $(OBJS) $(TARGET)

rebuild: clean all