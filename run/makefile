CXX = g++
CXXFLAGS = -fdiagnostics-color=always -g -pthread -O3
TARGET = ./hin_sscs
SRCS = ${wildcard ../src/*.cpp}

all: ${TARGET}

${TARGET}: ${SRCS}
	${CXX} ${CXXFLAGS} ${SRCS} -o ${TARGET}

clean:
	rm -f ${TARGET}