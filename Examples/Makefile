CXX_FLAGS=-g -fopenmp -I../include
CLANG_EXISTS:=$(shell which clang++)
ifeq "$(CLANG_EXISTS)" ""
	COMPILE_COMMAND=g++ $(CXX_FLAGS) -fmax-errors=2
else
	COMPILE_COMMAND=clang++ $(CXX_FLAGS) -ferror-limit=2
endif

all:	initialise runsIterator offsetRunsIterator draw2D crossIterator squareIterator iterator multiThreading

initialise:	Initialise/Initialise.cpp
		${COMPILE_COMMAND} -o build/Initialise Initialise/Initialise.cpp

runsIterator:	RunsIterator/RunsIterator.cpp
		${COMPILE_COMMAND} -o build/RunsIterator RunsIterator/RunsIterator.cpp

offsetRunsIterator:	OffsetRunsIterator/OffsetRunsIterator.cpp
		${COMPILE_COMMAND} -o build/OffsetRunsIterator OffsetRunsIterator/OffsetRunsIterator.cpp

draw2D:		Draw2D/Draw2D.cpp
		${COMPILE_COMMAND} -o build/Draw2D Draw2D/Draw2D.cpp

crossIterator:	CrossIterator/CrossIterator.cpp
		${COMPILE_COMMAND} -o build/CrossIterator CrossIterator/CrossIterator.cpp

squareIterator: SquareIterator/SquareIterator.cpp
		${COMPILE_COMMAND} -o build/SquareIterator SquareIterator/SquareIterator.cpp

iterator:	Iterator/Iterator.cpp
		${COMPILE_COMMAND} -o build/Iterator Iterator/Iterator.cpp

multiThreading:	MultiThreading/MultiThreading.cpp
		${COMPILE_COMMAND} -o build/MultiThreading MultiThreading/MultiThreading.cpp
