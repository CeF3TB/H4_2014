CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include -I../ -I./interface -I../CommonTools/interface

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lTMVA



all: 



makePedestalFile: makePedestalFile.cpp 
	$(CC) -Wall $(INCLUDES) -o makePedestalFile makePedestalFile.cpp $(ROOTFLAG) $(EXTRALIBS)

makeAnalysisTree: makeAnalysisTree.cpp RunHelper.o HodoCluster.o EnergyCalibration.o CalibrationUtility.o 
	$(CC) -Wall $(INCLUDES) -o makeAnalysisTree makeAnalysisTree.cpp RunHelper.o HodoCluster.o EnergyCalibration.o CalibrationUtility.o $(ROOTFLAG) $(EXTRALIBS)






HodoCluster.o: src/HodoCluster.cc
	$(CC) $(CFLAGS) -c src/HodoCluster.cc $(ROOTFLAG) $(EXTRALIBS)

DrawTools.o: ../CommonTools/src/DrawTools.cc
	$(CC) $(CFLAGS) -c ../CommonTools/src/DrawTools.cc $(ROOTFLAG) $(EXTRALIBS)

RunHelper.o: src/RunHelper.cc
	$(CC) $(CFLAGS) -c src/RunHelper.cc $(ROOTFLAG) $(EXTRALIBS)

CalibrationUtility.o: src/CalibrationUtility.cc
	$(CC) $(CFLAGS) -c src/CalibrationUtility.cc $(ROOTFLAG) $(EXTRALIBS)

EnergyCalibration.o: src/EnergyCalibration.cc
	$(CC) $(CFLAGS) -c src/EnergyCalibration.cc $(ROOTFLAG) $(EXTRALIBS)

FitTools.o: src/FitTools.cc
	$(CC) $(CFLAGS) -c src/FitTools.cc $(ROOTFLAG) $(EXTRALIBS)

PositionTools.o: src/PositionTools.cc
	$(CC) $(CFLAGS) -c src/PositionTools.cc $(ROOTFLAG) $(EXTRALIBS)




clean:
	rm *.o
