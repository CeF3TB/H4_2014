CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include -I../ -I./interface -I../CommonTools/interface

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lTMVA



all: makeAnalysisTree calibrateCeF3 drawEnergyScan alignTracking drawPositioningPlots




makeAnalysisTree: makeAnalysisTree.cpp RunHelper.o HodoCluster.o TagHelper.o EnergyCalibration.o AlignmentOfficer.o
	$(CC) -Wall $(INCLUDES) -o makeAnalysisTree makeAnalysisTree.cpp RunHelper.o HodoCluster.o TagHelper.o EnergyCalibration.o AlignmentOfficer.o $(ROOTFLAG) $(EXTRALIBS)

drawPositioningPlots: drawPositioningPlots.cpp DrawTools.o
	$(CC) -Wall $(INCLUDES) -o drawPositioningPlots drawPositioningPlots.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

positionAnalysis: positionAnalysis.cpp 
	$(CC) -Wall $(INCLUDES) -o positionAnalysis positionAnalysis.cpp $(ROOTFLAG) $(EXTRALIBS)

alignTracking: alignTracking.cpp  HodoCluster.o
	$(CC) -Wall $(INCLUDES) -o alignTracking alignTracking.cpp HodoCluster.o $(ROOTFLAG) $(EXTRALIBS)

calibrateCeF3: calibrateCeF3.cpp DrawTools.o
	$(CC) -Wall $(INCLUDES) -o calibrateCeF3 calibrateCeF3.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawEnergyScan: drawEnergyScan.cpp DrawTools.o
	$(CC) -Wall $(INCLUDES) -o drawEnergyScan drawEnergyScan.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)




HodoCluster.o: src/HodoCluster.cc
	$(CC) $(CFLAGS) -c src/HodoCluster.cc $(ROOTFLAG) $(EXTRALIBS)

DrawTools.o: ../CommonTools/src/DrawTools.cc
	$(CC) $(CFLAGS) -c ../CommonTools/src/DrawTools.cc $(ROOTFLAG) $(EXTRALIBS)

RunHelper.o: ../CommonTools/src/RunHelper.cc
	$(CC) $(CFLAGS) -c ../CommonTools/src/RunHelper.cc $(ROOTFLAG) $(EXTRALIBS)

TagHelper.o: src/TagHelper.cc
	$(CC) $(CFLAGS) -c src/TagHelper.cc $(ROOTFLAG) $(EXTRALIBS)

EnergyCalibration.o: src/EnergyCalibration.cc
	$(CC) $(CFLAGS) -c src/EnergyCalibration.cc $(ROOTFLAG) $(EXTRALIBS)

AlignmentOfficer.o: src/AlignmentOfficer.cc
	$(CC) $(CFLAGS) -c src/AlignmentOfficer.cc $(ROOTFLAG) $(EXTRALIBS)

FitTools.o: src/FitTools.cc
	$(CC) $(CFLAGS) -c src/FitTools.cc $(ROOTFLAG) $(EXTRALIBS)

PositionTools.o: src/PositionTools.cc
	$(CC) $(CFLAGS) -c src/PositionTools.cc $(ROOTFLAG) $(EXTRALIBS)




clean:
	rm *.o
