// Trim functions

#ifndef THRESHOLDMAP
#define THRESHOLDMAP

#include <TH2D.h>
#include "TestRoc.h"
#include "TestRange.h"

class ThresholdMap
{

public:
	ThresholdMap();
	
	TH2D* GetMap(char* mapName, TestRoc *roc, TestRange *testRange, int nTrig, int mode = -1);
// 	bool CheckMap();
	void MeasureMap(char* mapName, TestRoc *roc, TestRange *testRange, int nTrig);
	void SetParameters(char *mapName, int mode = -1);
	void SetCals();
	void SetXTalk();
	void SetDoubleWbc();
	void SetSingleWbc();
	void SetReverseMode();
		
protected:

	TH2D *histo;
	int dacReg;
	bool cals, reverseMode, xtalk, doubleWbc;

};


#endif

