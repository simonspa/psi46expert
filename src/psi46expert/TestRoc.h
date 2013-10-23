// Implementation of the tests at ROC level

#ifndef ROCTEST
#define ROCTEST

#include "BasePixel/Roc.h"
#include "TestDoubleColumn.h"
#include "TestParameters.h"

class TestRange;
class Test;
class TestRoc : public Roc
{

public:

    TestRoc(TBInterface * const aTBInterface, TestParameters * testParameters, const int aChipId, const int aHubId, const int aPortId, const int aoutChipPosition);
    void Execute(SysCommand &command);

    TestDoubleColumn * GetDoubleColumn(int column);
    TestPixel * GetPixel(int col, int row);
    TestPixel * GetTestPixel();
    TestParameters * GetTestParameters();
    TestRange * GetRange();

    // == Setting DACS ======================================

    int AdjustVana(double current0, double goalcurrent);
    void AdjustCalDelVthrComp();
    void AdjustCalDelVthrComp(int column, int row, int vcal, int belowNoise);
    int AdjustCalDel(int mode = 0);
    int AdjustCalDel(int col, int row, int mode = 0);
    int GetOptimalCalDel(int col, int row, int mode = 0);
    int AdjustPulseHeightRange(int range = 1);                                  ///< Adjusts the DACs to include all pulse height values in the ADC range
    void AdjustUltraBlackLevel(int ubLevel);

    // == Tests =============================================

    void ChipTest();
    void DoTest(Test * aTest);
    void DoIV(Test * aTest);
    void Test1();

    void PhError();
    void Rainbow();
    void Rainbow2();

    void TestM();
    void PowerOnTest(int nTests);
    void ADCSamplingTest();
    //  void DataTriggerLevelScan();
    double GetTemperature();
    void TrimVerification();

    void DoTrim();
    void DoTrimVcal();
    void DoTrimLow();
    void DoPhCalibration();
    void ThrMaps();
    void DoPulseShape();
    double DoPulseShape(int column, int row, int vcal);

    // == Helper routines ====================================

    int CountReadouts(int count);
    double Threshold(int sCurve[], int start, int sign, int step, double thrLevel);
    void SendSignals(int start, int stop, int step, int nTrig, char * dacName);
    void ReadSignals(int start, int stop, int step, int nTrig, int sCurve[]);
    int fitLeft(TH2D * his, int line, int maxEntries);
    int fitRight(TH2D * his, int line, int maxEntries);

    // == Histos =============================================

    TH1D * DACHisto();
    TH2D * TrimMap();

protected:

    TestParameters * testParameters;

};


#endif


