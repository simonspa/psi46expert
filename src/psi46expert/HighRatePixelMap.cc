#include <iostream>
using namespace std;

#include "HighRatePixelMap.h"

#include "TestModule.h"
#include "TestRoc.h"
#include "TestPixel.h"
#include "interface/Delay.h"
#include "BasePixel/RawPacketDecoder.h"
#include "DataFilter.h"
#include "interface/Log.h"

#include <TMath.h>
#include <TParameter.h>

HRPixelMap::HRPixelMap(TestRange * aTestRange, TestParameters * testParameters, TBInterface * aTBInterface)
{
    testRange = aTestRange;
    tbInterface = aTBInterface;
    this->testParameters = testParameters;
}

HRPixelMap::~HRPixelMap()
{

}

void HRPixelMap::ModuleAction(void)
{
    tbInterface->Flush();

    /* ??? */
    tbInterface->getCTestboard()->DataBlockSize(100);
    tbInterface->Flush();

    /* Unmask the ROC */
    int nroc = module->NRocs();
    psi::LogInfo() << "[HRPixelMap] Excluding masked pixels ... " << psi::endl;
    for (int i = 0; i < nroc; i++) {
        if (!testRange->IncludesRoc(i))
            continue;
        for (int col = 0; col < ROC_NUMCOLS; col++) {
            if (!testRange->IncludesColumn(i, col))
                continue;
            for (int row = 0; row < ROC_NUMROWS; row++) {
                if (testRange->IncludesPixel(i, col, row))
                    module->GetRoc(i)->EnablePixel(col, row);
            }
        }
    }
    tbInterface->Flush();

    /* Send a reset to the chip */
    tbInterface->getCTestboard()->Pg_SetCmd(0, PG_RESR);
    tbInterface->getCTestboard()->Pg_Single();
    tbInterface->Flush();

    /* Set clock stretch */
    //if (testParameters->HRPixelMapClockStretch > 1)
    //    tbInterface->SetClockStretch(STRETCH_AFTER_CAL, testParameters->HRPixelMapStretchDelay, testParameters->HRPixelMapClockStretch);

    /* Get the digital and analog voltages / currents */
    psi::LogInfo() << "[HRPixelMap] Measuring chip voltages and currents ..." << psi::endl;
    TParameter<float> vd("hr_pixelmap_digital_voltage", tbInterface->GetVD());
    TParameter<float> id("hr_pixelmap_digital_current", tbInterface->GetID());
    TParameter<float> va("hr_pixelmap_analog_voltage", tbInterface->GetVA());
    TParameter<float> ia("hr_pixelmap_analog_current", tbInterface->GetIA());
    vd.Write();
    id.Write();
    va.Write();
    ia.Write();

    float seconds = testParameters->HRPixelMapAquisitionTime;
    int repetitions = testParameters->HRPixelMapRepetitions;

    /* Data filters */
    RawData2RawEvent rs;
    RawEventDecoder ed(nroc, module->GetRoc(0)->has_analog_readout(), module->GetRoc(0)->has_row_address_inverted());
    HitMapper hm(nroc, seconds * repetitions);
    EventCounter count;
    MultiplicityHistogrammer mh;
    PulseHeightHistogrammer phh;
    ConfigParameters * configParameters = ConfigParameters::Singleton();
    phh.LoadCalibration(nroc, configParameters->directory);
    int tct = tbInterface->GetParameter("tct");
    int ttk = tbInterface->GetParameter("ttk");
    int deserAdjust = tbInterface->GetParameter("deserAdjust");

    /* Repeat measurements multiple times to collect statistics */
    for (int rep = 0; rep < repetitions; rep++) {
        if (testParameters->HRPixelMapRepetitions > 1)
            psi::LogInfo() << "[HRPixelMap] Masuring iteration " << rep + 1 << "/" << testParameters->HRPixelMapRepetitions << " ..." << psi::endl;
        /* Prepare the data aquisition (store to testboard RAM) */
        int memsize = tbInterface->getCTestboard()->Daq_Open(30000000);
        tbInterface->getCTestboard()->Daq_Select_Deser160(deserAdjust);

        /* Enable DMA (direct memory access) controller */
        tbInterface->getCTestboard()->Daq_Start();

        /* Set the trigger frequency (f = ???) */

        tbInterface->getCTestboard()->Pg_SetCmd(0, PG_RESR);
        tbInterface->getCTestboard()->Pg_Single();
        tbInterface->Flush();

        /* Issue continuous Reset-(Calibrate-)Trigger-Token pattern */
        //tbInterface->getCTestboard()->Pg_SetCmd(0, PG_CAL + tct);
        tbInterface->getCTestboard()->Pg_SetCmd(0, PG_TRG + ttk);
        tbInterface->getCTestboard()->Pg_SetCmd(1, PG_TOK);
        tbInterface->getCTestboard()->Pg_Loop(testParameters->HRPixelMapTriggerRate);

        tbInterface->Flush();


        for (float t = seconds; t >= 1; t--) {
            cout << "\r[HRPixelMap] Taking data (" << t << " seconds) ... ";
            cout.flush();
            gDelay->Mdelay(1000);
        }
        cout << "\r[HRPixelMap] Taking data (" << (seconds - (int)(seconds)) << " seconds) ... ";
        cout.flush();
        gDelay->Mdelay((int)((seconds - (int)(seconds)) * 1000));
        cout << "\r[HRPixelMap] Taking data (" << seconds << " seconds) ... done" << endl;

        /* Stop triggering */
        tbInterface->getCTestboard()->Pg_Stop();
        tbInterface->getCTestboard()->Pg_SetCmd(0, PG_RESR);
        tbInterface->getCTestboard()->Pg_Single();
        tbInterface->Flush();

        /* Wait for data aquisition to finish */
        gDelay->Mdelay(100);


        /* Disable data aquisition */
        tbInterface->getCTestboard()->Daq_Stop();
        tbInterface->Flush();

        /* Number of words stored in memory */
        //int nwords = (data_end - data_pointer) / 2;
        //psi::LogInfo() << "[HRPixelMap] Megabytes in RAM: " << nwords * 2. / 1024. / 1024. << psi::endl;

        /* Prepare data decoding */
        RAMRawDataReader rd(tbInterface->getCTestboard(), memsize);

        /* Decoding chain */
        rd >> rs >> ed >> hm >> count >> mh >> phh >> pipe_end;

        /* Free the memory in the RAM */
        tbInterface->getCTestboard()->Daq_Close();
    }

    /* Store histograms */
    int core_hits = 0;
    for (int i = -2; i < nroc; i++) {
        /* -2: Module with double sized edges, -1: Module */
        TH2I * map = (TH2I *) hm.getHitMap(i)->Clone();
        histograms->Add(map);

        /* Make dcol map and hit distribution */
        if (i >= 0) {
            for (int c = 3; c <= 50; c++) {
                for (int r = 1; r <= 79; r++) {
                    core_hits += map->GetBinContent(c, r);
                }
            }

            TH1I * dcol_map = new TH1I(Form("dcol_map_C%i", i), Form("DCol hit map ROC %i", i), 26, 0, 26);
            int x, y, z;
            map->GetMaximumBin(x, y, z);
            z = map->GetBinContent(x, y);
            TH1I * hit_dist = new TH1I(Form("hit_dist_C%i", i), Form("Hit distribution ROC %i", i), z > 100 ? 100 : z, 0, z);
            for (int dcol = 0; dcol < 26; dcol++) {
                int sum = 0;
                for (int row = 0; row < 80; row++) {
                    sum += map->GetBinContent(2 * dcol + 1, row + 1);
                    sum += map->GetBinContent(2 * dcol + 2, row + 1);
                    hit_dist->Fill(map->GetBinContent(2 * dcol + 1, row + 1));
                    hit_dist->Fill(map->GetBinContent(2 * dcol + 2, row + 1));
                }
                dcol_map->SetBinContent(dcol + 1, sum);
            }
            dcol_map->Sumw2();
            dcol_map->SetMinimum(0);
            histograms->Add(dcol_map);
            histograms->Add(hit_dist);

            TH1I * multi = (TH1I *) mh.getRocMultiplicity(i)->Clone();
            histograms->Add(multi);
        }
    }
    histograms->Add((TH1I *) hm.getHitsVsTimeDcol()->Clone());
    histograms->Add((TH1I *) hm.getHitsVsTimeRoc()->Clone());
    histograms->Add((TH1I *) phh.getPulseHeightDistribution()->Clone());
    histograms->Add((TH2F *) phh.getPulseHeightMap()->Clone());
    histograms->Add((TH2F *) phh.getPulseHeightWidthMap()->Clone());
    histograms->Add((TH1I *) phh.getCalPulseHeightDistribution()->Clone());
    histograms->Add((TH2F *) phh.getCalPulseHeightMap()->Clone());
    histograms->Add((TH2F *) phh.getCalPulseHeightWidthMap()->Clone());

    float active_area = nroc * (52 - 2 * 2) * (80 - 1) * 0.01 * 0.015; /* cm2 */
    float active_time = count.DataCounter * 25e-9; /* s */
    if (testParameters->HRPixelMapClockStretch > 1)
        active_time *= testParameters->HRPixelMapClockStretch;
    TH2I * map = (TH2I *) hm.getHitMap(-1);
    psi::LogInfo() << "[HRPixelMap] Number of triggers: " << count.DataCounter << psi::endl;
    psi::LogInfo() << "[HRPixelMap] Number of hits: " << map->GetEntries() << psi::endl;
    psi::LogInfo() << "[HRPixelMap] Rate: " << (core_hits / active_time / active_area / 1e6);
    psi::LogInfo() << " +/- " << (TMath::Sqrt(core_hits) / active_time / active_area / 1e6);
    psi::LogInfo() << " megahits / s / cm2" << psi::endl;
    psi::LogInfo() << "[HRPixelMap] Number of ROC sequence problems: " << count.RocSequenceErrorCounter << psi::endl;
    psi::LogInfo() << "[HRPixelMap] Number of decoding problems: " << ed.GetDecodingErrors() << psi::endl;

    TParameter<float> triggers("pixelmap_triggers", count.DataCounter);
    triggers.Write();

    /* Disable clock stretch */
    //if (testParameters->HRPixelMapClockStretch > 1)
    //    tbInterface->SetClockStretch(0, 0, 0);

    /* Reset the chip */
    tbInterface->getCTestboard()->Pg_Stop();
    tbInterface->getCTestboard()->Pg_SetCmd(0, PG_RESR);
    tbInterface->getCTestboard()->Pg_Single();
    tbInterface->Flush();
}
