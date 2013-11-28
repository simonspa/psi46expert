#include <iostream>
using namespace std;

#include "HighRateEfficiency.h"

#include "TestModule.h"
#include "TestRoc.h"
#include "TestPixel.h"
#include "interface/Delay.h"
#include "BasePixel/RawPacketDecoder.h"
#include "DataFilter.h"
#include "interface/Log.h"

#include <TMath.h>
#include <TParameter.h>

HREfficiency::HREfficiency(TestRange * aTestRange, TestParameters * testParameters, TBInterface * aTBInterface)
{
    testRange = aTestRange;
    tbInterface = aTBInterface;
    this->testParameters = testParameters;
}

HREfficiency::~HREfficiency()
{

}

void HREfficiency::ModuleAction(void)
{
    TBAnalogInterface * ai = (TBAnalogInterface *) tbInterface;
    ai->Flush();


    /* Unmask all ROCs */
    for (int i = 0; i < module->NRocs(); i++) {
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
    ai->Flush();

    /* Send a reset to the chip */
    ai->getCTestboard()->Pg_SetCmd(0, PG_RESR);
    ai->getCTestboard()->Pg_Single();
    ai->Flush();

    gDelay->Mdelay(10);

    /* Prepare the data aquisition (store to testboard RAM) */
    int memsize = ai->getCTestboard()->Daq_Open(30000000);
    ai->getCTestboard()->Daq_Select_Deser160(4);

    /* Enable DMA (direct memory access) controller */
    ai->getCTestboard()->Daq_Start();

    /* Get the digital and analog voltages / currents */
    psi::LogInfo() << "Measuring chip voltages and currents ..." << psi::endl;
    TParameter<float> vd("hr_efficiency_digital_Voltage", ai->GetVD());
    TParameter<float> id("hr_efficiency_digital_current", ai->GetID());
    TParameter<float> va("hr_efficiency_analog_voltage", ai->GetVA());
    TParameter<float> ia("hr_efficiency_analog_current", ai->GetIA());
    vd.Write();
    id.Write();
    va.Write();
    ia.Write();

    int tct = ai->GetParameter("tct");
    int ttk = ai->GetParameter("ttk");
    /* Set the number of triggers. Total triggers for all pixels: 4160 * ntrig */
    int ntrig = testParameters->HREfficiencyTriggers;

    /* iterate over columns and rows to get each pixel efficiency */
    for (int col = 0; col < 52; col++) {
        cout << "\rSending calibrate signals ... " << ((int)(100 * col / 51.)) << " % " << flush;
        for (int row = 0; row < 80; row++) {
            /* Arm the pixel */
            for (int i = 0; i < module->NRocs(); i++) {
                if (testRange->IncludesPixel(i, col, row))
                    module->GetRoc(i)->ArmPixel(col, row);
            }
            ai->CDelay(5000);
            ai->Flush();

            /* Send a reset to the chip */
            ai->getCTestboard()->Pg_SetCmd(0, PG_RESR);
            ai->getCTestboard()->Pg_Single();
            ai->CDelay(500);
            ai->Flush();

            /* send ntrig triggers with calibrates */
            ai->getCTestboard()->Pg_SetCmd(0, PG_CAL + tct);
            ai->getCTestboard()->Pg_SetCmd(1, PG_TRG + ttk);
            ai->getCTestboard()->Pg_SetCmd(2, PG_TOK);
            for (int t = 0; t < ntrig; t++) {
                ai->getCTestboard()->Pg_Single();
                ai->CDelay(500);
            }
            ai->Flush();

            /* Disarm the pixel */
            for (int i = 0; i < module->NRocs(); i++) {
                if (testRange->IncludesPixel(i, col, row))
                    module->GetRoc(i)->ClrCal();
            }
            ai->Flush();
        }
    }
    cout << endl;

    /* Stop triggering */
    ai->getCTestboard()->Pg_SetCmd(0, PG_RESR);
    ai->getCTestboard()->Pg_Single();
    ai->Flush();

    /* Wait for data aquisition to finish */
    gDelay->Mdelay(100);

    /* Disable data aquisition */
    ai->getCTestboard()->Daq_Stop();
    ai->Flush();


    /* Number of words in memory */
    //int nwords = (data_end - data_pointer) / 2;
    //psi::LogInfo() << "Megabytes in RAM: " << nwords * 2. / 1024. / 1024. << psi::endl;

    /* Prepare data decoding */
    int nroc = module->NRocs();
    RAMRawDataReader rd(ai->getCTestboard(), memsize);
    RawData2RawEvent rs;
    RawEventDecoder ed(nroc, module->GetRoc(0)->has_analog_readout(), module->GetRoc(0)->has_row_address_inverted());
    EfficiencyMapper em(nroc, ntrig);

    /* Decoding chain */
    rd >> rs >> ed >> em >> pipe_end;

    /* Store histograms */
    float background = 0;
    float background_core = 0;
    int background_core_pixels = 0;
    float efficiency = 0;
    float core_efficiency = 0;
    int core_efficiency_pixels = 0;
    for (int i = -1; i < nroc; i++) {
        TH2I * effmap = (TH2I *) em.getEfficiencyMap(i)->Clone();
        effmap->SetMinimum(0);
        effmap->SetMaximum(ntrig);
        histograms->Add(effmap);
        TH1I * effdist = (TH1I *) em.getEfficiencyDist(i)->Clone();
        histograms->Add(effdist);
        TH2I * bkgmap = (TH2I *) em.getBackgroundMap(i)->Clone();
        histograms->Add(bkgmap);
        if (i == -1) {
            efficiency = effdist->GetMean();
            background = bkgmap->GetEntries();
        }
        if (i >= 0) {
            for (int c = 2; c < 50; c++) {
                for (int r = 0; r < 80; r++) {
                    if (!testRange->IncludesPixel(i, c, r))
                        continue;
                    if (r < 79) {
                        background_core += bkgmap->GetBinContent(c + 1, r + 1);
                        background_core_pixels += 1;
                    }
                    core_efficiency += effmap->GetBinContent(c + 1, r + 1);
                    core_efficiency_pixels += 1;
                }
            }
        }
    }
    float eff_err = efficiency / 100.0;
    eff_err = TMath::Sqrt(eff_err * (1 - eff_err) * (4160 * ntrig)) * 100.0 / (4160 * ntrig);

    core_efficiency *= 1.0 / (ntrig * core_efficiency_pixels);
    float core_eff_err = TMath::Sqrt(core_efficiency * (1 - core_efficiency) * (core_efficiency_pixels * ntrig)) * 100.0 / (core_efficiency_pixels * ntrig);
    core_efficiency *= 100.0;

    /* Sensor area, edge double columns, top row, and double column under test excluded */
    float active_area = background_core_pixels * 0.01 * 0.015; /* cm2 */
    float active_time = ntrig * 4160 * 25e-9; /* s */
    float rate = background_core / active_area / active_time / 1e6;
    float rate_err = TMath::Sqrt(background_core) / active_area / active_time / 1e6;
    psi::LogInfo() << Form("%-19s %8i", "Number of triggers:", ntrig * 4160) << psi::endl;
    psi::LogInfo() << Form("%-19s %8i", "Number of hits:", background) << psi::endl;  
    psi::LogInfo() << Form("%-19s %8.3f +/- %.3f MHz / cm2", "Rate:", rate, rate_err) << psi::endl;;
    psi::LogInfo() << Form("%-19s %8.3f +/- %.3f %%", "Overall efficiency:", efficiency, eff_err) << psi::endl;
    psi::LogInfo() << Form("%-19s %8.3f +/- %.3f %%", "Core efficiency:", core_efficiency, core_eff_err) << psi::endl;
    psi::LogInfo() << Form("%-19s %8i", "Decoding problems:", ed.GetDecodingErrors()) << psi::endl;

    /* Free the memory in the RAM */
    ai->getCTestboard()->Daq_Close();

    /* Reset the chip */
    ai->Single(RES);
    ai->Flush();
}
