
void FormatCanvas(TCanvas* canvas);

void FormatHistLong(TH1D* hist, int colour, int style);

void FormatHist(TH1D* hist, int colour, int style);

void CalcPropLikeRatio(const int nProbBins = 25)
{

  const int nParticles = 5;

  TString testFileName =
    "/data/t2k/phsxvg/DUNEPid/data/4APA/EventSamples/train_muon_realflux_split_2.root";

  TString trainFileName[nParticles] = {
    "/data/t2k/phsxvg/DUNEPid/data/4APA/TrainMVA/mvaPlots_muon_all.root",
    "/data/t2k/phsxvg/DUNEPid/data/4APA/TrainMVA/mvaPlots_electron_all.root",
    "/data/t2k/phsxvg/DUNEPid/data/4APA/TrainMVA/mvaPlots_proton_all.root",
    "/data/t2k/phsxvg/DUNEPid/data/4APA/TrainMVA/mvaPlots_pich_all.root",
    "/data/t2k/phsxvg/DUNEPid/data/4APA/TrainMVA/mvaPlots_photon_all.root"};

  TString weightFileName[nParticles] = {
    "/data/t2k/phsxvg/DUNEPid/data/4APA/weights/muon_all_BDT.weights.xml",
    "/data/t2k/phsxvg/DUNEPid/data/4APA/weights/electron_all_BDT.weights.xml",
    "/data/t2k/phsxvg/DUNEPid/data/4APA/weights/proton_all_BDT.weights.xml",
    "/data/t2k/phsxvg/DUNEPid/data/4APA/weights/pich_all_BDT.weights.xml",
    "/data/t2k/phsxvg/DUNEPid/data/4APA/weights/photon_all_BDT.weights.xml"};

  TFile* fileForBinning = TFile::Open(trainFileName[0]);
  fileForBinning->cd("Method_BDT/BDT");

  int numbins = MVA_BDT_Train_S->GetNbinsX();
  const int nbins = numbins;

  fileForBinning->Close();

  double likeRatioS[nParticles][nbins] = {0};

  double binWidth[nParticles], MVAMin[nParticles], MVAMax[nParticles];

  for (int p = 0; p < nParticles; p++) {
    cout << endl;
    cout << "Calculating likelihood ratios from file " << trainFileName[p] << endl;
    cout << endl;

    TFile* trainFile = TFile::Open(trainFileName[p]);
    trainFile->cd("Method_BDT/BDT");

    if (MVA_BDT_Train_S->GetNbinsX() != nbins) {
      cerr << "MVA_BDT_Train_S->GetNbinsX() not equal to nbins in file " << trainFileName[p]
           << endl;
      exit(1);
    }
    if (MVA_BDT_Train_B->GetNbinsX() != nbins) {
      cerr << "MVA_BDT_Train_B->GetNbinsX() not equal to nbins in file " << trainFileName[p]
           << endl;
      exit(1);
    }

    binWidth[p] = MVA_BDT_Train_S->GetBinWidth(1);
    MVAMin[p] = MVA_BDT_Train_S->GetBinLowEdge(1);
    MVAMax[p] = MVA_BDT_Train_S->GetBinLowEdge(nbins) + MVA_BDT_Train_S->GetBinWidth(nbins);

    if (MVA_BDT_Train_B->GetBinWidth(1) != binWidth[p]) {
      cerr << "MVA_BDT_Train_B->GetBinWidth(1) not equal to binWidth in file " << trainFileName[p]
           << endl;
      exit(1);
    }
    if (MVA_BDT_Train_B->GetBinLowEdge(1) != MVAMin[p]) {
      cerr << "MVA_BDT_Train_B->GetBinLowEdge(1) not equal to MVAMin in file " << trainFileName[p]
           << endl;
      exit(1);
    }
    if (MVA_BDT_Train_B->GetBinLowEdge(nbins) + MVA_BDT_Train_B->GetBinWidth(nbins) != MVAMax[p]) {
      cerr << "MVA_BDT_Train_B->GetBinLowEdge(nbins) + MVA_BDT_Train_B->GetBinWidth(nbins) not "
              "equal to MVAMax in file "
           << trainFileName[p] << endl;
      exit(1);
    }

    for (int i = 1; i <= nbins; i++) {
      likeRatioS[p][i] =
        MVA_BDT_Train_S->GetBinContent(i) / TMath::Max(0.01, MVA_BDT_Train_B->GetBinContent(i));
    }
  }

  TFile* testfile = TFile::Open(testFileName);
  testfile->cd();
  TTree* PIDVar = (TTree*)testfile->Get("mvaTree");

  float evalRatio, coreHaloRatio, concentration, conicalness, dEdxStart, dEdxEnd, dEdxEndRatio;

  PIDVar->SetBranchAddress("evalRatio", &evalRatio);
  PIDVar->SetBranchAddress("coreHaloRatio", &coreHaloRatio);
  PIDVar->SetBranchAddress("concentration", &concentration);
  PIDVar->SetBranchAddress("conicalness", &conicalness);
  PIDVar->SetBranchAddress("dEdxStart", &dEdxStart);
  PIDVar->SetBranchAddress("dEdxEnd", &dEdxEnd);
  PIDVar->SetBranchAddress("dEdxEndRatio", &dEdxEndRatio);

  double MVAValue;
  int MVABin;

  float evalRatioR, coreHaloRatioR, concentrationR, conicalnessR, dEdxStartR, dEdxEndR,
    dEdxEndRatioR;

  vector<double> likeRatio[nParticles];
  vector<double>::iterator likeIter;

  TH1D* hMVAValue[nParticles] = {0};
  for (int p = 0; p < nParticles; p++)
    hMVAValue[p] = new TH1D(Form("hMVAValue%d", p), "", nbins, -0.5, 0.5);

  for (int p = 0; p < nParticles; p++) {
    cout << endl;
    cout << "Reading weights from file " << weightFileName[p] << endl;

    TMVA::Reader* reader = new TMVA::Reader();

    reader->AddVariable("evalRatio", &evalRatioR);
    reader->AddVariable("coreHaloRatio", &coreHaloRatioR);
    reader->AddVariable("concentration", &concentrationR);
    reader->AddVariable("conicalness", &conicalnessR);
    reader->AddVariable("dEdxStart", &dEdxStartR);
    reader->AddVariable("dEdxEnd", &dEdxEndR);
    reader->AddVariable("dEdxEndRatio", &dEdxEndRatioR);

    reader->BookMVA("BDT", weightFileName[p]);

    for (int i = 0; i < PIDVar->GetEntries(); ++i) {
      PIDVar->GetEntry(i);

      evalRatioR = evalRatio;
      coreHaloRatioR = coreHaloRatio;
      concentrationR = concentration;
      conicalnessR = conicalness;
      dEdxStartR = dEdxStart;
      dEdxEndR = dEdxEnd;
      dEdxEndRatioR = dEdxEndRatio;

      MVAValue = reader->EvaluateMVA("BDT");

      MVABin = 1 + int((MVAValue - MVAMin[p]) / binWidth[p]);

      hMVAValue[p]->Fill(MVAValue);

      likeRatio[p].push_back(likeRatioS[p][MVABin]);
    }

    delete reader;
  }

  TCanvas* MVAVal = new TCanvas("MVAVal", "MVAVal", 800, 600);
  FormatCanvas(MVAVal);
  hMVAValue[0]->GetXaxis()->SetTitle("BDT variable");
  hMVAValue[0]->GetYaxis()->SetRangeUser(0, 250);
  FormatHistLong(hMVAValue[0], 2, 1);
  hMVAValue[0]->Draw();
  FormatHist(hMVAValue[1], 4, 1);
  hMVAValue[1]->Draw("same");
  FormatHist(hMVAValue[2], 6, 1);
  hMVAValue[2]->Draw("same");
  FormatHist(hMVAValue[3], 8, 1);
  hMVAValue[3]->Draw("same");
  FormatHist(hMVAValue[4], 1, 1);
  hMVAValue[4]->Draw("same");

  TLegend* lM = new TLegend(0.54, 0.57, 0.84, 0.87);
  lM->AddEntry(hMVAValue[0], "Muon weights", "L");
  lM->AddEntry(hMVAValue[1], "Electron weights", "L");
  lM->AddEntry(hMVAValue[2], "Proton weights", "L");
  lM->AddEntry(hMVAValue[3], "#pi^{#pm} weights", "L");
  lM->AddEntry(hMVAValue[4], "Photon weights", "L");
  lM->SetBorderSize(0);
  lM->Draw();

  int likeVectorSize = likeRatio[0].size();
  for (int p = 1; p < nParticles; p++) {
    if (likeRatio[p].size() != likeVectorSize) {
      cerr << "Likelihood ratio vector sizes not equal" << endl;
      exit(1);
    }
  }

  int vectorIndex = 0;

  double particleLikeRatio[nParticles] = {0};
  double totalLikeRatio = 0.0;

  TH1D* hParticlePLR[nParticles] = {0};

  for (int p = 0; p < nParticles; p++)
    hParticlePLR[p] = new TH1D(Form("hParticlePLR%d", p), "", nProbBins, 0.0, 1.0);

  TH1D* hLikeRatio[nParticles] = {0};

  for (int p = 0; p < nParticles; p++)
    hLikeRatio[p] = new TH1D(Form("hLikeRatio%d", p), "", 40, -10.0, 10.0);

  while (vectorIndex < likeRatio[0].size()) {
    totalLikeRatio = 0.0;

    for (int p = 0; p < nParticles; p++) {
      likeIter = likeRatio[p].begin() + vectorIndex;
      particleLikeRatio[p] = (*likeIter);
      totalLikeRatio += (*likeIter);
    }

    for (int p = 0; p < nParticles; p++) {
      //if particleLikeRatio[p] / totalLikeRatio = 1.0, this goes to overflow bin of hParticlePLR[p]
      //for this reason, set it to be no more than 0.9999
      hParticlePLR[p]->Fill(TMath::Min(0.9999, particleLikeRatio[p] / totalLikeRatio));
    }

    ++vectorIndex;
  }

  TCanvas* ParticlePLR = new TCanvas("ParticlePLR", "ParticlePLR", 800, 600);
  FormatCanvas(ParticlePLR);
  hParticlePLR[0]->GetXaxis()->SetTitle("Proportion of likelihood ratio");
  hParticlePLR[0]->GetYaxis()->SetRangeUser(0, 800);
  FormatHistLong(hParticlePLR[0], 2, 1);
  hParticlePLR[0]->Draw();
  FormatHist(hParticlePLR[1], 4, 1);
  hParticlePLR[1]->Draw("same");
  FormatHist(hParticlePLR[2], 6, 1);
  hParticlePLR[2]->Draw("same");
  FormatHist(hParticlePLR[3], 8, 1);
  hParticlePLR[3]->Draw("same");
  FormatHist(hParticlePLR[4], 1, 1);
  hParticlePLR[4]->Draw("same");

  TLegend* lP = new TLegend(0.6, 0.55, 0.83, 0.85);
  lP->AddEntry(hParticlePLR[0], "Muon", "L");
  lP->AddEntry(hParticlePLR[1], "Electron", "L");
  lP->AddEntry(hParticlePLR[2], "Proton", "L");
  lP->AddEntry(hParticlePLR[3], "#pi^{#pm}", "L");
  lP->AddEntry(hParticlePLR[4], "Photon", "L");
  lP->SetBorderSize(0);
  lP->Draw();

} //end of CalcPropLikeRatio()

void FormatCanvas(TCanvas* canvas)
{
  canvas->Divide(1);
  canvas->cd(1);
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.15);
}

void FormatHistLong(TH1D* hist, int colour, int style)
{
  hist->SetTitle("");
  hist->SetLineColor(colour);
  hist->SetLineWidth(2);
  hist->SetLineStyle(style);
  hist->SetStats(0);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
}

void FormatHist(TH1D* hist, int colour, int style)
{
  hist->SetTitle("");
  hist->SetLineColor(colour);
  hist->SetLineWidth(2);
  hist->SetLineStyle(style);
  hist->SetStats(0);
}
