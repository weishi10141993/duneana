void PlotOpHitYZ(TString file_name, int event)
{
  TFile * file = new TFile(file_name,"READ");
  TTree * tree = (TTree *)file->Get("hitdumper/hitdumpertree");

  TCut c1 = Form("event == %d", event);

  //2D Histogram with SBND PDS dimensions
  TH2F * hPhotocatode = new TH2F("hPhotocatode", ";Z [cm];Y [cm]; #PE", 48, 0, 1390, 40, -640, 640);
  hPhotocatode->SetStats(0);
  tree->Draw("ophit_opdet_y:ophit_opdet_z >> hPhotocatode", (c1)*"ophit_pe", "COLZ");
}
