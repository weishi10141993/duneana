void PlotPhotonsYZ(TString file_name, int event)
{
  TFile * file = new TFile(file_name,"READ");
  TTree * tree = (TTree *)file->Get("opanalyzer/PhotonsPerOpDet");

  // Cuts to select only PMTs and particular event
 // TCut c1 = "isPMT";
  TCut c2 = Form("EventID == %d", event);
  
  //2D Histogram with SBND PDS dimensions
  TH2F * hPhotocatode = new TH2F("hPhotocathode", ";Z [cm];Y [cm]; #PE",48, 0, 1390, 40, -640, 640);
  hPhotocatode->SetStats(0);
  tree->Draw("OpDetY:OpDetZ >> hPhotocathode", (c2)*"CountAll", "COLZ");
}
