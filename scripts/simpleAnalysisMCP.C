{
  //TFile *_file0 = TFile::Open("plots/moduleCharacterization_step1_HPK528_unirr_52deg_T10C_Vov3.50_MCP.root");
  TFile *_file0 = TFile::Open("plots/moduleCharacterization_step1_FNAL2023_example.root");
  TTree* d1=(TTree*)_file0->Get("data_bar09L-R_Vov3.50_th15");

  //energy definitions and calibrations. Need to make landau fits
  d1->SetAlias("tE","(energyL+energyR)");
  d1->SetAlias("eCalib","tE/400");
  d1->SetAlias("eLCalib","energyL/200");
  d1->SetAlias("eRCalib","energyR/200");

  //phase defintions for TOFHIR from calibrated time
  d1->SetAlias("phi_L","(timeL-floor(timeL/6250.)*6250)/6250.");
  d1->SetAlias("phi_R","(timeR-floor(timeR/6250.)*6250)/6250.");

  //time walk & position corrected (interplay/2D) to be fully understood
  d1->SetAlias("phi_L_C","phi_L-tP_L_9-tW_L_9-0.5");
  d1->SetAlias("phi_R_C","phi_R-tP_R_9-tW_R_9-0.5");
  d1->SetAlias("phase_TOFHIR","(phi_L+phi_R)/2.-0.5");
  d1->SetAlias("phase_TOFHIR_C","(phi_L_C+phi_R_C)/2.");

  //phase defintions for H4DAQ
  d1->SetAlias("phase_H4","(t_MCP-t_CLK_P-0.849*6.25-floor((t_MCP-t_CLK_P)/6.25)*6.25)/6.25");
  d1->SetAlias("phase_H4_C","phase_H4-round(phase_H4)");

  //deltaT or deltaPhase
  d1->SetAlias("delta_phase","phase_TOFHIR-phase_H4_C");
  d1->SetAlias("delta_phase_C","phase_TOFHIR_C-phase_H4_C");

  //position correction (C is for L+R average, L or R to check separate corrections. Need to be applied to evaluate tWalk if computing L+R corrections)
  d1->SetAlias("tP_C_9","-0.000553625*x");
  d1->SetAlias("tP_L_9","-0.0417938+0.0013515*x");
  d1->SetAlias("tP_R_9","0.00882593-0.0024042*x");

  //time walk corrections (C is L+R average, L or R to check separate corrections)
  d1->SetAlias("tW_C_9","0.0546111-0.0953301*eCalib+0.0289362*eCalib*eCalib-0.00386627*eCalib*eCalib*eCalib");
  d1->SetAlias("tW_L_9","0.0723849-0.0909822*eLCalib+0.0244774*eLCalib*eLCalib-0.00243648*eLCalib*eLCalib*eLCalib");
  d1->SetAlias("tW_R_9","0.0631986-0.0648019*eRCalib+0.0102749*eRCalib*eRCalib-4.6009e-05*eRCalib*eRCalib*eRCalib");

  //  d1->SetAlias("tW_C_9","0.0638698-6.50564e-05*tE+7.68871e-09*tE*tE");
  // d1->SetAlias("tW_C_9","-1.48564e-01+2.17896e-01*tE/750-1.08096e-01*tE*tE/750/750+1.98548e-02*tE*tE*tE/750/750/750");

  //standard accepance cuts: ~>0.8 MIP & !saturated, 500<amp_MCP<1000, distance from MCP <10mm 
  d1->SetAlias("cut","(energyL+energyR)>300. && amp_MCP>150 && abs(delta_phase)<0.2 && (energyL+energyR)<1800 && amp_MCP<1000 && energyL<900 && energyR<900");
  //looser cut to check dependence on MCP amplitude
  // d1->SetAlias("cutLooseMCP","(energyL+energyR)>700. && amp_MCP>200 && abs(delta_phase)<0.2 && (energyL+energyR)<1800 && amp_MCP<1000 && energyL<900 && energyR<900 && sqrt((x-7.1)*(x-7.1)+(y-8.8)*(y-8.8))<=10");

  //series of various plots to evaluate resolution and/or corrections
  //    d1->Draw("phase_TOFHIR:phase_H4_C>>hh(2000,-2,2,2000,-2,2)","(energyL+energyR)>700. && amp_MCP>500 && abs(delta_phase)<0.2 && (energyL+energyR)<1800","COLZ");
  d1->Draw("(delta_phase):tE/400>>hh(300,0,3,2000,-2,2)","cut","COLZ");
  // d1->Draw("(delta_phase)-tW_C_9:phase_H4_C>>hh1(80,-2,2,2000,-2,2)","(energyL+energyR)>700. && amp_MCP>500 && abs(delta_phase)<0.2 && (energyL+energyR)<1800","COLZ");
  d1->Draw("(delta_phase)-tW_C_9-tP_C_9:phase_TOFHIR>>hh1(80,-2,2,2000,-2,2)","cut","PROF");
  //d1->Draw("(delta_phase)-tW_C_9>>hh(2000,-2,2)","cut");
  //  d1->Draw("((delta_phase)-tW_C_9):x>>hh2(200,-50,50,2000,-2,2)","cut","COLZ");
  //d1->Draw("(delta_phase)-tW_C_9-tP_C_9:amp_MCP>>hh2(20,0,1000,2000,-2,2)","cutLooseMCP","COLZ");
  d1->Draw("(delta_phase)-tW_C_9-tP_C_9>>hh(2000,-2,2)","cut");
  //  d1->Draw("(delta_phase)-tW_C_9-tP_C_9:phase_TOFHIR>>hh1(80,-2,2,2000,-2,2)","cut","COLZ");
  //d1->Draw("(delta_phase)-tW_C_9-tP_C_9:phase_TOFHIR>>hh1(80,-2,2,2000,-2,2)","cut","COLZ");
  //  d1->Draw("phi_L-0.5-phase_H4_C-tP_L_9-tW_L_9:energyL/500:x>>hh3(50,-50,50,50,0.5,2.)","cut && energyL>400","PROF2DCOLZ");
  //  d1->Draw("phi_L-0.5-phase_H4_C-tP_L_9-tW_L_9:x>>hh3(100,-50,50,2000,-2,2)","cut && energyL>400","COLZ");
  //  d1->Draw("phi_L-0.5-phase_H4_C-tP_L_9-tW_L_9:energyL/500>>hh3(60,0.5,2,2000,-2,2)","cut && energyL>400","COLZ");
  //  d1->Draw("phi_L-0.5-phase_H4_C-tP_L_9-tW_L_9>>hh4(2000,-2,2)","cut && energyL>400","COLZ");
  //d1->Draw("phi_L-0.5-phase_H4_C-tP_L_9-tW_L_9:phase_TOFHIR>>hh4(80,-2,2,2000,-2,2)","cut && energyL>400","COLZ");
  //  d1->Draw("phi_R-0.5-phase_H4_C:x>>hh3(100,-50,50,2000,-2,2)","cut && energyR>300","COLZ");
  //  d1->Draw("phi_R-0.5-phase_H4_C-tP_R_9-tW_R_9:energyR/350>>hh3(60,0.5,2,2000,-2,2)","cut && energyR>300","COLZ");
  //  d1->Draw("phi_R-0.5-phase_H4_C-tP_R_9-tW_R_9>>hh4(2000,-2,2)","cut && energyR>300","");
  //  d1->Draw("phi_R-0.5-phase_H4_C-tP_R_9-tW_R_9:phase_TOFHIR>>hh4(80,-2,2,2000,-2,2)","cut && energyR>300","COLZ");
  //  d1->Draw("phi_L-0.5-phase_H4_C-tP_L_9-tW_L_9:phase_TOFHIR>>hh5(80,-2,2,2000,-2,2)","cut && energyL>400","COLZ");
  //  d1->Draw("(delta_phase_C):phase_TOFHIR>>hh6(80,-2,2,2000,-2,2)","cut","COLZ");
  //  d1->Draw("(delta_phase_C):amp_MCP>>hh2(20,0,1000,4000,-2,2)","cutLooseMCP","COLZ");
  //  d1->Draw("(delta_phase_C)>>hh2(4000,-2,2)","cut","COLZ");
  //  d1->Draw("(delta_phase_C)>>hh2(4000,-2,2)","cut && abs(phase_TOFHIR)<0.02","COLZ");
  // d1->Draw("(delta_phase_C):x>>hh2(200,-50,50,2000,-2,2)","cut","COLZ");
}
