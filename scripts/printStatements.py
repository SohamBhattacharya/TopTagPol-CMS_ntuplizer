#!bin/sh python

import os


outfileName = "statements.txt"

toReplace = "$$$"

l_str = [
    "sv_pT_reco",
    "sv_eta_reco",
    "sv_phi_reco",
    "sv_m_reco",
    "sv_E_reco",
    "sv_etarel_reco",
    "sv_phirel_reco",
    "sv_deltaR_reco",
    "sv_ntracks_reco",
    "sv_chi2_reco",
    "sv_ndf_reco",
    "sv_normchi2_reco",
    "sv_dxy_reco",
    "sv_dxyerr_reco",
    "sv_dxysig_reco",
    "sv_d3d_reco",
    "sv_d3derr_reco",
    "sv_d3dsig_reco",
    "sv_costhetasvpv_reco",
    "sv_enratio_reco",
]


#templateStr = (
#    "sprintf(brName, \"jet_%s_{name}\", str_jetName.c_str());\n"
#    "tree->Branch(brName, &vv_jet_{name});"
#)
#linegaps = 1

templateStr = (
    "jetInfo->vv_jet_{name}.push_back(v_jet_{name});"
)
linegaps = 0


with open(outfileName, "w") as f :
    
    for iEntry in range(0, len(l_str)) :
        
        #temp_str = templateStr %(l_str[iEntry])
        temp_str = templateStr
        temp_str = temp_str.format(name = l_str[iEntry])
        
        print(temp_str + "\n" * linegaps)
        
        temp_str += "\n"
        temp_str += "\n" * linegaps
        
        f.write(temp_str)
        

