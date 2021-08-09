library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
#use Severe dataset - severe.rds / Moderate dataset -moderate.rds / Control dataset - control.rds
severe<-readRDS("severe.rds")

#comorbid genes list
cytokines<-c("ADIPOQ", "ADM", "ADM2", "AGRP", "AGT", "AMBN", "AMELX", "AMH", "ANGPTL5", "ANGPTL7", "APLN", "AREG", "ARMET", "ARMETL1", "ARTN", "AVP", "AZU1", "BDNF", "BMP1", "BMP10", "BMP15", "BMP2", "BMP3", "BMP4", "BMP5", "BMP6", "BMP7", "BMP8A", "BMP8B", "BTC", "C19orf10", "C3", "C5", "CALCA", "CALCB", "CAMP", "CAT", "CCK", "CCL1", "CCL11", "CCL13", "CCL14", "CCL14-CCL15", "CCL15", "CCL16", "CCL17", "CCL18", "CCL19", "CCL2", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", "CCL3", "CCL3L1", "CCL3L2", "CCL3L3", "CCL4", "CCL4L1", "CCL4L2", "CCL5", "CCL7", "CCL8", "CD320", "CD40LG", "CD70", "CECR1", "CER1", "CGA", "CGB", "CGB1", "CGB2", "CGB5", "CGB7", "CGB8", "CHGA", "CHGB", "CKLF", "CLCF1", "CLEC11A", "CMA1", "CMTM1", "CMTM2", "CMTM3", "CMTM4", "CMTM5", "CMTM6", "CMTM7", "CMTM8", "CNTF", "CORT", "CRH", "CSF1", "CSF2", "CSF3", "CSH1", "CSH2", "CSHL1", "CSPG5", "CTF1", "CTGF", "CTSG", "CX3CL1", "CXCL1", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL14", "CXCL16", "CXCL17", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL9", "CYR61", "DEFA1", "DEFA3", "DEFA5", "DEFB1", "DEFB103A", "DEFB104A", "DEFB4", "DKK1", "EBI3", "EDN1", "EDN2", "EDN3", "EGF", "EPGN", "EPO", "EREG", "ESM1", "FAM3B", "FAM3C", "FAM3D", "FASLG", "FGF1", "FGF10", "FGF11", "FGF12", "FGF13", "FGF14", "FGF16", "FGF17", "FGF18", "FGF19", "FGF2", "FGF20", "FGF21", "FGF22", "FGF23", "FGF3", "FGF4", "FGF5", "FGF6", "FGF7", "FGF8", "FGF9", "FIGF", "FIGNL2", "FLT3LG", "FSHB", "GAL", "GALP", "GAST", "GCG", "GDF1", "GDF10", "GDF11", "GDF15", "GDF2", "GDF3", "GDF5", "GDF6", "GDF7", "GDF9", "GDNF", "GH1", "GH2", "GHRH", "GHRL", "GIP", "GKN1", "GMFB", "GMFG", "GNRH1", "GNRH2", "GPHA2", "GPHB5", "GPI", "GREM1", "GREM2", "GRN", "GRP", "GUCA2A", "HAMP", "HBEGF", "HDGF", "HDGFRP3", "HGF", "HTN3", "IAPP", "IFNA1", "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA2", "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNB1", "IFNE", "IFNG", "IFNK", "IFNW1", "IGF1", "IGF2", "IL10", "IL11", "IL12A", "IL12B", "IL13", "IL15", "IL16", "IL17A", "IL17B", "IL17C", "IL17D", "IL17F", "IL18", "IL19", "IL1A", "IL1B", "IL1F10", "IL1F5", "IL1F6", "IL1F7", "IL1F8", "IL1F9", "IL1RN", "IL2", "IL20", "IL21", "IL22", "IL23A", "IL24", "IL25", "IL26", "IL27", "IL28A", "IL28B", "IL29", "IL3", "IL31", "IL32", "IL33", "IL34", "IL4", "IL5", "IL6", "IL6ST", "IL7", "IL8", "IL9", "INHA", "INHBA", "INHBB", "INHBC", "INHBE", "INS", "INS-IGF2", "INSL3", "INSL4", "INSL5", "INSL6", "JAG1", "JAG2", "KGFLP1", "KGFLP2", "KITLG", "KL", "LACRT", "LECT2", "LEFTY1", "LEFTY2", "LEP", "LHB", "LIF", "LRSAM1", "LTA", "LTB", "LTBP1", "LTBP2", "LTBP3", "LTBP4", "MDK", "MIA", "MIF", "MLN", "MSTN", "NAMPT", "NDP", "NENF", "NGF", "NMB", "NODAL", "NOV", "NPFF", "NPPA", "NPPB", "NPPC", "NPY", "NRG1", "NRG2", "NRG3", "NRG4", "NRTN", "NTF3", "NTF4", "NTS", "NUDT6", "OGN", "OSGIN1", "OSM", "OSTN", "OXT", "P11", "PDGFA", "PDGFB", "PDGFC", "PDGFD", "PDGFRA", "PDGFRB", "PDGFRL", "PDYN", "PENK", "PF4", "PF4V1", "PGF", "PLAU", "PMCH", "PNOC", "POMC", "PPBP", "PPBPL1", "PPBPL2", "PPY", "PRL", "PRLH", "PROK1", "PROK2", "PSPN", "PTH", "PTH2", "PTHLH", "PTN", "PYY", "QRFP", "RABEP1", "RABEP2", "REG1A", "RETN", "RETNLB", "RLN1", "RLN2", "RLN3", "RNASE2", "S100A6", "SAA1", "SAA2", "SBDS", "SCG2", "SCGB3A1", "SCT", "SCYE1", "SECTM1", "SEMA3A", "SEMA3B", "SEMA3C", "SEMA3D", "SEMA3E", "SEMA3F", "SEMA3G", "SEMA4A", "SEMA4B", "SEMA4C", "SEMA4D", "SEMA4F", "SEMA4G", "SEMA5A", "SEMA5B", "SEMA6A", "SEMA6B", "SEMA6C", "SEMA6D", "SEMA7A", "SLIT1", "SLIT2", "SLURP1", "SPP1", "SST", "STC1", "STC2", "TAC1", "TDGF1", "TDGF3", "TG", "TGFA", "TGFB1", "TGFB2", "TGFB3", "THPO", "TNC", "TNF", "TNFRSF11B", "TNFSF10", "TNFSF11", "TNFSF12", "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", "TNFSF18", "TNFSF4", "TNFSF8", "TNFSF9", "TOR2A", "TRH", "TSHB", "TSLP", "TXLNA", "TYMP", "UCN", "UCN2", "UCN3", "UTS2", "UTS2D", "VEGFA", "VEGFB", "VEGFC", "VGF", "VIP", "XCL1", "XCL2")
cyto_recep<-c("ACVR1B", "ACVR1C", "ACVR2A", "ACVR2B", "ACVRL1", "ADCYAP1R1", "ADIPOR1", "ADIPOR2", "ADRB1", "ADRB2", "AGTR1", "AGTR2", "AMHR2", "ANGPT1", "ANGPT4", "ANGPTL1", "ANGPTL2", "ANGPTL3", "ANGPTL4", "ANGPTL6", "APLNR", "AR", "AVPR1A", "AVPR1B", "AVPR2", "BMPR1A", "BMPR1B", "BMPR2", "BRD8", "C3AR1", "C5AR1", "CALCR", "CALCRL", "CCBP2", "CCR1", "CCR10", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCR8", "CCR9", "CCRL1", "CCRL2", "CD40", "CMKLR1", "CNTFR", "CRHR1", "CRHR2", "CRIM1", "CRLF1", "CRLF2", "CRLF3", "CSF1R", "CSF2RA", "CSF2RB", "CSF3R", "CX3CR1", "CXCR3", "CXCR4", "CXCR5", "CXCR6", "CXCR7", "CYSLTR1", "CYSLTR2", "DARC", "EDNRA", "EDNRB", "EGFR", "ENG", "EPOR", "ESR1", "ESR2", "ESRRA", "ESRRB", "ESRRG", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FGFRL1", "FLT1", "FLT3", "FLT4", "FPR1", "FPR2", "FPR2", "FSHR", "GALR2", "GALR3", "GCGR", "GHR", "GHRHR", "GHSR", "GIPR", "GLP1R", "GLP2R", "GNRHR", "GPER", "GPR17", "GPR32", "GPR33", "GPR44", "GPR77", "HNF4A", "HNF4G", "HTR3A", "HTR3B", "HTR3C", "HTR3D", "HTR3E", "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "IGF1R", "IGF2R", "IL10RA", "IL10RB", "IL11RA", "IL11RB", "IL12RB1", "IL12RB2", "IL13RA1", "IL13RA2", "IL15RA", "IL15RB", "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE", "IL18R1", "IL18RAP", "IL1R1", "IL1R2", "IL1RAP", "IL1RL1", "IL1RL2", "IL20RA", "IL20RB", "IL21R", "IL22RA1", "IL22RA2", "IL23R", "IL27RA", "IL28RA", "IL2RA", "IL2RB", "IL2RG", "IL31RA", "IL3RA", "IL4R", "IL5RA", "IL6R", "IL7R", "IL8RA", "IL8RB", "IL9R", "INSR", "KDR", "LEPR", "LGR4", "LGR5", "LGR6", "LHCGR", "LIFR", "LTB4R", "LTB4R2", "LTBR", "MC1R", "MC2R", "MC3R", "MC4R", "MCHR1", "MCHR2", "MET", "MLNR", "MPL", "MTNR1A", "MTNR1B", "NGFR", "NMBR", "NPR1", "NPR3", "NR0B1", "NR0B2", "NR1D1", "NR1D2", "NR1H2", "NR1H3", "NR1H4", "NR1I2", "NR1I3", "NR2C1", "NR2C2", "NR2E1", "NR2E3", "NR2F1", "NR2F2", "NR2F6", "NR3C1", "NR3C2", "NR4A1", "NR4A2", "NR4A3", "NR5A1", "NR5A2", "NR6A1", "NRP1", "NRP2", "OGFR", "OPRD1", "OPRK1", "OPRL1", "OPRM1", "OSMR", "OXTR", "PGR", "PGRMC2", "PLAUR", "PLXNA1", "PLXNA2", "PLXNA3", "PLXNA4", "PLXNB1", "PLXNB2", "PLXNB3", "PLXNC1", "PLXND1", "PPARA", "PPARD", "PPARG", "PRLHR", "PRLR", "PTAFR", "PTGDR", "PTGDS", "PTGER1", "PTGER2", "PTGER3", "PTGER4", "PTGFR", "PTH1R", "PTH2R", "RARA", "RARB", "RARG", "ROBO1", "ROBO2", "ROBO3", "RORA", "RORB", "RORC", "RXFP1", "RXFP2", "RXFP3", "RXRA", "RXRB", "RXRG", "S1PR1", "S1PR2", "SCTR", "SDC1", "SDC2", "SDC3", "SDC4", "SORT1", "SSTR1", "SSTR2", "SSTR5", "ST2", "TACR1", "TEK", "TGFBR1", "TGFBR2", "TGFBR3", "THRA", "THRB", "TIE1", "TNFRSF10A", "TNFRSF10B", "TNFRSF10C", "TNFRSF10D", "TNFRSF11A", "TNFRSF12A", "TNFRSF13B", "TNFRSF13C", "TNFRSF14", "TNFRSF17", "TNFRSF18", "TNFRSF19", "TNFRSF1A", "TNFRSF1B", "TNFRSF21", "TNFRSF25", "TNFRSF4", "TNFRSF6B", "TNFRSF8", "TNFRSF9", "TRHR", "TSHR", "TUBB3", "VDR", "VIPR1", "VIPR2", "XCR1")
rare_inf<-c("APOL1", "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9", "CARD9", "CCR5", "DARC", "EVER1", "EVER2", "CFD", "FUT2", "IFNGR1", "IFNGR2", "IL12B", "IL12RB1", "IL17F", "IL17RA", "IRAK4", "IRF7", "NEMO", "OX40", "PACRG", "PARK2", "PFC", "CFP", "SAP", "SH2D1A", "STAT1", "TMC6", "TMC8", "UNC93B1", "XIAP")
syndromic<-c("AHDC1", "ALG6", "ASH1L", "ASPM", "ATP1A1", "CSDE1", "DISC1", "DPYD", "HNRNPU", "HSD11B1", "KCND3", "KCNJ10", "KDM5B", "KIF14", "MIR137", "MTF1", "MTHFR", "MTOR", "NTNG1", "NTRK1", "POGZ", "POMGNT1", "RERE", "SLC45A1", "WDR26", "BCL11A", "CTNNA2", "CYP27A1", "DNMT3A", "FBXO11", "HDAC4", "HECW2", "KIF5C", "LNPK", "LRP2", "MBD5", "MYT1L", "NRXN1", "POU3F3", "SATB2", "SCN1A", "SCN2A", "TRIP12", "TTN", "ACY1", "CACNA1D", "CNTN4", "CTNNB1", "DHX30", "FOXP1", "MED12L", "PCCB", "SETD2", "SETD5", "SLC6A1", "STAG1", "TBC1D23", "TBL1XR1", "XPC", "ZBTB20", "CEP135", "GRID2", "NAA15", "NR3C2", "TBCK", "WDFY3", "CAMK2A", "KDM3B", "MEF2C", "NIPBL", "NR2F1", "NSD1", "PPP2CA", "TRIO", "ZSWIM6", "AHI1", "ALDH5A1", "ARID1B", "DLL1", "HIVEP2", "JARID2", "PHIP", "PPP2R5D", "SNX14", "SYNGAP1", "WASF1", "ZNF292", "ACTB", "ACTL6B", "AUTS2", "BRAF", "CAMK2B", "CDK13", "CEP41", "CNTNAP2", "CTTNBP2", "HOXA1", "INTS1", "KMT2C", "KMT2E", "RAC1", "RALA", "RELN", "RHEB", "TRRAP", "CHD7", "CLN8", "CYP11B1", "KAT6A", "PRKDC", "RAD21", "TRAPPC9", "VPS13B", "CACNA1B", "DOLK", "EHMT1", "GABBR2", "NFIB", "NTNG2", "NTRK2", "PRUNE2", "RORB", "SMARCA2", "STXBP1", "TSC1", "ZNF462", "EBF3", "KCNMA1", "PCDH15", "PTEN", "SMC3", "WAC", "ZMYND11", "BRSK2", "DEAF1", "DHCR7", "HEPACAM", "HRAS", "KMT2A", "KMT5B", "PACS1", "PAK1", "PAX6", "PHF21A", "SHANK2", "ANKS1B", "ARID2", "C12orf57", "CACNA1C", "CEP290", "CUX2", "GRIN2B", "MED13L", "PAH", "PRICKLE1", "PTPN11", "SCN8A", "SETD1B", "SMARCC2", "SOX5", "SYT1", "CDK8", "CHAMP1", "NBEA", "PCCA", "SETDB2", "CCNK", "CDC42BPB", "CHD8", "DYNC1H1", "ESR2", "FOXG1", "IRF2BPL", "PACS2", "PRKD1", "TRAPPC6B", "YY1", "ALDH1A3", "ARNT2", "BBS4", "CHD2", "GATM", "MAGEL2", "MEIS2", "NTRK3", "RORA", "SIN3A", "UBE3A", "ANKRD11", "CHMP1A", "CREBBP", "CTCF", "GRIN2A", "TRAF7", "TSC2", "USP7", "WWOX", "CACNA1G", "CHD3", "DLG4", "KANSL1", "KDM6B", "MED13", "NF1", "PPM1D", "PSMD12", "RAI1", "RNF135", "SGSH", "TANC2", "TAOK1", "TLK2", "VAMP2", "ASXL3", "SETBP1", "TCF4", "ATP1A3", "BRD4", "CNOT3", "DMPK", "KPTN", "MBOAT7", "NACC1", "NFIX", "PIK3R2", "PRR12", "SMARCA4", "UNC13A", "ADNP", "EEF1A2", "KCNB1", "KCNQ2", "DYRK1A", "SIK1", "SON", "ADSL", "CHKB", "DEPDC5", "EP300", "LZTR1", "PRODH", "SHANK3", "TBX1", "TCF20", "AFF2", "AP1S2", "ARHGEF9", "ARX", "ATRX", "BCORL1", "BRWD3", "CASK", "CDKL5", "CLCN4", "CNKSR2", "DDX3X", "DMD", "FMR1", "FRMPD4", "HCFC1", "HDAC8", "HNRNPH2", "HUWE1", "IQSEC2", "KDM5C", "KDM6A", "MAOA", "MECP2", "NEXMIF", "OCRL", "OPHN1", "PCDH19", "PHF8", "RLIM", "RPL10", "RPS6KA3", "SLC6A8", "SLC9A6", "SMC1A", "TAF1", "UPF3B", "USP9X")
copd<-c("ARHGAP42", "NPNT", "ADGRG6", "HTR4", "THSD4", "CHRNA5", "SFTPD", "SCGB1A1", "CHRNA3", "HTR4", "THSD4", "CDC123", "HCG20", "CFDP1", "LINC01807", "ADGRG6", "TET2", "ADGRG6", "FAM13A", "CASZ1", "C1GALT1", "AGER", "THSD4", "ADAM19", "ZDHHC20P2", "CHRNA3", "CHRNB4", "THSD4", "CHRNA3", "CHRNA5", "IP6K3", "LEMD2", "ARMC2", "EEFSEC", "FAM13A", "FAM13A", "FAM13A", "AHNAK", "LINC01679", "LINC00322", "FAM13A", "RIN3", "ID4", "CHRNA3", "HTR4", "TET2", "EEFSEC", "HERC1", "TMEM254", "CDH11", "HMGA2", "TNS1", "LINC01940", "HDAC4", "ZDHHC18", "PLB1", "FAM13A", "MTCL1", "HSF2BP", "HLA-DQA1", "HLA-DQB1", "BMP8A", "PPIEL", "ZFPM2", "CFDP1", "RARB", "LINC00886", "SLC30A10", "GTF2I", "ADAMTSL3", "ADAM19", "C1DP2", "C1DP3", "HDAC7", "CACNA2D3", "HYKK", "NPM1P35", "ANXA11", "FAM13A", "IREB2", "COL15A1", "FAM13A", "HLA-B", "STN1", "SYNPO2L", "PPARGC1A", "DHX15", "C1orf143", "TGFB2", "ASAP2", "VGLL4", "ARHGEF38", "DENND2D", "LRMDA", "FGFR3P1", "ZDHHC20P2", "ADCY5", "MFHAS1", "MICAL3", "PABPC4", "HEYL", "FAM13A", "AGER", "ITGA1", "SPATA9", "FAM13A", "CHRNA5", "SGF29", "MFAP2", "ZKSCAN1", "ARMC2", "ADAM19", "DOCK1", "EEFSEC", "ADGRG6", "AHNAK", "STAU1", "DDX27", "RSRC1", "RPL10P16", "ZNF490", "FGF18", "SPPL2C", "MAPT-AS1", "CASC15", "FAM13A", "PKD2L1", "TNS1", "LINC01807", "DMWD", "CASC17", "RNU7-155P", "RAB4B", "MIA-RAB4B", "RAB4B-EGLN2", "MMP12", "MMP3", "BTC", "DSP", "AMZ1", "C1orf143", "LINC02778", "LINC01748", "SFTPD", "RIN3", "CDH13", "PITPNA", "RIN3", "PSORS1C1", "MTCL1", "APIP", "TMEM219", "PABPC4", "HEYL", "CHRM3", "CHRM3-AS2", "ITGB8", "PSORS1C1", "TGFB2", "THRA", "SYN3", "DTWD1", "RFX6", "RPS29P13", "BTBD1", "LINC01804", "CDH13", "TNRC6A", "RSPH6A", "ZNF319", "TEPP", "CCDC69", "GM2A", "RREB1", "EFCAB5", "MIR4734", "EPOP", "CHRNA3", "IREB2", "IREB2", "LINC02252", "PDZD2", "COX10-AS1", "CDRT15P1", "CHRNA3", "SFTPD", "RARB", "LINC01876", "RBMS3", "SLMAP", "STN1", "HSPA4", "SFTPD", "FGD6", "ITPK1", "CYS1", "FUT8", "Metazoa_SRP", "TNPO1", "RNU6-234P", "ASB1", "ASTN2", "CCDC91", "FTO", "GPR65", "RNU6-835P", "CNTN4", "EEFSEC", "RSL24D1P6", "TSPY26P", "CFDP1", "TESK2", "TSPAN14", "Y_RNA", "MMP15", "CRACR2B", "IREB2", "BNIP3P1", "SGCD", "SLC35F3", "ZBTB38", "EML4", "EMP2", "TEKT5", "NUP42", "KLHL7", "PCDH15", "CAPRIN2", "KCNMA1", "SIX4", "SIX1", "RAB4B", "MIA-RAB4B", "RAB4B-EGLN2", "ZBTB9", "GGNBP1", "LINC02812", "SLC25A3P1", "RIN3", "SUZ12P1", "ADAM20P3", "DLG2", "LINC02252", "PLXNA4", "PCDH9", "CEP70", "USP24", "HSBP1")
cardio<-c("FGF5", "PRDM8", "APOE", "NOS3", "CABCOCO1", "INSR", "RGL3", "ARHGAP42", "NPR3", "KCNK3", "ZNF831", "CLCN6", "ATXN2", "SH2B3", "FES", "ZPR1", "ATP2B1", "FTO", "LINC02227", "PLCB1", "PPIAP55", "LDLR", "SMARCA4", "CASZ1", "CYP1A2", "CSK", "CELSR2", "LSP1", "ADRB1", "LPA", "PLCE1", "PDILT", "PLEKHG1", "CACNB2", "GUCY1A1", "ABCG8", "NT5C2", "SWAP70", "MYO9B", "RNU1-96P", "SOS2", "OPRL1", "SORCS3", "ADO", "NPNT", "HOXC4", "MARK2P11", "HOTTIP", "GOSR2", "AGT", "BANK1", "ABHD17C", "CDC25A", "WNT2B", "ADH1B", "CMIP", "ZFPM2", "MECOM", "INSR", "IGF2BP1", "B4GALNT2P1", "DPEP1", "GCKR", "C20orf181", "NAA38", "ZDHHC18", "ZC3HC1", "CHDH", "TSC22D2", "SLC9B1", "UNC13D", "PRKAG2", "CACNA1D", "ENPEP", "PRSS36", "PRSS8", "LINC01571", "C16orf97", "VEGFA", "USP38", "RXFP2", "TXNL4B", "HPR", "KLF14", "BCAS3", "HSPA4", "HGFAC", "FBN2", "LINC01679", "LINC00322", "PCSK9", "SOX6", "NFKBIA", "PLCE1", "PKHD1", "SVEP1", "SEPTIN9", "SCAF1", "HSD17B12", "THADA", "COL4A1", "LINC01478", "MAD1L1", "RORA", "ALG9", "LINC01169", "LFNG", "GRIFIN", "HIPK2", "FGF9", "RN7SL766P", "AXIN1", "MECOM", "EBF1", "MXD3", "HNF4G", "HNF1A-AS1", "MACROD1", "LIMA1", "JPH2", "ATG7", "PALM2-AKAP2", "LCORL", "CLUAP1", "NPC1", "NAV1", "IPO9-AS1", "PDE10A", "FAM133B", "CDK6", "CACNA1D", "MRAS", "TMEM72-AS1", "DBH-AS1", "DBH", "TBC1D19", "NPR3", "ZNF268", "ADCY9", "UBN1", "USP36", "KANK2", "VDR", "OR10AD1", "PREX1", "TERF2", "CFDP1", "LINC02513", "SLC7A1", "REXO1", "RUVBL1", "RNF43", "TSPOAP1-AS1", "MAFB", "NPR1", "INPP5A", "HNRNPFP1", "SLC9A3R2", "DDX23", "WNT3A", "SEPTIN14P17", "LINC02235", "NCALD", "PAQR5", "C1QTNF4", "GLYAT", "LINC02468", "USP34", "LINC02245", "LINC02576", "HIVEP2", "Y_RNA", "PRDM6", "C5orf56", "NDUFAF6", "FREM1", "CPEB2", "HOXB3", "HOXB-AS3", "HOXB5", "UBASH3B", "CCND2", "CCND2-AS1", "DNAJC27-AS1", "MTRF1", "RAC1P3", "MCF2L", "PCNX1", "PLPP3", "ZFAND2A", "RYK", "LINC02004", "LPL", "FGD5", "PAX2", "CERT1", "HMGCR", "RN7SL804P", "YEATS4", "UHRF1", "RTKN2", "ARID5B", "GML", "CYP11B2", "MECOM", "PKD2L1", "RPS6", "RPS23P2", "PCDH18", "UPF3A", "GMDS", "FOXC1", "HSF2BP", "ABCF3", "ARNTL", "DLEU1", "DLEU7", "BAHCC1", "RNU1-46P", "MEX3C", "NMT1", "DOCK7", "GTF2B", "RN7SL299P", "SIPA1L2", "GTF2I", "CASC15", "IRX1", "TANGO2", "PGPEP1", "SMARCE1", "TRIM48", "OR4A6P", "PHF13", "SRRM1", "PRKCE", "HBS1L", "MYB", "AKR1B1", "AKR1B10", "RPL23AP31", "PRKAG3", "CASC15", "ARL15", "LINC01344", "SPTBN5", "C22orf31", "MYEOV", "DCUN1D5", "PRDM16", "HS1BP3", "FOSL2", "HIBADH", "PNKD", "TMBIM1", "MXI1", "LINC01153", "RN7SKP167", "CDKAL1", "UBXN2B", "GIT2", "TCHP", "FBRSL1", "S1PR2", "DNMT1", "NFATC2", "WT1-AS", "LHFPL2", "TIMD4", "DPYSL2", "WRN", "FTH1P7", "MARK3", "CKB", "CDCA5", "DDX56", "NPC1L1", "WAPL", "YES1", "APOLD1", "FAM131C", "EPHA2", "CCDC148", "SORBS1", "ARMC4", "C20orf203", "NOL4L-DT", "RERE", "CCDC170", "CADPS2", "BICC1", "MPPED2-AS1", "CRTC1", "LINC01122", "ADK", "LINC01485", "CPEB4", "LINC00824", "ZCCHC7", "ATL1")
hyperten<-c("ZFAT", "APOE", "MTHFR", "CUX2", "PRDM8", "FGF5", "PRDM8", "FGF5", "NPR3", "CASZ1", "FGF5", "PRDM8", "WNT2B", "RN7SL222P", "PPIAP43", "CACNB2", "PRDM8", "FGF5", "PRDM8", "FGF5", "NT5C2", "NPR3", "PLCD3", "FES", "KCNK3", "WNT2B", "FGF5", "PRDM8", "CNNM2", "LINC02577", "PRPS1P1", "CYP17A1", "LSP1", "ARHGAP42", "RGL3", "ATP2B1", "PRDM8", "FGF5", "RGL3", "INSR", "ZNF831", "CYP1A2", "CSK", "RGL3", "FGF5", "PRDM8", "ZNF831", "KCNK3", "KCNK3", "ATXN2", "KCNK3", "CACNA1D", "FGD5", "MOV10", "PSMC3", "C17orf82", "NT5C2", "ST13P13", "ARHGAP42", "CNNM2", "CABCOCO1", "SBF2", "SLC39A8", "PHETA1", "CUX2", "FGF5", "NR1H3", "MTCH2", "FGD5", "RNU1-96P", "UMOD", "SEMA7A", "GPR39", "SLC44A4", "ULK4", "CAPZA1", "BAG6", "AMPD3", "EML6", "HFE", "MYO6", "NPR3", "BORCS7-ASMT", "BORCS7", "HOTTIP", "CERS5", "ARHGAP42", "HOTTIP", "FAM219B", "MPI", "CERS5", "HCP5", "USP8", "ARHGAP42", "NUP160", "BANK1", "LINC02571", "HLA-B", "HOTTIP", "PRDM16", "MTND1P22", "PTPRJ", "HOXB-AS3", "HOXB3", "CASZ1", "RPL35P4", "HNRNPA1P73", "HOXA-AS2", "HOXA3", "ARHGAP42", "ADRB1", "NHLRC2", "GUCY1B1", "GRM5", "CAPZA1", "NOS3", "FURIN", "KCNK3", "YES1", "RNU6-758P", "BCAR1", "CYP17A1", "LRRC10B", "PGPEP1", "L3MBTL4", "MACROD2", "HOXA11", "HOXA10", "OR4X1", "FOLH1", "FRMD3", "FURIN", "WDR92", "RPTOR", "CASZ1", "RPL35P4", "HNRNPA1P73", "SWAP70", "PLCE1", "FRMD3", "HOXA-AS2", "HOXA3", "LSP1", "HOXA11", "HOXA10", "LINC02163", "ATP2B1", "'-", "GPR20", "FGD4", "GML", "LSP1", "FGD5", "STN1", "LRRC7", "SLC12A9", "MYO9B", "CYP11B2", "GML", "ENPEP", "DLK1", "MEG3", "HMGCS2", "ATP2B1", "SWAP70", "LRRC4C", "PLEKHH2", "LINC00862", "CDKN2B-AS1", "CACNB2", "HECTD4", "NFATC2", "OR5B12", "UMOD", "PTPMT1", "NDUFS3", "LSP1", "FSTL4", "STK3", "SLC39A12", "CACNB2")
diab_mell<-c("HERPUD1", "CETP", "ZPR1", "GCKR", "PDILT", "UMOD", "THEMIS3P", "AKR1B1P6", "THEMIS3P", "AKR1B1P6", "CELSR2", "LPL", "LIPC-AS1", "LIPC", "ALDH1A2", "PPM1G", "NRBP1", "FADS1", "FADS2", "CUBN", "PRKAG2", "APOC1P1", "APOC1", "SEMA5B", "MPPED2-AS1", "APOB", "SEMA5B", "TRIB1", "BAZ1B", "TCF7L2", "TDRD15", "ITPK1", "ITPK1", "VPS37D", "MLXIPL", "FTO", "BCL3", "GCKR", "SHROOM3", "HBB", "LINC02570", "LINC00243", "RNU6-1150P", "AGPAT3", "PDILT", "LPA", "CPS1", "TRIM31-AS1", "TRIM31", "TSBP1-AS1", "HLA-DRA", "C2", "OR11A1", "UNCX", "ZFHX3", "SUGP1", "CDKN2B-AS1", "UMOD", "BSX", "RNU7-2P", "ZFHX3", "SLC30A3", "SIK3", "TCF7L2", "HCG14", "TRIM27", "PIP5K1B", "OR2B2", "NFATC1", "MX1", "TMPRSS2", "FAM238C", "ANKRD26", "CERT1", "PDILT", "UMOD", "HIST1H1PS2", "POLK", "PRKAG2", "TPRKB", "CARMIL1", "PLEKHH2", "MAML3", "BCAS3", "SLC39A8", "CDKN2B-AS1")
obesity<-c("FTO", "UGT1A4", "UGT1A1", "UGT1A8", "UGT1A6", "UGT1A7", "UGT1A5", "UGT1A9", "UGT1A10", "UGT1A3", "IL1RAP", "SAA2", "ST13P5", "AC093326.1", "TMEM18", "AL772161.2", "ABO", "AC090771.1", "RNU4-17P", "IL16", "AC108467.1", "THAP12P9", "MIA", "SNRPA", "CCL16", "AC011473.3", "KLK12", "THAP12P9", "AC108467.1", "RNU4-17P", "AC090771.1", "AC011511.5", "ICAM5", "PLA2G2A", "TMEM18", "AC093326.1", "REXO4", "AL136114.1", "AC103796.1", "BDNF", "TFAP2B", "CCL23", "CCL18", "ECM1", "AC131157.1", "BCDIN3D", "AL049693.1", "RPS17P5", "ADCY3", "DNAJC27", "NRXN3", "AL513166.1", "CA6", "BDNF-AS", "CFH", "AC008946.1", "CCL25", "TXNL4B", "HPR", "FCN2", "FAIM2", "DGKG", "GIPR", "SERPINA5", "SERPINA4", "TDGF1", "LRRC2", "AL591074.2", "SH2B1", "AC009159.4", "AC009159.1", "MAP2K5", "RNU4-17P", "AC090771.2", "POC5", "CPN1", "GPRC5B", "GPR139", "RPL5P33", "WFIKKN2", "GNAT2", "CNTN5", "ADCY3", "TRIM66", "TNNI3K", "FPGT-TNNI3K", "EXOC4 x OR10R1P - AL365440.1", "CFB", "AL645922.1", "CADM2", "AC062015.1", "MIR5702", "AZU1", "PRTN3", "AC090771.2", "NXPE1", "AC020549.3", "EXOC4 x OR10R1P", "AC007092.1", "AC005586.1", "AC005586.2", "EXOC4 x AL365440.1 - OR6Y1", "EXOC4 x AL365440.2", "OR10R3P", "GPR139", "GPRC5B", "LINGO2", "BCDIN3D-AS1", "KCNMA1", "PRKCH", "CUL4A", "RNU6-382P", "AC016925.2", "RNU6-1060P", "AL139231.1", "AC011029.1", "HNF4G", "AL139260.1", "MYCBP", "AL139260.2", "AL139260.3", "OR10R1P - AL365440.1 x EXOC4", "EXOC4 x OR6Y1", "LEPR", "PACS1", "TENM3-AS1", "GNPTAB", "EXOC4 x AL365440.2", "OR10R2", "AL356295.1", "LINC01122", "NMI", "TNFAIP6", "ETV5", "DGKG", "RMST", "TGM2", "HOXB5", "HOXB3", "HOXB-AS3", "POC5", "SLC25A5P9", "FPGT-TNNI3K", "TNNI3K", "C7", "RABEP2", "PAX5", "AC034111.2", "AC126283.2", "HS6ST3", "SLC41A2", "LINC00599", "TNKS", "RPTOR", "LINC01524", "SSR3", "TIPARP-AS1", "TDH", "AL445687.2", "RNU6-716P", "AL354741.1", "FAM155A", "AL356364.1", "LYPLAL1-AS1", "AC118549.1", "NCAM2", "TECTB", "ADCY9", "MRPS22", "GPR1", "GAS8", "URAHP", "AC011483.1", "KLK7", "LINC01615", "LINC02544", "AC008060.1", "HNRNPA1P31", "RNU6-61P", "AC124242.3", "CCDC77", "AC022690.1", "AC022690.2", "TCF7L2", "RNU7-165P", "RNU6-649P", "AC092169.1", "NPC1", "SLC35A1", "AL049697.1", "AL512603.2", "AL512603.1", "CBLIF", "MRPL16", "ANO5", "AC008115.2", "GPR19", "AL161910.1", "SEMA4D", "CHL1", "NEGR1", "TENT5A", "AL590824.1", "AC091857.1", "RPL36AP21", "INHBB", "LINC01101", "SLC8A1", "AC013262.1", "RN7SL201P", "PTDSS2", "TIPARP-AS1", "SSR3", "AP000894.1", "BOLA2P1", "LINC02654", "RNU6-1075P", "LINC02465", "RTP4", "CXCR5", "AP004609.1", "SUGP1", "PNPLA3", "ARG1", "SCNN1A", "CASC15", "TBCE", "FXN", "AL358113.1", "AC096639.1", "SOCS2-AS1", "SOX6", "OR2W5", "OR2B11", "NUP98", "FARS2", "SDCCAG8", "METTL15", "AL136962.1", "AL512484.1", "RAMP3", "AC073968.2", "LHFPL3", "TTC28", "TTC28-AS1", "PCDH9", "AC068688.1", "AC079052.1", "AC092047.1", "VIPR1", "WWOX", "AC004052.1", "LINC02503", "PTPRG", "FHIT", "MTRNR2L7", "AL132657.1")

#change vars=.... from above depending on the comorbidity being analyzed
#cluster 0-17 for severe; 0-16 for moderate; 0-18 for control 
sev0 <- as.matrix(FetchData(object=severe, vars = cytokines,cells=c(names(severe$seurat_clusters[severe$seurat_clusters==0]))))
values<-as.list(sev0)
h<-dim(sev0)[1]*dim(sev0)[2]
for(i in 1:h){
  nam[i]<-"cluster 0"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev0.csv")

sev1 <- as.matrix(FetchData(object=severe, vars = cytokines,cells=c(names(severe$seurat_clusters[severe$seurat_clusters==1]))))
values<-as.list(sev1)
h<-dim(sev1)[1]*dim(sev1)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 1"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev1.csv")


sev <- as.matrix(FetchData(object=severe, vars = cytokines,
  cells=c(names(severe$seurat_clusters[severe$seurat_clusters==2]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 2"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev2.csv")

sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==3]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 3"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev3.csv")

sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==4]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 4"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev4.csv")



sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==5]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 5"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev5.csv")



sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==6]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 6"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev6.csv")



sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==7]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 7"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev7.csv")



sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==8]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 8"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev8.csv")



sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==9]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 9"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev9.csv")



sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==10]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 10"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev10.csv")



sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==11]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 11"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev11.csv")




sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==12]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 12"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev12.csv")



sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==13]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 13"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev13.csv")



sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==14]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 14"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev14.csv")



sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==15]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 15"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev15.csv")




sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==16]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 16"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev16.csv")



sev <- as.matrix(FetchData(object=severe, vars = cytokines,
                           cells=c(names(severe$seurat_clusters[severe$seurat_clusters==17]))))
values<-as.list(sev)
h<-dim(sev)[1]*dim(sev)[2]
nam<-character(0)
for(i in 1:h){
  nam[i]<-"cluster 17"
}
df<-data.frame(expr=values)
df<-t(df)
write.csv(df,"sev17.csv")

#Combine sev 1 to sev 17 - sev.csv
#cluster expression
#cluster 1	1.977239278
#cluster 2	1.847175308
#cluster 15 1.372357886

g <- ggplot(expr, aes(x=sample, y=expression, fill=sample))
g + geom_boxplot()+geom_point()+ggtitle("Cytokines")


