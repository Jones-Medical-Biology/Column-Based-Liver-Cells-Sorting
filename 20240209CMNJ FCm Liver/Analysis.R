library(flowCore)
library(flowViz)
library(flowWorkspace)
library(ggcyto)

# get_fcs_resultsQC(QC_folder = "resultsQC/", 
#                   raw_folder = "Exp_20240209_1/") -> 
#   cs
setwd("/Users/brianschildt/OneDrive - University of Florida/Notebook Hoffman Brian/Projects/Column Based Liver Cells Sorting/20240209CMNJ FCm Liver")
fcsfiles <- dir("resultsQC", pattern = "*.fcs",full.names = TRUE) # list contents but only in directory with .fcs
fcsfiles.FS <- read.flowSet(fcsfiles)
flowPlot(fcsfiles.FS[[1]], plotParameters = c("FL8-A","SSC-A"))
fcsfiles.FS.GS <- GatingSet(fcsfiles.FS)
gs_get_pop_paths(fcsfiles.FS.GS)
rg_scatter <- rectangleGate("FSC-A"=c(200,4e6), "SSC-A"=c(250, 5e6), filterId="scatter")
gs_pop_add(fcsfiles.FS.GS, rg_scatter,parent = "root") #Good plot

rg_singlet <- rectangleGate("FSC-Width"=c(200,4e3), "FSC-H"=c(250, 1e6), filterId="singlet")
gs_pop_add(fcsfiles.FS.GS, rg_singlet, parent = "scatter") #Good plot

rg_WBC <- rectangleGate("FL8-A"=c(4e3,1e6), "SSC-A"=c(1e4, 5e6), filterId="CD45+")
gs_pop_add(fcsfiles.FS.GS, rg_WBC, parent = "singlet") #Good plot

rg_T <- rectangleGate("FL5-A"=c(8e3,1e6), "SSC-A"=c(1e4, 1e6), filterId="Tcells")
gs_pop_add(fcsfiles.FS.GS, rg_T, parent = "CD45+")

rg_B <- rectangleGate("FL19-A"=c(1e4,5e5), "SSC-A"=c(1e5, 5e6), filterId="Bcells")
gs_pop_add(fcsfiles.FS.GS, rg_B, parent = "CD45+")

rg_NK <- rectangleGate("FL13-A"=c(1e4,1e6), "SSC-A"=c(1e4, 5e6), filterId="NKcells")
gs_pop_add(fcsfiles.FS.GS, rg_NK, parent = "CD45+")

rg_NT <- rectangleGate("FL5-A"=c(1e2,6e3), "SSC-A"=c(1e4, 5e5), filterId="notT")
gs_pop_add(fcsfiles.FS.GS, rg_NT, parent = "CD45+")

rg_NB <- rectangleGate("FL19-A"=c(0,2e3), "SSC-A"=c(1e4, 5e5), filterId="notB")
gs_pop_add(fcsfiles.FS.GS, rg_NB, parent = "notT")

rg_K <- rectangleGate("FL13-A"=c(-3e2,0), "SSC-A"=c(1e4, 4e5), filterId="kupffer")
gs_pop_add(fcsfiles.FS.GS, rg_K, parent = "notB")

rg_Vsig <- rectangleGate("FL1-A"=c(4e2,2e3), "SSC-A"=c(1e4, 3e5), filterId="vsig4+")
gs_pop_add(fcsfiles.FS.GS, rg_Vsig, parent = "kupffer")

rg_VsigNK <- rectangleGate("FL1-A"=c(-5e2,3e2), "SSC-A"=c(1e4, 3e5), filterId="vsig4-")
gs_pop_add(fcsfiles.FS.GS, rg_VsigNK, parent = "kupffer")

rg_NBC <- rectangleGate("FL8-A"=c(-1000,3e3), "SSC-A"=c(1e3, 1e6), filterId="CD45-")
gs_pop_add(fcsfiles.FS.GS, rg_NBC, parent = "singlet")

rg_OC <- rectangleGate("FL3-A"=c(1e3,5e4), "SSC-A"=c(1e3, 1e6), filterId="CD146+")
gs_pop_add(fcsfiles.FS.GS, rg_OC, parent = "CD45-")

rg_lsec <- rectangleGate("FL17-A"=c(0,5e4), "SSC-A"=c(1e2, 5e4), filterId="lsec")
gs_pop_add(fcsfiles.FS.GS, rg_lsec, parent = "CD146+") #Noah will address later scales need to be fixed

rg_clec <- rectangleGate("FL17-A"=c(0,1e3), "SSC-A"=c(1e3, 1e6), filterId="CLEC4M")
gs_pop_add(fcsfiles.FS.GS, rg_clec, parent = "CD146+")

rg_clecm <- rectangleGate("FL17-A"=c(4e2,5e4), "SSC-A"=c(1e3, 1e6), filterId="CLEC4")
gs_pop_add(fcsfiles.FS.GS, rg_clecm, parent = "CD45-")

rg_HSC <- rectangleGate("FL3-A"=c(-1e3,1e3), "SSC-A"=c(1e3, 1e6), filterId="CD146-")
gs_pop_add(fcsfiles.FS.GS, rg_HSC, parent = "CD45-")

recompute(fcsfiles.FS.GS)
fcsfiles.FS.GS[[1]] |> autoplot(bins = 128)

# hello

ggcyto(fcsfiles.FS.GS[[1]], aes('FL13-A','SSC-A'), subset = "kupffer") +
  geom_hex(bins = 64) +
  xlim(-1e4,100)

ggcyto(fcsfiles.FS.GS[[5]], aes(x = `FL19-A`, y = `SSC-A`), subset = "notB") -> p
p + geom_hex(bins = 64) + scale_x_flowjo_biexp()

plot(fcsfiles.FS.GS)

# Part 1

gate_name <- "scatter"
ggcyto(fcsfiles.FS.GS[1],aes(`FSC-A`,`SSC-A`), subset = "root") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "singlet"
ggcyto(fcsfiles.FS.GS[1],aes(`FSC-Width`,`FSC-H`), subset = "singlet") +
  geom_hex(bins = 64) +
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7) +
  scale_y_log10() +
  geom_gate(gate_name)
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CD45+"
ggcyto(fcsfiles.FS.GS[1],aes(`FL8-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "Tcells"
ggcyto(fcsfiles.FS.GS[1],aes(`FL5-A`,`SSC-A`), subset = "CD45+") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "Bcells"
ggcyto(fcsfiles.FS.GS[1],aes(`FL19-A`,`SSC-A`), subset = "CD45+") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "NKcells"
ggcyto(fcsfiles.FS.GS[1],aes(`FL13-A`,`SSC-A`), subset = "CD45+") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "notT"
ggcyto(fcsfiles.FS.GS[1],aes(`FL5-A`,`SSC-A`), subset = "CD45+") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "notB"
ggcyto(fcsfiles.FS.GS[1],aes(`FL19-A`,`SSC-A`), subset = "notT") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "kupffer"
ggcyto(fcsfiles.FS.GS[1],aes(`FL13-A`,`SSC-A`), subset = "notB") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "vsig4+"
ggcyto(fcsfiles.FS.GS[1],aes(`FL1-A`,`SSC-A`), subset = "kupffer") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "vsig4-"
ggcyto(fcsfiles.FS.GS[1],aes(`FL1-A`,`SSC-A`), subset = "kupffer") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CD45-"
ggcyto(fcsfiles.FS.GS[1],aes(`FL8-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CD146-"
ggcyto(fcsfiles.FS.GS[1],aes(`FL3-A`,`SSC-A`), subset = "CD45-") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CD146+"
ggcyto(fcsfiles.FS.GS[1],aes(`FL3-A`,`SSC-A`), subset = "CD45-") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CLEC4"
ggcyto(fcsfiles.FS.GS[1],aes(`FL17-A`,`SSC-A`), subset = "CD146+") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CLEC4M"
ggcyto(fcsfiles.FS.GS[1],aes(`FL17-A`,`SSC-A`), subset = "CD146+") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

# Part 2

fcsfiles.FS <- read.flowSet(fcsfiles)
flowPlot(fcsfiles.FS[[1]], plotParameters = c("FL8-A","SSC-A"))
fcsfiles.FS.GS <- GatingSet(fcsfiles.FS)
gs_get_pop_paths(fcsfiles.FS.GS)
rg_scatter <- rectangleGate("FSC-A"=c(200,4e6), "SSC-A"=c(250, 5e6), filterId="scatter")
gs_pop_add(fcsfiles.FS.GS, rg_scatter,parent = "root") #Good plot

rg_singlet <- rectangleGate("FSC-Width"=c(200,4e3), "FSC-H"=c(250, 1e6), filterId="singlet")
gs_pop_add(fcsfiles.FS.GS, rg_singlet, parent = "scatter") #Good plot

rg_WBC <- rectangleGate("FL8-A"=c(4e3,1e6), "SSC-A"=c(1e4, 5e6), filterId="CD45+_scat")
gs_pop_add(fcsfiles.FS.GS, rg_WBC, parent = "singlet") #Good plot

rg_T <- rectangleGate("FL5-A"=c(8e3,1e6), "SSC-A"=c(1e4, 1e6), filterId="Tcells_scat")
gs_pop_add(fcsfiles.FS.GS, rg_T, parent = "singlet")

rg_B <- rectangleGate("FL19-A"=c(1e4,5e5), "SSC-A"=c(1e5, 5e6), filterId="Bcells_scat")
gs_pop_add(fcsfiles.FS.GS, rg_B, parent = "singlet")

rg_NK <- rectangleGate("FL13-A"=c(1e4,1e6), "SSC-A"=c(1e4, 5e6), filterId="NKcells_scat")
gs_pop_add(fcsfiles.FS.GS, rg_NK, parent = "singlet")

rg_NT <- rectangleGate("FL5-A"=c(1e2,6e3), "SSC-A"=c(1e4, 5e5), filterId="notT_scat")
gs_pop_add(fcsfiles.FS.GS, rg_NT, parent = "singlet")

rg_NB <- rectangleGate("FL19-A"=c(0,2e3), "SSC-A"=c(1e4, 5e5), filterId="notB_scat")
gs_pop_add(fcsfiles.FS.GS, rg_NB, parent = "singlet")

rg_K <- rectangleGate("FL13-A"=c(-3e2,0), "SSC-A"=c(1e4, 4e5), filterId="kupffer_scat")
gs_pop_add(fcsfiles.FS.GS, rg_K, parent = "singlet")

rg_Vsig <- rectangleGate("FL1-A"=c(4e2,2e3), "SSC-A"=c(1e4, 3e5), filterId="vsig4+_scat")
gs_pop_add(fcsfiles.FS.GS, rg_Vsig, parent = "singlet")

rg_VsigNK <- rectangleGate("FL1-A"=c(-5e2,3e2), "SSC-A"=c(1e4, 3e5), filterId="vsig4-_scat")
gs_pop_add(fcsfiles.FS.GS, rg_VsigNK, parent = "singlet")

rg_NBC <- rectangleGate("FL8-A"=c(-1000,3e3), "SSC-A"=c(1e3, 1e6), filterId="CD45-_scat")
gs_pop_add(fcsfiles.FS.GS, rg_NBC, parent = "singlet")

rg_OC <- rectangleGate("FL3-A"=c(1e3,5e4), "SSC-A"=c(1e3, 1e6), filterId="CD146+_scat")
gs_pop_add(fcsfiles.FS.GS, rg_OC, parent = "singlet")

rg_lsec <- rectangleGate("FL17-A"=c(0,5e4), "SSC-A"=c(1e2, 5e4), filterId="lsec_scat")
gs_pop_add(fcsfiles.FS.GS, rg_lsec, parent = "singlet") #Noah will address later scales need to be fixed

rg_clec <- rectangleGate("FL17-A"=c(0,1e3), "SSC-A"=c(1e3, 1e6), filterId="CLEC4M_scat")
gs_pop_add(fcsfiles.FS.GS, rg_clec, parent = "singlet")

rg_clecm <- rectangleGate("FL17-A"=c(4e2,5e4), "SSC-A"=c(1e3, 1e6), filterId="CLEC4_scat")
gs_pop_add(fcsfiles.FS.GS, rg_clecm, parent = "singlet")

rg_HSC <- rectangleGate("FL3-A"=c(-1e3,1e3), "SSC-A"=c(1e3, 1e6), filterId="CD146-_scat")
gs_pop_add(fcsfiles.FS.GS, rg_HSC, parent = "singlet")

recompute(fcsfiles.FS.GS)
gate_name <- "Tcells_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL5-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "Bcells_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL19-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "NKcells_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL13-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "notT_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL5-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "notB_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL19-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "kupffer_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL13-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "vsig4+_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL1-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "vsig4-_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL1-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CD45-_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL8-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CD146-_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL3-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CD146+_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL3-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CLEC4_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL17-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CLEC4M_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL17-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "CD146-_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL3-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

gate_name <- "Tcells_scat"
ggcyto(fcsfiles.FS.GS[1],aes(`FL5-A`,`SSC-A`), subset = "singlet") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p, width = 4, height = 4)

pdf("gating_tree.pdf")
plot(fcsfiles.FS.GS)
dev.off()

# 1. Complete the above set of plots according to the gating hierarchy (see
#    gating_tree.pdf)  
#
# 2. Plot everything again but this time change the subset to
#    "scatter" (i.e.: subset = "scatter")
# 3. Maybe I get credit on GitHub now
# 4. Reference the _scat plots to correct the gating parameters (just take some notes while looking at the pdfs, make a guess, and then go and plug them into the script)
# 5. Make smaller (complete/double-check)
