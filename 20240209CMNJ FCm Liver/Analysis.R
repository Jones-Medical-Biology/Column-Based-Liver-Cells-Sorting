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
gs_pop_add(fcsfiles.FS.GS, rg_scatter,parent = "root")

rg_WBC <- rectangleGate("FL8-A"=c(4e3,1e6), "SSC-A"=c(1e4, 5e6), filterId="CD45+")
gs_pop_add(fcsfiles.FS.GS, rg_WBC, parent = "scatter")

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
gs_pop_add(fcsfiles.FS.GS, rg_NBC, parent = "scatter")

rg_OC <- rectangleGate("FL3-A"=c(1e3,5e4), "SSC-A"=c(1e3, 1e6), filterId="CD146+")
gs_pop_add(fcsfiles.FS.GS, rg_OC, parent = "CD45-")

rg_lsec <- rectangleGate("FL17-A"=c(0,5e4), "SSC-A"=c(1e2, 5e4), filterId="lsec")
gs_pop_add(fcsfiles.FS.GS, rg_lsec, parent = "CD146+") #Noah will address later scales need to be fixed

rg_clec <- rectangleGate("FL17-A"=c(0,5e4), "SSC-A"=c(1e2, 5e4), filterId="CLEC4M")
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

gate_name <- "CD146+"
ggcyto(fcsfiles.FS.GS[1],aes(`FL17-A`,`SSC-A`), subset = "CD45-") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name)
ggsave(paste(gate_name, sep = "", ".pdf"), p)

gate_name <- "CLEC4"
ggcyto(fcsfiles.FS.GS[1],aes(`FL17-A`,`SSC-A`), subset = "CD45-") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p)

gate_name <- "CLEC4M"
ggcyto(fcsfiles.FS.GS[1],aes(`FL17-A`,`SSC-A`), subset = "CD146+") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p)

gate_name <- "CD146-"
ggcyto(fcsfiles.FS.GS[1],aes(`FL3-A`,`SSC-A`), subset = "CD45-") +
  geom_hex(bins = 64)+
  scale_x_flowjo_biexp(widthBasis = -1, maxValue = 1e7)+
  scale_y_log10()+
  geom_gate(gate_name) -> p
ggsave(paste(gate_name, sep = "", ".pdf"), p)

pdf("gating_tree.pdf")
plot(fcsfiles.FS.GS)
dev.off()

# 1. Complete the above set of plots according to the gating hierarchy (see
#    gating_tree.pdf) 
#
# 2. Plot everything again but this time change the subset to
#    "scatter" (i.e.: subset = "scatter")
