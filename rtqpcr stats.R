library(readxl)
rtqpcr_stats <- read_excel("~/Pinto/DO Transcriptomics Project Data/rtqpcr stats.xlsx")
summary(aov(copies ~do, data=rtqpcr_stats))

summary(aov(aob ~zone, data=rtqpcr_stats))
summary(aov(aob ~zone*do, data=rtqpcr_stats))
summary(aov(cmx ~zone*do, data=rtqpcr_stats))
summary(aov(nxr ~zone*do, data=rtqpcr_stats))
summary(aov(hzo ~zone*do, data=rtqpcr_stats))


ss_rt <-subset(rtqpcr_stats, zone == c("anox","aero"))
summary(aov(nxr ~zone*do, data=ss_rt))

ifas_rt <-subset(rtqpcr_stats, zone == "ifas")
summary(aov(cmx ~do, data=ifas_rt))
cmx1 <-aov(cmx ~do, data=ifas_rt)
