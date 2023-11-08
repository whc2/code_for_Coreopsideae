#!/usr/bin/env Rscript

#load package
library(ggplot2)
library(RColorBrewer)

#read data
args <- commandArgs(TRUE)
if (length(args) != 1) {
    cat("Rscript intact_LTR_time_density.R LTR_Insert_sum.tsv\n")
	quit()
}
dat <- read.table(args[1], sep="\t", header=TRUE)
#select the desired columns
#dat <- dat[,c(8, 10, 12)]
#colnames(dat) <- c("Identity", "SuperFamily", "Insertion_Time")
colnames(dat) <- c("SuperFamily", "Insertion_Time")

#draw density
dens.p <- ggplot(dat, aes(x=Insertion_Time, color=SuperFamily)) +
			geom_density(adjust=1) +
			xlim(0, 2000000) +
#			scale_color_brewer(palette="Set2") +
			theme(line=element_line(color="black", size=0.2), text=element_text(color="black", size=8),
		          axis.ticks=element_line(color="black", size=0.2), axis.text=element_text(color="black", size=8))#,
#		          panel.background=element_rect(fill="white", color="black"), panel.grid=element_blank())
#ggsave(paste(args[1], "LTR_Insertion_Time.pdf", sep="."), dens.p, device="svg", width=14, height=8, units="cm")
ggsave(paste(args[1], "_Insertion_Time.pdf", sep=""), dens.p, device="pdf", width=14, height=8, units="cm")

