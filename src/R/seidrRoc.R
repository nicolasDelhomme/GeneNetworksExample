# library
suppressPackageStartupMessages({
  library(here)
  library(pracma)
})

# Genie3
dat <- read.delim(
  here("analysis/gene_network/seidr/aggregated-genie3.roc"),
  header=FALSE,skip=1,nrows = 100,
  col.names = c("TP","FP","PR","ALGO"))

auc <- round(trapz(dat[,2],dat[,1]),digits=3)

plot(dat[,2],dat[,1],type="l",main=sprintf("%s (AUC = %s)","IRP",auc),
     xlab="False Positive Rate",ylab="True Positive Rate")

abline(0,1,lty=2)

# IRP
dat <- read.delim(
  here("analysis/gene_network/seidr/aggregated-irp.roc"),
  header=FALSE,skip=1,nrows = 100,
  col.names = c("TP","FP","PR","ALGO"))

auc <- round(trapz(dat[,2],dat[,1]),digits=3)

plot(dat[,2],dat[,1],type="l",main=sprintf("%s (AUC = %s)","IRP",auc),
     xlab="False Positive Rate",ylab="True Positive Rate")

abline(0,1,lty=2)

# Backbone
dat <- read.delim(
  here("analysis/gene_network/seidr/backbone-10-percent.roc"),
  header=FALSE,skip=1,nrows = 100,
  col.names = c("TP","FP","PR","ALGO"))

auc <- round(trapz(dat[,2],dat[,1]),digits=3)

plot(dat[,2],dat[,1],type="l",main=sprintf("%s (AUC = %s)","IRP",auc),
     xlab="False Positive Rate",ylab="True Positive Rate")

abline(0,1,lty=2)
