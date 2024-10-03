library(glue) ; g=glue ; len=length
library(gridExtra)
library(parallel)
library(magrittr)
library(ggplot2)
library(cowplot)
library(stringr)
library(rlist)
library(dplyr)
library(purrr)
library(qs)

fastqpath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/02_FASTQS/231108_VL00297_217_AAF52C3M5/outs/fastq_path/AAF52C3M5"
system("mkdir fastqs")
system(g("gsutil cp {fastqpath}/SI-NT-D9_*_R2_*.fastq.gz fastqs"))

puckpath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/pucks"
system("mkdir pucks")
system(g("gsutil cp {puckpath}/Puck_230525*.csv pucks"))

fastq = "fastqs/SI-NT-D9_S2_L001_R2_001.fastq.gz"

library(ShortRead) # sread, quality, id
R = fastq %>% FastqStreamer(n = 100000) %>% yield %>% sread %>% subseq(start=1, end=32) %>% as.data.frame
R %<>% tidyr::separate(x, into=c("sb1","up","sb2"), sep=c(8,26))
R %<>% filter(up=="TCTTCAGCGTTCCCGAGA") %>% transmute(sb=paste0(sb1,sb2))

pucks = map(list.files("pucks",full.names=T), function(puck){read.csv(puck)[[1]]}) %>% setNames(list.files("pucks"))

for (i in 1:len(pucks)) {
  R %<>% mutate(!!(names(pucks)[[i]]) := (sb%in%pucks[[i]]))
}

df = R %>% select(-1) %>% map_int(sum) %>% {data.frame(a=names(.),b=unname(.))} %>% mutate(b=round(b/nrow(R)*100,2) %>% paste0("%"))
View(df)
