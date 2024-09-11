library(dplyr)
library(magrittr)
library(GenomicRanges)
library(knitr)
library(ggplot2)
library(tidyr)

## remove hash to set the working unless you already did in console.
#setwd("path/to/deltafiles")

readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

## this is a function to run one or multiple delta files.

delta_files <- list.files(pattern = "*.delta")

for (file in delta_files){
mumgp = readDelta(file)
mumgp = mumgp %>% filter(strand == "+")

## remove hash if you also want negative sense strand.
#mumgp = mumgp %>% filter(strand == "-")
mumgp %>% head %>% kable

## this will help create a title name using the "qid" made my nucmer.
## unique() is used to work around all the rows created with the same id.

unique_id <- unique(mumgp$qid)

## this will piece together a title

title_text <- paste("Alignment Profile:", unique_id, "vs. Reference Genome")


## this creates simple mummer plot using data from .delta file

ggplot(mumgp, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) + 
  geom_segment() +
  geom_point(alpha=.5) +
  geom_point(aes(x = re, y = qe), alpha=.5) + 
  labs(title=title_text, x="Reference Sequence", 
       y="Consensus Sequence") +
  facet_grid(qid~., scales='free', space='free', switch='y') +
  theme_bw() + 
  theme(plot.title = element_text(hjust=0.5), 
        #strip.text.y=element_text(angle=180, size=5), ##changes the size of y-axis text
        legend.position=c(.99,.01), 
        legend.justification=c(1,0),
        strip.background=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_colour_brewer(palette='Set1') +
  scale_x_continuous(limits=c(0, 1000), breaks=seq(0, 1000, by=100)) +
  scale_y_continuous(limits=c(0, 1000), breaks=seq(0, 1000, by=100)) 

## this will save it as a high resolution plot
output_file <- paste0(tools::file_path_sans_ext(file), "_mummer.png")
ggsave(filename = output_file,
       width = 10,
       height = 6,
       dpi = 300)
}
