# ANKRD26

| chr | start | end | gene |
| :--- | :--- | :--- | :--- |
| chr10 | 26991914 | 27070472 | ANKRD26ex10+ |
| chr10 | 27070473 | 27100498 | ANKRD26ex1\_9 |
| chr10 | 28532493 | 28623112 | WAC |

make region-specific GTF for featureCounts

```text
bedtools intersect -b WAC.txt -a gencode.v33.annotation.gtf > WAC.gtf
bedtools intersect -b ANKRD26ex10 -a gencode.v33.annotation.gtf > ANKRD26ex10.gtf
bedtools intersect -b ANKRD26ex1_9.txt -a gencode.v33.annotation.gtf > ANKRD26ex1_9.gt
```

use featureCounts to get the number of reads mapping exons in each gene

input: case/control aligned, sorted and indexed BAMs and region-specific GTF

\(BAMs: after rsem index and alignment on raw paired-end reads\)

count reads in case-specific BAMs and control-specific BAMs for each region, and bin per gene. example command:

```text
cd /broad/hptmp/kgarg/ANKRD26/resources/featureCount/subread-2.0.0-Linux-x86_64/bin
./featureCounts -T 4 -a ../../../ANKRD26ex10.gtf -o ../../../ANKRD26ex10cases_out.txt -O -g gene_name /broad/hptmp/kgarg/ANKRD26/out/SH_S2_R1Aligned.out.bam /broad/hptmp/kgarg/ANKRD26/out/CC_S3_R1Aligned.out.bam /broad/hptmp/kgarg/ANKRD26/out/BAC_S1_R1Aligned.out.bam
# parse output
cat ANKRD26ex10cases_out.txt| awk '{print $1"\t"$7"\t"$8"\t"$9}'> ANKRD26ex10cases_sub.txt
```

read all files \(2 per region\) into R and merge

```text
setwd('/broad/hptmp/kgarg/ANKRD26/resources')
library(matrixStats)
library(dplyr)
library(ggplot2)

ANKRD26ex10cases=fread('ANKRD26ex10cases_sub.txt')
ANKRD26ex10controls=fread('ANKRD26ex10controls_sub.txt')
ANKRD26ex1_9cases=fread('ANKRD26ex1_9cases_sub.txt')
ANKRD26ex1_9controls=fread('ANKRD26ex1_9controls_sub.txt')
WAC_controls=fread('WACcontrols_sub.txt')
WAC_cases=fread('WACcases_sub.txt')

#merge, re-name some things and bind

ANK19 = merge(ANKRD26ex1_9cases, ANKRD26ex1_9controls,by='Geneid')
ANK10 = merge(ANKRD26ex10cases, ANKRD26ex10controls,by='Geneid')
WAC = merge(WAC_cases, WAC_controls,by='Geneid')
ANK10$Geneid[1] = 'ANKRD26_10'
ANK19$Geneid[1] = 'ANKRD26_19'
bind <- rbind(ANK10, ANK19, WAC)
bind_num <- subset(bind, select = -Geneid)
```

counts per million function from before.

do log2\(cpm\)

```text
cpm <- function(x, log = F ){
  out <- (x/sum(x))*1000000
  if(log == T){
    out <- log2(out)
  }
  out
}
bind_cpm <- cbind(subset(bind,select=Geneid),cpm(bind_num, log=T)) 
```

combine counts for cases and controls

```text
bind_cpm <- bind_cpm %>% 
  mutate(case = bind_cpm$SH + bind_cpm$CC + bind_cpm$BAC) %>%
  mutate(control = bind_cpm$VGS + bind_cpm$NH + bind_cpm$DH) %>%
  subset(., select=c(Geneid,case,control)) %>%
  filter(Geneid!='RNU6-490P') 
```

find mean and sd for each row

```text
bind_cpm <- cbind(bind_cpm, 
                  rowMeans(bind_cpm %>% subset(., select = -Geneid)))

bind_cpm <- cbind(bind_cpm,                  
                  rowSds(x=as.matrix(bind_cpm %>% subset(., select = -c(Geneid, V2)), na.rm=TRUE)))
bind_num = bind_cpm
bind_cpm <- bind_cpm %>%
  `colnames<-`(c('Geneid', 'case', 'control', 'mean', 'sd')) %>%
  gather(., key=tag,value=count,-Geneid, -mean, -sd) 
```

facet by gene

```text
ggplot(bind_cpm %>% group_by(Geneid), aes(x = tag, y =count,fill=as.factor(tag))) + 
  geom_bar(stat ='identity') + 
  facet_wrap(~Geneid, nrow=1) +
  scale_fill_manual(values = c('red', 'blue')) +
  geom_errorbar(data=bind_cpm,aes(x=tag, ymin=count-sd, ymax=count+sd), width=0.2, colour="black", alpha=0.8, size=1)+
  ylab('log2(cpm)')
```

![](.gitbook/assets/image%20%2823%29.png)

