# ChIP-seq

Data

There are four datasets \(paired-end reads\):

| Dataset | Description |
| :--- | :--- |
| BCL11A.1 | ChIP experiment, Replicate 1 |
| Input.1 | Input DNA, Replicate 1 |
| BCL11A.2 | ChIP experiment, Replicate 2 |
| Input.2 | Input DNA, Replicate 2 |

Workflow:

\(1\) Align reads to hg19 using BWA \(output: sorted BAM alignment files + index\)

```text
bwa aln -t 4 resources/hg19bwaidx $fastq > out/$out.fq.sai
bwa sampe resources/hg19bwaidx out/$out1.fq.sai out/$out2.fq.sai  $fastq1 $fastq2  > out/$out.sam
samtools view -b out/$out.sam > out/$out.bam
samtools sort -o -O out/$out.bam -T sorted out/st.$out.bam
samtools index -b out/st.$out.bam out/st.$out.bam.bai
```

\(2\) Filter out unmapped and multiple-mapped reads \(having multiple-mapped reads increases sensitivity of peak detection, and number of usable reads\). 

```text
samtools view -h -F 4 -b out/st.$out.bam > out/st.mapped.$out.bam
samtools index -b out/st.mapped.$out.bam out/st.mapped.$out.bam.bai
```

\(3\) Peak calling \(MACS2\) returns BED/GTF file

\(4\) Convert to bigWig for genome browsers

 

![Overlap Input and BCL11A tracks to see peaks](.gitbook/assets/image%20%2815%29.png)



