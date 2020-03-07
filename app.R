library(shiny)
library(shinyWidgets)
require(tidyverse)
require(annotables) #devtools
require(grImport)
require(BuenColors)  #devtools
require(cowplot)
require(edgeR) 
require(lobstr)
library("ggplot2")
library("magrittr")
library(shinythemes)
library(EnsDb.Hsapiens.v86)
library(ensembldb)
library(SNPlocs.Hsapiens.dbSNP149.GRCh38)
library(doParallel)
library(myvariant)
library(annotables)
library(rlist)
library(shinyjs)
library(BuenColors)
library(GenomicAlignments)
#options(repos = BiocInstaller::biocinstallRepos())
#options(repos = BiocManager::repositories())

genes <-annotables::grch38 %>%
  dplyr::select(symbol)

server <- function(input, output, session) {
  
  
  fix_names <- function(x){
    x %>%
      gsub("-", "", .) %>%
      str_squish(.) %>%
      gsub(" ", "_", .)
  }
  get_pos <- function(prot_conseq){
    str_extract_all(prot_conseq, "\\d+") %>%
      lapply(., function(x)x[[1]]) %>%
      unlist() %>%
      as.numeric()
  }
  wrapper <- function(df){
    voi <- c("stop_gained","frameshift_variant", "missense_variant",
             "splice_donor_variant", "splice_acceptor_variant", "start_lost")
    df %>%
      rename_all(~fix_names(names(df))) %>%
      dplyr::filter(Annotation %in% voi) %>%
      mutate(pos = get_pos(Protein.Consequence)) %>%
      dplyr::filter(!is.na(pos)) %>%
      mutate(variant = ifelse(Annotation == "missense_variant", "Missense", "LoF")) %>%
      dplyr::select(pos, Allele.Count, variant)
  }
  observeEvent(input$refresh, {
    shinyjs::reset("form")
    invalidateLater(1000000,session)
    v2_rawData<-c()
    v3_rawData<-c()
    get_domains() 
    get_v2_data()
  })
  v2_rawData <- eventReactive(input$v2, {
    read.csv(input$v2$datapath)
  })
  v3_rawData <- eventReactive(input$v3, {
    read.csv(input$v3$datapath)
  })
  update_choices <- function() {
    gene= input$genes[1]
    enst <- read_tsv('https://raw.githubusercontent.com/kopal-garg/resources/master/by_gene_enst.tsv')
    enst_row <- enst %>% dplyr::filter(enst$gene==input$genes[1])
    edb <- EnsDb.Hsapiens.v86
    txs <- transcripts(edb, filter = GeneNameFilter(input$genes[1]),
                       columns = c("protein_id", "tx_biotype",'uniprot_id'))
    up_txs<-data.frame(uniprot=txs$uniprot_id,tx_id=txs$tx_id,start=txs@ranges@start,end=txs@ranges@start+txs@ranges@width)
    gns <- up_txs %>% dplyr::filter(up_txs$tx_id == enst_row$transcript)
    gns <- subset(gns, select=-c(start,end))
    rel_json <- drawProteins::get_features(gns$uniprot[1])
    drawProteins::feature_to_dataframe(rel_json) -> rel_data
    Category <- lapply(rel_json[[1]][["features"]], function(x) x$category) %>% unlist()
    rel_data <- cbind(Category, rel_data)
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='VARIANT')
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='CONFLICT')
    rel_data <- rel_data %>% dplyr::filter(rel_data$description!='NONE')
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='MOD_RES')
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='COMPBIAS')
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='REGION')
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='CROSSLNK')
    rel_data_m <- rel_data %>% dplyr::filter(rel_data$Category == 'MOLECULE_PROCESSING')
    rel_data_d <- rel_data %>% dplyr::filter(rel_data$Category == 'DOMAINS_AND_SITES')
    if (nrow(rel_data_m)==0){rel_data <- rel_data_d}
    if (nrow(rel_data_d)==0){rel_data <- rel_data_m}
    if (nrow(rel_data_m)>0 && nrow(rel_data_d)>0){rel_data <- rbind(rel_data_m,rel_data_d)}
    #rel_data$type=str_replace(rel_data$type,'CHAIN','Scaffold')
    
    if(is.null(input$v3) | is.null(input$v2)) {
      transcript=enst_row$transcript
      q = paste0('dbnsfp.ensembl.transcriptid:', transcript)
      a_0_1000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=0)
      a_1001_2000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=1001)
      a_2001_3000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=2001)
      a_3001_4000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=3001)
      a_4001_5000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=4001)
      a_5001_6000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=5001)
      a_6001_7000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=6001)
      a_7001_8000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=7001)
      a_8001_9000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=8001)
      
      a_0_1000 <- subset(a_0_1000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_1001_2000 <- subset(a_1001_2000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_2001_3000 <- subset(a_2001_3000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_3001_4000 <- subset(a_3001_4000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_4001_5000 <- subset(a_4001_5000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_5001_6000 <- subset(a_5001_6000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_6001_7000 <- subset(a_6001_7000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_7001_8000 <- subset(a_7001_8000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_8001_9000 <- subset(a_8001_9000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      
      a <- rbind(a_0_1000,a_1001_2000, a_2001_3000, a_3001_4000, a_4001_5000,a_5001_6000,a_6001_7000,a_7001_8000,a_8001_9000)
      a <- a %>% dplyr::filter(a$gnomad_genome.ac.ac>0 | a$gnomad_exome.ac.ac>0)
      a$gnomad_exome.ac.ac[is.na(a$gnomad_exome.ac.ac)] <- 0
      a$gnomad_genome.ac.ac[is.na(a$gnomad_genome.ac.ac)] <-0
      a$Allele_Count <- a$gnomad_genome.ac.ac + a$gnomad_exome.ac.ac
      for (i in 1:nrow(a)){ 
        if(!is.null(a$snpeff.ann.protein.position[[i]][1])){
          a$pos[i] <-a$snpeff.ann.protein.position[[i]][1]}
        if(is.null(a$snpeff.ann.protein.position[[i]][1])){
          a$pos[i] <- 0}
        a$variant[i] = a$snpeff.ann.effect[[i]][1]
      }
      
      a <- a %>% dplyr::filter(pos!=0) 
      a <- subset(a, select= -c(`_score`,`snpeff._license`,`vcf.alt`,`vcf.position`,`vcf.ref`,`gnomad_exome._license`, `gnomad_genome._license`, `snpeff.ann.protein.length`,`snpeff.ann.protein.position`,`gnomad_exome.ac.ac`,`gnomad_genome.ac.ac`,`snpeff.ann.effect`))
      data_in <- data.frame(
        pos = a$pos,
        AC =a$Allele_Count,
        variant=a$variant
      )
      data_in$pos=as.numeric(levels(data_in$pos))[data_in$pos]
      data=data_in # in-built v2
    }
    else{
      v2=v2_rawData()
      v3=v3_rawData()
      v2 <- v2 %>% dplyr::filter(as.character(v2$Protein.Consequence)!="")
      v3 <- v3 %>% dplyr::filter(as.character(v3$Protein.Consequence)!="")
      v2 <- wrapper(v2)
      v3 <- wrapper(v3)
      
      data_uploaded <- full_join(v2, v3, by=c("pos","variant")) %>%
        mutate_all(., ~replace_na(., 0)) %>%
        mutate(AC = Allele.Count.x + Allele.Count.y) %>%
        group_by(pos)
      
      data=data_uploaded
      
    }
    updateSelectInput(session, 'Variant',choices=c(1:max(data$pos)),selected=1)
    updateSelectInput(session, "DomainChoices",choices = rel_data$type,selected=rel_data$type[1])
    val <- c(max(as.numeric(rel_data$begin)), max(as.numeric(rel_data$end)))
    updateSliderInput(session, 'CustomStart',min =1,value = 1,max  =max(val))
    updateSliderInput(session, 'CustomEnd',min =1 ,value = 1,max= max(val))

  }
  observe({
    update_choices()
    
  })
  get_transcript <- function() {
    gene= input$genes[1]
    enst <- read_tsv('https://raw.githubusercontent.com/kopal-garg/resources/master/by_gene_enst.tsv')
    enst_row <- enst %>% dplyr::filter(enst$gene==input$genes[1])
    enst_row$transcript
    
  }
  get_uniprot<- function() {
    canonical_transcript= get_transcript()
    edb <- EnsDb.Hsapiens.v86
    txs <- transcripts(edb, filter = GeneNameFilter(input$genes[1]),
                       columns = c("protein_id", "tx_biotype",'uniprot_id'))
    up_txs<-data.frame(uniprot=txs$uniprot_id,tx_id=txs$tx_id,start=txs@ranges@start,end=txs@ranges@start+txs@ranges@width)
    gns <- up_txs %>% dplyr::filter(up_txs$tx_id == canonical_transcript)
    gns <- subset(gns, select=-c(start,end))
    gns$uniprot[1]
    
  }
  output$v2_link <- renderText({
    canonical_transcript=get_transcript()
    paste0('https://gnomad.broadinstitute.org/transcript/', canonical_transcript,'?dataset=gnomad_r2_1')
    
  })
  output$v3_link <- renderText({
    canonical_transcript=get_transcript()
    paste0('https://gnomad.broadinstitute.org/transcript/',canonical_transcript, '?dataset=gnomad_r3')
  })
  get_domains<- function() {
    rel_json <- drawProteins::get_features(get_uniprot())
    drawProteins::feature_to_dataframe(rel_json) -> rel_data
    Category<- lapply(rel_json[[1]][["features"]], function(x) x$category) %>% unlist()
    rel_data <- cbind(rel_data, Category)
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='VARIANT')
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='CONFLICT')
    rel_data <- rel_data %>% dplyr::filter(rel_data$description!='NONE')
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='MOD_RES')
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='COMPBIAS')
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='REGION')
    rel_data <- rel_data %>% dplyr::filter(rel_data$type!='CROSSLNK')
    rel_data_m <- rel_data %>% dplyr::filter(rel_data$Category == 'MOLECULE_PROCESSING')
    rel_data_d <- rel_data %>% dplyr::filter(rel_data$Category == 'DOMAINS_AND_SITES')
    if (nrow(rel_data_m)==0){rel_data <- rel_data_d}
    if (nrow(rel_data_d)==0){rel_data <- rel_data_m}
    if (nrow(rel_data_m)>0 && nrow(rel_data_d)>0){rel_data <- rbind(rel_data_m,rel_data_d)}
    #rel_data$type=str_replace(rel_data$type,'CHAIN','Scaffold')
    # val <- c(max(as.numeric(rel_data$begin)), max(as.numeric(rel_data$end)))
    # updateSliderInput(session, 'CustomStart',min =1,value = 1,max  =max(val))
    # updateSliderInput(session, 'CustomEnd',min =1 ,value = 1,max= max(val))
    
    rel_data
  }
  output$DomainTable_out <- DT::renderDataTable({
    domain <- get_domains()
    domain <- subset(domain,select=-c(length,accession,entryName,taxid,order,Category))
    DT::datatable(domain, options = list(editable = 'row', searching = FALSE, paging= FALSE), selection = list(target = 'row'))
  }, selection = list(target = 'rows'))
  
  
  get_v2_data<- function(){
    if(is.null(input$v3) | is.null(input$v2)) {
      transcript=get_transcript()
      q = paste0('dbnsfp.ensembl.transcriptid:', transcript)
      a_0_1000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=0)
      a_1001_2000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=1001)
      a_2001_3000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=2001)
      a_3001_4000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=3001)
      a_4001_5000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=4001)
      a_5001_6000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=5001)
      a_6001_7000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=6001)
      a_7001_8000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=7001)
      a_8001_9000<-queryVariant(q=q,fields='vcf,gnomad_exome.ac.ac,gnomad_genome.ac.ac,snpeff.ann.effect,snpeff.ann.protein.position,snpeff.ann.protein.length,mutdb.uniprot_id,vcf',dotfield=T,size=1000,from=8001)
      a_0_1000 <- subset(a_0_1000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_1001_2000 <- subset(a_1001_2000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_2001_3000 <- subset(a_2001_3000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_3001_4000 <- subset(a_3001_4000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_4001_5000 <- subset(a_4001_5000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_5001_6000 <- subset(a_5001_6000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_6001_7000 <- subset(a_6001_7000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_7001_8000 <- subset(a_7001_8000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      a_8001_9000 <- subset(a_8001_9000$hits, select=c("_id","_score","snpeff._license",'snpeff.ann.effect','snpeff.ann.protein.length','snpeff.ann.protein.position','vcf.alt','vcf.ref','vcf.position','gnomad_exome._license','gnomad_exome.ac.ac','gnomad_genome.ac.ac','gnomad_exome._license','gnomad_genome._license'))
      
      a <- rbind(a_0_1000,a_1001_2000, a_2001_3000, a_3001_4000, a_4001_5000,a_5001_6000,a_6001_7000,a_7001_8000,a_8001_9000)
      
      a <- a %>% dplyr::filter(a$gnomad_genome.ac.ac>0 | a$gnomad_exome.ac.ac>0)
      a$gnomad_exome.ac.ac[is.na(a$gnomad_exome.ac.ac)] <- 0
      a$gnomad_genome.ac.ac[is.na(a$gnomad_genome.ac.ac)] <-0
      a$Allele_Count <- a$gnomad_genome.ac.ac + a$gnomad_exome.ac.ac
      for (i in 1:nrow(a)){ 
        if(!is.null(a$snpeff.ann.protein.position[[i]][1])){
          a$pos[i] <-a$snpeff.ann.protein.position[[i]][1]}
        if(is.null(a$snpeff.ann.protein.position[[i]][1])){
          a$pos[i] <- 0}
        a$variant[i] = a$snpeff.ann.effect[[i]][1]
      }
      
      a <- a %>% dplyr::filter(pos!=0)
      a <- subset(a, select= -c(`_score`,`snpeff._license`,`vcf.alt`,`vcf.position`,`vcf.ref`,`gnomad_exome._license`, `gnomad_genome._license`, `snpeff.ann.protein.length`,`snpeff.ann.protein.position`,`gnomad_exome.ac.ac`,`gnomad_genome.ac.ac`,`snpeff.ann.effect`))
      data_in <- data.frame(
        pos = a$pos,
        AC =a$Allele_Count,
        variant=a$variant
      )
      data_in$pos=as.numeric(levels(data_in$pos))[data_in$pos]
      data=data_in # in-built v2
    }
    else{
      v2=v2_rawData()
      v3=v3_rawData()
      v2 <- v2 %>% dplyr::filter(as.character(v2$Protein.Consequence)!="")
      v3 <- v3 %>% dplyr::filter(as.character(v3$Protein.Consequence)!="")
      v2 <- wrapper(v2)
      v3 <- wrapper(v3)
      
      data_uploaded <- full_join(v2, v3, by=c("pos","variant")) %>%
        mutate_all(., ~replace_na(., 0)) %>%
        mutate(AC = Allele.Count.x + Allele.Count.y) %>%
        group_by(pos)
      
      data=data_uploaded
    }
    data
  }
  output$Transcript <- renderText({
    canonical_transcript=get_transcript()
    paste0('Canonical Transcript: ', canonical_transcript)
  })
  output$SelectedGene <- renderText({
    paste0('Selected Gene: ', input$genes)
    
  })
  output$selected_rows = renderText({
    input$DomainTable_out_rows_selected
  })
  output$Uniprot <- renderText({
    uniprotID=get_uniprot()
    paste0('Uniprot ID: ', uniprotID)
  })
  output$DomainTable <- renderTable({
    domain <- get_domains()
    domain_new <- domain %>% dplyr::filter(domain$type %in% input$DomainChoices)
    domain_df <- data.frame(
      x1=domain_new$begin,
      x2=domain_new$end,
      color=domain_new$description,
      y1=c(0.999),
      y2=c(1.001)
    )
    domain_df
  })
  output$plot <- renderPlot({
    plot_fig()
  })
  plot_fig <- function() {
    domain <- get_domains()
    chain <- domain %>% dplyr::filter(domain$type == "CHAIN")
    domain_new <- domain %>% dplyr::filter(domain$type %in% input$DomainChoices)
    selected_rows <- input$DomainTable_out_rows_selected
    #selected <- cat(paste(shQuote(selected_rows, type="cmd"), collapse=", "))
    
    domain <- domain[selected_rows , ]
    
    domain_new <- rbind(chain, domain)
    domain_new <- unique(domain_new)
    #domain_new <- subset(domain_new, select = c('type', 'description', 'begin', 'end'))
    print(domain_new)
    gc()
    print(lobstr::mem_used())
    
    if (input$CustomDomain != ""){
      custom_domain <- c()
      custom_domain$type= input$CustomDomain
      custom_domain$description = 'Custom'
      custom_domain$begin = as.numeric(input$CustomStart)
      custom_domain$end = as.numeric(input$CustomEnd)
      custom_domain$length = 22
      custom_domain$accession = 'accession'
      custom_domain$entryName = 'entryName'
      custom_domain$taxid = 9606
      custom_domain$order  = 1
      custom_domain$Category = 'Category'
      
      custom_domain <- as.data.frame(custom_domain)
      domain_new <- rbind(domain_new, custom_domain)}
    
    
    domain_df <- data.frame(
      x1=domain_new$begin,
      x2=domain_new$end,
      color=domain_new$type,
      y1=c(0.999),
      y2=c(1.001)
    )
    variants=as.numeric(input$Variant)
    variants_df <- data.frame(pos = variants, y=1.0012)
    data = c()
    data=get_v2_data()
    palette <- as.character(jdb_palette("lawhoops"))
    palette <- list.append('#A9A9A9',palette, palette)
    
    palette_new <- palette[seq.int(1L, length(palette), 4)]
    
    # a <- c()
    # types=unique(domain_df$type)
    # types = types[types!= 'Scaffold']
    # for (i in 1:length(types)){
    #   a$domain_type[i] <- types[i]
    #   a$color[i] <- palette_new[i]
    # }
    p=ggplot(data, aes(x = pos, y = 1, size = AC)) +
      geom_rect(domain_df, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,
                                       fill=color), color="black",
                alpha=0.5, inherit.aes = FALSE) +
      geom_point(alpha = 0.5) +
      geom_point(data=variants_df, mapping=aes(x = pos, y=y), shape=25, fill="red",
                 inherit.aes = FALSE, size=3) +
      pretty_plot(fontsize = 9) + scale_size(trans = "log10") +
      scale_fill_manual(values =c(palette_new)) +
      L_border() +
      labs(x = "Canonical Protein Position", y = "", size = "gnomAD\nallele count",
           fill = "domain") + 
      theme(legend.position = "bottom")+
      theme(axis.title.x=element_blank(),
            axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    p
  }
  output$download <- downloadHandler(
    filename = paste0(input$genes,"_output.pdf"), 
    content = function(file) {
      ggsave(filename= file, width = 7, height=1.6)
      print(plot_fig())
      dev.off()
    })
} 

ui <- fluidPage(
  useShinyjs(),
  title = "gnomAD Gene Model",
  hr(),
  div(id = "form",
      downloadButton(outputId = 'download', label = 'Save'),
      actionButton("refresh", "Refresh"),
      plotOutput('plot'),
      hr(),
      fluidRow(
        column(4,
               fileInput("v2", "gnomAD V2", multiple=FALSE,accept = c("text/csv",
                                                                      "text/comma-separated-values,text/plain",
                                                                      ".csv")),
               fileInput("v3", "gnomad V3", multiple=FALSE,accept=  c("text/csv",
                                                                      "text/comma-separated-values,text/plain",
                                                                      ".csv")),
               textOutput('v2_link'),
               textOutput('v3_link')),
        column(4, 
               selectInput(inputId = 'genes',
                           label = 'Select Gene:',
                           choices = genes, selected = 'BCL11A'),
               #selectInput(inputId = 'DomainChoices', label='Select Domain: ',multiple = TRUE,selected=c('CHAIN'),choices = c('CHAIN')),
               selectInput(inputId = 'Variant', label='Select Positions: ', multiple=TRUE,choices = c(1,2,3,4,5,6,7,8,9,10),selected = 1),
        ),
        
        column(4,textInput('CustomDomain', label = 'Custom Domain Type: '),
               sliderInput('CustomStart', label  = 'Start: ', value = 0,min=0,max =1000),
               sliderInput('CustomEnd', label = 'End: ',value= 0,min=0,max =1000)),
        DT::dataTableOutput('DomainTable_out')
        #textOutput(outputId = 'selected_rows'),
      )))




shinyApp(ui = ui, server = server)
