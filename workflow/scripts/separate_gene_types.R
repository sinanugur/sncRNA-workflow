#!/usr/bin/env Rscript

####################################################################################################
## Script: separate_gene_types.R                                                              
## Date:   24 Feb 2017                                                                            
## Author: Sinan Ugur Umu (SUU), sinanugur@gmail.com                                              
## Latest Update: 16 Nov 2017                                                                     
##                                                                                                
## Description: Pull out some genes from the gencode table.                                                                                                
##                                                                                                
####################################################################################################


keep = c("lncRNA","misc_RNA","protein_coding","snoRNA","snRNA","scaRNA")

#setwd("/WORKING/projects/miRJanus/analyses/main_results/count_tables_for_all_analyses")


#gencode_combined_table=read.csv("/WORKING/projects/miRJanus/analyses/main_results/genome_mapping_analyses/gene.gencode.combined.table.csv",header = T,sep="\t")

arguments=commandArgs(TRUE)



if(length(arguments) >= 1) {

	  #output=tools::file_path_sans_ext(arguments[1])
	gencode_combined_table=read.csv(arguments[1],header = T,sep="\t")
  output_directory=paste0(dirname(arguments[1]),"/")
} else {
	gencode_combined_table=read.csv("exons.gencode.combined.table.csv",header = T,sep="\t")
  output_directory=""
}


rownames(gencode_combined_table)=gencode_combined_table[,1]

gencode_gene_table=read.csv("../../databases/gencode.gene-name.csv",header = T,sep="\t")


select_out_type=function(rna_type){
  
  selected_df=gencode_combined_table[as.character(gencode_gene_table[gencode_gene_table$Type == rna_type,]$ID),]
  
  write.table(selected_df,paste0(output_directory,rna_type,".tsv"),row.names = F,sep = "\t")
  
  return()
}

sapply(keep, select_out_type)


  


