#!/usr/bin/env Rscript

####################################################################################################
## Script: generic_table_creator.R                                                                ##
## Date:   13 Jan 2017                                                                            ##
## Author: Sinan Ugur Umu (SUU), sinanugur@gmail.com                                              ##
## Latest Update: 15 May 2017                                                                     ##
##                                                                                                ##
##                                                                                                ##
##                                                                                                ##
##                                                                                                ##
##                                                                                                ##
##                                                                                                ##
## Description: This code will generate a data table of all reads with a minimum                  ##
##            expression of 100 read counts for at least 1 sample.                                ##
##                                                                                                ##
##            It is very similar to "isomirna_table_creator.R". The input data was created by     ##
##            bowtie2 and rnaseq_tool.py                                                          ##
##                                                                                                ##
##                                                                                                ##
##                                                                                                ##
## Input:  *(argument[1]).txt                                                                     ##
## Output: argument[1].count.table.csv                                                            ##
##                                                                                                ##
##                                                                                                ##
##                                                                                                ##
##                                                                                                ##
##                                                                                                ##
##                                                                                                ##
##                                                                                                ##
####################################################################################################

arguments=commandArgs(TRUE)
require(dplyr)


#file_type=paste0("*.\\.",arguments[1],".txt$") change this to

file_type=paste0("*.\\.",arguments[1],"[.txt$|.tsv$]") #this for seqbuster 12/10/2017


if(length(arguments) > 1) {

  setwd(arguments[2])
  
}

read_file_list=list.files("./",file_type,full.names = T)  
print(read_file_list)






read_files=function(file_list)
{
  for (i in file_list){
    
    
    if(exists("file_df")){
      
      #file_df=merge(file_df,read.csv(i,sep = "\t", header = T),by=mergename,all=T)
      
      file_df_tmp=read.csv(i,sep = "\t", header = T)
      file_df_tmp[,1]=as.character(file_df_tmp[,1])
      file_df_tmp[,2]=as.integer(file_df_tmp[,2])
      
      file_df = file_df %>% full_join(file_df_tmp,by=mergename)
      
      
    }
    
    else {
      
      file_df=read.csv(i,sep = "\t", header = T,stringsAsFactors = F)
      
      mergename=names(file_df)[1]
      file_df[,1]=as.character(file_df[,1])
      file_df[,2]=as.integer(file_df[,2])
      
      
    }
    
    
    
  }
  
  return(file_df)
  
}


read_file_df=read_files(read_file_list)

read_file_df[is.na(read_file_df)]=0

write.table(read_file_df,paste0(arguments[1],".tsv"),row.names = F,sep = "\t")



