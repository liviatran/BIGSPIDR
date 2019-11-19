#atlasMaker 0.3
#by: Livia Tran
#11/17/19

###This function downloads nucleotide alignment sequences from the ANHIG/IMGTHLA Github Repository
#and finds exon boundaries within the nucleotide alignment to determine which amino acids in the 
#protein alignment belong to which exon


library(stringr)
library(DescTools)

atlasMaker<-function(loci){
  pep_start<-nuc<-nuc_df<-extract_ref<-nuc_extract<-start<-end<-pipe_split<-boundaries<-boundary_split<-atlas<-sapply(loci, function(x) NULL)
  
  for(i in 1:length(loci)){
    #download nucleotide alignment from ANHIG/IMGTHLA Github repository
    nuc[[loci[[i]]]]<-readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",paste(ifelse(loci[[i]] %in% c("DRB1", "DRB3", "DRB4", "DRB5"),"DRB",loci[[i]]),"_nuc.txt",sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)
    
    #reduces repeated whitespace between allele and nucleotide
    nuc[[loci[[i]]]]<-str_squish(nuc[[loci[[i]]]])
    
    #remove impertinent header/footer information
    nuc[[loci[[i]]]] <- head(nuc[[loci[[i]]]],-2)
    nuc[[loci[[i]]]] <- tail(nuc[[loci[[i]]]],-6)
    
    #removes all whitespace, except for the whitespace between the allele and nucleotide sequence
    nuc[[loci[[i]]]]<-paste(substr(nuc[[loci[i]]],1,regexpr(" ",text = nuc[[loci[i]]],fixed = TRUE)), gsub(" ","",substr(nuc[[loci[i]]],regexpr(" ",text = nuc[[loci[i]]],fixed = TRUE),nchar(nuc[[loci[i]]]))),sep = "")
    
    #splits at white spaces to yield allele and nucleotide sequences
    nuc[[loci[i]]]  <- strsplit(nuc[[loci[i]]]," ", fixed=T)
    
    #binds the previously split strings by row 
    nuc[[loci[i]]] <- do.call(rbind,nuc[[loci[i]]])
    
    #extracts beginning alignment enumeration
    pep_start[[loci[[i]]]]<-as.numeric(gsub("codon", "", nuc[[loci[[i]]]][[2,2]]))
    
    colnames(nuc[[loci[[i]]]])<-c(paste(loci[[i]], "alleles", sep="_"), "pepseq")
    
    #finds start and end of each "cDNA" block 
    start[[loci[i]]]<-as.numeric(grep("cDNA", nuc[[loci[i]]]))
    end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,nrow(nuc[[loci[i]]])))
    
    #due to ANHIG formatting, cases where an allele contains newly reference nucleotide sequences will not
    #contain the same number of rows as previous reference peptide blocks
    #this for loop is invoked to add "."for all other alleles for each character in the new reference nucleotide sequence
    #to preserve structural integrity
    for(k in 1:length(start[[loci[i]]])){
      if(nrow(nuc[[loci[[i]]]][start[[loci[i]]][k]:end[[loci[i]]][k],])!=nrow(nuc[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],])){
        x<-as.data.frame(nuc[[loci[i]]][,1][start[[loci[i]]][1]:end[[loci[i]]][1]][-c(1,2)], stringsAsFactors = F)
        colnames(x)<-paste(loci[[i]], "alleles", sep="_")
        x<-cbind.data.frame(x, pepseq=as.character(paste(rep(".", nchar(tail(nuc[[loci[i]]][,2], 1))), collapse = "")), stringsAsFactors=FALSE)
        y<-data.frame(tail(nuc[[loci[i]]], (nrow(nuc[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],][nrow(nuc[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],])!=nrow(nuc[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), stringsAsFactors = F)
        colnames(y)<-c(paste(loci[[i]], "alleles", sep="_"), "pepseq")
        x$pepseq[match(y[,1], x[,1])]<-y$pepseq
        nuc[[loci[i]]]<-as.matrix(rbind(head(nuc[[loci[i]]], -(nrow(nuc[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],][nrow(nuc[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],])!=nrow(nuc[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), x))
        start[[loci[i]]]<-as.numeric(grep("cDNA", nuc[[loci[i]]]))
        end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,nrow(nuc[[loci[i]]])))}
    }
    
    #subsets nucleotide alignment based on start and end blocks
    for(e in 1:length(start[[loci[[i]]]])){
      nuc_extract[[loci[i]]]<-cbind(nuc_extract[[loci[i]]], nuc[[loci[i]]][start[[loci[i]]][e]:end[[loci[i]]][e],])}
    
    
    #removes first two rows containing AA position and "Prot"
    nuc_extract[[loci[i]]] <- nuc_extract[[loci[i]]][-c(1,2,3),]
    
    #extracts only the nucleotide reference sequence
    extract_ref[[loci[i]]]<-rbind(nuc_extract[[loci[[i]]]][1,])
    
    #pastes every two columsn together to get full nucleotide reference sequence without allele repetition
    cols<-seq(0, ncol(extract_ref[[loci[i]]]), by=2)
    extract_ref[[loci[i]]]<-cbind(extract_ref[[loci[i]]][,1], apply(extract_ref[[loci[i]]][,cols, drop=F], 1 ,paste, collapse = ""))
    
    #split reference sequences at pipes
    pipe_split[[loci[[i]]]]<-strsplit(extract_ref[[loci[[i]]]][,2], "|", fixed = T)
    
    
    for(z in 1:(length(pipe_split[[loci[[i]]]][[1]])-1)){
      if(StrLeft(pipe_split[[loci[[i]]]][[1]][z], 2)!="..") {
        pipe_split[[loci[[i]]]][[1]][z]<-paste(pipe_split[[loci[[i]]]][[1]][z], StrLeft(pipe_split[[loci[[i]]]][[1]][z+1], 2), sep="")
      }
      
      if(StrLeft(pipe_split[[loci[[i]]]][[1]][z], 2)=="..") {
        pipe_split[[loci[[i]]]][[1]][z]<-pipe_split[[loci[[i]]]][[1]][z]
      }
      
      if(StrLeft(pipe_split[[loci[[i]]]][[1]][z+1], 2)=="..") {next}
      pipe_split[[loci[[i]]]][[1]][z+1]<-str_replace(pipe_split[[loci[[i]]]][[1]][z+1], substr(pipe_split[[loci[[i]]]][[1]][z+1],1, 2), "")
    } 
    #splits every nucleotide within each boundary, and removes InDels
    for(k in 1:length(pipe_split[[loci[[i]]]][[1]])){
      boundary_split[[loci[[i]]]][[k]]<-strsplit(pipe_split[[loci[[i]]]][[1]][[k]], "*")[[1]]
      if(boundary_split[[loci[[i]]]][[k]][[length(boundary_split[[loci[[i]]]][[k]])]]=="."){next}
      if(boundary_split[[loci[[i]]]][[k]][1]=="."){next}
      boundary_split[[loci[[i]]]][[k]]<-boundary_split[[loci[[i]]]][[k]][boundary_split[[loci[[i]]]][[k]]!="."]
    }
    
    
    #forms a blank dataframe with the number of rows equal to the number of exon boundaries 
    atlas[[loci[[i]]]]<-data.frame(matrix("", ncol=1, nrow =(length(boundary_split[[loci[[i]]]])-1)), stringsAsFactors = F)
    
    colnames(atlas[[loci[[i]]]])<-"Exon"
    
    #pastes together individual nucleotides to form peptides, counts number of peptides 
    #present between boundaries
    for(q in 1:(length(boundary_split[[loci[[i]]]])-1)){
      j <- seq.int(1L,length(boundary_split[[loci[[i]]]][[q]]),by = 3L)
      boundaries[[loci[[i]]]][[q]]<-length(paste0(boundary_split[[loci[[i]]]][[q]][j],boundary_split[[loci[[i]]]][[q]][j+1], boundary_split[[loci[[i]]]][[q]][j+2]))
    }
    
    #breaks out of for loop to add alignment start enumeration to actual start
    boundaries[[loci[[i]]]][[1]]<-boundaries[[loci[[i]]]][[1]]+pep_start[[loci[[i]]]]
    
    positions<-sapply(length(boundary_split[[loci[[i]]]])-1, function(x) NULL)
    
    #fills in atlas information
    #boundaries obtained by adding up cumulative lengths
    for(q in 1:(length(boundary_split[[loci[[i]]]])-1)){
      atlas[[loci[[i]]]][q,1]<-paste("exon", paste(seq(1, length(boundary_split[[loci[[i]]]]))[[q]], seq(1, length(boundary_split[[loci[[i]]]]))[[q+1]], sep=":"), sep="_")
      positions[[q]]<-as.numeric((0 + cumsum(boundaries[[loci[[i]]]])[[q]]))}
    atlas[[loci[[i]]]]$Boundary<-positions
  } 
  return(atlas)
}

#usage
AA_atlas<-atlasMaker(c('A','B','C','DMA','DMB','DOA','DOB','DPA1','DPB1','DQA1','DQA2','DQB1','DRA','DRB1','E','F','G','HFE','MICA','MICB','TAP1','TAP2'))

save(AA_atlas, file="AA_atlas.rda")



