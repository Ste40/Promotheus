
TSSparsing<- function(GSEdirectory, strandness, Acclist, GSE){
    
  setwd(GSEdirectory)
  
  if(length(list.files(pattern="*_starts_cvg.txt")) != length(Acclist)){
    
    bamfiles<- list.files(pattern= "*1.sorted.bam$", recursive=TRUE)
    bamfiles<- normalizePath(bamfiles)
    bamfiles<- sub(".*/", "", bamfiles)
    
    bampath<- bamfiles
    bampath<- gsub("E:\\", "/mnt/e/", bampath, fixed=TRUE)
    bampath<- gsub("\\", "/", bampath, fixed=TRUE)
    bampath<- gsub(" ", "\\\ ", bampath, fixed=TRUE)
    
    if(strandness[[1]][1] == "This is SingleEnd Data"){
      
      if(strandness[[1]][2] < 0.8){  #reversely single end
        
        message(paste("The inferred library for", GSE, "is Single End reversely stranded. I'm preparing files for TSS identification, please wait"))
        
        for(i in 1:length(bampath)){
          
          bamfile<- bampath[i]
          string_sorted<- bamfile
          
          string<- paste('samtools index', string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          negativo<- paste0(gsub("_sorted.bam", "", string_sorted),"_negative.bam")
          string<- paste('samtools view -F 16 -b -o', negativo, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', negativo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          positivo<- paste0(gsub("_sorted.bam", "", string_sorted),"_positive.bam")
          string<- paste('samtools view -f 16 -b -o', positivo, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', positivo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string_positive_cvg<- paste0(gsub(".bam", "", positivo),'.cvg.txt')
          string<- paste('samtools depth -aa -o', string_positive_cvg, positivo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string_negative_cvg<- paste0(gsub(".bam", "", negativo),'.cvg.txt')
          string<- paste('samtools depth -aa -o', string_negative_cvg, negativo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools view', positivo, '| cut -f 4 >', paste0(bamfile,"_positive.sam"))
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools view', negativo, '| cut -f 4 >', paste0(bamfile,"_negative.sam"))
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          tab<- read.csv(paste0(bamfiles[i],"_positive.sam"), sep="\t", header=F)
          
          POS<- tab$V1
          lengthref<- 1:nrow(read.csv(str_replace(bamfiles[i], ".bam", "_negative.cvg.txt")))*0
          reference<- c(0,lengthref)
          
          for(j in 1:length(POS)){
            reference[POS[j]]<- reference[POS[j]]+1
            
          }
          
          #write.table(reference, paste0(gsub(".bam", "", bamfiles[i]),"_positive_starts.txt"), sep="\t", row.names=TRUE, quote=FALSE)
          tabella<-data.frame(reference)
          tab<- read.csv(paste0(bamfiles[i],"_negative.sam"), sep="\t", header=F)
          
          POS<- tab$V1
          
          #string<- paste('samtools view ', negativo, '| awk \'{print length($10)}\' | head -1000 | sort -u')
          string<- paste('bedtools bamtobed -i', negativo ,'| awk \'{print $3}\' > negative_bed.txt')
          bash<- paste('bash -c', shQuote(string))
          system(bash, intern=TRUE)
          lengthreads<-read.csv("negative_bed.txt", sep="\t", header=F)
          lengthreads<-as.vector(lengthreads[,1])
          
          POS<- lengthreads
          #POS<- rev(POS)
          reference<- c(0,lengthref)
          
          for(j in 1:length(POS)){
            reference[POS[j]]<- reference[POS[j]]+1
            }
          
          #write.table(reference, paste0(gsub(".bam", "", bamfiles[i]),"_negative_starts.txt"), sep="\t", row.names=TRUE, quote=FALSE)
          #ref_revs<-rev(reference)
          cvg_positivo<-read.table(str_replace(bamfiles[i], ".bam", "_positive.cvg.txt"), sep="\t",header=FALSE)
          cvg_negativo<-read.table(str_replace(bamfiles[i], ".bam", "_negative.cvg.txt"),sep="\t",header=FALSE)
          pos<-c(1:length(cvg_positivo[,1]))
          tabella<-data.frame(pos, tabella, cvg_positivo[,3], reference,  cvg_negativo[,3])
          colnames(tabella)<-c("x","positive_starts", "positive_cvg", "negative_starts", "negative_cvg")
          write.table(tabella, paste0(gsub(".bam", "", bamfiles[i]),"_starts_cvg.txt"), sep="\t", row.name=FALSE, quote=FALSE)
        }
        
      } else { #single end stranded
        
        message(paste("The inferred for", GSE, "is Single End stranded/unstranded. I'm preparing files for TSS identification, please wait"))
        
        for(i in 1:length(bampath)){
          
          bamfile<- bampath[i]
          string_sorted<- bamfile
          
          string<- paste('samtools index', string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          positivo<- paste0(gsub("_sorted.bam", "", string_sorted),"_positive.bam")
          string<- paste('samtools view -F 16 -b -o', positivo, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', positivo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          negativo<- paste0(gsub("_sorted.bam", "", string_sorted),"_negative.bam")
          string<- paste('samtools view -f 16 -b -o', negativo, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', negativo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string_positive_cvg<- paste0(gsub(".bam", "", positivo),'.cvg.txt')
          string<- paste('samtools depth -aa -o', string_positive_cvg, positivo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string_negative_cvg<- paste0(gsub(".bam", "", negativo),'.cvg.txt')
          string<- paste('samtools depth -aa -o', string_negative_cvg, negativo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools view', positivo, '| cut -f 4 >', paste0(bamfile,"_positive.sam"))
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools view', negativo, '| cut -f 4 >', paste0(bamfile,"_negative.sam"))
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          tab<- read.csv(paste0(bamfiles[i],"_positive.sam"), sep="\t", header=F)
          
          POS<- tab$V1
          lengthref<- 1:nrow(read.csv(str_replace(bamfiles[i], ".bam", "_negative.cvg.txt")))*0
          reference<- c(0,lengthref)
          
          for(j in 1:length(POS)){
            
            reference[POS[j]]<- reference[POS[j]]+1
            
          }
          
          #write.table(reference, paste0(gsub(".bam", "", bamfiles[i]),"_positive_starts.txt"), sep="\t", row.names=TRUE, quote=FALSE)
          tabella<-data.frame(reference)
          tab<- read.csv(paste0(bamfiles[i],"_negative.sam"), sep="\t", header=F)
          
          POS<- tab$V1
          
          #string<- paste('samtools view ', negativo, '| awk \'{print length($10)}\' | head -1000 | sort -u')
          string<- paste('bedtools bamtobed -i', negativo ,'| awk \'{print $3}\' > negative_bed.txt')
          bash<- paste('bash -c', shQuote(string))
          system(bash, intern=TRUE)
          lengthreads<-read.csv("negative_bed.txt", sep="\t", header=F)
          lengthreads<-as.vector(lengthreads[,1])
          
          POS<- lengthreads
          #POS<- rev(POS)
          reference<- c(0,lengthref)
          
          
          for(j in 1:length(POS)){
            reference[POS[j]]<- reference[POS[j]]+1
            }
          
          #write.table(reference, paste0(gsub(".bam", "", bamfiles[i]),"_negative_starts.txt"), sep="\t", row.names=TRUE, quote=FALSE)
          #ref_revs<-rev(reference)
          cvg_positivo<-read.table(str_replace(bamfiles[i], ".bam", "_positive.cvg.txt"), sep="\t",header=FALSE)
          cvg_negativo<-read.table(str_replace(bamfiles[i], ".bam", "_negative.cvg.txt"),sep="\t",header=FALSE)
          pos<-c(1:length(cvg_positivo[,1]))
          tabella<-data.frame(pos, tabella, cvg_positivo[,3], reference,  cvg_negativo[,3])
          colnames(tabella)<-c("x","positive_starts", "positive_cvg", "negative_starts", "negative_cvg")
          write.table(tabella, paste0(gsub(".bam", "", bamfiles[i]),"_starts_cvg.txt"), sep="\t", row.name=FALSE, quote=FALSE)
        }
        
      }
    } 
    
    
    #Single end
    #######################################
    #######################################
    #Paired end
    
    
    else {
      
      if(strandness[[1]][2] < 0.8){ #reversely paired end
        
        message(paste("The inferred library for", GSE, "is Paired End reversely stranded, 
                      I'm preparing files for TSS identification, please wait"))
        
        for(i in 1:length(bampath)){
          
          bamfile<- bampath[i]
          string_sorted<- bamfile
          
          string<- paste('samtools index', string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          negativo1<- paste0(gsub("_sorted.bam", "", string_sorted),"_negative_1.bam")
          string<- paste('samtools view -f 144 -b -o', negativo1, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', negativo1)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          negativo2<- paste0(gsub("_sorted.bam", "", string_sorted),"_negative_2.bam")
          string<- paste('samtools view -f 64 -F 16 -b -o', negativo2, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', negativo2)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          negativo<- paste0(gsub("_sorted.bam", "", string_sorted),"_negativo.bam")
          string<- paste('samtools merge -f', negativo, negativo1, negativo2)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', negativo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          positivo1<- paste0(gsub("_sorted.bam", "", string_sorted),"_positive_1.bam")
          string<- paste('samtools view -f 128 -F16 -b -o', positivo1, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', positivo1)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          positivo2<- paste0(gsub("_sorted.bam", "", string_sorted),"_positive_2.bam")
          string<- paste('samtools view -f 80 -b -o', positivo2, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', positivo2)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          positivo<- paste0(gsub("_sorted.bam", "", string_sorted),"_positivo.bam")
          string<- paste('samtools merge -f', positivo, positivo1, positivo2)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', positivo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string_positive_cvg<- paste0(gsub(".bam", "", positivo),'.cvg.txt')
          string<- paste('samtools depth -aa -o', string_positive_cvg, positivo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string_negative_cvg<- paste0(gsub(".bam", "", negativo),'.cvg.txt')
          string<- paste('samtools depth -aa -o', string_negative_cvg, negativo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools view', positivo, '| cut -f 4 >', paste0(bamfile,"_positive.sam"))
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools view', negativo, '| cut -f 4 >', paste0(bamfile,"_negative.sam"))
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          tab<- read.csv(paste0(bamfiles[i],"_positive.sam"), sep="\t", header=F)
          
          POS<- tab$V1
          lengthref<- 1:nrow(read.csv(str_replace(bamfiles[i], ".bam", "_negativo.cvg.txt")))*0
          reference<- c(0,lengthref)
          
          for(j in 1:length(POS)){ 
            reference[POS[j]]<- reference[POS[j]]+1
            
          }
          
          #write.table(reference, paste0(gsub(".bam", "", bamfiles[i]),"_positive_starts.txt"), sep="\t", row.names=TRUE, quote=FALSE)
          tabella<-data.frame(reference)
          #tab<- read.csv(paste0(bamfiles[i],"_negative.sam"), sep="\t", header=F)
          
          #POS<- tab$V1
          
          #string<- paste('samtools view ', negativo, '| awk \'{print length($10)}\' | head -1000 | sort -u')
          string<- paste('bedtools bamtobed -i', negativo ,'| awk \'{print $3}\' > negative_bed.txt')
          bash<- paste('bash -c', shQuote(string))
          system(bash, intern=TRUE)
          lengthreads<-read.csv("negative_bed.txt", sep="\t", header=F)
          lengthreads<-as.vector(lengthreads[,1])
          
          POS<- lengthreads
          #POS<- rev(POS)
          reference<- c(0,lengthref)
          
          for(j in 1:length(POS)){
            reference[POS[j]]<- reference[POS[j]]+1
            }
          
          #write.table(reference, paste0(gsub(".bam", "", bamfiles[i]),"_negative_starts.txt"), sep="\t", row.names=TRUE, quote=FALSE)
          #ref_revs<-rev(reference)
          cvg_positivo<-read.table(str_replace(bamfiles[i], ".bam", "_positivo.cvg.txt"), sep="\t",header=FALSE)
          cvg_negativo<-read.table(str_replace(bamfiles[i], ".bam", "_negativo.cvg.txt"),sep="\t",header=FALSE)
          pos<-c(1:length(cvg_positivo[,1]))
          tabella<-data.frame(pos, tabella, cvg_positivo[,3], reference,  cvg_negativo[,3])
          colnames(tabella)<-c("x","positive_starts", "positive_cvg", "negative_starts", "negative_cvg")
          write.table(tabella, paste0(gsub(".bam", "", bamfiles[i]),"_starts_cvg.txt"), sep="\t", row.name=FALSE, quote=FALSE)
        }
        
      } else { 
        
        message(paste("The inferred for", GSE, "is Paired End stranded/unstranded. I'm preparing files for TSS identification, please wait"))
        
        for(i in 1:length(bampath)){
          
          bamfile<- bampath[i]
          string_sorted<- bamfile
          
          string<- paste('samtools index', string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          negativo1<- paste0(gsub("_sorted.bam", "", string_sorted),"_negative_1.bam")
          string<- paste('samtools view -f 144 -b -o', negativo1, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', negativo1)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          negativo2<- paste0(gsub("_sorted.bam", "", string_sorted),"_negative_2.bam")
          string<- paste('samtools view -f 64 -F 16 -b -o', negativo2, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', negativo2)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          negativo<- paste0(gsub("_sorted.bam", "", string_sorted),"_negativo.bam")
          string<- paste('samtools merge -f', negativo, negativo1, negativo2)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', negativo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          positivo1<- paste0(gsub("_sorted.bam", "", string_sorted),"_positive_1.bam")
          string<- paste('samtools view -f 128 -F16 -b -o', positivo1, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', positivo1)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          positivo2<- paste0(gsub("_sorted.bam", "", string_sorted),"_positivo_2.bam")
          string<- paste('samtools view -f 80 -b -o', positivo2, string_sorted)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', positivo2)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          positivo<- paste0(gsub("_sorted.bam", "", string_sorted),"_positivo.bam")
          string<- paste('samtools merge -f', positivo, positivo1, positivo2)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools index', positivo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string_positive_cvg<- paste0(gsub(".bam", "", positivo),'.cvg.txt')
          string<- paste('samtools depth -aa -o', string_positive_cvg, positivo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string_negative_cvg<- paste0(gsub(".bam", "", negativo),'.cvg.txt')
          string<- paste('samtools depth -aa -o', string_negative_cvg, negativo)
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools view', positivo, '| cut -f 4 >', paste0(bamfile,"_positive.sam"))
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          string<- paste('samtools view', negativo, '| cut -f 4 >', paste0(bamfile,"_negative.sam"))
          bash<- paste('bash -c', shQuote(string))
          system(bash)
          
          tab<- read.csv(paste0(bamfiles[i],"_positive.sam"), sep="\t", header=F)
          
          POS<- tab$V1
          lengthref<- 1:nrow(read.csv(str_replace(bamfiles[i], ".bam", "_negativo.cvg.txt")))*0
          reference<- c(0,lengthref)
          
          for(j in 1:length(POS)){ 
            
            reference[POS[j]]<- reference[POS[j]]+1
            
          }
          
          #write.table(reference, paste0(gsub(".bam", "", bamfiles[i]),"_positive_starts.txt"), sep="\t", row.names=TRUE, quote=FALSE)
          tabella<-data.frame(reference)
          tab<- read.csv(paste0(bamfiles[i],"_negative.sam"), sep="\t", header=F)
          
          #POS<- tab$V1
          
          #string<- paste('samtools view ', negativo, '| awk \'{print length($10)}\' | head -1000 | sort -u')
          string<- paste('bedtools bamtobed -i', negativo ,'| awk \'{print $3}\' > negative_bed.txt')
          bash<- paste('bash -c', shQuote(string))
          system(bash, intern=TRUE)
          lengthreads<-read.csv("negative_bed.txt", sep="\t", header=F)
          lengthreads<-as.vector(lengthreads[,1])
          
          POS<- lengthreads
          #POS<- rev(POS)
          reference<- c(0,lengthref)
          
          for(j in 1:length(POS)){
            reference[POS[j]]<- reference[POS[j]]+1
            }
          
          #write.table(reference, paste0(gsub(".bam", "", bamfiles[i]),"_negative_starts.txt"), sep="\t", row.names=TRUE, quote=FALSE)
          #ref_revs<-rev(reference)
          cvg_positivo<-read.table(str_replace(bamfiles[i], ".bam", "_positivo.cvg.txt"), sep="\t",header=FALSE)
          cvg_negativo<-read.table(str_replace(bamfiles[i], ".bam", "_negativo.cvg.txt"),sep="\t",header=FALSE)
          pos<-c(1:length(cvg_positivo[,1]))
          tabella<-data.frame(pos, tabella, cvg_positivo[,3], reference,  cvg_negativo[,3])
          colnames(tabella)<-c("x","positive_starts", "positive_cvg", "negative_starts", "negative_cvg")
          write.table(tabella, paste0(gsub(".bam", "", bamfiles[i]),"_starts_cvg.txt"), sep="\t", row.name=FALSE, quote=FALSE)
        }  
        
      }
      
    }
    
    if(length(list.files(pattern="*_starts_cvg.txt")) != length(Acclist)){
      token=FALSE
      
    } else {token=TRUE}
    
  } else { token=TRUE }

return(token)

}
