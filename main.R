#Tool for analysis of Genome Quebec sanger sequencing results:
#Kai Ellis 2020

#Minimum requirements: a folder named Ab1 with .ab1 sanger sequencing files in it
#Optional folders:
#plate_info with one or more genome quebec plate data .csv files
#   this file will result in ab1 files being automatically filtered based on Genome quebec quality control results

#References with one or more .txt files of reference genes in fasta format
#   This folder will enable optional alignment of reference sequences with the imported sanger seq results

#metadata with one or more .csv files of sample metadata
#   This folder will enable optional output of combined metadat/sequencing result csvs
#   NOTE: the metadata .csvs must include a column with your full sample names (the same names that were sent in the 
#   genome quebec request excel sheet)
s
#Functions:
testFile <- function(x) {
  out <- tryCatch({
    
    message("Testing input")
    
    setwd(paste0(x))
  },
  error=function(cond) {
    message(paste("Path does not seem to exist:", x))
    message("Here's the original error message:")
    message(cond)
    return(NA)
  }
  )
  return(out)
  
}

#Libraries
if("seqinr" %in% (.packages())){
  detach("package:seqinr", unload=TRUE)
}

if(!("BiocManager" %in% (.packages()))){
  library(BiocManager)
  library(sangerseqR)
  library(ape)
  library(dplyr)
}

#part 1 of the program, asks the user for a working directory and either sets it or allows the user to try a 
#diffrent input

P0 = T
while (P0 == T){
  P0 = F
  message("Current working directory is ", getwd(), "\n")
  current_files = list.files(getwd())
  
  message("file contains:")
  
  print(current_files)
  numfiles = length(list.files(path = paste0(getwd(),"\\Ab1"), pattern="*.\\.ab1", full.names = T, recursive = F))
  print(paste("With ", numfiles, " .ab1 files in the Ab1 folder" ))
  
  test <- toupper(readline(prompt = "Would you like to change it? (Y/N), enter [B] to quit: "))
  if(test == "Y"){
    P1 = T
    rm(test)
    next
  } 
  
  if(test == "N"){
    P1 = F
    rm(test)
    core_files <- list.files(getwd())
    next
  } 
  
  if(test == "B"){
    rm(test)
    break
  } else {
    message("invalid entry, please try again\n")
    P0 = T
  }
  
}

while (P1 == T){
  Path <- readline(prompt = "Insert file path, or enter [B] to quit: ")
  if (toupper(Path) == "B"){
    break
  }
  
  test <- testFile(Path)
  
  if (is.na(test)){
    next
  }
  
  core_files = list.files(paste(Path))
  
  message("file contains:")
  
  print(core_files)
  numfiles = length(list.files(path = paste0(Path,"\\Ab1"), pattern="*.\\.ab1", full.names = T, recursive = F))
  print(paste("With ", numfiles, " .ab1 files in the Ab1 folder" ))
  
  
  
  if(toupper(readline(prompt = "This is the correct file? (Y/N) ")) == "Y"){
    cat("Working directory set")
    P1 = F
    
  }
}

#Part 2 of the program, which asks the user which trimming methodology they want to use
#to alter their imported ab1 files and outputs it as a dataframe for either output or alignment

P2 = T
while (P2 == T) {
  
  message("Import all files in the .ab1 file with which type of trimming?")
  menu <- readline(prompt = cat(" [1] Bulk Trim \n [2] Manual Trim \n [3] No Trim \n [9] Help \n [B] Quit"))
  
  if (menu == "1"){
    Bulk = T
    
    T5 <- as.numeric(readline(prompt = "How many bases would you like trimmed from the 5' end?  "))
    T3 <- as.numeric(readline(prompt = "What position would you like the 3' end truncated at?  "))
    
    files <- list.files(path = "./Ab1", pattern="*.\\.ab1", full.names = T, recursive = F)
    class(paste0(T3))
    T3
    #get info from the plate_info file if provided
    if("plate_info" %in% core_files){
      message("plate_info file detected; only reading .ab1 files that passed GQ QC")
      #first we read all csvs in the plate_info file and bind them together
      temp <- list.files(path = "./plate_info", pattern="*.csv", full.names = T)
      fullDat <- lapply(temp, read.csv ) %>% bind_rows(.id = "Submission.lot")
      #then we create a test variable with all files that passed qc connected by the symbol for OR (|)
      qcPass <- fullDat %>%  filter(Ok=="true"|Ok==TRUE) 
      test <- paste(qcPass$Chromatogram, collapse = "|")
      
      #finally, we can read through the list of available files, reading in only those which passed QC
      PrimarySeq = lapply (files, function(x) {
        Pass = grepl(test, basename(x))
        
        if (Pass == TRUE) {
          ReadFile <- read.abif(x)
          Readable <- sangerseq(ReadFile)
          Seq1 = paste(Readable@primarySeq)
        }
      })
      
      #this returns the name of successful sequencing reactions, or null for unsuccessfull ones
      PrimaryNames = lapply (files, function(x) {
        Pass = grepl(test, basename(x))
        
        if (Pass == TRUE) {
          Seq1 = paste(basename(x))
        }
      })
    }
    
    
    #this returns the sequnce of successful sequencing reactions, or null for unsuccessfull ones
    
    #drop all nulls and bind into a data frame
    PrimaryNames <- PrimaryNames[-which(sapply(PrimaryNames, is.null))]
    PrimarySeq <- PrimarySeq[-which(sapply(PrimarySeq, is.null))]
    Sequences <- do.call(rbind, Map(data.frame, Label = PrimaryNames, Primer = PrimaryNames,    Sequence=PrimarySeq))
    
    numseq = nrow(Sequences)
    
    message(paste0(numseq, " .ab1files read"))
    
    #cut down to size and allocate primers and labels into a copy of the dataframe
    SequencesF <- data.frame(Sequences)
    SequencesF$Primer <- sub("_.*", "",(sub("^(.*?)_", "", Sequences$Primer)))
    SequencesF$Label <- sub("_.*", "", Sequences$Label)
    SequencesF$Sequence <- sapply(as.character(Sequences$Sequence), strtrim, T3)
    SequencesF$Sequence  <-  sub(paste0(".{",T5,"}"),"",SequencesF$Sequence)
    
    #ask user if they like the cut of the gib
    while(Bulk == T){
      message("Here are the first 6 sequences of the ", numseq, " read")
      
      print(head(SequencesF))
      
      #if not, retrim
      bulkmenu <- toupper(readline(prompt = "Trim differently? (Y/N), or enter [B] to return to previous trim menu  "))
      if(bulkmenu == "B"){
        rm(bulkmenu)
        break
      }
      
      if(bulkmenu == "Y"){
        T5 <- as.numeric(readline(prompt = "How many bases would you like trimmed from the 5' end?  "))
        T3 <- as.numeric(readline(prompt = "What position would you like the 3' end truncated at?  "))
        
        SequencesF$Sequence <- sapply(as.character(Sequences$Sequence), strtrim, T3)
        SequencesF$Sequence  <-  sub(paste0(".{",T5,"}"),"",SequencesF$Sequence)
        
      }
      
      if(bulkmenu == "N"){
        rm(Sequences)
        rm(bulkmenu)
        Bulk = F
        P2 = F
      }
      
      
    }
  }
  
  if (menu == "2"){
    files <- list.files(path = "./Ab1", pattern="*.\\.ab1", full.names = T, recursive = F)
    
    if("plate_info" %in% core_files){
      message("plate_info file detected; only reading .ab1 files that passed GQ QC")
      #first we read all csvs in the plate_info file and bind them together
      temp <- list.files(path = "./plate_info", pattern="*.csv", full.names = T)
      fullDat <- lapply(temp, read.csv ) %>% bind_rows(.id = "Submission.lot")
      #then we create a test variable with all files that passed qc connected by the symbol for OR (|)
      qcPass <- fullDat %>%  filter(Ok=="true"|Ok==TRUE) 
      test <- paste(qcPass$Chromatogram, collapse = "|")
      
      #test all files and replace with NA's
      for(i in 1:length(files)){
        pass <- grepl(test, basename(files[i]))
        if (pass==F){
          files[i] <- NA
        }
        files <- files[complete.cases(files)]
      }
      
    }
    
    
    
    df <- data.frame(Label = character(),
                     Primer = character(),
                     Sequence = character(),
                     fivetrim = numeric(),
                     threetrim = numeric(),
                     stringsAsFactors=FALSE)
    
    i=1
    noglobe = T
    readline(prompt= "please expand your plot viewer, or you will get an error.\n Press [Enter] to view the first chromatogram ")
    while (i < length(files)+1){
      #main while loop, generates a dataframe for all df. df that are skipped are NA'd. df are trimmed after
      #based upon trim values
      
      message(paste0("you are on sequence ", i, " of ", length(files)))
      
      tmp <- readsangerseq(files[i])
      
      if(noglobe == T){
        tr3 <- length(tmp@primarySeq)
        tr5 <- 0
      }
      
      V3 <- length(tmp@primarySeq)-tr3
      wid <- ceiling((length(tmp@primarySeq)-V3-tr5)/8)
      
      chromatogram(tmp, width = wid, trim3= V3, trim5 = tr5)
      
      manMenu <- readline(prompt = cat(" [1] Alter global filter \n [2] Trim sequence \n [3] Skip sequence \n [4] Skip remaining sequences \n [5] Help \n [9] Restart manual trim \n [B] Quit"))
      
      if(manMenu == "1"){
        tr5 <- as.numeric(readline(prompt = "Which position would you like all chromatograms to start at? "))-1
        tr3 <- as.numeric(readline(prompt = "Which position would you like all chromatograms to end at? "))
        noglobe = F
      }
      
      if(manMenu == "2"){
        df[i,1] <- paste(basename(files[i]))
        df[i,2] <- paste(basename(files[i]))
        df[i,3] <- paste(tmp@primarySeq)
        df[i,4] <- as.numeric(readline(prompt = "How many bases do you want trimmed from the 5' end? "))
        df[i,5] <- as.numeric(readline(prompt = "What base do you the 3' end truncated at? "))
        
        i = i+1
      }
      
      if(manMenu == "3"){
        df[i,1] <- NA 
        df[i,2] <- NA 
        df[i,3] <- NA 
        df[i,4] <- NA
        df[i,5] <- NA
        i = i+1
      }
      
      if(manMenu == "4"){
        df[i,1] <- NA 
        df[i,2] <- NA 
        df[i,3] <- NA 
        df[i,4] <- NA
        df[i,5] <- NA
        i = length(files)+1}
      
      if(manMenu == "5"){
        message(cat("Options for manual trim: \n\n
                    [1] alter global filter: adjusts the range of base numbers that will be viewed on the current, \n
                    and future chromatograms. Ex: inserting values of 1 and 500 will ensure all chromatograms only include \n
                    bases between the first and 500th. Values can be altered as many times as you want.\n\n
                    [2] Trim sequence: provide an initial 5' base, and a termial 3' base \n
                    any bases before the 5' and after the 3' will be trimmed.\n\n
                    [3] Skip sequence: This sequence will not be imported or trimmed. \n\n
                    [4] Skip remaining sequences: All remaining sequences will not be imported or trimmed\n\n
                    [9] Restart: removes all trimming and jumps back to the first provided ab1 file
        "))
        readline(prompt = "press [enter] to return to main menu")
      }
      
      if(manMenu == "9"){
        df <- data.frame(name = character(),
                         primer = character(),
                         seq = character(),
                         fivetrim = numeric(),
                         threetrim = numeric(),
                         stringsAsFactors=FALSE)
        i=1
        
      }
      
      if (toupper(manMenu) == "B"){
        
        rm(i,tr3,tr5,V3,wid,df)
        break
      }
      
      if( i >=length(files)){
        df <- df[complete.cases(df),]
        
        for (i in 1:nrow(df)){
          df[i,1] <- sub("_.*", "", df[i,1])
          df[i,2] <- sub("_.*", "",(sub("^(.*?)_", "", df[i,2])))
          df[i,3] <- strtrim(df[i,3],df[i,5]) 
          df[i,3]  <-  sub(paste0(".{",df[i,4],"}"),"",df[i,3])
        }
        SequencesF <- df
        P2 <- F
        
      }
    }
  }
  
  if (menu == "3"){
    files <- list.files(path = "./Ab1", pattern="*.\\.ab1", full.names = T, recursive = F)
    
    if("plate_info" %in% core_files){
      message("plate_info file detected; only reading .ab1 files that passed GQ QC")
      #first we read all csvs in the plate_info file and bind them together
      temp <- list.files(path = "./plate_info", pattern="*.csv", full.names = T)
      fullDat <- lapply(temp, read.csv ) %>% bind_rows(.id = "Submission.lot")
      #then we create a test variable with all files that passed qc connected by the symbol for OR (|)
      qcPass <- fullDat %>%  filter(Ok=="true"|Ok==TRUE) 
      test <- paste(qcPass$Chromatogram, collapse = "|")
      for(i in 1:length(files)){
        pass <- grepl(test, basename(files[i]))
        if (pass==F){
          files[i] <- NA
        }
      }
      
      
    }
    
    files <- files[complete.cases(files)]
    PrimarySeq = lapply (files, function(x) {
      ReadFile <- read.abif(x)
      Readable <- sangerseq(ReadFile)
      Seq1 = paste(Readable@primarySeq)
    })
    
    #this returns the name of successful sequencing reactions, or null for unsuccessfull ones
    PrimaryNames = lapply (files, function(x) {
      Seq1 = paste(basename(x))
    })
    
    #bind into a data frame
    SequencesF <- do.call(rbind, Map(data.frame, Label = PrimaryNames, Primer = PrimaryNames, Sequence=PrimarySeq))
    
    numseq = nrow(SequencesF)
    
    message(paste0(numseq, " .ab1files read"))
    
    #cut down to size and allocate primers and labels into a copy of the dataframe
    SequencesF$Primer <- sub("_.*", "",(sub("^(.*?)_", "", SequencesF$Primer)))
    SequencesF$Label <- sub("_.*", "", SequencesF$Label)
    
    message("Here are the first 6 sequences of the ", numseq, " read")
    
    print(head(SequencesF))
    
    #if not, retrim
    nomenu <- toupper(readline(prompt = "Trim sequences? (Y/N)? "))
    if(nomenu == "N"){
      P2=F
      next
    }
    if(nomenu == "Y"){
      rm(nomenu)
      next
    }
  }
  
  if (toupper(menu) == "B"){
    break
  }
  
  if (menu == "9"){
    cat("Options for trimming: \n\n
        Bulk will trim all imported sequences to the input 3' and 5' values.
        This means if you enter 50 and 300 for the  3' and 5' values all bases 
        before the 50th and after the 300th will be trimmed.Good for quick 
        and dirty analysis of large datasets or for if you don't require precise trimming\n
        \n.
        Manual will display the chromatogram for each sequence, 
        and ask how many bases you want trimmed off each end.
        Good for small datasets or if you want more precise control \n\n
        None will simply import the sequences with no trimming.
        ")
    readline(prompt = "press [enter] to continue")
  } 
}

if(P2 == F){
  P3 = T
}

while (P3 == T){
  message("Sequences imported, what next?")
  menu = toupper(readline(prompt=" [1] Align sequences \n [2] Output sequences \n [9] Help \n [B] Quit"))
  if (menu == "1"){
    
    bin <- sapply(as.character(SequencesF$Sequence), strsplit, split = "")
    
    names(bin)<- paste(1:length(bin), SequencesF$Label, sep = "_")
    
    if ("References" %in% core_files){
      message("References file detected, import reference sequences? ")
      refq <- toupper(readline(prompt ="(Y/N)"))
      if(refq == "Y"){
        temp <- list.files(path = "./References", pattern="*.txt", full.names = T)
        numtxt = grep("\\.txt", temp)
        if(sum(numtxt) >= 1){
          library(seqinr)
          metaDat <- lapply(temp, read.fasta ) %>% bind_rows(.id = "label")
          reference <- sapply(as.character(temp), strsplit, split="")
          names(reference) <- paste(1:length(temp), names(temp), sep="_")
          bin <- c(bin, reference)
        }else{
          readline(prompt = "No reference .*/.txt files were provided in reference folder, press [Enter] to move to alignment")
        }
      }
    }
    
    bin <- as.DNAbin(bin)
    
    Align <- muscle(bin, quiet = F)
    
    message("Sequences aligned, would you like to view alignment graphics?")
    viewMen <- toupper(readline(prompt = "(Y/N) \n or press [B] to return to previous menu" ))
    
    if(viewMen == "Y"){
      checkAlignment(Align)
      readline("press [Enter] to move on")
      P3 = F
      layout( matrix(1, ncol=1) )
      next
    }
    
    if(viewMen == "N"){
      P3 = F
      next
    }
    
    
    
  }
  
  if(menu == "2"){
    if(!("Outputs" %in% core_files)){
      dir.create("Outputs")
      message("No outputs folder detected, outputs folder created")
    }
    if("Metadata" %in% core_files){
      MetaMenu <- toupper(readline(prompt="Metadata folder detected, options for output: \n [1] Merge sequence data with metadata csv 
                                      \n [2] output trimmed sequences as fasta file
                                      \n [B] return to previous menu" ))
      if (MetaMenu == "1"){
        temp <- list.files(path = "./Metadata", pattern="*.csv", full.names = T)
        numcsv = grep("\\.csv", temp)
        if(sum(numcsv) >= 1){
          metaDat <- lapply(temp, read.csv ) %>% bind_rows(.id = "Csv")
          temp <- SequencesF
          y<- left_join(temp,metaDat, by = "Label")
          write.csv(y, file = "./Outputs/MetadataSeqs.csv")
          rm(temp,y)
          message("MetadataSeqs.csv can now be found in the outputs folder")
        }
      }
      if (MetaMenu == "2"){
        bin <- sapply(as.character(SequencesF$Sequence), strsplit, split = "")
        
        names(bin)<- paste(1:length(bin), SequencesF$Label, sep = "_")
        
        bin <- as.DNAbin(bin)
        
        write.dna(bin, file= "./Outputs/Trim.fas", format="fasta", nbcol=-1, colsep="")
        message("Trim.fas can now be found in the outputs folder")
      }else{
        next
      } 
    } else{write.dna(Align, file= "./Outputs/Align.fas", format="fasta", nbcol=-1, colsep="")
      message("Align.fas can now be found in the outputs folder")
    }
  }
  
  if (menu == "9"){
    message(cat("Options for imported sequences:\n\n
                 [1] Alignment: uses the MUSCLE algorithm to align imported DNA sequences. \n
                Additional Reference sequences can be imported here if a Reference folder is provided\n
                 [2] Output: A folder called outputs will be generated with your trimmed data file inside \n
                potential file types include: txt, fasta, and csv. \n
                Metadata can be included in the file generation if a Metadata folder is provided"))
    readline(prompt = "Press enter to return to main menu")
  }
  
  if (menu == "B"){
    break
  }
  
}
P4=F
if(P3 == F){
  P4 = T
}
while (P4 == T){
  message("Sequences aligned, what next?")
  finmenu = toupper(readline(prompt = cat(" [1] Generate exploratory phylogenetic tree \n [2] Output aligned data \n [B] Quit")))
  if(finmenu == "1"){
    library(phangorn)
    
    message("Generating NJ tree based upon a F81 distance matrix")
    
    phyObject <- phyDat(Align, type = "DNA")
    
    dna_dist <- dist.ml(phyObject, model="F81")
    dna_NJ <- NJ(dna_dist)
    plot(dna_NJ)
    
    message("Would you like to save this tree?")
    tsave <- toupper(readline(prompt = "(Y/N)"))
    if (tsave == "Y"){
      if(!("Outputs" %in% core_files)){
        dir.create("Outputs")
        message("No outputs file detected, file created")
      }
      
      pdf("./Outputs/prelimTree.pdf")
      
      plot(dna_NJ)
      
      dev.off()
      message("prelimTree.pdf can now be found in the outputs folder")
      
    }
  }
  if(finmenu == "2"){
    if(!("Outputs" %in% core_files)){
      dir.create("Outputs")
      message("No outputs folder detected, folder created")
    }
    if("Metadata" %in% core_files){
      MetaMenu <- toupper(readline(prompt=cat("Metadata folder detected, options for output: \n [1] Merge sequence data with metadata csv 
                                      \n [2] output aligned seqeuences as fasta file
                                      \n [B]return to previous menu" )))
      if (MetaMenu == "1"){
        temp <- list.files(path = "./Metadata", pattern="*.csv", full.names = T)
        numcsv = grep("\\.csv", temp)
        if(sum(numcsv) >= 1){
          metaDat <- lapply(temp, read.csv ) %>% bind_rows(.id = "Csv")
          temp <- SequencesF
          y<- left_join(temp,metaDat, by = "Label")
          write.csv(y, file = "./Outputs/MetadataSeqs.csv")
          rm(temp,y)
          message("MetadataSeqs.csv can now be found in the outputs folder")
        }
      }
      if (MetaMenu == "2"){
        write.dna(Align, file= "./Outputs/Align.fas", format="fasta", nbcol=-1, colsep="")
        message("Align.fas can now be found in the outputs folder")
      }else{
        next
      } 
    } else{write.dna(Align, file= "./Outputs/Align.fas", format="fasta", nbcol=-1, colsep="")
      message("Align.fas can now be found in the outputs folder")
    }
  }
  
  if(finmenu == "B"){
    break
  } 
}



