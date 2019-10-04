#!/usr/bin/env Rscript

#source("/home/bioinformatica/datos/03_Analysis/jgonzalez.69/scripts_chipseq/library_functions_ZCL_chip.R")
source("./functions_PeakZCL.R")

#change this path for your own

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="PATH to a chipseq file name", metavar="character"),
  make_option(c("-c", "--control"), type="character", default=NULL,
              help="PATH to a control file name", metavar="character"),
  make_option(c("-d", "--output_dir"), type="character", default=NULL,
              help="A directory for ZCL to create its output", metavar="character"),
  make_option(c("-n", "--normalization"), type="character",
              help="Select the signal normalization method.Options are SES,SDS an N (without).Default method is SES", metavar="character"),
  make_option(c("-w", "--Decimating"), type="numeric", default=50,
              help="Decimating signal value.If you decide not to decimate the original signal, the program will take a long time. If the value is high, you can lose a lot of information. Default is 50", metavar="character"),
  make_option(c("-l", "--length_signal_wavelet"), type="numeric", default=200000,
              help="Length signal wavelet.Default value is 2000000.If signal is longer than 4000000, wavelet function will not work.If the value is less than 100000, the peak calling function may cause errors", metavar="character"),
  make_option(c("-s", "--scale"), type="numeric", default=30,
              help="Number of scales of wavelet. Default is 30", metavar="character"),
  make_option(c("-g", "--Levth"), type="numeric", default=5,
              help="Zero-crossing lines leverage threshold. Default is 5", metavar="character"),
  make_option(c("-i", "--Threads"), type="numeric", default=0,
              help="Baseline intensity threshold. Default is 0", metavar="character"),
  make_option(c("-j", "--clustering"), type="numeric", default=500,
              help="Bp range between peaks to clustering purposes. Default is 500", metavar="character"),
  make_option(c("-a", "--area_chip"), type="numeric", default=1.6,
              help="Value of the area of the selected peaks. The value is given in log10(x). Default is 0.69", metavar="character"),
  make_option(c("-k", "--fold_change"), type="numeric", default=0.1,
              help="Difference between peak area on the chip and control. The value is given in log10(x). Default is 0.2", metavar="character"),
  make_option(c("-o", "--out_file"), type="character",default=NULL,
              help="Output BED file name.Default name will be the same at your experiment.", metavar="character"),
  make_option(c("-y", "--all_peaks"), type="character", default=FALSE,
              help="Show all peaks detected", metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if(is.null(opt$file)){
  print_help(opt_parser)
  stop("At least two argument must be supplied (input file and output directory)","\n",
  "Example of use: Rscript (PATH to ZCL_Rscript file in your system) -f (PATH to chipseq file) -c (PATH to control file), -d (PATH to output directory)", call.=FALSE)
}

fileName_chip<-basename(opt$file)
cat("loading chip file ",fileName_chip,"\n")
chipchrom<-fread(opt$file)
chipchrom<-as.data.frame(chipchrom)
chr_number<-as.character(unlist(chipchrom$V1[1]))

# cat("decimating chip signal...",fileName_chip,"\n")
decimated_signal_chip<-decimate(chipchrom$V3,opt$Decimating)
decimated_signal_chip[decimated_signal_chip<0]<-0

# cat("estimating CWT and peak calling...",fileName_chip,"\n")
start_end_peaks<-estimating_CWT_and_peak_calling(decimated_signal_chip,opt$length_signal_wavelet,opt$scale,opt$Levth,opt$Decimating)

# cat("clustering.....",fileName_chip,"\n")
peaks_clustered_no_filtered<-bining(start_end_peaks,opt$clustering/opt$Decimating)
peaks_clustered_no_filtered<-unique(peaks_clustered_no_filtered)

peaks_clustered_Ranges<-IRanges(start=peaks_clustered_no_filtered$start,end=peaks_clustered_no_filtered$end)
filtered_peaks<-peaks_clustered_Ranges[width(peaks_clustered_Ranges)<=10000/opt$Decimating] #ultimas listas con 1
df_filtered_peaks<-as.data.frame(filtered_peaks)
peaks_clustered<-data.frame("start"=df_filtered_peaks$start,"end"=df_filtered_peaks$end)

# cat("peaks in original signal....",fileName_chip,"\n")
peaks_in_position<-((peaks_clustered)-1)*opt$Decimating

####ahora debemos generar una "señal" que es subsetting de la diezmada, que son los rangos de los posiciones
# cat("selected signal from chip","\n")
signal_selected_chip<-c()
for (i in 1:nrow(peaks_in_position)){
  peaks_signal_ranges_chip<-chipchrom$V3[peaks_in_position$start[i]:peaks_in_position$end[i]]
  signal_selected_chip<-c(signal_selected_chip,peaks_signal_ranges_chip)
}

if(!is.null(opt$control)){
  fileName_control <-basename(opt$control)
  cat("loading input file ",fileName_control,"\n")
  controlchrom<-fread(opt$control)
  # cat("selected signal from input","\n")

  signal_selected_input<-c()
  for (i in 1:nrow(peaks_in_position)){
    peaks_signal_ranges_input<-controlchrom$V3[peaks_in_position$start[i]:peaks_in_position$end[i]]
    signal_selected_input<-c(signal_selected_input,peaks_signal_ranges_input)
  }


  if(opt$normalization=="SES"){
    # cat("normalizing signal....",fileName_control,"with ", opt$normalization,"method","\n")
    factor<-NormalizationChIP(signal_selected_chip,signal_selected_input,method="SES")
  }else if(opt$normalization=="SDS"){
    # cat("normalizing signal....",fileName_control,"with ", opt$normalization,"method","\n")
    factor<-NormalizationChIP(signal_selected_chip,signal_selected_input,method="SDS")
  }else if(opt$normalization=="N"){
    factor<-1
  }
  c_normalized<-(controlchrom$V3)*(factor)
  # cat("decimating input signal...",fileName_control,"\n")
  decimated_control_signal<-decimate(c_normalized,opt$Decimating)
  decimated_control_signal[decimated_control_signal<0]<-0
}



# cat("quantification...",fileName_chip,"\n")
quantified<-quantification(peaks_clustered,decimated_signal_chip)
if(!is.null(opt$control)){
  quantified_input<-quantification(peaks_clustered,decimated_control_signal)
}



if(!is.null(opt$control)){
  # cat("statistics...",fileName_chip,"\n")
  peaks_in_position<-statistics(peaks_in_position,opt$area_chip,opt$fold_change,quantified,quantified_input)
}else{
  peaks_in_position$area_chip<-log10(unlist(quantified))
}

cat("printing peaks file in BED format for ",fileName_chip,"\n")
options(scipen=100)
bed<-BED_files(peaks_in_position,chr_number,opt$control,opt$all_peaks)


if(is.null(opt$out_file)){
  fileName <- paste0(opt$output_dir,fileName_chip,".bed")
  write.table(bed,file=fileName,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

}else{
  write.table(bed,file=paste0(opt$output_dir,opt$out_file,".bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

}

cat(fileName_chip,"DONE!","\n")
