library(devtools)
if(!require("BiocManager",quietly = T)){
  install.packages("BiocManager")
}

BiocManager::install("GenomicRanges")
devtools::install_github("VanLoo-lab/ascat/ASCAT")
BiocManager::install("VariantAnnotation")
BiocManager::install("copynumber")
devtools::install_github("Wedge-Oxford/battenberg")


library(Battenberg)


setwd("/media/kovac/Resources1/MartinD")

# na ziskanie tumour a normal name z bam filu
# samtools view -H filename | grep "@RG*" | awk '{print $10}' | awk -F:  ' { print $2}' | head -n 1

install_version("readr", "1.3.1")

devtools::install_github("Wedge-Oxford/battenberg@HEAD")



impute_basedir = "/media/kovac/Rigel/Osteosarcoma.ICGC.NatComm.2017/"

#

#find -type f -name *.T.bam | samtools view -H *.T.bam | grep "@RG*" | awk '{print $10}' | awk -F:  ' { print $2}' | head -n 1

  # ziskanie vsetkych filov
files <- dir("/media/kovac/Rigel/Osteosarcoma.ICGC.NatComm.2017", recursive=TRUE, full.names=TRUE, pattern="\\.bam$")
files <- sub(impute_basedir,"",files)
files <- strsplit(files," ")
library(stringr)
library(stringi)
files <- str_sort(files, numeric = TRUE)
files
#tito pacienti uz zbehli
files <- files[-(1:35)]
#files <- files[-23]
files

#c_files <- read.delim("/media/kovac/Resources1/MartinD/cfiles.txt", header = FALSE)
#t_files <- read.delim("/media/kovac/Resources1/MartinD/tufiles.txt", header = FALSE)
#View(c_files)
#rg_sign <-paste("","@RG*","",sep="'")

#for (i in c_files) {
#  command_to_exec <- ""
#  command_to_exec <- paste("samtools view -H ",i)
#  grep_ap <- paste("| grep ",rg_sign," ")
#  command_to_exec <- paste(command_to_exec,grep_ap,"| awk '{print $10}' | awk -F:  ' { print $2}' | head -n 1")
#  print(command_to_exec)
#}

#for (i in t_files) {
#  command_to_exec <- ""
#  command_to_exec <- paste("samtools view -H ",i)
#  grep_ap <- paste("| grep ",rg_sign," ")
#  command_to_exec <- paste(command_to_exec,grep_ap,"| awk '{print $10}' | awk -F:  ' { print $2}' | head -n 1")
#  print(command_to_exec)
#}

normal_names <- read.delim("/media/kovac/Resources1/MartinD/normal_names.txt", header = FALSE)
tumour_names <- read.delim("/media/kovac/Resources1/MartinD/tumour_names.txt", header = FALSE)
View(normal_names)
View(tumour_names)
move_index <- 1
prev_index <- 0

for(i in 19:26){
  setwd("/media/kovac/Resources1/MartinD/")
  dir_to_create_name <- paste("P",i,"_batt_done",sep="")
  tryCatch({
    system(paste("mkdir ",dir_to_create_name))
    setwd(paste("/media/kovac/Resources1/MartinD/",dir_to_create_name,sep=""))
  })
  first_data <- paste(impute_basedir,files[move_index+prev_index],sep = "")
  second_data <- paste(impute_basedir,files[move_index+prev_index+1],sep = "")
  prev_index <- prev_index + 1
  normal_data <- ""
  tumour_data <- ""

  if(stri_sub(first_data,-5,-5) == "T"){
    normal_data <- second_data
    tumour_data <- first_data
  }else{
    normal_data <- first_data
    tumour_data <- second_data
  }




  tryCatch({
    battenberg(
      normalname = normal_names$V1[move_index],
      tumourname = tumour_names$V1[move_index],

      normal_data_file = normal_data,
      tumour_data_file = tumour_data,

      g1000prefix = "/home/kovac/Downloads/Battenberg/battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr",
      g1000allelesprefix ="/home/kovac/Downloads/Battenberg/battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr" ,
      gccorrectprefix = "/home/kovac/Downloads/Battenberg/battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_",
      repliccorrectprefix = "/home/kovac/Downloads/Battenberg/battenberg_wgs_replic_correction_1000g_v3/1000_genomes_replication_timing_chr_",
      imputeinfofile = "/home/kovac/Downloads/Battenberg/impute_info.txt",
      problemloci = "/home/kovac/Downloads/Battenberg/probloci_270415.txt.gz",
      skip_allele_counting = F,
      skip_preprocessing = F,
      skip_phasing = F,
      allelecounter_exe = "/media/kovac/Resources1/MartinD/dockerbattenberg/alleleCount/bin/alleleCounter",
      #impute_exe = "/home/kovac/miniconda3/bin/impute2",
      impute_exe = "/home/kovac/Downloads/impute_v2.3.2_x86_64_dynamic/impute2",
      ismale = T,
      write_battenberg_phasing = F,
      GENOMEBUILD = "hg19",
      data_type = "wgs",
      nthreads = 10,
      chrom_coord_file = "/home/kovac/Downloads/Battenberg/gcCorrect_chromosome_coordinates_hg19.txt"
    )},error = function(e) {
      next
    }
  )
  move_index <- move_index + 1
}

# sed -i 's/ /\t/g' f85ae2b7-cebf-17a2-e040-11ac0c48033a*_allHaplotypeInfo.txt