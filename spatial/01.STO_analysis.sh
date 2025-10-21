###download gtf and fasta from ensembl
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

###install singularit and download SAW from DockerHub
singularity build SAW_6.1.sif docker://stomics/saw:06.1.0

###build index
singularity exec SAW_6.1.sif mapping --runMode genomeGenerate \
    --genomeDir reference/STAR_SJ100 \
    --genomeFastaFiles reference/genome.fa \
    --sjdbGTFfile reference/genes.gtf \
    --sjdbOverhang 99 \
    --runThreadN 12


####mapping

bash ./SAW/stereoPipeline_v6.1.sh \
       -sif SAW_6.1.sif \
       -genomeSize 4 \
       -splitCount 1 \
       -maskFile SN.h5 \
       -fq1 ${path}/${lane1}/${lane1}_R1.fq.gz,...,${path}/${laneN}/${laneN}_R1.fq.gz \
       -fq2 ${path}/${lane1}/${lane1}_R1.fq.gz,...,${path}/${laneN}/${laneN}_R1.fq.gz \
       -speciesName <speciesName> \ ###human mouse
       -tissueType <tissueName> \   ####embryo
       -refIndex reference/STAR_SJ100 \
       -annotationFile reference/genes/Homo_sapiens.GRCh38.110.gtf \
       -threads 5 \
       -outDir ${outDir}/result \
       -imageRecordFile ${path}/image/<SN_date_time_version>.ipr \ # [optional] when image is given and has passed QC
       -imageCompressedFile ${path}/image/<SN_date_time_version>tar.gz \ # [optional] when image is given and has passed QC
       -doCellBin Y  # [optional] when you want to do the cell segmentation and get cell gene expression data
