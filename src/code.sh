#!/bin/bash
set -e -o pipefail

main() {

  #################
  # Download data #
  #################

  echo ""
  echo "=== Setup ==="
  echo "  [*] Downloading input files..." 
  dx-download-all-inputs --parallel > /dev/null

  ################
  # Housekeeping #
  ################

  echo "  [*] Performing some housekeeping..."

  # Path setup 

  export PATH=$PATH:/stjude/bin
  export PERL5LIB=/stjude/lib/perl:/stjude/bin
  #export CLASSPATH=/stjude/lib/java/
  for jar in $(ls /stjude/lib/java/*.jar); do
    export CLASSPATH=$CLASSPATH:$jar
  done
  # Link binary to where our pipeline looks for it
  ln -s /usr/bin/env /bin/env 

  echo "   [-] Classpath: $CLASSPATH"
  echo "   [-] Perl5lib:  $PERL5LIB"

  # Cloud config 
  export SJ_CONFIGS=/stjude/configs
  . import_config.sh genome "GRCh37-lite"

  # Linking files where cloud configs expects them
  ln -s $bam_path /home/dnanexus/${bam_name}
  ln -s $bam_index_path /home/dnanexus/${bam_index_name}

  ##########
  # RNApeg #
  ##########

  echo "=== RNApeg ==="

  echo " [*] SJ_CONFIGS=$SJ_CONFIGS"
  echo " [*] Running junction_extraction_wrapper.pl"
  if [ ${ref_name} == "GRCh37-lite" ]
  then
    junction_extraction_wrapper.pl -no-config -bam /home/dnanexus/${bam_name} -o /home/dnanexus -fasta /stjude/reference/Homo_sapiens/${ref_name}/FASTA/${ref_name}.fa -refflat /stjude/reference/Homo_sapiens/${ref_name}/mRNA/Combined/all_refFlats.txt -refgene /stjude/reference/Homo_sapiens/${ref_name}/mRNA/Combined/all_refFlats.txt -now
  else
    junction_extraction_wrapper.pl -no-config -bam /home/dnanexus/${bam_name} -o /home/dnanexus -fasta /stjude/reference/Homo_sapiens/${ref_name}/FASTA/${ref_name}.fa -refflat /stjude/reference/Homo_sapiens/${ref_name}/mRNA/RefSeq/refFlat-sharp.txt -refgene /stjude/reference/Homo_sapiens/${ref_name}/mRNA/RefSeq/refFlat-sharp.txt -now
  fi
  #echo " [*] Running RNApeg via Docker"
  #docker run -v /home/dnanexus:/data -v /stjude/reference:/reference mnedmonson/public:rnapeg RNApeg.sh -b /data/${bam_name} -r /reference/Homo_sapiens/GRCh37-lite/mRNA/Combined/all_refFlats.txt -f /reference/Homo_sapiens/GRCh37-lite/FASTA/GRCh37-lite.fa

  #mv /home/dnanexus/sample.bam.junctions.tab.shifted.tab /home/dnanexus/${bam_prefix}.junctions.tab.shifted.tab

  junctions=$(dx upload ${bam_name}.junctions.tab --brief)
  dx-jobutil-add-output junctions "$junctions" --class=file
  junctions_shifted=$(dx upload ${bam_name}.junctions.tab.shifted.tab --brief)
  dx-jobutil-add-output junctions_shifted "$junctions_shifted" --class=file
  junctions_shifted_bed=$(dx upload ${bam_name}.junctions.tab.shifted.bed --brief)
  dx-jobutil-add-output junctions_shifted_bed "$junctions_shifted_bed" --class=file
  junctions_annotated=$(dx upload ${bam_name}.junctions.tab.shifted.tab.annotated.tab --brief)
  dx-jobutil-add-output junctions_annotated "$junctions_annotated" --class=file

}
