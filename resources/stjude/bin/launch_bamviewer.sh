#!/bin/sh
# Need to add:
# -bam BAMFILE
# -name CHR
# EITHER:
#   -start START -end END
#   OR
#   -center CENTER
export CLASSPATH=/nfs_exports/apps/gnu-apps/NextGen/NewView/av.jar:/nfs_exports/apps/gnu-apps/NextGen/picard-tools-1.13/sam-1.13.jar:/nfs_exports/apps/gnu-apps/NextGen/edmonson_java/bin/mysql-connector-java-5.1.10-bin.jar
#/nfs_exports/apps/compilers/jdk1.6.0_15/bin/java -Xmx3072m Ace2.AceViewer -nib  /nfs_exports/apps/gnu-apps/NextGen/nextgensupport/WashU_hg18_nib -db-server SJMEMGB01.stjude.org -db-database hg18 -db-user pallas $*
/nfs_exports/apps/compilers/jdk1.6.0_15/bin/java -Xmx3072m Ace2.AceViewer -nib  /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/NIB -db-server SJMEMGB01.stjude.org -db-database GRCh37-lite -db-user pallas $*

