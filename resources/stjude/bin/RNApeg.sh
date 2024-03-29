#!/usr/bin/env sh
# - if building from Alpine, must use "sh" since Alpine doesn't support bash
# - if building from Ubuntu, must use "bash" as "sh" doesn't support all syntax

BAMFILE=""
# RNA BAM
FASTA=""
# reference genome
REFFLAT=""
# main junction file used in correction
REFGENE=""
# optional refFlat for final gene annotation step only, if not specified REFFLAT
# will be used.
OUTPUT_DIR=/results

usage() {
    echo "RNApeg.sh [-h] -b bamfile -f fasta -r refflat [-rg refflat]"
}

while [ ! -z "$1" ]; do
    case "$1" in
        -h) usage && exit 0; shift;;
        -b) BAMFILE=$2; shift;;
        -r) REFFLAT=$2; shift;;
        -f) FASTA=$2; shift;;
        -rg) REFGENE=$2; shift;;
	*) echo "ERROR: unrecognized parameter $1"; usage; exit 1; shift;;
    esac
    shift
done


#######################
### Validate inputs ###
#######################
if [[ ! -f "$BAMFILE" ]]; then
    >&2 echo "ERROR: Bamfile '$BAMFILE' DNE"
    usage
    exit 1
elif [[ ! -f "$BAMFILE.bai" ]]; then
    >&2 echo "ERROR: Bamfile index (.bai) '$BAMFILE.bai' DNE"
    usage
    exit 1
elif [[ ! -f $REFFLAT ]]; then
    >&2 echo "ERROR: refFlat file '$REFFLAT' DNE"
    usage
    exit 1
elif [[ ! -f $FASTA ]]; then
    >&2 echo "ERROR: FASTA file '$FASTA' DNE"
    usage
    exit 1
elif [[ ! -f $FASTA.fai ]]; then
    >&2 echo "ERROR: FASTA index file '$FASTA.fai' DNE; see samtools faidx command"
    usage
    exit 1
elif [[ ! -d "$OUTPUT_DIR" ]]; then
    usage
    >&2 echo "ERROR: Output directory '$OUTPUT_DIR' does not exist; need mountpoint?"
    exit 1
fi

if [[ -z "$REFGENE" ]]; then
    REFGENE=$REFFLAT
    echo set REFGENE to default $REFGENE
else
    echo user REFGENE is $REFGENE
fi

if [[ ! -f $REFGENE ]]; then
    >&2 echo "ERROR: refGene file '$REFGENE' DNE"
    usage
    exit 1
fi

##################
### Run RNApeg ###
##################

echo "[*] Running junction_extraction_wrapper.pl"
junction_extraction_wrapper.pl -no-config -bam $BAMFILE -o $OUTPUT_DIR -fasta $FASTA -refflat $REFFLAT -refgene $REFGENE -now


