usage="bash sraDownloader.sh ACCESSION OUTDIR"

# standard input
accession=$1
outdir=$2

# change accordingly to your system (require aspera connect to be installed)
asperaExec=$HOME/.aspera/connect/bin/ascp
asperaKey=$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh
sraOutDir=$HOME/ncbi/public/sra

# remove lock, if present
lockName=$sraOutDir/$accession.sra.lock
if [[ -e $lockName ]]; then rm -rf $lockName; fi

# download .sra file for $accession
clPrefetch="prefetch -t ascp -a \"$asperaExec|$asperaKey\" $accession"
echo "prefetch cl used: $clPrefetch"
eval $clPrefetch
sraOutFile=$sraOutDir/$accession.sra

# extract splitted .fastq reads from .sra
if [[ -e $outdir ]]; then rm -rf $outdir; fi
fastq-dump $sraOutFile --split-files -O $outdir --gzip
rm $sraOutFile
