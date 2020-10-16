# URA
Metagenomic classification of unmapped reads from targeted NGS sequencing of human DNA

Requirements:
python3, Krona (https://github.com/marbl/Krona), samtools

Example command:
python3 getUnmappedClassifications.py -k /somedirectory/Krona/KronaTools -b /somedirectory/sequences.fasta -f /somedirectory/mysample-T.bam -s mysample-T -o /somedirectory/output/ -v
