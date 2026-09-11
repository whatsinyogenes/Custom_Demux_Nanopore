# Custom Nanopore Demultiplexing

NOTE: This was for personal use and I thought it may be useful to others. 

A simple bash pipeline for demultiplexing Illumina derived unique dual index "IDT for Illumina". This pipeline can be customised for any index. There are numerous applications such as: sequencing QC before sending for illumina, pooled sequencing without the use of NBD114, plasmid sequencing, etc.

<p align="center">
  <img src="/Custom_demux.jpg" width="1000">
</p>

chatgpt^^

Usage:

  bash Nano_process_pipeline \
      --input FASTQ_OR_DIRECTORY \
      --indexes INDEX_CSV \
      --samples SEQUENCING_MANIFEST_TSV \
      --reference REFERENCE_FASTA \
      --adapters FASTPLONG_ADAPTER_FASTA \
      --output OUTPUT_DIRECTORY


Required arguments:

  --input PATH
        FASTQ file or directory containing FASTQ/FASTQ.GZ files.

  --indexes CSV
        Master index manifest CSV containing all UDP barcode sequences.
        
    Index,i7 in Adapter,i7 for Sample Sheet,i5 in Adapter,i5 for Sample Sheet in Forward Orientation,I5 Bases for Sample Sheet in Reverse Complement Orientation
        UDP0001,CGCTCAGTTC,GAACTGAGCG,TCGTGGAGCG,TCGTGGAGCG,CGCTCCACGA
        UDP0002,TATCTGACCT,AGGTCAGATA,CTACAAGATA,CTACAAGATA,TATCTTGTAG
        UDP0003,ATATGAGACG,CGTCTCATAT,TATAGTAGCT,TATAGTAGCT,AGCTACTATA

  --samples TSV
        Sequencing manifest specifying which UDPs were actually used.

        Format:

            sample_name    index_pair
            patient1       UDP0001
            patient2       UDP0003
            patient3       UDP0007

  --reference FASTA
        Reference genome FASTA for minimap2.

  --adapters FASTA
        Adapter FASTA for Fastplong to guide trimming.
      
    Format:

    >i7_front
    CAAGCAGAAGACGGCATACGAGAT

    >i7_rear
    GTCTCGTGGGCTCG

    >i5_front
    AATGATACGGCGACCACCGAGATCTACAC

    >i5_rear
    GTCGGCAGCGTC



  --output DIRECTORY
        Output directory.

General options:

  --threads N
        Number of threads.
        Default: 4

  --dorado PATH
        Dorado executable.
        Default: dorado

  --fastplong PATH
        Fastplong executable.
        Default: fastplong

  --minimap2 PATH
        minimap2 executable.
        Default: minimap2

  --samtools PATH
        samtools executable.
        Default: samtools

  --python PATH
        Python executable.
        Default: python3


Fastplong:

  --fastplong-passes N
        Number of Fastplong passes.
        Default: 3

  --distance-threshold VALUE
        Fastplong distance threshold.
        Default: 0.10

  --trimming-extension N
        Fastplong trimming extension.
        Default: 0

  --min-length N
        Minimum retained read length.
        Default: 50


Alignment:

  --minimap2-preset PRESET
        minimap2 preset.
        Default: sr


Dorado scoring:

  --max-barcode-penalty N
  --barcode-end-proximity N
  --min-barcode-penalty-dist N
  --min-separation-only-dist N
  --flank-left-pad N
  --flank-right-pad N
  --front-barcode-window N
  --rear-barcode-window N
  --midstrand-flank-score VALUE
  --rear-only-barcodes true|false


Barcode flanks:

  --i7-front SEQUENCE
  --i7-rear SEQUENCE
  --i5-front SEQUENCE
  --i5-rear SEQUENCE


Pipeline control:

  --resume
        Resume completed stages.
        Default: enabled

  --no-resume
        Do not use .done markers.

  --force
        Re-run stages even if .done markers exist.
        Existing stage output is cleaned before rerunning.

  --keep-intermediate
        Keep intermediate files.
        Default: enabled

  --no-keep-intermediate
        Remove intermediate files when possible.

  --skip-demux
        Skip Dorado demultiplexing.

  --skip-trimming
        Skip Fastplong.

  --skip-alignment
        Skip minimap2/samtools.


Help:

  -h
  --help
