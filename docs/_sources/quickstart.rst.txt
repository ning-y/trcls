Quickstart
==========

SAM Alignment File
------------------

*trcls* requires a pre-aligned SAM file. The program has been used
successfully with *bowtie2* and *HISAT2* aligned files. SAM alignment files
from other sequence alignment programs may not work if they represent alignments
(particularly, the CIGAR and MD fields) differently from bowtie2 or
HISAT2.

The SAM alignment file should be trimmed to the region of interest. This is most
commonly performed with *samtools*. If you intend to use the annotated SAM
alignment file downstream, which is almost always the case, be sure to include
the headers with the ``-h`` option.

For example, to extract the alignments mapping to FLNA to a file flna.sam
from the pre-aligned file pre-aligned.bam.

.. code-block:: bash

    samtools view -h pre-aligned.bam -o flna.sam chrX:154348524-154374638

Obtaining a GTF Annotation File
-------------------------------

trcls additionally requires a GTF annotation file of known transcripts from
the region of interest. Most typically, these annotations originate from a
single gene. There are some requirements that must be met with regard to the
contents of this file.

1. There must be feature lines with a *feature* field of value *'exon'*. These
   are, in fact, the only feature lines which trcls reads.
2. Each feature line with feature field of value 'exon' must have an
   *attribute* field with tag-value pair with tag *'transcript_id'* and value
   corresponding to a unique identifier (usually the accession number)
   corresponding to the transcript with which this exon belongs.

For examples of valid GTF files, please see the file `FLNA.gtf`_ from the
project repository. These files can be generated from the `UCSC Table Browser`_.
For GTF file specifications, please refer to `this page`_.

.. _FLNA.gtf: https://github.com/ningyuansg/trcls/blob/master/test/FLNA.gtf
.. _UCSC Table Browser: https://genome.ucsc.edu/cgi-bin/hgTables
.. _this page: https://ensembl.org/info/website/upload/gff.html

Output
------

The annotated SAM alignment is output via *stdout*. Log messages are sent to
*stderr*.

.. code-block:: bash

    python3 trcls flna.sam flna.gtf > flna-annotated.sam 2> flna.log


Annotations
-----------

Annotations are added to each alignment via an optional field tagged ``TR`` of
type string, i.e. ``TR:Z``. If an alignment fails to annotate, the value for
this field is set as a single asterisk character (``*``); otherwise, it is a
comma separated string of possible transcripts which this sequence read may have
originated from. Possible pre-mRNA is always represented by value *'pre-mRNA'*.

For example:

============================== ====================================
Tag                            Meaning
============================== ====================================
``TR:Z:*``                     No annotation information available.
``TR:Z:pre-mRNA``              Originates from pre-mRNA only.
``TR:Z:NM_001456``             Originates from NM_001456 only.
``TR:Z:NM_001110556,pre-mRNA`` Originates from either pre-mRNA,
                               or NM_0001110556.
============================== ====================================
