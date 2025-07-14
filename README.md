# zhangzhen/bulkfoodiepipeline

## Introduction

**zhangzhen/bulkfoodiepipeline** is a bioinformatics pipeline used for in-vivo bulk FOODIE sequencing data. It pre-processes raw data from FastQ inputs, aligns the reads, computes conversion rates of cytosines, and calls footprints of transcription factors.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
sample,fastq_1,fastq_2
b250417-MB-Tn5nwz-B2462-06-2h-CS3,b250417-MB-Tn5nwz-B2462-06-2h-CS3_L3_Q0060W0149.R1.fastq.gz,b250417-MB-Tn5nwz-B2462-06-2h-CS3_L3_Q0060W0149.R2.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run zhangzhen/bulkfoodiepipeline \
   -profile conda \
   --input samplesheet.csv \
   --depth <NUM> \
   --expected_ratio_file <EXPECTED_RATIO_FILE> \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

zhangzhen/bulkfoodiepipeline was originally written by Zhang Zhen, Shen Ke, Wang Quangui.

We thank the following people for their extensive assistance in the development of this pipeline:
- group members at Sunney Xie Lab

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use `bulk-foodie-pipeline` for your analysis, please cite the `FOODIE` article as follows:

> R. He, W. Dong, Z. Wang, C. Xie, L. Gao, W. Ma, K. Shen, D. Li, Y. Pang, F. Jian, J. Zhang, Y. Yuan, X. Wang, Z. Zhang, Y. Zheng, S. Liu, C. Luo, X. Chai, J. Ren, Z. Zhu, & X.S. Xie, **Genome-wide single-cell and single-molecule footprinting of transcription factors with deaminase**, _Proc. Natl. Acad. Sci. U.S.A._ 121 (52) e2423270121, [doi: 10.1073/pnas.2423270121](https://doi.org/10.1073/pnas.2423270121) (2024).

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use zhangzhen/bulkfoodiepipeline for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
