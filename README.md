Author: vaibhakv


# VCF Explorer Shiny App

![VCF Explorer Logo](https://github.com/vaibhakv/VCFexplorer/blob/main/images/Screenshot%20from%202023-12-15%2015-26-55.png)

## Overview
Interface for exploring and visualizing VCF files. Users can upload VCF files, choose a chromosome, and visualize the distribution of Single Nucleotide Variants (SNVs) and Insertion/Deletion variants (Indels) on the selected chromosome.

## Features

- **File Upload:** Choose a VCF file from your device.
- **Chromosome Selection:** Select the chromosome you want to explore.
- **Variant Counts:** Visualize the counts of SNVs and Indels in a bar plot.
- **Chromosome Visualization:** View the genomic coordinates of variants on an Ideogram.
- **Sample Information:** Display sample names and total counts for SNVs and Indels.

## How to Use

1. **Upload VCF File:** Click on "Choose a VCF File" and select your VCF file.
2. **Select Chromosome:** Choose a chromosome from the dropdown menu.
3. **Submit:** Click on the "Submit" button to generate visualizations.
4. **Explore:** Explore the bar plot and genomic coordinates on the Ideogram.

## Installation and Setup

### Prerequisites

- R
- Shiny library
- VariantAnnotation library
- Gviz library
- GenomicRanges library

### Run the App

```R
# Install required libraries (if not installed)
install.packages(c("shiny", "VariantAnnotation", "Gviz", "GenomicRanges"))

# Run the app
shiny::runGitHub("your-vcf-explorer-repo")
