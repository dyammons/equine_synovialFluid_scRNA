## Instructions to obtain processed Seurat objects  

The processed data are available on Zenodo and can be downloaded by visiting the [project repository web page](https://zenodo.org/doi/10.5281/zenodo.10680754).
Once on the web page scroll down and select download for the file(s) of interest.

Alternatively, use `wget` in a terminal to retrieve the data:
```sh
wget https://zenodo.org/record/11979927/files/eq_synovial_fluid_annotated.rds  # Full synovial fluid dataset
wget https://zenodo.org/record/11979927/files/eq_synovial_fluid_myeloid_annotated.rds  # Myeloid dataset
wget https://zenodo.org/record/11979927/files/eq_synovial_fluid_tcell_annotated.rds  # T cell dataset
wget https://zenodo.org/record/11979927/files/eq_synovium_annotated.rds  # Full synovium dataset
```

Prefer to use tools in python or R? Check out `zenodo_get` or `inborutils` to download within the respective software. 

<details><summary>Example usage of zenodo_get </summary>
<p>

Below is the code needed to install `zendo_get` using `pip` and the command to download the repositiry specific to this project (this should be completed in an environment with python3 installed).  

Visit the [`zendo_get`](https://github.com/dvolgyes/zenodo_get) page for most up to date instructions.

```sh
#install the python tool using pip
pip3 install zenodo_get

#download the Zenodo repository
zenodo_get 10.5281/zenodo.11979927
```

</p>
</details>

## Instructions to obtain count matrices from NCBI GEO  

To download via command line, navigate to directory in which you wish to download the data and run the following:
```sh
#pull down the data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254840/suppl/GSE254840_RAW.tar

#unpack the tar ball
tar -xvf GSE254840_RAW.tar
```

To download via NCBI webpage, navigate to the [GSE254840](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254840) page and download the zip folder containing the count matrices. Once downloaded, you will likely need to extract/unpack the data before using (should be able to do by right clicking on the `GSE254840_RAW.tar` file).

## Instructions to obtain raw data from SRA
Navigate to directory of interest and run for each file you wish to pull down. Then with SRA toolkit installed run:

NOTE: all raw data files are around 1TB of data
```sh
prefetch -v --max-size=55000000 SRR #smallest file for ex
fastq-dump --gzip --split-files SRR
```
From there you will have to modify the file names to be compatible with the Cell Ranger input (if using Cell Ranger). It expects:  
`[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz`  

Feel free to reach out (email or create an issue on GitHub) if you have trouble getting the data downloaded.

