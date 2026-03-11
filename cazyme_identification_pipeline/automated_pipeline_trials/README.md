# How to use CAZyme Domain annotation. 

## Step 1: Run Hmmer and check coverage
```*USAGE* bash run_01_hmmer.sh <proteome.fasta> ```

## Step 2: Resolve domains, i.e., keep the best hit domain if multiple are overlapping in the same region
```*USAGE*: python 02_run_resolved_domain.py output/dbcan_with_cov.tsv ```

### Step 2.5: Sanity Check 
In between 02 and 03 you can run the following command to check how many are filtered out due to poor coverage 

```
$
for enzyme in GH29 GH95 GH141 GH107 GH168; do ## inset domains of interest
    count1=$(awk -v e="$enzyme" '$1 ~ e {print $4}' prokaryotic_MAGs_vs_dbCAN.domtblout.tsv | sort -u | wc -l)
    count2=$(awk -v e="$enzyme" '$1 ~ e {print $4}' coverage_60_percent_domtblout.tsv | sort -u | wc -l)
    count3=$(awk -v e="$enzyme" '$1 ~ e {print $4}' domain_resolved_coverage_domtblout.tsv | sort -u | wc -l)
    echo "$enzyme  beforeanyfiltering: $count1  aftercoverage: $count2 coverage&overlapcleaned: $count3"
$
```

## Step 3: Run Interpro Scan ##
Interproscan takes time. So its better to run it on your domains of interest rather than running it on all the proteins. This script takes a list of of your enzymes of interest, creates a subset of step 2's output, extracts sequences and runs interpro scan.



*USAGE*: ```bash 03_run_interpro.sh domains_of_interest.txt proteome.fasta resolved_domains.tsv (output step2) ```

## Step 4: Obtain Domain Architecture files 
In this step, the domblout file and the interproscan file will be merged and produce an output. 
It generates a table, in this table a priority to given to CAZydb annotation over pfam, if there is 70% overlap


```*USAGE*: bash 04_run_domain_architecture.sh```


## Step 5:  

Performs mmseqs cluster analysis for each of your files. Geenrates unique cluster for each domain of interest. 
Also filters out those clusters that have <5 proteins per cluster. 

```*USAGE* bash 05_run_cluster_analysis.sh domains.txt proteome.faa merged_domain_architecture_filtered.tsv```

## Step 6: Generate esmfold predicted strucutres for proteins filtered in step 5 

```*USAGE*  bash 06_structure_analysis.sh```





