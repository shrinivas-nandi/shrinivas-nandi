cd output 

mkdir -p tmp

# move all extra files 
for f in *; do
    case "$f" in
        domains_subset.tsv|interpro_annotations.tsv|proteins_of_interest.fasta|tmp)
            continue
            ;;
        *)
            mv "$f" tmp/
            ;;
    esac
done


# Now with the saved files generate an architecture file.
# format: protein_name length component database start end description

awk -F'\t' 'BEGIN{OFS="\t"}
{
    component=$1
    sub(/\.hmm$/, "", component)
    print $4,$6,component,"dbCAN",$18,$19,"NA"
}' domains_subset.tsv > dbcan_for_domain_architecture.tsv


# Pfam extraction (safe for description text)
awk -F'\t' 'BEGIN{OFS="\t"}
$4=="Pfam" {
    print $1,$3,$5,$4,$7,$8,$6
}' interpro_annotations.tsv > pfam_for_domain_architecture.tsv


# Merge and order by protein then start
cat dbcan_for_domain_architecture.tsv pfam_for_domain_architecture.tsv \
| sort -t $'\t' -k1,1 -k5,5n \
> merged_domain_architecture.tsv


# Identify proteins that have Pfam but no dbCAN hit
awk -F'\t' '$4=="Pfam"{print $1}' merged_domain_architecture.tsv | sort -u > tmp/pfam_proteins.txt
awk -F'\t' '$4=="dbCAN"{print $1}' merged_domain_architecture.tsv | sort -u > tmp/dbcan_proteins.txt

awk -F'\t' 'BEGIN{OFS="\t"}
{
    prot=$1
    db=$4
    start=$5
    end=$6

    if(db=="dbCAN"){
        db_start[prot,++n[prot]]=start
        db_end[prot,n[prot]]=end
        print
        next
    }

    # check overlap with dbCAN
    keep=1
    for(i=1;i<=n[prot];i++){
        s=db_start[prot,i]
        e=db_end[prot,i]

        overlap_start=(start>s?start:s)
        overlap_end=(end<e?end:e)

        if(overlap_start<=overlap_end){
            overlap=overlap_end-overlap_start+1
            shorter=((end-start)<(e-s)?(end-start):(e-s))+1

            if(overlap/shorter >= 0.7){
                keep=0
                break
            }
        }
    }

    if(keep) print
}' merged_domain_architecture.tsv > merged_domain_architecture_filtered.tsv


comm -23 tmp/pfam_proteins.txt tmp/dbcan_proteins.txt > pfam_only_proteins.txt

rm dbcan_for_domain_architecture.tsv
rm merged_domain_architecture.tsv
rm pfam_for_domain_architecture.tsv









