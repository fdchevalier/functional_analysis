#!/bin/bash
# Title: func-table-mut.sh
# Version: 0.5
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2017-09-11
# Modified in: 2018-01-18
# Licence : GPL v3



#======#
# Aims #
#======#

aim="Generate mutation table from a VCF file for a given gene. The table contains gene and CDS positions of the mutation, protein mutation and allele and sample numbers."



#==========#
# Versions #
#==========#

# v0.5 - 2018-01-18: complex event (FreeBayes specific) column added
# v0.4 - 2018-01-16: warning when indels toward the end cannot be process / bug regarding diff output corrected
# v0.3 - 2017-10-10: output naming option added / genomic position added
# v0.2 - 2017-10-02: bug due to * alt allele (missing data due to upstream deletion) corrected
# v0.1 - 2017-09-28: bug induced by Haplotype Caller (SNP coded as InDel when InDel at the same site too) corrected / reference sequenced added in fasta
# v0.0 - 2017-09-01: creation

version=$(grep -i -m 1 "version" "$0" | cut -d ":" -f 2 | sed "s/^ *//g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -v|--vcf file -r|--ref file -g|-gff file -o|-out file -p|--pop file -h|--help

Aim: $aim

Version: $version

Options:
    -v, --vcf       VCF input file
    -r, --ref       reference genome
    -g, --gff       GFF file of the gene of interest
    -o, --out       name of the output report [default: gff-filename.tsv]
    -p, --pop       population file (tab separated values) [optional]. This file must be formated as follows:
                        - first column contains sample names
                        - second column contains corresponding population names
    -h, --help      this message
    "
}


# Info message
function info {
    if [[ -t 1 ]]
    then
        echo -e "\e[32mInfo:\e[00m $1"
    else
        echo -e "Info: $1"
    fi
}


# Warning message
function warning {
    if [[ -t 1 ]]
    then
        echo -e "\e[33mWarning:\e[00m $1"
    else
        echo -e "Warning: $1"
    fi
}


# Error message
## usage: error "message" exit_code
## exit code optional (no exit allowing downstream steps)
function error {
    if [[ -t 1 ]]
    then
        echo -e "\e[31mError:\e[00m $1"
    else
        echo -e "Error: $1"
    fi

    if [[ -n $2 ]]
    then
        exit $2
    fi
}


# Dependency test
function test_dep {
    which $1 &> /dev/null
    if [[ $? != 0 ]]
    then
        error "Package $1 is needed. Exiting..." 1
    fi
}


# Progress bar
## Usage: ProgressBar $mystep $myend
function ProgressBar {
    if [[ -t 1 ]]
    then
        # Process data
        let _progress=(${1}*100/${2}*100)/100
        let _done=(${_progress}*4)/10
        let _left=40-$_done
        # Build progressbar string lengths
        _fill=$(printf "%${_done}s")
        _empty=$(printf "%${_left}s")

        # Build progressbar strings and print the ProgressBar line
        # Output example:
        # Progress : [########################################] 100%
        #printf "\rProgress : [${_fill// /=}${_empty// / }] ${_progress}%%"
        printf "\r\e[32mProgress:\e[00m [${_fill// /=}${_empty// / }] ${_progress}%%"
        
        [[ ${_progress} == 100 ]] && echo ""
    fi
}


# Clean up function for trap command
## Usage: clean_up file1 file2 ...
function clean_up {
    rm -rf $@
    exit 1
}


# Reverse complement
## Usage: rev_cpt fasta_file
function rev_cpt {
    revseq -notag -sequence "$1" -filter | grep -v ">" | tr -d "\n" | sed "s/$/\n/" > "${1}_rev"
    mv "${1}_rev" "$1"
}


#==============#
# Dependencies #
#==============#

test_dep sed
test_dep bedtools
test_dep revseq     # EMBOSS package
test_dep transeq    # EMBOSS package



#===========#
# Variables #
#===========#

set -e

# Options
while [[ $# -gt 0 ]]
do
    case $1 in
        -v|--vcf    ) myvcf="$2"             ; shift 2 ;;
        -r|--ref    ) mygenome="$2"          ; shift 2 ;;
        -g|--gff    ) mygff="$2"             ; shift 2 ;;
        -o|--out    ) report="${2%.tsv}.tsv" ; shift 2 ;;
        -p|--pop    ) mypop_f="$2"           ; shift 2 ;;
        -h|--help   ) usage ; exit 0 ;;
        *           ) error "Invalid option: $1\n$(usage)" 1 ;;
    esac
done


# Check the existence of obligatory options
if [[ -z "$myvcf" ]]
then
    error "A VCF is required. Exiting...\n$(usage)" 1
elif [[ -z "$mygenome" ]]
then
    error "A reference genome is required. Exiting...\n$(usage)" 1
elif [[ -z "$mygff" ]]
then
    error "A GFF is required. Exiting...\n$(usage)" 1
fi

# Check nature of files
[[ ! -s "$myvcf" ]]    && error "VCF file does not exist or is empty. Exiting..." 1
[[ ! -s "$mygenome" ]] && error "Genome file does not exist or is empty. Exiting..." 1
[[ ! -s "$mygff" ]]    && error "GFF file does not exist or is empty. Exiting..." 1
[[ -n "$mypop_f" && ! -s "$mypop_f" ]] && error "Population file does not exist or is empty. Exiting..." 1


# GFF filename for output
fn=$(basename "${mygff%.*}")

# Output filename from GFF
if [[ -z "$report" ]]
then
    fn_out="$fn"
    report="${fn}.tsv"
else 
    fn_out=$(basename "${report%.*}")
fi
[[ -s "$report" ]] && error "Output file $report exist. Exiting..." 1

# Sequence output files
seq_cds_f="[na]_${fn_out}_cds.fa"
seq_aa_f="[aa]_${fn_out}.fa"
[[ -s "$seq_cds_f" || -s "$seq_aa_f" ]] && error "CDS or amino acid sequence files exist. Exiting..." 1

# Temporary files and folder
tmp_fd="tmp"
tmp="$tmp_fd/mygene"
report_tmp="$tmp_fd/report"

[[ -d "$tmp_fd" ]] && rm -R "$tmp_fd/"
mkdir "$tmp_fd"



#============#
# Processing #
#============#

# Trap
trap "clean_up "$tmp_fd" "$report" "$seq_cds_f" "$seq_aa_f"" SIGINT SIGTERM    # Clean_up function to remove tmp files
wait


#----------#
# GFF file #
#----------#

info "Extraction information from gff"

# Get features from the GFF
myfeatures=$(cut -f 3 "$mygff")

# Check gene feature 
[[ ! $(echo "$myfeatures" | grep -i "gene") ]]      && error "No gene detected in the gff file." 1
[[ $(echo "$myfeatures" | grep -c -i "gene") > 1 ]] && error "More than one gene detected in the gff file." 1

# Update features to get a list of unique items
myfeatures=$(echo "$myfeatures" | sort | uniq)


# Create bed files corresponding to each feature
for i in $myfeatures
do
    awk -v i=$i '$3 == i {print $1 "\t" $4-1 "\t" $5}' "$mygff" > "$tmp_fd/${fn}_${i}.bed"
done


# Store filename of the gene bed file
mygene_ext=$(echo "$myfeatures" | grep -i "gene")
mygene="$tmp_fd/${fn}_${mygene_ext}.bed"

# Store filename of the CDS bed file
mycds_ext=$(echo "$myfeatures" | egrep -i "cds|exon")
[[ $(echo "$mycds_ext" | wc -l) > 1 ]] && error "Exon and CDS features detected. Exiting..." 1
mycds="$tmp_fd/${fn}_${mycds_ext}.bed"


# Create bed file for intron if feature does not exist
if [[ ! $(echo "$myfeatures" | grep -i "intron") ]]
then
    i=$(echo "$myfeatures" | egrep -i "exon|cds")
    paste <(sort -k2n "$tmp_fd/${fn}_${i}.bed" | cut -f 3 | head -n -1) <(sort -k2n "$tmp_fd/${fn}_${i}.bed" | cut -f 2 | tail -n +2) > "$tmp_fd/${fn}_intron.bed"
    sed -i "s/^/$(cut -f 1 "$mygene")\t/g" "$tmp_fd/${fn}_intron.bed"
    myfeatures=$(echo -e "$myfeatures\nintron")
fi


# Update features
myfeatures=$(echo "$myfeatures" | egrep -i "exon|cds|intron|utr")


# Strand orientation
strand=$(grep -i "gene" "$mygff" | cut -f 7)


#--------------------#
# Reference sequence #
#--------------------#

info "Generating reference sequence"

# Extract reference sequence from genome
bedtools getfasta -fi "$mygenome" -bed "$mygene" -fo "$tmp.fa"
sed -ri "s/(>.*):.*$/\1/g" "$tmp.fa"

# Adjust CDS coordinates
mystart=$(cut -f 2 "$mygene")
awk -v i=$mystart '{print $1 "\t" $2-i "\t" $3-i}' "$mycds" > "${tmp}_cds.bed"

# Generate reference CDS
bedtools getfasta -fi "${tmp}.fa" -bed "${tmp}_cds.bed" -fo "${tmp}_cds.fa" &> /dev/null

# Reverse complement the sequence if needed
[[ $strand == - ]] && rev_cpt "${tmp}_cds.fa"

# Translate CDS into protein sequence
transeq -sequence "${tmp}_cds.fa" -outseq "${tmp}_aa.fa" &> /dev/null

# Real gene start
if [[ $strand == - ]]
then
    mystart_gene=$(( $(cut -f 3 "$mygene") - $mystart + 1 ))
fi

# Add reference sequences 
echo -e ">${fn}_ref\n$(cat "${tmp}_cds.fa")"                    > "$seq_cds_f"
echo -e ">${fn}_ref\n$(tail -n +2 "${tmp}_aa.fa" | tr -d "\n")" > "$seq_aa_f"


#-------------#
# Populations #
#-------------#

# Header of the table report
myhdr_pop="All_ref_allele_nb\tAll_alt_allele_nb\tAll_hmz_sample_nb\tAll_htz_sample_nb"

# Columns from VCF file to get genotype from
myclns="10-"

# If different populations
if [[ -n "$mypop_f" ]]
then
    # Population list
    mypops=$(cut -f 2 "$mypop_f" | sort | uniq)
    
    info "Identifying samples related to $(echo $mypops | sed "s/ /, /g") population(s)"
    
    # VCF column header (to find samples)
    myhdr_vcf=$(grep -m 1 "#C" "$myvcf")

    # Get corresponding VCF columns for each population
    for p in $mypops
    do
        myspl=$(awk -v p=$p '$2 == p {print $1}' "$mypop_f" | sed -r "s/(.*)/\^\1\$/g" | tr "\n" "|" | sed "s/|$//g")
        mycln=$(echo "$myhdr_vcf" | tr "\t" "\n" | egrep -n "$myspl" | cut -d ":" -f 1 | tr "\n" "," | sed "s/,$//g")
        myclns=(${myclns[@]} $mycln)

        # Update header of the table report
        myhdr_pop="${myhdr_pop}\t${p}_ref_allele_nb\t${p}_alt_allele_nb\t${p}_hmz_sample_nb\t${p}_htz_sample_nb"
    done
fi


#-----------------------#
# Variants in sequences #
#-----------------------#

# Intersect VCF with GFF to avoid any bad surprise
bedtools intersect -a "$myvcf" -b "$mygene" > "${tmp_fd}/$(basename ${myvcf%.vcf})_${fn}.vcf"
myvcf="${tmp_fd}/$(basename ${myvcf%.vcf})_${fn}.vcf"

# End of the loop (i.e. number of lines in the VCF)
mystep=1
myend=$(wc -l < "$myvcf")

# Determine GT column field
GT_f=$(sed "/#/d" "$myvcf" | head -1 | cut -f 9 | tr ":" "\n" | cut -d ":" -f 1 | grep -n "GT" | cut -d ":" -f 1)

# For each site of the VCF file
while read myline
do
    ProgressBar $mystep $myend ||:

    # Get coordinate (position) of the site
    mypos_genome=$(echo "$myline" | cut -f 2)
    
    # Adjust the coordinate
    mypos=$(( $mypos_genome - $mystart ))
    mypos_sed=$(( $mypos - 1 ))

    # In which feature the site falls in
    awk -v mypos=$mypos_genome '{print $1 "\t" mypos-1 "\t" mypos}' "$mygene" > "${tmp}_pos.bed" # Bed file of the variant
    myregion="-"   # Reset $myregion
    for r in $myfeatures
    do 
        if [[ $(bedtools intersect -a "${tmp}_pos.bed" -b "$tmp_fd/${fn}_${r}.bed") ]]
        then
            myregion=$r
            continue
        fi
    done
    
    # Position within the gene
    [[ $strand == - ]] && mypos_ig=$(( $mystart_gene - $mypos )) || mypos_ig=$mypos


    # Get alleles 
    myref=$(echo "$myline" | cut -f 4)
    myalt_all=$(echo "$myline" | cut -f 5)

    # For each alleles found in the ALT field
    for ((a=1 ; a <= $(echo "$myalt_all" | awk -F "," '{print NF}') ; a++))
    do
        # Select the alternative allele
        myalt=$(echo "$myalt_all" | cut -d "," -f $a)

        # Skip if allele codes for missing data
        [[ $myalt == "*" ]]  && continue

        # Adjust sequence if on the minus strand
        if [[ $strand == - ]]
        then
            myref_rev=$(revseq -notag -sequence <(echo $myref) -filter | tail -n +2 | tr -d "\n")
            myalt_rev=$(revseq -notag -sequence <(echo $myalt) -filter | tail -n +2 | tr -d "\n")
            mymut_tag="${myref_rev}>${myalt_rev}"
        else
            mymut_tag="${myref}>${myalt}"
        fi

        # Gene tag
        mygene_tag="g.${mypos_ig}${mymut_tag}"

        # Default is SNP
        mymut_type="SNP"
        tag=$a
        myalt_tag=$tag
        
        # If INDEL
        [[ $(echo $myref | grep -o "." | wc -l) != $(echo $myalt | grep -o "." | wc -l) && $(echo $myref | wc -m) != 2 ]] && mymut_type="Del" && myalt_tag=$(echo $myref | sed "s/./$tag/g") 
        [[ $(echo $myref | grep -o "." | wc -l) != $(echo $myalt | grep -o "." | wc -l) && $(echo $myalt | wc -m) != 2 ]] && mymut_type="Ins"

        # Test if the reference allele is in the sequence
        if [[ ! $(egrep "^.{$mypos_sed}$myref" <(tail -n +2 "$tmp.fa")) ]]
        then
            warnings=($warnings "Mutation $myalt at $mypos in the gene (genomic position: $mypos_genome) was not identified. This is likely due to out of coordinates mutation.")
            continue
        fi

        # If SNP but long tag
        [[ $mymut_type == "SNP" && $(echo $myref | wc -m) != 2 ]] && myalt_tag=$(echo $myref | sed "s/./$tag/g")

        # If SNP but complex event (FreeBayes specific)
        evt=$(diff <(echo $myref | grep -o .) <(echo $myalt | grep -o .) | grep -c "<" ||:)
        [[ $mymut_type == "SNP" && $evt > 1 ]] && cplx_evt="Yes ($evt)" || cplx_evt="No"

        # Mutate the reference seqence with a tag (the tag is critical for INDELs otherwise the coordinates change in bedtools can't extract the CDS)
        sed -n "1p" "$tmp.fa" > "${tmp}_${mypos_ig}${mymut_tag}.fa"
        sed -n "2p" "$tmp.fa" | sed -r "s/^(.{$mypos_sed})$myref/\1${myalt_tag}/" >> "${tmp}_${mypos_ig}${mymut_tag}.fa"
        
        # Isolate exons to have only CDS
        bedtools getfasta -fi "${tmp}_${mypos_ig}${mymut_tag}.fa" -bed "${tmp}_cds.bed" -fo "${tmp}_${mypos_ig}${mymut_tag}_cds.fa" &> /dev/null

        # Replace tag by the real mutation
        sed -i "s/${myalt_tag}/${myalt}/g" "${tmp}_${mypos_ig}${mymut_tag}_cds.fa"
        sed -i "s/${myalt_tag}/${myalt}/g" "${tmp}_${mypos_ig}${mymut_tag}.fa"

        # Reverse complement the sequence if needed
        [[ $strand == - ]] && rev_cpt "${tmp}_${mypos_ig}${mymut_tag}_cds.fa"
        
        # Translate CDS to get protein sequence
        transeq -sequence "${tmp}_${mypos_ig}${mymut_tag}_cds.fa" -outseq "${tmp}_${mypos_ig}${mymut_tag}_aa.fa" &> /dev/null


        #~~~~~~~~~~~~~#
        # Report line #
        #~~~~~~~~~~~~~#
        
        # Genotype and allele frequency for all and each population
        for ((c=0 ; c < ${#myclns[@]} ; c++))
        do
            
            # Genotype frequency
            alleles=$(echo "$myline" | cut -f ${myclns[$c]} | tr "\t" "\n" | cut -d ":" -f $GT_f)
            hmz=$(echo "$alleles" | egrep -c "$a/$a"     ||:)
            htz=$(echo "$alleles" | egrep -c "0/$a|$a/0" ||:)

            # Allele frequency
            alleles=$(echo  "$alleles" | tr "/" "\n")
            myref_nb=$(echo "$alleles" | grep -c "0"  ||:)
            myalt_nb=$(echo "$alleles" | grep -c "$a" ||:)

            [[ $c == 0 ]] && myGT_ln="${myref_nb}\t${myalt_nb}\t$hmz\t$htz"
            [[ $c > 0 ]]  && myGT_ln="${myGT_ln}\t${myref_nb}\t${myalt_nb}\t$hmz\t$htz"
        done


        # Mutation in protein
        mymut=""    # Reset variable
        if [[ $(echo $myregion | egrep -i "exon|cds") && "$mymut_type" == SNP ]]
        then
            mymut=$(paste \
                <(diff --unchanged-line-format="" --old-line-format="p.%L" --new-line-format="" <(tail -n +2 "${tmp}_aa.fa" | grep -o ".") <(tail -n +2 "${tmp}_${mypos_ig}${mymut_tag}_aa.fa" | grep -o ".")) \
                <(diff --unchanged-line-format="" --old-line-format="" --new-line-format="%dn%L" <(tail -n +2 "${tmp}_aa.fa" | grep -o ".") <(tail -n +2 "${tmp}_${mypos_ig}${mymut_tag}_aa.fa" | grep -o ".")) \
                | sed "s/\t//g" | tr "\n" "," | sed "s/,$//")
        elif [[ $(echo $myregion | egrep -i "exon|cds") && "$mymut_type" =~ "Del"|"Ins" ]]
        then
            mydiff=$(diff <(tail -n +2 "${tmp}_aa.fa" | grep -o ".") <(tail -n +2 "${tmp}_${mypos_ig}${mymut_tag}_aa.fa" | grep -o ".") ||:)

            if [[ $(echo "$mydiff" | head -n 1 | grep "[0-9]*d[0-9]*") ]]
            then
                mypos_del=$(echo "$mydiff" | head -n 1 | cut -d "d" -f 1 ||:)
                myaa=$(echo "$mydiff" | tail -n +2 |  cut -d " " -f 2)
                
                if [[ $(echo "$mypos_del" | grep ",") ]]
                then
                    myaa1=$(echo "$myaa" | head -n 1)
                    myaa2=$(echo "$myaa" | tail -n 1)
                    mypos_del=$(echo "$mypos_del" | sed "s/,/_$myaa2/")
                    myaa=$myaa1
                fi

                mymut="p.${myaa}${mypos_del}del"
            else
                myaa_ref=$(echo "$mydiff" | grep -m 1 "<" | cut -d " " -f 2)
                myaa_alt=$(echo "$mydiff" | grep -m 1 ">" | cut -d " " -f 2)
                myaa_pos=$(echo "$mydiff" | head -n 1 | cut -d "," -f 1 | cut -d "c" -f 1 ||:)
                mystop_pos=$(tail -n +2 "${tmp}_${mypos_ig}${mymut_tag}_aa.fa" | tr -d "\n" | grep -o "." | grep -m 1 -n "*" | cut -d ":" -f 1 ||:)
                mystop_pos=$(( $mystop_pos - $myaa_pos ))
                mymut="p.${myaa_ref}${myaa_pos}${myaa_alt}fsX${mystop_pos}"
            fi
        fi

        
        # Mutation in CDS: if in an exon, is mutation synonymous or not (i.e. if there is a mutation in the protein)
        if [[ $(echo $myregion | egrep -i "exon|cds") ]]
        then
            if [[ "$mymut_type" == "SNP" ]]
            then
                [[ -z "$mymut" ]] && mymut_snp_type="Syn"
                [[ -n "$mymut" ]] && mymut_snp_type="Non_syn" 
            else
                mymut_snp_type="-"
            fi

            # Position in coding sequence
            [[ $strand == - ]] && f=3 || f=2

            mypos_cds=$(bedtools intersect -a "$tmp_fd/${fn}_${myregion}.bed" -b "${tmp}_pos.bed" -wa | cut -f $f)
            mybed_ln=$(awk -v i=$mypos_cds -v f=$f ' $f == i {print NR-1}' "$tmp_fd/${fn}_${myregion}.bed")     # This line is adjusted to avoid the first one
            if [[ $mybed_ln > 0 ]]
            then
                mypos_add=$(sed -n "1,${mybed_ln}p" "$tmp_fd/${fn}_${myregion}.bed" | awk '{SUM += $3-$2} END {print SUM}')
            else
                mypos_add=0
            fi

            [[ $strand != - ]] && mypos_cds=$(( $mypos_genome - $mypos_cds + $mypos_add + 1))
            [[ $strand == - ]] && mypos_cds=$(( $mypos_cds - $mypos_genome + $mypos_add + 1))
            mycds_tag="c.${mypos_cds}${mymut_tag}"
        else
            mymut_snp_type="-"
            mycds_tag="-"
        fi

        # Set empty mymut to a value
        [[ -z "$mymut" ]] && mymut="-"

        # Write the report line
        echo -e "${mypos_genome}\t${mygene_tag}\t${mycds_tag}\t$myregion\t${mymut_type}\t${cplx_evt}\t${mymut_snp_type}\t$mymut\t${myGT_ln}" >> "$report_tmp"


        #~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Mutated sequence storage #
        #~~~~~~~~~~~~~~~~~~~~~~~~~~#

        if [[ $(echo $myregion | egrep -i "exon|cds") ]]
        then
           seq_cds_na=$(echo -e "${seq_cds_na}\n>${fn}_${mycds_tag}|$(cat "${tmp}_${mypos_ig}${mymut_tag}_cds.fa")")
           [[ "$mymut_snp_type" == "Non_syn" ]] && seq_aa=$(echo -e "${seq_aa}\n>${fn}_${mymut}|$(tail -n +2 "${tmp}_${mypos_ig}${mymut_tag}_aa.fa" | tr -d "\n")")
       fi


    done

    # Increase counter for ProgressBar
    ((mystep++))
    
done < "$myvcf"


# Print warnings if any
if [[ -n "$warnings" ]]
then
    for i in "${warnings[@]}"
    do
        warning "$i"
    done
fi


# Header of the table report
echo -e "Genomic_pos\tGene_pos\tCDS_pos\tGene_region\tMutation_type\tComplex_event\tSNP_type\tProtein_mutation\t${myhdr_pop}" > "$report"

[[ $strand != - ]] && cat "$report_tmp" >> "$report"
[[ $strand == - ]] && tac "$report_tmp" >> "$report"

# Generate multifasta files
if [[ $strand == - ]]
then
    echo -e "$seq_cds_na" | tac | sed "/^$/d" | tr "|" "\n"  >> "$seq_cds_f"
    echo -e "$seq_aa"     | tac | sed "/^$/d" | tr "|" "\n"  >> "$seq_aa_f"
else
    echo -e "$seq_cds_na" | sed "/^$/d" | tr "|" "\n" >> "$seq_cds_f"
    echo -e "$seq_aa"     | sed "/^$/d" | tr "|" "\n" >> "$seq_aa_f"
fi

exit 0
