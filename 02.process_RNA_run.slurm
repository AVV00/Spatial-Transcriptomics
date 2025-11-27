#!/bin/bash

#SBATCH --job-name=SpatialOmicsCoord
#SBATCH --output=/public/home/miaozhu_gibh/project/spatial_full_ychen/log/SpatialOmicsCoord.out 
#SBATCH --error=/public/home/miaozhu_gibh/project/spatial_full_ychen/log/SpatialOmicsCoord.out 
#SBATCH --open-mode=append
 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB 
#SBATCH --partition=fat
 
source /public/home/miaozhu_gibh/miniconda3/etc/profile.d/conda.sh 
source /public/home/miaozhu_gibh/miniconda3/bin/activate salus 

main=/public/home/miaozhu_gibh/project/spatial_full_ychen
sampleid=SL08-E13.5

data=$main/data/0221-SL08-E13.5
out=$main/result/$sampleid
species=mouse
result=$main/result
out=$result/$sampleid

if [ ! -d $out ];then
  mkdir -p $out
fi

######################################################################################################
# DATA

if [[ $species == 'mouse' ]];then
  index='/public/home/miaozhu_gibh/reference/refdata-cellranger-arc-mm10-2020-A-2.0.0/star_index'
  tss_bed=/public/home/miaozhu_gibh/reference/refdata-cellranger-arc-mm10-2020-A-2.0.0/regions/tss.bed
  bedFile=/public/home/miaozhu_gibh/reference/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.bed12
else
  index='/public/home/miaozhu_gibh/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star_index'
  tss_bed=/public/home/miaozhu_gibh/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/regions/tss.bed
  bedFile=/public/home/miaozhu_gibh/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.bed12
fi

Exact_bool=TRUE
if [[ $Exact_bool == "TRUE" ]];then 
    SpatialCoord_file=${main}/data/BC_fq/BC_0221/whitelist_zone/correction/Lane02_0221_A1_SpatialCoord.starcode.txt.gz
    whitelist=${main}/data/BC_fq/BC_0221/whitelist_zone/correction/Lane02_0221_A1_whitelist.starcode.txt
else
    SpatialCoord_file=${main}/data/BC_fq/BC_0221/whitelist_zone/Lane02_0221_A1_SpatialCoord.txt.gz
    whitelist=${main}/data/BC_fq/BC_0221/whitelist_zone/Lane02_0221_A1_whitelist.txt
fi

segment_anything=$main/scripts/salus/segment-anything/sam_vit_h_4b8939.pth

######################################################################################################
# MAIN

function main {
    cd $out

    mergefastq \
        $data $data "R1" "R2" "fastq.fq.gz" "fastq_R2.fq.gz"

    fastp_qc \
        $data $data

    # (Optinal)
    cutadapt_RNA \
        $data/cleanFastq $data/cleanFastq
    
    [[ $Exact_bool == "TRUE" ]] && {
        split_fastq \
            $data/cleanFastq $data/cleanFastq
        correct_fastq_barcode_seq \
            $data/cleanFastq $data/cleanFastq/correction
        star_prefix=star_outs_Exact
        Starsolo_mapping \
            $data/cleanFastq/correction $out 8 $star_prefix
    } || {
        star_prefix=star_outs_Simple
        Starsolo_mapping \
            $data/cleanFastq $out 8 $star_prefix
    }

    # append_SpatialCoord \
    #     $out/$star_prefix/Solo.out/GeneFull_Ex50pAS/raw $out/$star_prefix/Solo.out/GeneFull_Ex50pAS/raw

    binGrid \
        $out/$star_prefix/Solo.out/GeneFull_Ex50pAS/raw $out/$star_prefix/simple_grids/${sampleid}_b100 100

    SegAny \
        $out/$star_prefix/simple_grids/${sampleid}_b100 $out/$star_prefix/segmentation 100

}

######################################################################################################

function mergefastq () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: merge fastq ------\n"

    local in_dir=$1
    local out_dir=$2
    local num=("${@:3:2}")
    local match_pattern=("${@:5:2}")

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir $out_dir

    {
        echo $(date)
        for idx in ${!num[@]}; do
            [ -s  ${out_dir}/${sampleid}_S1_L001_${num[idx]}_001.fastq.gz ] && {
                echo "NOTE: $out_dir/${sampleid}_S1_L001_${num[idx]}_001.fastq.gz exists!!!"
                continue
            }

            for file in `ls ${in_dir}/*/*_${match_pattern[idx]}` ; do
                [ -s $file ] && {
                    echo ${file}
                    cat ${file} >> ${out_dir}/${sampleid}_S1_L001_${num[idx]}_001.fastq.gz
                } || {
                    echo "ERROR: $file isn't exists!!!"
                    exit 1
                }
            done
            echo "merge file: "${out_dir}/${sampleid}_S1_L001_${num[idx]}_001.fastq.gz
        done
    } 2>&1 

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: merge fastq ------\n"
}

function fastp_qc () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: fastp ------\n"

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir $out_dir

    {
        echo $(date)
        [ -e  $out_dir/cleanFastq ] && {
            echo "NOTE: $out_dir/cleanFastq exists!!!"
            exit 1
        }
    
        [ -s $in_dir/*R1*.gz -a -s $in_dir/*R2*.gz ] && {
            local cmd="bash $main/scripts/salus/fastp.sh $sampleid $in_dir $out_dir"
            echo $cmd
            eval $cmd
        } || {
            echo "ERROR: $in_dir/*R1*.gz or $in_dir/*R2*.gz isn't exists!!!"
            exit 1
        }

    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: fastp ------\n"
}

function cutadapt_RNA () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: cutadapt ------\n"

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir $out_dir

    {
        echo $(date)
        [ -e $out_dir/${sampleid}.cutadapt.fq.gz ] && {
            echo "NOTE: ${sampleid}.cutadapt.fq.gz exists!!!"
            exit 1
        }
    
        [ -s $in_dir/${sampleid}*R1*clean.fq.gz -a -s $in_dir/${sampleid}*R2*clean.fq.gz ] && {
            local cmd="cutadapt -a 'AAAAAAAAAAAAAAANNNNNNNNNNGCGAGAGTCGAGGGTGCTGTAGTCACAAGA;anywhere' \
                -a 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT;anywhere' \
                -a 'TATGACGCCAGACATGTTGCGAACGGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT;anywhere' \
                -g 'GTTCGTCTTCTGCCGTATGCTCTACAGAGCACCCGAGCCTCTACACATATTCTCTGTC;anywhere' \
                --error-rate 0.11 --overlap 24 -Z --minimum-length 30 \
                --action trim --report full --times 8 --poly-a -j 4\
                -o $out_dir/${sampleid}.R1.cutadapt.fq.gz -p $out_dir/${sampleid}.R2.cutadapt.fq.gz\
                --json $out_dir/${sampleid}.cutadapt.json $in_dir/${sampleid}*R1*clean.fq.gz $in_dir/${sampleid}*R2*clean.fq.gz"
            echo $cmd
            eval $cmd
        } || {
            echo "ERROR: $in_dir/${sampleid}*clean.fq.gz isn't exists!!!"
            exit 1
        }

    } 2>&1 | tee $out_dir/${sampleid}.cutadapt.qc.txt

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: cutadapt ------\n"
}

function Starsolo_mapping () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: STARsolo ------\n"

    local in_dir=$1
    local out_dir=$2
    local thread=$3
    local prefix=$4

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir -p $out_dir

    {
        echo $(date)
        [ -s  $out_dir/$prefix/Aligned.sortedByCoord.out.bam ] && {
            echo "NOTE: $out_dir/$prefix/Aligned.sortedByCoord.out.bam exists!!!"
            exit 1
        }
    
        [ -s ${in_dir}/*R1*.gz -a -s ${in_dir}/*R2*.gz ] && {
            # [[ $Exact_bool == "TRUE" ]] && {
            #     local matchWLtype="Exact"
            # } || {
            #     local matchWLtype="1MM_multi_Nbase_pseudocounts"
            # }
            local matchWLtype="1MM_multi_Nbase_pseudocounts"

            local cmd="/public/home/miaozhu_gibh/software/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
                --runThreadN $thread \
                --genomeDir $index \
                --readFilesCommand zcat \
                --outFileNamePrefix $prefix/ \
                --readFilesIn $in_dir/${sampleid}*R2*.gz $in_dir/${sampleid}*R1*.gz \
                --soloType CB_UMI_Simple \
                --soloCBstart 1 --soloCBlen 30 \
                --soloUMIstart 31 --soloUMIlen 10 \
                --soloBarcodeReadLength 0 \
                --soloCBwhitelist $whitelist \
                --soloCBmatchWLtype $matchWLtype \
                --soloUMIdedup 1MM_CR \
                --soloMultiMappers Unique \
                --clip3pAdapterSeq AAAAAAAAAAA \
                --clip3pAdapterMMp 0.1 \
                --clipAdapterType Hamming \
                --soloUMIfiltering MultiGeneUMI_CR \
                --soloCellFilter None \
                --soloStrand Forward \
                --soloFeatures GeneFull_Ex50pAS \
                --outFilterMultimapNmax 20 \
                --limitOutSJcollapsed 2000000 \
                --outFilterScoreMinOverLread 0 \
                --outFilterMatchNminOverLread 0 \
                --outFilterMatchNmin 25 \
                --outReadsUnmapped Fastx \
                --outSAMattributes CR UR CY UY CB UB NH HI AS nM MD GX GN \
                --outSAMtype BAM SortedByCoordinate"
            echo $cmd
            eval $cmd
        } || {
            echo "ERROR: $in_dir/*R1*.gz or $in_dir/*R2*.gz isn't exists!!!"
            exit 1
        }            
    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: STARsolo ------\n"
}

function split_fastq () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` '\n------ START: Split fastqs to small files ------\n' 

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir -p $out_dir

    if [ ! -d $out_dir/smallfastqs ];then
        mkdir -p $out_dir/smallfastqs
    fi

    dosplitfull() {
        seqkit split2 $1 -s 8000000 -O $2/smallfastqs -e .gz
    }

    export -f dosplitfull

    local cmd="parallel --delay 1 dosplitfull ::: ${in_dir}/*R1*.fq.gz  ::: ${out_dir}"
    echo $cmd
    eval $cmd

    ls -d $out_dir/smallfastqs/* | grep "R1.*part" > $out_dir/filesr1.xls

    echo -e `date +"%Y-%m-%d %H:%M:%S"` '\n------ END: Split fastqs to small files ------\n'  
}

function correct_fastq_barcode_seq () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` '\n------ START: add corrected barcode for R1 fastq ------\n' 

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir -p $out_dir

    process_fastq() {
        base=$(basename $1)
        out_dir=$(dirname $1)
        prefix=${base%%.fq.gz}

        bash ${3}/scripts/salus/correct_barcode_from_fastq.sh \
            ${2} \
            false \
            $1 \
            ${out_dir}/${prefix}_correct_barcode.fq.gz \
            3
    }

    export -f process_fastq

    local cmd="parallel --jobs 10 --colsep '\t' process_fastq :::: $in_dir/filesr1.xls ::: ${whitelist} ::: ${main}"
    echo $cmd
    eval $cmd

    local cmd="python $main/scripts/salus/correct_barcode_summary.py \"${in_dir}/smallfastqs/*_bc_stats.tsv\" ${in_dir}/corrected_bc_stats.tsv "
    echo $cmd
    eval $cmd

    parallel 'if [ -f '${out_dir}'/{}_add_corrected_barcode_R1.fq.gz ]; then echo "found {}_add_corrected_barcode_R1.fq.gz or {} folder";
    else ls '${in_dir}'/smallfastqs/{}*_correct_barcode.fq.gz | xargs cat > '${out_dir}'/{}_add_corrected_barcode_R1.fq.gz && echo "Generated {}_add_corrected_barcode_R1.fq.gz"; fi' ::: ${sampleid}

    ln -s ${in_dir}/*R2*.fq.gz ${out_dir}/${sampleid}_add_corrected_barcode_R2.fq.gz

    echo -e `date +"%Y-%m-%d %H:%M:%S"` '\n------ END: add corrected barcode for R1 fastq ------\n'  
}

function correct_fastq_barcode () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: correct_fastq_barcode ------\n"

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir -p $out_dir

    {
        echo $(date)
        # [ -e $out_dir/corrected ] && {
        #     echo "NOTE: $out_dir/corrected exists!!!"
        #     exit 1
        # }
    
        [ -s $whitelist -a -e $in_dir ] && {
            local cmd="python3 $main/scripts/salus/correct_barcode.py \
                --whitelist $whitelist \
                --read1 $in_dir/${sampleid}*R1*.gz \
                --read2 $in_dir/${sampleid}*R2*.gz \
                --outname $out_dir/$sampleid \
                --distance 3"
            echo $cmd
            eval $cmd
        } || {
            echo "ERROR: $whitelist or $in_dir isn't exists!!!"
            exit 1
        }

    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: correct_fastq_barcode ------\n"
}

function append_SpatialCoord () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: append_SpatialCoord ------\n"

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir $out_dir

    {
        echo $(date)
        [ -s  $out_dir/spatial.txt ] && {
            echo "NOTE: $out_dir/spatial.txt exists!!!"
            exit 1
        }
    
        [ -s $in_dir/barcodes.tsv -a -s $SpatialCoord_file ] && {
            zcat $SpatialCoord_file | awk 'FNR==NR {{gsub("N","A",$1);a[$1]=$2" "$3; next}} $1 in a {{print $1, a[$1]}}' - $in_dir/barcodes.tsv > $out_dir/spatial.txt
        } || {
            echo "ERROR: $SpatialCoord_file or $in_dir/barcodes.tsv isn't exists!!!"
            exit 1
        }

    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: append_SpatialCoord ------\n"
}

function binGrid () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: binGrid ------\n"

    local in_dir=$1
    local out_dir=$2
    local binsize=$3

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir -p $out_dir

    {
        echo $(date)
        [ -s  $out_dir/matrix.mtx.gz ] && {
            echo "NOTE: $out_dir/matrix.mtx.gz exists!!!"
            exit 1
        }
    
        [ -s $SpatialCoord_file -a -e $in_dir ] && {
            local cmd="python3 $main/scripts/salus/splanegridict.py \
                -d \
                -b $binsize \
                -i $SpatialCoord_file \
                -r $in_dir \
                -o $out_dir"
            echo $cmd
            eval $cmd
        } || {
            echo "ERROR: $SpatialCoord_file or $in_dir isn't exists!!!"
            exit 1
        }

    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: binGrid ------\n"
}

function SegAny () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: SegAny ------\n"

    local in_dir=$1
    local out_dir=$2
    local binsize=$3

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir -p $out_dir

    {
        echo $(date)
        [ -s $out_dir/${sampleid}_b${binsize}/mask.npy ] && {
            echo "NOTE: $out_dir/${sampleid}_b${binsize}/mask.npy exists!!!"
            exit 1
        }
    
        [ -s $segment_anything -a -e $in_dir ] && {
            local cmd="python3 $main/scripts/salus/SegAny.py \
                -p $in_dir \
                -o $out_dir/pre_${sampleid}_b${binsize} \
                -d $segment_anything \
                -b $binsize"
            echo $cmd
            eval $cmd

            [ -s $out_dir/pre_${sampleid}_b${binsize}/mask.npy ] && {
                local cmd="python3 $main/scripts/salus/SegAny.py \
                    -p $in_dir \
                    -o $out_dir/${sampleid}_b${binsize} \
                    -d $segment_anything \
                    -b $binsize \
                    -k $out_dir/pre_${sampleid}_b${binsize}/mask.npy"
                echo $cmd
                eval $cmd
            } || {
                echo "ERROR: $out_dir/pre_${sampleid}_b${binsize}/mask.npy isn't exists!!!"
                exit 1
            }
        } || {
            echo "ERROR: $segment_anything or $in_dir isn't exists!!!"
            exit 1
        }

    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: SegAny ------\n"
}

######################################################################################################
# RUN
main 2>&1 >> ${main}/log/${sampleid}_RNA.log