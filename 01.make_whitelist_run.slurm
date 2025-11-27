#!/bin/bash

#SBATCH --job-name=SpatialOmicsCoord
#SBATCH --output=/public/home/miaozhu_gibh/project/spatial_full_ychen/log/SpatialOmicsCoord.out 
#SBATCH --error=/public/home/miaozhu_gibh/project/spatial_full_ychen/log/SpatialOmicsCoord.out 
#SBATCH --open-mode=append
 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB 
#SBATCH --partition=fat 
 
source /public/home/miaozhu_gibh/miniconda3/etc/profile.d/conda.sh 
source /public/home/miaozhu_gibh/miniconda3/bin/activate salus 

main=/public/home/miaozhu_gibh/project/spatial_full_ychen
sampleid=Lane02_0221
data=$main/data/BC_fq/BC_0221
out=$data

######################################################################################################
# MAIN

function main {
    cd $out

    mregefastq \
        $data $data
    
    fastp_qc \
        $data $data

    # (Optional)
    cutadapt_BC \
        $data/cleanFastq $data/cleanFastq

    convert_SpatialCoord \
        $data/cleanFastq $data/cleanFastq

    fastqc_qc \
        $data/cleanFastq $data/cleanFastq

    extract_whitelist_zone \
        $data/cleanFastq $data/whitelist_zone

    starcode_correction \
        $data/whitelist_zone $data/whitelist_zone/correction

    estimate_tissueBoundary \
        $data/whitelist_zone $data/estimate_tissueBoundary

    remove_duplicate \
        $data/whitelist_zone $data/whitelist_zone 
}

######################################################################################################

function mregefastq () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: merge fastq ------\n"

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir $out_dir

    {
        echo $(date)
        [ -s  ${out_dir}/${sampleid}.fq.gz ] && {
            echo "NOTE: ${out_dir}/${sampleid}.fq.gz exists!!!"
            exit 1
        }

        for file in `ls ${in_dir}/*${sampleid}*.{fastq,fq}.gz` ; do
            [ -s $file ] && {
                echo ${file}
                cat ${file} >> ${out_dir}/${sampleid}.fq.gz
            } || {
                echo "ERROR: $file isn't exists!!!"
                exit 1
            }
        done
        echo "merge file: "${out_dir}/${sampleid}.fq.gz
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
        # [ -e  $out_dir/cleanFastq ] && {
        #     echo "NOTE: $out_dir/cleanFastq exists!!!"
        #     exit 1
        # }
    
        [ -s $in_dir/${sampleid}*.fastq.gz -o -s $in_dir/${sampleid}*.fq.gz ] && {
            local cmd="bash $main/scripts/salus/fastp_whitelist.sh $sampleid $in_dir $out_dir"
            echo $cmd
            eval $cmd
        } || {
            echo "ERROR: $in_dir/${sampleid}*.{fastq,fq}.gz isn't exists!!!"
            exit 1
        }

    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: fastp ------\n"
}

function cutadapt_BC () {
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
    
        [ -s $in_dir/${sampleid}*clean.fq.gz ] && {
            local cmd="cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
                --error-rate 0.11 --overlap 5 -Z --action lowercase --report full --times 1 -j 4\
                -o $out_dir/${sampleid}.cutadapt.fq.gz \
                --json $out_dir/${sampleid}.cutadapt.json $in_dir/${sampleid}*clean.fq.gz"
            echo $cmd
            eval $cmd
        } || {
            echo "ERROR: $in_dir/${sampleid}*clean.fq.gz isn't exists!!!"
            exit 1
        }

    } 2>&1 | tee $out_dir/${sampleid}.cutadapt.qc.txt

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: cutadapt ------\n"
}

function convert_SpatialCoord () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: SpatialOmicsCoord ------\n"

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir $out_dir

    {
        echo $(date)
        [ -s  $out_dir/newCoord_${sampleid}*.gz ] && {
            echo "NOTE: $out_dir/newCoord_${sampleid}.clean.fq.gz exists!!!"
            exit 1
        }
    
        [ -s $in_dir/*${sampleid}*clean.fq.gz ] && {
            local cmd="python3 $main/scripts/salus/SpatialOmicsCoord_0221.py \
                $in_dir/$sampleid*clean.fq.gz \
                $out_dir \
                2160 4096 -r1 0.3 -r2 0.2 -r1s 0.25"
            echo $cmd
            eval $cmd
        } || {
            echo "ERROR: $in_dir/${sampleid}.clean.fq.gz isn't exists!!!"
            exit 1
        }

    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: SpatialOmicsCoord ------\n"
}

function fastqc_qc () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: fastqc ------\n"

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir $out_dir

    {
        echo $(date)
        [ -s  $in_dir/newCoord_${sampleid}*fastqc.html ] && {
            echo "NOTE: $in_dir/newCoord_${sampleid}*fastqc.html exists!!!"
            exit 1
        }
    
        [ -s $in_dir/newCoord_${sampleid}*.gz ] && {
            local cmd="fastqc -t 4 -f fastq  \
                -o $out_dir \
                $in_dir/newCoord_${sampleid}*.gz"
            echo $cmd
            eval $cmd

            local cmd="multiqc $out_dir/newCoord* -n ${sampleid}_multiqc -o $out_dir"
            echo $cmd
            eval $cmd
        } || {
            echo "ERROR: $in_dir/newCoord_${sampleid}.clean.fq.gz isn't exists!!!"
            exit 1
        }

    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: fastqc ------\n"
}

function extract_whitelist_zone () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: extract each zone SpatialCoord ------\n"

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir $out_dir

    {
        echo $(date)
        [ -s $out_dir/${sampleid}_*_whitelist.txt -o -s $out_dir/${sampleid}_*_SpatialCoord.txt.gz ] && {
            echo "NOTE: zone whitelsit exists!!!"
            exit 1
        }

        [ -s $in_dir/newCoord_${sampleid}*clean.fq.gz ] && {
            local cmd="python3 $main/scripts/salus/ChipZone_whitelist.py \
                            $in_dir/newCoord_$sampleid*clean.fq.gz \
                            $out_dir \
                            $sampleid \
                            1"
            echo $cmd
            eval $cmd
        } || {
            echo "ERROR: $in_dir/newCoord_${sampleid}.clean.fq.gz isn't exists!!!"
            exit 1        
        }
    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: extract each zone SpatialCoord ------\n"
}

function estimate_tissueBoundary () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: estimate tissueBoundary ------\n"

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir $out_dir

    {
        echo $(date)
        for file in $in_dir/${sampleid}_*_SpatialCoord.txt.gz ; do
            local f=$(basename $file)
            f=${f/_SpatialCoord.txt.gz/}
            [ -s $out_dir/${f}  ] && {
                echo "NOTE: $out_dir/${f} exists!!!"
                continue
            }

            [ -s $file ] && {
                local cmd="python3 $main/scripts/salus/estBound_makeBins.py \
                            -f $file \
                            -b 100 \
                            -o $out_dir/$f"
                echo $cmd
                eval $cmd
            } || {
                echo "ERROR: $file isn't exists!!!"
                exit 1        
            }
        done
    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: estimate tissueBoundary ------\n"
}

function starcode_correction () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: starcode_correction ------\n"

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir $out_dir

    {
        echo $(date)
        for file in $in_dir/${sampleid}_A1_whitelist.txt ; do
            local f=$(basename $file)
            out_file=${f/.txt/.starcode.txt}
            # [ -s $out_dir/${out_file}  ] && {
            #     echo "NOTE: $out_dir/${out_file} exists!!!"
            #     continue
            # }

            [ -s $file ] && {
                # local cmd="/public/home/miaozhu_gibh/software/starcode-1.4 \
                #     -t 20 -d 3 -s --non-redundant \
                #     -i $file -o $out_dir/${out_file}
                # "
                # echo $cmd
                # eval $cmd

                SpatialCoord=${f/whitelist.txt/SpatialCoord.txt.gz}
                corrected_SpatialCoord=${SpatialCoord/txt.gz/starcode.txt}

                zcat  $in_dir/$SpatialCoord | awk 'FNR==NR {{a[$1]=$1; next}} $1 in a {{print $0}}' $out_dir/${out_file} - > $out_dir/${corrected_SpatialCoord}
                bgzip -@ 4 $out_dir/${corrected_SpatialCoord}
            } || {
                echo "ERROR: $file isn't exists!!!"
                exit 1        
            }
        done
    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: starcode_correction ------\n"
}

function remove_duplicate () {
    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ START: remove whitelist duplicate ------\n"

    local in_dir=$1
    local out_dir=$2

    [ -e $in_dir ] || {
        echo "ERROR: no input"
        exit 1
    }
    [ -e $out_dir ] || mkdir $out_dir

    {
        echo $(date)
        for file in $in_dir/${sampleid}_*_whitelist.txt ; do
            local f=$(basename $file)
            f=${f/.txt/.rmdup.txt}
            [ -s $out_dir/${f}  ] && {
                echo "NOTE: $out_dir/${f} exists!!!"
                continue
            }

            [ -s $file ] && {
                local cmd="sort $file | uniq -u > $out_dir/${f}"
                echo $cmd
                eval $cmd
                [ -s $out_dir/${f} ] && {
                    rm $file
                }
            } || {
                echo "ERROR: $file isn't exists!!!"
                exit 1        
            }
        done
    } 2>&1

    echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ END: remove whitelist duplicate ------\n"
}

######################################################################################################
# RUN
main 2>&1 >> ${main}/log/${sampleid}_whitelist.log