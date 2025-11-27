#!/bin/bash
#
# Purpose:
#   Read index 2 FASTQ file with raw barcode reads and correct them
#   according to the provided whitelist.
#
#   For each FASTQ record, return:
#     - read name
#     - raw barcode sequences (CR:Z:bc_raw_seq)
#     - raw barcode qualities (CY:Z:bc_raw_qual)
#     - corrected barcode sequences (CB:z:bc_corrected_seq)
#       (if raw barcode sequence was correctable.)

#### modify based on correct_barcode_from_fastq.sh, remove the parameter `bc_suffix`

set -e
set -o pipefail

decompress_fastq_cat_cmd='cat';
decompress_fastq_zcat_cmd='zcat';
decompress_fastq_igzip_cmd='igzip -c -d';


# Number of threads to use to compress each FASTQ output file.
compress_fastq_threads="${compress_fastq_threads:-4}";

# Gzip compression level for bgzip, pigz and gzip.
compress_fastq_level="${compress_fastq_level:-6}";
# Gzip compression level for igzip (3 is maximum).
compress_fastq_igzip_level="3";

compress_fastq_bgzip_cmd="bgzip -@ ${compress_fastq_threads} -l ${compress_fastq_level} -c";
compress_fastq_pigz_cmd="pigz -p ${compress_fastq_threads} -${compress_fastq_level} -c";
compress_fastq_igzip_cmd="igzip -${compress_fastq_igzip_level} -c";
compress_fastq_gzip_cmd="gzip -${compress_fastq_level} -c";

correct_barcode_from_fastq () {
    local bc_whitelist_filename="${1}";
    local bc_remapping_filename="${2}";
    local fastq_with_raw_bc_filename="${3}";
    local corrected_bc_filename="${4}";
    local max_mismatches="${5:-1}";
    local min_frac_bcs_to_find="${6:-0.5}";

    if [ ${#@} -lt 4 ] ; then
        printf 'Usage:\n';
        printf '    correct_barcode_from_fastq \\\n';
        printf '        bc_whitelist_file bc_remapping_file fastq_with_raw_bc_file \\\n';
        printf '        corrected_bc_file [max_mismatches] [min_frac_bcs_to_find]\n\n';
        printf 'Purpose:\n';
        printf '    Read index 2 FASTQ file with raw barcode reads and correct them\n';
        printf '    according to the provided whitelist.\n\n';
        printf '    For each FASTQ record, return:\n';
        printf '      - read name\n';
        printf '      - raw barcode sequences (CR:Z:bc_raw_seq)\n';
        printf '      - raw barcode qualities (CY:Z:bc_raw_qual)\n';
        printf '      - corrected barcode sequences (CB:z:bc_corrected_seq)\n';
        printf '        (if raw barcode sequence was correctable.)\n\n';
        printf 'Arguments:\n';
        printf '    bc_whitelist_file:\n';
        printf '        File with barcode whitelist to use to correct raw barcode sequences.\n';
        printf '    bc_remapping_file:\n';
        printf '        File with barcodes to use for mapping corrected barcodes to other\n';
        printf '        barcodes (e.g. map 10x multiome ATAC barcodes to 10X multiome RNA\n';
        printf '        barcodes). Set to "false" or "none" to disable remapping.\n';
        printf '    fastq_with_raw_bc_file:\n';
        printf '        FASTQ index 2 file with raw barcode reads to correct.\n';
        printf '    corrected_bc_file:\n';
        printf '        Output file with read name, raw barcode sequence, raw barcode quality\n';
        printf '        and corrected barcode sequence (if correctable) for each FASTQ record\n';
        printf '        in FASTQ index 2 file.\n';
        printf '    max_mismatches\n';
        printf '        Maximum amount of mismatches allowed between raw barcode and whitelist.\n';
        printf '        Default: 1\n';
        printf '    min_frac_bcs_to_find\n';
        printf '        Minimum fraction of reads that need to have a barcode that matches the\n';
        printf '        whitelist.\n';
        printf '        Default: 0.5\n';
        return 1;
    fi

    if type igzip > /dev/null 2>&1 ; then
        # Decompress gzipped FASTQ files with igzip if installed (6x faster than gzip).
        local decompress_fastq_gzipped_cmd="${decompress_fastq_igzip_cmd}";
    else
        # Decompress gzipped FASTQ files with gzip.
        local decompress_fastq_gzipped_cmd="${decompress_fastq_zcat_cmd}";
    fi

    # Detect if input FASTQ files are compressed (gzip/zstd) or not.
    if [ "${fastq_with_raw_bc_filename}" != "${fastq_with_raw_bc_filename%.gz}" ] ; then
        local decompress_fastq_with_raw_bc_cmd="${decompress_fastq_gzipped_cmd}";
    elif [ "${fastq_with_raw_bc_filename}" != "${fastq_with_raw_bc_filename%.zst}" ] ; then
        local decompress_fastq_with_raw_bc_cmd="${decompress_fastq_zstd_cmd}";
    else
        local decompress_fastq_with_raw_bc_cmd="${decompress_fastq_cat_cmd}";
    fi

    if [ ! -e "${bc_whitelist_filename}" ] ; then
        printf 'Error: Barcode whitelist file "%s" could not be found.\n' "${bc_whitelist_filename}" >&2;
        return 1;
    fi

    case "${bc_remapping_filename,,}" in
        ""|false|none)
            ;;
        *)
            if [ ! -e "${bc_remapping_filename}" ] ; then
                printf 'Error: Barcode remapping file "%s" could not be found.\n' "${bc_remapping_filename}" >&2;
                return 1;
            fi;;
    esac

    if [ ! -e "${fastq_with_raw_bc_filename}" ] ; then
        printf 'Error: FASTQ file with raw barcodes "%s" could not be found.\n' "${fastq_with_raw_bc_filename}" >&2;
        return 1;
    fi

    local first_barcode='';

    # Read first barcode from barcode whitelist file.
    if [ "${bc_whitelist_filename%.gz}" == "${bc_whitelist_filename}" ] ; then
        # Uncompressed file.
        read -r first_barcode < "${bc_whitelist_filename}";
    else
        # Unset pipefail.
        set +o pipefail

        # Gzip compressed file.
        first_barcode=$(zcat "${bc_whitelist_filename}" | head -n 1);

        # Set pipefail.
        set -o pipefail
    fi

    # Get length of first barcode.
    local -i bc_length="${#first_barcode}";

    # Check if the first barcode is not empty and only contains ACGT and no other characters.
    if [[ ${bc_length} -eq 0 || "${first_barcode}" != "${first_barcode//[^ACGT]/}" ]] ; then
        printf 'Error: The first line of "%s" does not contain a valid barcode.\n' "${bc_whitelist_filename}";
        return 1;
    fi

    # Get script dir.
    local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";

    /public/home/miaozhu_gibh/software/seq-deploy/bin/seqc run \
        -D bc_length=${bc_length} \
        -release \
        "${script_dir}/correct_barcode_from_fastq.seq" \
            "${bc_whitelist_filename}" \
            "${bc_remapping_filename}" \
            "${fastq_with_raw_bc_filename}" \
            "/dev/stdout" \
            "${corrected_bc_filename}.corrected_bc_stats.tsv" \
            "${max_mismatches}" \
            "${min_frac_bcs_to_find}" \
      | pigz -p 4 \
      > "${corrected_bc_filename}";

    return $?
}

correct_barcode_from_fastq "${@}";
