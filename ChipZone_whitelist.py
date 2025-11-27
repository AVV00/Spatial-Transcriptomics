import argparse
import gzip
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Argument parsing
parser = argparse.ArgumentParser(description='Process fastq.gz file based on defined ChipZones')
parser.add_argument('input_file', type=str, help='Input directory path')
parser.add_argument('out_dir', type=str, help='Output directory path')
parser.add_argument('sampleid', type=str, help='Sample ID')
parser.add_argument('IDy', type=str, help='IDy parameter (default: "1")')
args = parser.parse_args()

# def reverse_complement(sequence):
#     complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
#     return ''.join(complement.get(base, base) for base in reversed(sequence))

ChipZoneEdges = {}
# for IDy in ('1', '2', ''):
IDy = args.IDy
if IDy == '':
    xTuple = map(chr, range(65, 69))
    yRange = (0, 45000)  # maxY might be 44524
else:
    xTuple = map(chr, range(65, 70))
    if IDy != '':
        yRange = (1856, 42560)
        yRange = (1800, 42599)
    for IDx in xTuple:
        ZoneID = IDx + IDy
        if IDx == 'A' and IDy != '':
            xRange = (3891, 44596)
            xRange = (3800, 44599)
        elif IDx == 'B' and IDy != '':
            xRange = (48521, 89226)
            xRange = (48500, 89299)
        elif IDx == 'C' and IDy != '':
            xRange = (93151, 133856)
            xRange = (93100, 133899)
        elif IDx == 'D' and IDy != '':
            xRange = (137781, 178486)
            xRange = (137700, 178499)
        elif IDx == 'E' and IDy != '':
            xRange = (182410, 223115)
            xRange = (182400, 223199)
        elif IDx == 'A' and IDy == '':
            xRange = (3534, 55502)
        elif IDx == 'B' and IDy == '':
            xRange = (59524, 111492)
        elif IDx == 'C' and IDy == '':
            xRange = (115515, 167483)
        elif IDx == 'D' and IDy == '':
            xRange = (171505, 223473)
        else:
            raise ValueError("[x]It can only be memory error...")
        
        ChipZoneEdges[ZoneID] = (xRange, yRange)

# # Process data for each zone and save the split files
# for zone_id, (x_range, y_range) in ChipZoneEdges.items():
#     with gzip.open(f'{args.out_dir}/{args.sampleid}_{zone_id}_SpatialCoord.txt.gz', 'w') as spatialcoord_file, \
#          open(f'{args.out_dir}/{args.sampleid}_{zone_id}_whitelist.txt', 'w') as whitelist_file:
         
#         with gzip.open(args.input_file, 'rt') as input_file :
#             for title, seq, qual in FastqGeneralIterator(input_file):
#                 parts = title.split('_')
#                 x, y = map(float, parts[-2:])
#                 rc_sequence = str(Seq(seq).reverse_complement())

#                 if x_range[0] <= x <= x_range[1] and y_range[0] <= y <= y_range[1]:
#                     spatialcoord_file.write(f"{rc_sequence}\t{x}\t{y}\t{zone_id}\n".encode())
#                     whitelist_file.write(f"{rc_sequence}\n")

# Open files outside the loop
file_handles = {}
for zone_id, (x_range, y_range) in ChipZoneEdges.items():
    spatialcoord_file = gzip.open(f'{args.out_dir}/{args.sampleid}_{zone_id}_SpatialCoord.txt.gz', 'at')
    whitelist_file = open(f'{args.out_dir}/{args.sampleid}_{zone_id}_whitelist.txt', 'a')
    file_handles[zone_id] = (spatialcoord_file, whitelist_file)

with gzip.open(args.input_file, 'rt') as input_file:
    for title, seq, qual in FastqGeneralIterator(input_file):
        parts = title.split('_')
        x, y = map(float, parts[-2:])
        rc_sequence = str(Seq(seq).reverse_complement())
        
        for zone_id, (x_range, y_range) in ChipZoneEdges.items():
            if x_range[0] <= x <= x_range[1] and y_range[0] <= y <= y_range[1]:
                spatialcoord_file, whitelist_file = file_handles[zone_id]
                spatialcoord_file.write(f"{rc_sequence}\t{x}\t{y}\t{zone_id}\n")
                whitelist_file.write(f"{rc_sequence}\n")

# Close all file handles at the end
for spatial_file, whitelist in file_handles.values():
    spatial_file.close()
    whitelist.close()
