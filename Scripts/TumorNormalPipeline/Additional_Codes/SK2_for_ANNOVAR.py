import gzip
import sys

# Specify the input and output VCF files
input_vcf_file = sys.argv[1]
output_vcf_file = sys.argv[2]

def add_gt_to_vcf(input_file, output_file):
    with gzip.open(input_file, 'rt') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                # Write header lines without modification
                outfile.write(line)
            else:
                # Split the line into columns
                columns = line.strip().split('\t')

                if columns[6] != "PASS":
                    continue

                # Update the format field (FORMAT column)
                format_field = columns[8]  # Assuming GT is the first field in the FORMAT column
                columns[8] = 'GT:' + format_field

                # Update the sample field (Sample columns)
                for i in range(9, len(columns)):
                    sample_info = columns[i].split(':')
                    if i == 9:
                        sample_info.insert(0, '0/0')
                    else:
                        sample_info.insert(0, '0/1')

                    columns[i] = ':'.join(sample_info)

                # Write the modified line to the output file
                outfile.write('\t'.join(columns) + '\n')

# Call the function to add the GT field
add_gt_to_vcf(input_vcf_file, output_vcf_file)
