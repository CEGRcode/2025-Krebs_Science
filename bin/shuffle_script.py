import random
import sys

def shuffle_row(row):
    """Shuffle the numbers in a row while preserving the first two columns."""
    header = row[:2]  # Keep the first two columns unchanged
    numbers = row[2:]  # The columns to shuffle
    random.shuffle(numbers)  # Shuffle the numbers
    return header + numbers

def process_file(input_filename, output_filename):
    with open(input_filename, 'r') as infile:
        lines = infile.readlines()
    
    # Ensure the file has at least 1 header line
    if len(lines) < 1:
        print("Error: The input file must contain at least one line.")
        sys.exit(1)
    
    # Extract the header line
    header_line = lines[0].strip()
    header = header_line.split()
    
    # Process the remaining rows
    rows = [line.strip().split() for line in lines[1:]]
    shuffled_rows = [shuffle_row(row) for row in rows]
    
    with open(output_filename, 'w') as outfile:
        outfile.write('\t'.join(header) + '\n')
        for row in shuffled_rows:
            outfile.write('\t'.join(row) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python shuffle_script.py <input_file> <output_file>")
        sys.exit(1)

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    process_file(input_filename, output_filename)
    print(f'Processed {input_filename} and saved the output to {output_filename}.')
