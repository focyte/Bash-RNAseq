#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <index_path> <splice_sites_file> <gtf_file> <read_type> <data_type>"
    exit 1
fi

# Assign arguments to variables
INPUT_DIR="$1"
OUTPUT_DIR="$2"
INDEX_PATH="$3"
SPLICE_SITES="$4"
GTF_FILE="$5"
READ_TYPE="$6"
DATA_TYPE="$7"


# Run individual scripts

# Choose map script based on DATA_TYPE and READ_TYPE
if [ "$DATA_TYPE" = "Processed" ]; then
    if [ "$READ_TYPE" = "Unpaired" ]; then
        bash mapPU.sh "$INPUT_DIR" "$OUTPUT_DIR" "$INDEX_PATH" "$SPLICE_SITES"
    elif [ "$READ_TYPE" = "Paired" ]; then
        bash mapPP.sh "$INPUT_DIR" "$OUTPUT_DIR" "$INDEX_PATH" "$SPLICE_SITES"
    else
        echo "Invalid read type. Supported types: Unpaired, Paired."
        exit 1
    fi
elif [ "$DATA_TYPE" = "Raw" ]; then
    if [ "$READ_TYPE" = "Paired" ]; then
        bash mapRP.sh "$INPUT_DIR" "$OUTPUT_DIR" "$INDEX_PATH" "$SPLICE_SITES"
    else
        bash mapRU.sh "$INPUT_DIR" "$OUTPUT_DIR" "$INDEX_PATH" "$SPLICE_SITES"
    fi
else
    echo "Invalid data type. Supported types: Processed, Raw."
    exit 1
fi

bash samToBam.sh "$INPUT_DIR" "$OUTPUT_DIR"
bash indexBam.sh "$INPUT_DIR" "$OUTPUT_DIR"
bash featureCount.sh "$INPUT_DIR" "$OUTPUT_DIR" "$GTF_FILE"