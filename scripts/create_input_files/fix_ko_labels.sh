#!/bin/bash

INPUT=$1
OUTPUT=$2

awk 'BEGIN{FS=OFS="\t"}
{
    $1 = sprintf("ko%05d", $1)
    print
}' "$INPUT" > "$OUTPUT"