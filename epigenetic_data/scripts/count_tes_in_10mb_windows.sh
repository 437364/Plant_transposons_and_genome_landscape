#!/bin/bash

# Enhanced script to count LTR-TEs in 1Mb genomic windows with summary statistics
# Outpu used in notebooks/element_density_Lsyl_Jeff.ipynb for density comparisons between Lsyl and Jeff, boxplots per family per element type, and centromeric element proportion analysis.
# Create output directories
mkdir -p data/local_patterns/Jeff/TE_counts_by_window
mkdir -p data/local_patterns/Lsyl/TE_counts_by_window

# Function to process BED files and generate summaries
count_TEs_in_windows() {
    local species=$1
    local window_file=$2
    local te_bed_dir=$3
    local output_dir=$4
    
    echo "Processing ${species} TEs..."
    
    # Create summary file
    summary_file="${output_dir}/summary_statistics.txt"
    echo -e "TE_Family\tCategory\tTotal_TEs\tWindows_with_TEs\tMean_TEs_per_window\tMax_TEs_per_window" > "${summary_file}"
    
    # Find all BED files for this species
    find "${te_bed_dir}" -name "*.bed" -type f | sort | while read te_file; do
        # Extract category and family info from path
        category=$(echo "${te_file}" | grep -oP '(?<=_bed_files/)[^/]+')
        family=$(basename "${te_file}" .bed)
        
        # Get relative path for output naming
        rel_path=$(echo "${te_file}" | sed "s|${te_bed_dir}/||" | sed 's|/|_|g' | sed 's|\.bed$||')
        output_file="${output_dir}/${rel_path}_counts.bed"
        
        echo "  Processing: ${category}/${family}"
        
        # Count overlaps
        bedtools intersect -a "${window_file}" -b "${te_file}" -c > "${output_file}"
        
        # Calculate summary statistics
        total_elements=$(wc -l < "${te_file}")
        windows_with_tes=$(awk '$5 > 0' "${output_file}" | wc -l)
        mean_per_window=$(awk '{sum+=$5; count++} END {if(count>0) print sum/count; else print 0}' "${output_file}")
        max_per_window=$(awk 'BEGIN{max=0} {if($5>max) max=$5} END{print max}' "${output_file}")
        
        # Append to summary
        echo -e "${family}\t${category}\t${total_elements}\t${windows_with_tes}\t${mean_per_window}\t${max_per_window}" >> "${summary_file}"
    done
    
    echo "Summary statistics saved to: ${summary_file}"
}

# Process Jeff TEs
echo "================================"
echo "Processing Jeff genome"
echo "================================"
count_TEs_in_windows \
    "Jeff" \
    "data/local_patterns/Jeff_bins_10Mb.bed" \
    "data/local_patterns/Jeff/transposon_bed_files" \
    "data/local_patterns/Jeff/TE_counts_by_window"

echo ""
echo "================================"
echo "Processing Lsyl genome"
echo "================================"
count_TEs_in_windows \
    "Lsyl" \
    "data/local_patterns/Lsyl_bins_10Mb.bed" \
    "data/local_patterns/Lsyl/transposon_bed_files" \
    "data/local_patterns/Lsyl/TE_counts_by_window"

echo ""
echo "================================"
echo "All done!"
echo "================================"
echo "Results are in:"
echo "  - data/local_patterns/Jeff/TE_counts_by_window/"
echo "  - data/local_patterns/Lsyl/TE_counts_by_window/"
echo ""
echo "Summary statistics:"
echo "  - data/local_patterns/Jeff/TE_counts_by_window/summary_statistics.txt"
echo "  - data/local_patterns/Lsyl/TE_counts_by_window/summary_statistics.txt"
