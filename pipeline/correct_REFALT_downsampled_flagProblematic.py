#!/usr/bin/env python3
"""
Batch correct REF/ALT - SLURM Array Job Version
"""

import pandas as pd
import pysam
from pathlib import Path
import gzip
from tqdm import tqdm
import logging
import sys
from datetime import datetime
import shutil
import os

# Get SLURM array task ID
SLURM_ARRAY_TASK_ID = int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
SLURM_ARRAY_TASK_COUNT = int(os.environ.get('SLURM_ARRAY_TASK_COUNT', 1))

# Setup logging (same as before)
def setup_logging(output_base_dir, task_id):
    """Setup logging with task ID"""
    log_dir = Path(output_base_dir) / "logs"
    log_dir.mkdir(exist_ok=True, parents=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"ref_correction_task{task_id}_{timestamp}.log"
    
    file_formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_formatter = logging.Formatter(
        '%(levelname)s - %(message)s'
    )
    
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_formatter)
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(console_formatter)
    
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return log_file


def correct_single_file(input_file, reference_fasta, output_dir):
    """
    Correct a single file with is_complex flag
    """
    try:
        sample_id = input_file.stem.replace('.tumour.tso500_sparse_downsample.txt', '')
        if sample_id.endswith('.txt'):
            sample_id = sample_id[:-4]
        
        output_file = output_dir / f"{sample_id}.tumour.tso500_sparse_downsample_corrected.txt.gz"
        
        if output_file.exists():
            logging.info(f"Skipping {sample_id} - already processed")
            return {'sample_id': sample_id, 'status': 'skipped', 'n_variants': 0, 'n_swapped': 0, 'n_complex': 0}
        
        # Load data
        if str(input_file).endswith('.gz'):
            df = pd.read_csv(input_file, sep='\t', compression='gzip')
        else:
            df = pd.read_csv(input_file, sep='\t')
        
        # Load reference genome
        fasta = pysam.FastaFile(reference_fasta)
        
        corrected_rows = []
        n_swapped = 0
        n_complex = 0
        
        for idx, row in df.iterrows():
            chrom = row['CHROM']
            pos = row['POS']
            current_ref = row['REF']
            current_alt = str(row['ALT']).split(',')[0]
            ref_count = row['REFcounts']
            alt_count = row['ALTcounts']
            total_dp = row['DP']
            
            # Get true reference
            try:
                true_ref = fasta.fetch(chrom, pos - 1, pos).upper()
            except:
                # If can't fetch reference, keep original and flag as complex
                corrected_rows.append({
                    'CHROM': chrom,
                    'POS': pos,
                    'REF': current_ref,
                    'ALT': current_alt,
                    'DP': total_dp,
                    'REFcounts': ref_count,
                    'ALTcounts': alt_count,
                    'BAF_ALT': alt_count / total_dp if total_dp > 0 else 0,
                    'GT_inferred': './.',
                    'is_complex': True,  # ← FLAG
                    'correction_type': 'fetch_failed'
                })
                n_complex += 1
                continue
            
            # ============================================================
            # CASE 1: REF is correct
            # ============================================================
            if current_ref == true_ref:
                new_ref, new_alt = current_ref, current_alt
                new_ref_count, new_alt_count = ref_count, alt_count
                is_complex = False
                correction_type = 'correct'
            
            # ============================================================
            # CASE 2: REF and ALT are swapped
            # ============================================================
            elif current_alt == true_ref:
                new_ref, new_alt = current_alt, current_ref
                new_ref_count, new_alt_count = alt_count, ref_count
                is_complex = False
                correction_type = 'swapped'
                n_swapped += 1
            
            # ============================================================
            # CASE 3: COMPLEX - neither matches (FLAGGED)
            # ============================================================
            else:
                new_ref = true_ref
                if ref_count >= alt_count:
                    new_alt, new_ref_count, new_alt_count = current_ref, 0, ref_count
                else:
                    new_alt, new_ref_count, new_alt_count = current_alt, 0, alt_count
                is_complex = True  # ← FLAG AS COMPLEX
                correction_type = 'complex_mismatch'
                n_complex += 1
            
            # Calculate BAF
            baf_alt = new_alt_count / total_dp if total_dp > 0 else 0
            
            # Infer genotype
            if total_dp == 0:
                gt_inferred = './.'
            elif new_alt_count == 0:
                gt_inferred = '0/0'
            elif new_ref_count == 0:
                gt_inferred = '1/1'
            else:
                alt_fraction = new_alt_count / total_dp
                if 0.2 <= alt_fraction <= 0.8:
                    gt_inferred = '0/1'
                elif alt_fraction < 0.2:
                    gt_inferred = '0/0'
                else:
                    gt_inferred = '1/1'
            
            # ============================================================
            # BUILD OUTPUT ROW WITH FLAGS
            # ============================================================
            corrected_rows.append({
                'CHROM': chrom,
                'POS': pos,
                'REF': new_ref,
                'ALT': new_alt,
                'DP': total_dp,
                'REFcounts': new_ref_count,
                'ALTcounts': new_alt_count,
                'BAF_ALT': baf_alt,
                'GT_inferred': gt_inferred,
                'is_complex': is_complex,           # ← NEW COLUMN
                'correction_type': correction_type   # ← NEW COLUMN (optional, for diagnostics)
            })
        
        fasta.close()
        
        # Create output dataframe
        corrected_df = pd.DataFrame(corrected_rows)
        
        # Save with compression
        corrected_df.to_csv(output_file, sep='\t', index=False, compression='gzip')
        
        logging.info(f"✓ {sample_id}: {len(corrected_df)} variants, {n_swapped} swapped, {n_complex} complex")
        
        return {
            'sample_id': sample_id,
            'status': 'success',
            'n_variants': len(corrected_df),
            'n_swapped': n_swapped,
            'n_complex': n_complex,
            'input_file': str(input_file),
            'output_file': str(output_file)
        }
        
    except Exception as e:
        logging.error(f"✗ {input_file.name}: {str(e)}")
        return {
            'sample_id': input_file.stem,
            'status': 'failed',
            'error': str(e),
            'input_file': str(input_file),
            'output_file': ''
        }


def process_batch(file_list, reference_fasta, output_dir, task_id):
    """
    Process a batch of files assigned to this task
    """
    results = []
    
    logging.info(f"Task {task_id}: Processing {len(file_list)} files")
    
    for input_file in tqdm(file_list, desc=f"Task {task_id}"):
        result = correct_single_file(input_file, reference_fasta, output_dir)
        results.append(result)
    
    return results


def main():
    """
    Main execution function for SLURM array job
    """
    
    # ========================================================================
    # CONFIGURATION
    # ========================================================================
    

    FOLDER1 = "/storage/scratch01/groups/co/pipeline_testing/alleleSpecific_downsampling/downsampled_JBLAB231"
    FOLDER2 = "/storage/scratch01/groups/co/pipeline_testing/alleleSpecific_downsampling/downsampled_JBLAB17041"

    REFERENCE_FASTA = "/home/_groups/co/reference/hg19/hg19.fa"

    OUTPUT_FOLDER1 = "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/corrected_downsampled_JBLAB231_flagProblematic"
    OUTPUT_FOLDER2 = "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/corrected_downsampled_JBLAB17041_flagProblematic"
    
    # ========================================================================
    
    # Setup paths
    folder1 = Path(FOLDER1)
    folder2 = Path(FOLDER2)
    output_folder1 = Path(OUTPUT_FOLDER1)
    output_folder2 = Path(OUTPUT_FOLDER2)
    
    output_folder1.mkdir(exist_ok=True, parents=True)
    output_folder2.mkdir(exist_ok=True, parents=True)
    
    output_base_dir = output_folder1.parent
    
    # Setup logging
    log_file = setup_logging(output_base_dir, SLURM_ARRAY_TASK_ID)
    
    logging.info("="*80)
    logging.info(f"REF/ALT CORRECTION - SLURM ARRAY TASK {SLURM_ARRAY_TASK_ID}")
    logging.info("="*80)
    logging.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Log file: {log_file}")
    
    # Find all files
    files_folder1 = sorted(list(folder1.glob("*.tumour.tso500_sparse_downsample.txt.gz")))
    files_folder2 = sorted(list(folder2.glob("*.tumour.tso500_sparse_downsample.txt.gz")))
    
    all_files = files_folder1 + files_folder2
    all_output_dirs = [output_folder1] * len(files_folder1) + [output_folder2] * len(files_folder2)
    
    total_files = len(all_files)
    
    logging.info(f"Total files found: {total_files}")
    logging.info(f"  Folder 1: {len(files_folder1)}")
    logging.info(f"  Folder 2: {len(files_folder2)}")
    
    # Divide files among tasks
    files_per_task = total_files // SLURM_ARRAY_TASK_COUNT
    remainder = total_files % SLURM_ARRAY_TASK_COUNT
    
    # Calculate start and end indices for this task
    start_idx = SLURM_ARRAY_TASK_ID * files_per_task + min(SLURM_ARRAY_TASK_ID, remainder)
    end_idx = start_idx + files_per_task + (1 if SLURM_ARRAY_TASK_ID < remainder else 0)
    
    # Get files for this task
    task_files = all_files[start_idx:end_idx]
    task_output_dirs = all_output_dirs[start_idx:end_idx]
    
    logging.info(f"Task {SLURM_ARRAY_TASK_ID}: Processing files {start_idx} to {end_idx-1}")
    logging.info(f"Task {SLURM_ARRAY_TASK_ID}: {len(task_files)} files assigned")
    
    # Process files
    results = []
    for input_file, output_dir in zip(task_files, task_output_dirs):
        result = correct_single_file(input_file, REFERENCE_FASTA, output_dir)
        results.append(result)
    
    # Save task results
    results_df = pd.DataFrame(results)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    task_summary_file = output_base_dir / "logs" / f"task_{SLURM_ARRAY_TASK_ID}_summary_{timestamp}.csv"
    results_df.to_csv(task_summary_file, index=False)
    
    # Summary
    n_success = (results_df['status'] == 'success').sum()
    n_failed = (results_df['status'] == 'failed').sum()
    n_skipped = (results_df['status'] == 'skipped').sum()
    
    logging.info("\n" + "="*80)
    logging.info(f"TASK {SLURM_ARRAY_TASK_ID} SUMMARY")
    logging.info("="*80)
    logging.info(f"Files processed: {len(results_df)}")
    logging.info(f"Success: {n_success}")
    logging.info(f"Failed: {n_failed}")
    logging.info(f"Skipped: {n_skipped}")
    
    if n_success > 0:
        success_results = results_df[results_df['status'] == 'success']
        total_variants = success_results['n_variants'].sum()
        total_swapped = success_results['n_swapped'].sum()
        
        logging.info(f"Total variants: {total_variants:,}")
        logging.info(f"Total swapped: {total_swapped:,}")
    
    logging.info(f"Task summary saved to: {task_summary_file}")
    logging.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if n_failed > 0:
        logging.warning(f"Task {SLURM_ARRAY_TASK_ID} completed with {n_failed} failures")
        sys.exit(1)
    else:
        logging.info(f"Task {SLURM_ARRAY_TASK_ID} completed successfully")


if __name__ == "__main__":
    main()

########

    

### # Or as background job (recommended for long runs)
### nohup python correct_REFALT_downsampled.py > correction_run.out 2>&1 &

# Check it's running
### ps aux | grep correct_ref_alt_batch.py

'''
# Watch log in real-time
tail -f /path/to/output/parent/logs/ref_correction_*.log

# Check number of processed files
ls -1 /path/to/output_folder1/*.corrected.txt.gz | wc -l
ls -1 /path/to/output_folder2/*.corrected.txt.gz | wc -l

# Check for errors in log
grep "ERROR" /path/to/output/parent/logs/ref_correction_*.log
grep "✗" /path/to/output/parent/logs/ref_correction_*.log
```

## What Gets Created

After running, you'll have this structure:
```
/your/output/parent/
├── folder1_corrected/
│   ├── sample1.tumour.tso500_sparse_downsample_corrected.txt.gz
│   ├── sample2.tumour.tso500_sparse_downsample_corrected.txt.gz
│   └── ...
├── folder2_corrected/
│   ├── sample3.tumour.tso500_sparse_downsample_corrected.txt.gz
│   └── ...
└── logs/
    ├── ref_correction_20260318_143022.log           # Main log file
    ├── correction_script_20260318_143022.py         # Copy of script
    ├── correction_summary_20260318_143022.csv       # Summary table
    └── failed_files_20260318_143022.csv            # Failed files (if any)
'''