#!/usr/bin/env python3
"""
Performance Analysis Script with Automated Log Parsing (Multi-Run Support)

This script analyzes the performance of the carbohydrate analysis module
by parsing actual log files and extracting timing and classification metrics.
Supports analysis of multiple file processing runs within a single log file.
"""

import re
from datetime import datetime
from typing import Dict, List, Optional
from pathlib import Path


class LogAnalyzer:
    """
    Analyzes performance metrics from log files.
    
    This class parses log files to extract timing information for each stage
    of the processing pipeline and provides performance analysis. Supports
    multi-run analysis where multiple files are processed in a single log.
    """
    
    def __init__(self, log_file_path: str):
        """
        Initialize the LogAnalyzer.
        
        Parameters:
        -----------
        log_file_path : str
            Path to the log file to analyze
        """
        self.log_file_path = Path(log_file_path)
        self.log_lines = []
        self.runs = []  # List of runs (one per file processed)
        
        # Compile regex patterns for log parsing
        self._compile_patterns()
    
    def _compile_patterns(self):
        """Compile all regex patterns for log parsing."""
        # Timestamp pattern
        self.timestamp_pattern = re.compile(r'^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})')
        
        # File processing boundary marker
        self.file_processing_pattern = re.compile(r'Processing file (\d+)/(\d+): (.+)')
        
        # Stage 1: InChIKey to CID Resolution
        self.stage1_resolving_pattern = re.compile(r'Resolving (inchikey|smiles): (.+?)\.\.\.')
        self.stage1_resolved_pattern = re.compile(r'Resolved (.+?)\.\.\. -> CID (\d+)')
        self.stage1_no_cid_pattern = re.compile(r'No CID found for (.+?)\.\.\.')
        
        # Stage 2: Batch Property Fetching
        self.stage2_start_pattern = re.compile(r'Fetching properties for (\d+) CIDs in chunks of (\d+)')
        self.stage2_chunk_pattern = re.compile(r'Processing chunk (\d+)/(\d+) \((\d+) CIDs\)')
        self.stage2_complete_pattern = re.compile(r'Fetched properties for (\d+)/(\d+) CIDs')
        
        # Stage 3: Classification
        self.stage3_start_pattern = re.compile(r'Starting batch processing of (\d+) compounds')
        self.stage3_fast_path_pattern = re.compile(r'Classifying via ChEBI ancestry')
        self.stage3_fallback_pattern = re.compile(r'falling back to PubChem classification')
        self.stage3_complete_pattern = re.compile(r'Batch processing complete: (\d+)/(\d+) compounds processed successfully')
    
    def load_log_file(self):
        """Load the log file into memory and split into runs."""
        print(f"Loading log file: {self.log_file_path}")
        with open(self.log_file_path, 'r', encoding='utf-8') as f:
            self.log_lines = f.readlines()
        print(f"Loaded {len(self.log_lines)} log lines")
        
        # Split log into runs based on file processing boundaries
        self._split_into_runs()
    
    def _split_into_runs(self):
        """Split log lines into separate runs (one per file processed)."""
        self.runs = []
        current_run = {
            'file_num': None,
            'total_files': None,
            'filename': None,
            'start_line': 0,
            'end_line': 0,
            'lines': []
        }
        
        for i, line in enumerate(self.log_lines):
            match = self.file_processing_pattern.search(line)
            if match:
                # Save previous run if it has content
                if current_run['lines']:
                    current_run['end_line'] = i - 1
                    self.runs.append(current_run)
                
                # Start new run
                current_run = {
                    'file_num': int(match.group(1)),
                    'total_files': int(match.group(2)),
                    'filename': match.group(3),
                    'start_line': i,
                    'end_line': 0,
                    'lines': []
                }
            
            current_run['lines'].append(line)
        
        # Don't forget the last run
        if current_run['lines']:
            current_run['end_line'] = len(self.log_lines) - 1
            self.runs.append(current_run)
        
        # Filter out empty runs (runs without filenames - typically metadata sections)
        self.runs = [run for run in self.runs if run['filename'] is not None]
        
        print(f"Found {len(self.runs)} processing run(s)")
    
    def extract_timestamp(self, line: str) -> Optional[datetime]:
        """Extract timestamp from a log line."""
        match = self.timestamp_pattern.match(line)
        if match:
            return datetime.strptime(match.group(1), '%Y-%m-%d %H:%M:%S')
        return None
    
    def analyze_stage1(self, lines: List[str]) -> Dict:
        """Analyze Stage 1: InChIKey to CID Resolution."""
        resolving_events = []
        resolved_events = []
        failed_events = []
        
        for line in lines:
            timestamp = self.extract_timestamp(line)
            if not timestamp:
                continue
            
            if self.stage1_resolving_pattern.search(line):
                resolving_events.append(timestamp)
            if self.stage1_resolved_pattern.search(line):
                resolved_events.append(timestamp)
            if self.stage1_no_cid_pattern.search(line):
                failed_events.append(timestamp)
        
        total_identifiers = len(resolving_events)
        resolved_cids = len(resolved_events)
        failed_ids = len(failed_events)
        
        if resolving_events and resolved_events:
            start_time = min(resolving_events)
            end_time = max(resolved_events + failed_events) if failed_events else max(resolved_events)
            duration = (end_time - start_time).total_seconds()
            rate = duration / total_identifiers if total_identifiers > 0 else 0
        else:
            start_time = None
            end_time = None
            duration = 0
            rate = 0
        
        return {
            'total_identifiers': total_identifiers,
            'resolved_cids': resolved_cids,
            'failed_ids': failed_ids,
            'success_rate': (resolved_cids / total_identifiers * 100) if total_identifiers > 0 else 0,
            'start_time': start_time,
            'end_time': end_time,
            'duration_seconds': duration,
            'seconds_per_identifier': rate
        }
    
    def analyze_stage2(self, lines: List[str]) -> Dict:
        """Analyze Stage 2: Batch Property Fetching."""
        start_event = None
        end_event = None
        total_cids = 0
        fetched_cids = 0
        chunk_events = []
        
        for line in lines:
            timestamp = self.extract_timestamp(line)
            if not timestamp:
                continue
            
            match_start = self.stage2_start_pattern.search(line)
            if match_start:
                start_event = timestamp
                total_cids = int(match_start.group(1))
                continue
            
            if self.stage2_chunk_pattern.search(line):
                chunk_events.append(timestamp)
                continue
            
            match_complete = self.stage2_complete_pattern.search(line)
            if match_complete:
                end_event = timestamp
                fetched_cids = int(match_complete.group(1))
                continue
        
        if start_event and end_event:
            duration = (end_event - start_event).total_seconds()
            rate = duration / total_cids if total_cids > 0 else 0
        else:
            duration = 0
            rate = 0
        
        return {
            'total_cids': total_cids,
            'fetched_cids': fetched_cids,
            'chunks': len(chunk_events),
            'start_time': start_event,
            'end_time': end_event,
            'duration_seconds': duration,
            'seconds_per_cid': rate
        }
    
    def analyze_stage3(self, lines: List[str]) -> Dict:
        """Analyze Stage 3: Classification."""
        start_event = None
        end_event = None
        total_compounds = 0
        processed_compounds = 0
        fast_path_count = 0
        fallback_count = 0
        
        for line in lines:
            timestamp = self.extract_timestamp(line)
            if not timestamp:
                continue
            
            match_start = self.stage3_start_pattern.search(line)
            if match_start:
                start_event = timestamp
                total_compounds = int(match_start.group(1))
                continue
            
            if self.stage3_fast_path_pattern.search(line):
                fast_path_count += 1
                continue
            
            if self.stage3_fallback_pattern.search(line):
                fallback_count += 1
                continue
            
            match_complete = self.stage3_complete_pattern.search(line)
            if match_complete:
                end_event = timestamp
                processed_compounds = int(match_complete.group(1))
                continue
        
        if start_event and end_event:
            duration = (end_event - start_event).total_seconds()
            rate = duration / processed_compounds if processed_compounds > 0 else 0
        else:
            duration = 0
            rate = 0
        
        total_classified = fast_path_count + fallback_count
        fast_path_pct = (fast_path_count / total_classified * 100) if total_classified > 0 else 0
        fallback_pct = (fallback_count / total_classified * 100) if total_classified > 0 else 0
        
        return {
            'total_compounds': total_compounds,
            'processed_compounds': processed_compounds,
            'fast_path_count': fast_path_count,
            'fallback_count': fallback_count,
            'total_classified': total_classified,
            'fast_path_percentage': fast_path_pct,
            'fallback_percentage': fallback_pct,
            'start_time': start_event,
            'end_time': end_event,
            'duration_seconds': duration,
            'seconds_per_compound': rate
        }
    
    def analyze_run(self, run: Dict) -> Dict:
        """Analyze a single processing run."""
        stage1 = self.analyze_stage1(run['lines'])
        stage2 = self.analyze_stage2(run['lines'])
        stage3 = self.analyze_stage3(run['lines'])
        
        total_duration = (
            stage1.get('duration_seconds', 0) +
            stage2.get('duration_seconds', 0) +
            stage3.get('duration_seconds', 0)
        )
        
        return {
            'file_num': run['file_num'],
            'filename': run['filename'],
            'stage1': stage1,
            'stage2': stage2,
            'stage3': stage3,
            'total_duration': total_duration
        }
    
    def analyze_all(self) -> Dict:
        """Run complete analysis on all runs."""
        self.load_log_file()
        
        if not self.runs:
            print("No runs found in log file")
            return {}
        
        print(f"\nAnalyzing {len(self.runs)} run(s)...")
        
        # Analyze each run separately
        run_results = []
        for i, run in enumerate(self.runs):
            print(f"\n--- Run {i+1}/{len(self.runs)}: {run['filename']} ---")
            result = self.analyze_run(run)
            run_results.append(result)
            
            # Print summary for this run
            s1 = result['stage1']
            s2 = result['stage2']
            s3 = result['stage3']
            print(f"  Stage 1: {s1['total_identifiers']} InChIKeys, {s1['resolved_cids']} resolved ({s1['success_rate']:.1f}%), {s1['duration_seconds']:.1f}s, {s1['seconds_per_identifier']:.3f}s/id")
            print(f"  Stage 2: {s2['total_cids']} CIDs, {s2['duration_seconds']:.1f}s, {s2['seconds_per_cid']:.4f}s/CID")
            print(f"  Stage 3: {s3['processed_compounds']} compounds, {s3['fast_path_count']} fast ({s3['fast_path_percentage']:.1f}%), {s3['duration_seconds']:.1f}s, {s3['seconds_per_compound']:.3f}s/compound")
            print(f"  Total: {result['total_duration']:.1f}s ({result['total_duration']/60:.1f} min)")
        
        # Calculate aggregate statistics
        aggregate = self._calculate_aggregate_stats(run_results)
        
        return {
            'runs': run_results,
            'aggregate': aggregate,
            'num_runs': len(self.runs)
        }
    
    def _calculate_aggregate_stats(self, run_results: List[Dict]) -> Dict:
        """Calculate aggregate statistics across all runs."""
        if not run_results:
            return {}
        
        # Collect metrics from all runs
        total_identifiers = sum(r['stage1']['total_identifiers'] for r in run_results)
        total_resolved = sum(r['stage1']['resolved_cids'] for r in run_results)
        total_failed = sum(r['stage1']['failed_ids'] for r in run_results)
        
        total_cids_stage2 = sum(r['stage2']['total_cids'] for r in run_results)
        total_chunks = sum(r['stage2']['chunks'] for r in run_results)
        
        total_compounds = sum(r['stage3']['processed_compounds'] for r in run_results)
        total_fast_path = sum(r['stage3']['fast_path_count'] for r in run_results)
        total_fallback = sum(r['stage3']['fallback_count'] for r in run_results)
        
        total_duration_s1 = sum(r['stage1']['duration_seconds'] for r in run_results)
        total_duration_s2 = sum(r['stage2']['duration_seconds'] for r in run_results)
        total_duration_s3 = sum(r['stage3']['duration_seconds'] for r in run_results)
        total_duration = total_duration_s1 + total_duration_s2 + total_duration_s3
        
        # Calculate rates
        avg_rate_s1 = total_duration_s1 / total_identifiers if total_identifiers > 0 else 0
        avg_rate_s2 = total_duration_s2 / total_cids_stage2 if total_cids_stage2 > 0 else 0
        avg_rate_s3 = total_duration_s3 / total_compounds if total_compounds > 0 else 0
        
        success_rate = (total_resolved / total_identifiers * 100) if total_identifiers > 0 else 0
        fast_path_pct = (total_fast_path / (total_fast_path + total_fallback) * 100) if (total_fast_path + total_fallback) > 0 else 0
        
        return {
            'stage1': {
                'total_identifiers': total_identifiers,
                'resolved_cids': total_resolved,
                'failed_ids': total_failed,
                'success_rate': success_rate,
                'duration_seconds': total_duration_s1,
                'seconds_per_identifier': avg_rate_s1
            },
            'stage2': {
                'total_cids': total_cids_stage2,
                'chunks': total_chunks,
                'duration_seconds': total_duration_s2,
                'seconds_per_cid': avg_rate_s2
            },
            'stage3': {
                'processed_compounds': total_compounds,
                'fast_path_count': total_fast_path,
                'fallback_count': total_fallback,
                'fast_path_percentage': fast_path_pct,
                'fallback_percentage': 100 - fast_path_pct,
                'duration_seconds': total_duration_s3,
                'seconds_per_compound': avg_rate_s3
            },
            'total_duration': total_duration
        }
    
    def generate_formulas(self, results: Dict) -> Dict:
        """Generate runtime estimation formulas."""
        if not results or 'aggregate' not in results:
            return {}
        
        aggregate = results['aggregate']
        stage1 = aggregate['stage1']
        stage2 = aggregate['stage2']
        stage3 = aggregate['stage3']
        
        t1_rate = stage1.get('seconds_per_identifier', 0)
        t2_rate = stage2.get('seconds_per_cid', 0)
        t3_rate = stage3.get('seconds_per_compound', 0)
        success_rate = stage1.get('success_rate', 0) / 100.0
        
        # Calculate total rate per input InChIKey (accounting for success rate)
        total_rate = t1_rate + (t2_rate * success_rate) + (t3_rate * success_rate)
        
        return {
            'stage1_rate': t1_rate,
            'stage2_rate': t2_rate,
            'stage3_rate': t3_rate,
            'success_rate': success_rate,
            'total_rate': total_rate,
            'formula_stage1': f"T1(N) = {t1_rate:.3f} x N seconds",
            'formula_stage2': f"T2(N) = {t2_rate:.4f} x N seconds",
            'formula_stage3': f"T3(N) = {t3_rate:.3f} x N seconds",
            'formula_total': f"T_total(N) = {total_rate:.3f} x N seconds"
        }
    
    def generate_report(self, results: Dict, formulas: Dict):
        """Generate comprehensive performance report."""
        print("\n" + "=" * 80)
        print("PERFORMANCE ANALYSIS REPORT")
        print("=" * 80)
        
        if not results:
            print("\nNo results to report")
            return
        
        aggregate = results.get('aggregate', {})
        if not aggregate:
            print("\nNo aggregate data available")
            return
        
        # Overall summary
        print(f"\nANALYSIS SUMMARY")
        print("-" * 80)
        print(f"Number of runs: {results['num_runs']}")
        print(f"Total duration: {aggregate['total_duration']:.1f}s ({aggregate['total_duration']/60:.1f} minutes, {aggregate['total_duration']/3600:.1f} hours)")
        
        # Stage breakdown
        stage1 = aggregate['stage1']
        stage2 = aggregate['stage2']
        stage3 = aggregate['stage3']
        
        stage1_pct = (stage1['duration_seconds'] / aggregate['total_duration'] * 100) if aggregate['total_duration'] > 0 else 0
        stage2_pct = (stage2['duration_seconds'] / aggregate['total_duration'] * 100) if aggregate['total_duration'] > 0 else 0
        stage3_pct = (stage3['duration_seconds'] / aggregate['total_duration'] * 100) if aggregate['total_duration'] > 0 else 0
        
        print(f"\nStage 1 (CID Resolution): {stage1['duration_seconds']:.1f}s ({stage1_pct:.1f}%)")
        print(f"Stage 2 (Property Fetch): {stage2['duration_seconds']:.1f}s ({stage2_pct:.1f}%)")
        print(f"Stage 3 (Classification): {stage3['duration_seconds']:.1f}s ({stage3_pct:.1f}%)")
        
        # Detailed stage information
        print("\n" + "-" * 80)
        print("AGGREGATE STAGE 1: InChIKey to CID Resolution")
        print("-" * 80)
        print(f"Total identifiers (all runs): {stage1['total_identifiers']}")
        print(f"Resolved CIDs: {stage1['resolved_cids']}")
        print(f"Failed: {stage1['failed_ids']}")
        print(f"Success rate: {stage1['success_rate']:.1f}%")
        print(f"Average time per identifier: {stage1['seconds_per_identifier']:.3f}s")
        print(f"\n{formulas['formula_stage1']}")
        
        print("\n" + "-" * 80)
        print("AGGREGATE STAGE 2: Batch Property Fetching")
        print("-" * 80)
        print(f"Total CIDs (all runs): {stage2['total_cids']}")
        print(f"Chunks processed: {stage2['chunks']}")
        print(f"Average time per CID: {stage2['seconds_per_cid']:.4f}s")
        print(f"\n{formulas['formula_stage2']}")
        print("Note: Property fetching is highly efficient due to batch API")
        
        print("\n" + "-" * 80)
        print("AGGREGATE STAGE 3: Classification")
        print("-" * 80)
        print(f"Total compounds (all runs): {stage3['processed_compounds']}")
        print(f"\nPath usage (across all runs):")
        print(f"  Fast path (ChEBI ancestry): {stage3['fast_path_count']} ({stage3['fast_path_percentage']:.1f}%)")
        print(f"  Fallback (PubChem): {stage3['fallback_count']} ({stage3['fallback_percentage']:.1f}%)")
        print(f"\nAverage time per compound: {stage3['seconds_per_compound']:.3f}s")
        print(f"\n{formulas['formula_stage3']}")
        
        # Per-run summary
        print("\n" + "-" * 80)
        print("PER-RUN BREAKDOWN")
        print("-" * 80)
        for i, run in enumerate(results['runs']):
            print(f"\nRun {i+1}: {run['filename']}")
            print(f"  InChIKeys: {run['stage1']['total_identifiers']} -> {run['stage1']['resolved_cids']} CIDs ({run['stage1']['success_rate']:.1f}%)")
            print(f"  Stage 1: {run['stage1']['duration_seconds']:.1f}s ({run['stage1']['seconds_per_identifier']:.3f}s/id)")
            print(f"  Stage 2: {run['stage2']['duration_seconds']:.1f}s ({run['stage2']['seconds_per_cid']:.4f}s/CID)")
            print(f"  Stage 3: {run['stage3']['duration_seconds']:.1f}s ({run['stage3']['seconds_per_compound']:.3f}s/compound)")
            print(f"           Fast: {run['stage3']['fast_path_count']} ({run['stage3']['fast_path_percentage']:.1f}%), Fallback: {run['stage3']['fallback_count']} ({run['stage3']['fallback_percentage']:.1f}%)")
            print(f"  Total: {run['total_duration']:.1f}s ({run['total_duration']/60:.1f} min)")
        
        # Total runtime formula
        print("\n" + "-" * 80)
        print("RUNTIME ESTIMATION FORMULA")
        print("-" * 80)
        print(f"{formulas['formula_total']}")
        print(f"or approximately: T_total(N) = {formulas['total_rate']/60:.3f} x N minutes")
        
        # Example calculations
        print("\nExamples:")
        total_rate = formulas['total_rate']
        for n in [100, 500, 1000, 2000]:
            time_seconds = total_rate * n
            time_minutes = time_seconds / 60
            time_hours = time_minutes / 60
            if time_hours >= 1:
                print(f"  N = {n:4d}: ~{time_hours:.1f} hours")
            else:
                print(f"  N = {n:4d}: ~{time_minutes:.1f} minutes")
        
        print("\n" + "=" * 80)


def main():
    """Main function to run the analysis."""
    import sys
    
    # Get log file path from command line or use default
    if len(sys.argv) > 1:
        log_file = sys.argv[1]
    else:
        log_file = 'logs/20260117-2.log'
    
    # Initialize analyzer
    analyzer = LogAnalyzer(log_file)
    
    # Run analysis
    results = analyzer.analyze_all()
    
    # Generate formulas
    formulas = analyzer.generate_formulas(results)
    
    # Generate report
    analyzer.generate_report(results, formulas)


if __name__ == '__main__':
    main()
