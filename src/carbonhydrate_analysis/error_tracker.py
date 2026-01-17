"""
Error Tracker

This module provides error tracking and adaptive behavior for handling
server errors and failures.
"""

import time
from typing import List, Dict, Any
from collections import deque
from pathlib import Path
import json
from loguru import logger
from . import config


class ServerErrorTracker:
    """
    Tracks server errors and implements adaptive sleep behavior.
    
    Monitors error frequency and triggers longer delays when
    error rate exceeds threshold to reduce server load.
    """
    
    def __init__(
        self,
        window_size: int = 100,
        error_window: float = 60.0,
        error_threshold: int = 10,
        sleep_time: float = 180.0
    ):
        """
        Initialize error tracker.
        
        Parameters:
        -----------
        window_size : int
            Maximum number of errors to track
        error_window : float
            Time window in seconds to consider for error rate
        error_threshold : int
            Number of errors in window to trigger sleep
        sleep_time : float
            Time to sleep when threshold exceeded
        """
        self.error_times: deque = deque(maxlen=window_size)
        self.error_window = error_window
        self.error_threshold = error_threshold
        self.sleep_time = sleep_time
    
    def record_error(self) -> bool:
        """
        Record a server error and check if sleep is needed.
        
        Returns:
        --------
        bool
            True if sleep threshold was exceeded and sleep occurred
        """
        current_time = time.time()
        self.error_times.append(current_time)
        
        # Check if we have too many errors in the recent window
        recent_errors = sum(
            1 for t in self.error_times 
            if current_time - t <= self.error_window
        )
        
        if recent_errors >= self.error_threshold:
            logger.warning(
                f"High error rate detected: {recent_errors} errors in {self.error_window}s. "
                f"Sleeping for {self.sleep_time}s to reduce server load..."
            )
            time.sleep(self.sleep_time)
            self.reset()  # Clear error history after sleep
            return True
        
        return False
    
    def reset(self) -> None:
        """Clear error history."""
        self.error_times.clear()
    
    def get_recent_error_count(self) -> int:
        """
        Get count of errors in recent window.
        
        Returns:
        --------
        int
            Number of errors in the configured time window
        """
        current_time = time.time()
        return sum(
            1 for t in self.error_times 
            if current_time - t <= self.error_window
        )


class FailedIdentifierTracker:
    """
    Tracks failed identifiers and CIDs for later retry.
    
    Maintains separate lists for different identifier types
    and provides persistence to disk.
    """
    
    def __init__(self, cache_dir: Path):
        """
        Initialize failed identifier tracker.
        
        Parameters:
        -----------
        cache_dir : Path
            Directory to store failed identifier files
        """
        self.cache_dir = cache_dir
        self.failed_identifiers: Dict[str, List[str]] = {}
        self.failed_cids: List[int] = []
        
        # Ensure cache directory exists
        self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    def add_failed_identifier(self, identifier_type: str, identifier: str) -> None:
        """
        Record a failed identifier.
        
        Parameters:
        -----------
        identifier_type : str
            Type of identifier (e.g., 'inchikey', 'smiles')
        identifier : str
            The identifier that failed
        """
        if identifier_type not in self.failed_identifiers:
            self.failed_identifiers[identifier_type] = []
        
        if identifier not in self.failed_identifiers[identifier_type]:
            self.failed_identifiers[identifier_type].append(identifier)
            logger.debug(f"Recorded failed {identifier_type}: {identifier[:30]}...")
    
    def add_failed_cid(self, cid: int) -> None:
        """
        Record a failed CID.
        
        Parameters:
        -----------
        cid : int
            The CID that failed
        """
        if cid not in self.failed_cids:
            self.failed_cids.append(cid)
            logger.debug(f"Recorded failed CID: {cid}")
    
    def save_failed_identifiers(self, identifier_type: str) -> None:
        """
        Save failed identifiers to disk.
        
        Parameters:
        -----------
        identifier_type : str
            Type of identifier to save
        """
        if identifier_type not in self.failed_identifiers:
            return
        
        if not self.failed_identifiers[identifier_type]:
            return
        
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        filename = f'failed_{identifier_type}_{timestamp}.json'
        filepath = self.cache_dir / filename
        
        try:
            data = {
                'timestamp': timestamp,
                'identifier_type': identifier_type,
                'count': len(self.failed_identifiers[identifier_type]),
                'identifiers': self.failed_identifiers[identifier_type]
            }
            
            with open(filepath, 'w') as f:
                json.dump(data, f, indent=2)
            
            logger.info(
                f"Saved {len(self.failed_identifiers[identifier_type])} "
                f"failed {identifier_type}s to {filename}"
            )
        except Exception as e:
            logger.error(f"Failed to save failed identifiers: {e}")
    
    def save_failed_cids(self) -> None:
        """Save failed CIDs to disk."""
        if not self.failed_cids:
            return
        
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        filename = f'failed_cids_{timestamp}.json'
        filepath = self.cache_dir / filename
        
        try:
            data = {
                'timestamp': timestamp,
                'count': len(self.failed_cids),
                'cids': self.failed_cids
            }
            
            with open(filepath, 'w') as f:
                json.dump(data, f, indent=2)
            
            logger.info(f"Saved {len(self.failed_cids)} failed CIDs to {filename}")
        except Exception as e:
            logger.error(f"Failed to save failed CIDs: {e}")
    
    def get_failed_count(self, identifier_type: str = None) -> int:
        """
        Get count of failed identifiers.
        
        Parameters:
        -----------
        identifier_type : str, optional
            Specific type to count, or None for total
        
        Returns:
        --------
        int
            Count of failed identifiers
        """
        if identifier_type:
            return len(self.failed_identifiers.get(identifier_type, []))
        else:
            return sum(len(ids) for ids in self.failed_identifiers.values())
    
    def clear(self) -> None:
        """Clear all tracked failures."""
        self.failed_identifiers.clear()
        self.failed_cids.clear()
