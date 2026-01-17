"""
Rate Limiter

This module provides rate limiting functionality to prevent API rate limit violations.
"""

import time
from typing import Optional
from loguru import logger


class RateLimiter:
    """
    Simple rate limiter with configurable delay between calls.
    
    Ensures a minimum time delay between consecutive operations
    to comply with API rate limits.
    """
    
    def __init__(self, delay: float = 0.2):
        """
        Initialize rate limiter.
        
        Parameters:
        -----------
        delay : float
            Minimum delay between calls in seconds
        """
        self.delay = delay
        self.last_call_time: Optional[float] = None
    
    def wait(self) -> None:
        """
        Wait if necessary to maintain rate limit.
        
        Calculates time since last call and sleeps if needed
        to maintain the minimum delay.
        """
        if self.last_call_time is not None:
            elapsed = time.time() - self.last_call_time
            if elapsed < self.delay:
                sleep_time = self.delay - elapsed
                logger.debug(f"Rate limiting: sleeping for {sleep_time:.3f}s")
                time.sleep(sleep_time)
        
        self.last_call_time = time.time()
    
    def reset(self) -> None:
        """Reset rate limiter state."""
        self.last_call_time = None
    
    def set_delay(self, delay: float) -> None:
        """
        Update rate limit delay.
        
        Parameters:
        -----------
        delay : float
            New delay in seconds
        """
        self.delay = delay
        logger.debug(f"Rate limit delay updated to {delay}s")


class AdaptiveRateLimiter(RateLimiter):
    """
    Rate limiter with adaptive delay based on error rates.
    
    Automatically increases delay when errors are detected
    and gradually decreases it during stable operation.
    """
    
    def __init__(
        self,
        initial_delay: float = 0.2,
        min_delay: float = 0.1,
        max_delay: float = 2.0,
        increase_factor: float = 2.0,
        decrease_factor: float = 0.9
    ):
        """
        Initialize adaptive rate limiter.
        
        Parameters:
        -----------
        initial_delay : float
            Initial delay in seconds
        min_delay : float
            Minimum delay in seconds
        max_delay : float
            Maximum delay in seconds
        increase_factor : float
            Factor to multiply delay by on errors
        decrease_factor : float
            Factor to multiply delay by on success
        """
        super().__init__(delay=initial_delay)
        self.min_delay = min_delay
        self.max_delay = max_delay
        self.increase_factor = increase_factor
        self.decrease_factor = decrease_factor
        self.consecutive_successes = 0
        self.success_threshold = 10  # Decrease delay after this many successes
    
    def on_error(self) -> None:
        """
        Handle error event by increasing delay.
        
        Resets success counter and increases delay to reduce load.
        """
        old_delay = self.delay
        self.delay = min(self.delay * self.increase_factor, self.max_delay)
        self.consecutive_successes = 0
        
        if self.delay != old_delay:
            logger.info(f"Rate limit increased: {old_delay:.3f}s -> {self.delay:.3f}s")
    
    def on_success(self) -> None:
        """
        Handle success event by potentially decreasing delay.
        
        Gradually decreases delay after sustained success to improve throughput.
        """
        self.consecutive_successes += 1
        
        if self.consecutive_successes >= self.success_threshold:
            old_delay = self.delay
            self.delay = max(self.delay * self.decrease_factor, self.min_delay)
            self.consecutive_successes = 0
            
            if self.delay != old_delay:
                logger.debug(f"Rate limit decreased: {old_delay:.3f}s -> {self.delay:.3f}s")
    
    def reset(self) -> None:
        """Reset rate limiter to initial state."""
        super().reset()
        self.consecutive_successes = 0
