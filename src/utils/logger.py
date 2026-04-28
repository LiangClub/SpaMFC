"""
SpaMFC 日志模块

提供统一的日志功能，支持：
1. 控制台和文件双输出
2. 时间戳日志文件
3. 进度跟踪
4. 错误追踪
"""

import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional

LOG_LEVELS = {
    'DEBUG': '调试',
    'INFO': '信息',
    'WARNING': '警告',
    'ERROR': '错误',
    'CRITICAL': '严重'
}


class ChineseFormatter(logging.Formatter):
    def format(self, record):
        record.levelname = LOG_LEVELS.get(record.levelname, record.levelname)
        return super().format(record)


def setup_logger(
    name: str = "SpaMFC",
    log_dir: Optional[str] = None,
    level: int = logging.INFO,
    verbose: bool = True,
    save_log: bool = True
) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG if verbose else level)
    logger.handlers.clear()
    
    console_format = ChineseFormatter(
        '%(asctime)s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)
    
    if save_log and log_dir:
        log_path = Path(log_dir)
        log_path.mkdir(parents=True, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = log_path / f"spamfc_{timestamp}.log"
        
        file_format = ChineseFormatter(
            '%(asctime)s | %(levelname)-8s | %(name)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)
        
        logger.debug(f"日志文件: {log_file}")
    
    return logger


def get_logger(name: str = "SpaMFC") -> logging.Logger:
    logger = logging.getLogger(name)
    
    if not logger.handlers:
        return setup_logger(name, save_log=False)
    
    return logger


class ProgressLogger:
    def __init__(self, logger: logging.Logger, total: int, desc: str = "处理中"):
        self.logger = logger
        self.total = total
        self.desc = desc
        self.current = 0
        self.start_time = datetime.now()
    
    def update(self, n: int = 1):
        self.current += n
        elapsed = (datetime.now() - self.start_time).total_seconds()
        if self.total > 0:
            pct = self.current / self.total * 100
            self.logger.debug(f"{self.desc}: {self.current}/{self.total} ({pct:.1f}%) - {elapsed:.1f}秒")
    
    def close(self):
        elapsed = (datetime.now() - self.start_time).total_seconds()
        self.logger.info(f"{self.desc}完成: {self.current}项, 耗时{elapsed:.2f}秒")


def log_step(logger: logging.Logger, step: int, total_steps: int, message: str):
    logger.info(f"[步骤 {step}/{total_steps}] {message}")


def log_section(logger: logging.Logger, title: str, char: str = "=", width: int = 60):
    logger.info("")
    logger.info(char * width)
    logger.info(title)
    logger.info(char * width)


def log_subsection(logger: logging.Logger, title: str, char: str = "-", width: int = 50):
    logger.info("")
    logger.info(char * width)
    logger.info(title)
    logger.info(char * width)


def log_error(logger: logging.Logger, error: Exception, context: str = ""):
    import traceback
    logger.error(f"错误 [{context}]: {str(error)}")
    logger.debug(traceback.format_exc())