import sys

from loguru import logger
from typing_extensions import Final

_format: Final = "<green>{time}</green> | {function} | <level>{message}</level>"

logger.remove()
logger.add(sys.stderr, colorize=True, format=_format)
logger.add("ensep.log", colorize=True, format=_format)