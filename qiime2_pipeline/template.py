import subprocess
from typing import Any
from datetime import datetime


class Settings:

    workdir: str
    outdir: str
    threads: int
    debug: bool
    mock: bool
    for_publication: bool

    def __init__(
            self,
            workdir: str,
            outdir: str,
            threads: int,
            debug: bool,
            mock: bool,
            for_publication: bool):

        self.workdir = workdir
        self.outdir = outdir
        self.threads = threads
        self.debug = debug
        self.mock = mock
        self.for_publication = for_publication


class Logger:

    WARNING = 'WARNING'
    INFO = 'INFO'
    DEBUG = 'DEBUG'

    name: str
    level: str

    def __init__(self, name: str, level: str):
        self.name = name
        assert level in [self.INFO, self.DEBUG]
        self.level = level

    def warning(self, msg: Any):
        self.__print(level=self.WARNING, msg=msg)

    def info(self, msg: Any):
        if self.level in [self.WARNING]:
            return
        self.__print(level=self.INFO, msg=msg)

    def debug(self, msg: Any):
        if self.level in [self.INFO, self.WARNING]:
            return
        self.__print(level=self.DEBUG, msg=msg)

    def __print(self, level: str, msg: Any):
        print(f'{self.name}\t{level}\t{datetime.now()}', flush=True)
        print(f'{msg}\n', flush=True)


class Processor:

    CMD_LINEBREAK = ' \\\n  '
    MAX_TRY = 3

    settings: Settings
    workdir: str
    outdir: str
    threads: int
    debug: bool
    mock: bool

    logger: Logger

    def __init__(self, settings: Settings):

        self.settings = settings
        self.workdir = settings.workdir
        self.outdir = settings.outdir
        self.threads = settings.threads
        self.debug = settings.debug
        self.mock = self.settings.mock

        self.logger = Logger(
            name=self.__class__.__name__,
            level=Logger.DEBUG if self.debug else Logger.INFO
        )

    def call(self, cmd: str):
        self.logger.info(cmd)
        if self.mock:
            return

        tried = 0
        while True:
            try:
                subprocess.check_call(cmd, shell=True)
                break  # succeed and break from the loop
            except Exception as e:
                self.logger.info(f'Failed: {e}')
                tried += 1

            if tried >= self.MAX_TRY:
                raise Exception('Failed too many times')
