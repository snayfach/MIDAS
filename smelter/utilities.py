#!/usr/bin/env python3

import os
from os.path import abspath
import sys
import time
import threading
import subprocess
import multiprocessing


# Thread-safe and timestamped prints.
tslock = multiprocessing.RLock()


def timestamp(t):
    # We do not use "{:.3f}".format(time.time()) because its result may be
    # up to 0.5 ms in the future due to rounding.  We want flooring here.
    s = str(int(t * 10))
    return s[:-1] + "." + s[-1:]


def tsfmt(msg):
    ts = timestamp(time.time()) + " "
    msg = ts + msg.replace("\n", "\n" + ts)
    return msg


def tsout(msg):
    with tslock:
        sys.stdout.write(str(msg))
        sys.stdout.write("\n")


def tserr(msg):
    with tslock:
        sys.stderr.write(str(msg))
        sys.stderr.write("\n")


def tsprint(msg):
    tserr(tsfmt(msg))


def check_output(command, quiet=False):
    assert sys.version_info >= (3, 6), "Please run this script with Python version >= 3.6."
    shell = isinstance(command, str)
    if not quiet:
        command_str = command if shell else " ".join(command)
        tsprint(repr(command_str))
    return subprocess.check_output(command, shell=shell).decode('utf-8')


def backtick(command):
    return check_output(command, quiet=True).strip()


def makedirs(newdir, exist_ok):
    try:
        os.makedirs(newdir, exist_ok=exist_ok)
    except Exception as e:
        if not exist_ok:
            e.help_text = f"If directory '{abspath(newdir)}' exists, please rename or remove it, then try again."
        raise


# A simple rows -> hashes converter.
# Credit: github.com/snayfatch/MIDAS
def parse_table(rows):
    headers = next(rows)  # pylint: disable=stop-iteration-return
    for values in rows:
        assert len(headers) == len(values)
        yield dict(zip(headers, values))


def tsv_rows(path):
    # TODO:  Support s3 and compressed files.
    with open(path, "r") as stream:
        for line in stream:
            yield line.rstrip("\n").split("\t")


# A simple progress tracker class.
# Credit:  github.com/chanzuckerberg/idseq-bench
class ProgressTracker:  # pylint: disable=too-few-public-methods

    def __init__(self, target):
        self.target = target
        self.current = 0
        self.t_start = time.time()
        self.t_last_print = 0

    def advance(self, amount):
        PESSIMISM = 1.25
        self.current += amount
        now = time.time()
        if amount == 0 or now - self.t_last_print > 30:  # update every 30 seconds or when forced with 0 amount
            self.t_last_print = now
            t_elapsed = now - self.t_start
            t_remaining = (t_elapsed / self.current) * self.target - t_elapsed
            t_remaining *= PESSIMISM
            t_eta = self.t_start + t_elapsed + t_remaining
            t_eta_str = time.strftime("%H:%M:%S", time.localtime(t_eta))
            tsprint(f"*** {self.current/self.target*100:3.1f} percent done, {t_elapsed/60:3.1f} minutes elapsed, {t_remaining/60:3.1f} minutes remaining, ETA {t_eta_str} ***")


if __name__ == "__main__":
    tsprint(f"Hello from {backtick('pwd')}.  Put tests here.")
