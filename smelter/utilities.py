#!/usr/bin/env python3

import os
import sys
import time
import threading
import subprocess


# Thread-safe and timestamped prints.
tslock = threading.RLock()


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
    # Assuming python >= 3.6
    shell = isinstance(command, str)
    if not quiet:
        command_str = command if shell else " ".join(command)
        print(repr(command_str))
    return subprocess.check_output(command, shell=shell).decode('utf-8')


def backtick(command):
    return check_output(command, quiet=True).strip()


def makedirs(newdir, exist_ok):
    try:
        os.makedirs(newdir, exist_ok=exist_ok)
    except Exception as e:
        if not exist_ok:
            e.help_text = f"If directory '{os.path.abspath(newdir)}' exists, please rename or remove it, then try again."
        raise


if __name__ == "__main__":
    tsprint(f"Hello from {backtick('pwd')}.  Put tests here.")
