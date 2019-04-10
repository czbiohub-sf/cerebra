import os
import subprocess


# Set input language USA unicode encoding setting
# Necessary because click assumes ascii input unless otherwise specified
# https://click.palletsprojects.com/en/7.x/python3/
unicode_usa = 'en_US.utf-8'
os.environ['LC_LANG'] = unicode_usa
os.environ["LC_ALL"] = unicode_usa

HOME = os.path.expanduser("~")


def sanitize_path(path):
    """Make sure path is absolute and user-expanded"""
    return os.path.abspath(os.path.expanduser(path))


def maybe_add_slash(path):
    """Add a final trailing slash if it wasn't there already"""
    with_trailing_slash = path if path.endswith('/') else path + '/'
    return with_trailing_slash


def decode(stderr_or_stdout):
    decoded = stderr_or_stdout.decode("utf-8").splitlines()
    lines = [x.strip() for x in decoded]
    return lines


def get_stdout_from_command(command):
    """Run a program on the command line, and save the stdout

    Parameters
    ----------
    command : list
        list of strings to submit to subprocess.run

    Returns
    -------
    lines : list
        Newline-separated strings from output of command
    """
    result = subprocess.run(command, stdout=subprocess.PIPE)
    lines = decode(result.stdout)
    return lines


def get_stdout_stderr_from_command(command):
    """Run a program on the command line, and save the stdout

    Parameters
    ----------
    command : list
        list of strings to submit to subprocess.run

    Returns
    -------
    lines : list
        Newline-separated strings from output of command
    """
    result = subprocess.run(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout = decode(result.stdout)
    stderr = decode(result.stderr)
    return stdout, stderr
