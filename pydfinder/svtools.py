#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

from functools import wraps
from docopt import docopt
from schema import Schema, And, Or, Use, SchemaError

__version__ = '0.1.0'


def argparsed(func):
    @wraps(func)
    def wrapped(argv):
        args = docopt(func.__doc__, argv=argv)
        return func(args)
    return wrapped


@argparsed
def sim(args):
    """
Usage: svtools sim -r =<file0> <tool1> <file1> <tool2> <file2>

Dump audio meta data of the <files>.

Options:
  -r=<file0>      Give the file of true variants
    """
    print "here..."


@argparsed
def diff(args):
    """
Usage: svtools diff -r <file0> <tool1> <file1> <tool2> <file2>

Show the difference between results obtained by two tools.

Options:
  -r <file0>    The file of true variants.
  <tool1>       Tool name should be sprites, lumpy or pindel.
  <file1>       The file of results obtained by <tool1>.
  <tool2>       Tool name should be sprites, lumpy or pindel.
  <file2>       The file of results obtained by <tool2>.
    """

    tools = ('sprites', 'lumpy', 'pindel')
    schema = Schema({
        '-r': os.path.isfile,
        '<tool1>': lambda s: s in tools,
        '<file1>': os.path.isfile,
        '<tool2>': lambda s: s in tools,
        '<file2>': os.path.isfile
    })

    try:
        args = schema.validate(args)
    except SchemaError as e:
        exit(e)

    print args


def help(argv):
    if len(argv) > 1:
        cmd = argv[-1]
        try:
            print(globals()[cmd].__doc__)
        except KeyError:
            exit("%r is not a svtools command. See 'svtools help'." % cmd)
    else:
        docopt(main.__doc__, argv='-h')


def main(argv=None):
    """
A toolbox for research on structural variation detection.

Usage: svtools <command> [<options>...]

General Options:
  -h, --help      Show help.
  --version       Show version and exit.

Commands:
  sim             Generate an artificial genome.
  diff            Show the difference between calls of two tools

See 'svtools help <command>' for more information on a specific command.
    """

    args = docopt(
        main.__doc__,
        version='svtools version %s' % __version__,
        options_first=True,
        argv=argv or sys.argv[1:]
    )

    cmd = args['<command>']
    try:
        method = globals()[cmd]
        assert callable(method)
    except (KeyError, AssertionError):
        exit("%r is not a svtools command. See 'svtools help'." % cmd)

    argv = [args['<command>']] + args['<options>']
    return method(argv)

if __name__ == "__main__":
    main()
