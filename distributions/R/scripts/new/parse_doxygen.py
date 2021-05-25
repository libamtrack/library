import re
from pprint import pprint
from CppHeaderParser import CppHeader
from sys import argv


param_regex = re.compile(
    r"""
    [@\\]param                       # PROPER parameter classifier
    (?P<mode>(\[(in|out|in,out)])?)  # parameter may be in, out, both or unspecified (def. in)
    \s*
    (?P<name>[a-z_][a-z0-9_]*)       # won't match C99's unicode names
    .*                               # skip anything that's not relevant right now
    (\(array\s*of\s*(length|size)\s*
    (?P<lname>[a-z_][a-z0-9_]*)\))?  # if param is an array, capture it's length
    """, re.VERBOSE | re.MULTILINE | re.IGNORECASE
)


if __name__ == '__main__':
    all_functions = []
    for file in argv[1:]:
        all_functions += CppHeader(file).functions
    pprint(all_functions)
