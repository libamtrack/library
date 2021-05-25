from ctypes import *
from collections import defaultdict


mapping = defaultdict(lambda: 'error.type', {
    c_bool: 'logical',
    c_char: 'character',
    POINTER(c_char): 'character',
    c_ubyte: 'raw',
    c_int16: 'integer',
    c_int32: 'integer',
    c_int64: 'integer',
    c_uint16: 'integer',
    c_uint32: 'integer',
    c_uint64: 'integer',
    c_float: 'single',
    c_double: 'double'
})

for ctype in list(mapping.keys()):
    mapping[POINTER(ctype)] = mapping[ctype]
