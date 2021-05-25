import re
import argparse
import ctypes  # used in eval
import warnings
from sys import stderr
from pathlib import Path

from CppHeaderParser import CppHeader, CppMethod, CppVariable

from ctype_translator import mapping


param_regex = re.compile(
    r"""
    [@\\]param                       # PROPER parameter classifier
    (?P<mode>(\[(in|out|in,out)])?)  # parameter may be in, out, both or unspecified (def. in)
    \s*
    (?P<name>[a-z_][a-z0-9_]*)       # won't match C99's unicode names
    .*$                              # also capture the description
    """, re.VERBOSE | re.MULTILINE | re.IGNORECASE
)
array_regex = re.compile(r'array\s*of\s*(length|size)\s*(?P<lname>[a-z0-9_]+)')


def extract_functions_from_file(path):
    try:
        return CppHeader(path).functions
    except IOError:
        raise
    except (TypeError, ValueError, AttributeError) as e:
        raise ValueError(f'{path} is not a valid C/C++ header file') from e


class Parameter:
    def __init__(self, desc: CppVariable, ordinal: int, mode: str, size: str = None):
        self.ord = ordinal
        self.type = mapping[eval(desc['ctypes_type'])]
        self.name = desc['name']
        self.mode = mode[1:-1].lower() if mode != '' else 'in'
        if size is not None:
            try:
                if size.startswith('0x'):
                    self.size = int(size, 16)
                elif size[0] == '0':
                    self.size = int(size, 8)
                else:
                    self.size = int(size)
            except ValueError:
                if size[0].isnumeric():
                    raise ValueError(f'{size} is neither a variable name nor a numeric literal')
                else:
                    self.size = size
        else:
            self.size = None
        self.targets = []  # only used with size parameters

    @property
    def conversion(self):
        return f'\t{self.name} <- as.{self.type}({self.name})'


def create_wrapper_for_function(fun: CppMethod):
    func_name = fun['name'].replace('_', '.')
    signature = [f'{func_name} <- function(']
    center = ['){']
    ending = ['\treturn(AUTO__RETVAL)', '}']
    parameter_list = []
    call_proper = f'\tAUTO___RET__PARAMS <- .C("{fun["name"]}", '
    call_params = []
    before_call = []
    after_call = []

    parameters: dict[str, Parameter] = {}
    size_parameters = set()
    for param_desc, param_doxy, i in zip(
            fun['parameters'],
            param_regex.finditer(fun['doxygen']),
            range(len(fun['parameters']))):
        param_mode = param_doxy.group('mode')
        array_info = array_regex.search(param_doxy.group(0))
        param_size = array_info.group('lname') if array_info else None
        try:
            param = parameters[param_desc['name']] = Parameter(param_desc, i, param_mode, param_size)
        except ValueError as e:
            raise ValueError(
                'Cannot create wrapper for function '
                f'{fun["name"]} (defined in {fun["filename"]}'
            ) from e
        if isinstance(param.size, str):
            size_parameters.add(param.size)

    list_parameters = sorted(parameters.values(), key=lambda x: x.ord)
    out_params = [p for p in list_parameters if p.mode in ['out', 'in,out']]
    if len(out_params) == 0:
        warnings.warn(f'{fun["name"]} appears to have no output parameters (defined in {fun["filename"]}')
        return ''  # no output parameters, no point in creating a wrapper
    for param in size_parameters:
        try:
            if parameters[param].mode == 'in':
                parameters[param].mode = 'size'
        except KeyError as e:
            raise ValueError(
                f'{param} is not a parameter of '
                f'{fun["name"]} (defined in {fun["filename"]}'
            ) from e

    for param in parameters.values():
        if isinstance(param.size, str):
            if param.mode in ['in', 'in,out']:
                parameters[param.size].targets.append(param)

    for param in size_parameters:
        before_call.append(
            f'\t{param} <- min({",".join([f"length({target.name})" for target in parameters[param].targets])})')

    for param in list_parameters:
        call_params.append(param.name)
        if param.mode in ['in', 'in,out']:
            parameter_list.append('\t' + param.name)
            before_call.append(param.conversion)
            if isinstance(param.size, int):
                before_call.append(f'\tstopifnot(length({param.name}) >= {param.size})')
        if param.mode in ['out', 'in,out']:
            if param.mode == 'out':
                before_call.append(
                    f'\t{param.name} <- vector(mode = "'
                    f'{"double" if param.type == "single" else param.type}'
                    f'", length = {param.size if param.size else 1})'
                )
            after_call.append(f'\t{param.name} <- AUTO___RET__PARAMS${param.name}')

    call_proper = call_proper + ', '.join(call_params) + ')'
    if len(out_params) == 1:
        after_call.append(f'\tAUTO__RETVAL <- {out_params[0].name}')
    else:
        after_call.append(f'''\tAUTO__RETVAL <- list({",".join([f'"{p.name}" = {p.name}' for p in out_params])})''')

    return '\n'.join(signature
                     + [',\n'.join(parameter_list)]
                     + center
                     + before_call
                     + [call_proper]
                     + after_call
                     + ending) + '\n\n\n'


def create_wrappers_for_header_file(path: str, out_dir: str):
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_path = Path(out_dir + '/' + Path(path).name.replace('.h', '.R'))
    with open(out_path, 'w') as fout:
        for func in extract_functions_from_file(path):
            try:
                print(create_wrapper_for_function(func), file=fout, end='')
            except Exception as e:
                print(e, file=stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Create R wrapper for functions declared in a header file(s)'
    )
    parser.add_argument('infile', nargs='+',
                        help='C header files to be parsed')
    parser.add_argument('-o', dest='out_dir', metavar='output_dir',
                        help='directory to write wrappers to, (default ./out)',
                        default='./out')
    args = parser.parse_args()
    for infile in args.infile:
        create_wrappers_for_header_file(infile, args.out_dir)


if __name__ == '__main__':
    main()
