"""
Adaptation of Pyviz' nbsite paramdoc: https://github.com/pyviz-dev/nbsite/blob/master/nbsite/paramdoc.py
See Also: https://github.com/Jhsmit/param_docs


"""
import inspect
from functools import partial

try:
    from numpydoc.numpydoc import DEDUPLICATION_TAG
except ImportError:
    DEDUPLICATION_TAG = ''

import param

from param.parameterized import label_formatter

param.parameterized.docstring_signature = False
param.parameterized.docstring_describe_params = False

# Parameter attributes which are never shown
IGNORED_ATTRS = [
    'precedence', 'check_on_set', 'instantiate', 'pickle_default_value',
    'watchers', 'compute_default_fn', 'doc', 'owner', 'per_instance',
    'constant', 'is_instance', 'name', 'allow_None', 'time_fn',
    'time_dependent', 'label', 'inclusive_bounds'
]

# Default parameter attribute values (value not shown if it matches defaults)
DEFAULT_VALUES = {'allow_None': False, 'readonly': False}


def print_lines(app, what, name, obj, options, lines):
    print(lines)


def param_format_basic(app, what, name, obj, options, lines):
    # if what == 'module':
    #     lines = ["start"]

    if what == 'class' and isinstance(obj, param.parameterized.ParameterizedMetaclass):
        lines += ['..', DEDUPLICATION_TAG]  # Prevent numpydoc from mangling the docstring

        try:
            lines.insert(0, '')
            lines.insert(0, f'**{obj.header}**')
        except AttributeError:
            pass

        parameters = ['name']
        mro = obj.mro()[::-1]
        inherited = []
        for cls in mro[:-1]:
            if not issubclass(cls, param.Parameterized) or cls is param.Parameterized:
                continue
            cls_params = [p for p in cls.param if p not in parameters and
                          cls.param[p] == obj.param[p]]
            if not cls_params:
                continue
            parameters += cls_params
            cname = cls.__name__
            module = cls.__module__
            params = ', '.join([f'**{p}**' for p in cls_params])
            inherited.extend(['', '    :class:`{module}.{name}`: {params}'.format(
                name=cname, module=module, params=params)
            ])

        params = [p for p in obj.param if p not in parameters]
        for child in params:
            pobj = obj.param[child]
            if child in ["print_level", "name"]:
                continue
            if pobj.precedence and pobj.precedence < 0:
                continue

            label = label_formatter(pobj.name)
            doc = pobj.doc or ""
            members = inspect.getmembers(pobj)
            params_str = ""
            for m in members:
                # This formats default='apple', objects=['apple', 'pear', 'banana'] etc
                if (m[0][0] != "_" and m[0] not in IGNORED_ATTRS and
                    not inspect.ismethod(m[1]) and not inspect.isfunction(m[1]) and
                    m[1] is not None and DEFAULT_VALUES.get(m[0]) != m[1] and
                    (m[0] != 'label' or pobj.label != label)):
                    subtitutions = {'objects': 'options'}
                    key = subtitutions.get(m[0], m[0])

                    params_str += "%s=%s, " % (key, repr(m[1]))
            params_str = params_str[:-2]
            ptype = pobj.__class__.__name__

            display_name = pobj.label or ' '.join(s.capitalize() for s in child.split('_'))
            if params_str.lstrip():
                lines.extend([f"| **{display_name}** (*{ptype}*, {params_str})"])
            else:
                lines.extend([f"| **{display_name}** (*{ptype}*)"])
            lines.extend(["| " + doc])
            lines.append('')

        if inherited:
            lines.extend(["Additional GUI elements on: "]+inherited)

        lines.append('|')


def param_skip(app, what, name, obj, skip, options):
    if what == 'class' and not skip:
        return (
            getattr(obj, '__qualname__', '').startswith('Parameters.deprecate') or
            (isinstance(obj, partial) and obj.args and isinstance(obj.args[0], param.Parameterized)) or
            (getattr(obj, '__qualname__', '').startswith('Parameterized.') and
             getattr(obj, '__class__', str).__name__ == 'function')
        )
