"""
takes requirements from setup.cfg and makes req-<selection>.txt file
which can be used with `conda install --file req-<selection>.txt

selections are:
all
base
docs
pdf
web

"""

from configparser import ConfigParser
from pathlib import Path

conversions = {
    'torch': 'pytorch'
}

def convert(req_list):
    """
    Apply conversion dictionary and remove null items from requirement list

    Parameters
    ----------
    req_list

    Returns
    -------

    """
    return [conversions.get(item, item) for item in req_list if item]

def make_requirements_file(*extras):
    cp = ConfigParser()
    cp.read_string(Path('setup.cfg').read_text())
    base = convert(cp.get('options', 'install_requires').split('\n'))
    Path('_req-base.txt').write_text('\n'.join(base))

    for extra in extras:
        out = convert(cp.get('options.extras_require', extra).split('\n'))
        Path(f'_req-{extra}.txt').write_text('\n'.join(out))

        base += out

    Path(f'_req-all.txt').write_text('\n'.join(base))

# cp = ConfigParser()
# cp.read_string(Path('setup.cfg').read_text())
#
# web = cp.get('options.extras_require', 'web')
# print(web)

#combined = [i for i in base if i]

#print(combined)

if __name__ == '__main__':
    make_requirements_file('web', 'pdf', 'docs')