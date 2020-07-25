from pyhdx import VERSION_STRING
from pbr import git
from pbr import version
v = git.get_git_short_sha(r'C:\Users\jhsmi\pp\PyHDX\.git')

print(v)

print(VERSION_STRING)
print(git._run_git_functions())

git_dir = 'C:\\Users\\jhsmi\\pp\\PyHDX\\.git'
tagged = git._run_git_command(
    ['describe'], git_dir,
    throw_on_error=True).replace('-', '.')
print(tagged)

target_version = version.SemanticVersion.from_pip_string(tagged)
print(target_version)

print('whoo')
