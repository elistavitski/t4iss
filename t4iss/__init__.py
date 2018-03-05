import os.path

print("Imported t4iss.")
print("Please set dbroot, feff_cmd, and data_dir")

global_cache = dict()

def set_global(name, val):
    ''' set a global variable.'''
    global global_cache
    global_cache[name] = val

def print_globals():
    for key, val in global_cache.items():
        print("{} : {}".format(key, val))


here = os.path.dirname(os.path.realpath(__file__)) + "/.."
dbroot = os.path.join(here, 'data', 'XANES')

global_cache['dbroot'] = dbroot
global_cache['feff_cmd'] = os.path.join(here,'lib','feff_cmd.sh')
global_cache['data_dir'] = os.path.expanduser("~/data")
