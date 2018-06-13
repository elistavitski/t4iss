import os
import shutil



def setup_engine():

    here = os.path.dirname(os.path.realpath(__file__))
    user_home = os.path.expanduser('~')

    # defaults
    t4iss_defaults = dict()
    t4iss_defaults['t4iss_data'] = os.path.join(user_home, '/Documents/Work/databases/MP')
    t4iss_defaults['t4iss_scratch'] = os.path.join(user_home, '.t4iss', 'scratch')
    t4iss_defaults['mcr_path'] = os.path.join(here, 'mcr')
    t4iss_defaults['scripts_path'] = os.path.join(here, 'scripts')
    t4iss_defaults['octave_path'] = shutil.which('octave')
    t4iss_defaults['matlab_path'] = shutil.which('matlab')

    print("\nImported t4iss with defaults:\n")


    if (t4iss_defaults['matlab_path'] is not None) and (t4iss_defaults['octave_path'] is not None):
        print('Error: MATLAB/OCTAVE is  cannot be found...')
    else:
        return t4iss_defaults