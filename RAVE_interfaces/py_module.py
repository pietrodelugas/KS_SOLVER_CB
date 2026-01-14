# setup.py
from setuptools import setup, Extension
import os
current_dir = os.path.dirname(__file__)

# Define the extension module
rave_user_events_module = Extension(
    'rave_user_events',                    # Module name
    sources=[os.path.join(current_dir, 'rave_user_events_wrapper.c')],  # Source files
    #libraries=['rave_user_events'],
    extra_compile_args=['-O3', '-I'+os.path.join(current_dir,'.')],
    #extra_link_args=['-L'+os.path.join(current_dir,'../lib'), '-Wl,-rpath='+os.path.join(current_dir,'../lib')]
)

# Setup function to build the module
setup(
    name='rave_user_events',
    version='1.0',
    description='Python interface for the rave_user_events C library',
    ext_modules=[rave_user_events_module]
)

