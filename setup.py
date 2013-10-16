from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup, Extension


dispatcher_ext = Extension('jip.dispatcher',
                           ['jip/dispatcher/jip_binding.c',
                            'jip/dispatcher/jip_dispatcher.c'])

setup(
    name='pyjip',
    version="0.2",
    description='JIP pipeline library',
    author_email='thasso.griebel@gmail.com',
    url='',
    license="BSD",
    long_description='''This is yet another pipeline library''',
    packages=['jip', 'jip.cli', 'jip.vendor', 'jip.scripts'],
    package_data={
        'jip.scripts': ['*.jip']
    },
    install_requires=["sqlalchemy==0.8.2",
                      "jinja2==2.7",
                      ],
    ext_modules=[dispatcher_ext],
    entry_points={
        "console_scripts": [
            'jip = jip.cli.jip_main:main'
        ]
    }
)
