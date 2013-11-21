import os
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

name = 'pyjip'
version = '0.3'
description = 'JIP pipeline library'
author_email = "thasso.griebel@gmail.com"
url = ""
packages = ['jip', 'jip.cli', 'jip.vendor', 'jip.scripts']
try:
    with open('README.rst') as rf:
        readme = rf.read()
except:
    readme = ''


def default_install():
    try:
        from setuptools import setup, Extension
    except:
        from distribute_setup import use_setuptools
        use_setuptools()
        from setuptools import setup, Extension
    dispatcher_ext = Extension('jip.dispatcher',
                               ['jip/dispatcher/jip_binding.c',
                                'jip/dispatcher/jip_dispatcher.c'])
    setup(
        name=name,
        version=version,
        description=description,
        author_email=author_email,
        url=url,
        license="BSD",
        long_description=readme,
        packages=packages,
        package_data={
            'jip.scripts': ['*.jip']
        },
        install_requires=["sqlalchemy==0.8.2",
                          "jinja2==2.7",
                          "argparse"
                          ],
        ext_modules=[dispatcher_ext],
        entry_points={
            "console_scripts": [
                'jip = jip.cli.jip_main:main'
            ]
        }
    )


def rtd_install():
    from distutils import setup, Extension

    dispatcher_ext = Extension('jip.dispatcher',
                               ['jip/dispatcher/jip_binding.c',
                                'jip/dispatcher/jip_dispatcher.c'])
    setup(
        name=name,
        version=version,
        description=description,
        author_email=author_email,
        url=url,
        license="BSD",
        long_description=readme,
        packages=packages,
        package_data={
            'jip.scripts': ['*.jip']
        },
        ext_modules=[dispatcher_ext],
    )

if on_rtd:
    rtd_install()
else:
    default_install()
