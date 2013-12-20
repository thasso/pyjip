try:
    from setuptools import setup, Extension
except:
    from distribute_setup import use_setuptools
    use_setuptools()
    from setuptools import setup, Extension

name = 'pyjip'
version = '0.4'
description = 'JIP pipeline library'
author_email = "thasso.griebel@gmail.com"
url = ""
packages = ['jip', 'jip.cli', 'jip.vendor', 'jip.scripts', 'jip.dispatcher']
try:
    with open('Readme.rst') as rf:
        readme = rf.read()
except:
    readme = ''


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
    install_requires=["sqlalchemy>=0.8.2",
                      "jinja2>=2.7",
                      "argparse"
                      ],
    ext_modules=[dispatcher_ext],
    entry_points={
        "console_scripts": [
            'jip = jip.cli.jip_main:main'
        ]
    }
)
