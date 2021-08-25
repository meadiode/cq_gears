'''
CQ_Gears - CadQuery based involute profile gear generator

Copyright 2021 meadiode@github

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''


from setuptools import setup, find_packages

version = '0.45'

setup(
    name='cq_gears',
    version=version,
    
    url='https://github.com/meadiode/cq_gears',
    license='Apache Public License 2.0',
    
    author='meadiode@github',
    author_email='meadiode@protonmail.com',
    
    description='Involute profile gear generator',
    
    long_description=open('./README.md').read(),
    packages=find_packages(),
    install_requires=(['cadquery', 'numpy']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    test_suite='tests',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Development Status :: 3 - Alpha',
        'Environment :: Plugins',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Information Technology',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Multimedia :: Graphics :: 3D Modeling',
        'Topic :: Scientific/Engineering',
    ],
)
