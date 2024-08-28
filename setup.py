from setuptools import setup, find_packages

setup(
    name='GlycoTools',
    version='0.0.1',
    license='MIT',
    description='ML tools for carbohydrate chemistry',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Natasha Videcrantz Faurschou',
    author_email='natashavidecrantz@gmail.com',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    package_dir={"": "GlycoPredict"},
    packages=setuptools.find_packages(where="GlycoPredict")
    ,
    include_package_data=True,
    extras_require={
        'dev': ['jupyterlab'],
    },
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'GlycoPredict=run_scripts.predict:main',
        ],
    },
)