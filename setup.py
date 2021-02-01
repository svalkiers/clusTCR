from setuptools import setup, find_packages
import versioneer

# conda package requirements
requirements = [
    'numpy',
    'pandas',
    'networkx',
    'scikit-learn',
    'markov_clustering',
    'pyteomics',
    'parmap',
    'faiss'
    # 'tcrdist3'
]

setup(
    name='clusTCR',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="a Python interface for rapid clustering of large sets of CDR3 sequences",
    license="MIT",
    author="Sebastiaan Valkiers & Max Van Houcke",
    author_email='sebastiaan.valkiers@uantwerpen.be',
    url='https://github.com/svalkiers/clustcr',
    packages=find_packages(exclude=('analysis',)),
    package_data={
        'clustcr': ['../data/*', '../data/*/*', 'input/adaptive_imgt_mapping.csv']
    },
    include_package_data=True,
    install_requires=requirements,
    keywords='https://github.com/svalkiers/clustcr',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
