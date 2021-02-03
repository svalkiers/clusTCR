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
        'clustcr': ['input/adaptive_imgt_mapping.csv', 
                    'input/vdjdb_trb.tsv', 
                    'modules/olga/default_models/*/*',
                    'analysis/cq_classifier.pkl']
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
