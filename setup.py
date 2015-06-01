try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'My Project',
    'author': 'Zhen Zhang',
    'url': 'https://github.com/zhangzhen/dfinder2',
    'download_url': 'Where to download it.',
    'author_email': 'zhangz@csu.edu.cn',
    'version': '0.1',
    'install_requires': ['nose', 'pysam'],
    'packages': ['dfinder2'],
    'scripts': [],
    'name': 'dfinder2'
}

setup(**config)
