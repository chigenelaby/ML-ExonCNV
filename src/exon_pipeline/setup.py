from setuptools import setup, find_packages

setup(
    name="exon_pipeline",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        # 在这里列出你的依赖包
    ],
    entry_points={
        'console_scripts': [
            # 如果有可执行脚本，可以在这里定义
        ],
    },
)
